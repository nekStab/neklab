      submodule(neklab_linops) resolvent
         use LightKrylov, only: Id_rdp, axpby_linop_rdp
         use LightKrylov, only: gmres, fgmres, gmres_dp_opts, gmres_dp_metadata
         use LightKrylov_Logger
         use stdlib_logger, only: information_level
         use neklab_nek_forcing, only: set_neklab_forcing, zero_neklab_forcing
         implicit none
      
         complex(kind=dp), parameter :: zero_cdp = cmplx(0.0_dp, 0.0_dp, kind=dp)
         complex(kind=dp), parameter :: one_cdp = cmplx(1.0_dp, 0.0_dp, kind=dp)
         complex(kind=dp), parameter :: im_cdp = cmplx(0.0_dp, 1.0_dp, kind=dp)
         real(dp), parameter :: pi_dp = 4.0_dp*atan(1.0_dp)
         type(nek_zvector) :: resolvent_forcing
      
      contains
      
         module procedure resolvent_matvec
         type(exptA_linop) :: exptA
         type(nek_dvector) :: b
         real(kind=dp) :: tau
         character(len=*), parameter :: this_procedure = 'resolvent_matvec'
      ! Integration time.
         tau = merge(1.0_dp, 2*pi_dp/abs(self%omega), self%omega == 0.0_dp)
      ! Initialize exponential propagator.
         exptA = exptA_linop(tau, self%baseflow); call exptA%init()
      ! Resolvent computation.
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
      ! Compute forced-response with zero initial condition.
               b = evaluate_rhs(vec_in, self%omega, ifadj)
      ! Solve the linear system for the real part.
               vec_out%re = solve_resolvent_real_part(exptA, b, ifadj)
      ! Evaluate the imaginary part.
               exptA%tau = tau/4.0_dp; call exptA%init()
               vec_out%im = evaluate_imaginary_part(vec_in, self%omega, vec_out%re, ifadj)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module, this_procedure)
         end select
         end procedure
      
         module procedure resolvent_rmatvec
         type(exptA_linop) :: exptA
         type(nek_dvector) :: b
         real(kind=dp) :: tau
         character(len=*), parameter :: this_procedure = 'resolvent_rmatvec'
      ! Force adjoint mode.
      ! Integration time.
         tau = merge(1.0_dp, 2*pi_dp/abs(self%omega), self%omega == 0.0_dp)
      ! Initialize exponential propagator.
         exptA = exptA_linop(tau, self%baseflow); call exptA%init()
      ! Resolvent computation.
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
      ! Compute forced-response with zero initial condition.
               b = evaluate_rhs(vec_in, self%omega, .true.)
      ! Solve the linear system for the real part.
               vec_out%re = solve_resolvent_real_part(exptA, b, .true.)
      ! Evaluate the imaginary part.
               exptA%tau = tau/4.0_dp; call exptA%init()
               vec_out%im = evaluate_imaginary_part(vec_in, self%omega, vec_out%re, .true.)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module, this_procedure)
         end select
         end procedure
      
      !-------------------------------------
      !-----     Utility functions     -----
      !-------------------------------------
      
         function evaluate_rhs(forcing, omega, adjoint) result(b)
            type(nek_zvector), intent(in) :: forcing
            real(dp), intent(in) :: omega
            logical, intent(in) :: adjoint
            type(nek_dvector) :: b
            complex(dp) :: alpha
            real(dp) :: sign
            if (adjoint) then
               ifadj = .true.; call bcast(ifadj, lsize)
            end if
      ! Set the sign (adjoint or direct).
            sign = merge(-1.0_dp, 1.0_dp, adjoint)
      ! Time integration.
            call opzero(vxp, vyp, vzp)
            time = 0.0_dp; lastep = 0; fintim = param(10)
            do istep = 1, nsteps
      ! Update the resolvent forcing for the corresponding time.
               alpha = exp(sign*im_cdp*omega*time)
               call resolvent_forcing%axpby(zero_cdp, forcing, alpha)
               associate (force => resolvent_forcing%re)
                  call set_neklab_forcing(force%vx, force%vy, force%vz, ipert=1)
               end associate
      ! Nek5000 step.
               call nek_advance()
            end do
            call nek2vec(b, vxp, vyp, vzp, prp, tp)
      ! Cleanup
            call zero_neklab_forcing()
            if (adjoint) then
               ifadj = .false.; call bcast(ifadj, lsize)
            end if
         end function
      
         function solve_resolvent_real_part(exptA, b, adjoint) result(x)
            type(nek_dvector), intent(inout) :: b
            type(exptA_linop), intent(in) :: exptA
            logical, intent(in) :: adjoint
            type(nek_dvector) :: x
      ! Local variables.
            type(gmres_dp_opts) :: opts
            integer :: info
            type(axpby_linop_rdp) :: S
            real(dp), parameter :: tol = 1.0e-6_dp
            real(dp) :: alpha
      ! Initialize S = I - exptA.
            S = axpby_linop_rdp(Id_rdp(), exptA, 1.0_dp, -1.0_dp)
            call x%zero()
      ! GMRES call.
            call logger_setup(nio=0, log_level=information_level, log_stdout=.true., log_timestamp=.true.)
            opts = gmres_dp_opts(kdim=64)
            call gmres(S, b, x, info, options=opts, atol=1.0e-12_dp, rtol=tol, transpose=adjoint)
         end function
      
         function evaluate_imaginary_part(forcing, omega, x, adjoint) result(y)
            type(nek_zvector), intent(in) :: forcing
            real(dp), intent(in) :: omega
            type(nek_dvector), intent(in) :: x
            logical, intent(in) :: adjoint
            type(nek_dvector) :: y
            complex(dp) :: alpha
            real(dp) :: sign
            if (adjoint) then
               ifadj = .true.; call bcast(ifadj, lsize)
            end if
      ! Set the sign (adjoint or direct).
            sign = merge(-1.0_dp, 1.0_dp, adjoint)
      ! Initial condition.
            call vec2nek(vxp, vyp, vzp, prp, tp, x)
      ! Time-integration.
            time = 0.0_dp; lastep = 0; fintim = param(10)
            do istep = 1, nsteps
      ! Update resolvent forcing.
               alpha = exp(sign*im_cdp*omega*time)
               call resolvent_forcing%axpby(zero_cdp, forcing, alpha)
               associate (force => resolvent_forcing%re)
                  call set_neklab_forcing(force%vx, force%vy, force%vz, ipert=1)
               end associate
      ! Nek5000 step.
               call nek_advance()
            end do
      ! Cleanup
            call zero_neklab_forcing()
            call nek2vec(y, vxp, vyp, vzp, prp, tp)
            if (adjoint) then
               ifadj = .false.; call bcast(ifadj, lsize)
            end if
         end function
      end submodule
