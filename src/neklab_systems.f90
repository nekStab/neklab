      module neklab_systems
         ! Fortran standard Library
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_optval, only: optval
         use stdlib_linalg, only: svd, diag
         !---------------------------------------
         !-----     LightKrylov Imports     -----
         !---------------------------------------
         ! Default real kind.
         use LightKrylov, only: dp, qr, apply_inverse_permutation_matrix
         ! Abstract types for real-valued vectors and utilities
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: orthonormalize_basis, orthogonalize_against_basis, zero_basis, copy, rand_basis, linear_combination
         use LightKrylov, only: abstract_system_rdp, abstract_jacobian_linop_rdp
         use LightKrylov_Logger
         use LightKrylov_utils, only: assert_shape
         ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter :: this_module = 'neklab_systems'
         
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         !! Local number of grid points for the pressure mesh.

         !-------------------------------------------
         !-----     LIGHTKRYLOV SYSTEM TYPE   -------
         !-------------------------------------------

         type, extends(abstract_system_rdp), public :: nek_system
         contains
            private
            procedure, pass(self), public :: init => init_DNS
            procedure, pass(self), public :: eval => nonlinear_map
         end type nek_system

         type, extends(abstract_system_rdp), public :: nek_system_upo
         contains
            private
            procedure, pass(self), public :: init => init_DNS_upo
            procedure, pass(self), public :: eval => nonlinear_map_upo
         end type nek_system_upo

         !---------------------------------------------
         !-----     LIGHTKRYLOV JACOBIAN TYPE     -----
         !---------------------------------------------

         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian
         contains
            private
            procedure, pass(self), public :: init    => init_jac_exptA
            procedure, pass(self), public :: matvec  => jac_exptA_matvec
            procedure, pass(self), public :: rmatvec => jac_exptA_rmatvec
         end type nek_jacobian

         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_upo
         contains
            private
            procedure, pass(self), public :: init    => init_jac_map
            procedure, pass(self), public :: matvec  => jac_direct_map
            procedure, pass(self), public :: rmatvec => jac_adjoint_map
         end type nek_jacobian_upo

      contains

         !-------------------------------------------------------
         !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM    -----
         !-------------------------------------------------------

         subroutine init_DNS(self)
            class(nek_system), intent(in) :: self
            call nekgsync()
         
            ! Setup Nek5000 logical flags for perturbative solver.
            ifpert = .false.; ifbase = .false.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
         
            ! Deactivate OIFS.
            if (ifchar) ifchar = .false.
         
            ! Force CFL to 0.5
            if (param(26) > 0.4) then
            if (nid == 0) then
               write (6, *) "WARNING : Target CFL is larger than 0.5"
               write (6, *) "          Forcing it to 0.5"
            end if
            param(26) = 0.4_dp
            end if
         
            ! Compute appropriate step size.
            call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
            dt = param(10)/nsteps; param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0; fintim = nsteps*dt
         
            ! Force contant time step.
            param(12) = -abs(param(12))
         
            ! Broadcast parameters.
            call bcast(param, 200*wdsize)
         
            return
         end subroutine init_DNS

         subroutine nonlinear_map(self, vec_in, vec_out, atol)
            ! Dynamical system.
            class(nek_system),          intent(in)  :: self
            ! Input vector.
            class(abstract_vector_rdp), intent(in)  :: vec_in
            ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! Solver tolerances if needed
            real(dp),                   intent(in)  :: atol

            ! Deactivate linear solver
            ifpert = .false.; call bcast(ifpert, lsize)
            ifbase = .true.; call bcast(ifbase, lsize)

            call self%init()

            ! Set parameters
            param(22) = atol/100         ! pressure tolerance
            param(21) = atol/100         ! velocity tolerance
            param(31) = 0; npert = 0 ! number of perturbations
            lastep = 0
            ! Broadcast parameters.
            call bcast(param, 200*wdsize)

            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  ! Set the initial condition
                  call vec2nek(vx, vy, vz, pr, t, vec_in)

                  ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do

                  ! Extract the final solution to vector.
                  call nek2vec(vec_out, vx, vy, vz, pr, t)
            
                  ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select

            return
         end subroutine nonlinear_map

         !---------------------------------------------------------
         !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN    -----
         !---------------------------------------------------------

         subroutine init_jac_exptA(self)
            class(nek_jacobian), intent(in) :: self
            call nekgsync()
      
            ! Setup Nek5000 logical flags for perturbative solver.
            ifpert = .true.; ifbase = .false.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
      
            ! Force single perturbation mode.
            if (param(31) > 1) then
               if (nid == 0) write (6, *) "neklab does not support (yet) npert > 1."
               call nek_end()
            else
               param(31) = 1; npert = 1
            end if
      
            ! Deactivate OIFS.
            if (ifchar) then
            if (nid == 0) then
               write (6, *) "WARNING : OIFS is not available for linearized solver."
               write (6, *) "          Turning it off."
            end if
            ifchar = .false.
            end if
      
            ! Force CFL to 0.5
            if (param(26) > 0.5) then
            if (nid == 0) then
               write (6, *) "WARNING : Target CFL is larger than 0.5"
               write (6, *) "          Forcing it to 0.5"
            end if
            param(26) = 0.5_dp
            end if
      
            ! Compute appropriate step size.
            call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
            call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
            dt = param(10)/nsteps; param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0; fintim = nsteps*dt
      
            ! Force contant time step.
            param(12) = -abs(param(12))
      
            ! Broadcast parameters.
            call bcast(param, 200*wdsize)

            if (nid .eq. 0) then
               print *, 'init_jac_exptA:'
               print *, '     param(12) = ', param(12)
               print *, '     param(10) = ', param(10)
               print *, '     dt = ', dt
               print *, '     npert = ', npert
            end if
      
            return
         end subroutine init_jac_exptA

         subroutine jac_exptA_matvec(self, vec_in, vec_out)
            class(nek_jacobian),         intent(in)  :: self
            class(abstract_vector_rdp), intent(in)  :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            ! Nek-related setup.
            ifadj = .false.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
            ! Reinitialize Jacobian to incorporate changes in the baseflow
            call self%init()

            ! Set the baseflow initial condition
            call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  ! Set the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
                  !call outpost_dnek(vec_in, 'pr0')
      
                  ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
            
                  ! Extract the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
                  !call outpost_dnek(vec_out, 'pr1')

                  ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
                  !call outpost_dnek(vec_out, 'pr2')
               end select
            end select

            !STOP 9
      
            return
         end subroutine jac_exptA_matvec
      
         subroutine jac_exptA_rmatvec(self, vec_in, vec_out)
            class(nek_jacobian), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! Nek-related setup.
            ifadj = .true.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
            ! Reinitialize Jacobian to incorporate changes in the baseflow
            call self%init()

            ! Set the baseflow initial condition
            call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  ! Set the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
               
                  ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
            
                  ! Extract the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               
                  ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            ifadj = .false.
            return
         end subroutine jac_exptA_rmatvec

         !-----------------------------------------------------------
         !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM_UPO    -----
         !-----------------------------------------------------------

         subroutine init_DNS_upo(self)
            class(nek_system_upo), intent(in) :: self
            call nekgsync()
         
            ! Setup Nek5000 logical flags for perturbative solver.
            ifpert = .false.; ifbase = .false.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
         
            ! Deactivate OIFS.
            if (ifchar) ifchar = .false.
         
            ! Force CFL to 0.5
            if (param(26) > 0.4) then
            if (nid == 0) then
               write (6, *) "WARNING : Target CFL is larger than 0.5"
               write (6, *) "          Forcing it to 0.5"
            end if
            param(26) = 0.4_dp
            end if
         
            ! Compute appropriate step size.
            call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
            dt = param(10)/nsteps; param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0; fintim = nsteps*dt
         
            ! Force contant time step.
            param(12) = -abs(param(12))
         
            ! Broadcast parameters.
            call bcast(param, 200*wdsize)
         
            return
         end subroutine init_DNS_upo

         subroutine nonlinear_map_upo(self, vec_in, vec_out, atol)
            ! Dynamical system.
            class(nek_system_upo),      intent(in)  :: self
            ! Input vector.
            class(abstract_vector_rdp), intent(in)  :: vec_in
            ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! Solver tolerances if needed
            real(dp),                   intent(in)  :: atol

            ! Deactivate linear solver
            ifpert = .false.; call bcast(ifpert, lsize)
            ifbase = .false.; call bcast(ifbase, lsize)

            call self%init()

            ! Set parameters
            param(21) = atol         ! pressure tolerance
            param(21) = atol         ! velocity tolerance
            param(31) = 0; npert = 0 ! number of perturbations
            lastep = 0
            ! Broadcast parameters.
            call bcast(param, 200*wdsize)

            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
                  ! Set the initial condition
                  call ext_vec2nek(vx, vy, vz, pr, t, vec_in)
            
                  ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do

                  ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vx, vy, vz, pr, t)
                  
                  ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select

            return
         end subroutine nonlinear_map_upo

         !-------------------------------------------------------------
         !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN_UPO    -----
         !-------------------------------------------------------------

         subroutine init_jac_map(self)
            class(nek_jacobian_upo), intent(in) :: self
            ! internal
            real(dp), parameter :: cfl_limit = 0.4_dp

            call nekgsync()
      
            ! Setup Nek5000 flags for coupled solver
            ifpert = .true.; ifbase = .true.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
      
            ! Force single perturbation mode.
            if (param(31) > 1) then
               if (nid == 0) write (6, *) "neklab does not (yet) support npert > 1."
               call nek_end()
            else
               param(31) = 1; npert = 1
            end if
      
            ! Deactivate OIFS.
            if (ifchar) then
               if (nid == 0) then
                  write (6, *) "WARNING : OIFS is not available for linearized solver."
                  write (6, *) "          Turning it off."
               end if
               ifchar = .false.
            end if
      
            ! Force CFL to cfl_limit
            ! we are a bit conservative because the baseflow will evolve
            if (param(26) > cfl_limit) then
               if (nid == 0) then
                  write (6, *) "WARNING : Target CFL is larger than", cfl_limit
                  write (6, *) "          Forcing it to", cfl_limit
               end if
               param(26) = cfl_limit
            end if
      
            ! Compute appropriate step size based on the baseflow
            call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
            call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
            dt = param(10)/nsteps; param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0
            
            ! Set integration time
            fintim = nsteps*dt
      
            ! Force contant time step.
            param(12) = -abs(param(12))
      
            ! Broadcast parameters.
            call bcast(param, 200*wdsize)
      
            return
         end subroutine init_jac_map
         
         subroutine jac_direct_map(self, vec_in, vec_out)
            ! Dynamical system.
            class(nek_jacobian_upo),    intent(in)  :: self
            ! Input vector.
            class(abstract_vector_rdp), intent(in)  :: vec_in
            ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            type(nek_ext_dvector) :: vec

            ! Set the baseflow initial condition
            call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
            
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
                  ! Set the perturbation initial condition
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

                  ! Set integration time
                  param(10) = vec_in%T
                  fintim = param(10)

                  ! Reinitialize Jacobian to incorporate changes in the baseflow
                  call self%init()
            
                  ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
            
                  ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)

                  ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
            
                  ! Evaluate f'(X(T), T) * dT and add it to the position residual
                  call compute_fdot(vec)
                  call vec_out%axpby(1.0_dp, vec, vec_in%T)
            
                  ! Evaluate f'(X(0), 0).T @ dx and add phase condition
                  ! Get the initial point of the orbit
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
                  call compute_fdot(vec) 
                  vec_out%T = vec_in%dot(vec)
               end select
            end select
      
            return
         end subroutine jac_direct_map

         subroutine jac_adjoint_map(self, vec_in, vec_out)
            class(nek_jacobian_upo), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            type(nek_ext_dvector) :: vec

            ! Nek-related setup.
            ifadj = .true.; lastep = 0; 
            call bcast(ifadj, lsize)

            ! Set the baseflow initial condition
            call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      
            
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
                  ! Sets the initial condition for Nek5000's linearized solver.
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

                  ! Set integration time
                  param(10) = vec_in%T
                  fintim = param(10)

                  ! Reinitialize Jacobian to incorporate changes in the baseflow & period
                  call self%init()
           
                  ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
                  ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)

                  ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)

                  ! Evaluate f'(X(T), T) * dT and add it to the position residual
                  param(10) = vec_out%T
                  call compute_fdot(vec)
                  call vec_out%axpby(1.0_dp, vec, vec_in%T)

                  ! Evaluate f'(X(0), 0).T @ dx and add phase condition
                  ! Set the initial point of the orbit
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
                  call compute_fdot(vec) 
                  vec_out%T = vec_in%dot(vec)
               end select
            end select         
      
            ifadj = .false.
            return
         end subroutine jac_adjoint_map

         subroutine compute_fdot(vec)
            class(nek_ext_dvector), intent(out) :: vec
            ! internal
            type(nek_ext_dvector) :: vec_in
            ! copy initial condition
            call nek2ext_vec(vec_in, vx, vy, vz, pr, t)
            ! Integrate the nonlinear equations forward in time.
            time = 0.0_dp
            do istep = 1, 1
               call nek_advance()
            end do
            ! Extract f(X(t0+dt))
            call nek2ext_vec(vec, vx, vy, vz, pr, t)
            ! Approximate derivative at t = t0:
            !
            !    f'(X(t0)) ~ ( f(X(t0+dt)) - f(X(t0)) ) / dt
            !
            call vec%sub(vec_in)
            call vec%scal(1.0/dt)
            vec%T = 0.0_dp ! ensure that the period shift vec_out%T is zero
            return
         end subroutine compute_fdot
      
      end module neklab_systems
