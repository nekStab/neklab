      submodule(neklab_systems) periodic_orbit
         implicit none
         ! Phase condition F_T = < f(x_ref), u - u_ref >: x_ref defines a fixed Poincare section.
         type(nek_ext_dvector) :: xref_phase, fref_phase
         logical :: has_phase_ref = .false.
      contains
         module procedure init_nek_upo_system
         character(len=128) :: msg
         character(len=*), parameter :: this_procedure = 'init_nek_upo_system'
         associate (opts => self%opts)
            ! save minimal tolerances from .par file
            opts%ptol_min_par = param(21)
            opts%vtol_min_par = param(22)
            write(msg,'(A,2(1X,ES15.8))') "Record minimal values for ptol/vtol from .par: ", opts%ptol_min_par, opts%vtol_min_par
            call nek_log_message(msg, this_module, this_procedure)

            if (present(cfl_limit)) then
               opts%cfl_limit = cfl_limit
               param(26) = cfl_limit
               call nek_log_message("Set CFL limit", this_module, this_procedure)
            else if (opts%cfl_limit == 0.0_dp) then
               opts%cfl_limit = param(26)
               write(msg,'(A,F15.8)') "Set CFL limit from .par: ", opts%cfl_limit
               call nek_log_message(msg, this_module, this_procedure)
            end if

            if (present(p2v_factor)) then
               opts%p2v_factor = p2v_factor
               call nek_log_message("Set vtol/ptol", this_module, this_procedure)
            else if (opts%p2v_factor == 0.0_dp) then
               opts%p2v_factor = param(22) / param(21)
               write(msg,'(A,ES15.8)') "Set vtol/ptol from .par: ", opts%p2v_factor
               call nek_log_message(msg, this_module, this_procedure)
            end if
            
            if (present(newton_tol_factor)) then
               opts%newton_tol_factor = newton_tol_factor
               call nek_log_message("Set newton_tol_factor", this_module, this_procedure)
            else if (opts%newton_tol_factor == 0.0_dp) then
               opts%newton_tol_factor = 0.1_dp
               write(msg,'(A,ES15.8)') "Set newton_tol_factor to default: ", opts%newton_tol_factor
               call nek_log_message(msg, this_module, this_procedure)
            end if
            opts%is_set = .true.
         end associate
         has_phase_ref = .false.
         call nek_log_message("Phase reference cleared (set at next call of eval)", this_module, this_procedure)
         end procedure 

         module procedure sync_nek_upo_jacobian
            character(len=*), parameter :: this_procedure = 'sync_nek_upo_jacobian'
            select type (jac => self%jacobian)
            type is (nek_upo_jacobian)
               jac%opts = self%opts
            class default
               call nek_stop_error('nek_temp_system expects a nek_upo_jacobian.', this_module, this_procedure)
            end select
         end procedure

         module procedure nonlinear_map_UPO
         character(len=*), parameter :: this_procedure = 'nonlinear_map_UPO'
         character(len=128) :: msg
         real(dp) :: ptol, vtol
         type(nek_ext_dvector) :: du
         if (.not. self%opts%is_set) call nek_stop_error('System used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)

               ptol = max(self%opts%ptol_min_par, atol * self%opts%newton_tol_factor)
               vtol = max(self%opts%vtol_min_par, atol * self%opts%newton_tol_factor * self%opts%p2v_factor)

      ! Set the initial condition for the nonlinear solver.
               call ext_vec2nek(vx, vy, vz, pr, t, vec_in)

      ! Set appropriate tolerances and Nek status.
               call setup_nonlinear_solver(recompute_dt = .true., 
     &                                     endtime      = vec_in%T,
     &                                     cfl_limit    = self%opts%cfl_limit,
     &                                     ptol         = ptol,
     &                                     vtol         = vtol)
               write (msg, '(A,F9.6)') 'Current period estimate, T = ', vec_in%T
               call nek_log_message(msg, this_module, this_procedure)

      ! Integrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

               end do

      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vx, vy, vz, pr, t)

      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)

      ! Phase condition: F_T = < f(x_ref), u - u_ref >
               if (.not. has_phase_ref) call set_phase_reference(vec_in)
               du = vec_in
               call du%sub(xref_phase)
               vec_out%T = du%dot(fref_phase)

            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_ext_dvector','IN',this_module,this_procedure)
         end select
      ! sync options with jacobian
         call self%sync_jacobian()
         end procedure nonlinear_map_UPO
      
         module procedure jac_direct_map
         character(len=*), parameter :: this_procedure = 'jac_direct_map'
         integer :: nrst
         real(dp) :: atol, vtol, ptol
         type(nek_ext_dvector) :: vec
         if (.not. self%opts%is_set) call nek_stop_error('Jacobian used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)
               ptol = max(self%opts%ptol_min_par, atol * self%opts%newton_tol_factor)
               vtol = max(self%opts%vtol_min_par, atol * self%opts%newton_tol_factor * self%opts%p2v_factor)

      ! Set baseflow.
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Ensure correct nek status -> set end time.
               call setup_linear_solver(solve_baseflow = .true.,
     &                                  transpose      = .false.,
     &                                  recompute_dt   = .true.,
     &                                  endtime        = get_period_abs(self%X),
     &                                  cfl_limit      = self%opts%cfl_limit,
     &                                  ptol           = ptol,
     &                                  vtol           = vtol)

      ! Set the initial condition for the linearized solver.
               call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Intgrate the coupled equations forward
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Evaluate f'(X(T), T) * dT and add it to the position residual
      ! Here we assume that vx,vy,vz contains the endpoint of the nonlinear trajectory
               call compute_fdot(vec)
               call vec_out%axpby(vec_in%T, vec, 1.0_dp)

      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the nonlinear trajectory
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
               call compute_fdot(vec)
      ! Phase row: < f(x_ref), du >  (derivative of F_T)
               if (.not. has_phase_ref) call set_phase_reference(self%X)
               vec_out%T = vec_in%dot(fref_phase)

      ! Reset tolerances
               param(22) = atol
               param(21) = atol

            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_ext_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_direct_map
      
         module procedure jac_adjoint_map
         character(len=*), parameter :: this_procedure = 'jac_adjoint_map'
         integer :: nrst
         real(dp) :: atol, ptol, vtol
         type(nek_ext_dvector) :: vec
         if (.not. self%opts%is_set) call nek_stop_error('Jacobian used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)
               ptol = max(self%opts%ptol_min_par, atol * self%opts%newton_tol_factor)
               vtol = max(self%opts%vtol_min_par, atol * self%opts%newton_tol_factor * self%opts%p2v_factor)
               
      ! Set baseflow.
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Ensure correct nek status -> set end time
               call setup_linear_solver(solve_baseflow = .true.,
     &                                  transpose      = .true.,
     &                                  recompute_dt   = .true.,
     &                                  endtime        = get_period_abs(self%X),
     &                                  cfl_limit      = self%opts%cfl_limit,
     &                                  ptol           = ptol,
     &                                  vtol           = vtol)

      ! Set the initial condition for the linearized solver.
               call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Evaluate f'(X(T), T) * dT and add it to the position residual
               call compute_fdot(vec)
               call vec_out%axpby(vec_in%T, vec, 1.0_dp)

      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the orbit
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
               call compute_fdot(vec)
      
      ! Transpose of the border:
      !   column: f(x_ref) * s        (s = vec_in%T)
      !   row:    < f(X(T)), v >      (vx, vy, vz contain the end point X(T))
               call compute_fdot(vec)
               if (.not. has_phase_ref) call set_phase_reference(self%X)
               call vec_out%axpby(vec_in%T, fref_phase, 1.0_dp)
               vec_out%T = vec_in%dot(vec)

      ! Reset tolerances
               param(22) = atol
               param(21) = atol
            
            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_ext_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_adjoint_map

         module procedure jac_compute_rst
            ! internal
            character(len=*), parameter :: this_procedure = 'jac_compute_rst'
            type(nek_ext_dvector) :: vec_rst
            character(len=128) :: msg
            integer :: irst, itmp
            real(dp) :: rtmp
            select type(vec_out)
            type is (nek_ext_dvector)
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, this_procedure)
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt, silent = .true.)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst
                  call nek_advance()
                  irst = istep - nsteps
                  call nek2ext_vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp
               call setup_linear_solver(endtime = time, silent = .true.)
            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module, this_procedure)
            end select
         end procedure

         module procedure jac_get_rst
            character(len=*), parameter :: this_procedure = 'jac_get_rst'
            type(nek_ext_dvector) :: vec_rst
            character(len=128) :: msg
            select type(vec_in)
            type is (nek_ext_dvector)
               if (vec_in%has_rst_fields()) then
                  call vec_in%get_rst(vec_rst, istep)
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
               end if
            class default
               call type_error('vec_in','nek_ext_dvector','IN',this_module, this_procedure)
            end select
         end procedure

         subroutine set_phase_reference(x)
         ! Store x_ref and f(x_ref) for the phase condition.
         ! Overwrites vx, vy, vz, pr, t (compute_fdot advances one step).
            class(abstract_vector_rdp), intent(in) :: x
            character(len=*), parameter :: this_procedure = 'set_phase_reference'
            select type (x)
            type is (nek_ext_dvector)
               xref_phase = x
               call abs_ext_vec2nek(vx, vy, vz, pr, t, x)
               call compute_fdot(fref_phase)          ! sets fref_phase%T = 0
               has_phase_ref = .true.
               call nek_log_message("Phase reference set", this_module, this_procedure)
            class default
               call type_error('x','nek_ext_dvector','IN',this_module, this_procedure)
            end select
         end subroutine set_phase_reference
      end submodule
