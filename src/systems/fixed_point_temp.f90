      submodule(neklab_systems) fixed_point_temp
         implicit none
      contains
         module procedure init_nek_temp_system
         character(len=128) :: msg
         character(len=*), parameter :: this_procedure = 'init_nek_temp_system'
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
         end procedure

         module procedure sync_nek_temp_jacobian
            character(len=*), parameter :: this_procedure = 'sync_nek_temp_jacobian'
            select type (jac => self%jacobian)
            type is (nek_temp_jacobian)
               jac%opts = self%opts
            class default
               call nek_stop_error('nek_temp_system expects a nek_temp_jacobian.', this_module, this_procedure)
            end select
         end procedure

         module procedure nonlinear_map_temp
         real(dp) :: vtol, ptol
         character(len=*), parameter :: this_procedure = 'nonlinear_map_temp'
         if (.not. self%opts%is_set) call nek_stop_error('System used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
         
               ptol = max(self%opts%ptol_min_par, atol * self%opts%newton_tol_factor)
               vtol = max(self%opts%vtol_min_par, atol * self%opts%newton_tol_factor * self%opts%p2v_factor)

      ! Set the initial condition for the nonlinear solver.
               call vec2nek(vx, vy, vz, pr, t, vec_in)

      ! Set appropriate tolerances.
               call setup_nonlinear_solver(recompute_dt = .true., 
     &                                     cfl_limit    = self%opts%cfl_limit,
     &                                     ptol         = ptol,
     &                                     vtol         = vtol)

      ! Integrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)

      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
      ! sync options with jacobian
         call self%sync_jacobian()
         end procedure nonlinear_map_temp
      
         module procedure jac_exptA_temp_matvec
         character(len=*), parameter :: this_procedure = 'jac_exptA_temp_matvec'
         integer :: nrst
         real(dp) :: atol, ptol, vtol
         if (.not. self%opts%is_set) call nek_stop_error('Jacobian used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)
               ptol = max(self%opts%ptol_min_par, atol * self%opts%newton_tol_factor)
               vtol = max(self%opts%vtol_min_par, atol * self%opts%newton_tol_factor * self%opts%p2v_factor)

      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = self%opts%cfl_limit,
     &                                  ptol           = ptol,
     &                                  vtol           = vtol)

      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Extract the final solution to vector.   
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Reset tolerances.
               param(22) = atol
            
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_exptA_temp_matvec
      
         module procedure jac_exptA_temp_rmatvec
      ! internal
         character(len=*), parameter :: this_procedure = 'jac_exptA_temp_rmatvec'
         integer :: nrst
         real(dp) :: atol, ptol, vtol
         if (.not. self%opts%is_set) call nek_stop_error('Jacobian used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)
               ptol = max(self%opts%ptol_min_par, atol * self%opts%newton_tol_factor)
               vtol = max(self%opts%vtol_min_par, atol * self%opts%newton_tol_factor * self%opts%p2v_factor)

      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status.
               call setup_linear_solver(transpose      = .true., 
     &                                  solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = self%opts%cfl_limit,
     &                                  ptol           = ptol,
     &                                  vtol           = vtol)

      ! Set the initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Reset tolerances.
               param(22) = atol
               
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_exptA_temp_rmatvec

         module procedure jac_exptA_temp_compute_rst
            ! internal
            character(len=*), parameter :: this_procedure = 'jac_exptA_temp_compute_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            integer :: irst, itmp
            real(dp) :: rtmp
            select type(vec_out)
            type is (nek_dvector)
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
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp
               call setup_linear_solver(endtime = time, silent = .true.)
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         end procedure

         module procedure jac_exptA_temp_get_rst
            character(len=*), parameter :: this_procedure = 'jac_exptA_temp_get_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            select type(vec_in)
            type is (nek_dvector)
               if (vec_in%has_rst_fields()) then
                  call vec_in%get_rst(vec_rst, istep)
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
               end if
            class default
               call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
            end select
         end procedure
      end submodule
