      submodule(neklab_linops) exponential_propagator_temp
         implicit none
      contains
         module procedure init_exptA_temp
         logical :: silent_
         character(len=128) :: msg
         character(len=*), parameter :: this_procedure = 'exptA_temp_init'
         silent_ = optval(silent, .false.)
         associate (opts => self%opts)
            ! save minimal tolerances from .par file
            opts%ptol_min_par = param(21)
            opts%vtol_min_par = param(22)
            write(msg,'(A,2(1X,ES15.8))') "Record minimal values for ptol/vtol from .par: ", opts%ptol_min_par, opts%vtol_min_par
            call nek_log_message(msg, this_module, this_procedure)
            
            if (present(bf)) then
               self%baseflow = bf
               call nek_log_message("Update baseflow", this_module, this_procedure)
            end if

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

            if (present(ptol)) then
               opts%ptol = ptol
               param(21) = ptol
               call nek_log_message("Set pressure tolerance", this_module, this_procedure)
               opts%vtol = opts%ptol * opts%p2v_factor
               param(22) = opts%ptol * opts%p2v_factor
               call nek_log_message("Set velocity tolerance", this_module, this_procedure)
            else if (opts%ptol == 0.0_dp) then
               opts%ptol = param(21)
               write(msg,'(A,ES15.8)') "Set pressure tolerance from .par: ", opts%ptol
               opts%vtol = opts%ptol * opts%p2v_factor
               param(22) = opts%ptol * opts%p2v_factor
               write(msg,'(A,ES15.8)') "Set pressure tolerance from .par: ", opts%vtol
               call nek_log_message(msg, this_module, this_procedure)
            end if
            
         ! For the baseflow field for dt/nsteps/cfl computation.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
            call nek_log_information("Set self%baseflow -> vx, vy, vz, pr, t", this_module, this_procedure)
         ! Setup Nek5000 for perturbation solver.
            call setup_linear_solver(solve_baseflow    = .false.,
     &                               solve_temperature = .true.,
     &                               silent            = .true.,
     &                               endtime           = self%tau,
     &                               recompute_dt      = .true.,
     &                               cfl_limit         = opts%cfl_limit,
     &                               vtol              = opts%vtol,
     &                               ptol              = opts%ptol)
            if (.not. silent_) call nek_status(.true.)

            opts%is_set = .true.
         end associate
         end procedure

         module procedure exptA_temp_matvec
         character(len=*), parameter :: this_procedure = 'exptA_temp_matvec'
         integer :: nrst
         if (.not. self%opts%is_set) call nek_stop_error('Operator used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Set nek configuration (after the v[xzy] and v[xyz]p fields are updated)
               call setup_linear_solver(transpose    = .false.,
     &                                  silent       = .true.,
     &                                  endtime      = self%tau,
     &                                  recompute_dt = .true.,
     &                                  cfl_limit    = self%opts%cfl_limit)

      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do
         
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure

         module procedure exptA_temp_rmatvec
         character(len=*), parameter :: this_procedure = 'exptA_temp_rmatvec'
         integer :: nrst
         if (.not. self%opts%is_set) call nek_stop_error('Operator used uninitialized.', this_module, this_procedure)
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Set nek configuration (after the v[xzy] and v[xyz]p fields are updated)
               call setup_linear_solver(transpose    = .true.,
     &                                  silent       = .true.,
     &                                  endtime      = self%tau,
     &                                  recompute_dt = .true.,
     &                                  cfl_limit    = self%opts%cfl_limit)

      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure

         module procedure exptA_temp_compute_rst
            character(len=*), parameter :: this_procedure = 'exptA_temp_compute_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            integer :: irst, itmp
            real(dp) :: rtmp
            select type(vec_out)
            type is (nek_dvector)
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_information(msg, this_module, this_procedure)
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt, silent = .true.)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst
                  call nek_advance()
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, istep - nsteps)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp
               call setup_linear_solver(endtime = time, silent = .true.)
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         end procedure

         module procedure exptA_temp_get_rst
            character(len=*), parameter :: this_procedure = 'exptA_temp_get_rst'
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
