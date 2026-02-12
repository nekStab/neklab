      submodule(neklab_linops) exponential_propagator_temp
         implicit none
      contains
         module procedure init_exptA_temp
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_information("Set self%baseflow -> vx, vy, vz, pr, t", this_module, "init_exptA_temp")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         end procedure

         module procedure exptA_temp_matvec
         character(len=*), parameter :: this_procedure = 'exptA_temp_matvec'
         integer :: nrst
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
     &                                  cfl_limit    = 0.5_dp)

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
     &                                  cfl_limit    = 0.5_dp)

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
               call nek_log_debug(msg, this_module, this_procedure)
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt)
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
