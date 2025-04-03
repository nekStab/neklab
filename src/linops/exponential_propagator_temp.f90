      submodule(neklab_linops) exponential_propagator_temp
         implicit none
      contains
         module procedure init_exptA_temp
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_message("Set self%baseflow -> vx, vy, vz, pr, t", this_module, "init_exptA_temp")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         end procedure

         module procedure exptA_matvec_temp
         integer :: nrst, itmp
         real(dp) :: rtmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               call setup_linear_solver(transpose     = .false.,
     &                                  silent        = .false.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)

      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst .and. vec_in%has_rst_fields()) then
                     call vec_in%get_rst(vec_rst, istep)
                     call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
                  end if

               end do

      ! Compute restart fields.
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_information(msg, this_module, 'exptA_matvec_temp')
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst

                  call nek_advance()
                  
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp

      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'", 
     & this_module, 'exptA_matvec_temp')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'", 
     & this_module, 'exptA_matvec_temp')
         end select
         end procedure

         module procedure exptA_rmatvec_temp
         integer :: nrst, itmp
         real(dp) :: rtmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               call setup_linear_solver(transpose     = .true.,
     &                                  silent        = .false.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
         ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)

         ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

         ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst .and. vec_in%has_rst_fields()) then
                     call vec_in%get_rst(vec_rst, istep)
                     call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
                  end if

               end do

         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

         ! Compute restart fields.
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_information(msg, this_module, 'exptA_rmatvec_temp')
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst

                  call nek_advance()
                  
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp

            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'", 
     & this_module, 'exptA_rmatvec_temp')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'", 
     & this_module, 'exptA_rmatvec_temp')
         end select
         end procedure
      end submodule
