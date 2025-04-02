      submodule(neklab_linops) exponential_propagator
         implicit none
      contains
         module procedure init_exptA
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_message("Set self%baseflow -> vx, vy, vz, pr, t", this_module, "init_exptA")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false., 
     &                            endtime        = self%tau, 
     &                            recompute_dt   = .true., 
     &                            cfl_limit      = 0.5_dp)
         end procedure

         module procedure exptA_matvec
         integer :: nrst, itmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               nrst = abs(param(27)) - 1
               call setup_linear_solver(transpose    = .false., 
     &                                  silent       = .true., 
     &                                  endtime      = self%tau, 
     &                                  recompute_dt = .true., 
     &                                  cfl_limit    = 0.5_dp)
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
      ! Set restart fields if present
                  if (istep <= nrst.and.vec_in%has_rst_fields()) then
                     call get_restart(vec_in, istep)
                  end if
               end do
      ! Copy the final solution to vector and compute restart fields
               call save_restart(vec_out, nrst)
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'exptA_matvec')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'exptA_matvec')
         end select
         end procedure

         module procedure exptA_rmatvec
         integer :: nrst, itmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               nrst = abs(param(27)) - 1
               call setup_linear_solver(transpose    = .true., 
     &                                  silent       = .true.,  
     &                                  endtime      = self%tau,  
     &                                  recompute_dt = .true.,  
     &                                  cfl_limit    = 0.5_dp)
      ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
      ! Set restart fields if present
                  if (istep <= nrst.and.vec_in%has_rst_fields()) then
                     call get_restart(vec_in, istep)
                  end if
               end do
      ! Copy the final solution to vector and compute restart fields
               call save_restart(vec_out, nrst)
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'exptA_rmatvec')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'exptA_rmatvec')
         end select
         end procedure
      
      end submodule
