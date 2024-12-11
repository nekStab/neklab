      submodule(neklab_linops) exponential_propagator
         implicit none
      contains
         module procedure init_exptA
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call logger%log_message("Set self%baseflow -> vx, vy, vz, pr, t", module=this_module, procedure="init_exptA")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow=.false., endtime=self%tau, recompute_dt=.true., cfl_limit=0.5_dp)
         end procedure
      
         module procedure exptA_matvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               call setup_linear_solver(transpose=.false., silent=.false.)
      ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call stop_error('Output must be a nek_dvector', module=this_module, procedure='exptA_matvec')
            end select
         class default
            call stop_error('Input must be a nek_dvector', module=this_module, procedure='exptA_matvec')
         end select
         end procedure
      
         module procedure exptA_rmatvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               call setup_linear_solver(transpose=.true., silent=.true.)
      ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call stop_error('Output must be a nek_dvector', module=this_module, procedure='exptA_rmatvec')
            end select
         class default
            call stop_error('Input must be a nek_dvector', module=this_module, procedure='exptA_rmatvec')
         end select
         end procedure
      end submodule
