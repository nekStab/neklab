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
            call setup_linear_solver(transpose=.false., silent=.true.)
            ! Force baseflow.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)

            ! Set initial condition for the linearized solver.
            select type(vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select

            ! Integrate the equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            enddo

            ! Copy the final solution to vector.
            select type(vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end procedure

         module procedure exptA_rmatvec
            call setup_linear_solver(transpose=.true., silent=.true.)
            ! Force baseflow.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)

            ! Set initial condition for the linearized solver.
            select type(vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select

            ! Integrate the equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            enddo

            ! Copy the final solution to vector.
            select type(vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end procedure
      end submodule
