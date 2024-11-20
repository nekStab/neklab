      submodule(neklab_systems) fixed_point
         implicit none
      contains
         module procedure nonlinear_map
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Set the initial condition
               call vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances
               call setup_nonlinear_solver(recompute_dt=.true., cfl_limit=0.4_dp,
     $   vtol = atol/10.0, ptol = atol/10.0)
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
         end procedure nonlinear_map
      
         module procedure jac_exptA_matvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow=.false.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp)
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
         end procedure jac_exptA_matvec
      
         module procedure jac_exptA_rmatvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
               call setup_linear_solver(transpose=.true., solve_baseflow=.false.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp)
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
         end procedure jac_exptA_rmatvec
      end submodule
