      submodule(neklab_systems) fixed_point
         implicit none
      contains
         module procedure nonlinear_map
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

      ! Set appropriate tolerances.
               call setup_nonlinear_solver(recompute_dt = .true., 
     &                                     cfl_limit    = 0.4_dp,
     &                                     vtol         = atol*0.1,
     &                                     ptol         = atol*0.1)
            
      ! Set the initial condition for the nonlinear solver.
               call vec2nek(vx, vy, vz, pr, t, vec_in)
      
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
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'nonlinear_map')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'nonlinear_map')
         end select
         end procedure nonlinear_map
      
         module procedure jac_exptA_matvec
      ! internal
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

      ! Ensure correct nek status.
               atol = param(22)
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)
            
      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Set the initial condition for the linearized solver.
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

      ! Reset tolerances.
               param(22) = atol

            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'jac_exptA_matvec')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'jac_exptA_matvec')
         end select
         end procedure jac_exptA_matvec
      
         module procedure jac_exptA_rmatvec
      ! internal
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

      ! Ensure correct nek status.
               atol = param(22)
               call setup_linear_solver(transpose      = .true., 
     &                                  solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)

      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Set the initial condition for the linearized solver.

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

      ! Reset tolerances.
               param(22) = atol

            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'jac_exptA_rmatvec')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'jac_exptA_rmatvec')
         end select
         end procedure jac_exptA_rmatvec
      end submodule
