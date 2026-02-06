      submodule(neklab_systems) fixed_point_temp
         implicit none
      contains
         module procedure nonlinear_map_temp
         character(len=*), parameter :: this_procedure = 'nonlinear_map_temp'
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
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure nonlinear_map_temp
      
         module procedure jac_exptA_temp_matvec
      ! internal
         character(len=*), parameter :: this_procedure = 'jac_exptA_temp_matvec'
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
            
      ! Ensure correct nek status
               atol = param(22)
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)
      
      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

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
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_exptA_temp_rmatvec
      end submodule
