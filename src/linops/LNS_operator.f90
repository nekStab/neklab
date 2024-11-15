      submodule(neklab_linops) LNS_operator
         implicit none
      contains
         module procedure init_LNS
      ! Force the baseflow field for dt/nsteps/clf computation
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call logger%log_message('Set self%baseflow -> vx, vy, vz, pr, t',
     $   module = this_module, procedure = 'init_LNS')
      ! Setup Nek5000 for perturbation solver
         call setup_linear_solver(solve_baseflow=.false., recompute_dt=.true.,
     $   cfl_limit = 0.5_dp)
         end procedure init_LNS
      
         module procedure LNS_matvec
         real(dp), dimension(lv, 1) :: Lux, Luy, Luz, dv1, dv2, dv3
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Force the baseflow field.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Ensure correct settings for Nek
               call setup_linear_solver(transpose=.false., silent=.true.)
      ! Sets the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Apply LNS operator
               call apply_Lv(Lux, Luy, Luz, vxp, vyp, vzp, trans=.false.)
      ! Compute divergence of velocity dp = D^T @ u
               call opdiv(prp, Lux, Luy, Luz)
      ! Solve for pressure correction to enforce continuity D^T @ D @ dp = D^T @ u
               call project_perturbation(prp) ! in-place
      ! Compute corresponding velocity correction
               call opgradt(dv1, dv2, dv3, prp)
      ! Compute final velocity
               call opsub3(vxp, vyp, vzp, Lux, Luy, Luz, dv1, dv2, dv3)
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end select
         end procedure LNS_matvec
      
         module procedure LNS_rmatvec
         real(dp), dimension(lv, 1) :: Lux, Luy, Luz, dv1, dv2, dv3
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Force the baseflow field.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Ensure the appropriate solver state
               call setup_linear_solver(transpose=.true., silent=.true.)
      ! Sets the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Apply adjoint LNS operator
               call apply_Lv(Lux, Luy, Luz, vxp, vyp, vzp, trans=.true.)
      ! Compute divergence of velocity dp = D^T @ u
               call opdiv(prp, Lux, Luy, Luz)
      ! Solve for pressure correction to enforce continuity D^T @ D @ dp = D^T @ u
               call project_perturbation(prp) ! in-place
      ! Compute corresponding velocity correction
               call opgradt(dv1, dv2, dv3, prp)
      ! Compute final velocity
               call opsub3(vxp, vyp, vzp, Lux, Luy, Luz, dv1, dv2, dv3)
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end select
         end procedure LNS_rmatvec
      end submodule LNS_operator
