      submodule (neklab_linops) pressure_projection
         implicit none
      contains
         module procedure apply_DTD
            real(dp), dimension(lv) :: vxtmp, vytmp, vztmp
            ifield = 1
            select type (vec_in)
            type is (nek_pr_dvector)
               select type (vec_out)
               type is (nek_pr_dvector)
                  ! compute wp = D^T @ D @ dpr
                  call compute_LNS_gradp(vxtmp, vytmp, vztmp, vec_in%pr)
                  !call opgradt(vxtmp, vytmp, vztmp, vec_in%pr)
                  call bcdirvc(vxtmp, vytmp, vztmp, v1mask, v2mask, v3mask)
                  call opdiv(vec_out%pr, vxtmp, vytmp, vztmp)
               end select
            end select
         end procedure apply_DTD

         module procedure project_div0
            real(dp), dimension(lv) :: dv1, dv2, dv3
            type(nek_pr_dvector)   :: dpr, x
            type(cg_dp_opts)     :: opts
            type(cg_dp_metadata) :: meta
            integer :: info
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Copy into to nek
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Compute divergence of velocity dp = D^T @ u
                  call opdiv(prp, vxp, vyp, vzp)
      ! Copy from nek to nek_pdvector
                  call nek2pr_vec(dpr, prp)
      ! Solve linear system D^T @ D x = dpr
                  opts = cg_dp_opts(maxiter = 100, if_print_metadata = .true.)
                  meta = cg_dp_metadata()
                  call cg(self%DTD, dpr, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, meta=meta)
      ! Copy solution back to nek
                  call pr_vec2nek(prp, x)
      ! Compute velocity correction
                  call opgradt(dv1, dv2, dv3, prp)
      ! Compute output
                  call opsub2(vxp, vyp, vzp, dv1, dv2, dv3)
      ! Zero out pressure
                  call rzero(prp, lp)
      ! Copty output to vec_out
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
         end procedure project_div0
      end submodule pressure_projection