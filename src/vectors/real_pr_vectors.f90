      submodule(neklab_vectors) real_pr_vectors
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------
      
         module procedure construct_nek_pr_dvector
      ! Pressure.
         if (present(pr)) then
            out%pr = pr
         else
            out%pr = 0.0_dp
         end if
         end procedure
      
      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nek_pr_dzero
         call self%scal(0.0_dp)
         end procedure
      
         module procedure nek_pr_drand
         logical :: normalize
         integer :: i, ieg, iel
         real(kind=dp) :: xl(ldim), fcoeff(3), alpha
         normalize = optval(ifnorm, .false.)
         do i = 1, lp
            ieg = lglel(iel)
            xl(1) = xm2(i, 1, 1, 1)
            xl(2) = ym2(i, 1, 1, 1)
            if (if3D) xl(ldim) = zm2(i, 1, 1, 1)
            call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
            self%pr(i) = self%pr(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
         end do
         if (normalize) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
         end if
         end procedure
      
         module procedure nek_pr_dscal
         call dscal(lp, alpha, self%pr, 1)
         end procedure
      
         module procedure nek_pr_daxpby
         call self%scal(alpha)
         select type (vec)
         type is (nek_pr_dvector)
            call daxpy(lp, beta, vec%pr, 1, self%pr, 1)
         end select
         end procedure
      
         module procedure nek_pr_ddot
         real(kind=dp), external :: glsc2
         select type (vec)
         type is (nek_pr_dvector)
            !alpha = glsc3(self%pr, vec%pr, bm2, lp)
            !alpha = glsc3(self%pr, bm2inv, vec%pr, lp)/volvm2    ! rnorm from convprn in navier1.f
            alpha = glsc2(self%pr, vec%pr, lp)                    ! convprn in navier1.f
         end select
         end procedure
      
         module procedure nek_pr_dsize
         N = lp
         end procedure
      
      end submodule
