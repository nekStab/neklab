      submodule(neklab_vectors) complex_vectors
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------
      
         module procedure construct_nek_zvector
      ! Definition of the constructor cause an Internal Compiler Error.
      ! Implementation is delayed for now.
         out%re%vx = vx(:, 1); out%im%vx = vx(:, 2)
         out%re%vy = vy(:, 1); out%im%vy = vy(:, 2)
         if (present(vz)) then
            out%re%vz = vz(:, 1); out%im%vz = vz(:, 2)
         else
            out%re%vz = 0.0_dp; out%im%vz = 0.0_dp
         end if
         if (present(pr)) then
            out%re%pr = pr(:, 1); out%im%pr = pr(:, 2)
         else
            out%re%pr = 0.0_dp; out%im%pr = 0.0_dp
         end if
         if (present(theta)) then
            out%re%theta = theta(:, :, 1); out%im%theta = theta(:, :, 2)
         else
            out%re%theta = 0.0_dp; out%im%theta = 0.0_dp
         end if
         end procedure
      
      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nek_zzero
         call self%scal(zero_cdp)
         end procedure
      
         module procedure nek_zrand
         real(dp) :: alpha
         call nek_drand(self%re); call nek_drand(self%im)
         if (optval(ifnorm, .false.)) then
            alpha = self%norm(); call self%scal(one_cdp/alpha)
         end if
         end procedure
      
         module procedure nek_zscal
         type(nek_zvector) :: wrk
      ! Scratch array.
         wrk = self
      ! Scale complex vector.
         call nek_daxpby(self%re, alpha%re, wrk%im, -alpha%im)
         call nek_daxpby(self%im, alpha%re, wrk%re, alpha%im)
         end procedure
      
         module procedure nek_zaxpby
         type(nek_zvector) :: wrk
      ! Scratch array.
         select type (vec)
         type is (nek_zvector)
            wrk = vec
      ! Scale vectors before addition.
            call self%scal(alpha); call wrk%scal(beta)
      ! Vector addition.
            call nek_daxpby(self%re, 1.0_dp, wrk%re, 1.0_dp)
            call nek_daxpby(self%im, 1.0_dp, wrk%im, 1.0_dp)
         end select
         end procedure
      
         module procedure nek_zdot
         real(kind=dp) :: alpha_r, alpha_i
         select type (vec)
         type is (nek_zvector)
            alpha_r = self%re%dot(vec%re) + self%im%dot(vec%im)
            alpha_i = self%re%dot(vec%im) - self%im%dot(vec%re)
            alpha = cmplx(alpha_r, alpha_i, kind=dp)
         end select
         end procedure
      
         module procedure nek_zsize
         integer :: i
         n = 2*lv + lp
         if (if3d) n = n + lv
         if (ifto) n = n + lv
         if (ldimt > 1) then
         do i = 2, ldimt
            if (ifpsco(i - 1)) n = n + lv
         end do
         end if
         end procedure
      
      end submodule
