      submodule(neklab_vectors) real_vectors
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------
      
         module procedure construct_nek_dvector
      ! Velocity arrays.
         out%vx = vx; out%vy = vy
         if (present(vz)) then
            out%vz = vz
         else
            out%vz = 0.0_dp
         end if
      
      ! Pressure.
         if (present(pr)) then
            out%pr = pr
         else
            out%pr = 0.0_dp
         end if
      
      ! Temperature and passive scalars.
         if (present(theta)) then
            out%theta = theta
         else
            out%theta = 0.0_dp
         end if
         end procedure
      
      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nek_dzero
         call self%scal(0.0_dp)
         end procedure
      
         module procedure nek_drand
         logical :: normalize
         integer :: i, n, ieg, iel
         real(kind=dp) :: xl(ldim), fcoeff(3), alpha
      
         n = nx1*ny1*nz1*nelv
      
         call random_number(self%vx)
         call random_number(self%vy)
         if (if3d) call random_number(self%vz)
      
         call opdssum(self%vx, self%vy, self%vz)
         call opcolv(self%vx, self%vy, self%vz, vmult)
         call dsavg(self%vx); call dsavg(self%vy)
         if (if3d) call dsavg(self%vz)
      
         call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
         if (optval(ifnorm, .false.)) then
            alpha = self%norm(); call self%scal(1.0_dp/alpha)
         end if
         end procedure
      
         module procedure nek_dscal
         call dscal(lv, alpha, self%vx, 1)
         call dscal(lv, alpha, self%vy, 1)
         if (if3d) call dscal(lv, alpha, self%vz, 1)
         call dscal(lp, alpha, self%pr, 1)
         if (ifto) call dscal(lv, alpha, self%theta(:, 1), 1)
         end procedure
      
         module procedure nek_daxpby
         call self%scal(alpha)
         select type (vec)
         type is (nek_dvector)
            call daxpy(lv, beta, vec%vx, 1, self%vx, 1)
            call daxpy(lv, beta, vec%vy, 1, self%vy, 1)
            if (if3d) call daxpy(lv, beta, vec%vz, 1, self%vz, 1)
            call daxpy(lp, beta, vec%pr, 1, self%pr, 1)
            if (ifto) call daxpy(lv, beta, vec%theta(:, 1), 1, self%theta(:, 1), 1)
         end select
         end procedure
      
         module procedure nek_ddot
         real(kind=dp), external :: op_glsc2_wt, glsc3
         integer :: i
      
         ifield = 1
         select type (vec)
         type is (nek_dvector)
            alpha = op_glsc2_wt(self%vx, self%vy, self%vz, vec%vx, vec%vy, vec%vz, bm1)
            if (ifto) then
               alpha = alpha + glsc3(self%theta(:, 1), vec%theta(:, 1), bm1, lv)
            end if
            if (ldimt > 1) then
            do i = 2, ldimt
               if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), vec%theta(:, i), bm1, lv)
            end do
            end if
         end select
         end procedure
      
         module procedure nek_dsize
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
