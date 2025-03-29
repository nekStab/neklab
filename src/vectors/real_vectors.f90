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
         integer :: ix, iy, iz, iel, ieg, ijke
         integer :: iface, kx1, kx2, ky1, ky2, kz1, kz2
         real(kind=dp) :: xl(ldim), fcoeff(3), alpha
         normalize = optval(ifnorm, .false.)
      
         ifield = 1 ! for bcdirvc

         do iel = 1, nelv
         do iz = 1, lz1
         do iy = 1, ly1
         do ix = 1, lx1
            ieg = lglel(iel)
            xl(1) = xm1(ix, iy, iz, iel)
            xl(2) = ym1(ix, iy, iz, iel)
            if (if3d) xl(3) = zm1(ix, iy, iz, iel)
            ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(iel-1)))
      
            call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
            self%vx(ijke) = self%vx(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
      
            call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
            self%vy(ijke) = self%vy(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)

            if (if3d) then
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%vz(ijke) = self%vz(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
            end if
         end do
         end do
         end do
         end do
      
      ! Face averaging.
         call opdssum(self%vx, self%vy, self%vz)
         call opcolv(self%vx, self%vy, self%vz, vmult)
         call dsavg(self%vx)
         call dsavg(self%vy)
         if (if3d) call dsavg(self%vz)
         call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
         if (normalize) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
         end if
         end procedure
      
         module procedure nek_dscal
         integer :: n1, n2
         n1 = nx1*ny1*nz1*nelv
         n2 = nx2*ny2*nz2*nelv
         call cmult(self%vx, alpha, n1)
         call cmult(self%vy, alpha, n1)
         if (if3d) call cmult(self%vz, alpha, n1)
         call cmult(self%pr, alpha, n2)
         if (ifto) call cmult(self%theta(:, 1), alpha, n1)
         end procedure
      
         module procedure nek_daxpby
         integer :: n1, n2
         n1 = nx1*ny1*nz1*nelv
         n2 = nx2*ny2*nz2*nelv
         call self%scal(beta)
         select type (vec)
         type is (nek_dvector)
            call add2s2(self%vx, vec%vx, alpha, n1)
            call add2s2(self%vy, vec%vy, alpha, n1)
            if (if3d) call add2s2(self%vz, vec%vz, alpha, n1)
            call add2s2(self%pr, vec%pr, alpha, n2)
            if (ifto) call add2s2(self%theta(:, 1), vec%theta(:, 1), alpha, n1)
         class default
            call stop_error("The intent [IN] argument 'vec' must be of type 'nek_dvector'",
     & this_module, 'nek_daxpby')
         end select

         end procedure
      
         module procedure nek_ddot
         real(kind=dp), external :: glsc3
         integer :: i, n
         n = nx1*ny1*nz1*nelv
         select type (vec)
         type is (nek_dvector)
            alpha =         glsc3(self%vx, vec%vx, bm1, n)
            alpha = alpha + glsc3(self%vy, vec%vy, bm1, n)
            if (if3d) alpha = alpha + glsc3(self%vz, vec%vz, bm1, n)
            if (ifto) then
               alpha = alpha + glsc3(self%theta(:, 1), vec%theta(:, 1), bm1, n)
            end if
            if (ldimt > 1) then
            do i = 2, ldimt
               if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), vec%theta(:, i), bm1, n)
            end do
            end if
         class default
            call stop_error("The intent [IN] argument 'vec' must be of type 'nek_dvector'",
     & this_module, 'nek_ddot')
         end select

         end procedure
      
         module procedure nek_dsize
         integer :: i, n1
         n1 = nx1*ny1*nz1*nelv
         n = 2*n1 + nx2*ny2*nz2*nelv
         if (if3d) n = n + n1
         if (ifto) n = n + n1
         if (ldimt > 1) then
         do i = 2, ldimt
            if (ifpsco(i - 1)) n = n + n1
         end do
         end if
         end procedure
      
      end submodule
