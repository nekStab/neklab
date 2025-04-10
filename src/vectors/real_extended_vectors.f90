      submodule(neklab_vectors) real_extended_vectors
         implicit none
      contains
      
         module procedure construct_nek_ext_dvector
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
      
      ! Time.
         if (present(time)) then
            out%T = time
         else
            out%T = 0.0_dp
         end if
         end procedure
      
         module procedure nek_ext_dzero
         call self%scal(0.0_dp)
      ! clear restart fields if present
         self%nrst = 0
         end procedure
      
         module procedure nek_ext_drand
         logical :: normalize
         integer :: ix, iy, iz, iel, ieg, ijke, m
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

            if (ifto) then
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%theta(ijke, 1) = self%theta(ijke, 1) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
            end if

            if (ldimt > 1) then
               do m = 2, ldimt
                  call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
                  if (ifpsco(m - 1)) self%theta(ijke, m) = self%theta(ijke, m) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
               end do
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
      
         if (ifto) call dsavg(self%theta(:, 1))
         if (ifto) call bcdirsc(self%theta(:, 1))
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) call dsavg(self%theta(:, m))
               if (ifpsco(m - 1)) call bcdirsc(self%theta(:, m))
            end do
         end if

         call random_number(self%T)
         
         if (normalize) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
         end if
         
      ! clear restart fields if present
         self%nrst = 0
         end procedure
      
         module procedure nek_ext_dscal
         integer :: irst, m, lv, lp

         lv = nx1*ny1*nz1*nelv
         lp = nx2*ny2*nz2*nelv

         call           cmult(self%vx,          alpha, lv)
         call           cmult(self%vy,          alpha, lv)
         if (if3d) call cmult(self%vz,          alpha, lv)
         call           cmult(self%pr,          alpha, lp)
         if (ifto) call cmult(self%theta(:, 1), alpha, lv)
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) then
                  call  cmult(self%theta(:, m), alpha, lv)
               end if
            end do
         end if

         self%T = alpha*self%T

         do irst = 1, self%nrst

            call           cmult(self%vxrst(:, irst),       alpha, lv)
            call           cmult(self%vyrst(:, irst),       alpha, lv)
            if (if3d) call cmult(self%vzrst(:, irst),       alpha, lv)
            call           cmult(self%prrst(:, irst),       alpha, lp)
            if (ifto) call cmult(self%thetarst(:, irst, 1), alpha, lv)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call  cmult(self%thetarst(:, irst, m), alpha, lv)
                  end if
               end do
            end if

            self%Trst(irst) = alpha*self%Trst(irst)

         end do
         end procedure
      
         module procedure nek_ext_daxpby
         integer :: irst, m, lv, lp

         lv = nx1*ny1*nz1*nelv
         lp = nx2*ny2*nz2*nelv

         call self%scal(beta)

         select type (vec)
         type is (nek_ext_dvector)

            call           add2s2(self%vx,          vec%vx,          alpha, lv)
            call           add2s2(self%vy,          vec%vy,          alpha, lv)
            if (if3d) call add2s2(self%vz,          vec%vz,          alpha, lv)
            call           add2s2(self%pr,          vec%pr,          alpha, lp)
            if (ifto) call add2s2(self%theta(:, 1), vec%theta(:, 1), alpha, lv)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call  add2s2(self%theta(:, m), vec%theta(:, m), alpha, lv)
                  end if
               end do
            end if

            self%T = beta*self%T + alpha*vec%T

            do irst = 1, self%nrst

               call           add2s2(self%vxrst(:, irst),       vec%vx,          alpha, lv)
               call           add2s2(self%vyrst(:, irst),       vec%vy,          alpha, lv)
               if (if3d) call add2s2(self%vzrst(:, irst),       vec%vz,          alpha, lv)
               call           add2s2(self%prrst(:, irst),       vec%pr,          alpha, lp)
               if (ifto) call add2s2(self%thetarst(:, irst, 1), vec%theta(:, 1), alpha, lv)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) then
                        call  add2s2(self%thetarst(:, irst, m), vec%theta(:, m), alpha, lv)
                     end if
                  end do
               end if

               self%Trst(irst) = beta*self%Trst(irst) + alpha*vec%Trst(irst)

            end do

            self%nrst = max(self%nrst, vec%nrst)

         class default
            call type_error('vec','nek_dvector','IN',this_module,'nek_ext_daxpby')
         end select

         end procedure
      
         module procedure nek_ext_ddot
         real(kind=dp), external :: glsc3
         integer :: n, m

         n = nx1*ny1*nz1*nelv

         select type (vec)
         type is (nek_ext_dvector)

            alpha =                   glsc3(self%vx,          vec%vx,          bm1, n)
            alpha =           alpha + glsc3(self%vy,          vec%vy,          bm1, n)
            if (if3d) alpha = alpha + glsc3(self%vz,          vec%vz,          bm1, n)
            if (ifto) alpha = alpha + glsc3(self%theta(:, 1), vec%theta(:, 1), bm1, n)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     alpha =  alpha + glsc3(self%theta(:, m), vec%theta(:, m), bm1, n)
                  end if
               end do
            end if

            alpha = alpha + self%T*vec%T

         class default
            call type_error('vec','nek_dvector','IN',this_module,'nek_ext_ddot')
         end select

         end procedure
      
         module procedure nek_ext_dsize
         integer :: lv, m
         lv = nx1*ny1*nz1*nelv
         n = 2*lv + nx2*ny2*nz2*nelv + 1
         if (if3d) n = n + lv
         if (ifto) n = n + lv
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) n = n + lv
            end do
         end if
         end procedure

         module procedure ext_dsave_rst
         integer :: m, lv, lp, torder
         character(len=128) :: msg

         lv = nx1*ny1*nz1*nelv
         lp = nx2*ny2*nz2*nelv
         torder = abs(param(27)) ! integration order in time

         ! sanity checks
         if (irst == torder) then
            write(msg,'(2(A,I0),A)') 'Cannot save rst fields ', torder, ' for a simulation of temporal order ', torder, '.'
            if (nid == 0) print "('ERROR: ',A)", msg
            call log_error(msg, this_module, 'dsave_rst')
         else
            write(msg,'(A,I0)') 'Saving rst fields: ', irst
            if (nid == 0) print "('INFO: ',A)", msg
            call log_information(msg, this_module, 'dsave_rst')
         end if

         select type (vec_rst)
         type is (nek_ext_dvector)

            call           copy(self%vxrst(:, irst),       vec_rst%vx,          lv)
            call           copy(self%vyrst(:, irst),       vec_rst%vy,          lv)
            if (if3d) call copy(self%vzrst(:, irst),       vec_rst%vz,          lv)
            call           copy(self%prrst(:, irst),       vec_rst%pr,          lp)
            if (ifto) call copy(self%thetarst(:, irst, 1), vec_rst%theta(:, 1), lv)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call  copy(self%thetarst(:, irst, m), vec_rst%theta(:, m), lv)
                  end if
               end do
            end if

            self%Trst(irst) = vec_rst%T

         class default
            call type_error('vec_rst','nek_dvector','IN',this_module,'ext_dsave_rst')
         end select
         end procedure

         module procedure ext_dget_rst
         integer :: m, lv, lp
         character(len=128) :: msg

         lv = nx1*ny1*nz1*nelv
         lp = nx2*ny2*nz2*nelv

         ! sanity checks
         if (irst < 1) then
            write(msg,'(A,I0)') 'Invalid input for irst: ', irst
            if (nid == 0) print "('ERROR: ',A)", msg
            call log_error(msg, this_module, 'ext_dget_rst')
         else if (irst > self%nrst) then
            write(msg,'(A,I0)') 'No rst field to retrieve: ', irst
            if (nid == 0) print "('WARN: ',A)", msg
            call log_warning(msg, this_module, 'ext_dget_rst')
         else
            write(msg,'(A,I0)') 'Retrieving rst fields: ', irst
            if (nid == 0) print "('INFO: ',A)", msg
            call log_debug(msg, this_module, 'ext_dget_rst')
         end if

         select type (vec_rst)
         type is (nek_ext_dvector)

            call copy(vec_rst%vx, self%vxrst(:,irst), lv)
            call copy(vec_rst%vy, self%vyrst(:,irst), lv)
            if (if3d) call copy(vec_rst%vz, self%vzrst(:, irst), lv)
            call copy(vec_rst%pr, self%prrst(:, irst), lp)
            if (ifto) call copy(vec_rst%theta(:, 1), self%thetarst(:, irst, 1), lv)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) call copy(vec_rst%theta(:, m), self%thetarst(:, irst, m), lv)
               end do
            end if

            vec_rst%T = self%Trst(irst)

         class default
            call type_error('vec_rst','nek_dvector','OUT',this_module,'ext_dget_rst')
         end select
         end procedure

         module procedure ext_dhas_rst_fields
         has_rst_fields = .false.
         if (self%nrst > 0) has_rst_fields = .true.
         end procedure  

         module procedure ext_dclear_rst_fields
         if (self%nrst == 0) then
            if (nid == 0) print "('INFO: ',A)", 'No rst fields to clear'
            call log_debug('No rst fields to clear', this_module, 'ext_dclear_rst_fields')
         end if
         self%nrst = 0
         end procedure
      
      end submodule
