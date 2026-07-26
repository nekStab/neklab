      submodule(neklab_vectors) real_extended_vectors
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------
         
         module procedure construct_nek_ext_dvector
            character(len=*), parameter :: this_procedure = 'construct_nek_ext_dvector'
            call out%init() ! by default without restart storage
            out%vx = vx
            out%vy = vy
            if (present(vz)) then
               if (if3d) then
                  out%vz = vz
               else
                  call stop_error('vz passed but if3d = .false.', this_module, this_procedure)
               end if
            end if
            if (present(pr)) out%pr = pr
            if (present(theta)) then
               if (ifto) then
                  out%theta = theta
               else
                  call stop_error('theta passed but ifto = .false.', this_module, this_procedure)
               end if
            end if
            if (present(T)) out%T = T
         end procedure

      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nek_ext_dzero
            character(len=*), parameter :: this_procedure = 'nek_ext_dzero'
            call assert_ext_initialized(self, this_procedure)
            self%vx = 0.0_dp
            self%vy = 0.0_dp
            if (if3d) self%vz = 0.0_dp
            self%pr = 0.0_dp
            if (ifto) self%theta = 0.0_dp
            self%T = 0.0_dp
            ! clear restart fields
            self%nrst = 0
            if (self%has_rst) then
               self%vxrst = 0.0_dp
               self%vyrst = 0.0_dp
               if (if3d) self%vzrst = 0.0_dp
               self%prrst = 0.0_dp
               if (ifto) self%thetarst = 0.0_dp
               self%Trst = 0.0_dp
            end if
         end procedure
      
         module procedure nek_ext_drand
            character(len=*), parameter :: this_procedure = 'nek_ext_drand'
            logical :: normalize
            integer :: ix, iy, iz, iel, ieg, ijke, m
            integer :: iface, kx1, kx2, ky1, ky2, kz1, kz2
            real(kind=dp) :: xl(ldim), fcoeff(3), alpha
            call assert_ext_initialized(self, this_procedure)
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
            character(len=*), parameter :: this_procedure = 'nek_ext_dscal'
            integer :: irst, m, lvn, lpn

            call assert_ext_initialized(self, this_procedure)
            lvn = nx1*ny1*nz1*nelv
            lpn = nx2*ny2*nz2*nelv

            call           cmult(self%vx,          alpha, lvn)
            call           cmult(self%vy,          alpha, lvn)
            if (if3d) call cmult(self%vz,          alpha, lvn)
            call           cmult(self%pr,          alpha, lpn)
            if (ifto) call cmult(self%theta(:, 1), alpha, lvn)
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call  cmult(self%theta(:, m), alpha, lvn)
                  end if
               end do
            end if

            self%T = alpha * self%T

            if (self%nrst > 0) call assert_ext_rst_storage(self, this_procedure)
            do irst = 1, self%nrst

               call           cmult(self%vxrst(:, irst),       alpha, lvn)
               call           cmult(self%vyrst(:, irst),       alpha, lvn)
               if (if3d) call cmult(self%vzrst(:, irst),       alpha, lvn)
               call           cmult(self%prrst(:, irst),       alpha, lpn)
               if (ifto) call cmult(self%thetarst(:, irst, 1), alpha, lvn)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) then
                        call  cmult(self%thetarst(:, irst, m), alpha, lvn)
                     end if
                  end do
               end if

               self%Trst(irst) = alpha * self%Trst(irst)

            end do
         end procedure
      
         module procedure nek_ext_daxpby
            character(len=*), parameter :: this_procedure = 'nek_ext_daxpby'
            integer :: irst, m, lvn, lpn

            call assert_ext_initialized(self, this_procedure)
            lvn = nx1*ny1*nz1*nelv
            lpn = nx2*ny2*nz2*nelv

            call self%scal(beta)

            select type (vec)
            type is (nek_ext_dvector)
               call assert_ext_initialized(vec, this_procedure)

               call           add2s2(self%vx,          vec%vx,          alpha, lvn)
               call           add2s2(self%vy,          vec%vy,          alpha, lvn)
               if (if3d) call add2s2(self%vz,          vec%vz,          alpha, lvn)
               call           add2s2(self%pr,          vec%pr,          alpha, lpn)
               if (ifto) call add2s2(self%theta(:, 1), vec%theta(:, 1), alpha, lvn)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) then
                        call  add2s2(self%theta(:, m), vec%theta(:, m), alpha, lvn)
                     end if
                  end do
               end if

               self%T = self%T + alpha*vec%T

               if (self%nrst > 0) call assert_ext_rst_storage(vec, this_procedure)
               do irst = 1, self%nrst

                  call           add2s2(self%vxrst(:, irst),       vec%vxrst(:, irst),       alpha, lvn)
                  call           add2s2(self%vyrst(:, irst),       vec%vyrst(:, irst),       alpha, lvn)
                  if (if3d) call add2s2(self%vzrst(:, irst),       vec%vzrst(:, irst),       alpha, lvn)
                  call           add2s2(self%prrst(:, irst),       vec%prrst(:, irst),       alpha, lpn)
                  if (ifto) call add2s2(self%thetarst(:, irst, 1), vec%thetarst(:, irst, 1), alpha, lvn)
                  if (ldimt > 1) then
                     do m = 2, ldimt
                        if (ifpsco(m - 1)) then
                           call  add2s2(self%thetarst(:, irst, m), vec%thetarst(:, irst, m), alpha, lvn)
                        end if
                     end do
                  end if

                  self%Trst(irst) = self%Trst(irst) + alpha*vec%Trst(irst)

               end do

               self%nrst = max(self%nrst, vec%nrst)

            class default
               call type_error('vec','nek_dvector','IN',this_module,this_procedure)
            end select

         end procedure
      
         module procedure nek_ext_ddot
            character(len=*), parameter :: this_procedure = 'nek_ext_ddot'
            real(kind=dp), external :: glsc3
            integer :: lvn, m

            call assert_ext_initialized(self, this_procedure)
            lvn = nx1*ny1*nz1*nelv

            select type (vec)
            type is (nek_ext_dvector)
               call assert_ext_initialized(vec, this_procedure)

               alpha =                   glsc3(self%vx,          vec%vx,          bm1, lvn)
               alpha =           alpha + glsc3(self%vy,          vec%vy,          bm1, lvn)
               if (if3d) alpha = alpha + glsc3(self%vz,          vec%vz,          bm1, lvn)
               if (ifto) alpha = alpha + glsc3(self%theta(:, 1), vec%theta(:, 1), bm1, lvn)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) then
                        alpha =  alpha + glsc3(self%theta(:, m), vec%theta(:, m), bm1, lvn)
                     end if
                  end do
               end if

               alpha = alpha + self%T * vec%T

            class default
               call type_error('vec','nek_dvector','IN',this_module,this_procedure)
            end select

         end procedure
      
         module procedure nek_ext_dsize
            integer :: lvn, lpn, m
            lvn = nx1*ny1*nz1*nelv
            lpn = nx2*ny2*nz2*nelv
            n = 2*lv + lpn + 1
            if (if3d) n = n + lvn
            if (ifto) n = n + lvn
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) n = n + lvn
               end do
            end if
         end procedure

         ! initialization

         module procedure nek_ext_dinit
            character(len=*), parameter :: this_procedure = 'nek_ext_init'
            logical :: with_rst
            character(len=128) :: msg

            with_rst = optval(alloc_rst, .false.)

            ! --- primary fields ---
            if (allocated(self%vx))    deallocate(self%vx)
            if (allocated(self%vy))    deallocate(self%vy)
            if (allocated(self%vz))    deallocate(self%vz)
            if (allocated(self%pr))    deallocate(self%pr)
            if (allocated(self%theta)) deallocate(self%theta)

            allocate(self%vx(lv), source=0.0_dp)
            allocate(self%vy(lv), source=0.0_dp)
            if (if3d) allocate(self%vz(lv), source=0.0_dp)
            allocate(self%pr(lp), source=0.0_dp)
            if (ifto) allocate(self%theta(lv, ldimt), source=0.0_dp)
            self%T = 0.0_dp

            ! --- restart storage ---
            if (allocated(self%vxrst))    deallocate(self%vxrst)
            if (allocated(self%vyrst))    deallocate(self%vyrst)
            if (allocated(self%vzrst))    deallocate(self%vzrst)
            if (allocated(self%prrst))    deallocate(self%prrst)
            if (allocated(self%thetarst)) deallocate(self%thetarst)
            if (allocated(self%Trst))     deallocate(self%Trst)

            if (with_rst) then
               allocate(self%vxrst(lv, lorder-1), source=0.0_dp)
               allocate(self%vyrst(lv, lorder-1), source=0.0_dp)
               if (if3d) allocate(self%vzrst(lv, lorder-1), source=0.0_dp)
               allocate(self%prrst(lp, lorder-1), source=0.0_dp)
               if (ifto) allocate(self%thetarst(lv, lorder-1, ldimt), source=0.0_dp)
               allocate(self%Trst(lorder-1), source=0.0_dp)
               self%has_rst = .true.
            else
               self%has_rst = .false.
            end if

            self%nrst = 0
            self%initialized = .true.

            write(msg,'(A,L1)') 'nek_ext_dvector initialized, restart storage = ', with_rst
            call log_debug(msg, this_module, this_procedure)
         end procedure

         module procedure nek_ext_dis_initialized
            ok = self%initialized
         end procedure

         module procedure ext_dsave_rst
            character(len=*), parameter :: this_procedure = 'ext_dsave_rst'
            integer :: m, lvn, lpn, torder
            character(len=128) :: msg

            call assert_ext_rst_storage(self, this_procedure)
            lvn = nx1*ny1*nz1*nelv
            lpn = nx2*ny2*nz2*nelv
            torder = abs(param(27)) ! integration order in time
            if (irst < 1 .or. irst > lorder-1) then
               write(msg,'(A,I0)') 'irst out of range: ', irst
               call stop_error(msg, this_module, this_procedure)
            end if

            ! sanity checks
            if (irst == torder) then
               write(msg,'(2(A,I0),A)') 'Cannot save rst fields ', torder, ' for a simulation of temporal order ', torder, '.'
               if (nid == 0) print "('ERROR: ',A)", msg
               call stop_error(msg, this_module, this_procedure)
            else
               write(msg,'(A,I0)') 'Saving rst fields: ', irst
               if (nid == 0) print "('INFO: ',A)", msg
               call log_information(msg, this_module, this_procedure)
            end if

            select type (vec_rst)
            type is (nek_ext_dvector)
               call assert_ext_initialized(vec_rst, this_procedure)

               call           copy(self%vxrst(:, irst),       vec_rst%vx,          lvn)
               call           copy(self%vyrst(:, irst),       vec_rst%vy,          lvn)
               if (if3d) call copy(self%vzrst(:, irst),       vec_rst%vz,          lvn)
               call           copy(self%prrst(:, irst),       vec_rst%pr,          lpn)
               if (ifto) call copy(self%thetarst(:, irst, 1), vec_rst%theta(:, 1), lvn)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) then
                        call  copy(self%thetarst(:, irst, m), vec_rst%theta(:, m), lvn)
                     end if
                  end do
               end if

               self%Trst(irst) = vec_rst%T

            class default
               call type_error('vec_rst','nek_dvector','IN',this_module, this_procedure)
            end select

            self%nrst = max(self%nrst, irst)
         end procedure

         module procedure ext_dget_rst
            character(len=*), parameter :: this_procedure = 'ext_dget_rst'
            integer :: m, lvn, lpn
            character(len=128) :: msg

            call assert_ext_rst_storage(self, this_procedure)
            lvn = nx1*ny1*nz1*nelv
            lpn = nx2*ny2*nz2*nelv

            ! sanity checks
            if (irst < 1) then
               write(msg,'(A,I0)') 'Invalid input for irst: ', irst
               if (nid == 0) print "('ERROR: ',A)", msg
               call stop_error(msg, this_module, this_procedure)
            else if (irst > self%nrst) then
               write(msg,'(A,I0)') 'No rst field to retrieve: ', irst
               if (nid == 0) print "('WARN: ',A)", msg
               call log_warning(msg, this_module, this_procedure)
               return
            else
               write(msg,'(A,I0)') 'Retrieving rst fields: ', irst
               if (nid == 0) print "('INFO: ',A)", msg
               call log_debug(msg, this_module, this_procedure)
            end if

            select type (vec_rst)
            type is (nek_ext_dvector)
               call assert_ext_initialized(vec_rst, this_procedure)

               call           copy(vec_rst%vx,          self%vxrst(:, irst),       lvn)
               call           copy(vec_rst%vy,          self%vyrst(:, irst),       lvn)
               if (if3d) call copy(vec_rst%vz,          self%vzrst(:, irst),       lvn)
               call           copy(vec_rst%pr,          self%prrst(:, irst),       lpn)
               if (ifto) call copy(vec_rst%theta(:, 1), self%thetarst(:, irst, 1), lvn)
               if (ldimt > 1) then
                  do m = 2, ldimt
                     if (ifpsco(m - 1)) call copy(vec_rst%theta(:, m), self%thetarst(:, irst, m), lvn)
                  end do
               end if

               vec_rst%T = self%Trst(irst)

            class default
               call type_error('vec_rst','nek_dvector', 'OUT', this_module, this_procedure)
            end select
         end procedure

         module procedure ext_dhas_rst_fields
            has = merge(.true., .false., self%nrst > 0)
         end procedure

         module procedure ext_dhas_rst_storage
            has = self%has_rst
         end procedure

         module procedure ext_dclear_rst_fields
            character(len=*), parameter :: this_procedure = 'ext_dclear_rst_fields'
            logical :: deallocate_rst
            deallocate_rst = optval(dealloc_rst, .false.)
            if (self%nrst == 0) then
               if (nid == 0) print "('INFO: ',A)", 'No rst fields to clear'
               call log_debug('No rst fields to clear', this_module, this_procedure)
            end if
            self%nrst = 0
            if (deallocate_rst) then
               if (.not. self%has_rst) then
                  if (nid == 0) print "('INFO: ',A)", 'No rst storage to deallocate'
                  call log_debug('No rst storage to deallocate', this_module, this_procedure)
               end if
               if (allocated(self%vxrst))    deallocate(self%vxrst)
               if (allocated(self%vyrst))    deallocate(self%vyrst)
               if (allocated(self%vzrst))    deallocate(self%vzrst)
               if (allocated(self%prrst))    deallocate(self%prrst)
               if (allocated(self%thetarst)) deallocate(self%thetarst)
               if (allocated(self%Trst))     deallocate(self%Trst)
               self%has_rst = .false.
            end if
         end procedure

         ! --- private helpers ---

         subroutine assert_ext_initialized(vec, procedure_name)
            class(nek_ext_dvector), intent(in) :: vec
            character(len=*), intent(in) :: procedure_name
            if (.not. vec%initialized) then
               call stop_error('nek_ext_dvector used before init', this_module, procedure_name)
            end if
         end subroutine

         subroutine assert_ext_rst_storage(vec, procedure_name)
            class(nek_ext_dvector), intent(in) :: vec
            character(len=*), intent(in) :: procedure_name
            if (.not. vec%initialized) then
               call stop_error('nek_ext_dvector used before init', this_module, procedure_name)
            else if (.not. vec%has_rst) then
               call stop_error('nek_ext_dvector initialized without rst storage', this_module, procedure_name)
            end if
         end subroutine
      
      end submodule
