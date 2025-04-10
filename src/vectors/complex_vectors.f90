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
      ! clear restart fields if present
         self%nrst = 0
         end procedure
      
         module procedure nek_zrand
         real(dp) :: alpha
         call nek_drand(self%re); call nek_drand(self%im)
         if (optval(ifnorm, .false.)) then
            alpha = self%norm(); call self%scal(one_cdp/alpha)
         end if
      ! clear restart fields already cleared in nek_drand
         end procedure
      
         module procedure nek_zscal
         integer :: i
         type(nek_zvector) :: wrk
      ! Scratch array.
         wrk = self
      ! Scale complex vector.
         call nek_daxpby(alpha%re, wrk%im, -alpha%im, self%re)
         call nek_daxpby(alpha%re, wrk%re,  alpha%im, self%im)
         do i = 1, self%nrst
            call nek_daxpby(alpha%re, wrk%im_rst(i), -alpha%im, self%re_rst(i))
            call nek_daxpby(alpha%re, wrk%re_rst(i),  alpha%im, self%im_rst(i))
         end do
         end procedure
      
         module procedure nek_zaxpby
         integer :: i
         type(nek_zvector) :: wrk
      ! Scratch array.
         select type (vec)
         type is (nek_zvector)
         
            wrk = vec
      ! Scale vectors before addition.
            call self%scal(beta); call wrk%scal(alpha)
      ! Vector addition.
            call nek_daxpby(1.0_dp, wrk%re, 1.0_dp, self%re)
            call nek_daxpby(1.0_dp, wrk%im, 1.0_dp, self%im)

            do i = 1, self%nrst
               call nek_daxpby(1.0_dp, wrk%re, 1.0_dp, self%re_rst(i))
               call nek_daxpby(1.0_dp, wrk%im, 1.0_dp, self%im_rst(i))
            end do

         class default
            call type_error('vec','nek_zvector','IN',this_module,'nek_zaxpby')
         end select
         end procedure
      
         module procedure nek_zdot
         real(kind=dp) :: alpha_r, alpha_i
         select type (vec)
         type is (nek_zvector)
            alpha_r = self%re%dot(vec%re) + self%im%dot(vec%im)
            alpha_i = self%re%dot(vec%im) - self%im%dot(vec%re)
            alpha = cmplx(alpha_r, alpha_i, kind=dp)
         class default
            call type_error('vec','nek_zvector','IN',this_module,'nek_zdot')
         end select
         end procedure
      
         module procedure nek_zsize
         integer :: lv, m
         lv = nx1*ny1*nz1*nelv
         n = 2*lv + nx2*ny2*nz2*nelv
         if (if3d) n = n + lv
         if (ifto) n = n + lv
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) n = n + lv
            end do
         end if
         end procedure

         module procedure zsave_rst
         integer :: m, lv, lp, torder
         character(len=128) :: msg

         lv = nx1*ny1*nz1*nelv
         lp = nx2*ny2*nz2*nelv
         torder = abs(param(27)) ! integration order in time

         ! sanity checks
         if (irst == torder) then
            write(msg,'(2(A,I0),A)') 'Cannot save rst fields ', torder, ' for a simulation of temporal order ', torder, '.'
            call log_error(msg, this_module, 'zsave_rst')
         else
            write(msg,'(A,I0)') 'Saving rst fields: ', irst
            call log_debug(msg, this_module, 'zsave_rst')
         end if

         select type (vec_rst)
         type is (nek_zvector)
            ! associate?
            call copy(self%re_rst(irst)%vx, vec_rst%re%vx, lv)
            call copy(self%im_rst(irst)%vx, vec_rst%im%vx, lv)
            call copy(self%re_rst(irst)%vy, vec_rst%re%vy, lv)
            call copy(self%im_rst(irst)%vy, vec_rst%im%vy, lv)
            if (if3d) then
               call copy(self%re_rst(irst)%vz, vec_rst%re%vz, lv)
               call copy(self%im_rst(irst)%vz, vec_rst%im%vz, lv)
            end if
            call copy(self%re_rst(irst)%pr, vec_rst%re%pr, lp)
            call copy(self%im_rst(irst)%pr, vec_rst%im%pr, lp)

            if (ifto) then
               call copy(self%re_rst(irst)%theta(:, 1), vec_rst%re%theta(:, 1), lv)
               call copy(self%im_rst(irst)%theta(:, 1), vec_rst%im%theta(:, 1), lv)
            end if
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call copy(self%re_rst(irst)%theta(:, m), vec_rst%re%theta(:, m), lv)
                     call copy(self%im_rst(irst)%theta(:, m), vec_rst%im%theta(:, m), lv)
                  end if
               end do
            end if

         class default
            call type_error('vec_rst','nek_zvector','IN',this_module,'zsave_rst')
         end select
         end procedure
   
         module procedure zget_rst
         integer :: m, lv, lp
         character(len=128) :: msg

         lv = nx1*ny1*nz1*nelv
         lp = nx2*ny2*nz2*nelv

         ! sanity checks
         if (irst < 1) then
            write(msg,'(A,I0)') 'Invalid input for irst: ', irst
            call log_error(msg, this_module, 'zget_rst')
         else if (irst > self%nrst) then
            write(msg,'(A,I0)') 'No rst field to retrieve: ', irst
            call log_warning(msg, this_module, 'zget_rst')
         else
            write(msg,'(A,I0)') 'Retrieving rst fields: ', irst
            call log_information(msg, this_module, 'zget_rst')
         end if

         select type (vec_rst)
         type is (nek_zvector)

            call copy(vec_rst%re%vx, self%re_rst(irst)%vx, lv)
            call copy(vec_rst%im%vx, self%im_rst(irst)%vx, lv)
            call copy(vec_rst%re%vy, self%re_rst(irst)%vy, lv)
            call copy(vec_rst%im%vy, self%im_rst(irst)%vy, lv)
            if (if3d) then
               call copy(vec_rst%re%vz, self%re_rst(irst)%vz, lv)
               call copy(vec_rst%im%vz, self%im_rst(irst)%vz, lv)
            end if
            call copy(vec_rst%re%pr, self%re_rst(irst)%pr, lp)
            call copy(vec_rst%im%pr, self%im_rst(irst)%pr, lp)
            if (ifto) then
               call copy(vec_rst%re%theta(:, 1), self%re_rst(irst)%theta(:, 1), lv)
               call copy(vec_rst%im%theta(:, 1), self%im_rst(irst)%theta(:, 1), lv)
            end if
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call copy(vec_rst%re%theta(:, m), self%re_rst(irst)%theta(:, m), lv)
                     call copy(vec_rst%im%theta(:, m), self%im_rst(irst)%theta(:, m), lv)
                  end if
               end do
            end if

         class default
            call type_error('vec_rst','nek_zvector','OUT',this_module,'zget_rst')
         end select
         end procedure

         module procedure zhas_rst_fields
         has_rst_fields = .false.
         if (self%nrst > 0) has_rst_fields = .true.
         end procedure 

         module procedure zclear_rst_fields
         if (self%nrst == 0) then
            call log_debug('No rst fields to clear', this_module, 'zclear_rst_fields')
         end if
         self%nrst = 0
         end procedure
      
      end submodule
