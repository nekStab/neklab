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
               call nek_daxpby(1.0_dp, wrk%re_rst(i), 1.0_dp, self%re_rst(i))
               call nek_daxpby(1.0_dp, wrk%im_rst(i), 1.0_dp, self%im_rst(i))
            end do
            self%nrst = max(self%nrst, wrk%nrst)
         class default
            call stop_error("The intent [IN] argument 'vec' must be of type 'nek_zvector'",
     & this_module, 'nek_zaxpby')
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
            call stop_error("The intent [IN] argument 'vec' must be of type 'nek_zvector'",
     & this_module, 'nek_zdot')
         end select
         end procedure
      
         module procedure nek_zsize
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

         module procedure zsave_rst
         integer :: i, n1, n2
         integer :: irst, torder
         character(len=128) :: msg
         n1 = nx1*ny1*nz1*nelv
         n2 = nx2*ny2*nz2*nelv
         torder = abs(param(27)) ! integration order in time
         ! sanity check
         if (self%nrst == torder - 1) then
            write(msg,'(2(A,I0),A)') 'Cannot save rst fields ', torder, ' for a simulation of temporal order ', torder, '.'
            call logger%log_error(msg, this_module, 'zsave_rst')
         else
            self%nrst = self%nrst + 1
            write(msg,'(A,I0)') 'Saving rst fields: ', self%nrst
            call logger%log_information(msg, this_module, 'zsave_rst')
         end if
         irst = self%nrst
         select type (vec_rst)
         type is (nek_zvector)
            ! associate?
            call copy(self%re_rst(irst)%vx, vec_rst%re%vx, n1)
            call copy(self%im_rst(irst)%vx, vec_rst%im%vx, n1)
            call copy(self%re_rst(irst)%vy, vec_rst%re%vy, n1)
            call copy(self%im_rst(irst)%vy, vec_rst%im%vy, n1)
            if (if3d) then
               call copy(self%re_rst(irst)%vz, vec_rst%re%vz, n1)
               call copy(self%im_rst(irst)%vz, vec_rst%im%vz, n1)
            end if
            call copy(self%re_rst(irst)%pr, vec_rst%re%pr, n2)
            call copy(self%im_rst(irst)%pr, vec_rst%im%pr, n2)
            if (ifto) then
               call copy(self%re_rst(irst)%theta(:, 1), vec_rst%re%theta(:, 1), n1)
               call copy(self%im_rst(irst)%theta(:, 1), vec_rst%im%theta(:, 1), n1)
               if (ldimt > 1) then
                  do i = 2, ldimt
                     if (ifpsco(i - 1)) then
                        call copy(self%re_rst(irst)%theta(:, i), vec_rst%re%theta(:, i), n1)
                        call copy(self%im_rst(irst)%theta(:, i), vec_rst%im%theta(:, i), n1)
                     end if
                  end do
               end if
            end if
         class default
            call stop_error("The intent [IN] argument 'vec_rst' must be of type 'nek_zvector'",
     & this_module, 'zsave_rst')
         end select
         end procedure
   
         module procedure zget_rst
         integer :: i, n1, n2
         character(len=128) :: msg
         n1 = nx1*ny1*nz1*nelv
         n2 = nx2*ny2*nz2*nelv
         ! sanity check
         if (irst < 1) then
            write(msg,'(A,I0)') 'Invalid input for irst: ', irst
            call logger%log_error(msg, this_module, 'zget_rst')
         else if (irst > self%nrst) then
            write(msg,'(A,I0)') 'No rst field to retrieve: ', irst
            call logger%log_warning(msg, this_module, 'zget_rst')
         else
            write(msg,'(A,I0)') 'Retrieving rst fields: ', irst
            call logger%log_information(msg, this_module, 'zget_rst')
         end if
         select type (vec_rst)
         type is (nek_zvector)
            call copy(vec_rst%re%vx, self%re_rst(irst)%vx, n1)
            call copy(vec_rst%im%vx, self%im_rst(irst)%vx, n1)
            call copy(vec_rst%re%vy, self%re_rst(irst)%vy, n1)
            call copy(vec_rst%im%vy, self%im_rst(irst)%vy, n1)
            if (if3d) then
               call copy(vec_rst%re%vz, self%re_rst(irst)%vz, n1)
               call copy(vec_rst%im%vz, self%im_rst(irst)%vz, n1)
            end if
            call copy(vec_rst%re%pr, self%re_rst(irst)%pr, n2)
            call copy(vec_rst%im%pr, self%im_rst(irst)%pr, n2)
            if (ifto) then
               call copy(vec_rst%re%theta(:, 1), self%re_rst(irst)%theta(:, 1), n1)
               call copy(vec_rst%im%theta(:, 1), self%im_rst(irst)%theta(:, 1), n1)
            end if
            if (ldimt > 1) then
               do i = 2, ldimt
                  if (ifpsco(i - 1)) then
                     call copy(vec_rst%re%theta(:, i), self%re_rst(irst)%theta(:, i), n1)
                     call copy(vec_rst%im%theta(:, i), self%im_rst(irst)%theta(:, i), n1)
                  end if
               end do
            end if
         class default
            call stop_error("The intent [OUT] argument 'vec_rst' must be of type 'nek_zvector'",
     & this_module, 'get_rst')
         end select
         end procedure

         module procedure zhas_rst_fields
         has_rst_fields = .false.
         if (self%nrst > 0) has_rst_fields = .true.
         end procedure 

         module procedure zclear_rst_fields
         if (self%nrst == 0) then
            call logger%log_debug('No rst fields to clear', this_module, 'zclear_rst_fields')
         end if
         self%nrst = 0
         end procedure
      
      end submodule
