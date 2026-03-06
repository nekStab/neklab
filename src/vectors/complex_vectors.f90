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
         integer :: i
         ! this causes an internal compiler error, switch to manual
         !call self%re%zero()
         !call self%im%zero()
         self%re%vx = 0.0_dp
         self%re%vy = 0.0_dp
         self%re%vz = 0.0_dp
         self%re%pr = 0.0_dp
         self%re%theta = 0.0_dp
         self%re%vyrst = 0.0_dp
         self%re%vxrst = 0.0_dp
         self%re%vzrst = 0.0_dp
         self%re%prrst = 0.0_dp
         self%re%thetarst = 0.0_dp
         self%im%vx = 0.0_dp
         self%im%vy = 0.0_dp
         self%im%vz = 0.0_dp
         self%im%pr = 0.0_dp
         self%im%theta = 0.0_dp
         self%im%vyrst = 0.0_dp
         self%im%vxrst = 0.0_dp
         self%im%vzrst = 0.0_dp
         self%im%prrst = 0.0_dp
         self%im%thetarst = 0.0_dp
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
         integer :: lvn, lpn, m
         lvn = nx1*ny1*nz1*nelv
         lpn = nx2*ny2*nz2*nelv
         n = 2*lvn + lpn
         if (if3d) n = n + lvn
         if (ifto) n = n + lvn
         if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) n = n + lvn
            end do
         end if
         end procedure

         module procedure zsave_rst
         integer :: m, lvn, lpn, torder
         character(len=*), parameter :: this_procedure = 'zsave_rst'
         character(len=128) :: msg

         lvn = nx1*ny1*nz1*nelv
         lpn = nx2*ny2*nz2*nelv
         torder = abs(param(27)) ! integration order in time

         ! sanity checks
         if (irst == torder) then
            write(msg,'(2(A,I0),A)') 'Cannot save rst fields ', torder, ' for a simulation of temporal order ', torder, '.'
            call log_error(msg, this_module, this_procedure)
         else
            write(msg,'(A,I0)') 'Saving rst fields: ', irst
            call log_debug(msg, this_module, this_procedure)
         end if

         select type (vec_rst)
         type is (nek_zvector)
            ! associate?
            call copy(self%re%vxrst(:, irst), vec_rst%re%vx, lvn)
            call copy(self%im%vxrst(:, irst), vec_rst%im%vx, lvn)
            call copy(self%re%vyrst(:, irst), vec_rst%re%vy, lvn)
            call copy(self%im%vyrst(:, irst), vec_rst%im%vy, lvn)
            if (if3d) then
               call copy(self%re%vzrst(:, irst), vec_rst%re%vz, lvn)
               call copy(self%im%vzrst(:, irst), vec_rst%im%vz, lvn)
            end if
            call copy(self%re%prrst(:, irst), vec_rst%re%pr, lpn)
            call copy(self%im%prrst(:, irst), vec_rst%im%pr, lpn)

            if (ifto) then
               call copy(self%re%thetarst(:, irst, 1), vec_rst%re%theta(:, 1), lvn)
               call copy(self%im%thetarst(:, irst, 1), vec_rst%im%theta(:, 1), lvn)
            end if
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call copy(self%re%thetarst(:, irst, m), vec_rst%re%theta(:, m), lvn)
                     call copy(self%im%thetarst(:, irst, m), vec_rst%im%theta(:, m), lvn)
                  end if
               end do
            end if

         class default
            call type_error('vec_rst','nek_zvector','IN', this_module, this_procedure)
         end select
         end procedure
   
         module procedure zget_rst
         integer :: m, lvn, lpn
         character(len=*), parameter :: this_procedure = 'zget_rst'
         character(len=128) :: msg

         lvn = nx1*ny1*nz1*nelv
         lpn = nx2*ny2*nz2*nelv

         ! sanity checks
         if (irst < 1) then
            write(msg,'(A,I0)') 'Invalid input for irst: ', irst
            call log_error(msg, this_module, this_procedure)
         else if (irst > self%nrst) then
            write(msg,'(A,I0)') 'No rst field to retrieve: ', irst
            call log_warning(msg, this_module, this_procedure)
         else
            write(msg,'(A,I0)') 'Retrieving rst fields: ', irst
            call log_information(msg, this_module, this_procedure)
         end if

         select type (vec_rst)
         type is (nek_zvector)

            call copy(vec_rst%re%vx, self%re%vxrst(:, irst), lvn)
            call copy(vec_rst%im%vx, self%im%vxrst(:, irst), lvn)
            call copy(vec_rst%re%vy, self%re%vyrst(:, irst), lvn)
            call copy(vec_rst%im%vy, self%im%vyrst(:, irst), lvn)
            if (if3d) then
               call copy(vec_rst%re%vz, self%re%vzrst(:, irst), lvn)
               call copy(vec_rst%im%vz, self%im%vzrst(:, irst), lvn)
            end if
            call copy(vec_rst%re%pr, self%re%prrst(:, irst), lpn)
            call copy(vec_rst%im%pr, self%im%prrst(:, irst), lpn)
            if (ifto) then
               call copy(vec_rst%re%theta(:, 1), self%re%thetarst(:, irst, 1), lvn)
               call copy(vec_rst%im%theta(:, 1), self%im%thetarst(:, irst, 1), lvn)
            end if
            if (ldimt > 1) then
               do m = 2, ldimt
                  if (ifpsco(m - 1)) then
                     call copy(vec_rst%re%theta(:, m), self%re%thetarst(:, irst, m), lvn)
                     call copy(vec_rst%im%theta(:, m), self%im%thetarst(:, irst, m), lvn)
                  end if
               end do
            end if

         class default
            call type_error('vec_rst','nek_zvector','OUT',this_module, this_procedure)
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
