      submodule(neklab_vectors) complex_vectors
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------

         module procedure construct_nek_zvector
            character(len=*), parameter :: this_procedure = 'construct_nek_zvector'
            call out%init() ! re/im initialized without rst storage
            out%re%vx = vx(:, 1); out%im%vx = vx(:, 2)
            out%re%vy = vy(:, 1); out%im%vy = vy(:, 2)
            if (present(vz)) then
               if (if3d) then
                  out%re%vz = vz(:, 1); out%im%vz = vz(:, 2)
               else
                  call stop_error('vz passed but if3d = .false.', this_module, this_procedure)
               end if
            end if
            if (present(pr)) then
               out%re%pr = pr(:, 1); out%im%pr = pr(:, 2)
            end if
            if (present(theta)) then
               if (ifto) then
                  out%re%theta = theta(:, :, 1); out%im%theta = theta(:, :, 2)
               else
                  call stop_error('theta passed but ifto = .false.', this_module, this_procedure)
               end if
            end if
         end procedure
      
      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nek_zzero
            character(len=*), parameter :: this_procedure = 'nek_zzero'
            call assert_zinitialized(self, this_procedure)
            call nek_dzero(self%re)
            call nek_dzero(self%im)
            ! this causes an internal compiler error, switch to manual
            !call self%re%zero()
            !call self%im%zero()
         end procedure
      
         module procedure nek_zrand
            character(len=*), parameter :: this_procedure = 'nek_zrand'
            real(dp) :: alpha
            call assert_zinitialized(self, this_procedure)
            call nek_drand(self%re)
            call nek_drand(self%im)
            if (optval(ifnorm, .false.)) then
               alpha = self%norm()
               call nek_zscal(self, one_cdp/alpha)
            end if
         end procedure
      
         module procedure nek_zscal
            character(len=*), parameter :: this_procedure = 'nek_zscal'
            type(nek_dvector) :: re_old, im_old
            call assert_zinitialized(self, this_procedure)

            ! Save components of complex vector before scaling.
            re_old = self%re
            im_old = self%im

            ! re_new = ar*re_old - ai*im_old
            call nek_dzero(self%re)
            call nek_daxpby(alpha%re,  re_old, 1.0_dp, self%re)
            call nek_daxpby(-alpha%im, im_old, 1.0_dp, self%re)

            ! im_new = ar*im_old + ai*re_old
            call nek_dzero(self%im)
            call nek_daxpby(alpha%re, im_old, 1.0_dp, self%im)
            call nek_daxpby(alpha%im, re_old, 1.0_dp, self%im)
         end procedure
         
         module procedure nek_zaxpby
            character(len=*), parameter :: this_procedure = 'nek_zaxpby'
            type(nek_zvector) :: wrk
            call assert_zinitialized(self, this_procedure)
         
            select type (vec)
            type is (nek_zvector)
            call assert_zinitialized(vec, this_procedure)

               wrk = vec                       ! deep copy (re/im allocatables)
               call nek_zscal(self, beta)
               call nek_zscal(wrk, alpha)
               call nek_daxpby(1.0_dp, wrk%re, 1.0_dp, self%re)
               call nek_daxpby(1.0_dp, wrk%im, 1.0_dp, self%im)

            class default
               call type_error('vec','nek_zvector','IN',this_module,this_procedure)
            end select

         end procedure
      
         module procedure nek_zdot
            character(len=*), parameter :: this_procedure = 'nek_zdot'
            real(kind=dp) :: alpha_r, alpha_i
            call assert_zinitialized(self, this_procedure)

            select type (vec)
            type is (nek_zvector)
               call assert_zinitialized(vec, this_procedure)
               ! <self, vec> = <re,re> + <im,im> + i(<re,im> - <im,re>)
               alpha_r = nek_ddot(self%re, vec%re) + nek_ddot(self%im, vec%im)
               alpha_i = nek_ddot(self%re, vec%im) - nek_ddot(self%im, vec%re)
               alpha = cmplx(alpha_r, alpha_i, kind=dp)

            class default
               call type_error('vec','nek_zvector','IN',this_module,this_procedure)
            end select
         end procedure
      
         module procedure nek_zsize
            ! delegate: real + imaginary linear-algebra dimensions
            n = nek_dsize(self%re) + nek_dsize(self%im)
         end procedure

         ! initialization

         module procedure nek_zinit
            character(len=*), parameter :: this_procedure = 'nek_zinit'
            ! re and im are always initialized together with identical layout,
            ! so their state flags (initialized / has_rst / nrst) stay in lockstep.
            call nek_dinit(self%re, alloc_rst)
            call nek_dinit(self%im, alloc_rst)
         end procedure

         module procedure nek_zis_initialized
            ok = .false.
            if (nek_dis_initialized(self%re) .and. nek_dis_initialized(self%im)) ok = .true.
         end procedure

         module procedure zsave_rst
            character(len=*), parameter :: this_procedure = 'zsave_rst'
            call assert_zrst_storage(self, this_procedure)
            select type (vec_rst)
            type is (nek_zvector)
               call dsave_rst(self%re, vec_rst%re, irst)
               call dsave_rst(self%im, vec_rst%im, irst)
            class default
               call type_error('vec_rst','nek_zvector','IN', this_module, this_procedure)
            end select
         end procedure

         module procedure zget_rst
            character(len=*), parameter :: this_procedure = 'zget_rst'
            call assert_zrst_storage(self, this_procedure)
            select type (vec_rst)
            type is (nek_zvector)
               call dget_rst(self%re, vec_rst%re, irst)
               call dget_rst(self%im, vec_rst%im, irst)
            class default
               call type_error('vec_rst','nek_zvector','OUT', this_module, this_procedure)
            end select
         end procedure

         module procedure zhas_rst_fields
            has = dhas_rst_fields(self%re)
         end procedure

         module procedure zhas_rst_storage
            has = dhas_rst_storage(self%re)
         end procedure

         module procedure zclear_rst_fields
            call dclear_rst_fields(self%re, dealloc_rst)
            call dclear_rst_fields(self%im, dealloc_rst)
         end procedure

         ! --- private helpers ---

         subroutine assert_zinitialized(vec, procedure_name)
            class(nek_zvector), intent(in) :: vec
            character(len=*), intent(in) :: procedure_name
            if (.not. (nek_dis_initialized(vec%re) .and. nek_dis_initialized(vec%im))) then
               call stop_error('nek_zvector used before init', this_module, procedure_name)
            end if
         end subroutine

         subroutine assert_zrst_storage(vec, procedure_name)
            class(nek_zvector), intent(in) :: vec
            character(len=*), intent(in) :: procedure_name
            if (.not. (nek_dis_initialized(vec%re) .and. nek_dis_initialized(vec%im))) then
               call stop_error('nek_zvector used before init', this_module, procedure_name)
            else if (.not. dhas_rst_storage(vec%re)) then
               call stop_error('nek_zvector initialized without rst storage', this_module, procedure_name)
            end if
         end subroutine

      end submodule