      module neklab_otd
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_sorting, only: sort, sort_index
         use stdlib_linalg, only: eig
         use stdlib_optval, only: optval
      ! Default real kind.
         use LightKrylov, only: dp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov_Utils, only: abstract_opts
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek, neklab_forcing
         use neklab_nek_setup, only: setup_linear_solver
         use neklab_linops, only: apply_L
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
         type, extends(abstract_linop_rdp), public :: nek_otd
            integer :: r = lpert
            type(nek_dvector), allocatable :: baseflow
            type(nek_dvector), allocatable :: basis(:)
      ! Forcing components
            real(dp), dimension(lv, lpert) :: OTDfx = 0.0_dp
            real(dp), dimension(lv, lpert) :: OTDfy = 0.0_dp
            real(dp), dimension(lv, lpert) :: OTDfz = 0.0_dp
         contains
            private
            procedure, pass(self), public :: init => init_OTD
            procedure, pass(self), public :: matvec => apply_LNS
            procedure, pass(self), public :: rmatvec => apply_adjLNS
            procedure, pass(self), public :: spectral_analysis
            procedure, pass(self), public :: outpost_OTDmodes
            procedure, pass(self), public :: generate_forcing
            procedure, pass(self), public :: set_forcing
         end type nek_otd
      
         type, extends(abstract_opts), public :: otd_opts
            integer :: startstep = 1
      !! Timestep at which to start OTD computation
            integer :: printstep = 5
      !! Output frequency for spectral analysis
            integer :: orthostep = 10
      !! Reorthogonalization frequency
            integer :: iostep = 100
      !! Output frequency for the projected basis vectors
            integer :: iorststep = 100
      !! Output frequency for the basis (for restarts)
            integer :: n_usrIC = 0
      !! Number of user-defined initial conditions for the perturbations
            logical :: trans = .false.
      !! Direct of adjoint?
            logical :: solve_baseflow = .true.
      !! Solve nonlinear equations alongside the linear equations?
            logical :: if_output_initial_conditions = .false.
      !! Outpost initial state (OTDinit)
            character(len=128) :: OTDIC_basename = 'OTDIC_'
      !! Base filename for initial conditions
         end type
      
      ! Module level pointer to the current instance of nek_otd
         type(nek_otd), pointer, public :: otd_instance => null()
      
      contains
      
         subroutine apply_LNS(self, vec_in, vec_out)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  call apply_L(vec_out%vx, vec_out%vy, vec_out%vz,
     $   vec_in%vx, vec_in%vy, vec_in%vz, vec_in%pr, trans = .false.)
               end select
            end select
            return
         end subroutine apply_LNS
      
         subroutine apply_adjLNS(self, vec_in, vec_out)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  call apply_L(vec_out%vx, vec_out%vy, vec_out%vz,
     $   vec_in%vx, vec_in%vy, vec_in%vz, vec_in%pr, trans = .true.)
               end select
            end select
            return
         end subroutine apply_adjLNS
      
         subroutine init_OTD(self, opts)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            type(otd_opts), intent(in) :: opts
      ! internal
            character(len=128) :: msg
            character(len=128) :: ifile
            logical :: loadIC, exist_IC
            integer :: i, r
      
      ! number of OTD modes
            r = self%r
      
      ! Associate the current instance and the neklab forcing
            call set_neklab_forcing(self)
      
      ! allocate variables
            allocate (self%basis(r), source=self%baseflow)
            call zero_basis(self%basis)
      
            loadIC = .false.
            if (opts%n_usrIC < 0) then
               write (msg, *) 'Incorrect number of IC fields to load. nIC=', opts%n_usrIC
               call logger%log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else if (opts%n_usrIC > r) then
               write (msg, *) 'Inconsistent number of IC fields to load. nIC=', opts%n_usrIC, ' r=', self%r
               call logger%log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               loadIC = .true.
            end if
      
            call rand_basis(self%basis)
      
            if (loadIC) then
            do i = 1, opts%n_usrIC
               write (ifile, '(A,I2.2,".fld")') trim(opts%OTDIC_basename), i
               inquire (file=ifile, exist=exist_IC)
               if (exist_IC) then
                  write (msg, *) 'Loading IC file: ', trim(ifile)
                  call logger%log_message(msg, module=this_module, procedure='init_OTD')
                  call load_fld(ifile)
                  call nek2vec(self%basis(i), vx, vy, vz, pr, t)
               else
                  write (msg, *) 'Cannot find IC file: ', trim(ifile)
                  call logger%log_message(msg, module=this_module, procedure='init_OTD')
                  if (nid == 0) print *, trim(msg)
                  call nek_end()
               end if
            end do
            end if
      
      ! orthonormalize
            call orthonormalize_basis(self%basis)
      
      ! force baseflow
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! set perturbations
            do i = 1, r
               call vec2nek(vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i), self%basis(i))
            end do
      
      ! outpost initial conditions (diagnostics)
            if (opts%if_output_initial_conditions) then
               call outpost(vx, vy, vz, pr, t, 'bfi')
               do i = 1, r
                  call outpost(vxp(1, i), vyp(1, i), vzp(1, i), prp(1, i), tp(1, :, i), 'pri')
               end do
            end if
      
      ! setup nek5000
            call setup_linear_solver(recompute_dt=.true., cfl_limit=0.4_dp, solve_baseflow=opts%solve_baseflow)
      
            return
         end subroutine init_OTD
      
         subroutine spectral_analysis(self, Lr, sigma, svec, lambda, eigvec, ifprint)
      ! Linear Operator.
            class(nek_otd), intent(in) :: self
            real(dp), dimension(self%r, self%r), intent(out) :: Lr
            real(dp), dimension(self%r), intent(out) :: sigma
            real(dp), dimension(self%r, self%r), intent(out) :: svec
            complex(dp), dimension(self%r), intent(out) :: lambda
            complex(dp), dimension(self%r, self%r), intent(out) :: eigvec
            logical, intent(in) :: ifprint
      ! internals
            real(dp), dimension(self%r) :: s
            real(dp), dimension(self%r, self%r) :: Lsym
            complex(dp), dimension(self%r) :: l
            complex(dp), dimension(self%r, self%r) :: v
            integer :: i, r
            integer :: idx(self%r)
            character(len=128) :: msg
      
            r = self%r
      ! compute eigenvalues of the symmetrized operator
            Lsym = 0.5_dp*(Lr + transpose(Lr))
            call eig(Lsym, l, right=v)
            s = real(l)
            call sort_index(s, idx, reverse=.true.)
            do i = 1, r
               sigma(i) = l(idx(i))
               svec(:, i) = v(:, idx(i))
            end do
            call sort(sigma, reverse=.true.)
            if (ifprint) then
               write (msg, '(I10,1X,F15.8,*(1X,E15.8))') istep, time, sigma
               call logger%log_message(msg, module=this_module, procedure='OTD Ls')
            end if
      ! outpost projected modes?
      
      ! compute eigenvalues of Lr
            call eig(Lr, l, right=v)
            s = real(l)
            call sort_index(s, idx, reverse=.true.)
            do i = 1, r
               lambda(i) = l(idx(i))
               eigvec(:, i) = v(:, idx(i))
            end do
            if (ifprint) then
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, real(lambda)
               call logger%log_message(msg, module=this_module, procedure='OTD Lr%Re')
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, aimag(lambda)
               call logger%log_message(msg, module=this_module, procedure='OTD Lr%Im')
            end if
            return
         end subroutine spectral_analysis
      
         subroutine outpost_OTDmodes(self, eigvec)
      !! We assume the OTD basis has already been copied to v[xyz]p
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            complex(dp), dimension(self%r, self%r), intent(in) :: eigvec
      ! internal
            complex(dp), dimension(lv, self%r) :: vxm
            complex(dp), dimension(lv, self%r) :: vym
            complex(dp), dimension(lv, self%r) :: vzm
            integer :: i, r
            character(len=3) :: file_prefix
      
            r = self%r
      
      ! project modes (real part)
            call mxm(vxp, lv, real(eigvec), r, vxm%re, r)
            call mxm(vyp, lv, real(eigvec), r, vym%re, r)
            if (if3d) then
               call mxm(vzp, lv, real(eigvec), r, vzm%re, r)
            end if
      !! project modes (imaginary part)
      !      call mxm(vxp, lv,aimag(eigvec), r, vxm%im, r)
      !      call mxm(vyp, lv,aimag(eigvec), r, vym%im, r)
      !      if (if3d) then
      !         call mxm(vzp,lv,aimag(eigvec), r, vzm%im, r)
      !      end if
            do i = 1, self%r
               write (file_prefix, '(A,I2.2)') 'm', i
               call outpost(vxm(1, i)%re, vym(1, i)%re, vzm(1, i)%re,
     $   prp(1, i), tp(1, :, i), file_prefix)
            end do
            return
         end subroutine outpost_OTDmodes
      
         subroutine generate_forcing(self, Lr, Phi)
      !! We assume the OTD basis has already been copied to v[xyz]p
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            real(dp), dimension(self%r, self%r), intent(in) :: Lr
            real(dp), dimension(self%r, self%r), intent(in) :: Phi
      ! internal
            integer :: r
            r = self%r
            call mxm(vxp, lv, Lr - Phi, r, self%OTDfx, r)
            call mxm(vyp, lv, Lr - Phi, r, self%OTDfy, r)
            call mxm(vzp, lv, Lr - Phi, r, self%OTDfz, r)
            return
         end subroutine generate_forcing
      
         subroutine set_forcing(self, ffx, ffy, ffz, ix, iy, iz, ieg, ipert)
      ! Linear Operator.
            class(nek_otd), intent(in) :: self
            real(dp), intent(inout) :: ffx
            real(dp), intent(inout) :: ffy
            real(dp), intent(inout) :: ffz
            integer, intent(in) :: ix
            integer, intent(in) :: iy
            integer, intent(in) :: iz
            integer, intent(in) :: ieg
            integer, intent(in) :: ipert
      ! internal
            integer :: e, ijke
            e = gllel(ieg)
            ijke = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            if (ipert /= 0) then
               ffx = ffx + self%OTDfx(ijke, ipert)
               ffy = ffy + self%OTDfy(ijke, ipert)
               ffz = ffz + self%OTDfz(ijke, ipert)
            end if
         end subroutine set_forcing
      
         subroutine set_neklab_forcing(instance)
            type(nek_otd), target, intent(inout) :: instance
      ! set pointer for current instance
            otd_instance => instance
            if (associated(otd_instance)) then
               call logger%log_debug('OTD instance associated!', module=this_module, procedure='set_neklab_forcing')
            else
               call logger%log_message('OTD instance not associated!', module=this_module, procedure='set_neklab_forcing')
               call nek_end()
            end if
      ! set neklab forcing pointer to otd_forcing function
            neklab_forcing => wrapper_neklab_forcing
            if (associated(neklab_forcing)) then
               call logger%log_debug('neklab_forcing set to OTD forcing!', module=this_module, procedure='set_neklab_forcing')
            else
               call logger%log_message('neklab_forcing not associated!', module=this_module, procedure='set_neklab_forcing')
               call nek_end()
            end if
            return
         end subroutine set_neklab_forcing
      
         subroutine wrapper_neklab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg, ipert)
            real(dp), intent(inout) :: ffx
            real(dp), intent(inout) :: ffy
            real(dp), intent(inout) :: ffz
            integer, intent(in) :: ix
            integer, intent(in) :: iy
            integer, intent(in) :: iz
            integer, intent(in) :: ieg
            integer, intent(in) :: ipert
            if (associated(otd_instance)) then
               call otd_instance%set_forcing(ffx, ffy, ffz, ix, iy, iz, ieg, ipert)
            else
               call logger%log_message('otd_instance not associated.', module=this_module, procedure='wrapper_neklab_forcing')
               call nek_end()
            end if
            return
         end subroutine wrapper_neklab_forcing
      
      end module neklab_otd
