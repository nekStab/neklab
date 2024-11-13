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
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek
         use neklab_nek_setup, only: setup_linear_solver
         use neklab_linops, only: apply_L
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
         character(len=*), parameter, private :: OTDIC_basename = 'OTDIC_'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
         type, extends(abstract_linop_rdp), public :: nek_otd
      ! to be set at definition
            type(nek_dvector), allocatable :: baseflow
            integer :: r = lpert
      ! defined in init
            integer :: startstep = 1
            integer :: printstep = 5
            integer :: orthostep = 10
            integer :: iostep = 100
            integer :: iorststep = 100
            real(dp) :: OTDrsttime = 0.0_dp
            logical :: trans = .false.
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
         end type nek_otd
      
      contains
      
         subroutine apply_LNS(self, vec_in, vec_out)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_dvector) :: veci, veco
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  call copy(veci, vec_in)
                  call apply_L(veco%vx, veco%vy, veco%vz,
     $   veci%vx, veci%vy, veci%vz, veci%pr, trans = .false.)
                  call copy(vec_out, veco)
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
      ! internal
            type(nek_dvector) :: veci, veco
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  call copy(veci, vec_in)
                  call apply_L(veco%vx, veco%vy, veco%vz,
     $   veci%vx, veci%vy, veci%vz, veci%pr, trans = .true.)
                  call copy(vec_out, veco)
               end select
            end select
            return
         end subroutine apply_adjLNS
      
         subroutine init_OTD(self, startstep, nIC)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! When to start computing OTD modes
            integer, optional, intent(in) :: startstep
            integer, optional, intent(in) :: nIC
            integer :: nICflds
      ! internal
            character(len=128) :: msg
            character(len=128) :: ifile
            logical :: loadIC, exist_IC
            integer :: i, r
      
            r = self%r
      
      ! set fields to read
            nICflds = optval(nIC, 0)
      
      ! set startstep
            self%startstep = optval(startstep, 1)
      
      ! set number of modes and allocate data
            call zero_basis(self%basis)
      
            loadIC = .false.
            if (nICflds < 0) then
               write (msg, *) 'Incorrect number of IC fields to load. nIC=', nICflds
               call logger%log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else if (nICflds > r) then
               write (msg, *) 'Inconsistent number of IC fields to load. nIC=', nICflds, ' r=', self%r
               call logger%log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               loadIC = .true.
            end if
      
            call rand_basis(self%basis)
      
            if (loadIC) then
            do i = 1, nICflds
               write (ifile, '(A,I2.2,".fld")') OTDIC_basename, i
               inquire (file=ifile, exist=exist_IC)
               if (exist_IC) then
                  write (msg, *) 'Loading IC file', ifile
                  call logger%log_message(msg, module=this_module, procedure='init_OTD')
                  call load_fld(ifile)
                  call nek2vec(self%basis(i), vx, vy, vz, pr, t)
                  self%OTDrsttime = time
               else
                  write (msg, *) 'Cannot find IC file', ifile
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
      
      ! setup nek5000
            call setup_linear_solver(recompute_dt=.true., cfl_limit=0.4_dp, solve_baseflow=.true.)
      
            return
         end subroutine init_OTD
      
         subroutine spectral_analysis(self, Lr, sigma, svec, lambda, eigvec)
      ! Linear Operator.
            class(nek_otd), intent(in) :: self
            real(dp), dimension(self%r, self%r), intent(out) :: Lr
            real(dp), dimension(self%r), intent(out) :: sigma
            real(dp), dimension(self%r, self%r), intent(out) :: svec
            complex(dp), dimension(self%r), intent(out) :: lambda
            complex(dp), dimension(self%r, self%r), intent(out) :: eigvec
      ! internals
            real(dp), dimension(self%r) :: s
            real(dp), dimension(self%r, self%r) :: Lsym
            complex(dp), dimension(:), allocatable :: l
            complex(dp), dimension(:, :), allocatable :: v
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
            if (mod(istep, self%printstep) == 0) then
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, sigma
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
            if (mod(istep, self%printstep) == 0) then
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, real(lambda)
               call logger%log_message(msg, module=this_module, procedure='OTD Lr Re')
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, aimag(lambda)
               call logger%log_message(msg, module=this_module, procedure='OTD Lr Im')
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
      ! project modes (imaginary part)
      !call mxm(vxp, lv,aimag(eigvec), r, vxm%im, r)
      !call mxm(vyp, lv,aimag(eigvec), r, vym%im, r)
      !if (if3d) then
      !   call mxm(vzp,lv,aimag(eigvec), r, vzm%im, r)
      !end if
      
            do i = 1, self%r
            if (self%r > 9) then
               write (file_prefix, '(A,I2.2)') 'o', i
            else
               write (file_prefix, '(A,I1)') 'ot', i
            end if
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
      
      end module neklab_otd
