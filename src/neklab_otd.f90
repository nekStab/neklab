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
         use neklab_nek_forcing, only: neklab_forcing, set_neklab_forcing
         use neklab_nek_setup, only: nek_stop_error, nek_log_debug, nek_log_message, nek_log_information
         use neklab_utils, only: nek2vec, vec2nek
         use neklab_nek_setup, only: setup_linear_solver
         use neklab_linops, only: apply_L
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
         character(len=*), parameter, private :: logfile_Ls = 'Ls.dat'
         character(len=*), parameter, private :: logfile_Lr = 'Lr.dat'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
         type, extends(abstract_linop_rdp), public :: nek_otd
            integer :: r = lpert
            type(nek_dvector), allocatable :: baseflow
            type(nek_dvector), allocatable :: basis(:)
         contains
            private
            procedure, pass(self), public :: init => init_OTD
            procedure, pass(self), public :: matvec => apply_LNS
            procedure, pass(self), public :: rmatvec => apply_adjLNS
            procedure, pass(self), public :: spectral_analysis
            procedure, pass(self), public :: outpost_OTDmodes
            procedure, pass(self), public :: generate_forcing
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
     & vec_in%vx, vec_in%vy, vec_in%vz, vec_in%pr, trans = .false.)
                  class default
                  call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'apply_LNS')
               end select
            class default
               call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'apply_LNS')
            end select
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
     & vec_in%vx, vec_in%vy, vec_in%vz, vec_in%pr, trans = .true.)
               class default
                  call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'apply_adjLNS')
               end select
            class default
               call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'apply_adjLNS')
            end select
         end subroutine apply_adjLNS
      
         subroutine init_OTD(self, opts)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            type(otd_opts), intent(in) :: opts
      ! internal
            character(len=128) :: msg
            character(len=128) :: ifile
            logical :: loadIC, exist_IC
            integer :: i, j, r
      
      ! number of OTD modes
            r = self%r
      
      ! allocate variables
            allocate (self%basis(r), source=self%baseflow)
            call zero_basis(self%basis)
      
            loadIC = .false.
            if (opts%n_usrIC < 0) then
               write (msg, *) 'Incorrect number of IC fields to load. nIC=', opts%n_usrIC
               call nek_log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else if (opts%n_usrIC > r) then
               write (msg, *) 'Inconsistent number of IC fields to load. nIC=', opts%n_usrIC, ' r=', self%r
               call nek_log_message(msg, module=this_module, procedure='init_OTD')
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
                  call nek_log_message(msg, module=this_module, procedure='init_OTD')
                  call load_fld(ifile)
                  call nek2vec(self%basis(i), vx, vy, vz, pr, t)
               else
                  write (msg, *) 'Cannot find IC file: ', trim(ifile)
                  call nek_log_message(msg, module=this_module, procedure='init_OTD')
                  if (nid == 0) print *, trim(msg)
                  call nek_end()
               end if
            end do
            end if
      
      ! Prepare logfiles
            open (1234, file=logfile_Ls, status='replace', action='write'); close (1234)
            open (1234, file=logfile_Lr, status='replace', action='write'); close (1234)
      
      ! orthonormalize
            write (msg, '(A,*(1X,E10.3))') 'IC: norm.  err pre: ', (self%basis(i)%dot(self%basis(i)) - 1.0_dp, i=1, r)
            call nek_log_information(msg, module=this_module, procedure='OTD init')
            write (msg, '(A,*(1X,E10.3))') 'IC: ortho. err pre: ', ((self%basis(i)%dot(self%basis(j)), j=i + 1, r), i=1, r)
            call nek_log_information(msg, module=this_module, procedure='OTD init')
      
            call orthonormalize_basis(self%basis)
      
            write (msg, '(A,*(1X,E10.3))') 'IC: norm.  err post:', (self%basis(i)%dot(self%basis(i)) - 1.0_dp, i=1, r)
            call nek_log_debug(msg, module=this_module, procedure='OTD init')
            write (msg, '(A,*(1X,E10.3))') 'IC: ortho. err post:', ((self%basis(i)%dot(self%basis(j)), j=i + 1, r), i=1, r)
            call nek_log_debug(msg, module=this_module, procedure='OTD init')
      
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
            character(len=128) :: msg, fmt_Lr
      
            r = self%r
            write (fmt_Lr, '("(I8,1X,F15.8,A,",I0,"(1X,E15.8),A,",I0,"(1X,E15.8))")') r, r
      
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
               call nek_log_message(msg, module=this_module, procedure='OTD Ls')
      ! stamp logfile
               open (1234, file=logfile_Ls, status='old', action='write', position='append')
               write (1234, '(I8,1X,F15.8,A,*(1X,E15.8))') istep, time, ' Ls ', sigma
               close (1234)
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
               call nek_log_message(msg, module=this_module, procedure='OTD Lr%Re')
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, aimag(lambda)
               call nek_log_message(msg, module=this_module, procedure='OTD Lr%Im')
      ! stamp logfile
               open (1234, file=logfile_Lr, status='old', action='write', position='append')
               write (1234, fmt_Lr) istep, time, ' Lr%Re ', real(lambda), ' Lr%Im ', aimag(lambda)
               close (1234)
            end if
         end subroutine spectral_analysis
      
         subroutine outpost_OTDmodes(self, eigvec)
      !! We assume the OTD basis has already been copied to v[xyz]p
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            complex(dp), dimension(self%r, self%r), intent(in) :: eigvec
      ! internal
            real(dp), dimension(lv, self%r) :: vxr
            real(dp), dimension(lv, self%r) :: vyr
            real(dp), dimension(lv, self%r) :: vzr
      !real(dp), dimension(lv, self%r) :: vxi
      !real(dp), dimension(lv, self%r) :: vyi
      !real(dp), dimension(lv, self%r) :: vzi
      
            integer :: i, r
            character(len=3) :: file_prefix
      
            r = self%r
      ! project modes (real part)
            call mxm(vxp, lv, real(eigvec), r, vxr, r)
            call mxm(vyp, lv, real(eigvec), r, vyr, r)
            if (if3d) then
               call mxm(vzp, lv, real(eigvec), r, vzr, r)
            end if
      !! project modes (imaginary part)
      !      call mxm(vxp, lv,aimag(eigvec), r, vxi, r)
      !      call mxm(vyp, lv,aimag(eigvec), r, vyi, r)
      !      if (if3d) then
      !         call mxm(vzp,lv,aimag(eigvec), r, vzi, r)
      !      end if
            do i = 1, self%r
               write (file_prefix, '(A,I2.2)') 'm', i
               call outpost(vxr(1, i), vyr(1, i), vzr(1, i), prp(1, i), tp(1, :, i), file_prefix)
            end do
         end subroutine outpost_OTDmodes
      
         subroutine generate_forcing(self, Lr, Phi)
      !! We assume the OTD basis has already been copied to v[xyz]p
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            real(dp), dimension(self%r, self%r), intent(in) :: Lr
            real(dp), dimension(self%r, self%r), intent(in) :: Phi
      ! internal
            integer :: r, i
            real(dp), dimension(lv, self%r) :: OTDfx
            real(dp), dimension(lv, self%r) :: OTDfy
            real(dp), dimension(lv, self%r) :: OTDfz
            r = self%r
            call mxm(vxp, lv, Lr - Phi, r, OTDfx, r)
            call mxm(vyp, lv, Lr - Phi, r, OTDfy, r)
            call mxm(vzp, lv, Lr - Phi, r, OTDfz, r)
            do i = 1, r
               call set_neklab_forcing(OTDfx(:, i), OTDfy(:, i), OTDfy(:, i), ipert=i)
            end do
         end subroutine generate_forcing
      
      end module neklab_otd
