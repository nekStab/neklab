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
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek, setup_nonlinear_solver, setup_linear_solver
         use neklab_linops, only: apply_A
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
            integer :: r = 1
            integer :: startstep = 1
            integer :: printstep = 5
            integer :: orthostep = 10
            integer :: iostep    = 100
            real(dp), allocatable :: sigma(:)
            real(dp), allocatable :: svec(:,:)
            complex(dp), allocatable :: lambda(:)
            complex(dp), allocatable :: eigvec(:,:)
            type(nek_dvector), allocatable :: basis(:)
            type(nek_dvector), allocatable :: modes_re(:)
            type(nek_dvector), allocatable :: modes_im(:)
            type(nek_dvector), allocatable :: Lu(:)
            real(dp), allocatable :: Lr(:,:)
            real(dp), allocatable :: Phi(:,:)
            real(dp) :: OTDrsttime = 0.0_dp
         contains
            private
            procedure, pass(self), public ::  matvec => apply_LNS
            procedure, pass(self), public :: rmatvec => apply_adjLNS
            procedure, pass(self), public :: init    => init_OTD
            procedure, pass(self), public :: generate_Lu
            procedure, pass(self), public :: compute_Lr
            procedure, pass(self), public :: compute_Phi
            procedure, pass(self), public :: project_OTDbasis
            procedure, pass(self), public :: generate_OTD_forcing
         end type nek_otd

      contains

         subroutine apply_LNS(self, vec_in, vec_out)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! Input vector.
            class(abstract_vector_rdp), intent(in)  :: vec_in
            ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            type(nek_dvector) :: vec
            select type(vec_in)
            type is (nek_dvector)
               select type(vec_out)
               type is (nek_dvector)
                  call copy(vec, vec_in)
                  call apply_L(vec%vx, vec%vy, vec%vz, trans=.false.)
                  call copy(vec_out, vec)
               end select
            end select
            return
         end subroutine apply_LNS

         subroutine apply_adjLNS(self, vec_in, vec_out)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! Input vector.
            class(abstract_vector_rdp), intent(in)  :: vec_in
            ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            type(nek_dvector) :: vec
            select type(vec_in)
            type is (nek_dvector)
               select type(vec_out)
               type is (nek_dvector)
                  call copy(vec, vec_in)
                  call apply_L(vec%vx, vec%vy, vec%vz, trans=.true.)
                  call copy(vec_out, vec)
               end select
            end select
            return
         end subroutine apply_adjLNS

         subroutine init_OTD(self, bf, startstep, nIC)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! Baseflow
            type(nek_dvector), intent(in) :: bf
            ! When to start computing OTD modes
            integer, optional, intent(in) :: startstep
            integer, optional, intent(in) :: nIC
            integer :: nICflds
            ! internal
            read(dp), dimension(lv) :: bvx, bvy, bvz
            character(len=128) :: msg
            character(len=128) :: ifile
            logical :: loadIC, exist_IC
            integer :: i

            nICflds = optval(nIC, 0)

            ! set startstep
            self%startstep = optval(startstep, 1)

            ! set number of modes and allocate data
            self%r = npert
            allocate(self%basis(self%r), source=bf);    call zero_basis(self%basis)
            allocate(self%modes_re(self%r), source=bf); call zero_basis(self%modes_re)
            allocate(self%modes_im(self%r), source=bf); call zero_basis(self%modes_im)
            allocate(self%Lu(self%r), source=bf);       call zero_basis(self%Lu)
            allocate(self%Lr(self%r, self%r));  self%Lr = 0.0_dp
            allocate(self%Phi(self%r, self%r)); self%Phi = 0.0_dp

            loadIC = .false.
            if (nICflds < 0) then
               write(msg,*) 'Incorrect number of IC fields to load. nIC=', nICflds
               call logger%log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else if (nICflds > self%r) then
               write(msg,*) 'Inconsistent number of IC fields to load. nIC=', nICflds, ' r=', self%r
               call logger%log_message(msg, module=this_module, procedure='init_OTD')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               loadIC = .true.
            end if

            call rand_basis(self%basis)

            if (loadIC) then
               do i = 1, nICflds
                  write(ifile,'(A,I2.2,".fld")') OTDIC_basename, i
                  inquire (file=ifile,exist=exist_IC)
                  if (exist_IC) then
                     write(msg,*) 'Loading IC file', ifile
                     call logger%log_message(msg, module=this_module, procedure='init_OTD')
                     call load_fld(ifile)
                     call nek2vec(sefl%basis(i), vx, vy, vz, pr, t)
                     self%OTDrsttime = time
                  else
                     write(msg,*) 'Cannot find IC file', ifile
                     call logger%log_message(msg, module=this_module, procedure='init_OTD')
                     if (nid == 0) print *, trim(msg)
                     call nek_end()
                  end if
               end do
            end if

            ! orthonormalize
            call orthonormalize_basis(self%basis)

            return
         end subroutine init_OTD

         subroutine generate_Lu(self, trans)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! transpose ?
            logical, optional, intent(in) :: trans
            logical :: transpose
            ! internal
            integer :: i
            transpose = optval(trans, .false.)
            do i = 1, self%r
               if (tranpose) then
                  call self%apply_rmatvec(self%basis(i),self%Lu(i))
               else
                  call self%apply_matvec(self%basis(i),self%Lu(i))
               end if
            end do
            return
         end subroutine generate_Lu

         subroutine compute_Lr(self)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            call innerprod(self%Lr, self%basis, self%Lu)
            return
         end subroutine compute_Lr

         subroutine compute_Phi(self,)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! internal
            integer i, j
            call zero_basis(self%Phi)
            do i = 1, self%r
               do j = i+1, self%rand
                  self%Phi(i,j) =  self%Lr(i,j)
                  self%Phi(j,1) = -self%Lr(i,j)
               end do
            end do
            return
         end subroutine compute_Phi

         subroutine spectral_analysis(self)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! internals
            real(dp), dimension(self%r)           :: s
            complex(dp), dimension(self%r)        :: l
            complex(dp), dimension(self%r,self%r) :: v
            integer :: i
            integer :: idx(self%r)
            character(len=128) :: msg

            ! compute eigenvalues of the symmetrized operator
            call eig(0.5_dp*(self%Lr + transpose(self%Lr)), l, right=self%svec)
            self%sigma = real(l)
            call sort(self%sigma, reverse=.true.)
            if (mod(istep,self%printstep) == 0) then
               write(msg,'(I7,1X,F15.8,*(1X,E15.8))') istep, time, self%sigma
               call logger%log_message(msg, module=this_module, procedure='OTD Ls')
            end if
            ! outpost projected modes?

            ! compute eigenvalues of Lr
            call eig(self%Lr, l, right=v)
            s = real(l)
            allocate(idx, mold=s)
            call sort_index(s, idx, reverse=.true.)
            do i = 1, self%r
               self%lambda(i) = l(idx(i))
               self%eigvec(:,i) = v(:,idx(i))
            end do
            if (mod(istep,self%printstep) == 0) then
               write(msg,'(I7,1X,F15.8,*(1X,E15.8))') istep, time, real(self%lambda)
               call logger%log_message(msg, module=this_module, procedure='OTD Lr Re')
               write(msg,'(I7,1X,F15.8,*(1X,E15.8))') istep, time, aimag(self%lambda)
               call logger%log_message(msg, module=this_module, procedure='OTD Lr Im')
            end if           
            
            return
         subroutine spectral_analysis

         subroutine outpost_OTDmodes(self)
            !! We assume the OTD basis has already been copied to v[xyz]p
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! internal
            real(dp), dimension(lv,self%r) :: vxm_re, vxm_im
            real(dp), dimension(lv,self%r) :: vym_re, vym_im
            real(dp), dimension(lv,self%r) :: vzm_re, vzm_im
            integer :: i
            character*3 :: prefix

            if (mod(istep, self%printstep) /= 0) then
               self%spectral_analysis()
            end if

            ! project modes
            call mxm(vxp,lv,real(self%lambda),self%r,vxm_re,self%r)
            !call mxm(vxp,lv,aimag(self%lambda),self%r,vxm_im,self%r)
            call mxm(vyp,lv,real(self%lambda),self%r,vym_re,self%r)
            !call mxm(vyp,lv,aimag(self%lambda),self%r,vym_im,self%r)
            if (if3d) then
              call mxm(vzp,lv,real(self%lambda),self%r,vzm_re,self%r)
              !call mxm(vzp,lv,aimag(self%lambda),self%r,vzm_im,self%r)
            end if

            do i = 1, self%r
               if (self%r > 9) then
                  write(prefix,'(A,I2.2)') 'o', i
               else
                  write(prefix,'(A,I1)') 'ot', i
               end if
               call outpost(vxm_re(1,i),vym_re(1,i),vzm_re(1,i),prp(1,i),tp,prefix)
            end do
            return
         end subroutine outpost_OTDmodes

         subroutine generate_OTD_forcing(self, OTDfx, OTDfy, OTDfz)
            ! Linear Operator.
            class(nek_otd), intent(in) :: self
            ! Forcing components
            real(dp), dimension(lv, self%r), intent(out) :: OTDfx
            real(dp), dimension(lv, self%r), intent(out) :: OTDfy
            real(dp), dimension(lv, self%r), intent(out) :: OTDfz

            call mxm(VXP, lv, self%Lr - self%Phi, self%r, OTDfx, self%r)
            call mxm(VYP, lv, self%Lr - self%Phi, self%r, OTDfy, self%r)
            call mxm(VZP, lv, self%Lr - self%Phi, self%r, OTDfz, self%r)
            return
         end subroutine generate_OTD_forcing   
      
      end module neklab_otd
