      module neklab_vectors
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use LightKrylov, only: abstract_vector_rdp, abstract_vector_cdp
         use LightKrylov_Logger
      
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 'neklab_vectors'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         complex(kind=dp), parameter :: zero_cdp = cmplx(0.0_dp, 0.0_dp, kind=dp)
         complex(kind=dp), parameter :: one_cdp  = cmplx(1.0_dp, 0.0_dp, kind=dp)
         complex(kind=dp), parameter :: im_cdp   = cmplx(0.0_dp, 1.0_dp, kind=dp)
      
         public :: mth_rand
      
      !----------------------------------------
      !-----     NEK REAL VECTOR TYPE     -----
      !----------------------------------------
      
      ! --> Type.
         type, extends(abstract_vector_rdp), public :: nek_dvector
            ! public data
            real(kind=dp), allocatable :: vx(:), vy(:), vz(:)
            real(kind=dp), allocatable :: pr(:)
            real(kind=dp), allocatable :: theta(:,:)
            ! private state and restart storage
            real(kind=dp), allocatable, private :: vxrst(:,:), vyrst(:,:), vzrst(:,:)
            real(kind=dp), allocatable, private :: prrst(:,:)
            real(kind=dp), allocatable, private :: thetarst(:,:,:)
            integer, private :: nrst = 0
            logical, private :: has_rst = .false.
            logical, private :: initialized = .false.
         contains
            private
            ! required basic type-bound procedures
            procedure, pass(self), public :: zero     => nek_dzero
            procedure, pass(self), public :: rand     => nek_drand
            procedure, pass(self), public :: scal     => nek_dscal
            procedure, pass(self), public :: axpby    => nek_daxpby
            procedure, pass(self), public :: dot      => nek_ddot
            procedure, pass(self), public :: get_size => nek_dsize
            ! initialization
            procedure, pass(self), public :: init           => nek_dinit
            procedure, pass(self), public :: is_initialized => nek_dis_initialized
            ! additional type-bound procedures
            procedure, pass(self), public :: save_rst         => dsave_rst
            procedure, pass(self), public :: get_rst          => dget_rst
            procedure, pass(self), public :: has_rst_fields   => dhas_rst_fields
            procedure, pass(self), public :: has_rst_storage  => dhas_rst_storage
            procedure, pass(self), public :: clear_rst_fields => dclear_rst_fields
         end type nek_dvector
      
      ! --> Constructor.
         interface nek_dvector
            module function construct_nek_dvector(vx, vy, vz, pr, theta) result(out)
               real(kind=dp), dimension(lv), intent(in) :: vx, vy
               real(kind=dp), dimension(lv), optional, intent(in) :: vz
               real(kind=dp), dimension(lp), optional, intent(in) :: pr
               real(kind=dp), dimension(lv, ldimt), optional, intent(in) :: theta
               type(nek_dvector) :: out
            end function
         end interface
      
      ! --> Type-bound procedures.
         interface

            ! required by LightKrylov abstract_vector_rdp
            module subroutine nek_dzero(self)
               class(nek_dvector), intent(inout) :: self
            end subroutine
      
            module subroutine nek_drand(self, ifnorm)
               class(nek_dvector), intent(inout) :: self
               logical, optional, intent(in) :: ifnorm
            end subroutine
      
            module subroutine nek_dscal(self, alpha)
               class(nek_dvector), intent(inout) :: self
               real(kind=dp), intent(in) :: alpha
            end subroutine
      
            module subroutine nek_daxpby(alpha, vec, beta, self)
               class(nek_dvector), intent(inout) :: self
               real(kind=dp), intent(in) :: alpha
               class(abstract_vector_rdp), intent(in) :: vec
               real(kind=dp), intent(in) :: beta
            end subroutine
      
            real(kind=dp) module function nek_ddot(self, vec) result(alpha)
               class(nek_dvector), intent(in) :: self
               class(abstract_vector_rdp), intent(in) :: vec
            end function
      
            integer pure module function nek_dsize(self) result(n)
               class(nek_dvector), intent(in) :: self
            end function

            ! initialization
            module subroutine nek_dinit(self, alloc_rst)
               class(nek_dvector), intent(inout) :: self
               logical, optional, intent(in) :: alloc_rst ! default is .false. (no restart storage)
            end subroutine

            logical pure module function nek_dis_initialized(self) result(ok)
               class(nek_dvector), intent(in) :: self
            end function

            ! restart storage and organization
            module subroutine dsave_rst(self, vec_rst, irst)
               class(nek_dvector), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_rst
               integer, intent(in) :: irst
            end subroutine
            
            module subroutine dget_rst(self, vec_rst, irst)
               class(nek_dvector), intent(in) :: self
               class(abstract_vector_rdp), intent(inout) :: vec_rst
               integer, intent(in) :: irst
            end subroutine
            
            logical pure module function dhas_rst_fields(self) result(has) ! nrst > 0
               class(nek_dvector), intent(in) :: self
            end function

            logical pure module function dhas_rst_storage(self) result(has)
               class(nek_dvector), intent(in) :: self
            end function
            
            module subroutine dclear_rst_fields(self, dealloc_rst)
               class(nek_dvector), intent(inout) :: self
               logical, optional, intent(in) :: dealloc_rst   ! default is .false. .true. => free stacks, mrst = 0
            end subroutine
         end interface
      
      !-------------------------------------------------
      !-----     NEK EXTENDED REAL VECTOR TYPE     -----
      !-------------------------------------------------
      
      ! --> Type.
         type, extends(abstract_vector_rdp), public :: nek_ext_dvector
            ! public data
            real(kind=dp), allocatable :: vx(:), vy(:), vz(:)
            real(kind=dp), allocatable :: pr(:)
            real(kind=dp), allocatable :: theta(:,:)
            real(kind=dp) :: T
            ! private state and restart storage
            real(kind=dp), allocatable, private :: vxrst(:,:), vyrst(:,:), vzrst(:,:)
            real(kind=dp), allocatable, private :: prrst(:,:)
            real(kind=dp), allocatable, private :: thetarst(:,:,:)
            real(kind=dp), allocatable, private :: Trst(:)
            integer, private :: nrst = 0
            logical, private :: has_rst = .false.
            logical, private :: initialized = .false.
         contains
            private
            ! required basic type-bound procedures
            procedure, pass(self), public :: zero     => nek_ext_dzero
            procedure, pass(self), public :: rand     => nek_ext_drand
            procedure, pass(self), public :: scal     => nek_ext_dscal
            procedure, pass(self), public :: axpby    => nek_ext_daxpby
            procedure, pass(self), public :: dot      => nek_ext_ddot
            procedure, pass(self), public :: get_size => nek_ext_dsize
            ! initialization
            procedure, pass(self), public :: init           => nek_ext_dinit
            procedure, pass(self), public :: is_initialized => nek_ext_dis_initialized
            ! additional type-bound procedures
            procedure, pass(self), public :: save_rst         => ext_dsave_rst
            procedure, pass(self), public :: get_rst          => ext_dget_rst
            procedure, pass(self), public :: has_rst_fields   => ext_dhas_rst_fields
            procedure, pass(self), public :: has_rst_storage  => ext_dhas_rst_storage
            procedure, pass(self), public :: clear_rst_fields => ext_dclear_rst_fields
         end type nek_ext_dvector
      
      ! --> Constructor.
         interface nek_ext_dvector
            module function construct_nek_ext_dvector(vx, vy, vz, pr, theta, T) result(out)
               real(kind=dp), dimension(lv), intent(in) :: vx, vy
               real(kind=dp), dimension(lv), optional, intent(in) :: vz
               real(kind=dp), dimension(lp), optional, intent(in) :: pr
               real(kind=dp), dimension(lv, ldimt), optional, intent(in) :: theta
               real(kind=dp), optional, intent(in) :: T
               type(nek_ext_dvector) :: out
            end function
         end interface
      
      ! --> Type-bound procedures.
         interface

            ! required by LightKrylov abstract_vector_rdp
            module subroutine nek_ext_dzero(self)
               class(nek_ext_dvector), intent(inout) :: self
            end subroutine
      
            module subroutine nek_ext_drand(self, ifnorm)
               class(nek_ext_dvector), intent(inout) :: self
               logical, optional, intent(in) :: ifnorm
            end subroutine
      
            module subroutine nek_ext_dscal(self, alpha)
               class(nek_ext_dvector), intent(inout) :: self
               real(kind=dp), intent(in) :: alpha
            end subroutine
      
            module subroutine nek_ext_daxpby(alpha, vec, beta, self)
               class(nek_ext_dvector), intent(inout) :: self
               real(kind=dp), intent(in) :: alpha
               class(abstract_vector_rdp), intent(in) :: vec
               real(kind=dp), intent(in) :: beta
            end subroutine
      
            real(kind=dp) module function nek_ext_ddot(self, vec) result(alpha)
               class(nek_ext_dvector), intent(in) :: self
               class(abstract_vector_rdp), intent(in) :: vec
            end function
      
            integer pure module function nek_ext_dsize(self) result(n)
               class(nek_ext_dvector), intent(in) :: self
            end function

            ! initialization
            module subroutine nek_ext_dinit(self, alloc_rst)
               class(nek_ext_dvector), intent(inout) :: self
               logical, optional, intent(in) :: alloc_rst ! default is .false. (no restart storage)
            end subroutine

            logical pure module function nek_ext_dis_initialized(self) result(ok)
               class(nek_ext_dvector), intent(in) :: self
            end function

            ! restart storage and organization
            module subroutine ext_dsave_rst(self, vec_rst, irst)
               class(nek_ext_dvector), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_rst
               integer, intent(in) :: irst
            end subroutine

            module subroutine ext_dget_rst(self, vec_rst, irst)
               class(nek_ext_dvector), intent(in) :: self
               class(abstract_vector_rdp), intent(inout) :: vec_rst
               integer, intent(in) :: irst
            end subroutine
            
            logical pure module function ext_dhas_rst_fields(self) result(has) ! nrst > 0
               class(nek_ext_dvector), intent(in) :: self
            end function

            logical pure module function ext_dhas_rst_storage(self) result(has)
               class(nek_ext_dvector), intent(in) :: self
            end function
            
            module subroutine ext_dclear_rst_fields(self, dealloc_rst)
               class(nek_ext_dvector), intent(inout) :: self
               logical, optional, intent(in) :: dealloc_rst   ! default is .false. .true. => free stacks, mrst = 0
            end subroutine
         end interface
      
      !-------------------------------------------
      !-----     NEK COMPLEX VECTOR TYPE     -----
      !-------------------------------------------
      
      ! --> Type.
         type, extends(abstract_vector_cdp), public :: nek_zvector
            type(nek_dvector) :: re
            type(nek_dvector) :: im
         contains
            private
            ! required basic type-bound procedures
            ! required basic type-bound procedures
            procedure, pass(self), public :: zero     => nek_zzero
            procedure, pass(self), public :: rand     => nek_zrand
            procedure, pass(self), public :: scal     => nek_zscal
            procedure, pass(self), public :: axpby    => nek_zaxpby
            procedure, pass(self), public :: dot      => nek_zdot
            procedure, pass(self), public :: get_size => nek_zsize
            ! initialization
            procedure, pass(self), public :: init           => nek_zinit
            procedure, pass(self), public :: is_initialized => nek_zis_initialized
            ! additional type-bound procedures
            procedure, pass(self), public :: save_rst         => zsave_rst
            procedure, pass(self), public :: get_rst          => zget_rst
            procedure, pass(self), public :: has_rst_fields   => zhas_rst_fields
            procedure, pass(self), public :: has_rst_storage  => zhas_rst_storage
            procedure, pass(self), public :: clear_rst_fields => zclear_rst_fields
         end type nek_zvector
      
      ! --> Constructor.
         interface nek_zvector
            module function construct_nek_zvector(vx, vy, vz, pr, theta) result(out)
               real(kind=dp), dimension(lv, 2), intent(in) :: vx, vy
               real(kind=dp), dimension(lv, 2), optional, intent(in) :: vz
               real(kind=dp), dimension(lp, 2), optional, intent(in) :: pr
               real(kind=dp), dimension(lv, ldimt, 2), optional, intent(in) :: theta
               type(nek_zvector) :: out
            end function
         end interface
      
      ! --> Type-bound procedures.
         interface

            ! required by LightKrylov abstract_vector_cdp
            module subroutine nek_zzero(self)
               class(nek_zvector), intent(inout) :: self
            end subroutine
      
            module subroutine nek_zrand(self, ifnorm)
               class(nek_zvector), intent(inout) :: self
               logical, optional, intent(in) :: ifnorm
            end subroutine
      
            module subroutine nek_zscal(self, alpha)
               class(nek_zvector), intent(inout) :: self
               complex(kind=dp), intent(in) :: alpha
            end subroutine
      
            module subroutine nek_zaxpby(alpha, vec, beta, self)
               class(nek_zvector), intent(inout) :: self
               complex(kind=dp), intent(in) :: alpha
               class(abstract_vector_cdp), intent(in) :: vec
               complex(kind=dp), intent(in) :: beta
            end subroutine
      
            complex(kind=dp) module function nek_zdot(self, vec) result(alpha)
               class(nek_zvector), intent(in) :: self
               class(abstract_vector_cdp), intent(in) :: vec
            end function
      
            integer pure module function nek_zsize(self) result(n)
               class(nek_zvector), intent(in) :: self
            end function

            ! initialization
            module subroutine nek_zinit(self, alloc_rst)
               class(nek_zvector), intent(inout) :: self
               logical, optional, intent(in) :: alloc_rst ! default is .false. (no restart storage)
            end subroutine

            logical pure module function nek_zis_initialized(self) result(ok)
               class(nek_zvector), intent(in) :: self
            end function

            ! restart storage and organization
            module subroutine zsave_rst(self, vec_rst, irst)
               class(nek_zvector), intent(inout) :: self
               class(abstract_vector_cdp), intent(in) :: vec_rst
               integer, intent(in) :: irst
            end subroutine

            module subroutine zget_rst(self, vec_rst, irst)
               class(nek_zvector), intent(in) :: self
               class(abstract_vector_cdp), intent(inout) :: vec_rst
               integer, intent(in) :: irst
            end subroutine

            logical pure module function zhas_rst_fields(self) result(has) ! nrst > 0
               class(nek_zvector), intent(in) :: self
            end function

            logical pure module function zhas_rst_storage(self) result(has)
               class(nek_zvector), intent(in) :: self
            end function

            module subroutine zclear_rst_fields(self, dealloc_rst)
               class(nek_zvector), intent(inout) :: self
               logical, optional, intent(in) :: dealloc_rst ! default is .false. .true. => free stacks
            end subroutine
         end interface
      
      contains
      
         real(kind=dp) function mth_rand(ix, iy, iz, ieg, xl, fcoeff) !generate random number
            include 'INPUT'           ! IF3D
            integer ix, iy, iz, ieg
            real(kind=dp) xl(LDIM), fcoeff(3)
            mth_rand = fcoeff(1)*(ieg + xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy + fcoeff(3)*ix
            if (IF3D) mth_rand = fcoeff(1)*(ieg + xl(NDIM)*sin(mth_rand)) + fcoeff(2)*iz*ix + fcoeff(3)*iz
            mth_rand = 1.0e3_dp*sin(mth_rand)
            mth_rand = 1.0e3_dp*sin(mth_rand)
            mth_rand = cos(mth_rand)
         end function mth_rand
      
      end module neklab_vectors
