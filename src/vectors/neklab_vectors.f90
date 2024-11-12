      module neklab_vectors
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use LightKrylov, only: abstract_vector_rdp
      
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 'neklab_vectors'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      
         !----------------------------------------
         !-----     NEK REAL VECTOR TYPE     -----
         !----------------------------------------
      
         ! --> Type.
         type, extends(abstract_vector_rdp), public :: nek_dvector
            real(kind=dp), dimension(lv) :: vx, vy, vz
            real(kind=dp), dimension(lp) :: pr
            real(kind=dp), dimension(lv, ldimt) :: theta
         contains
            private
            procedure, pass(self), public :: zero => nek_dzero
            procedure, pass(self), public :: rand => nek_drand
            procedure, pass(self), public :: scal => nek_dscal
            procedure, pass(self), public :: axpby => nek_daxpby
            procedure, pass(self), public :: dot => nek_ddot
            procedure, pass(self), public :: get_size => nek_dsize
         end type nek_dvector

         ! --> Constructor.
         interface nek_dvector
            pure module function construct_nek_dvector(vx, vy, vz, pr, theta) result(out)
               real(dp), dimension(lv), intent(in) :: vx, vy
               real(dp), dimension(lv), optional, intent(in) :: vz
               real(dp), dimension(lp), optional, intent(in) :: pr
               real(dp), dimension(lv, ldimt), optional, intent(in) :: theta
               type(nek_dvector) :: out
            end function
         end interface
      
         ! --> Type-bound procedures.
         interface
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
      
            module subroutine nek_daxpby(self, alpha, vec, beta)
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
         end interface

         !-------------------------------------------------
         !-----     NEK EXTENDED REAL VECTOR TYPE     -----
         !-------------------------------------------------

         ! --> Type.
         type, extends(abstract_vector_rdp), public :: nek_ext_dvector
            real(kind=dp), dimension(lv) :: vx, vy, vz
            real(kind=dp), dimension(lp) :: pr
            real(kind=dp), dimension(lv, ldimt) :: theta
            real(kind=dp) :: T
         contains
            private
            procedure, pass(self), public :: zero => nek_ext_dzero
            procedure, pass(self), public :: rand => nek_ext_drand
            procedure, pass(self), public :: scal => nek_ext_dscal
            procedure, pass(self), public :: axpby => nek_ext_daxpby
            procedure, pass(self), public :: dot => nek_ext_ddot
            procedure, pass(self), public :: get_size => nek_ext_dsize
         end type nek_ext_dvector

         ! --> Constructor.
         interface nek_ext_dvector
            pure module function construct_nek_ext_dvector(vx, vy, vz, pr, theta, time) result(out)
               real(dp), dimension(lv), intent(in) :: vx, vy
               real(dp), dimension(lv), optional, intent(in) :: vz
               real(dp), dimension(lp), optional, intent(in) :: pr
               real(dp), dimension(lv, ldimt), optional, intent(in) :: theta
               real(dp), optional, intent(in) :: time
               type(nek_ext_dvector) :: out
            end function
         end interface

         ! --> Type-bound procedures.
         interface
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
      
            module subroutine nek_ext_daxpby(self, alpha, vec, beta)
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
            return
         end function mth_rand
      
      end module neklab_vectors
