      module neklab_pvectors
         use stdlib_optval, only: optval
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
      ! Abstract types for real-valued vectors.
         use LightKrylov, only: abstract_vector_rdp
      
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 'neklab_pvectors'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
      !----------------------------------------
      !-----     NEK REAL VECTOR TYPE     -----
      !----------------------------------------
      
         type, extends(abstract_vector_rdp), public :: nek_pdvector
            real(kind=dp), dimension(lp) :: pr
      !! Pressure components of the state vector.
         contains
            private
            procedure, pass(self), public :: zero => nek_dzero
      !! Sets a vector to zero.
            procedure, pass(self), public :: rand => nek_drand
      !! Create a random vector.
            procedure, pass(self), public :: scal => nek_dscal
      !! Scale a vector such that \( \mathbf{x} = \alpha \mathbf{x}\) with \(\alpha \in \mathbb{R} \).
            procedure, pass(self), public :: axpby => nek_daxpby
      !! Add (in-place) two vectors such that \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \) with \( \alpha \) and \( \beta \in \mathbb{R} \).
            procedure, pass(self), public :: dot => nek_ddot
      !! Compute the \( \ell_2 \) inner-product between two vectors.
            procedure, pass(self), public :: get_size => nek_dsize
      !! Return the size of the vector
         end type nek_pdvector
      
      contains
      
      !----------------------------------------------------------------------------
      !-----                                                                  -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES FOR REAL VECTORS     -----
      !-----                                                                  -----
      !----------------------------------------------------------------------------
      
         subroutine nek_dzero(self)
            class(nek_pdvector), intent(inout) :: self
      !! Vector to be zeroed-out.
            call self%scal(0.0_dp)
            return
         end subroutine nek_dzero
      
         subroutine nek_drand(self, ifnorm)
            class(nek_pdvector), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
            logical :: normalize
            integer :: i, n, ieg, iel
            real(kind=dp) :: xl(ldim), fcoeff(ldim), alpha
      
            normalize = optval(ifnorm, .false.)
      
            do i = 1, lp
               ieg = lglel(iel)
               xl(1) = xm2(i, 1, 1, 1)
               xl(2) = ym2(i, 1, 1, 1)
               if (if3D) xl(ldim) = zm2(i, 1, 1, 1)
         
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%pr(i) = self%pr(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
            end do
      
            if (normalize) then
               alpha = self%norm()
               call self%scal(1.0_dp/alpha)
            end if
      
            return
         end subroutine nek_drand
      
         subroutine nek_dscal(self, alpha)
            class(nek_pdvector), intent(inout) :: self
            real(kind=dp), intent(in) :: alpha
            call dscal(lp, alpha, self%pr, 1)
            return
         end subroutine nek_dscal
      
         subroutine nek_daxpby(self, alpha, vec, beta)
            class(nek_pdvector), intent(inout) :: self
            real(kind=dp), intent(in) :: alpha
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), intent(in) :: beta
            call self%scal(alpha)
            select type (vec)
            type is (nek_pdvector)
               call daxpy(lp, beta, vec%pr, 1, self%pr, 1)
            end select
            return
         end subroutine nek_daxpby
      
         real(kind=dp) function nek_ddot(self, vec) result(alpha)
            class(nek_pdvector), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), external :: glsc2
      
            select type (vec)
            type is (nek_pdvector)
               !alpha = glsc3(self%pr, vec%pr, bm2, lp)
               !alpha = glsc3(self%pr, bm2inv, vec%pr, lp)/volvm2    ! rnorm from convprn in navier1.f
               alpha = glsc2(self%pr, vec%pr, lp)                    ! convprn in navier1.f
            end select
      
            return
         end function nek_ddot
      
         integer pure function nek_dsize(self) result(N)
            class(nek_pdvector), intent(in) :: self
            N = lp
            return
         end function nek_dsize
      
      !---------------------------------------
      !
      !       UTILITY for random vectors
      !
      !---------------------------------------
      
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
      
      end module neklab_pvectors
