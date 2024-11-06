      module neklab_vectors
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
         character(len=*), parameter :: this_module = 'neklab_vectors'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
      !----------------------------------------
      !-----     NEK REAL VECTOR TYPE     -----
      !----------------------------------------
      
         type, extends(abstract_vector_rdp), public :: nek_dvector
      !! Type definition for Nek5000 state-vector (real).
            real(kind=dp), dimension(lv) :: vx, vy, vz
      !! Velocity components of the state vector.
            real(kind=dp), dimension(lp) :: pr
      !! Pressure components of the state vector.
            real(kind=dp), dimension(lv, ldimt) :: theta
      !! Temperature + passive scalars components.
         contains
            private
            procedure, pass(self), public :: zero => nek_dzero
      !! Sets a vector to zero.
            procedure, pass(self), public :: rand => nek_drand
      !! Create a random vector.
            procedure, pass(self), public :: scal => nek_dscal
      !! Scale a vector such that \( \mathbf{x} = \alpha \mathbf{x}$ with     $         \alpha \in \mathbb{R} \).
            procedure, pass(self), public :: axpby => nek_daxpby
      !! Add (in-place) two vectors such that \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \) with \( \alpha \) and \( \beta \in \mathbb{R} \).
            procedure, pass(self), public :: dot => nek_ddot
      !! Compute the \( \ell_2 \) inner-product between two vectors.
            procedure, pass(self), public :: get_size => nek_dsize
      !! Return the size of the vector
         end type nek_dvector

         type, extends(abstract_vector_rdp), public :: nek_ext_dvector
      !! Type definition for Nek5000 state-vector (real).
            real(kind=dp), dimension(lv) :: vx, vy, vz
      !! Velocity components of the state vector.
            real(kind=dp), dimension(lp) :: pr
      !! Pressure components of the state vector.
            real(kind=dp), dimension(lv, ldimt) :: theta
      !! Temperature + passive scalars components.
            real(kind=dp) :: T
      !! Period of the orbit
         contains
            private
            procedure, pass(self), public :: zero => nek_ext_dzero
      !! Sets a vector to zero.
            procedure, pass(self), public :: rand => nek_ext_drand
      !! Create a random vector.
            procedure, pass(self), public :: scal => nek_ext_dscal
      !! Scale a vector such that \( \mathbf{x} = \alpha \mathbf{x}$ with     $         \alpha \in \mathbb{R} \).
            procedure, pass(self), public :: axpby => nek_ext_daxpby
      !! Add (in-place) two vectors such that \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \) with \( \alpha \) and \( \beta \in \mathbb{R} \).
            procedure, pass(self), public :: dot => nek_ext_ddot
      !! Compute the \( \ell_2 \) inner-product between two vectors.
            procedure, pass(self), public :: get_size => nek_ext_dsize
      !! Return the size of the vector
         end type nek_ext_dvector
      
      contains
      
      !----------------------------------------------------------------------------
      !-----                                                                  -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES FOR REAL VECTORS     -----
      !-----                                                                  -----
      !----------------------------------------------------------------------------
      
         subroutine nek_dzero(self)
            class(nek_dvector), intent(inout) :: self
      !! Vector to be zeroed-out.
            call self%scal(0.0_dp)
            return
         end subroutine nek_dzero
      
         subroutine nek_drand(self, ifnorm)
            class(nek_dvector), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
            logical :: normalize
            integer :: i, n, ieg, iel
            real(kind=dp) :: xl(ldim), fcoeff(3), alpha
      
            n = nx1*ny1*nz1*nelv
            normalize = optval(ifnorm, .false.)
      
            do i = 1, n
               ieg = lglel(iel)
               xl(1) = xm1(i, 1, 1, 1)
               xl(2) = ym1(i, 1, 1, 1)
               if (if3d) xl(3) = zm1(i, 1, 1, 1)
      
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%vx(i) = self%vx(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%vy(i) = self%vy(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
            end do
      
      ! Face averaging.
            call opdssum(self%vx, self%vy, self%vz)
            call opcolv(self%vx, self%vy, self%vz, vmult)
            call dsavg(self%vx)
            call dsavg(self%vy)
            if (if3d) call dsavg(self%vz)
            call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
            if (normalize) then
               alpha = self%norm()
               call self%scal(1.0_dp/alpha)
            end if
      
            return
         end subroutine nek_drand
      
         subroutine nek_dscal(self, alpha)
            class(nek_dvector), intent(inout) :: self
            real(kind=dp), intent(in) :: alpha
            call dscal(lv, alpha, self%vx, 1)
            call dscal(lv, alpha, self%vy, 1)
            if (if3d) call dscal(lv, alpha, self%vz, 1)
            call dscal(lp, alpha, self%pr, 1)
            if (ifto) call dscal(lv, alpha, self%theta(:, 1), 1)
            return
         end subroutine nek_dscal
      
         subroutine nek_daxpby(self, alpha, vec, beta)
            class(nek_dvector), intent(inout) :: self
            real(kind=dp), intent(in) :: alpha
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), intent(in) :: beta
            call self%scal(alpha)
            select type (vec)
            type is (nek_dvector)
               call daxpy(lv, beta, vec%vx, 1, self%vx, 1)
               call daxpy(lv, beta, vec%vy, 1, self%vy, 1)
               if (if3d) call daxpy(lv, beta, vec%vz, 1, self%vz, 1)
               call daxpy(lp, beta, vec%pr, 1, self%pr, 1)
               if (ifto) call daxpy(lv, beta, vec%theta(:, 1), 1, self%theta(:, 1), 1)
            end select
            return
         end subroutine nek_daxpby
      
         real(kind=dp) function nek_ddot(self, vec) result(alpha)
            class(nek_dvector), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), external :: glsc3
            integer :: i
      
            select type (vec)
            type is (nek_dvector)
            ! Kinetic energy contribution.
               alpha = glsc3(self%vx, bm1, vec%vx, lv) + glsc3(self%vy, bm1, vec%vy, lv)
               if (if3d) then
                  alpha = alpha + glsc3(self%vz, bm1, vec%vz)
               end if
      
            ! Thermal energy contribution.
               if (ifto) then
                  alpha = alpha + glsc3(self%theta(:, 1), bm1, vec%theta(:, 1), lv)
               end if
      
            ! Whatever contribution from additional scalars.
               if (ldimt > 1) then
               do i = 2, ldimt
                  if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), bm1, vec%theta(:, i), lv)
               end do
               end if
            end select
      
            return
         end function nek_ddot

         integer pure function nek_dsize(self) result(N)
            class(nek_dvector), intent(in) :: self
            integer :: i
            N = 2*lv + lp
            if (if3d) N = N + lv
            if (ifto) N = N + lv
            if (ldimt > 1) then
               do i = 2, ldimt
                  if (ifpsco(i - 1)) N = N + lv
               end do
            end if
            return
         end function nek_dsize

      !-------------------------------------------------------------------------------------
      !-----                                                                           -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES FOR EXTENDED REAL VECTORS     -----
      !-----                                                                           -----
      !-------------------------------------------------------------------------------------
      
         subroutine nek_ext_dzero(self)
            class(nek_ext_dvector), intent(inout) :: self
      !! Vector to be zeroed-out.
            call self%scal(0.0_dp)
            return
         end subroutine nek_ext_dzero
      
         subroutine nek_ext_drand(self, ifnorm)
            class(nek_ext_dvector), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
            logical :: normalize
            integer :: i, n, ieg, iel
            real(kind=dp) :: xl(ldim), fcoeff(3), alpha
      
            n = nx1*ny1*nz1*nelv
            normalize = optval(ifnorm, .false.)
      
            do i = 1, n
               ieg = lglel(iel)
               xl(1) = xm1(i, 1, 1, 1)
               xl(2) = ym1(i, 1, 1, 1)
               if (if3d) xl(3) = zm1(i, 1, 1, 1)
      
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%vx(i) = self%vx(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%vy(i) = self%vy(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
            end do
      
      ! Face averaging.
            call opdssum(self%vx, self%vy, self%vz)
            call opcolv(self%vx, self%vy, self%vz, vmult)
            call dsavg(self%vx)
            call dsavg(self%vy)
            if (if3d) call dsavg(self%vz)
            call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
            call random_number(self%T)

            if (normalize) then
               alpha = self%norm()
               call self%scal(1.0_dp/alpha)
            end if
      
            return
         end subroutine nek_ext_drand
      
         subroutine nek_ext_dscal(self, alpha)
            class(nek_ext_dvector), intent(inout) :: self
            real(kind=dp), intent(in) :: alpha
            call dscal(lv, alpha, self%vx, 1)
            call dscal(lv, alpha, self%vy, 1)
            if (if3d) call dscal(lv, alpha, self%vz, 1)
            call dscal(lp, alpha, self%pr, 1)
            if (ifto) call dscal(lv, alpha, self%theta(:, 1), 1)
            self%T = alpha*self%T
            return
         end subroutine nek_ext_dscal
      
         subroutine nek_ext_daxpby(self, alpha, vec, beta)
            class(nek_ext_dvector), intent(inout) :: self
            real(kind=dp), intent(in) :: alpha
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), intent(in) :: beta
            call self%scal(alpha)
            select type (vec)
            type is (nek_ext_dvector)
               call daxpy(lv, beta, vec%vx, 1, self%vx, 1)
               call daxpy(lv, beta, vec%vy, 1, self%vy, 1)
               if (if3d) call daxpy(lv, beta, vec%vz, 1, self%vz, 1)
               call daxpy(lp, beta, vec%pr, 1, self%pr, 1)
               if (ifto) call daxpy(lv, beta, vec%theta(:, 1), 1, self%theta(:, 1), 1)
               self%T = alpha*self%T + vec%T
            end select
            return
         end subroutine nek_ext_daxpby
      
         real(kind=dp) function nek_ext_ddot(self, vec) result(alpha)
            class(nek_ext_dvector), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), external :: glsc3
            integer :: i
      
            select type (vec)
            type is (nek_ext_dvector)
            ! Kinetic energy contribution.
               alpha = glsc3(self%vx, bm1, vec%vx, lv) + glsc3(self%vy, bm1, vec%vy, lv)
               if (if3d) then
                  alpha = alpha + glsc3(self%vz, bm1, vec%vz)
               end if
      
            ! Thermal energy contribution.
               if (ifto) then
                  alpha = alpha + glsc3(self%theta(:, 1), bm1, vec%theta(:, 1), lv)
               end if
      
            ! Whatever contribution from additional scalars.
               if (ldimt > 1) then
               do i = 2, ldimt
                  if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), bm1, vec%theta(:, i), lv)
               end do
               end if
               alpha = alpha + self%T*vec%T
            end select
      
            return
         end function nek_ext_ddot

         integer pure function nek_ext_dsize(self) result(N)
            class(nek_ext_dvector), intent(in) :: self
            integer :: i
            N = 2*lv + lp + 1
            if (if3d) N = N + lv
            if (ifto) N = N + lv
            if (ldimt > 1) then
               do i = 2, ldimt
                  if (ifpsco(i - 1)) N = N + lv
               end do
            end if
            return
         end function nek_ext_dsize

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
      
      end module neklab_vectors