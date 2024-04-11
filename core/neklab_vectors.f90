      module neklab_vectors
      use LightKrylov, only: abstract_vector, wp
      implicit none
      include "SIZE"
      include "TOTAL"
      
      private
      
      public :: nopcopy
      
      integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
      integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
      type, extends(abstract_vector), public :: nek_dvector
      !! Type definition for Nek5000 state-vector (real).
      real(kind=wp), dimension(lv) :: vx, vy, vz
      !! Velocity components of the state vector.
      real(kind=wp), dimension(lp) :: pr
      !! Pressure components of the state vector.
      real(kind=wp), dimension(lv, ldimt) :: theta
      !! Temperature + passive scalars components.
      contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: dot
      end type nek_dvector
      
      contains
      
      !-----------------------------------------------------------
      !-----                                                 -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES     -----
      !-----                                                 -----
      !-----------------------------------------------------------
      
      subroutine zero(self)
      class(nek_dvector), intent(inout) :: self
      !! Vector to be zeroed-out.
      self%vx = 0.0_wp
      self%vy = 0.0_wp
      self%vz = 0.0_wp
      self%pr = 0.0_wp
      self%theta = 0.0_wp
      return
      end subroutine zero
      
      subroutine rand(self, ifnorm)
      class(nek_dvector), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      integer :: i, n, ieg, iel
      real(kind=wp) :: xl(ldim), fcoeff(3), alpha
      
      n = nx1*ny1*nz1*nelv
      do i = 1, n
      ieg = lglel(iel)
      xl(1) = xm1(i, 1, 1, 1)
      xl(2) = ym1(i, 1, 1, 1)
      if (if3d) xl(3) = zm1(i, 1, 1, 1)
      
      fcoeff(1) = 3.0e4_wp; fcoeff(2) = -1.5e3_wp ; fcoeff(3) = 0.5e5_wp
      self%vx(i) = self%vx(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      
      fcoeff(1) = 2.3e4_wp; fcoeff(2) = 2.3e3_wp; fcoeff(3) = -2.0e5_wp
      self%vy(i) = self%vy(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      enddo
      ! Face averaging.
      call opdssum(self%vx, self%vy, self%vz)
      call opcolv(self%vx, self%vy, self%vz, vmult)
      call dsavg(self%vx)
      call dsavg(self%vy)
      if (if3d) call dsavg(self%vz)
      call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
      return
      end subroutine rand
      
      subroutine scal(self, alpha)
      class(nek_dvector), intent(inout) :: self
      real(kind=wp), intent(in) :: alpha
      self%vx = alpha * self%vx
      self%vy = alpha * self%vy
      self%vz = alpha * self%vz
      self%pr = alpha * self%pr
      self%theta = alpha * self%theta
      return
      end subroutine scal
      
      subroutine axpby(self, alpha, vec, beta)
      class(nek_dvector), intent(inout) :: self
      real(kind=wp), intent(in) :: alpha
      class(abstract_vector), intent(in) :: vec
      real(kind=wp), intent(in) :: beta
      select type(vec)
      type is(nek_dvector)
      self%vx = alpha*self%vx + beta*vec%vx
      self%vy = alpha*self%vy + beta*vec%vy
      self%vz = alpha*self%vz + beta*vec%vz
      self%pr = alpha*self%pr + beta*vec%pr
      self%theta = alpha*self%theta + beta*vec%theta
      end select
      return
      end subroutine axpby
      
      real(kind=wp) function dot(self, vec) result(alpha)
      class(nek_dvector), intent(in) :: self
      class(abstract_vector), intent(in) :: vec
      real(kind=wp), external :: glsc3
      integer :: i, n
      
      n = nx1*ny1*nz1*nelv
      select type (vec)
      type is (nek_dvector)
      ! Kinetic energy contribution.
      alpha = glsc3(self%vx, bm1, vec%vx, n) + glsc3(self%vy, bm1, vec%vy, n)
      if (if3d) then
      alpha = alpha + glsc3(self%vz, bm1, vec%vz)
      end if
      
      ! Thermal energy contribution.
      if (ifto) then
      alpha = alpha + glsc3(self%theta(:, 1), bm1, vec%theta(:, 1), n)
      end if
      
      ! Whatever contribution from additional scalars.
      if (ldimt > 1) then
      do i = 2, ldimt
      if (ifpsco(i-1)) alpha = alpha + glsc3(self%theta(:, i), bm1, vec%theta(:, i), n)
      enddo
      endif
      end select
      
      return
      end function dot
      
      !-----------------------------------------
      !-----                               -----
      !-----     NEK-RELATED UTILITIES     -----
      !-----                               -----
      !-----------------------------------------
      
      subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      integer :: n, k
      real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
      real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
      n = nx1*ny1*nz1*nelv
      call copy(a1, b1, n)
      call copy(a2, b2, n)
      if (if3d) call copy(a3, b3, n)
      if (ifpo) call copy(a4, b4, nx2*ny2*nz2*nelv)
      if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
      do k = 1, npscal
      if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
      end do
      end if
      return
      end subroutine nopcopy
      
      real(kind=wp) function mth_rand(ix,iy,iz,ieg,xl,fcoeff) !generate random number
      include 'INPUT'           ! IF3D
      integer ix,iy,iz,ieg
      real xl(LDIM), fcoeff(3)
      mth_rand = fcoeff(1)*(ieg+xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy+fcoeff(3)*ix
      if (IF3D) mth_rand = fcoeff(1)*(ieg +xl(NDIM)*sin(mth_rand))+fcoeff(2)*iz*ix+fcoeff(3)*iz
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = cos(mth_rand)
      return
      end function mth_rand
      end module neklab_vectors
