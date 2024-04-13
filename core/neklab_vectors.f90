      module neklab_vectors
      use LightKrylov, only: abstract_vector, wp
      implicit none
      include "SIZE"
      include "TOTAL"
      
      private
      
      public :: nek2vec, vec2nek
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
      procedure, pass(self), public :: zero  => nek_dzero
      procedure, pass(self), public :: rand  => nek_drand
      procedure, pass(self), public :: scal  => nek_dscal
      procedure, pass(self), public :: axpby => nek_daxpby
      procedure, pass(self), public :: dot   => nek_ddot
      end type nek_dvector
      
      interface nek2vec
      module procedure nek2vec_prt
      module procedure nek2vec_std
      end interface
      
      interface vec2nek
      module procedure vec2nek_std
      module procedure vec2nek_prt
      end interface
      
      contains
      
      !-----------------------------------------------------------
      !-----                                                 -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES     -----
      !-----                                                 -----
      !-----------------------------------------------------------
      
      subroutine nek_dzero(self)
      class(nek_dvector), intent(inout) :: self
      !! Vector to be zeroed-out.
      self%vx = 0.0_wp
      self%vy = 0.0_wp
      self%vz = 0.0_wp
      self%pr = 0.0_wp
      self%theta = 0.0_wp
      return
      end subroutine nek_dzero
      
      subroutine nek_drand(self, ifnorm)
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
      
      fcoeff(1) = 3.0e4_wp; fcoeff(2) = -1.5e3_wp; fcoeff(3) = 0.5e5_wp
      self%vx(i) = self%vx(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      
      fcoeff(1) = 2.3e4_wp; fcoeff(2) = 2.3e3_wp; fcoeff(3) = -2.0e5_wp
      self%vy(i) = self%vy(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      end do
      ! Face averaging.
      call opdssum(self%vx, self%vy, self%vz)
      call opcolv(self%vx, self%vy, self%vz, vmult)
      call dsavg(self%vx)
      call dsavg(self%vy)
      if (if3d) call dsavg(self%vz)
      call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
      return
      end subroutine nek_drand
      
      subroutine nek_dscal(self, alpha)
      class(nek_dvector), intent(inout) :: self
      real(kind=wp), intent(in) :: alpha
      call dscal(lv, alpha, self%vx, 1)
      call dscal(lv, alpha, self%vy, 1)
      if (if3d) call dscal(lv, alpha, self%vz, 1)
      call dscal(lp, alpha, self%pr, 1)
      if (ifto) call dscal(lv, alpha, self%theta(:, 1), 1)
      return
      end subroutine nek_dscal
      
      subroutine nek_daxpby(self, alpha, vec, beta)
      class(nek_dvector), intent(inout) :: self
      real(kind=wp), intent(in) :: alpha
      class(abstract_vector), intent(in) :: vec
      real(kind=wp), intent(in) :: beta
      select type (vec)
      type is (nek_dvector)
      call daxpby(lv, beta, vec%vx, 1, alpha, self%vx, 1)
      call daxpby(lv, beta, vec%vy, 1, alpha, self%vy, 1)
      if (if3d) call daxpby(lv, beta, vec%vz, 1, alpha, self%vz, 1)
      call daxpby(lp, beta, vec%pr, 1, alpha, self%pr, 1)
      if (ifto) call daxpby(lv, beta, vec%theta(:, 1), 1, alpha, self%theta(:, 1), 1)
      end select
      return
      end subroutine nek_daxpby
      
      real(kind=wp) function nek_ddot(self, vec) result(alpha)
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
      if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), bm1, vec%theta(:, i), n)
      end do
      end if
      end select
      
      return
      end function nek_ddot
      
      !-----------------------------------------
      !-----                               -----
      !-----     NEK-RELATED UTILITIES     -----
      !-----                               -----
      !-----------------------------------------
      subroutine nek2vec_prt(vec, vx_, vy_, vz_, pr_, t_)
      include "SIZE"
      type(nek_dvector), intent(out) :: vec
      real(kind=wp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vx_
      real(kind=wp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vy_
      real(kind=wp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vz_
      real(kind=wp), dimension(lx2*ly2*lz2*lelv, lpert), intent(in) :: pr_
      real(kind=wp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
      call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
      return
      end subroutine nek2vec_prt
      
      subroutine nek2vec_std(vec, vx_, vy_, vz_, pr_, t_)
      include "SIZE"
      type(nek_dvector), intent(out) :: vec
      real(kind=wp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
      real(kind=wp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
      real(kind=wp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
      real(kind=wp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
      real(kind=wp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
     
      call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
      
      return
      end subroutine nek2vec_std
      
      subroutine vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
      include "SIZE"
      type(nek_dvector), intent(in) :: vec
      real(kind=wp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
      real(kind=wp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
      real(kind=wp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
      real(kind=wp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
      real(kind=wp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
      call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta) 
      
      return
      end subroutine vec2nek_std
      
      subroutine vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
      include "SIZE"
      type(nek_dvector), intent(in) :: vec
      real(kind=wp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vx_
      real(kind=wp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vy_
      real(kind=wp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vz_
      real(kind=wp), dimension(lx2*ly2*lz2*lelv, lpert), intent(out) :: pr_
      real(kind=wp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
      
      call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
      return
      end subroutine vec2nek_prt
      
      subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer :: n, k
      real, intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
      real, intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
      n = nx1*ny1*nz1*nelv
      call copy(a1, b1, n)
      call copy(a2, b2, n)
      if (if3D) call copy(a3, b3, n)
      if (ifpo) call copy(a4, b4, nx2*ny2*nz2*nelv)
      if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
      do k = 1, npscal
      if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
      end do
      end if
      return
      end subroutine nopcopy
      
      real(kind=wp) function mth_rand(ix, iy, iz, ieg, xl, fcoeff) !generate random number
      include 'INPUT'           ! IF3D
      integer ix, iy, iz, ieg
      real xl(LDIM), fcoeff(3)
      mth_rand = fcoeff(1)*(ieg + xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy + fcoeff(3)*ix
      if (IF3D) mth_rand = fcoeff(1)*(ieg + xl(NDIM)*sin(mth_rand)) + fcoeff(2)*iz*ix + fcoeff(3)*iz
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = 1.e3*sin(mth_rand)
      mth_rand = cos(mth_rand)
      return
      end function mth_rand
      end module neklab_vectors
