      module neklab_LNS_utils
         use stdlib_optval, only: optval
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Logger
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_LNS_utils'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
         public :: compute_LNS_conv
         public :: compute_LNS_gradp
         public :: compute_LNS_laplacian
         public :: apply_L
      
      contains
      
         subroutine compute_LNS_conv(c_x, c_y, c_z, uxp, uyp, uzp, trans)
            ! Compute the linearized convective terms
            ! We assume that the baseflow is set in vx, vy, vz
            real(dp), dimension(lv,1), intent(out) :: c_x
            real(dp), dimension(lv,1), intent(out) :: c_y
            real(dp), dimension(lv,1), intent(out) :: c_z
            real(dp), dimension(lv,1), intent(in)  :: uxp
            real(dp), dimension(lv,1), intent(in)  :: uyp
            real(dp), dimension(lv,1), intent(in)  :: uzp
            logical, optional :: trans
            ! internals
            integer :: i
            logical :: transpose
            ! scratch arrays
            real(dp), dimension(lv) :: ta1, ta2, ta3, tb1, tb2, tb3
            transpose = optval(trans, .false.)
            ! (u.grad) Ub
            call opcopy       (tb1,tb2,tb3,vx,vy,vz)          ! Save velocity
            call opcopy       (vx,vy,vz,uxp,uyp,uzp)          ! U <-- u
            if (transpose) then
               call convop_adj(ta1,ta2,ta3,tb1,tb2,tb3,vx,vy,vz)
            else
               call convop    (ta1,tb1)
               call convop    (ta2,tb2)
               if (if3d) then
                  call convop (ta3,tb3)
               else
                  call rzero  (ta3,lv)
               end if
            end if
            call opcopy       (c_x,c_y,c_z,ta1,ta2,ta3) ! copy to output
            ! (Ub.grad) u
            call opcopy       (vx,vy,vz,tb1,tb2,tb3)          ! Restore velocity
            if (transpose) then
               call convop_adj(ta1,ta2,ta3,tb1,tb2,tb3,vx,vy,vz)
            else
               call convop    (tb1,uxp)
               call convop    (tb2,uyp)
               if (if3d) then
                  call convop (tb3,uzp)
               else
                  call rzero  (tb3,lv)
               end if
            end if
            call opadd2       (c_x,c_y,c_z,tb1,tb2,tb3) ! add to output
            return
         end subroutine compute_LNS_conv

         subroutine compute_LNS_laplacian(d_x, d_y, d_z, uxp, uyp, uzp)
            ! compute the laplacian of the input field
            real(dp), dimension(lv,1), intent(out) :: d_x
            real(dp), dimension(lv,1), intent(out) :: d_y
            real(dp), dimension(lv,1), intent(out) :: d_z
            real(dp), dimension(lv,1), intent(in)  :: uxp
            real(dp), dimension(lv,1), intent(in)  :: uyp
            real(dp), dimension(lv,1), intent(in)  :: uzp
            ! internals
            real(dp), dimension(lv) :: h1, h2
            ifield = 1
            call copy(h1, vdiff(1,1,1,1,ifield), lv)
            call rzero(h2, lv)
            ! and apply to the velocity field to compute the diffusion term
            call ophx(d_x, d_y, d_z, uxp, uyp, uzp, h1, h2)
            return
         end subroutine compute_LNS_laplacian

         subroutine compute_LNS_gradp(gp_x, gp_y, gp_z, pp)
            ! compute the laplacian of the input field
            real(dp), dimension(lv,1), intent(out) :: gp_x
            real(dp), dimension(lv,1), intent(out) :: gp_y
            real(dp), dimension(lv,1), intent(out) :: gp_z
            real(dp), dimension(lp,1), intent(in)  :: pp
            ! internals
            real(dp), dimension(lv) :: ta1, ta2, wrk
            ! map the perturbation pressure to the velocity mesh
            call mappr(wrk,pp,ta1,ta2)
            ! compute the gradient on the velocity mesh direclty
            call gradm1(gp_x,gp_y,gp_z,wrk)
            return
         end subroutine compute_LNS_gradp

         subroutine apply_L(Lux, Luy, Luz, ux, uy, uz, pres, trans)
            !! Apply LNS operator (before the projection onto the divergence-free space)
            !! This function assumes that the input vector has already been loaded into v[xyz]p
            real(dp), dimension(lv,1), intent(out) :: Lux
            real(dp), dimension(lv,1), intent(out) :: Luy
            real(dp), dimension(lv,1), intent(out) :: Luz
            real(dp), dimension(lv,1), intent(in)  :: ux
            real(dp), dimension(lv,1), intent(in)  :: uy
            real(dp), dimension(lv,1), intent(in)  :: uz
            real(dp), dimension(lp,1), intent(in)  :: pres
            logical, optional, intent(in) :: trans
            !! adjoint?
            ! internal
            real(dp), dimension(lv) :: utmpx, utmpy, utmpz
      
            ifield = 1
            ! apply BCs
            call bcdirvc(ux,uy,uz,v1mask,v2mask,v3mask)
            
            ! Diffusion term
            call logger%log_information('diffusion term', module=this_module, procedure='compute_LNS')
            call compute_LNS_laplacian(Lux,Luy,Luz,ux,uy,uz)

            ! Pressure gradient
            call logger%log_information('pressure gradient', module=this_module, procedure='compute_LNS')
            call opgradt(utmpx,utmpy,utmpz,pres)

            ! subtract from output terms
            call opsub2(Lux,Luy,Luz,utmpx,utmpy,utmpz)

            ! Convective terms
            call logger%log_information('convective term', module=this_module, procedure='compute_LNS')
            call compute_LNS_conv(utmpx,utmpy,utmpz,ux,uy,uz,trans)
            
            ! subtract from output terms
            call opsub2(Lux,Luy,Luz,utmpx,utmpy,utmpz)
            
            return
         end subroutine apply_L
      
      end module neklab_LNS_utils
