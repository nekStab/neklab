      module neklab_linops
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp, rtol_dp
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: abstract_linop_cdp, abstract_vector_cdp
         use LightKrylov, only: cg, cg_dp_opts, cg_dp_metadata
         use LightKrylov_Logger
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek
         use neklab_nek_setup, only: setup_nonlinear_solver, setup_linear_solver
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_log_information, nek_stop_error
         use neklab_nek_forcing, only: set_neklab_forcing
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         integer, parameter :: lxyz = lx1*ly1*lz1
      
         public :: apply_exptA
         public :: compute_LNS_conv
         public :: compute_LNS_gradp
         public :: compute_LNS_laplacian
         public :: apply_Lv, apply_L
         public :: get_restart, save_restart
      
      !------------------------------------------
      !-----     EXPONENTIAL PROPAGATOR     -----
      !------------------------------------------
      
      ! --> Type.
         type, extends(abstract_linop_rdp), public :: exptA_linop
            real(kind=dp) :: tau
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: init => init_exptA
            procedure, pass(self), public :: matvec => exptA_matvec
            procedure, pass(self), public :: rmatvec => exptA_rmatvec
         end type exptA_linop
      
      ! --> Type-bound procedures: exponential_propagator.f90
         interface
            module subroutine init_exptA(self)
               class(exptA_linop), intent(in) :: self
            end subroutine
      
            module subroutine exptA_matvec(self, vec_in, vec_out)
               class(exptA_linop), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine
      
            module subroutine exptA_rmatvec(self, vec_in, vec_out)
               class(exptA_linop), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine
         end interface
      
      !--------------------------------------
      !-----     RESOLVENT OPERATOR     -----
      !--------------------------------------
      
      ! --> Type.
         type, extends(abstract_linop_cdp), public :: resolvent_linop
            real(kind=dp) :: omega
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: matvec => resolvent_matvec
            procedure, pass(self), public :: rmatvec => resolvent_rmatvec
         end type
      
      ! --> Type-bound procedures: resolvent.f90
         interface
            module subroutine resolvent_matvec(self, vec_in, vec_out)
               class(resolvent_linop), intent(inout) :: self
               class(abstract_vector_cdp), intent(in) :: vec_in
               class(abstract_vector_cdp), intent(out) :: vec_out
            end subroutine
      
            module subroutine resolvent_rmatvec(self, vec_in, vec_out)
               class(resolvent_linop), intent(inout) :: self
               class(abstract_vector_cdp), intent(in) :: vec_in
               class(abstract_vector_cdp), intent(out) :: vec_out
            end subroutine
         end interface
      
      contains
      
         subroutine apply_exptA(vec_out, A, vec_in, tau, info, trans)
      !! Subroutine for the exponential propagator that conforms with the abstract interface
      !! defined in expmlib.f90
            class(abstract_vector_rdp), intent(out) :: vec_out
      !! Output vector
            class(abstract_linop_rdp), intent(inout) :: A
      !! Linear operator
            class(abstract_vector_rdp), intent(in) :: vec_in
      !! Input vector.
            real(dp), intent(in) :: tau
      !! Integration horizon
            integer, intent(out) :: info
      !! Information flag
            logical, optional, intent(in) :: trans
            logical :: transpose
      !! Direct or Adjoint?
      
      ! optional argument
            transpose = optval(trans, .false.)
      
      ! time integrator
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  select type (A)
                  type is (exptA_linop)
                     ! set integration time
                     A%tau = tau
                     if (transpose) then
                        call A%rmatvec(vec_in, vec_out)
                     else
                        call A%matvec(vec_in, vec_out)
                     end if
                  class default
                     call nek_stop_error("The intent [INOUT] argument 'A' must be of type 'exptA_linop', "//
     & "'exptA_linop_proj', or 'exptA_linop_frc'", this_module, 'apply_exptA')
                  end select
               class default
                  call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'", 
     & this_module, 'apply_exptA')
               end select
            class default
               call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'", 
     & this_module, 'apply_exptA')
            end select
         end subroutine apply_exptA
      
         subroutine compute_LNS_conv(c_x, c_y, c_z, uxp, uyp, uzp, trans)
      ! Compute the linearized convective terms
      ! We assume that the baseflow is set in vx, vy, vz
            real(dp), dimension(lv, 1), intent(out) :: c_x
            real(dp), dimension(lv, 1), intent(out) :: c_y
            real(dp), dimension(lv, 1), intent(out) :: c_z
            real(dp), dimension(lv, 1), intent(in) :: uxp
            real(dp), dimension(lv, 1), intent(in) :: uyp
            real(dp), dimension(lv, 1), intent(in) :: uzp
            logical, optional :: trans
      ! internals
            integer :: i
            logical :: transpose
      ! scratch arrays
            real(dp), dimension(lv) :: ta1, ta2, ta3, tb1, tb2, tb3
            transpose = optval(trans, .false.)
      ! (u.grad) Ub
            call opcopy(tb1, tb2, tb3, vx, vy, vz)          ! Save velocity
            call opcopy(vx, vy, vz, uxp, uyp, uzp)          ! U <-- u
            if (transpose) then
               call convop_adj(ta1, ta2, ta3, tb1, tb2, tb3, vx, vy, vz)
            else
               call convop(ta1, tb1)
               call convop(ta2, tb2)
               if (if3d) then
                  call convop(ta3, tb3)
               else
                  call rzero(ta3, lv)
               end if
            end if
            call opcopy(c_x, c_y, c_z, ta1, ta2, ta3) ! copy to output
      ! (Ub.grad) u
            call opcopy(vx, vy, vz, tb1, tb2, tb3)          ! Restore velocity
            if (transpose) then
               call convop_adj(ta1, ta2, ta3, tb1, tb2, tb3, vx, vy, vz)
            else
               call convop(tb1, uxp)
               call convop(tb2, uyp)
               if (if3d) then
                  call convop(tb3, uzp)
               else
                  call rzero(tb3, lv)
               end if
            end if
            call opadd2(c_x, c_y, c_z, tb1, tb2, tb3) ! add to output
            return
         end subroutine compute_LNS_conv
      
         subroutine compute_LNS_laplacian(d_x, d_y, d_z, uxp, uyp, uzp)
      ! compute the laplacian of the input field
            real(dp), dimension(lv, 1), intent(out) :: d_x
            real(dp), dimension(lv, 1), intent(out) :: d_y
            real(dp), dimension(lv, 1), intent(out) :: d_z
            real(dp), dimension(lv, 1), intent(in) :: uxp
            real(dp), dimension(lv, 1), intent(in) :: uyp
            real(dp), dimension(lv, 1), intent(in) :: uzp
            call lap_1D(d_x, uxp)
            call lap_1D(d_y, uyp)
            if (if3d) call lap_1D(d_z, uzp)
      ! multiply by 1/Re
            call col2(d_x, vdiff, lv)
            call col2(d_y, vdiff, lv)
            if (if3d) call col2(d_z, vdiff, lv)
            return
         end subroutine compute_LNS_laplacian
      
         subroutine lap_1D(nabla_u, u)
            real(dp), dimension(lv, 1), intent(in) :: u       ! perturbation velocity component
            real(dp), dimension(lxyz, lelv), intent(out) :: nabla_u
      ! internals
            real(dp), dimension(lxyz, lelv) :: ux, uy, uz
            real(dp), dimension(lxyz) :: us, ur, ut
            integer e, i, nel
            nel = lx1 - 1
            call gradm1(ux, uy, uz, u)
            do e = 1, lelv
            if (if3d) then
               call local_grad3(ur, us, ut, ux, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  nabla_u(i, e) = jacmi(i, e)*(ur(i)*rxm1(i, 1, 1, e) + us(i)*sxm1(i, 1, 1, e) + ut(i)*txm1(i, 1, 1, e))
               end do
               call local_grad3(ur, us, ut, uy, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  nabla_u(i, e) = nabla_u(i, e) + jacmi(i, e)*(ur(i)*rym1(i, 1, 1, e) + us(i)*sym1(i, 1, 1, e) + ut(i)*tym1(i, 1, 1, e))
               end do
               call local_grad3(ur, us, ut, uz, nel, e, dxm1, dxtm1)
               do i = 1, lxyz
                  nabla_u(i, e) = nabla_u(i, e) + jacmi(i, e)*(ur(i)*rzm1(i, 1, 1, e) + us(i)*szm1(i, 1, 1, e) + ut(i)*tzm1(i, 1, 1, e))
               end do
            else ! 2D
               call local_grad2(ur, us, ux, nel, e, dxm1, dytm1)
               do i = 1, lxyz
                  nabla_u(i, e) = jacmi(i, e)*(ur(i)*rxm1(i, 1, 1, e) + us(i)*sxm1(i, 1, 1, e))
               end do
               call local_grad2(ur, us, uy, nel, e, dxm1, dytm1)
               do i = 1, lxyz
                  nabla_u(i, e) = nabla_u(i, e) + jacmi(i, e)*(ur(i)*rym1(i, 1, 1, e) + us(i)*sym1(i, 1, 1, e))
               end do
            end if ! if3d
            end do
            return
         end subroutine lap_1D
      
         subroutine compute_LNS_gradp(gp_x, gp_y, gp_z, pp)
      ! compute the laplacian of the input field
            real(dp), dimension(lv, 1), intent(out) :: gp_x
            real(dp), dimension(lv, 1), intent(out) :: gp_y
            real(dp), dimension(lv, 1), intent(out) :: gp_z
            real(dp), dimension(lp, 1), intent(in) :: pp
      ! internals
            real(dp), dimension(lv) :: ta1, ta2, wrk
      ! map the perturbation pressure to the velocity mesh
            call mappr(wrk, pp, ta1, ta2)
      ! compute the gradient on the velocity mesh direclty
            call gradm1(gp_x, gp_y, gp_z, wrk)
            return
         end subroutine compute_LNS_gradp
      
         subroutine apply_L(Lux, Luy, Luz, ux, uy, uz, pres, trans)
      !! Apply LNS operator including the subtraction of the pressure gradient
      !! (but without the projection onto the divergence-free space)
            real(dp), dimension(lv, 1), intent(out) :: Lux
            real(dp), dimension(lv, 1), intent(out) :: Luy
            real(dp), dimension(lv, 1), intent(out) :: Luz
            real(dp), dimension(lv, 1), intent(in) :: ux
            real(dp), dimension(lv, 1), intent(in) :: uy
            real(dp), dimension(lv, 1), intent(in) :: uz
            real(dp), dimension(lp, 1), intent(in) :: pres
            logical, optional, intent(in) :: trans
      !! adjoint?
      ! internal
            real(dp), dimension(lv) :: utmpx, utmpy, utmpz
      ! Apply the linear operator to the velocity components
            call apply_Lv(Lux, Luy, Luz, ux, uy, uz, trans)
      ! and subtract the pressure gradient term
            call logger%log_debug(' pressure gradient', module=this_module, procedure='compute_L')
            call compute_LNS_gradp(utmpx, utmpy, utmpz, pres)
            call opsub2(Lux, Luy, Luz, utmpx, utmpy, utmpz)
            return
         end subroutine apply_L
      
         subroutine apply_Lv(Lux, Luy, Luz, ux, uy, uz, trans)
      !! Apply the convective and diffusive terms of the LNS operator to the velocity perturbation
            real(dp), dimension(lv, 1), intent(out) :: Lux
            real(dp), dimension(lv, 1), intent(out) :: Luy
            real(dp), dimension(lv, 1), intent(out) :: Luz
            real(dp), dimension(lv, 1), intent(in) :: ux
            real(dp), dimension(lv, 1), intent(in) :: uy
            real(dp), dimension(lv, 1), intent(in) :: uz
            logical, optional, intent(in) :: trans
      !! adjoint?
      ! internal
            real(dp), dimension(lv) :: utmpx, utmpy, utmpz
      ! apply BCs
            call bcdirvc(ux, uy, uz, v1mask, v2mask, v3mask)
      ! Diffusion term
            call logger%log_debug('diffusion term', module=this_module, procedure='compute_Lv')
            call compute_LNS_laplacian(Lux, Luy, Luz, ux, uy, uz)
      ! Convective terms
            call logger%log_debug('convective term', module=this_module, procedure='compute_Lv')
            call compute_LNS_conv(utmpx, utmpy, utmpz, ux, uy, uz, trans)
      ! subtract from output terms
            call opsub2(Lux, Luy, Luz, utmpx, utmpy, utmpz)
            return
         end subroutine apply_Lv

         subroutine get_restart(vec_in, istep)
            type(nek_dvector), intent(in) :: vec_in
            integer, intent(in) :: istep
            ! internal
            type(nek_dvector) :: vec_rst
            call vec_in%get_rst(vec_rst, istep)
            call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
         end subroutine get_restart

         subroutine save_restart(vec_out, nrst)
            type(nek_dvector), intent(out) :: vec_out
            integer, intent(in) :: nrst
            ! internal
            type(nek_dvector) :: vec_rst
            integer :: itmp
            real(dp) :: rtmp
            character(len=128) :: msg
            ! Copy the final solution to vector.
            call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            ! Fill up the restart fields
            write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
            call nek_log_information(msg, this_module, 'exptA_matvec')
            itmp = nsteps
            rtmp = time
            call setup_linear_solver(endtime = time + nrst*dt)
            nsteps = itmp
            do istep = nsteps + 1, nsteps + nrst
               call nek_advance()
               call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
               call vec_out%save_rst(vec_rst)
            end do
            ! Reset iteration count and time
            istep = itmp
            time  = rtmp
         end subroutine save_restart
      
      end module neklab_linops
