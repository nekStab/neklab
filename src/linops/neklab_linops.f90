      module neklab_linops
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp, rtol_dp
         use LightKrylov, only: abstract_linop_rdp, abstract_sym_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: cg, cg_dp_opts, cg_dp_metadata
         use LightKrylov_Logger
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek, nek2pr_vec, pr_vec2nek
         use neklab_utils, only: setup_nonlinear_solver, setup_linear_solver
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
      
         integer, parameter :: lv    = lx1*ly1*lz1*lelv
         integer, parameter :: lp    = lx2*ly2*lz2*lelv
         integer, parameter :: lvxyz = lx1*ly1*lz1
      
         public :: apply_exptA
         public :: compute_LNS_conv
         public :: compute_LNS_gradp
         public :: compute_LNS_laplacian
         public :: apply_L
         public :: project_perturbation
      
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

         !-----------------------------------------------------
         !-----     LINEARIZED NAVIER-STOKES OPERATOR     -----
         !-----------------------------------------------------
      
         ! --> Type.
         type, extends(abstract_linop_rdp), public :: LNS_linop
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: init => init_LNS
            procedure, pass(self), public :: matvec => LNS_matvec
            procedure, pass(self), public :: rmatvec => LNS_rmatvec
         end type LNS_linop

         ! --> Type-bound procedures: LNS_operator.f90
         interface
            module subroutine init_LNS(self)
               class(LNS_linop), intent(in) :: self
            end subroutine
      
            module subroutine LNS_matvec(self, vec_in, vec_out)
               class(LNS_linop), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine

            module subroutine LNS_rmatvec(self, vec_in, vec_out)
               class(LNS_linop), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine
         end interface

         !---------------------------
         !-----     D^T @ D     -----
         !---------------------------
      
         ! --> Type.
         type, extends(abstract_sym_linop_rdp), public :: DTD_linop
         contains
            private
            procedure, pass(self), public :: matvec  => apply_DTD
            procedure, pass(self), public :: rmatvec => apply_DTD
         end type DTD_linop

         ! --> Type-bound procedures: pressure_projection.f90
         interface
            module subroutine apply_DTD(self, vec_in, vec_out)
               class(DTD_linop), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine
         end interface

         !--------------------------------------------------
         !-----     Projection onto div-free space     -----
         !--------------------------------------------------
      
         ! --> Type.
         type, extends(abstract_sym_linop_rdp), public :: P_div0_linop
            type(DTD_linop) :: DTD
         contains
            private
            procedure, pass(self), public :: matvec  => project_div0
            procedure, pass(self), public :: rmatvec => project_div0
         end type P_div0_linop

         ! --> Type-bound procedures: pressure_projection.f90
         interface
            module subroutine project_div0(self, vec_in, vec_out)
               class(P_div0_linop), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
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
                  end select
               end select
            end select
            return
         end subroutine apply_exptA

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
            ifield = 1
            call laplacian(d_x,uxp)
            call laplacian(d_y,uyp)
            if (if3d) call laplacian(d_z,uzp) 
            ! multiply by 1/Re
            call col2(d_x,vdiff(1,1,1,1,ifield),lv) 
            call col2(d_y,vdiff(1,1,1,1,ifield),lv) 
            if (if3d) call col2(d_z,vdiff(1,1,1,1,ifield),lv) 
            return
         end subroutine compute_LNS_laplacian

         subroutine laplacian(nabla_u,u)
   	      real(dp), dimension(lv,      1), intent(in)  :: u ! perturbation velocity component
   		   real(dp), dimension(lvxyz,lelv), intent(out) :: nabla_u
   		   ! internals
   		   real(dp), dimension(lv)    :: ux, uy, uz
   		   real(dp), dimension(lvxyz) :: ur, us, ut	      !	   
   		   integer :: e, i, nel
   		   nel = nx1-1
   		   call gradm1(ux,uy,uz,u)
   		   do e = 1, lelv
   		      if (if3d) then
   		   	   call local_grad3(ur,us,ut,ux,nel,e,dxm1,dxtm1)
   		   	   do i=1, lvxyz
   		   		   nabla_u(i,e) = jacmi(i,e)*(  ur(i)*rxm1(i,1,1,e)
     $	                                       + us(i)*sxm1(i,1,1,e)
     $	                                       + ut(i)*txm1(i,1,1,e) )
   		   	   enddo
   		   	   call local_grad3(ur,us,ut,uy,nel,e,dxm1,dxtm1)
   		   	   do i = 1, lvxyz
   		   	   	nabla_u(i,e) = nabla_u(i,e) + jacmi(i,e)*(  ur(i)*rym1(i,1,1,e)
     $	                                                      + us(i)*sym1(i,1,1,e)
     $	                                                      + ut(i)*tym1(i,1,1,e) )
   		   	   enddo
   		   	   call local_grad3(ur,us,ut,uz,nel,e,dxm1,dxtm1)
   		   	   do i = 1, lvxyz   
   		   	   	nabla_u(i,e) = nabla_u(i,e) + jacmi(i,e)*(  ur(i)*rzm1(i,1,1,e)
     $	                                                      + us(i)*szm1(i,1,1,e)
     $	                                                      + ut(i)*tzm1(i,1,1,e) )
   		   	   enddo
   		      else ! 2D
   		   	   call local_grad2(ur,us,ux,nel,e,dxm1,dytm1)
   		   	   do i = 1, lvxyz
   		   	      nabla_u(i,e) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $	                                     + us(i)*sxm1(i,1,1,e) )
   		   	   enddo
   		   	   call local_grad2(ur,us,uy,nel,e,dxm1,dytm1)
   		   	   do i = 1, lvxyz
   		   	   	nabla_u(i,e) = nabla_u(i,e) + jacmi(i,e)*(  ur(i)*rym1(i,1,1,e)
     $	                                                      + us(i)*sym1(i,1,1,e) )
   		   	   enddo
   		      endif ! if3d
   		   enddo
   		end subroutine laplacian

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
            ! apply BCs
            call bcdirvc(ux,uy,uz,v1mask,v2mask,v3mask)
            
            ! Diffusion term
            call logger%log_information('diffusion term', module=this_module, procedure='compute_LNS')
            call compute_LNS_laplacian(Lux,Luy,Luz,ux,uy,uz)

            ! Pressure gradient
            call logger%log_information('pressure gradient', module=this_module, procedure='compute_LNS')
            call compute_LNS_gradp(utmpx,utmpy,utmpz,pres)
            call opsub2(Lux,Luy,Luz,utmpx,utmpy,utmpz)

            ! Convective terms
            call logger%log_information('convective term', module=this_module, procedure='compute_LNS')
            call compute_LNS_conv(utmpx,utmpy,utmpz,ux,uy,uz,trans)
            call opsub2(Lux,Luy,Luz,utmpx,utmpy,utmpz)
            
            return
         end subroutine apply_L

         subroutine project_perturbation(dpr)
      !! Project perturbation velocity fields onto closest solenoidal space
      !! with div u = 0
      !! For this, we solve
      !!        D^T @ D @ dpr = D^T @ u
      !! using CG iteration to the recover the corresponding velocity correction
      !!        du = D @ p
      !! that we add to the previously computed velocity field.
            real(dp), dimension(lp, 1), intent(inout) :: dpr
      ! internal
            real(dp), dimension(lx1, ly1, lz1, lelv) :: vxtmp, vytmp, vztmp
            real(dp), dimension(lx1, ly1, lz1, lelv) :: h1, h2
            real(dp), dimension(lx2, ly2, lz2, lelv) :: wp, xCG, rpCG
            logical :: iconv
            integer :: iter
            real(dp) :: tolps, tolpss
            real(dp) :: rrp1, rrp2
            real(dp) :: alpha, beta
            real(dp) :: div0, rnorm, rnrm1, ratio, pDTDp
            real(dp), external :: glsc2 ! inner product
            common/scruz/wp, xCG, rpCG
            common/scrvh/h1, h2
      
            call ortho(dpr) ! project out null-space
      
            iter = 0
      
            call chktCG2(tolps, dpr, iconv) ! copute tolerance for pressure solver
            if (param(21) > 0 .and. tolps > abs(param(21))) tolps = abs(param(21))
      ! param(21) is the pressure tolerance
      
      !call uzprec(rpCG,dpr,h1,h2,0,wp)
      ! let's see whether we can use the preconditioner, this is what it does
      !call col3        (wp,dpr,h1m2,lp)     ! wp = dpr .* h1m1
      !call col3        (rpCG,wp,BM2inv,lp) ! rpCG = wp ./ bm2
      ! rpCG = dpr .* h1m2./bm2
      !rrp1 = glsc2 (rpCG,dpr,lp)                 ! < p , p >_(diag(h1m2/bm2)) = < r , z >
      !call copy    (dpr,rpCG,lp)
      
            rrp1 = glsc2(dpr, dpr, lp)
            call rzero(xCG, lp)                     ! xCG = 0
      
            beta = 0.0_dp
            div0 = 0.0_dp
      
            call convprn(iconv, rnorm, rrp1, dpr, rpCG, tolpss) ! computes rnorm, rrp1
      
            tolpss = tolps
            if (nid == 0) print *, 'rrp1', rrp1, tolpss
            cg_loop: do iter = 1, nmxp
               if (nid == 0) print *, 'Step', iter
      
               pDTDp = glsc2(dpr, wp, lp)        ! < p , p >_(D^T @ D)
      
               if (nid == 0) print *, 'pDTDp', pDTDp
      
               alpha = rrp1/pDTDp             ! step size
               call add2s2(xCG, dpr, alpha, lp)  ! x = x + alpha*dpr
               call add2s2(dpr, wp, -alpha, lp)  ! r = r - alpha*DTDdpr
      
               call ortho(dpr) ! project out null-space
      
      ! preconditioner
      !call uzprec  (rpCG,dpr,h1,h2,0,wp)
               call copy(rpCG, dpr)
      
      ! check for convergence
               call convprn(iconv, rnorm, rrp2, dpr, rpCG, tolpss) ! Convergence test for the pressure step;  < r , z >  ! rrp2 is output
               if (iter == 1) div0 = rnorm ! reference is the initial rnorm
               if (param(21) < 0) tolpss = abs(param(21))*div0 ! relative norm
               ratio = rnorm/div0
               if (nio == 0) write (6, '(I5,1P4E12.5," Divergence")') iter, tolpss, rnorm, div0, ratio
               if (iconv .and. iter > 1) exit cg_loop
      
               beta = rrp2/rrp1
               call add2s1(dpr, rpCG, beta, lp) ! dpr = beta*dpr + rpCG
      
               rrp1 = rrp2
            end do cg_loop
      
            iter = iter - 1
      
            if (iter > 0) call copy(dpr, xCG, lp) ! copy result into input/output vector
            call ortho(dpr)
      
            return
         end subroutine project_perturbation
      
      end module neklab_linops
