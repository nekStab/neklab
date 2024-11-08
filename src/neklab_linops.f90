      module neklab_linops
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_optval, only: optval
      ! Default real kind.
         use LightKrylov, only: dp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek, setup_nonlinear_solver, setup_linear_solver, nek_status
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
         public :: apply_exptA, apply_L
      
      !------------------------------------------
      !-----     EXPONENTIAL PROPAGATOR     -----
      !------------------------------------------
      
         type, extends(abstract_linop_rdp), public :: exptA_linop
            real(kind=dp) :: tau
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: init => init_exptA
            procedure, pass(self), public :: matvec => exptA_matvec
            procedure, pass(self), public :: rmatvec => exptA_rmatvec
         end type exptA_linop

      !-----------------------------------------------------
      !-----     LINEARIZED NAVIER-STOKES OPERATOR     -----
      !-----------------------------------------------------
      
         type, extends(abstract_linop_rdp), public :: LNS_linop
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: init => init_LNS
            procedure, pass(self), public :: matvec => LNS_matvec
            procedure, pass(self), public :: rmatvec => adjLNS_matvec
         end type LNS_linop
      
      contains
      
      !---------------------------------------------------
      !-----     TYPE-BOUND PROCEDURES FOR exptA     -----
      !---------------------------------------------------
      
         subroutine init_exptA(self)
            class(exptA_linop), intent(in) :: self
      
      ! Force the baseflow field for dt/nsteps/clf computation
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
            !call logger%log_message('Set self%baseflow -> vx, vy, vz, pr, t', module=this_module, procedure='init_exptA')
      
      ! Setup Nek5000 for perturbation solver
            call setup_linear_solver(solve_baseflow=.false., recompute_dt=.true., 
     $                     endtime = self%tau, cfl_limit=0.5_dp)
      
            return
         end subroutine init_exptA
      
         subroutine exptA_matvec(self, vec_in, vec_out)
            class(exptA_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Force the baseflow field.
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)

      ! Ensure the appropriate solver state
                  call setup_linear_solver(transpose=.false., silent=.true.)
      
      ! Sets the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
                 
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
   
      ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
      
            return
         end subroutine exptA_matvec
      
         subroutine exptA_rmatvec(self, vec_in, vec_out)
            class(exptA_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Force the baseflow field.
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)

      ! Ensure the appropriate solver state
                  call setup_linear_solver(transpose=.true., silent=.true.)
      
      ! Sets the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
                 
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
   
      ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
      
            return
         end subroutine exptA_rmatvec

      !-------------------------------------------------
      !-----     TYPE-BOUND PROCEDURES FOR LNS     -----
      !-------------------------------------------------
      
         subroutine init_LNS(self)
            class(LNS_linop), intent(in) :: self
      
      ! Force the baseflow field for dt/nsteps/clf computation
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
            !call logger%log_message('Set self%baseflow -> vx, vy, vz, pr, t', module=this_module, procedure='init_LNS')
      
      ! Setup Nek5000 for perturbation solver
            call setup_linear_solver(solve_baseflow=.false., recompute_dt=.true., cfl_limit=0.5_dp)
      
            return
         end subroutine init_LNS
      
         subroutine LNS_matvec(self, vec_in, vec_out)
            class(LNS_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            real(dp), dimension(lv, 1) :: dv1, dv2, dv3
            common /scrns/ dv1, dv2, dv3
      !
      ! NOTE: only for a single perturbation vector !
      !
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Force the baseflow field.
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)
            
      ! Ensure correct settings for Nek
                  call setup_linear_solver(transpose=.false., silent=.true.)

      ! Sets the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            
      ! Apply LNS operator in place
                  call apply_L(vxp, vyp, vzp, trans=.false.)

      ! Compute divergence of velocity dp = D^T @ u
                  call opdiv(prp, vxp, vyp, vzp)

                  !call outpost(vxp, vyp, vzp, prp, tp, 'tst')

      ! Solve for pressure correction to enforce continuity D^T @ D @ dp = D^T @ u
                  call project_perturbation(prp) ! in-place
      
      ! Compute corresponding velocity correction
                  call opgradt(dv1, dv2, dv3, prp)

      ! Compute final velocity
                  call opadd2(vxp, vyp, vzp, dv1, dv2, dv3)
            
      ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
      
         end subroutine LNS_matvec
      
         subroutine adjLNS_matvec(self, vec_in, vec_out)
            class(LNS_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            real(dp), dimension(lv, 1) :: dv1, dv2, dv3
            common /scrns/ dv1, dv2, dv3
      !
      ! NOTE: only for a single perturbation vector !
      !
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Force the baseflow field.
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)
               
      ! Ensure the appropriate solver state
                  call setup_linear_solver(transpose=.true., silent=.true.)
            
      ! Sets the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            
      ! Apply adjoint LNS operator in place
                  call apply_L(vxp, vyp, vzp, trans=.true.)

      ! Compute divergence of velocity dp = D^T @ u
                  call opdiv(prp, vxp, vyp, vzp)

                  !call outpost(vxp, vyp, vzp, prp, tp, 'tst')

      ! Solve for pressure correction to enforce continuity D^T @ D @ dp = D^T @ u
                  call project_perturbation(prp) ! in-place
      
      ! Compute corresponding velocity correction
                  call opgradt(dv1, dv2, dv3, prp)

      ! Compute final velocity
                  call opadd2(vxp, vyp, vzp, dv1, dv2, dv3)
      
      ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
      
         end subroutine adjLNS_matvec

         subroutine apply_L(vxp_in, vyp_in, vzp_in, trans)
      !! Apply LNS operator (before the projection onto the divergence-free space)
      !! This function assumes that the input vector has already been loaded into v[xyz]p
            real(dp), dimension(lv, 1), intent(inout) :: vxp_in
            real(dp), dimension(lv, 1), intent(inout) :: vyp_in
            real(dp), dimension(lv, 1), intent(inout) :: vzp_in
            logical, intent(in) :: trans
            !! adjoint?
            ! internal
            real(dp), dimension(lv, 1) :: resv1, resv2, resv3, dv1, dv2, dv3, h1, h2
            common /scrns/ resv1, resv2, resv3, dv1, dv2, dv3
            common /scrvh/ h1, h2

            ifield = 1
      ! apply BCs
            call bcdirvc(vxp, vyp, vzp, v1mask, v2mask, v3mask)
      
      !---------------------
      ! Pressure gradient
      !---------------------
            call opgradt(resv1, resv2, resv3, prp)
            
      !---------------------
      ! Convective term
      !---------------------
            if (trans) then
      ! construct convective terms and add them to bf[xyz]p
               call advabp_adjoint
            else  
      ! construct convective terms and add them to bf[xyz]p
               call advabp
            end if
      ! add the convective terms bf[xyz]p to the pressure gradient
            call opadd2(resv1, resv2, resv3, bfxp, bfyp, bfzp)
      
      !---------------------
      ! Diffusion term
      !---------------------
      ! set factors for the Helmholtz operator  (H = h1*A + h2*B)
            call copy(h1, vdiff(1,1,1,1,ifield), lv)
            call rzero(h2, lv)
      ! add a diagonal term to the operator for the ON/on/O/o boundary conditions    
      !      if (trans) call bc_out_adj(h2) ! not needed since h2 = 0
      ! and apply to the velocity field to compute the diffusion term
            call ophx(dv1, dv2, dv3, vxp, vyp, vzp, h1, h2)
      
      ! substract result form rest of rhs and put into v[xyz]p
            call opsub3(vxp, vyp, vzp, resv1, resv2, resv3, dv1, dv2, dv3)
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
            common /scruz/ wp, xCG, rpCG
            common /scrvh/ h1, h2

            call ortho(dpr) ! project out null-space

            iter = 0

            CALL chktCG2 (tolps,dpr,iconv) ! copute tolerance for pressure solver
            if (param(21) > 0 .and. tolps > abs(param(21))) tolps = abs(param(21))
            ! param(21) is the pressure tolerance

            !call uzprec(rpCG,dpr,h1,h2,0,wp)
            ! let's see whether we can use the preconditioner, this is what it does
            !call col3        (wp,dpr,h1m2,lp)     ! wp = dpr .* h1m1
            !call col3        (rpCG,wp,BM2inv,lp) ! rpCG = wp ./ bm2
            ! rpCG = dpr .* h1m2./bm2
            !rrp1 = glsc2 (rpCG,dpr,lp)                 ! < p , p >_(diag(h1m2/bm2)) = < r , z >
            !call copy    (dpr,rpCG,lp)

            rrp1 = glsc2(dpr,dpr,lp)
            call rzero   (xCG,lp)                     ! xCG = 0

            beta = 0.0_dp
            div0 = 0.0_dp

            call convprn (iconv,rnorm,rrp1,dpr,rpCG,tolpss) ! computes rnorm, rrp1

            tolpss = tolps
            if (nid == 0) print *, 'rrp1', rrp1, tolpss
            cg_loop: do iter = 1, nmxp
               if (nid == 0) print *, 'Step', iter
               
               pDTDp = glsc2(dpr,wp,lp)        ! < p , p >_(D^T @ D)

               if (nid == 0) print *, 'pDTDp', pDTDp
      
               alpha = rrp1/pDTDp             ! step size
               call add2s2 (xCG,dpr,alpha,lp)  ! x = x + alpha*dpr
               call add2s2 (dpr,wp,-alpha,lp)  ! r = r - alpha*DTDdpr
            
               call ortho(dpr) ! project out null-space
      
               ! preconditioner
               !call uzprec  (rpCG,dpr,h1,h2,0,wp)
               call copy(rpCG, dpr)

               ! check for convergence  
               call convprn (iconv,rnorm,rrp2,dpr,rpCG,tolpss) ! Convergence test for the pressure step;  < r , z >  ! rrp2 is output
               if (iter == 1) div0 = rnorm ! reference is the initial rnorm
               if (param(21) < 0) tolpss = abs(param(21))*div0 ! relative norm
               ratio = rnorm/div0
               if (nio == 0) write(6,'(I5,1P4E12.5," Divergence")') iter,tolpss,rnorm,div0,ratio
               if (iconv .and. iter > 1) exit cg_loop

               beta = rrp2/rrp1
               call add2s1 (dpr,rpCG,beta,lp) ! dpr = beta*dpr + rpCG

               rrp1 = rrp2
            end do cg_loop

            iter  = iter-1

            if (iter > 0) call copy (dpr,xCG,lp) ! copy result into input/output vector
            call ortho(dpr)
               
            return
         end subroutine project_perturbation
      
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
                     time = 0.0_dp
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

!         real(dp), dimension(lp, 1), intent(inout) :: dp
!            ! internal
!            real(dp), dimension(lx1, ly1, lz1, lelv) :: vxtmp, vytmp, vztmp
!            real(dp), dimension(lx1, ly1, lz1, lelv) :: h1 !, h2, h2inv
!            real(dp), dimension(lx2, ly2, lz2, lelv) :: wp, xCG, rpCG
!            logical :: iconv
!            integer :: iter
!            real(dp) :: tolps, tolpss
!            real(dp) :: rrp1, rrp2 
!            real(dp) :: alpha, beta
!            real(dp) :: div0, rnorm, rnrm1, ratio, pDTDp
!            real(dp), external :: glsc2 ! inner product
!            common /scruz/ wp, xCG, pCG, rpCG
!            common /scrvh/ h1 !, h2, h2inv
!            COMMON /CTOLPR/ DIVEX
!
!            call ortho(dp) ! project out null-space
!
!            iter = 0
!            divex = 0.0_dp
!
!            CALL chktCG2 (tolps,dp,iconv) ! copute tolerance for pressure solver
!            if (param(21) > 0 .and. tolps > abs(param(21))) tolps = abs(param(21))
!            ! param(21) is the pressure tolerance
!
!            call uzprec(rpCG,dp,h1,h2,0,wp)
!            ! let's see whether we can use the preconditioner, this is what it does
!            !call col3        (wp,dp,h1m2,lp)     ! wp = dp .* h1m1
!            !call col3        (rpCG,wp,BM2inv,lp) ! rpCG = wp ./ bm2
!            ! rpCG = dp .* h1m2./bm2
!            rrp1 = glsc2 (rpCG,dp,lp)                 ! < p , p >_(diag(h1m2/bm2)) = < r , z >
!            call copy    (dp,rpCG,lp)
!            call rzero   (xCG,lp)                     ! xCG = 0
!
!            beta = 0.0_dp
!            div0 = 0.0_dp
!
!            tolpss = tolps
!            cg_iter: do iter = 1, nmxp
!
!               call convpr  (dp,tolpss,iconv,rnorm)           ! convergence test for the pressure step    
!               call convprn (iconv,rnorm,rrp1,dp,rpCG,tolpss) ! Convergence test for the pressure step;  < r , z >
!
!               if (iter == 1) div0 = rnorm
!
!               if (param(21) < 0) tolpss = abs(param(21))*div0 ! relative norm
!
!               ratio = rnorm/div0
!
!               if (nio == 0) then 
!                  write(6,'(I5,1P4E12.5," Divergence")') iter,tolpss,rnorm,div0,ratio
!               end if
!
!               if (iconv .and. iter > 1) exit cg_iter
!
!               if (iter /= 1) then
!                  beta = rrp1/rrp2
!                  call add2s1 (dp,rpCG,beta,lp) ! dp = beta*dp + rpCG
!               endif
!
!               ! compute wp = D^T @ D @ dp
!               call opgradt(vxtmp, vytmp, vztmp, dp)
!               call opdiv(wp, vxtmp, vytmp, vztmp)
!
!               pDTDp = glsc2(dp,wp,lp)        ! < p , p >_(D^T @ D)
!      
!               alpha = rrp1/pDTDp             ! step size
!               call add2s2 (xCG,dp,alpha,lp)  ! x = x + alpha*dp
!               call add2s2 (dp,wp,-alpha,lp)  ! r = r - alpha*DTDdp
!      
!               if (iter == -1) then
!                  call convprn (iconv,rnrm1,rrpx,dp,rpCG,tolpss)
!                  if (iconv) then
!                     rnorm = rnrm1
!                     ratio = rnrm1/div0
!                     if (nio == 0) write(6,'(I5,1P4E12.5," Divergence")') iter,tolpss,rnorm,div0,ratio
!                  endif
!               endif
!      
!               call ortho(dp) ! project out null space
!      
!               rrp2 = rrp1
!               call uzprec  (rpCG,dp,h1,h2,0,wp)
!            end do
!
!            divex = rnorm
!            iter  = iter-1
!
!            if (iter > 0) call copy (dp,xCG,lp) ! copy result into input/output vector
!            call ortho(dp)
      
      end module neklab_linops
