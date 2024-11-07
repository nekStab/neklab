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
      
         public :: apply_exptA
      
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
                  call apply_L(trans=.false.)
            
      ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
      
         end subroutine LNS_matvec
      
         subroutine adjLNS_matvec(self, vec_in, vec_out)
            class(LNS_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
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
                  call apply_L(trans=.true.)
      
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
      
         end subroutine adjLNS_matvec

         subroutine apply_L(trans)
      !! Apply LNS operator (before the projection onto the divergence-free space)
      !! This function assumes that the input vector has already been loaded into v[xyz]p
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
      ! set factors for the Helmholtz operator
            call rzero(h1, lv)
            call copy(h2, vtrans(1,1,1,1,ifield), lv)
      ! add a diagonal term to the operator for the ON/on/O/o boundary conditions    
            if (trans) call bc_out_adj(h2)
      ! and apply to the velocity field to compute the diffusion term
            call ophx(dv1, dv2, dv3, vxp, vyp, vzp, h1, h2)
      
      ! substract result form rest of rhs and put into v[xyz]p
            call opsub3(vxp, vyp, vzp, resv1, resv2, resv3, dv1, dv2, dv3)
            return
         end subroutine apply_L

         subroutine project_perturbation()
      !! Project perturbation velocity fields onto closest solenoidal space
      !! with div u = 0
      !! For this, we solve
      !!        D^T @ D @ dp = D^T @ u
      !! using CG iteration to the recover the corresponding velocity correction
      !!        du = D @ p
      !! that we add to the previously computed velocity field.

      ! compute divergence of velocity dp = D^T @ u
            call opdiv(prp, vxp, vyp, vzp)
      ! project out null-space
            call ortho(prp)

!            subroutine uzawa (rcg,h1,h2,h2inv,intype,iter)
!               C-----------------------------------------------------------------------
!               C
!               C     Solve the pressure equation by (nested) preconditioned 
!               C     conjugate gradient iteration.
!               C     INTYPE =  0  (steady)
!               C     INTYPE =  1  (explicit)
!               C     INTYPE = -1  (implicit)
!               C
!               C-----------------------------------------------------------------------
!                     include 'SIZE'
!                     include 'TOTAL'
!                     COMMON  /CTOLPR/ DIVEX
!                     COMMON  /CPRINT/ IFPRINT
!                     LOGICAL          IFPRINT
!                     REAL             RCG  (LX2,LY2,LZ2,LELV)
!                     REAL             H1   (LX1,LY1,LZ1,LELV)
!                     REAL             H2   (LX1,LY1,LZ1,LELV)
!                     REAL             H2INV(LX1,LY1,LZ1,LELV)
!                     COMMON /SCRUZ/   WP   (LX2,LY2,LZ2,LELV)
!                    $ ,               XCG  (LX2,LY2,LZ2,LELV)
!                    $ ,               PCG  (LX2,LY2,LZ2,LELV) 
!                    $ ,               RPCG (LX2,LY2,LZ2,LELV)
!                
!                     real*8 etime1,dnekclock
!                     integer*8 ntotg,nxyz2
!               
!               
!                     etime1 = dnekclock()
!                     DIVEX = 0.
!                     ITER  = 0
!               c
!                     CALL CHKTCG2 (TOLPS,RCG,ICONV)
!                     if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
!                    $   TOLPS = abs(param(21))
!               C
!               c      IF (ICONV.EQ.1) THEN
!               c         IF (NID.EQ.0) WRITE(6,9999) ITER,DIVEX,TOLPS
!               c         return
!               c      ENDIF
!               
!                     nxyz2 = lx2*ly2*lz2
!                     ntot2 = nxyz2*nelv
!                     ntotg = nxyz2*nelgv
!               
!                     CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
!                     RRP1 = GLSC2 (RPCG,RCG,NTOT2)
!                     CALL COPY    (PCG,RPCG,NTOT2)
!                     CALL RZERO   (XCG,NTOT2)
!                     if (rrp1.eq.0) return
!                     BETA = 0.
!                     div0=0.
!               C
!                     tolpss = tolps
!                     DO 1000 ITER=1,NMXP
!               C
!               C        CALL CONVPR  (RCG,tolpss,ICONV,RNORM)
!                        call convprn (iconv,rnorm,rrp1,rcg,rpcg,tolpss)
!               
!                        if (iter.eq.1)      div0   = rnorm
!                        if (param(21).lt.0) tolpss = abs(param(21))*div0
!               
!                        ratio = rnorm/div0
!                        IF (IFPRINT.AND.NIO.EQ.0) 
!                    $   WRITE (6,66) iter,tolpss,rnorm,div0,ratio,istep
!                  66    format(i5,1p4e12.5,i8,' Divergence')
!               c
!                        IF (ICONV.EQ.1.and.iter.gt.1) GOTO 9000
!               c        IF (ICONV.EQ.1.and.(iter.gt.1.or.istep.le.2)) GOTO 9000
!               c        IF (ICONV.EQ.1) GOTO 9000
!               c        if (ratio.le.1.e-5) goto 9000
!               
!               
!                        IF (ITER .NE. 1) THEN
!                           BETA = RRP1/RRP2
!                           CALL ADD2S1 (PCG,RPCG,BETA,NTOT2)
!                        ENDIF
!               
!                        CALL CDABDTP  (WP,PCG,H1,H2,H2INV,INTYPE)
!                        PAP   = GLSC2 (PCG,WP,NTOT2)
!               
!                        IF (PAP.NE.0.) THEN
!                           ALPHA = RRP1/PAP
!                        ELSE
!                           pcgmx = glamax(pcg,ntot2)
!                           wp_mx = glamax(wp ,ntot2)
!                           ntot1 = lx1*ly1*lz1*nelv
!                           h1_mx = glamax(h1 ,ntot1)
!                           h2_mx = glamax(h2 ,ntot1)
!                           if (nid.eq.0) write(6,*) 'ERROR: pap=0 in uzawa.'
!                    $      ,iter,pcgmx,wp_mx,h1_mx,h2_mx
!                           call exitt
!                        ENDIF
!                        CALL ADD2S2 (XCG,PCG,ALPHA,NTOT2)
!                        CALL ADD2S2 (RCG,WP,-ALPHA,NTOT2)
!               
!                        if (iter.eq.-1) then
!                           call convprn (iconv,rnrm1,rrpx,rcg,rpcg,tolpss)
!                           if (iconv.eq.1) then
!                              rnorm = rnrm1
!                              ratio = rnrm1/div0
!                              if (nio.eq.0) 
!                    $         write (6,66) iter,tolpss,rnrm1,div0,ratio,istep
!                              goto 9000
!                           endif
!                        endif
!               
!                        call ortho(rcg)
!               
!                        RRP2 = RRP1
!                        CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
!               c        RRP1 = GLSC2 (RPCG,RCG,NTOT2)
!               
!                1000 CONTINUE
!                     if (nid.eq.0) WRITE (6,3001) ITER,RNORM,tolpss
!               c     if (istep.gt.20) CALL EMERXIT
!                3001 FORMAT(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
!                9000 CONTINUE
!               
!                     divex = rnorm
!                     iter  = iter-1
!               
!                     if (iter.gt.0) call copy (rcg,xcg,ntot2)
!                     call ortho(rcg)
!               
!                     etime1 = dnekclock()-etime1
!                     IF (NIO.EQ.0) WRITE(6,9999) ISTEP, '  U-Press std. ',
!                    &                            ITER,DIVEX,div0,tolpss,etime1
!                9999 FORMAT(I11,a,I7,1p4E13.4)
!               19999 FORMAT(I11,'  U-Press 1.e-5: ',I7,1p4E13.4)
!               C
!               C
!                     return
!                     END

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

!         ! project onto closest solenoidal space (~ incomprp)
!      
!      ! compute divergence of velocity dp = D^T @ u
!         call opdiv(prp, vxp, vyp, vzp)
!         call chsign(prp, lp)
!! call ortho(dp) ! project out null-space (is done already in esolver)
!! solve D @ A^(-1) @ D^T @ x = dp using the uzawa splitting
!! set factors for the Helmholtz operator
!         call rzero(h1, lv)
!         call copy(h2, vtrans(1, 1, 1, 1, ifield), lv)
!         call invers2(h2inv, h2, lv)
!         call esolver(prp, h1, h2, h2inv, intype)  ! --> overwrites prp with solution
!! compute gradient of pressure correction --> velocity correction
!         call opgradt(dv1, dv2, dv3, prp)
!! and add it to velocity
!         call opadd2(vxp, vyp, vzp, dv1, dv2, dv3)
      
      end module neklab_linops
