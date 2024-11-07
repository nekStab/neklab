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
         use neklab_utils, only: nek2vec, vec2nek, setup_nonlinear_solver, setup_linear_solver
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
      
      contains
      
      !-------------------------------------------------
      !-----     TYPE-BOUND PROCEDURES FOR LNS     -----
      !-------------------------------------------------
      
         subroutine init_LNS(self)
            class(LNS_linop), intent(in) :: self
      
      ! Force the baseflow field for dt/nsteps/clf computation
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
            call logger%log_message('Set self%baseflow -> vx, vy, vz, pr, t', module=this_module, procedure='init_LNS')
      
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
            integer, parameter :: ifield = 1 ! which material properties to use (vtrans)
            integer, parameter :: intype = 0 ! integration type (steady)
            integer, parameter :: imesh = 1 ! which mesh to use mesh
            integer :: iter
            real(dp), dimension(lx1*ly1*lz1*lelv, 1) :: resv1, resv2, resv3, dv1, dv2, dv3, h1, h2, h2inv
            common/scrns/resv1, resv2, resv3, dv1, dv2, dv3
            common/scrvh/h1, h2
            common/scrhi/h2inv
      
            call setup_linear_solver(transpose=.false., silent=.true.)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      !-------------------
      ! Apply LNS operator
      !-------------------
      
      ! apply BCs
            call bcdirvc(vxp, vyp, vzp, v1mask, v2mask, v3mask)
      
      ! compute pressure gradient
            call opgradt(resv1, resv2, resv3, prp)
      
      ! construct convective terms and add them to bf[xyz]p
            call advabp
      
      ! add the convective terms bf[xyz]p to the pressure gradient
            call opadd2(resv1, resv2, resv3, bfxp, bfyp, bfzp)
      
      ! set factors for the Helmholtz operator
            call rzero(h1, lv)
            call copy(h2, vtrans(1, 1, 1, 1, ifield), lv)
      ! and apply to the velocity field
            call ophx(dv1, dv2, dv3, vxp, vyp, vzp, h1, h2)
      
      ! substract result form rest of rhs and put into v[xyz]p
            call opsub3(vxp, vyp, vzp, resv1, resv2, resv3, dv1, dv2, dv3)
      
      ! project onto closest solenoidal space (~ incomprp)
      ! compute divergence of velocity dp = D^T @ u
            call opdiv(prp, vxp, vyp, vzp)
            call chsign(prp, lp)
            call ortho(prp) ! project out null-space (is done already in esolver)
      !! solve D @ D^T @ x = dp
            call rzero(h1, lv)
            call copy(h2, vtrans(1, 1, 1, 1, ifield), lv)
            iter = 200
            call hmh_flex_cg(prp, h1, h2, vmult, iter)
      !call hmh_gmres(prp, h1, h2, vmult, iter)
      !call rzero   (h1, lv)
      !call copy    (h2, vtrans(1,1,1,1,ifield), lv)
      !call invers2 (h2inv, h2, lv)
      !call esolver (prp, h1, h2, h2inv, intype)  ! --> overwrites prp with solution
      !! compute gradient of pressure correction --> velocity correction
            call opgradt(dv1, dv2, dv3, prp)
      !call opbinv  (resv1, resv2, resv3, dv1, dv2, dv3, h2inv)
      !! and add it to velocity
            call opadd2(vxp, vyp, vzp, dv1, dv2, dv3)
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
         end subroutine LNS_matvec
      
         subroutine adjLNS_matvec(self, vec_in, vec_out)
            class(LNS_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      !
      ! NOTE: only for a single perturbation vector !
      !
            integer, parameter :: ifield = 1 ! which material properties to use (vtrans)
            integer, parameter :: intype = 0 ! integration type (steady)
            integer, parameter :: imesh = 1 ! which mesh to use mesh
            real(dp), dimension(lx1, ly1, lz1, lelv) :: resv1, resv2, resv3, dv1, dv2, dv3, h1, h2, h2inv
            common/scrns/resv1, resv2, resv3, dv1, dv2, dv3
            common/scrvh/h1, h2
            common/scrhi/h2inv
      
            call setup_linear_solver(transpose=.true., silent=.true.)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      !---------------------------
      ! Apply adjoint LNS operator
      !---------------------------
      
      ! apply BCs
            call bcdirvc(vxp, vyp, vzp, v1mask, v2mask, v3mask)
      
      ! compute pressure gradient
            call opgradt(resv1, resv2, resv3, prp)
      
      ! construct convective terms and add them to bf[xyz]p
            call advabp_adjoint
      ! additional diagonal term for ON/on/O/o BCs
            call bc_out_adj(h2)
      
      ! add the convective terms bf[xyz]p to the pressure gradient
            call opadd2(resv1, resv2, resv3, bfxp, bfyp, bfzp)
      
      ! set factors for the Helmholtz operator
            call copy(h1, vtrans(1, 1, 1, 1, ifield), lv)
            call rzero(h2, lv)
      ! and apply to the velocity field
            call ophx(dv1, dv2, dv3, vxp, vyp, vzp, h1, h2)
      
      ! substract result form rest of rhs and put into v[xyz]p
            call opsub3(vxp, vyp, vzp, resv1, resv2, resv3, dv1, dv2, dv3)
      
      ! project onto closest solenoidal space (~ incomprp)
      
      ! compute divergence of velocity dp = D^T @ u
            call opdiv(prp, vxp, vyp, vzp)
            call chsign(prp, lp)
      ! call ortho(dp) ! project out null-space (is done already in esolver)
      ! solve D @ A^(-1) @ D^T @ x = dp using the uzawa splitting
      ! set factors for the Helmholtz operator
            call rzero(h1, lv)
            call copy(h2, vtrans(1, 1, 1, 1, ifield), lv)
            call invers2(h2inv, h2, lv)
            call esolver(prp, h1, h2, h2inv, intype)  ! --> overwrites prp with solution
      ! compute gradient of pressure correction --> velocity correction
            call opgradt(dv1, dv2, dv3, prp)
      ! and add it to velocity
            call opadd2(vxp, vyp, vzp, dv1, dv2, dv3)
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
         end subroutine adjLNS_matvec
      
      !---------------------------------------------------
      !-----     TYPE-BOUND PROCEDURES FOR exptA     -----
      !---------------------------------------------------
      
         subroutine init_exptA(self)
            class(exptA_linop), intent(in) :: self
      
      ! Force the baseflow field for dt/nsteps/clf computation
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
            call logger%log_message('Set self%baseflow -> vx, vy, vz, pr, t', module=this_module, procedure='init_exptA')
      
      ! Setup Nek5000 for perturbation solver
            call setup_linear_solver(solve_baseflow=.false., recompute_dt=.true., cfl_limit=0.5_dp)
      
            return
         end subroutine init_exptA
      
         subroutine exptA_matvec(self, vec_in, vec_out)
            class(exptA_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            call setup_linear_solver(transpose=.false., silent=.true.)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
            return
         end subroutine exptA_matvec
      
         subroutine exptA_rmatvec(self, vec_in, vec_out)
            class(exptA_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            call setup_linear_solver(transpose=.true., silent=.true.)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
            return
         end subroutine exptA_rmatvec
      
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
      
      end module neklab_linops
