      module neklab_plinops
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_optval, only: optval
      ! Default real kind.
         use LightKrylov, only: dp, rtol_dp, atol_dp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_sym_linop_rdp, abstract_vector_rdp
      ! Logging
         use LightKrylov_Logger
      ! Linear solvers
         use LightKrylov, only: cg, cg_dp_opts, cg_dp_metadata
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_pvectors
         use neklab_utils, only: nek2vec, vec2nek, nek2pvec, pvec2nek
         use neklab_utils, only: setup_nonlinear_solver, setup_linear_solver, nek_status
         use neklab_LNS_utils
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_plinops'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
      !---------------------------------------------------------
      !-----     Linearized NS Operator (velocity only)    -----
      !---------------------------------------------------------
      
         type, extends(abstract_sym_linop_rdp), public :: LNS_v_linop
         contains
            private
            procedure, pass(self), public :: matvec  => apply_L_direct
            procedure, pass(self), public :: rmatvec => apply_L_adjoint
         end type LNS_v_linop

      !---------------------------
      !-----     D^T @ D     -----
      !---------------------------
      
         type, extends(abstract_sym_linop_rdp), public :: DTD_linop
         contains
            private
            procedure, pass(self), public :: matvec  => apply_DTD
            procedure, pass(self), public :: rmatvec => apply_DTD
         end type DTD_linop

      !--------------------------------------------------
      !-----     Projection onto div-free space     -----
      !--------------------------------------------------
      
         type, extends(abstract_linop_rdp), public :: P_div0
            type(DTD_linop) :: DTD
         contains
            private
            procedure, pass(self), public :: matvec  => project_div0
            procedure, pass(self), public :: rmatvec => project_div0
         end type P_div0

      !---------------------------------------------------------------------
      !-----     LINEARIZED NAVIER-STOKES OPERATOR with projection     -----
      !---------------------------------------------------------------------
      
         type, extends(abstract_linop_rdp), public :: LNS_div0_linop
            type(nek_dvector) :: baseflow
            type(LNS_v_linop) :: L
            type(P_div0)      :: P
         contains
            private
            procedure, pass(self), public :: matvec  =>    LNS_matvec
            procedure, pass(self), public :: rmatvec => adjLNS_matvec
         end type LNS_div0_linop

      contains

      !-------------------------------------------------
      !-----     TYPE-BOUND PROCEDURES FOR DTD     -----
      !-------------------------------------------------
      
         subroutine apply_DTD(self, vec_in, vec_out)
            class(DTD_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            real(dp), dimension(lv) :: vxtmp, vytmp, vztmp
            ifield = 1
            select type (vec_in)
            type is (nek_pdvector)
               select type (vec_out)
               type is (nek_pdvector)
                  ! compute wp = D^T @ D @ dpr
                  !call compute_LNS_gradp(vxtmp, vytmp, vztmp, vec_in%pr)
                  call opgradt(vxtmp, vytmp, vztmp, vec_in%pr)
                  call bcdirvc(vxtmp, vytmp, vztmp, v1mask, v2mask, v3mask)
                  call opdiv(vec_out%pr, vxtmp, vytmp, vztmp)
               end select
            end select
            return
         end subroutine apply_DTD

      !----------------------------------------------------
      !-----     TYPE-BOUND PROCEDURES FOR P_div0     -----
      !----------------------------------------------------
      
         subroutine project_div0(self, vec_in, vec_out)
            class(P_div0), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            real(dp), dimension(lv) :: dv1, dv2, dv3
            type(nek_pdvector)   :: dpr, x
            type(cg_dp_opts)     :: opts
            type(cg_dp_metadata) :: meta
            integer :: info
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Copy into to nek
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Compute divergence of velocity dp = D^T @ u
                  call opdiv(prp, vxp, vyp, vzp)
      ! Copy from nek to nek_pdvector
                  call nek2pvec(dpr, prp)
      ! Solve linear system D^T @ D x = dpr
                  opts = cg_dp_opts(maxiter = 100, if_print_metadata = .true.)
                  meta = cg_dp_metadata()
                  call cg(self%DTD, dpr, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, meta=meta)
      ! Copy solution back to nek
                  call pvec2nek(prp, x)
      ! Compute velocity correction
                  call opgradt(dv1, dv2, dv3, prp)
      ! Compute output
                  call opsub2(vxp, vyp, vzp, dv1, dv2, dv3)
      ! Zero out pressure
                  call rzero(prp, lp)
      ! Copty output to vec_out
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
            return
         end subroutine project_div0

         subroutine apply_L_direct(self, vec_in, vec_out)
            class(LNS_v_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            real(dp), dimension(lv, 1) :: Lux, Luy, Luz
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Copy into to nek
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Apply LNS operator
                  call apply_L(Lux, Luy, Luz, vxp, vyp, vzp, prp, trans=.false.)
      ! Copty output to vec_out
                  call nek2vec(vec_out, Lux, Luy, Luz, prp, tp)
               end select
            end select
            return
         end subroutine apply_L_direct

         subroutine apply_L_adjoint(self, vec_in, vec_out)
            class(LNS_v_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            real(dp), dimension(lv, 1) :: Lux, Luy, Luz
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Copy into to nek
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Apply adjoint LNS operator
                  call apply_L(Lux, Luy, Luz, vxp, vyp, vzp, prp, trans=.true.)
      ! Copty output to vec_out
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               end select
            end select
            return
         end subroutine apply_L_adjoint

         subroutine LNS_matvec(self, vec_in, vec_out)
            class(LNS_div0_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            type(nek_dvector) :: v0, v1
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Force the baseflow field.
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Load input vector into temporary vector
                  call copy(v0, vec_in)
      ! Apply LNS operator in place
                  call self%L%matvec(v0, v1)
      ! Project onto div-free space
                  call self%P%matvec(v1, v0)
      ! Copy to output
                  call copy(vec_out, v0)
               end select
            end select
            return
         end subroutine LNS_matvec
            
         subroutine adjLNS_matvec(self, vec_in, vec_out)
            class(LNS_div0_linop), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
            ! internal
            type(nek_dvector) :: v0, v1
            ifadj = .true.
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Force the baseflow field.
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Load input vector into temporary vector
                  call copy(v0, vec_in)
      ! Apply LNS operator in place
                  call self%L%rmatvec(v0, v1)
      ! Project onto div-free space
                  call self%P%matvec(v1, v0)
      ! Copy to output
                  call copy(vec_out, v0)
               end select
            end select
            return
         end subroutine adjLNS_matvec
   
      end module neklab_plinops