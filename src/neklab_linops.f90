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
