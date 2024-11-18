      module neklab_linops
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov_Logger
         use neklab_vectors
         use neklab_utils, only: nek2vec, vec2nek
         use neklab_nek_setup, only: setup_nonlinear_solver, setup_linear_solver
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      
         public :: apply_exptA
      
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
         
         ! --> Type-bound procedures.
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
      
      end module neklab_linops
