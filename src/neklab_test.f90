      module neklab_test
    !! Unit tests for neklab vector types and linear operators.
    !!
    !! This module is compiled into the neklab library ONLY when the build is
    !! configured with the NEKLAB_TESTS switch (see makeneklab). Production
    !! builds neither compile nor link it.
    !!
    !! Two public entry points are provided:
    !!
    !!   * test_vector(x, name, nfail)
    !!         Runs LightKrylov's verify_vector_axioms on a neklab vector type.
    !!
    !!   * test_operator(A, mold, name, nfail [, test_adjoint])   [generic]
    !!         Runs LightKrylov's verify_linop_axioms on a neklab linear
    !!         operator. Resolves to the real (rdp) or complex (cdp) body by the
    !!         type of A. This checks linearity (additivity + homogeneity) of
    !!         matvec/rmatvec and, unless test_adjoint=.false., also verifies
    !!         adjoint consistency  <A u, v> = <u, A^H v>.
    !!
    !! The test case / baseflow is NOT passed in: each entry point is called
    !! from the corresponding case .usr file, so the case is already obvious
    !! from context.
    !!
    !! Adding coverage for a new type is one line at the relevant call site in
    !! the case .usr (istep==0 block): construct the vector/operator, then call
    !! test_vector / test_operator.
    !!
    !! A non-zero nfail, and any 'NEKLAB_TEST_FAIL' sentinel line in the log,
    !! signals failure. Each case .usr should `error stop 1` when nfail > 0 so
    !! the process exit code also reflects the result.

       use LightKrylov, only: dp
       use LightKrylov, only: abstract_vector_rdp, abstract_vector_cdp
       use LightKrylov, only: abstract_linop_rdp, abstract_linop_cdp
       use LightKrylov, only: verify_vector_axioms
       use LightKrylov, only: verify_linop_axioms
       use LightKrylov_Logger

       use neklab_vectors, only: nek_dvector, nek_ext_dvector, nek_zvector

       implicit none
       private

       character(len=*), parameter :: this_module = 'neklab_test'

    !> Tolerance for the operator axiom checks. Loose enough to absorb the
    !> roundoff of independent time integrations (matvec is a propagation, so
    !> A(a*x+b*y) and a*A(x)+b*A(y) differ at O(eps*||A||*nsteps)), tight
    !> enough that a genuine linearity/adjoint bug (stray nonlinear term, sign
    !> error, forgotten baseflow) fails by orders of magnitude.
       real(dp), parameter, public :: linop_rtol = 1.0e-9_dp

    !> Tolerance passed to verify_vector_axioms. The neklab dot product runs a
    !> global reduction (gop), so associativity across ranks does not hold at
    !> machine precision; use a loosened value rather than the default.
       real(dp), parameter, public :: vector_axioms_tol = 1.0e-11_dp

       public :: test_vector
       public :: test_operator

       interface test_vector
          module procedure test_vector_rdp
          module procedure test_vector_cdp
       end interface test_vector

       interface test_operator
          module procedure test_operator_rdp
          module procedure test_operator_cdp
       end interface test_operator

      contains

    !---------------------------------------------------------------------
    !-----                      VECTOR AXIOMS                        -----
    !---------------------------------------------------------------------

       subroutine test_vector_rdp(x, name, ntrials)
    !! Run the LightKrylov vector-space axiom suite on a real neklab vector.
          class(abstract_vector_rdp), intent(in) :: x
    !! A (zeroed) instance of the vector type to test; used as a mold from
    !! which the axiom suite draws random vectors internally.
          character(len=*), intent(in) :: name
    !! Human-readable name of the type, for logging.
          integer, optional, intent(in) :: ntrials
    !! Number of trials to run; default is 100.

          logical :: success

          success = verify_vector_axioms(x, ntrials=ntrials, tolerance=vector_axioms_tol)
          call report(name, 'vector axioms', success)

       end subroutine test_vector_rdp

       subroutine test_vector_cdp(x, name, ntrials)
    !! Run the LightKrylov vector-space axiom suite on a complex neklab vector.
          class(abstract_vector_cdp), intent(in) :: x
          character(len=*), intent(in) :: name
          integer, optional, intent(in) :: ntrials

          logical :: success

          success = verify_vector_axioms(x, ntrials=ntrials, tolerance=vector_axioms_tol)
          call report(name, 'vector axioms', success)

       end subroutine test_vector_cdp

    !---------------------------------------------------------------------
    !-----             OPERATOR AXIOMS  (REAL, rdp)                 -----
    !---------------------------------------------------------------------

       subroutine test_operator_rdp(A, mold, name, ntrials, test_zero_map, test_adjoint)
    !! Verify linearity (and, by default, adjoint consistency) of a real
    !! linear operator via LightKrylov's verify_linop_axioms.
          class(abstract_linop_rdp), intent(inout) :: A
    !! The operator to test. Must already be constructed and initialised,
    !! with its baseflow set, by the caller.
          class(abstract_vector_rdp), intent(in) :: mold
    !! An instance of the vector type the operator acts on (e.g. a nek_dvector).
    !! Used only as a mold to generate random inputs; its contents are
    !! irrelevant and are overwritten.
          character(len=*), intent(in) :: name
    !! Human-readable name of the operator type, for logging.
          integer, optional, intent(in) :: ntrials
    !! How many trials to run?
          logical, optional, intent(in) :: test_zero_map
    !! Also verify explicitly that A(0) = 0? (Default: .true.)
          logical, optional, intent(in) :: test_adjoint
    !! Also verify adjoint consistency <A u, v> = <u, A^H v>? (Default: .true.)

          logical :: success

          success = verify_linop_axioms(A, mold, ntrials=ntrials, tolerance=linop_rtol, test_zero_map=test_zero_map, test_adjoint=test_adjoint)
          call report(name, 'linop axioms', success)

       end subroutine test_operator_rdp

    !---------------------------------------------------------------------
    !-----           OPERATOR AXIOMS  (COMPLEX, cdp)               -----
    !---------------------------------------------------------------------

       subroutine test_operator_cdp(A, mold, name, ntrials, test_zero_map, test_adjoint)
    !! Verify linearity (and, by default, adjoint consistency) of a complex
    !! linear operator via LightKrylov's verify_linop_axioms.
          class(abstract_linop_cdp), intent(inout) :: A
          class(abstract_vector_cdp), intent(in) :: mold
    !! An instance of the vector type the operator acts on (e.g. a nek_zvector).
          character(len=*), intent(in) :: name
          integer, optional, intent(in) :: ntrials
          logical, optional, intent(in) :: test_zero_map    
          logical, optional, intent(in) :: test_adjoint

          logical :: success

          success = verify_linop_axioms(A, mold, ntrials=ntrials, tolerance=linop_rtol, test_zero_map=test_zero_map, test_adjoint=test_adjoint)
          call report(name, 'linop axioms', success)

       end subroutine test_operator_cdp

    !---------------------------------------------------------------------
    !-----                      LOG HELPER                          -----
    !---------------------------------------------------------------------

       subroutine report(name, test_name, success)
          character(len=*), intent(in) :: name
    !! Human-readable name of the tested type/operator.
          character(len=*), intent(in) :: test_name
    !! What was tested, e.g. 'vector axioms' or 'linop axioms'.
          logical, intent(in) :: success
          character(len=256) :: msg

          character(len=*), parameter :: this_procedure = 'report'

          if (success) then
             write (msg, '(3a)') trim(test_name), ' PASS : ', trim(name)
             call log_message(msg, this_module, this_procedure)
          else
             write (msg, '(4a)') 'NEKLAB_TEST_FAIL ', trim(test_name), ' : ', trim(name)
             call log_warning(msg, this_module, this_procedure)
          end if
       end subroutine report

      end module neklab_test