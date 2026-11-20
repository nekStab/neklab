c-----------------------------------------------------------------------
c     unit_tests.f -- drivers and helpers for the neklab unit tests
c                     (2D cylinder)
c
c     NOT compiled on its own. It is pulled into the .usr with
c
c           include 'unit_tests.f'
c
c     placed ABOVE subroutine userchk (a module must come before the
c     code that uses it). Fixed-form source, like the .usr itself.
c
c        1  fixed point     nek_dvector, nek_zvector, exptA_linop,
c                           nek_system + nek_jacobian           on BF.fld
c        2  resolvent       resolvent_linop                     on BF.fld
c        3  periodic orbit  nek_ext_dvector,
c                           nek_upo_system + nek_upo_jacobian   on snapshot.fld
c-----------------------------------------------------------------------
      module unit_tests_cyl
      use LightKrylov
      use neklab
      implicit none
      private
      character(len=128), parameter :: this_module = 'unit_tests_cyl'

      public :: test_fixed_point, test_resolvent, test_periodic_orbit
      public :: tests_completed
      public :: check, report, check_linop, check_taylor

c     ---- Settings (defaults; can be overwritten in userchk before
c          calling a driver, e.g.  tau = 0.2d0) -------------------------
      character(len=128), public :: base_flow = 'BF.fld'
      character(len=128), public :: snapshot  = 'snapshot.fld'
      logical,  public :: test_zero_map    = .true.     ! test zero invariant map
      logical,  public :: test_adjoint     = .true.     ! test adjoint compatibility
      logical,  public :: verb             = .true.     ! verbose output
      real(dp), public :: tau              = 0.1_dp     ! integration time of exptA
      real(dp), public :: omega            = 8.0_dp*atan(1.0_dp)/0.1_dp ! frequency
      real(dp), public :: period           = 0.1_dp     ! guessed period
      real(dp), public :: cfl_limit        = 0.4_dp
      real(dp), public :: ptol             = 1.0e-08_dp ! pressure tolerance
      real(dp), public :: p2v_factor       = 0.01_dp    ! vtol = ptol*p2v_factor
      integer,  public :: ntrials_vec      = 10         ! vector axioms (cheap)
      integer,  public :: ntrials_lin      = 2          ! operator axioms (time integrations)
      real(dp), public :: vec_tol          = 1.0e-11_dp ! vector axioms
      real(dp), public :: lin_tol          = 1.0e-06_dp ! operator axioms (100*ptol)
      real(dp), public :: lin_tol_res      = 1.0e-4_dp  ! resolvent (GMRES tol.)
      real(dp), public :: adj_tol          = 1.0e-2_dp  ! operator axioms
      real(dp), public :: eps0             = 1.0_dp     ! initial perturbation for Taylor test
      real(dp), public :: taylor_min_order = 1.8_dp
c     -------------------------------------------------------------------

      contains

c=======================================================================
c     DRIVERS
c=======================================================================

      subroutine test_fixed_point()
      include 'SIZE'
      include 'TOTAL'
      type(nek_dvector), allocatable :: dvec, bf
      type(nek_zvector), allocatable :: zvec
      type(exptA_linop), allocatable :: exptA
      type(nek_system),  allocatable :: sys
      real(dp), allocatable :: orders(:)
      real(dp) :: maxerr(4)
      logical  :: ok
      character(len=128) :: this_test

      allocate(dvec, zvec, bf)

      ! Base flow
      call load_fld(trim(base_flow))
      call nek2vec(bf, vx, vy, vz, pr, t)

      ! ---- nek_dvector
      this_test = 'nek_dvector.axioms'
      call dvec%zero()
      ok = verify_vector_axioms(dvec, 
     &                      ntrials          = ntrials_vec,
     &                      tolerance        = vec_tol,
     &                      verbose          = verb)
      call check(this_test, ok)

      ! ---- nek_zvector
      this_test = 'nek_zvector.axioms'
      call zvec%zero()
      ok = verify_vector_axioms(zvec, 
     &                      ntrials          = ntrials_vec,
     &                      tolerance        = vec_tol,
     &                      verbose          = verb)
      call check(this_test, ok)

      ! ---- exponential propagator
      exptA = exptA_linop(tau = tau, baseflow = bf)
      call exptA%init(cfl_limit=cfl_limit, ptol=ptol, p2v_factor=p2v_factor)

      this_test = 'exptA_linop.axioms'
      ok = verify_linop_axioms(exptA, bf,
     &                      ntrials          = ntrials_lin,
     &                      tolerance        = lin_tol,
     &                      adj_tolerance    = adj_tol,
     &                      test_zero_map    = test_zero_map,
     &                      test_adjoint     = test_adjoint,
     &                      verbose          = verb,
     &                      stop_on_failure  = .false.,
     &                      max_errors       = maxerr)
      call check_linop(this_test, ok, maxerr, .true., lin_tol, adj_tol)

      ! ---- nek_system
      sys = nek_system()
      sys%jacobian = nek_jacobian()
      sys%jacobian%X = bf
      call sys%init(cfl_limit=cfl_limit, p2v_factor=p2v_factor, newton_tol_factor=0.1_dp)

      this_test = 'nek_system.linearization'
      call verify_system_linearization(sys,
     &                      nsteps           = 6,
     &                      eps0             = eps0,
     &                      stepratio        = sqrt(0.1_dp),
     &                      tolerance        = 1.0e-8_dp,
     &                      verbose          = verb,
     &                      orders           = orders)
      call check_taylor(this_test, orders)

      this_test = 'nek_system%jacobian.axioms'
      ok = verify_linop_axioms(sys%jacobian, bf,
     &                      ntrials          = ntrials_lin,
     &                      tolerance        = lin_tol,
     &                      adj_tolerance    = adj_tol,
     &                      test_zero_map    = test_zero_map,
     &                      test_adjoint     = test_adjoint,
     &                      verbose          = verb,
     &                      stop_on_failure  = .false.,
     &                      max_errors       = maxerr)
      call check_linop(this_test, ok, maxerr, .true., lin_tol, adj_tol)
      end subroutine test_fixed_point

c-----------------------------------------------------------------------
      subroutine test_resolvent()
      include 'SIZE'
      include 'TOTAL'
      type(nek_dvector),     allocatable :: bf
      type(nek_zvector),     allocatable :: zvec
      type(resolvent_linop), allocatable :: resolvent
      real(dp) :: maxerr(4)
      logical  :: ok
      character(len=128) :: this_test

      call nek_log_message('Currently not available. Exiting.', this_module, 'test_resolvent')
      call nek_end()

      allocate(bf, zvec)

      ! Base flow
      call load_fld(trim(base_flow))
      call nek2vec(bf, vx, vy, vz, pr, t)

      ! ---- nek_zvector
      this_test = 'nek_zvector.axioms'
      call zvec%zero()
      ok = verify_vector_axioms(zvec,
     &                      ntrials          = ntrials_vec,
     &                      tolerance        = vec_tol,
     &                      verbose          = verb)
      call check(this_test, ok)

      ! ---- resolvent operator
      resolvent = resolvent_linop(omega=omega, baseflow=bf)
      call resolvent%init(bf=bf, omega=omega, cfl_limit=cfl_limit,
     &                    ptol=ptol, p2v_factor=p2v_factor)

      this_test = 'resolvent_linop.axioms'
      ok = verify_linop_axioms(resolvent, zvec,
     &                      ntrials          = 1,
     &                      tolerance        = lin_tol_res,
     &                      adj_tolerance    = adj_tol,
     &                      test_zero_map    = test_zero_map,
     &                      test_adjoint     = test_adjoint,
     &                      verbose          = verb,
     &                      stop_on_failure  = .false.,
     &                      max_errors       = maxerr)
      call check_linop(this_test, ok, maxerr, .true., lin_tol_res,
     &                 adj_tol)
      end subroutine test_resolvent

c-----------------------------------------------------------------------
      subroutine test_periodic_orbit()
      include 'SIZE'
      include 'TOTAL'
      type(nek_ext_dvector), allocatable :: evec, x0
      type(nek_upo_system),  allocatable :: sys_upo
      real(dp), allocatable :: orders(:)
      real(dp) :: maxerr(4)
      logical  :: ok
      character(len=128) :: this_test

      allocate(x0, evec)

      ! Snapshot + guessed period
      call load_fld(trim(snapshot))
      call nek2ext_vec(x0, vx, vy, vz, pr, t)
      x0%T = period

      ! ---- nek_ext_dvector
      this_test = 'nek_ext_dvector.axioms'
      call evec%zero()
      ok = verify_vector_axioms(evec, 
     &                      ntrials          = ntrials_vec,
     &                      tolerance        = vec_tol,
     &                      verbose          = verb)
      call check(this_test, ok)

      ! ---- nek_upo_system
      sys_upo = nek_upo_system()
      sys_upo%jacobian = nek_upo_jacobian()
      sys_upo%jacobian%X = x0
      call sys_upo%init(cfl_limit=cfl_limit, p2v_factor=p2v_factor, newton_tol_factor=0.1_dp)

      this_test = 'nek_upo_system.linearization'
      call verify_system_linearization(sys_upo,
     &                      nsteps           = 6,
     &                      eps0             = eps0,
     &                      stepratio        = sqrt(0.1_dp),
     &                      tolerance        = 1.0e-8_dp,
     &                      verbose          = verb,
     &                      orders           = orders)
      call check_taylor(this_test, orders)

      this_test = 'nek_upo_system%jacobian.axioms'
      ok = verify_linop_axioms(sys_upo%jacobian, x0,
     &                      ntrials          = ntrials_lin,
     &                      tolerance        = lin_tol,
     &                      adj_tolerance    = adj_tol,
     &                      test_zero_map    = test_zero_map,
     &                      test_adjoint     = test_adjoint,
     &                      verbose          = verb,
     &                      stop_on_failure  = .false.,
     &                      max_errors       = maxerr)
      call check_linop(this_test, ok, maxerr, .true., lin_tol, adj_tol)
      end subroutine test_periodic_orbit

c-----------------------------------------------------------------------
      subroutine tests_completed()
      call nek_log_message('############################', procedure='main')
      call nek_log_message('#     tests completed      #', procedure='main')
      call nek_log_message('############################', procedure='main')
      call nek_end()
      end subroutine tests_completed

c=======================================================================
c     HELPERS
c=======================================================================

      subroutine check(name, ok, tol, info)
      use stdlib_strings, only : padr
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      real(dp), intent(in), optional :: tol
      character(len=*), intent(in), optional :: info
      ! internal
      character(len=*), parameter :: esc   = achar(27)
      character(len=*), parameter :: green = esc//'[32m'
      character(len=*), parameter :: red   = esc//'[31m'
      character(len=*), parameter :: reset = esc//'[0m'
      character(len=*), parameter :: pfx = 'NEKLAB_CHECK '
      character(len=*), parameter :: p   = green//'PASSED'//reset
      character(len=*), parameter :: f   = red//'FAILED'//reset
      character(len=512) :: msg
      character(len=45) :: ntr
      character(len=20) :: tstr
      ntr = padr(trim(name),45)
      call nek_log_message('')
      call nek_log_message('#####################################################', procedure='main')
      call nek_log_message('#     test '//trim(name), procedure='main')
      call nek_log_message('#####################################################', procedure='main')
      call nek_log_message('')
      tstr = ''
      if (present(tol)) then
         write(tstr, '(A6,ES10.2,4X)') 'tol = ', tol
      end if
      if (ok) then
         write(msg, '(A,A45,A20,A)')      pfx, ntr, tstr, p
      else if (present(info)) then
         write(msg, '(A,A45,A20,A,1X,A)') pfx, ntr, tstr, f, trim(info)
      else
         write(msg, '(A,A45,A20,A)')      pfx, ntr, tstr, f
      end if
      call nek_log_message(msg, procedure='main')
      call nek_log_message('')
      end subroutine check

c-----------------------------------------------------------------------
      subroutine report(key, value)
      use stdlib_strings, only : padr
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: value
      character(len=512) :: msg
      write(msg, '(A,A40,2X,ES15.5)') 'NEKLAB_RESULT ',
     &                                padr(trim(key), 40), value
      call nek_log_message(msg, procedure='main')
      end subroutine report

c-----------------------------------------------------------------------
      subroutine check_linop(name, ok, maxerr, with_adjoint, tol, a_tol)
c     verify_linop_axioms result + the largest error of each axiom
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      real(dp), intent(in) :: maxerr(4)
      logical, intent(in) :: with_adjoint
      real(dp), intent(in) :: tol
      real(dp), intent(in) :: a_tol
      logical :: chk
      character(len=128) :: msg
      call report(trim(name)//'.max_err_additivity ', maxerr(1))
      call report(trim(name)//'.max_err_homogeneity', maxerr(2))
      call report(trim(name)//'.max_err_zero_map   ', maxerr(3))
      write(msg, '(A,ES9.2,A,ES9.2)') 'max error ', maxval(maxerr(1:3)),' > lin_tol ', tol
      chk = all(maxerr(1:3) <= tol)
      call check(trim(name)//'.linearity', chk, tol=tol, info=msg)
      if (with_adjoint) then
         call report(trim(name)//'.max_err_adjoint', maxerr(4))
         write(msg, '(A,ES9.2,A,ES9.2)') 'error ', maxerr(4), ' > adj_tol ', a_tol
         chk = maxerr(4) <= a_tol
         call check(trim(name)//'.adjoint', chk, tol=a_tol, info=msg)
      end if
      call check('GLOBAL : '//name, ok)
      call nek_log_message('')
      end subroutine check_linop

c-----------------------------------------------------------------------
      subroutine check_taylor(name, orders)
c     pass if at least half of the observed orders reach taylor_min_order
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: orders(:)
      character(len=128) :: key, msg
      integer :: k, ngood
      logical :: chk
      call nek_log_message('')
      do k = 1, size(orders)
         write(key, '(2A,I0)') trim(name), '.order_', k
         call report(key, orders(k))
      end do
      ngood = count(orders >= taylor_min_order)
      write(msg, '(I0,A,I0,A,F4.2)') ngood, ' of ', size(orders), ' orders >= ', taylor_min_order
      chk = (size(orders) > 0 .and. ngood >= size(orders)/2)
      call check(name, chk, info=msg)
      call nek_log_message('')
      end subroutine check_taylor

      end module unit_tests_cyl
