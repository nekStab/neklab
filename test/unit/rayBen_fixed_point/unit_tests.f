c-----------------------------------------------------------------------
c     unit_tests.f -- drivers and settings for the neklab unit tests
c                     (Rayleigh-Benard convection, with temperature)
c
c     NOT compiled on its own. It is pulled into the .usr with
c
c           include 'unit_tests.f'
c
c     test_fixed_point_temp: nek_dvector with temperature,
c        exptA_temp_linop, nek_temp_system + nek_temp_jacobian,
c        on the conduction state u = 0, T = 1 - y (or on base_flow).
c-----------------------------------------------------------------------
      include 'neklab_test_helpers.f'
c-----------------------------------------------------------------------
      module unit_tests_rbc
      use LightKrylov
      use neklab
      use neklab_test_helpers
      implicit none
      private
      character(len=128), parameter :: this_module = 'unit_tests_rbc'

      public :: test_fixed_point_temp
      public :: tests_completed, log_settings, run_tests
      public :: check, report, check_linop, check_taylor

c     ---- Settings (defaults; can be overwritten in userchk before
c          calling a driver, or from test/suites/unit.py) ---------------
      logical,  public :: test_zero_map    = .true.     ! test zero invariant map
      logical,  public :: test_adjoint     = .true.     ! test adjoint compatibility
      logical,  public :: verb             = .true.     ! verbose output
      real(dp), public :: tau              = 0.1_dp     ! integration time of exptA
      real(dp), public :: cfl_limit        = 0.4_dp
      real(dp), public :: ptol             = 1.0e-08_dp ! pressure tolerance
      real(dp), public :: p2v_factor       = 0.01_dp    ! vtol = ptol*p2v_factor
      integer,  public :: ntrials_vec      = 10         ! vector axioms (cheap)
      integer,  public :: ntrials_lin      = 2          ! operator axioms (time integrations)
      real(dp), public :: vec_tol          = 1.0e-10_dp ! vector axioms
      real(dp), public :: lin_tol          = 1.0e-06_dp ! operator axioms (100*ptol)
      real(dp), public :: adj_tol          = 1.0e-2_dp  ! operator axioms
      real(dp), public :: eps0             = 1.0_dp     ! initial perturbation for Taylor test
      real(dp), public :: taylor_tol       = 1.0e-8_dp  ! solver atol in the Taylor test
      real(dp), public :: taylor_min_order = 1.8_dp
c     -------------------------------------------------------------------

      contains

c=======================================================================
c     DRIVERS
c=======================================================================

      subroutine test_fixed_point_temp()
      include 'SIZE'
      include 'TOTAL'
      include 'ADJOINT'
      type(nek_dvector),      allocatable :: dvec, bf
      type(exptA_temp_linop), allocatable :: exptA
      type(nek_temp_system),  allocatable :: sys
      real(dp), allocatable :: orders(:)
      real(dp) :: maxerr(4)
      logical  :: ok
      integer  :: i
      character(len=128) :: this_test
      integer, parameter :: lv = lx1*ly1*lz1*lelv
      integer, parameter :: lp = lx2*ly2*lz2*lelv

      allocate(dvec, bf)

      ! ---- base flow (analytical)
      call rzero(vx, lv)
      call rzero(vy, lv)
      call rzero(vz, lv)
      call rzero(pr, lp)
      do i = 1, lv
         t(i,1,1,1,1) = 1.0_dp - ym1(i,1,1,1)
      end do
      call nek2vec(bf, vx, vy, vz, pr, t)
      
      ! reset baseflow temperature gradient for adjoint computation
      if (test_adjoint) call gradm1(dTdx,dTdy,dTdz,t)

      ! ---- nek_dvector (with temperature)
      this_test = 'nek_dvector_temp.axioms'
      call dvec%zero()
      ok = verify_vector_axioms(dvec,
     &                      ntrials          = ntrials_vec,
     &                      tolerance        = vec_tol,
     &                      verbose          = verb)
      call check(this_test, ok)

      ! ---- exponential propagator with temperature
      exptA = exptA_temp_linop(tau = tau, baseflow = bf)
      call exptA%init(cfl_limit=cfl_limit, ptol=ptol, p2v_factor=p2v_factor, recompute_dt=.false.)

      this_test = 'exptA_temp_linop.axioms'
      ok = verify_linop_axioms(exptA, bf,
     &                      ntrials          = ntrials_lin,
     &                      tolerance        = lin_tol,
     &                      adj_tolerance    = adj_tol,
     &                      test_zero_map    = test_zero_map,
     &                      test_adjoint     = test_adjoint,
     &                      verbose          = verb,
     &                      stop_on_failure  = .false.,
     &                      max_errors       = maxerr)
      call check_linop(this_test, ok, maxerr, test_adjoint, lin_tol, adj_tol)

      ! ---- nek_temp_system (integration time = endTime of the .par)
      sys = nek_temp_system()
      sys%jacobian = nek_temp_jacobian()
      sys%jacobian%X = bf
      call sys%init(cfl_limit=cfl_limit, p2v_factor=p2v_factor, newton_tol_factor=0.1_dp, recompute_dt=.false.)

      this_test = 'nek_temp_system.linearization'
      call verify_system_linearization(sys,
     &                      nsteps           = 6,
     &                      eps0             = eps0,
     &                      stepratio        = sqrt(0.1_dp),
     &                      tolerance        = taylor_tol,
     &                      verbose          = verb,
     &                      orders           = orders)
      call check_taylor(this_test, orders, taylor_min_order)

      this_test = 'nek_temp_system%jacobian.axioms'
      ok = verify_linop_axioms(sys%jacobian, bf,
     &                      ntrials          = ntrials_lin,
     &                      tolerance        = lin_tol,
     &                      adj_tolerance    = adj_tol,
     &                      test_zero_map    = test_zero_map,
     &                      test_adjoint     = test_adjoint,
     &                      verbose          = verb,
     &                      stop_on_failure  = .false.,
     &                      max_errors       = maxerr)
      call check_linop(this_test, ok, maxerr, test_adjoint, lin_tol, adj_tol)
      end subroutine test_fixed_point_temp

c-----------------------------------------------------------------------
      subroutine run_tests(group)
c     Run one group of checks; the test harness sets the group in the .par
c     (userParam01). To add a test: write a driver above and add a case.
      integer, intent(in) :: group
      character(len=64) :: msg
      select case (group)
      case (1)
         call test_fixed_point_temp()
      case default
         write(msg, '(A,I0,A)') 'unknown group ', group,
     &                          ' (userParam01 in the .par: 1)'
         call check('setup.group', .false., info=msg)
      end select
      end subroutine run_tests
c-----------------------------------------------------------------------
      subroutine log_settings()
c     Write every setting and the Nek5000 setup to the log, so the log
c     records what this run used (the test harness edits the defaults
c     above in its copy of this file).
      call log_section('unit test settings (unit_tests_rbc)')
      call log_setting('test_zero_map', test_zero_map)
      call log_setting('test_adjoint', test_adjoint)
      call log_setting('verb', verb)
      call log_setting('tau', tau)
      call log_setting('cfl_limit', cfl_limit)
      call log_setting('ptol', ptol)
      call log_setting('p2v_factor', p2v_factor)
      call log_setting('ntrials_vec', ntrials_vec)
      call log_setting('ntrials_lin', ntrials_lin)
      call log_setting('vec_tol', vec_tol)
      call log_setting('lin_tol', lin_tol)
      call log_setting('adj_tol', adj_tol)
      call log_setting('eps0', eps0)
      call log_setting('taylor_tol', taylor_tol)
      call log_setting('taylor_min_order', taylor_min_order)
      end subroutine log_settings

      end module unit_tests_rbc
