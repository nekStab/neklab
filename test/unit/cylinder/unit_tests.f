c-----------------------------------------------------------------------
c     unit_tests.f -- drivers and settings for the neklab unit tests (2D cylinder)
c
c     NOT compiled on its own. It is pulled into the .usr with
c
c           include 'unit_tests.f'
c
c        1  fixed point     nek_dvector, nek_zvector, exptA_linop,
c                           nek_system + nek_jacobian           on BF.fld
c        2  resolvent       resolvent_linop                     on BF.fld
c        3  periodic orbit  nek_ext_dvector,
c                           nek_upo_system + nek_upo_jacobian   on snapshot.fld
c-----------------------------------------------------------------------
      include 'neklab_test_helpers.f'
c-----------------------------------------------------------------------
      module unit_tests_cyl
      use LightKrylov
      use neklab
      use neklab_test_helpers
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
      call check_linop(this_test, ok, maxerr, test_adjoint, lin_tol, adj_tol)

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
      call check_taylor(this_test, orders, taylor_min_order)

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
      call check_linop(this_test, ok, maxerr, test_adjoint, lin_tol, adj_tol)
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
      call check_linop(this_test, ok, maxerr, test_adjoint, lin_tol_res,
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
      call check_taylor(this_test, orders, taylor_min_order)

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
      call check_linop(this_test, ok, maxerr, test_adjoint, lin_tol,
     &                 adj_tol)
      end subroutine test_periodic_orbit

      end module unit_tests_cyl
