c-----------------------------------------------------------------------
c     neklab_test_helpers.f -- reporting helpers shared by all tests
c
c     The .usr of every test case pulls it in with its first line
c
c           include 'neklab_test_helpers.f'
c
c     (each case folder holds a symlink to test/common/). Fixed-form
c     source, like the .usr itself.
c
c     check(name, ok)      NEKLAB_CHECK <name> [tol = x] PASSED|FAILED [info]
c                          (names of at most 45 characters)
c     report(key, value)   NEKLAB_RESULT <key> <value>, 17 significant digits
c                          (keys of at most 50 characters)
c     tests_completed()    the end marker the harness waits for
c
c     Besides the check/report helpers, it logs what a run used, so every
c     log (lightkrylov.log and the Nek logfile) records its own setup:
c        log_section('title'), log_setting(name, value)  -- one per setting
c-----------------------------------------------------------------------
      module neklab_test_helpers
      use LightKrylov, only: dp
      use neklab
      implicit none
      private

      public :: check, report, check_linop, check_taylor
      public :: tests_completed
      public :: log_section, log_setting

      interface log_setting
         module procedure log_setting_r, log_setting_i
         module procedure log_setting_l, log_setting_c
      end interface

      contains

c-----------------------------------------------------------------------
      subroutine check(name, ok, tol, info)
c     one PASSED/FAILED line; info is shown only when the check fails
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
c     one value for measured_values.json, in full precision: the harness
c     compares it with the references to tolerances down to 1e-6 and below
      use stdlib_strings, only : padr
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: value
      character(len=512) :: msg
      write(msg, '(A,A50,2X,ES24.16)') 'NEKLAB_RESULT ', padr(trim(key), 50), value
      call nek_log_message(msg, procedure='main')
      end subroutine report

c-----------------------------------------------------------------------
      subroutine check_linop(name, ok, maxerr, with_adjoint, tol, a_tol)
c     verify_linop_axioms result + the largest error of each axiom:
c       linearity (additivity, homogeneity, zero map) against tol,
c       adjoint consistency against a_tol (only if with_adjoint)
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
      call nek_log_message('')
      end subroutine check_linop

c-----------------------------------------------------------------------
      subroutine check_taylor(name, orders, min_order)
c     Taylor test: pass if at least half of the observed orders reach
c     min_order (the smallest steps may level off at the solver tolerance)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: orders(:)
      real(dp), intent(in) :: min_order
      character(len=128) :: key, msg
      integer :: k, ngood
      logical :: chk
      call nek_log_message('')
      do k = 1, size(orders)
         write(key, '(2A,I0)') trim(name), '.order_', k
         call report(key, orders(k))
      end do
      ngood = count(orders >= min_order)
      write(msg, '(I0,A,I0,A,F4.2)') ngood, ' of ', size(orders), ' orders >= ', min_order
      chk = (size(orders) > 0 .and. 2*ngood >= size(orders))
      call check(name, chk, info=msg)
      call nek_log_message('')
      end subroutine check_taylor

c-----------------------------------------------------------------------
      subroutine tests_completed()
c     end marker the harness waits for, then stop Nek
      call nek_log_message('############################', procedure='main')
      call nek_log_message('#     tests completed      #', procedure='main')
      call nek_log_message('############################', procedure='main')
      call nek_end()
      end subroutine tests_completed

c=======================================================================
c     LOGGING OF THE SETUP
c=======================================================================

      subroutine log_section(title)
      character(len=*), intent(in) :: title
      call nek_log_message('')
      call nek_log_message('---- '//trim(title)//' ----', procedure='setup')
      end subroutine log_section

      subroutine log_setting_r(name, value)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      character(len=256) :: msg
      write(msg, '(2X,A,T30,"= ",ES14.6)') name, value
      call nek_log_message(msg, procedure='setup')
      end subroutine log_setting_r

      subroutine log_setting_i(name, value)
      character(len=*), intent(in) :: name
      integer, intent(in) :: value
      character(len=256) :: msg
      write(msg, '(2X,A,T30,"= ",I14)') name, value
      call nek_log_message(msg, procedure='setup')
      end subroutine log_setting_i

      subroutine log_setting_l(name, value)
      character(len=*), intent(in) :: name
      logical, intent(in) :: value
      character(len=256) :: msg
      write(msg, '(2X,A,T30,"= ",L14)') name, value
      call nek_log_message(msg, procedure='setup')
      end subroutine log_setting_l

      subroutine log_setting_c(name, value)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: value
      character(len=256) :: msg
      write(msg, '(2X,A,T30,"= ",A14)') name, "'"//trim(value)//"'"
      call nek_log_message(msg, procedure='setup')
      end subroutine log_setting_c

      end module neklab_test_helpers
