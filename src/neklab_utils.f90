      module neklab_utils
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Logger
      ! Abstract types for real-valued vectors.
         use LightKrylov, only: abstract_vector_rdp
         use LightKrylov_Utils, only: abstract_opts
      ! Neklab vectors
         use neklab_vectors
         use neklab_nek_setup, only: setup_nek, nek_stop_error, nek_log_message
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_utils'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
         integer, parameter :: lt = lx1*ly1*lz1*lelt
      !! Local number of grid points for the temperature/passive scalar mesh.
      
      ! utilities for regular nek vectors
         public :: nek2vec, vec2nek, abs_vec2nek, outpost_dnek, outpost_dnek_abs_vector
      ! utilities for extended nek vectors
         public :: nek2ext_vec, ext_vec2nek, abs_ext_vec2nek, outpost_ext_dnek
         public :: get_period, get_period_abs
      ! miscellaneous
         public :: set_nek_opts, nopcopy

		! nek5000 setup options
         type, extends(abstract_opts), public :: abstract_nek_opts
            !! General Nek5000 solver options
				logical :: recompute_dt = .true.
				! Recompute timestep based on baseflow and CFL limit (default: .true.).
				logical :: variable_dt = .false.
				! Allow for variable timestep during integration (default: .false.).
				real(dp) :: endtime = 0.0
				! Integration time (default: taken from .par file)
				real(dp) :: ptol = 0.0_dp
				! tolerance setting for the pressure solves (default: taken from .par file)
				real(dp) :: vtol = 0.0_dp
				! tolerance setting for the velocity solves (default: taken from .par file)
				real(dp) :: cfl_limit = 0.5_dp	
				! CFL limit used to determine maximum dt (default: taken from .par file)
				logical, private :: initialized = .false.
         end type

         type, extends(abstract_nek_opts), public :: nek_opts_std
         contains
            procedure :: init           => init_nek_opts_std
            procedure :: is_initialized => is_initialized_nek_opts_std
            procedure :: print          => print_nek_opts_std
         end type
            
         type, extends(abstract_nek_opts), public :: nek_opts_prt
            !! Specific Nek5000 solver options for linear simulations
				logical :: solve_baseflow = .false.
				! Solve baseflow and perturbations simultaneously (default: .false.).
         contains
            procedure :: init           => init_nek_opts_prt
            procedure :: is_initialized => is_initialized_nek_opts_prt
            procedure :: print          => print_nek_opts_prt
         end type
      
      ! Nek vector utilities
         interface nek2vec
            module procedure nek2vec_std
            module procedure nek2vec_prt
         end interface
      
         interface vec2nek
            module procedure vec2nek_std
            module procedure vec2nek_prt
         end interface
      
         interface abs_vec2nek
            module procedure abstract_vec2nek_std
            module procedure abstract_vec2nek_prt
         end interface
      
         interface outpost_dnek
            module procedure outpost_dnek_vector
            module procedure outpost_dnek_basis
         end interface
      
      ! Extended nek vector utilities
         interface nek2ext_vec
            module procedure nek2ext_vec_std
            module procedure nek2ext_vec_prt
         end interface
      
         interface ext_vec2nek
            module procedure ext_vec2nek_std
            module procedure ext_vec2nek_prt
         end interface
      
         interface abs_ext_vec2nek
            module procedure abstract_ext_vec2nek_std
            module procedure abstract_ext_vec2nek_prt
         end interface
      
         interface outpost_ext_dnek
            module procedure outpost_ext_dnek_vector
            module procedure outpost_ext_dnek_basis
         end interface

		! Set nek opts
			interface set_nek_opts
            module procedure set_nek_opts_std
            module procedure set_nek_opts_prt
         end interface

      contains

      !---------------------------------------------------------------------
      ! Type bound procedures for opts
      !---------------------------------------------------------------------

			subroutine init_nek_opts_std(self, recompute_dt, variable_dt, endtime, ptol, vtol, cfl_limit)
				class(nek_opts_std), intent(inout) :: self
            logical, optional, intent(in) :: recompute_dt
				logical, optional, intent(in) :: variable_dt
				real(dp), optional, intent(in) :: endtime
				real(dp), optional, intent(in) :: ptol
				real(dp), optional, intent(in) :: vtol
				real(dp), optional, intent(in) :: cfl_limit
            ! internal
            character(len=*), parameter :: this_procedure = 'init_nek_opts_std'
            
            ! Defaults
            self%recompute_dt = optval(recompute_dt, .true.)
            self%variable_dt  = optval(variable_dt,  .false.)
            ! involving param
				self%endtime      = optval(endtime,      param(10))
				self%ptol         = optval(ptol,         param(21))
				self%vtol         = optval(vtol,         param(22))
				self%cfl_limit    = optval(cfl_limit,    param(26))
            
            if (.not. self%initialized) then
               call nek_log_message('options structure initialized.', this_module, this_procedure)
            end if

				self%initialized  = .true.
			end subroutine init_nek_opts_std

         pure logical function is_initialized_nek_opts_std(self) result(is_initialized)
            class(nek_opts_std), intent(in) :: self
            is_initialized = self%initialized
         end function is_initialized_nek_opts_std

         subroutine print_nek_opts_std(self)
            class(nek_opts_std), intent(in) :: self
            ! internal
            character(len=128) :: msg
            ! print info
            call nek_log_message('##  NEK OPTS STD ##', module=this_module)
            write (msg, '(A,L15)') padl('recompute_dt:', 20), self%recompute_dt
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,L15)') padl('variable_dt:', 20), self%variable_dt
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('endtime:', 20), self%endtime
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('cfl_limit:', 20), self%cfl_limit
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')            
            write (msg, '(A,E15.8)') padl('vtol:', 20), self%vtol
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,E15.8)') padl('ptol:', 20), self%ptol
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            call nek_log_message('##  NEK OPTS STD ##', module=this_module)
         end subroutine print_nek_opts_std

         subroutine init_nek_opts_prt(self, solve_baseflow, recompute_dt, variable_dt, endtime, ptol, vtol, cfl_limit)
				class(nek_opts_prt), intent(inout) :: self
            logical, optional, intent(in) :: solve_baseflow
            logical, optional, intent(in) :: recompute_dt
				logical, optional, intent(in) :: variable_dt
				real(dp), optional, intent(in) :: endtime
				real(dp), optional, intent(in) :: ptol
				real(dp), optional, intent(in) :: vtol
				real(dp), optional, intent(in) :: cfl_limit
            ! internal
            character(len=*), parameter :: this_procedure = 'init_nek_opts_prt'
            character(len=128) :: msg
            ! Defaults
            self%solve_baseflow = optval(solve_baseflow, .false.)
            self%recompute_dt   = optval(recompute_dt, .true.)
            self%variable_dt    = optval(variable_dt, .false.)
            ! involving param
				self%endtime        = optval(endtime,      param(10))
				self%ptol           = optval(ptol,         param(21))
				self%vtol           = optval(vtol,         param(22))
				self%cfl_limit      = optval(cfl_limit,    param(26))

            if (.not. self%initialized) then
               call nek_log_message('options structure initialized.', this_module, this_procedure)
            end if

				self%initialized = .true.
			end subroutine init_nek_opts_prt

         pure logical function is_initialized_nek_opts_prt(self) result(is_initialized)
            class(nek_opts_prt), intent(in) :: self
            is_initialized = self%initialized
         end function is_initialized_nek_opts_prt

         subroutine print_nek_opts_prt(self)
            class(nek_opts_prt), intent(in) :: self
            ! internal
            character(len=128) :: msg
            ! print info
            call nek_log_message('##  NEK OPTS PRT ##', module=this_module)
            write (msg, '(A,L15)') padl('solve_baseflow:', 20), self%solve_baseflow
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,L15)') padl('recompute_dt:', 20), self%recompute_dt
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,L15)') padl('variable_dt:', 20), self%variable_dt
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('endtime:', 20), self%endtime
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('cfl_limit:', 20), self%cfl_limit
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')            
            write (msg, '(A,E15.8)') padl('vtol:', 20), self%vtol
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,E15.8)') padl('ptol:', 20), self%ptol
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            call nek_log_message('##  NEK OPTS PRT ##', module=this_module)
         end subroutine print_nek_opts_prt

         subroutine set_nek_opts_std(opts, silent)
				type(nek_opts_std), intent(inout) :: opts
				logical, optional, intent(in) :: silent
            ! internal
            logical :: silent_
            
            ! optional arguments
            silent_ = optval(silent, .false.)

            ! safety check
            if (.not. opts%is_initialized()) call opts%init()
            if (.not. silent_) call opts%print()
            
            ! set opts
				call setup_nek(LNS = .false., recompute_dt = opts%recompute_dt, variable_dt = opts%variable_dt, endtime = opts%endtime, vtol = opts%vtol, ptol = opts%ptol, cfl_limit = opts%cfl_limit, silent = silent)

			end subroutine set_nek_opts_std

			subroutine set_nek_opts_prt(opts, transpose, silent)
				type(nek_opts_prt), intent(inout) :: opts
				logical, optional, intent(in) :: transpose
				logical, optional, intent(in) :: silent
            ! internal
            logical :: silent_
            
            ! optional arguments
            silent_ = optval(silent, .false.)

            ! safety check
            if (.not. opts%is_initialized()) call opts%init()
            if (.not. silent_) call opts%print()

            ! set opts
				call setup_nek(LNS = .true., transpose = transpose, solve_baseflow = opts%solve_baseflow, recompute_dt = opts%recompute_dt, variable_dt = opts%variable_dt, endtime = opts%endtime,vtol = opts%vtol, ptol = opts%ptol, cfl_limit = opts%cfl_limit, silent = silent)

			end subroutine set_nek_opts_prt

         subroutine nek2vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_

            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))

         end subroutine nek2vec_prt
      
         subroutine nek2vec_std(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_

            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)

         end subroutine nek2vec_std
      
         subroutine vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
            call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine vec2nek_std
      
         subroutine vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
      
            call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine vec2nek_prt
      
         subroutine abstract_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
            select type (vec)
            type is (nek_dvector)

               call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_dvector','IN',this_module,'abstract_vec2nek_std')
            end select
         end subroutine abstract_vec2nek_std
      
         subroutine abstract_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
            select type (vec)
            type is (nek_dvector)

               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_dvector','IN',this_module,'abstract_vec2nek_prt')
            end select
         end subroutine abstract_vec2nek_prt
      
      ! EXTENDED
      
         subroutine nek2ext_vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, 1), intent(in) :: vx_
            real(kind=dp), dimension(lv, 1), intent(in) :: vy_
            real(kind=dp), dimension(lv, 1), intent(in) :: vz_
            real(kind=dp), dimension(lp, 1), intent(in) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
         end subroutine nek2ext_vec_prt
      
         subroutine nek2ext_vec_std(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
      
         end subroutine nek2ext_vec_std
      
         subroutine ext_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
            call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine ext_vec2nek_std
      
         subroutine ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
      
            call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine ext_vec2nek_prt
      
         subroutine abstract_ext_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
            select type (vec)
            type is (nek_ext_dvector)

               call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_ext_dvector','IN',this_module,'abstract_ext_vec2nek_std')
            end select
         end subroutine abstract_ext_vec2nek_std
      
         subroutine abstract_ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
            select type (vec)
            type is (nek_ext_dvector)

               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_ext_dvector','IN',this_module,'abstract_ext_vec2nek_prt')
            end select
         end subroutine abstract_ext_vec2nek_prt
      
         real(dp) function get_period_abs(vec) result(period)
            class(abstract_vector_rdp), intent(in) :: vec
            select type (vec)
            type is (nek_ext_dvector)

               period = vec%T

            class default
               call type_error('vec','nek_ext_dvector','IN',this_module,'get_period_abs')
            end select
         end function get_period_abs
      
         pure real(dp) function get_period(vec) result(period)
            class(nek_ext_dvector), intent(in) :: vec
            period = vec%T
         end function get_period
      
         subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
            implicit none
            include 'SIZE'
            include 'TOTAL'
            integer :: n, k
            real(kind=dp), intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
            real(kind=dp), intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
            n = nx1*ny1*nz1*nelv

            call copy(a1, b1, n)
            call copy(a2, b2, n)
            if (if3D) call copy(a3, b3, n)

            if (ifpo) call copy(a4, b4, nx2*ny2*nz2*nelv)

            if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))

            if (ldimt > 1) then
            do k = 1, npscal
               if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
            end do
            end if
         end subroutine nopcopy
      
         subroutine outpost_dnek_vector(vec, prefix)
            type(nek_dvector), intent(in) :: vec
            character(len=3), intent(in) :: prefix

            call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)

         end subroutine outpost_dnek_vector
      
         subroutine outpost_dnek_abs_vector(vec, prefix)
            class(abstract_vector_rdp), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            select type (vec)
            type is (nek_dvector)

               call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
               
            class default
               call type_error('vec','nek_dvector','IN',this_module,'outpost_dnek_abs_vector')
            end select
         end subroutine outpost_dnek_abs_vector
      
         subroutine outpost_dnek_basis(vec, prefix)
            type(nek_dvector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            integer :: i
            do i = 1, size(vec)
               call outpost_dnek_vector(vec(i), prefix)
            end do
         end subroutine outpost_dnek_basis
      
         subroutine outpost_ext_dnek_vector(vec, prefix)
            type(nek_ext_dvector), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
         end subroutine outpost_ext_dnek_vector
      
         subroutine outpost_ext_dnek_basis(vec, prefix)
            type(nek_ext_dvector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            integer :: i
            do i = 1, size(vec)
               call outpost_ext_dnek_vector(vec(i), prefix)
            end do
         end subroutine outpost_ext_dnek_basis

      end module neklab_utils
