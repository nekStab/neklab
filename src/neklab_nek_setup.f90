      module neklab_nek_setup
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
         use stdlib_logger, only: debug_level, warning_level, information_level
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Logger
      ! Abstract types for real-valued vectors.
         use LightKrylov, only: abstract_vector_rdp
      ! Neklab vectors
         use neklab_vectors
 
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_nek_setup'
 
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
         integer, parameter :: lt = lx1*ly1*lz1*lelt
      !! Local number of grid points for the temperature/passive scalar mesh.
 
      ! Set up solver
         public :: setup_nek, setup_nonlinear_solver, setup_linear_solver, nek_status
      ! Utilities for logging
         public :: nek_log_message, nek_log_warning, nek_log_information, nek_log_debug, nek_stop_error
      ! Utilities for file handling
         public :: set_fldindex
 
      contains
 
         subroutine setup_nek(LNS, transpose, solve_baseflow, recompute_dt, variable_dt, endtime, vtol, ptol, cfl_limit, solve_temperature, silent)
            logical, intent(in) :: LNS
            logical, optional, intent(in) :: transpose
            logical :: transpose_
            logical, optional, intent(in) :: solve_baseflow
            logical :: solve_baseflow_
            logical, optional, intent(in) :: recompute_dt
            logical :: recompute_dt_
            logical, optional, intent(in) :: variable_dt
            logical :: variable_dt_
            real(dp), optional, intent(in) :: endtime
            real(dp) :: endtime_
            real(dp), optional, intent(in) :: vtol
            real(dp) :: vtol_
            real(dp), optional, intent(in) :: ptol
            real(dp) :: ptol_
            real(dp), optional, intent(in) :: cfl_limit
            real(dp) :: cfl_limit_
            logical, optional, intent(in) :: solve_temperature
            logical :: solve_temperature_
            logical, optional, intent(in) :: silent
            logical :: silent_
      ! internal
            character(len=*), parameter :: nekfmt = '(5X,A)'
            real(dp) :: dt_old
            integer :: idt, n
            integer, parameter :: max_idt = 200
            character(len=128) :: msg
            logical :: full_summary
            logical :: iffxdt
            common /FIXDT/ iffxdt
            integer, parameter :: pad = 30
            logical :: if_vel_zero, if_t_zero, if_vel_tiny, if_t_tiny
            real(dp), parameter :: tiny_tol = 1e-6_dp 
            real, external :: glamax
 
      ! Only print summary if we switch from linear to nonlinear solvers or vice versa
            full_summary = .false.
 
            transpose_ = optval(transpose, .false.)
            solve_baseflow_    = optval(solve_baseflow,    .false.)
            solve_temperature_ = optval(solve_temperature, .false.)
            endtime_ = optval(endtime, param(10))
            ptol_ = optval(ptol, param(21))
            vtol_ = optval(vtol, param(22))
            cfl_limit_ = optval(cfl_limit, param(26))
            silent_ = optval(silent, .false.)
 
            n = lx1*ly1*lz1*nelv
            if_vel_zero = (glamax(vx, n) == 0.0_dp  .and. glamax(vy, n) == 0.0_dp)
            if_t_zero   = (glamax(t,  n) == 0.0_dp)
            if_vel_tiny = (glamax(vx, n) < tiny_tol .and. glamax(vy, n) < tiny_tol)
            if_t_tiny   = (glamax(t,  n) < tiny_tol)
            if (if3D) then
               if_vel_zero = if_vel_zero .and. glamax(vz, n) == 0.0_dp
               if_vel_tiny = if_vel_tiny .and. glamax(vz, n) < tiny_tol
            end if
            if (if_vel_zero) then
               recompute_dt_ = .false.
               write(msg, '(A,ES8.1,A,F8.6,A)') 'Velocities are zero. dt =', dt, ' (.par)'
               if (ifaxis .and. if_t_zero) then
                  call nek_log_warning(msg, this_module, 'setup_nek')
               else
                  call nek_log_debug(msg, this_module, 'setup_nek')
               end if
            elseif (if_vel_tiny) then
               recompute_dt_  = .false.
               write(msg, '(A,ES8.1,A,F8.6,A)') 'Velocity amplitudes are tiny (max < ', tiny_tol, '). dt =', dt, ' (.par)'
               if (ifaxis .and. if_t_tiny) then
                  call nek_log_warning(msg, this_module, 'setup_nek')
               else
                  call nek_log_debug(msg, this_module, 'setup_nek')
               end if
            else
               recompute_dt_ = optval(recompute_dt, .false.)
            end if
            if (.not. recompute_dt_ .and. dt <= 0.0_dp) then
               write(msg, '(A,ES8.1,A,F8.6,A)') 'dt =', dt, ' (.par) and recompute_dt =', recompute_dt_, '.'
               call nek_stop_error(msg, this_module, 'setup_nek')
            end if
 
            variable_dt_ = optval(variable_dt, .false.)
 
            call nekgsync()
 
            if (nid == 0 .and. .not. silent_) then
               print *, ''
               print '("neklab ",A)', '################## SETUP NEK ###################'
               print *, ''
            end if
      
      ! general settings
            lastep = 0;
            if (LNS) then
               if (.not. ifpert) full_summary = .true.
               ifpert = .true.; call bcast(ifpert, lsize)
               if (transpose_) then
                  ifadj = .true.; call bcast(ifadj, lsize)
               else
                  ifadj = .false.; call bcast(ifadj, lsize)
               end if
               if (solve_baseflow_) then
                  ifbase = .true.; call bcast(ifbase, lsize)
               else
                  ifbase = .false.; call bcast(ifbase, lsize)
               end if
      ! Force single perturbation mode.
               npert = param(31)
               if (param(31) == 1) then
                  write (msg, '(A,I0,A,I0)') "Neklab single perturbation mode. lpert = ", lpert, ", npert = ", npert
               else
                  write (msg, '(A,I0,A,I0)') "Neklab multi-perturbation mode. lpert =", lpert, ", npert =", npert
               end if
               call nek_log_debug(msg, this_module, 'setup_nek')
               if (lpert /= npert) then
                  param(31) = lpert
                  npert = lpert
                  write (msg, '(A)') "Neklab requires lpert (SIZE) = npert (.par) to work reliably. Forcing npert=lpert."
                  call nek_log_message(msg, this_module, 'setup_nek')
               end if
      ! Deactivate OIFS.
               if (ifchar) then
                  write (msg, '(A)') "OIFS is not available for linearized solver. Turning it off."
                  call nek_log_warning(msg, this_module, 'setup_nek')
                  ifchar = .false.
               end if
            else
               if (ifpert) full_summary = .true.
               ifpert = .false.; call bcast(ifpert, lsize)
            end if
            if (solve_temperature_) then
               ifheat = .true.; call bcast(ifheat, lsize)
            end if
 
      ! Set integration time
            if (endtime_ <= 0.0_dp) then
               write (msg, *) 'Invalid endtime specified. Endtime =', endtime_
               call nek_stop_error(msg, this_module, 'setup_nek')
            end if
            param(10) = endtime_
            write (msg, '(A,F15.8)') padl('Set integration time: ', 30), param(10)
            if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
 
      ! Force CFL to chosen limit
            if (cfl_limit_ < 0.0_dp .or. cfl_limit_ > 0.5_dp) then
               write (msg, *) "Invalid target CFL. CLF_target =", cfl_limit_
               call nek_log_warning(msg, this_module, 'setup_nek')
               write (msg, *) "          Forcing it to", 0.5_dp
               call nek_log_warning(msg, this_module, 'setup_nek')
               cfl_limit_ = 0.5_dp
            end if
            param(26) = cfl_limit_
            write (msg, '(A,F15.8)') padl('Set CFL limit: ', 30), param(26)
            if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
 
      ! Force constant timestep if requested
            if (variable_dt_) then
               iffxdt = .false.
               ! always recompute dt and flush the internal dtold variable
               call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
               dt = param(26)/ctarg
               call compute_cfl(ctarg, vx, vy, vz, dt)
               if (dt <= 0.0_dp) then
                  call nek_stop_error('dt = 0 and a zero velocity field: no CFL to take dt from. '//
     &            'Set dt > 0 in the .par or start from a non-zero field.', this_module, 'setup_nek')
               end if
               param(12) = abs(dt)
               dtinit = param(12)
               ! flush dtold in subs1.f
               fintim = huge(1.0_dp) ! this avoids to jump to lastep = 1 irrespective of the value of time
               call setdt            ! resets dt based on cfl and a maximum increase of 20% compared to dtold ...
               idt = 0
               do while (abs(param(12) - dt) > 0.0_dp) ! we need to repeat this until dt_in = dt_out
                  idt = idt + 1
                  if (idt > max_idt) then
                     call nek_log_debug('max_idt reached.', this_module, 'setup_nek')
                     exit
                  end if
                  dt = param(12)
                  call setdt
               end do
               call nek_log_debug('Flushed dtold (subs1.f/setdt).', this_module, 'setup_nek')
               ! set correct fintim
               fintim = param(10)
               if (fintim <= 0.0_dp) then
                  write (msg, '(5X,A,F15.8)') padl('fintim (param(10)) = ', pad), fintim
                  call nek_log_message(msg, this_module, 'setup_nek')
                  write (msg, *) "When using variable dt fintim (param(10)) must be set to a valid value."
                  call nek_stop_error(msg, this_module, 'setup_nek')
               end if
               nsteps = 0
               if (LNS) then
                  write (msg, '(A,A)') padl('Set timestep: ', 30), 'read from baseflow'
                  if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
               else
                  write (msg, '(A,A)') padl('Set timestep: ', 30), 'variable dt'
                  if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
                  write (msg, '(A,ES15.8)') padl('Set fintim: ', 30), fintim
                  if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
               end if
            else
               iffxdt = .true.
               if (dt <= 0.0_dp) then
                  if (if_vel_zero) call nek_stop_error('dt = 0 and a zero velocity field: '//
     &               'no CFL to take dt from. Set dt > 0 in the .par or start from a non-zero field.',
     &               this_module, 'setup_nek')
                  recompute_dt_ = .true.
               end if
               ! Recompute dt
               if (recompute_dt_) then
                  dt_old = dt
                  call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
                  dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
                  dt = param(10)/nsteps; param(12) = dt
                  call compute_cfl(ctarg, vx, vy, vz, dt)
                  if (dt /= dt_old) then
                     write (msg, '(5X,A)') 'Recomputing dt/nsteps/cfl from target_cfl and current baseflow.'
                     call nek_log_information(msg, this_module, 'setup_nek')
                     write (msg, '(5X,A,F15.8)') padl('effective CFL = ', pad), ctarg
                     call nek_log_information(msg, this_module, 'setup_nek')
                  else
                     write (msg, '(5X,A)') 'Recomputing dt/nsteps/cfl from target_cfl and current baseflow.'
                     call nek_log_debug(msg, this_module, 'setup_nek')
                     write (msg, '(5X,A,F15.8)') padl('effective CFL = ', pad), ctarg
                     call nek_log_debug(msg, this_module, 'setup_nek')
                  end if
               else
      ! dt may have been set to endtime/nsteps by an earlier call: guard the
      ! ceiling against round-off, so that the step count does not change.
                  nsteps = ceiling(param(10)/dt*(1.0_dp - 1.0e-12_dp))
                  param(12) = dt
               end if
               if (nsteps < 1) then
                  write (msg, '(A,ES15.8,A,ES15.8)') 'invalid step count: dt = ', dt, ', endtime = ', param(10)
                  call nek_stop_error(msg, this_module, 'setup_nek')
               end if
               fintim = nsteps*dt
               param(12) = -abs(param(12))
               write (msg, '(A,ES15.8)') padl('Force constant timestep: ', 30), -param(12)
               if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
               write (msg, '(A,ES15.8)') padl('Set fintim: ', 30), fintim
               if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
               write (msg, '(A,I15)') padl('Set nsteps: ', 30), nsteps
               if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
            end if
 
      ! Set tolerances if requested
            param(21) = ptol_; TOLPDF = param(21); call bcast(TOLPDF,wdsize)
            param(22) = vtol_; TOLHDF = param(22); call bcast(TOLHDF,wdsize)
            restol(:) = param(22); call bcast(restol, (ldimt1+1)*wdsize)
            atol(:) = param(22); call bcast(atol, (ldimt1+1)*wdsize)
            write (msg, '(A,ES15.8)') padl('Set pressure tol: ', 30), param(21)
            if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
            write (msg, '(A,ES15.8)') padl('Set velocity tol: ', 30), param(22)
            if (nid == 0 .and. .not. silent_) print nekfmt, trim(msg)
 
      ! Broadcast parameters
            call bcast(param, 200*wdsize)
 
            if (nid == 0 .and. .not. silent_) then
               print *, ''
               print '("neklab ",A)', '############### SETUP COMPLETED ################'
               print *, ''
            end if
 
      ! Print status
            if (.not. silent_) call nek_status(full_summary)
 
         end subroutine setup_nek
 
         subroutine setup_nonlinear_solver(recompute_dt, variable_dt, endtime, vtol, ptol, cfl_limit, solve_temperature, silent)
            logical, optional, intent(in) :: recompute_dt
            logical, optional, intent(in) :: variable_dt
            real(dp), optional, intent(in) :: endtime
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            real(dp), optional, intent(in) :: cfl_limit
            logical, optional, intent(in) :: solve_temperature
            logical, optional, intent(in) :: silent
            call setup_nek(LNS=.false., recompute_dt=recompute_dt, variable_dt=variable_dt,
     &   endtime = endtime, vtol = vtol, ptol = ptol, cfl_limit = cfl_limit, solve_temperature=solve_temperature,
     &   silent = silent)
         end subroutine setup_nonlinear_solver
 
         subroutine setup_linear_solver(transpose, solve_baseflow, recompute_dt, variable_dt, endtime, vtol, ptol, cfl_limit, solve_temperature, silent)
            logical, optional, intent(in) :: transpose
            logical, optional, intent(in) :: solve_baseflow
            logical, optional, intent(in) :: recompute_dt
            logical, optional, intent(in) :: variable_dt
            real(dp), optional, intent(in) :: endtime
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            real(dp), optional, intent(in) :: cfl_limit
            logical, optional, intent(in) :: solve_temperature
            logical, optional, intent(in) :: silent
            call setup_nek(LNS=.true., transpose=transpose, solve_baseflow=solve_baseflow, recompute_dt=recompute_dt,
     &   variable_dt=variable_dt, endtime = endtime, vtol = vtol, ptol = ptol, cfl_limit = cfl_limit, solve_temperature=solve_temperature,
     &   silent = silent)
         end subroutine setup_linear_solver
 
         subroutine nek_status(full_summary)
            character(len=128) :: msg
            logical, optional, intent(in) :: full_summary
            logical :: full_summary_
            character(len=*), parameter :: nekfmt = '(5X,A)'
            integer, parameter :: pad = 20
            full_summary_ = optval(full_summary, .false.)
      ! overview
            if (nid == 0) then
               print *, ''
               print '("neklab ",A)', '################## NEK STATUS ##################'
               print *, ''
            end if
 
            if (ifpert) then
               if (full_summary_) then
                  call nek_log_message('')
                  call nek_log_message('LINEAR MODE:', this_module, 'nek_status', nekfmt)
                  write (msg, '(A,L15)') padl('ifpert: ', pad), ifpert
                  call nek_log_message(msg, this_module, 'nek_status', nekfmt)
                  write (msg, '(A,I15)') padl('npert: ', pad), npert
                  call nek_log_message(msg, this_module, 'nek_status', nekfmt)
                  if (ifadj) then
                     write (msg, '(A,L15)') padl('adjoint mode: ', pad), ifadj
                     call nek_log_message(msg, this_module, 'nek_status', nekfmt)
                  end if
                  if (ifbase) then
                     write (msg, '(A,L15)') padl('solve for baseflow: ', pad), ifbase
                     call nek_log_message(msg, this_module, 'nek_status', nekfmt)
                  end if
                  write (msg, '(A,L15)') padl('OIFS: ', pad), ifchar
                  call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               else
                  call nek_log_debug('LINEAR MODE', this_module, 'nek_status', nekfmt)
               end if
            else
               if (full_summary_) then
                  call nek_log_message('', this_module, 'nek_status', nekfmt)
                  call nek_log_message('NONLINEAR MODE:', this_module, 'nek_status', nekfmt)
                  write (msg, '(A,L15)') padl('OIFS: ', pad), ifchar
                  call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               else
                  call nek_log_debug('NONLINEAR MODE', this_module, 'nek_status', nekfmt)
               end if
            end if
            if (full_summary_) then
      ! params
               call nek_log_message('PARAMETERS:', this_module, 'nek_status', nekfmt)
               write (msg, '(A,F15.8)') padl('endtime: ', pad), param(10)
               call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               if (param(12) >= 0.0_dp) then
                  if (ifpert) then
                     write (msg, '(A,A)') padl('variable dt: ', pad), 'read from bflow'
                  else
                     write (msg, '(A,A)') padl('variable dt: ', pad), ' set at runtime'
                  end if
               else
                  write (msg, '(A,ES15.8)') padl('dt: ', pad), abs(param(12))
               end if
               call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               if (nsteps > 0) then
                  write (msg, '(A,I15)') padl('nsteps: ', pad), nsteps
                  call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               end if
               write (msg, '(A,F15.4)') padl('target CFL: ', pad), param(26)
               call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               write (msg, '(A,ES15.4)') padl('pressure tol: ', pad), param(21)
               call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               write (msg, '(A,ES15.4)') padl('velocity tol: ', pad), param(22)
               call nek_log_message(msg, this_module, 'nek_status', nekfmt)
               call nek_log_message('')
            end if
            if (nid == 0) then
               print *, ''
               print '("neklab ",A)', '################## NEK STATUS ##################'
               print *, ''
            end if
         end subroutine nek_status
 
         subroutine nek_log_message(msg, module, procedure, fmt)
            character(len=*), intent(in) :: msg
            character(len=*), optional, intent(in) :: module, procedure, fmt
            ! internal
            character(len=256) :: fmt_
            if (present(fmt)) then
               fmt_ = fmt
            else
               fmt_ = default_fmt('', module, procedure)
            end if
            call log_message(msg, module, procedure)          ! LK: I/O rank only
            if (nid == 0 .and. .not. logger_writes_to_stdout()) print fmt_, trim(msg)
         end subroutine nek_log_message
 
         subroutine nek_log_warning(msg, module, procedure, fmt)
            character(len=*), intent(in) :: msg
            character(len=*), optional, intent(in) :: module, procedure, fmt
            ! internal
            character(len=256) :: fmt_
            if (present(fmt)) then
               fmt_ = fmt
            else
               fmt_ = default_fmt('WARNING:', module, procedure)
            end if
            call log_warning(msg, module, procedure)
            if (nid == 0 .and. .not. logger_writes_to_stdout()) print fmt_, trim(msg)
         end subroutine nek_log_warning
 
         subroutine nek_log_debug(msg, module, procedure, fmt)
            character(len=*), intent(in) :: msg
            character(len=*), optional, intent(in) :: module
            character(len=*), optional, intent(in) :: procedure
            character(len=*), optional, intent(in) :: fmt
            ! internal
            character(len=128) :: fmt_
            if (present(fmt)) then
               fmt_ = fmt
            else
               fmt_ = default_fmt('DEBUG:', module, procedure)
            end if
            call log_debug(msg, module, procedure)
            if (nid == 0 .and. .not. logger_writes_to_stdout()) print fmt_, trim(msg)
         end subroutine nek_log_debug
 
         subroutine nek_log_information(msg, module, procedure, fmt)
            character(len=*), intent(in) :: msg
            character(len=*), optional, intent(in) :: module
            character(len=*), optional, intent(in) :: procedure
            character(len=*), optional, intent(in) :: fmt
            ! internal
            character(len=128) :: fmt_
            if (present(fmt)) then
               fmt_ = fmt
            else
               fmt_ = default_fmt('INFO:', module, procedure)
            end if
            call log_debug(msg, module, procedure)
            if (nid == 0 .and. .not. logger_writes_to_stdout()) print fmt_, trim(msg)
         end subroutine nek_log_information
 
         subroutine nek_stop_error(msg, module, procedure, fmt)
            character(len=*), intent(in) :: msg
            character(len=*), optional, intent(in) :: module
            character(len=*), optional, intent(in) :: procedure
            character(len=*), optional, intent(in) :: fmt
            call stop_error(msg, module=module, procedure=procedure)
         end subroutine nek_stop_error
 
         function default_fmt(prefix, module, procedure) result(fmt)
            character(len=*), optional, intent(in) :: prefix, module, procedure
            character(len=256) :: fmt
            character(len=200) :: s
            s = ''
            if (present(prefix)) s = trim(prefix)
            if (present(module)) then
               s = trim(s)//' '//trim(module)
               if (present(procedure)) s = trim(s)//' %'
            end if
            if (present(procedure)) s = trim(s)//' '//trim(procedure)//' :'
            write (fmt, '("(",A,",1X,A)")') '"'//trim(adjustl(s))//'"'
         end function default_fmt
 
         subroutine set_fldindex(file_prefix, fldindex) !file numbering suffix counter
            character(len=3) :: file_prefix
            integer :: fldindex
            ! internal
            integer :: iprefix
            integer :: nopen(99,2)
            common /RES_WANT/ nopen
            ! function
            integer :: i_find_prefix
      !     change prepost.f line 1094 from "save nopen" to "common /RES_WANT/ nopen"
            iprefix          = i_find_prefix(file_prefix,99)
            nopen(iprefix,1) = fldindex - 1
         end subroutine set_fldindex
      end module neklab_nek_setup
