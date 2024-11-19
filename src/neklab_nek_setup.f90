      module neklab_nek_setup
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
      ! Neklab vectors
         use neklab_vectors
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_utils'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         integer, parameter :: lt = lx1*ly1*lz1*lelt
      
      ! Set up solver
         public :: setup_nek, setup_nonlinear_solver, setup_linear_solver, nek_status
      
      contains
      
         subroutine setup_nek(LNS, transpose, solve_baseflow, recompute_dt, endtime, vtol, ptol, cfl_limit, silent)
            logical, intent(in) :: LNS
            logical, optional, intent(in) :: transpose
            logical :: transpose_
            logical, optional, intent(in) :: solve_baseflow
            logical :: solve_baseflow_
            logical, optional, intent(in) :: recompute_dt
            logical :: recompute_dt_
            real(dp), optional, intent(in) :: endtime
            real(dp) :: endtime_
            real(dp), optional, intent(in) :: vtol
            real(dp) :: vtol_
            real(dp), optional, intent(in) :: ptol
            real(dp) :: ptol_
            real(dp), optional, intent(in) :: cfl_limit
            real(dp) :: cfl_limit_
            logical, optional, intent(in) :: silent
            logical :: silent_
      ! internal
            character(len=128) :: msg
            logical :: full_summary
      
      ! Only print summary if we switch from linear to nonlinear solvers or vice versa
            full_summary = .false.
      
            transpose_ = optval(transpose, .false.)
            solve_baseflow_ = optval(solve_baseflow, .false.)
            endtime_ = optval(endtime, param(10))
            ptol_ = optval(ptol, param(21))
            vtol_ = optval(vtol, param(22))
            cfl_limit_ = optval(cfl_limit, param(26))
            silent_ = optval(silent, .false.)
      
            if (maxval(abs(vx)) == 0.0_dp .and. maxval(abs(vy)) == 0.0_dp .and. maxval(abs(vz)) == 0.0_dp) then
               recompute_dt_ = .false.
            else
               recompute_dt_ = optval(recompute_dt, .false.)
            end if
      
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
               if (param(31) > 1) then
                  write (msg, *) "Neklab does not (yet) support npert > 1."
                  call logger%log_message(msg, module=this_module, procedure='setup_nek')
                  if (nid == 0) print *, trim(msg)
                  call nek_end()
               else
                  param(31) = 1; npert = 1
               end if
      ! Deactivate OIFS.
               if (ifchar) then
                  write (msg, *) "OIFS is not available for linearized solver. Turning it off."
                  call logger%log_warning(msg, module=this_module, procedure='setup_nek')
                  if (nid == 0) print *, "WARNING :", trim(msg)
                  ifchar = .false.
               end if
            else
               if (ifpert) full_summary = .true.
               ifpert = .false.; call bcast(ifpert, lsize)
               param(31) = 0; npert = 0
            end if
      
      ! Set integration time
            if (endtime_ <= 0.0_dp) then
               write (msg, *) 'Invalid endtime specified. Endtime =', endtime_
               call logger%log_message(msg, module=this_module, procedure='setup_nek')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            end if
            param(10) = endtime_
            if (nid == 0 .and. .not. silent_) print '(5X,A)', 'Set integration time.'
      
      ! Force CFL to chosen limit
            if (cfl_limit_ < 0.0_dp .or. cfl_limit_ > 0.5_dp) then
               write (msg, *) "Invalid target CFL. CLF_target =", cfl_limit_
               call logger%log_warning(msg, module=this_module, procedure='setup_nek')
               if (nid == 0) print *, "WARNING :", trim(msg)
               write (msg, *) "          Forcing it to", 0.5_dp
               call logger%log_warning(msg, module=this_module, procedure='setup_nek')
               if (nid == 0) print *, trim(msg)
               cfl_limit_ = 0.5_dp
            end if
            param(26) = cfl_limit_
            if (nid == 0 .and. .not. silent_) print '(5X,A)', 'Set CFL limit.'
      
      ! Recompute dt
            if (recompute_dt_) then
               write (msg, '(5X,A)') 'Recomputing dt/nsteps/cfl from target_cfl and current baseflow.'
               call logger%log_information(msg, module=this_module, procedure='setup_nek')
               if (nid == 0) print '(A)', trim(msg)
               call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
               dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
               dt = param(10)/nsteps; param(12) = dt
               call compute_cfl(ctarg, vx, vy, vz, dt)
               write (msg, '(5X,A,F15.8)') padl('effective CFL = ', 20), ctarg
               call logger%log_information(msg, module=this_module, procedure='setup_nek')
               if (nid == 0) print '(A)', trim(msg)
            else
               nsteps = ceiling(param(10)/dt)
               param(12) = dt
            end if
            fintim = nsteps*dt
      
      ! Set tolerances if requested
            param(21) = ptol_
            param(22) = vtol_
            if (nid == 0 .and. .not. silent_) print '(5X,A)', 'Set velocity and pressure solver tolerances.'
      
      ! Force constant timestep
            param(12) = -abs(param(12))
            if (nid == 0 .and. .not. silent_) print '(5X,A)', 'Force a constant timestep.'
      
      ! Broadcast parameters
            call bcast(param, 200*wdsize)
      
            if (nid == 0 .and. .not. silent_) then
               print *, ''
               print '("neklab ",A)', '############### SETUP COMPLETED ################'
               print *, ''
            end if
      
      ! Print status
            if (.not. silent_) call nek_status(full_summary)
      
            return
         end subroutine setup_nek
      
         subroutine setup_nonlinear_solver(recompute_dt, endtime, vtol, ptol, cfl_limit, silent)
            logical, optional, intent(in) :: recompute_dt
            real(dp), optional, intent(in) :: endtime
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            real(dp), optional, intent(in) :: cfl_limit
            logical, optional, intent(in) :: silent
            call setup_nek(LNS=.false., recompute_dt=recompute_dt,
     $   endtime = endtime, vtol = vtol, ptol = ptol, cfl_limit = cfl_limit, silent = silent)
            return
         end subroutine setup_nonlinear_solver
      
         subroutine setup_linear_solver(transpose, solve_baseflow, recompute_dt, endtime, vtol, ptol, cfl_limit, silent)
            logical, optional, intent(in) :: transpose
            logical, optional, intent(in) :: solve_baseflow
            logical, optional, intent(in) :: recompute_dt
            real(dp), optional, intent(in) :: endtime
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            real(dp), optional, intent(in) :: cfl_limit
            logical, optional, intent(in) :: silent
            call setup_nek(LNS=.true., transpose=transpose, solve_baseflow=solve_baseflow, recompute_dt=recompute_dt,
     $   endtime = endtime, vtol = vtol, ptol = ptol, cfl_limit = cfl_limit, silent = silent)
            return
         end subroutine setup_linear_solver
      
         subroutine nek_status(full_summary)
            character(len=128) :: msg
            logical, optional, intent(in) :: full_summary
            logical :: full_summary_
            character(len=*), parameter :: nekfmt = '(5X,A)'
            full_summary_ = optval(full_summary, .false.)
      ! overview
            if (nid == 0) then
               print *, ''
               print '("neklab ",A)', '################## NEK STATUS ##################'
               print *, ''
            end if
      
            if (ifpert) then
            if (full_summary_) then
               call logger%log_message('LINEAR MODE:', module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, 'LINEAR MODE:'
               write (msg, '(A,L8)') padl('ifpert: ', 20), ifpert
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               write (msg, '(A,I8)') padl('npert: ', 20), npert
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               if (ifadj) then
                  write (msg, '(A,L8)') padl('adjoint mode: ', 20), ifadj
                  call logger%log_message(msg, module=this_module, procedure='nek_status')
                  if (nid == 0) print nekfmt, trim(msg)
               end if
               if (ifbase) then
                  write (msg, '(A,L8)') padl('solve for baseflow: ', 20), ifbase
                  call logger%log_message(msg, module=this_module, procedure='nek_status')
                  if (nid == 0) print nekfmt, trim(msg)
               end if
               write (msg, '(A,L8)') padl('OIFS: ', 20), ifchar
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
            else
               write (msg, '(A,A,L8,A,I8)') 'LINEAR MODE: ', padl('ifpert: ', 10), ifpert, padl('npert: ', 10), npert
               call logger%log_message('NEK5000 '//msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
            end if
            else
            if (full_summary_) then
               call logger%log_message('NONLINEAR MODE:', module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, 'NONLINEAR MODE:'
               write (msg, '(A,L8)') padl('OIFS: ', 20), ifchar
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
            else
               write (msg, '(A)') 'NONLINEAR MODE'
               call logger%log_message('NEK5000 '//msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
            end if
            end if
            if (full_summary_) then
      ! params
               call logger%log_message('PARAMETERS:', module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, 'PARAMETERS:'
               write (msg, '(A,F15.8)') padl('endtime: ', 20), param(10)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               write (msg, '(A,E15.8)') padl('dt: ', 20), abs(param(12))
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               write (msg, '(A,I8)') padl('nsteps: ', 20), nsteps
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               write (msg, '(A,F15.8)') padl('Target CFL: ', 20), param(26)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               write (msg, '(A,E15.8)') padl('pressure tol: ', 20), param(21)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
               write (msg, '(A,E15.8)') padl('velocity tol: ', 20), param(22)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (nid == 0) print nekfmt, trim(msg)
            end if
            if (nid == 0) then
               print *, ''
               print '("neklab ",A)', '################## NEK STATUS ##################'
               print *, ''
            end if
            return
         end subroutine nek_status
      
      end module neklab_nek_setup
