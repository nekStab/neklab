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
      ! Neklab vectors
         use neklab_vectors
      
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
         public :: nek2ext_vec, ext_vec2nek, abs_ext_vec2nek, outpost_ext_dnek, get_period, get_period_abs
      ! Set up solver
         public :: setup_nek, setup_nonlinear_solver, setup_linear_solver, nek_status
      ! miscellaneous
         public :: nopcopy
      
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
      
      contains
      
         subroutine nek2vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
            return
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
      
            return
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
      
            return
         end subroutine vec2nek_std
      
         subroutine vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
      
            call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
            return
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
            end select
            return
         end subroutine abstract_vec2nek_std
      
         subroutine abstract_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
            select type (vec)
            type is (nek_dvector)
               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
            end select
            return
         end subroutine abstract_vec2nek_prt
      
      ! EXTENDED
      
         subroutine nek2ext_vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
            return
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
      
            return
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
      
            return
         end subroutine ext_vec2nek_std
      
         subroutine ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
      
            call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
            return
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
            end select
            return
         end subroutine abstract_ext_vec2nek_std
      
         subroutine abstract_ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
            select type (vec)
            type is (nek_ext_dvector)
               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
            end select
            return
         end subroutine abstract_ext_vec2nek_prt
      
         pure real(dp) function get_period_abs(vec) result(period)
            class(abstract_vector_rdp), intent(in) :: vec
            select type (vec)
            type is (nek_ext_dvector)
               period = vec%T
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
            return
         end subroutine nopcopy
      
         subroutine outpost_dnek_vector(vec, prefix)
            type(nek_dvector), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
            return
         end subroutine outpost_dnek_vector
      
         subroutine outpost_dnek_abs_vector(vec, prefix)
            class(abstract_vector_rdp), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            select type (vec)
            type is (nek_dvector)
               call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
            end select
            return
         end subroutine outpost_dnek_abs_vector
      
         subroutine outpost_dnek_basis(vec, prefix)
            type(nek_dvector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            integer :: i
            do i = 1, size(vec)
               call outpost_dnek_vector(vec(i), prefix)
            end do
            return
         end subroutine outpost_dnek_basis
      
         subroutine outpost_ext_dnek_vector(vec, prefix)
            type(nek_ext_dvector), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
            return
         end subroutine outpost_ext_dnek_vector
      
         subroutine outpost_ext_dnek_basis(vec, prefix)
            type(nek_ext_dvector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            integer :: i
            do i = 1, size(vec)
               call outpost_ext_dnek_vector(vec(i), prefix)
            end do
            return
         end subroutine outpost_ext_dnek_basis
      
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
      
      end module neklab_utils
