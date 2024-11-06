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
         character(len=*), parameter :: this_module = 'neklab_utils'
          
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.

         ! utilities for regular nek vectors
         public :: nek2vec, vec2nek, abs_vec2nek, outpost_dnek
         ! utilities for extended nek vectors
         public :: nek2ext_vec, ext_vec2nek, abs_ext_vec2nek, outpost_ext_dnek
         ! Set up solver
         public :: setup_nek, setup_nonlinear_solver, setup_linear_solver, nek_status
         ! miscellaneous
         public :: nopcopy
      
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
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(in) :: pr_
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
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
      
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
            select type(vec)
            type is (nek_dvector)
               call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
            end select
            return
         end subroutine abstract_vec2nek_std
      
         subroutine abstract_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
            select type(vec)
            type is (nek_dvector)
               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
            end select
            return
         end subroutine abstract_vec2nek_prt


         ! EXTENDED


         subroutine nek2ext_vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
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
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
      
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
            select type(vec)
            type is (nek_ext_dvector)
               call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
            end select
            return
         end subroutine abstract_ext_vec2nek_std
      
         subroutine abstract_ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
            select type(vec)
            type is (nek_ext_dvector)
               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
            end select
            return
         end subroutine abstract_ext_vec2nek_prt
      
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

         subroutine setup_nek(if_LNS, if_adjoint, if_solve_baseflow, endtime, vtol, ptol, cfl_limit, full_summary, silent)
            logical,            intent(in) :: if_LNS
            logical, optional,  intent(in) :: if_adjoint
            logical                        :: if_adjoint_
            logical, optional,  intent(in) :: if_solve_baseflow
            logical                        :: if_solve_baseflow_
            real(dp), optional, intent(in) :: endtime
            real(dp)                       :: endtime_
            real(dp), optional, intent(in) :: vtol
            real(dp)                       :: vtol_
            real(dp), optional, intent(in) :: ptol
            real(dp)                       :: ptol_
            real(dp), optional, intent(in) :: cfl_limit
            real(dp)                       :: cfl_limit_
            logical, optional,  intent(in) :: full_summary
            logical                        :: full_summary_
            logical, optional,  intent(in) :: silent
            logical                        :: silent_
  
            if_adjoint_        = optval(if_adjoint,        .false.)
            if_solve_baseflow_ = optval(if_solve_baseflow, .false.)
            endtime_           = optval(endtime, param(10))
            ptol_              = optval(ptol, param(21))
            vtol_              = optval(vtol, param(22))
            cfl_limit_         = optval(cfl_limit, param(26))
            full_summary_      = optval(full_summary,      .false.)
            silent_            = optval(silent,            .false.)

            call nekgsync()

            ! general settings
            lastep = 0;

            ! Force contant time step.
            param(12) = -abs(param(12))
            
            if (if_LNS) then
               ifpert = .true.; call bcast(ifpert, lsize)
               if (if_adjoint_) then
                  ifadj = .true.; call bcast(ifadj, lsize)
               else
                  ifadj = .false.; call bcast(ifadj, lsize)
               end if
               if (if_solve_baseflow_) then
                  ifbase = .true.; call bcast(ifbase, lsize)
               else
                  ifbase = .false.; call bcast(ifbase, lsize)
               end if
               ! Force single perturbation mode.
               if (param(31) > 1) then
                  if (nid == 0) print *, "Neklab does not (yet) support npert > 1."
                  call nek_end()
               else
                  param(31) = 1; npert = 1
               end if
               ! Deactivate OIFS.
               if (ifchar) then
                  if (nid == 0) then
                     print *, "WARNING : OIFS is not available for linearized solver."
                     print *, "          Turning it off."
                  end if
                  ifchar = .false.
               end if
            else
               ifpert = .false.; call bcast(ifpert, lsize)
               param(31) = 0; npert = 0
            end if

            ! Force CFL to chosen limit
            if (cfl_limit_ < 0.0_dp .or. cfl_limit_ > 0.5_dp) then
               if (nid == 0) then
                  print *, "WARNING : Invalid target CFL. ", cfl_limit_
                  print *, "          Forcing it to", 0.5_dp
               end if
               cfl_limit_ = 0.5_dp
            end if
            param(26) = cfl_limit_

            ! Set integration time
            if (endtime_ <= 0.0_dp) then
               if (nid == 0) print *, 'Invalid endtime specified. Endtime =', endtime_
               call nek_end()
            end if
            param(10) = endtime_

            call compute_cfl(ctarg, vx, vy, vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(param(10)/dt)
            dt = param(10)/nsteps; param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt)
            if (nid == 0) print *, "INFO : Current CFL and target CFL :", ctarg, param(26)
            fintim = nsteps*dt

            ! Broadcast parameters
            call bcast(param, 200*wdsize)

            ! Print status
            if (.not. silent_) call nek_status(full_summary=full_summary_)

            return
         end subroutine setup_nek

         subroutine setup_nonlinear_solver(endtime, vtol, ptol, cfl_limit, full_summary, silent)
            real(dp), optional, intent(in) :: endtime
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            real(dp), optional, intent(in) :: cfl_limit
            logical, optional,  intent(in) :: full_summary
            logical, optional,  intent(in) :: silent
            call setup_nek(if_LNS=.false.,
     &                  endtime=endtime, vtol=vtol, ptol=ptol, cfl_limit=cfl_limit, full_summary=full_summary, silent=silent)
            return
         end subroutine setup_nonlinear_solver

         subroutine setup_linear_solver(if_adjoint, if_solve_baseflow, endtime, vtol, ptol, cfl_limit, full_summary, silent)
            logical, optional,  intent(in) :: if_adjoint
            logical, optional,  intent(in) :: if_solve_baseflow
            real(dp), optional, intent(in) :: endtime
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            real(dp), optional, intent(in) :: cfl_limit
            logical, optional,  intent(in) :: full_summary
            logical, optional,  intent(in) :: silent
            call setup_nek(if_LNS=.true., if_adjoint=if_adjoint, if_solve_baseflow=if_solve_baseflow,
     &                  endtime=endtime, vtol=vtol, ptol=ptol, cfl_limit=cfl_limit, full_summary=full_summary, silent=silent)
            return
         end subroutine setup_linear_solver

         subroutine nek_status(full_summary)
            character(len=:), allocatable :: msg
            logical, optional, intent(in) :: full_summary
            logical :: full_summary_
            full_summary_ = optval(full_summary, .false.)
            ! overview
            if (ifpert) then
               call logger%log_message('LINEAR MODE ACTIVE', module=this_module, procedure='nek_status')
               write(msg, '(A, L8)') padl('npert: ', 20), param(31)
               if (ifadj) then
                  write(msg, '(A, L8)') padl('adjoint mode: ', 20), ifadj
                  call logger%log_message(msg, module=this_module, procedure='nek_status')
               end if
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               if (ifbase) then
                  write(msg, '(A, L8)') padl('solve for baseflow: ', 20), ifbase
                  call logger%log_message(msg, module=this_module, procedure='nek_status')
               end if
            else
               call logger%log_message('  NONLINEAR MODE ACTIVE', module=this_module, procedure='nek_status')
            end if
            if (full_summary_) then
               write(msg, '(A, L8)') padl('OIFS: ', 20), ifchar
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               ! params
               call logger%log_message('PARAMETERS:', module=this_module, procedure='nek_status')
               write(msg, '(A, F15.8)') padl('endTime: ', 20), param(10)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               write(msg, '(A, E15.8)') padl('dt: ', 20), param(12)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               write(msg, '(A, I8)') padl('nsteps: ', 20), nsteps
               write(msg, '(A, F15.3)') padl('CFL: ', 20), param(26)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               write(msg, '(A, E15.8)') padl('pressure tol: ', 20), param(21)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
               write(msg, '(A, E15.8)') padl('velocity tol: ', 20), param(22)
               call logger%log_message(msg, module=this_module, procedure='nek_status')
            end if
         end subroutine nek_status
      
      end module neklab_utils