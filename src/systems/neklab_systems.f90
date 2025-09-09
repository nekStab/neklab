      module neklab_systems
      ! Fortran standard Library
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_optval, only: optval
         use stdlib_linalg, only: svd, diag
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp, qr
      ! Abstract types for real-valued vectors and utilities
         use LightKrylov, only: atol_dp
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: orthonormalize_basis, orthogonalize_against_basis, zero_basis, copy, rand_basis, linear_combination
         use LightKrylov, only: abstract_system_rdp, abstract_jacobian_linop_rdp
         use LightKrylov_Logger
         use LightKrylov_utils, only: assert_shape
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_systems'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      
         public :: compute_fdot
         public :: nek_constant_tol
         public :: nek_dynamic_tol
      
      !--------------------------------------------------
      !-----     NEKLAB SYSTEM FOR FIXED-POINTS   -------
      !--------------------------------------------------
      
      ! --> Type: nek_system
         type, extends(abstract_system_rdp), public :: nek_system
            type(nek_opts_std) :: nek_opts
            real(dp) :: tol_ratio_vel = 0.1_dp
            real(dp) :: tol_ratio_pr  = 0.1_dp
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map
         end type nek_system
      
      ! --> Type: nek_jacobian
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian
            type(nek_opts_prt) :: nek_opts
            real(dp) :: tol_ratio_vel = 0.1_dp
            real(dp) :: tol_ratio_pr  = 0.1_dp
         contains
            private
            procedure, pass(self), public :: matvec => jac_exptA_matvec
            procedure, pass(self), public :: rmatvec => jac_exptA_rmatvec
         end type nek_jacobian
      
      ! --> Type-bound procedures for nek_system & nek_jacobian
         interface
            module subroutine nonlinear_map(self, vec_in, vec_out, atol)
               class(nek_system), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
               real(dp), intent(in) :: atol
            end subroutine nonlinear_map

            module subroutine jac_exptA_matvec(self, vec_in, vec_out)
               class(nek_jacobian), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine jac_exptA_matvec
      
            module subroutine jac_exptA_rmatvec(self, vec_in, vec_out)
               class(nek_jacobian), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine jac_exptA_rmatvec
         end interface

      !--------------------------------------------------------------------
      !-----     NEKLAB SYSTEM FOR FIXED-POINTS with TEMP field     -------
      !--------------------------------------------------------------------
      
      ! --> Type: nek_system
         type, extends(abstract_system_rdp), public :: nek_system_temp
            type(nek_opts_std) :: nek_opts
            real(dp) :: tol_ratio_vel = 0.1_dp
            real(dp) :: tol_ratio_pr  = 0.1_dp
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map_temp
         end type nek_system_temp
      
      ! --> Type: nek_jacobian
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_temp
            type(nek_opts_prt) :: nek_opts
            real(dp) :: tol_ratio_vel = 0.1_dp
            real(dp) :: tol_ratio_pr  = 0.1_dp
         contains
            private
            procedure, pass(self), public :: matvec => jac_exptA_temp_matvec
            procedure, pass(self), public :: rmatvec => jac_exptA_temp_rmatvec
         end type nek_jacobian_temp
      
      ! --> Type-bound procedures for nek_system & nek_jacobian
         interface
            module subroutine nonlinear_map_temp(self, vec_in, vec_out, atol)
               class(nek_system_temp), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
               real(dp), intent(in) :: atol
            end subroutine nonlinear_map_temp

            module subroutine jac_exptA_temp_matvec(self, vec_in, vec_out)
               class(nek_jacobian_temp), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine jac_exptA_temp_matvec
      
            module subroutine jac_exptA_temp_rmatvec(self, vec_in, vec_out)
               class(nek_jacobian_temp), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine jac_exptA_temp_rmatvec
         end interface

      !-----------------------------------------------------
      !-----     NEKLAB SYSTEM FOR PERIODIC ORBITS   -------
      !-----------------------------------------------------
      
         type, extends(abstract_system_rdp), public :: nek_system_upo
            type(nek_opts_std) :: nek_opts
            real(dp) :: tol_ratio_vel = 0.1_dp
            real(dp) :: tol_ratio_pr  = 0.1_dp
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map_upo
         end type nek_system_upo
      
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_upo
            type(nek_opts_prt) :: nek_opts
            real(dp) :: tol_ratio_vel = 0.1_dp
            real(dp) :: tol_ratio_pr  = 0.1_dp
         contains
            private
            procedure, pass(self), public :: matvec => jac_direct_map
            procedure, pass(self), public :: rmatvec => jac_adjoint_map
         end type nek_jacobian_upo
      
      ! --> Type-bound procedures for nek_system_upo & nek_jacobian_upo
         interface
            module subroutine nonlinear_map_upo(self, vec_in, vec_out, atol)
               class(nek_system_upo), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
               real(dp), intent(in) :: atol
            end subroutine nonlinear_map_upo
      
            module subroutine jac_direct_map(self, vec_in, vec_out)
               class(nek_jacobian_upo), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine jac_direct_map
      
            module subroutine jac_adjoint_map(self, vec_in, vec_out)
               class(nek_jacobian_upo), intent(inout) :: self
               class(abstract_vector_rdp), intent(in) :: vec_in
               class(abstract_vector_rdp), intent(out) :: vec_out
            end subroutine jac_adjoint_map
         end interface
      
      contains
      
         subroutine compute_fdot(vec)
            class(nek_ext_dvector), intent(out) :: vec
      ! internal
            type(nek_ext_dvector) :: vec_in
      ! copy initial condition
            call nek2ext_vec(vec_in, vx, vy, vz, pr, t)
      ! Integrate the nonlinear equations forward in time.
            time = 0.0_dp
            do istep = 1, 1
               call nek_advance()
            end do
      ! Extract f(X(t0+dt))
            call nek2ext_vec(vec, vx, vy, vz, pr, t)
      ! Approximate derivative at t = t0:
      !
      !    f'(X(t0)) ~ ( f(X(t0+dt)) - f(X(t0)) ) / dt
      !
            call vec%sub(vec_in)
            call vec%scal(1.0/dt)
            vec%T = 0.0_dp ! ensure that the period shift vec_out%T is zero
            return
         end subroutine compute_fdot

         !---------------------------------------------------------------------
         !-----     Definition of two tolerance schedulers for Nek5000    -----
         !---------------------------------------------------------------------
      
         subroutine nek_constant_tol(tol, target_tol, rnorm, iter, info)
            !! Constant tolerance scheduler for the Newton iteration
            real(dp), intent(out) :: tol
            !! Tolerance to be used
            real(dp), intent(in) :: target_tol
            !! Target tolerance
            real(dp), intent(in)  :: rnorm
            !! Norm of the residual of the current iterate
            integer,  intent(in)  :: iter
            !! Newton iteration count
            integer,  intent(out)  :: info
            !! Information flag

            ! internals
            character(len=256) :: msg
            real(dp), parameter :: mintol = 10.0*atol_dp! minimum acceptable solver tolerance

            ! Sanity checks for the target tolerance
            if (target_tol < mintol) then
               ! Ensure target is above mintol
               tol = mintol
               write(msg,'(A,E11.4)') 'Input tolerance below minimum tolerance! Resetting solver tolerance to mintol= ', tol
               call nek_log_warning(msg, this_module, 'nek_constant_tol')
            else
               ! Ensure target is below maxtol
               tol = target_tol
               write(msg,'(A,E11.4)') 'Nek velocity and pressure tolerances set to tol= ', tol
               call nek_log_information(msg, this_module, 'nek_constant_tol')
            end if

            ! Update nek status (in case they have changed in between calls)
            param(21) = tol; TOLPDF = param(21); call bcast(TOLPDF,wdsize)
            param(22) = tol; TOLHDF = param(22); call bcast(TOLHDF,wdsize)
            restol(:) = param(22); call bcast(restol, (ldimt1+1)*wdsize)
            atol(:) = param(22); call bcast(atol, (ldimt1+1)*wdsize)
         end subroutine nek_constant_tol
      
         subroutine nek_dynamic_tol(tol, target_tol, rnorm, iter, info)
            !! Dynamic tolerance scheduler for the Newton iteration that sets the tolerance based on the current residual
            real(dp), intent(out) :: tol
            !! Tolerance to be used
            real(dp), intent(in) :: target_tol
            !! Target tolerance
            real(dp), intent(in)  :: rnorm
            !! Norm of the residual of the current iterate
            integer,  intent(in)  :: iter
            !! Newton iteration count
            integer,  intent(out)  :: info
            !! Information flag

            ! internals
            real(dp), parameter :: maxtol = 1.0e-04_dp   ! maximum acceptable solver tolerance
            real(dp), parameter :: mintol = 10.0*atol_dp ! minimum acceptable solver tolerance
            real(dp) :: tol_old, target_tol_
            character(len=256) :: msg

            ! Sanity checks for the target tolerance
            if (target_tol < mintol) then
               ! Ensure target is above mintol
               write(msg,'(A,E11.4)') 'Input target tolerance below minimum tolerance! Resetting target to mintol= ', mintol
               call nek_log_warning(msg, this_module, 'nek_dynamic_tol')
            end if
            target_tol_ = max(target_tol, mintol)
            if (target_tol > maxtol) then
               ! Ensure target is below maxtol
               write(msg,'(A,E11.4)') 'Input target tolerance above maximum tolerance! Resetting target to maxtol= ', maxtol
               call nek_log_warning(msg, this_module, 'nek_dynamic_tol')
            end if
            target_tol_ = min(target_tol_, maxtol)

            ! Compute new tolerance based on residual and target
            tol_old = tol
            tol = max(0.1_dp*rnorm, target_tol_)

            if (tol < 10*target_tol_) then
               ! If tolerance is close to target, make a single larger reduction step
               write(msg,'(A,E11.4)') 'Residual is close to target. Setting tolerance to input target= ', target_tol_
               call nek_log_information(msg, this_module, 'nek_dynamic_tol')
               tol = target_tol_
            end if
            if (tol > maxtol) then
               ! If tolerance is too large, cap it at maxtol
               write(msg,'(A,E11.4)') 'Residual is large. Setting tolerance to tol= ', maxtol
               call nek_log_information(msg, this_module, 'nek_dynamic_tol')
            end if
            tol = min(tol, maxtol)
      
            if (tol /= tol_old) then
               ! If tolerance is changed, stamp logs and update nek status
               if (tol == target_tol_) then
                  write(msg,'(A,E11.4)') 'Nek solver tolerance set to input target. tol= ', tol
               else
                  write(msg,'(A,E11.4)') 'Nek solver tolerance set to tol= ', tol
               end if
               call nek_log_information(msg, this_module, 'nek_dynamic_tol')
               param(21) = tol; TOLPDF = param(21); call bcast(TOLPDF,wdsize)
               param(22) = tol; TOLHDF = param(22); call bcast(TOLHDF,wdsize)
               restol(:) = param(22); call bcast(restol, (ldimt1+1)*wdsize)
               atol(:) = param(22); call bcast(atol, (ldimt1+1)*wdsize)
            else
               ! Do nothing
               write(msg,'(A,E11.4)') 'Nek solver tolerances unchanged at tol= ', tol_old
               call nek_log_information(msg, this_module, 'nek_dynamic_tol')
            end if
         end subroutine nek_dynamic_tol
                  
      end module neklab_systems
