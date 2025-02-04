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
      
      !--------------------------------------------------
      !-----     NEKLAB SYSTEM FOR FIXED-POINTS   -------
      !--------------------------------------------------
      
      ! --> Type: nek_system
         type, extends(abstract_system_rdp), public :: nek_system
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map
         end type nek_system
      
      ! --> Type: nek_jacobian
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian
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
      
      !-----------------------------------------------------
      !-----     NEKLAB SYSTEM FOR PERIODIC ORBITS   -------
      !-----------------------------------------------------
      
         type, extends(abstract_system_rdp), public :: nek_system_upo
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map_upo
         end type nek_system_upo
      
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_upo
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
      
      end module neklab_systems
