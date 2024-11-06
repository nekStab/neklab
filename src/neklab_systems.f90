      module neklab_systems
      ! Fortran standard Library
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_optval, only: optval
         use stdlib_linalg, only: svd, diag
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp, qr, apply_inverse_permutation_matrix
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
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_systems'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
      !-------------------------------------------
      !-----     LIGHTKRYLOV SYSTEM TYPE   -------
      !-------------------------------------------
      
         type, extends(abstract_system_rdp), public :: nek_system
         contains
            private
            procedure, pass(self), public :: init => init_DNS
            procedure, pass(self), public :: eval => nonlinear_map
         end type nek_system
      
         type, extends(abstract_system_rdp), public :: nek_system_upo
         contains
            private
            procedure, pass(self), public :: init => init_DNS_upo
            procedure, pass(self), public :: eval => nonlinear_map_upo
         end type nek_system_upo
      
      !---------------------------------------------
      !-----     LIGHTKRYLOV JACOBIAN TYPE     -----
      !---------------------------------------------
      
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian
         contains
            private
            procedure, pass(self), public :: init => init_jac_exptA
            procedure, pass(self), public :: matvec => jac_exptA_matvec
            procedure, pass(self), public :: rmatvec => jac_exptA_rmatvec
         end type nek_jacobian
      
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_upo
         contains
            private
            procedure, pass(self), public :: init => init_jac_map
            procedure, pass(self), public :: matvec => jac_direct_map
            procedure, pass(self), public :: rmatvec => jac_adjoint_map
         end type nek_jacobian_upo
      
      contains
      
      !-------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM    -----
      !-------------------------------------------------------
      
         subroutine init_DNS(self)
            class(nek_system), intent(in) :: self
      
      ! We assume the baseflow has been set already!
      
      ! Setup Nek5000 for perturbation solver
            call setup_nonlinear_solver(recompute_dt=.true., cfl_limit=0.4_dp, full_summary=.true.)
      
            return
         end subroutine init_DNS
      
         subroutine nonlinear_map(self, vec_in, vec_out, atol)
      ! Dynamical system.
            class(nek_system), intent(in) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Solver tolerances if needed
            real(dp), intent(in) :: atol
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Set appropriate tolerances
                  call setup_nonlinear_solver(vtol=atol/100, ptol=atol/100)
      
      ! Ensure correct nek status and reset dt based on new baseflow
                  call self%init()
      
      ! Set the initial condition
                  call vec2nek(vx, vy, vz, pr, t, vec_in)
      
      ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Extract the final solution to vector.
                  call nek2vec(vec_out, vx, vy, vz, pr, t)
      
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine nonlinear_map
      
      !---------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN    -----
      !---------------------------------------------------------
      
         subroutine init_jac_exptA(self)
            class(nek_jacobian), intent(in) :: self
      
      ! Force the baseflow field for dt/nsteps/clf computation
            call abs_vec2nek(vx, vy, vz, pr, t, self%X)
            call logger%log_message('Set self%X -> vx, vy, vz, pr, t', module=this_module, procedure='init_jac_exptA')
      
      ! Setup Nek5000 for perturbation solver
            call setup_linear_solver(if_solve_baseflow=.false., recompute_dt=.true., cfl_limit=0.5_dp, full_summary=.true.)
      
            return
         end subroutine init_jac_exptA
      
         subroutine jac_exptA_matvec(self, vec_in, vec_out)
            class(nek_jacobian), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Ensure correct nek status
                  call setup_linear_solver(if_adjoint=.false., silent=.true.)
      
      ! Set the baseflow initial condition
                  call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Set the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Extract the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine jac_exptA_matvec
      
         subroutine jac_exptA_rmatvec(self, vec_in, vec_out)
            class(nek_jacobian), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Ensure correct nek status
                  call setup_linear_solver(if_adjoint=.true., silent=.true.)
      
      ! Set the baseflow initial condition
                  call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Set the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Extract the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            ifadj = .false.
            return
         end subroutine jac_exptA_rmatvec
      
      !-----------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM_UPO    -----
      !-----------------------------------------------------------
      
         subroutine init_DNS_UPO(self)
            class(nek_system_upo), intent(in) :: self
      
      ! We assume the baseflow has been set already!
      
      ! Setup Nek5000 for perturbation solver
            call setup_nonlinear_solver(recompute_dt=.true., cfl_limit=0.4_dp, full_summary=.true.)
      
            return
         end subroutine init_DNS_UPO
      
         subroutine nonlinear_map_UPO(self, vec_in, vec_out, atol)
      ! Dynamical system.
            class(nek_system_upo), intent(in) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Solver tolerances if needed
            real(dp), intent(in) :: atol
      
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
      ! Set appropriate tolerances
                  call setup_nonlinear_solver(vtol=atol/100, ptol=atol/100)
      
      ! Ensure correct nek status and reset dt based on new baseflow
                  call self%init()
      
      ! Set the initial condition
                  call ext_vec2nek(vx, vy, vz, pr, t, vec_in)
      
      ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vx, vy, vz, pr, t)
      
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine nonlinear_map_UPO
      
      !-------------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN_UPO    -----
      !-------------------------------------------------------------
      
         subroutine init_jac_map(self)
            class(nek_jacobian_upo), intent(in) :: self
      ! Force the initial condition for dt/nsteps/clf computation
            call abs_vec2nek(vx, vy, vz, pr, t, self%X)
            call logger%log_message('Set self%X -> vx, vy, vz, pr, t', module=this_module, procedure='init_jac_map')
      
      ! Setup Nek5000 for perturbation solver
            call setup_linear_solver(if_solve_baseflow=.true., recompute_dt=.true.,
     $   cfl_limit = 0.4_dp, full_summary = .true.)
      
            return
         end subroutine init_jac_map
      
         subroutine jac_direct_map(self, vec_in, vec_out)
      ! Dynamical system.
            class(nek_jacobian_upo), intent(in) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_ext_dvector) :: vec
      
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
      ! Set the baseflow initial condition
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(if_adjoint=.false., endtime=vec_in%T, silent=.true.)
      
      ! Set the perturbation initial condition
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
      
      ! Evaluate f'(X(T), T) * dT and add it to the position residual
                  call compute_fdot(vec)
                  call vec_out%axpby(1.0_dp, vec, vec_in%T)
      
      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Get the initial point of the orbit
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
                  call compute_fdot(vec)
                  vec_out%T = vec_in%dot(vec)
               end select
            end select
      
            return
         end subroutine jac_direct_map
      
         subroutine jac_adjoint_map(self, vec_in, vec_out)
            class(nek_jacobian_upo), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_ext_dvector) :: vec
      
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
      ! Set the baseflow initial condition
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(if_adjoint=.true., endtime=vec_in%T, silent=.true.)
      
      ! Set the perturbation initial condition
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
      
      ! Evaluate f'(X(T), T) * dT and add it to the position residual
                  param(10) = vec_out%T
                  call compute_fdot(vec)
                  call vec_out%axpby(1.0_dp, vec, vec_in%T)
      
      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the orbit
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
                  call compute_fdot(vec)
                  vec_out%T = vec_in%dot(vec)
               end select
            end select
      
            ifadj = .false.
            return
         end subroutine jac_adjoint_map
      
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
