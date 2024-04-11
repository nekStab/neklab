      module neklab_linops
      use LightKrylov, only: abstract_linop, abstract_vector, wp
      use neklab_vectors
      implicit none
      include "SIZE"
      include "TOTAL"
      include "ADJOINT"
      
      private
      
      !------------------------------------------
      !-----                                -----
      !-----     EXPONENTIAL PROPAGATOR     -----
      !-----                                -----
      !------------------------------------------
      
      type, extends(abstract_linop), public :: exponential_propagator
      real(kind=wp) :: tau
      type(nek_dvector) :: baseflow
      contains
      private
      procedure, pass(self), public :: matvec => exptA_matvec
      procedure, pass(self), public :: rmatvec => exptA_rmatvec
      end type exponential_propagator
      
      contains
      
      !--------------------------------------------------
      !-----     TYPE-BOUND PROCEDURE FOR exptA     -----
      !--------------------------------------------------
      
      subroutine exptA_matvec(self, vec_in, vec_out)
      class(exponential_propagator), intent(in) :: self
      class(abstract_vector), intent(in) :: vec_in
      class(abstract_vector), intent(out) :: vec_out
      
      logical, save :: init
      data init/.false./
      
      if (.not. init) then
      ! Sets up the linearized solver.
      call setup_nek_linear_solver(self)
      ! Sets up the baseflow for Nek5000 solver.
      call nopcopy(vx, vy, vz, pr, t, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, self%baseflow%pr, self%baseflow%theta)
      ! Initialization flag .true.
      init = .true.
      endif
      
      ! Nek-related setup.
      lastep = 0 ; fintim = param(10)
      
      ! Sets the initial condition for Nek5000's linearized solver.
      select type(vec_in)
      type is(nek_dvector)
      call nopcopy(vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1), vec_in%vx, vec_in%vy, vec_in%vz, vec_in%pr, vec_in%theta)
      end select
      
      ! Integrate the equations forward in time.
      time = 0.0_wp
      do istep = 1, nsteps
      call nek_advance()
      enddo
      
      call outpost(vxp, vyp, vzp, prp, tp, "prt")
      
      ! Copy the final solution to vector.
      select type(vec_out)
      type is(nek_dvector)
      call nopcopy(vec_out%vx, vec_out%vy, vec_out%vz, vec_out%pr, vec_out%theta, vxp(:, 1), vyp(:, 1), vzp(:, 1), prp(:, 1), tp(:, :, 1))
      end select
      return
      contains
      subroutine setup_nek_linear_solver(linop)
      type(exponential_propagator), intent(in) :: linop
      
      call nekgsync()
      
      ! Setup Nek5000 logical flags for perturbative solver.
      ifpert = .true. ; ifadj = .false. ; ifbase = .false.
      call bcast(ifpert, lsize); call bcast(ifadj, lsize); call bcast(ifbase, lsize)
      
      ! Force single perturbation mode.
      if (param(31) > 1) then
      if (nid.eq.0) write(*, *) "neklab does not support (yet) npert > 1."
      call nek_end()
      else
      param(31) = 1; npert = 1
      endif
      
      ! Deactivate OIFS.
      if (ifchar) then
      if (nid.eq.0) then
      write(*, *) "WARNING : OIFS is not available for linearized solver."
      write(*, *) "          Turning it off."
      endif
      ifchar = .false.
      endif
      
      ! Force CFL to 0.5
      if (param(26) > 0.5) then
      if (nid.eq.0) then
      write(*, *) "WARNING : Target CFL is larger than 0.5"
      write(*, *) "          Forcing it to 0.5"
      endif
      param(26) = 0.5_wp
      endif
      
      ! Compute appropriate step size.
      param(10) = linop%tau
      call compute_cfl(ctarg, linop%baseflow%vx, linop%baseflow%vy, linop%baseflow%vz, 1.0_wp)
      dt = param(26)/ctarg ; nsteps = ceiling(linop%tau/dt)
      dt = linop%tau/nsteps ; param(12) = dt
      call compute_cfl(ctarg, linop%baseflow%vx, linop%baseflow%vy, linop%baseflow%vz, dt)
      if (nid.eq.0) write(*, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
      lastep = 0 ; fintim = nsteps*dt
      
      ! Force sontant time step.
      param(12) = -abs(param(12))
      
      ! Broadcast parameters.
      call bcast(param, 200*wdsize)
      return
      end subroutine setup_nek_linear_solver
      end subroutine exptA_matvec
      
      subroutine exptA_rmatvec(self, vec_in, vec_out)
      class(exponential_propagator), intent(in) :: self
      class(abstract_vector), intent(in) :: vec_in
      class(abstract_vector), intent(out) :: vec_out
      return
      end subroutine exptA_rmatvec
      
      end module neklab_linops
