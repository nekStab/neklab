      module neklab_linops
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
      
      !------------------------------------------
      !-----     EXPONENTIAL PROPAGATOR     -----
      !------------------------------------------
      
         type, extends(abstract_linop_rdp), public :: exptA_linop
            real(kind=dp) :: tau
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: matvec => exptA_matvec
            procedure, pass(self), public :: rmatvec => exptA_rmatvec
            procedure, pass(self), public :: init => init_exptA
         end type exptA_linop
      
      contains
      
      !--------------------------------------------------
      !-----     TYPE-BOUND PROCEDURE FOR exptA     -----
      !--------------------------------------------------
      
         subroutine init_exptA(self)
            class(exptA_linop), intent(in) :: self
            call nekgsync()
      
      ! Setup Nek5000 logical flags for perturbative solver.
            ifpert = .true.; ifbase = .false.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
      
      ! Force single perturbation mode.
            if (param(31) > 1) then
               if (nid == 0) write (6, *) "neklab does not support (yet) npert > 1."
               call nek_end()
            else
               param(31) = 1; npert = 1
            end if
      
      ! Deactivate OIFS.
            if (ifchar) then
            if (nid == 0) then
               write (6, *) "WARNING : OIFS is not available for linearized solver."
               write (6, *) "          Turning it off."
            end if
            ifchar = .false.
            end if
      
      ! Force CFL to 0.5
            if (param(26) > 0.5) then
            if (nid == 0) then
               write (6, *) "WARNING : Target CFL is larger than 0.5"
               write (6, *) "          Forcing it to 0.5"
            end if
            param(26) = 0.5_dp
            end if
      
      ! Compute appropriate step size.
            param(10) = self%tau
            call compute_cfl(ctarg, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(self%tau/dt)
            dt = self%tau/nsteps; param(12) = dt
            call compute_cfl(ctarg, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0; fintim = nsteps*dt
      
      ! Force contant time step.
            param(12) = -abs(param(12))
      
      ! Broadcast parameters.
            call bcast(param, 200*wdsize)
      
            return
         end subroutine init_exptA
      
         subroutine exptA_matvec(self, vec_in, vec_out)
            class(exptA_linop), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
      ! Nek-related setup.
            ifadj = .false.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
            return
         end subroutine exptA_matvec
      
         subroutine exptA_rmatvec(self, vec_in, vec_out)
            class(exptA_linop), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Nek-related setup.
            ifadj = .true.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
            ifadj = .false.
            return
         end subroutine exptA_rmatvec
      
      end module neklab_linops
