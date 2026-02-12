      submodule(neklab_systems) fixed_point
         implicit none
      contains
         module procedure nonlinear_map
         character(len=*), parameter :: this_procedure = 'nonlinear_map'
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

      ! Set the initial condition for the nonlinear solver.
               call vec2nek(vx, vy, vz, pr, t, vec_in)

      ! Set appropriate tolerances.
               call setup_nonlinear_solver(recompute_dt = .true., 
     &                                     cfl_limit    = 0.4_dp,
     &                                     vtol         = atol*0.1,
     &                                     ptol         = atol*0.1)      
      
      ! Integrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)

      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure nonlinear_map
      
         module procedure jac_exptA_matvec
         character(len=*), parameter :: this_procedure = 'jac_exptA_matvec'
         integer :: nrst
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)

      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Ensure correct nek status.
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)

      ! Set the initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Extract the final solution to vector.   
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Reset tolerances.
               param(22) = atol

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_exptA_matvec
      
         module procedure jac_exptA_rmatvec
         character(len=*), parameter :: this_procedure = 'jac_exptA_rmatvec'
         integer :: nrst
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)
            
      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
            
      ! Ensure correct nek status.
               call setup_linear_solver(transpose      = .true., 
     &                                  solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)

      ! Set the initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call self%get_rst(vec_in, istep)

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Reset tolerances.
               param(22) = atol

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_exptA_rmatvec

         module procedure jac_exptA_compute_rst
            character(len=*), parameter :: this_procedure = 'jac_exptA_compute_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            integer :: irst, itmp
            real(dp) :: rtmp
            select type(vec_out)
            type is (nek_dvector)
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, this_procedure)
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst
                  call nek_advance()
                  irst = istep - nsteps
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         end procedure

         module procedure jac_exptA_get_rst
            character(len=*), parameter :: this_procedure = 'jac_exptA_get_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            select type(vec_in)
            type is (nek_dvector)
               if (vec_in%has_rst_fields()) then
                  call vec_in%get_rst(vec_rst, istep)
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
               end if
            class default
               call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
            end select
         end procedure
      end submodule
