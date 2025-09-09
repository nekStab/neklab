      submodule(neklab_systems) periodic_orbit
         implicit none
      contains
         module procedure nonlinear_map_UPO
      ! internal
         real(dp) :: atol_v, atol_p
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)

      ! Set appropriate tolerances and Nek status.
               if (.not. self%nek_opts%is_initialized()) call self%nek_opts%init()
               write (msg, '(A,F9.6)') 'Current period estimate, T = ', vec_in%T
               call nek_log_message(msg, this_module, 'nonlinear_map_UPO')
               self%nek_opts%endtime = vec_in%T                                           ! set current period estimate
               self%nek_opts%vtol    = atol*self%tol_ratio_vel                            ! set velocity tolerance based on input value
               self%nek_opts%ptol    = atol*self%tol_ratio_pr*self%nek_opts%tol_ratio_pv  ! set pressure tolerance based on input value
               call set_nek_opts(self%nek_opts, stamp_log=.true., print_summary=.true.)

      ! Set the initial condition for the nonlinear solver.
               call ext_vec2nek(vx, vy, vz, pr, t, vec_in)

      ! Integrate the nonlinear equations forward.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

               end do

      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vx, vy, vz, pr, t)
               vec_out%T = vec_in%T

      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)

      ! Revert to input atol.
               self%nek_opts%vtol = atol
               self%nek_opts%ptol = atol*self%nek_opts%tol_ratio_pv
               call set_nek_opts(self%nek_opts)

            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module,'nonlinear_map_UPO')
            end select
         class default
            call type_error('vec_in','nek_ext_dvector','IN',this_module,'nonlinear_map_UPO')
         end select
         end procedure nonlinear_map_UPO
      
         module procedure jac_direct_map
      ! internal
         real(dp) :: atol_v, atol_p
         integer :: nrst, itmp, irst
         real(dp) :: rtmp
         type(nek_ext_dvector) :: vec_rst, vec
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)
               nrst = abs(param(27)) - 1

      ! Set baseflow.
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Set the initial condition for the linearized solver.
               call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Ensure correct nek status.
               atol_v = param(22)                                    ! save current velocity tolerance (set during nonlinear run)
               atol_p = param(21)                                    ! save current pressure tolerance (set during nonlinear run)
               self%nek_opts%vtol    = atol_v*self%tol_ratio_vel     ! set velocity tolerance based on current value
               self%nek_opts%ptol    = atol_p*self%tol_ratio_pr      ! set pressure tolerance based on current value
               self%nek_opts%endtime = get_period_abs(self%X)        ! set endtime
               call set_nek_opts(self%nek_opts, stamp_log=.true.)

      ! Intgrate the coupled equations forward.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst .and. vec_in%has_rst_fields()) then
                     call vec_in%get_rst(vec_rst, istep)
                     call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
                  end if

               end do

      ! Compute restart fields.
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, 'exptA_matvec')
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               self%nek_opts%endtime = time + nrst*dt
               call set_nek_opts(self%nek_opts)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst
                  
                  call nek_advance()
                  
                  irst = istep - nsteps
                  call nek2ext_vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time.
               istep = itmp
               time  = rtmp

      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Evaluate f'(X(T), T) * dT and add it to the position residual
      ! Here we assume that vx,vy,vz contains the endpoint of the nonlinear trajectory
               call compute_fdot(vec)
               call vec_out%axpby(vec_in%T, vec, 1.0_dp)

      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the nonlinear trajectory
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
               call compute_fdot(vec)
               vec_out%T = vec_in%dot(vec)

      ! Reset tolerances and endtime.
               self%nek_opts%vtol = atol_v
               self%nek_opts%ptol = atol_p
               self%nek_opts%endtime = time
               call set_nek_opts(self%nek_opts)

            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module,'jac_direct_map')
            end select
         class default
            call type_error('vec_in','nek_ext_dvector','IN',this_module,'jac_direct_map')
         end select
         end procedure jac_direct_map
      
         module procedure jac_adjoint_map
      ! internal
         real(dp) :: atol_v, atol_p
         integer :: nrst, itmp, irst
         real(dp) :: rtmp
         type(nek_ext_dvector) :: vec_rst, vec
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)

      ! Ensure correct nek status.
               atol_v = param(22)                                 ! save current velocity tolerance (set during nonlinear run)
               atol_p = param(21)                                 ! save current pressure tolerance (set during nonlinear run)
               self%nek_opts%vtol    = atol_v*self%tol_ratio_vel  ! set velocity tolerance based on current value
               self%nek_opts%ptol    = atol_p*self%tol_ratio_pr   ! set pressure tolerance based on current value
               self%nek_opts%endtime = get_period_abs(self%X)     ! set endtime
               call set_nek_opts(self%nek_opts, transpose = .true., stamp_log=.true.)

      ! Set baseflow.
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Set the initial condition for the linearized solver.
               call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst .and. vec_in%has_rst_fields()) then
                     call vec_in%get_rst(vec_rst, istep)
                     call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
                  end if

               end do

      ! Compute restart fields.
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, 'exptA_matvec')
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               self%nek_opts%endtime = time + nrst*dt
               call set_nek_opts(self%nek_opts, transpose= .true.)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst
                  
                  call nek_advance()
                  
                  irst = istep - nsteps
                  call nek2ext_vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time.
               istep = itmp
               time  = rtmp

      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Evaluate f'(X(T), T) * dT and add it to the position residual
               call compute_fdot(vec)
               call vec_out%axpby(vec_in%T, vec, 1.0_dp)

      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the orbit
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
               call compute_fdot(vec)
               vec_out%T = vec_in%dot(vec)

      ! Reset tolerances and endtime.
               self%nek_opts%vtol = atol_v
               self%nek_opts%ptol = atol_p
               self%nek_opts%endtime = time
               call set_nek_opts(self%nek_opts, transpose= .true.)
            
            class default
               call type_error('vec_out','nek_ext_dvector','OUT',this_module,'jac_adjoint_map')
            end select
         class default
            call type_error('vec_in','nek_ext_dvector','IN',this_module,'jac_adjoint_map')
         end select
         end procedure jac_adjoint_map
      end submodule
