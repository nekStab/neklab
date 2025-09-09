      submodule(neklab_systems) fixed_point
         implicit none
      contains
         module procedure nonlinear_map
         ! internal
         real(dp) :: atol_v, atol_p
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

      ! Set the initial condition for the nonlinear solver.
               call vec2nek(vx, vy, vz, pr, t, vec_in)

      ! Set appropriate tolerances.
               if (.not. self%nek_opts%is_initialized()) call self%nek_opts%init()     ! read from param initially
               self%nek_opts%vtol = atol*self%tol_ratio_vel                            ! set velocity tolerance based on input value
               self%nek_opts%ptol = atol*self%tol_ratio_pr*self%nek_opts%tol_ratio_pv  ! set pressure tolerance based on input value
               call set_nek_opts(self%nek_opts, stamp_log=.true., print_summary=.true.)
      
      ! Integrate the nonlinear equations forward.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)

      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)

      ! Revert to input atol.
               self%nek_opts%vtol = atol
               self%nek_opts%ptol = atol*self%nek_opts%tol_ratio_pv
               call set_nek_opts(self%nek_opts)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module,'nonlinear_map')
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module,'nonlinear_map')
         end select
         end procedure nonlinear_map

         module procedure jac_exptA_matvec
      ! internal
         real(dp) :: atol_v, atol_p
         integer :: nrst, itmp, irst
         real(dp) :: rtmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               nrst = abs(param(27)) - 1
            
      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Set the initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Ensure correct nek status.
               atol_v = param(22)                                    ! save current velocity tolerance (set during nonlinear run)
               atol_p = param(21)                                    ! save current pressure tolerance (set during nonlinear run)
               self%nek_opts%vtol = atol_v*self%tol_ratio_vel        ! set velocity tolerance based on current value
               self%nek_opts%ptol = atol_p*self%tol_ratio_pr         ! set pressure tolerance based on current value
               call set_nek_opts(self%nek_opts, stamp_log=.true.) 

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst .and. vec_in%has_rst_fields()) then
                     call vec_in%get_rst(vec_rst, istep)
                     call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
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
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time.
               istep = itmp
               time  = rtmp

      ! Extract the final solution to vector.   
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Reset tolerances and endtime.
               self%nek_opts%vtol = atol_v
               self%nek_opts%ptol = atol_p
               self%nek_opts%endtime = time
               call set_nek_opts(self%nek_opts)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module,'jac_exptA_matvec')
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module,'jac_exptA_matvec')
         end select
         end procedure jac_exptA_matvec
      
         module procedure jac_exptA_rmatvec
      ! internal
         real(dp) :: atol_v, atol_p
         integer :: nrst, itmp, irst
         real(dp) :: rtmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               nrst = abs(param(27)) - 1
         
      ! Set baseflow.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Set the initial condition for the linearized solver.
      
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Ensure correct nek status.
               atol_v = param(22)                                    ! save current velocity tolerance (set during nonlinear run)
               atol_p = param(21)                                    ! save current pressure tolerance (set during nonlinear run)
               self%nek_opts%vtol = atol_v*self%tol_ratio_vel        ! set velocity tolerance based on current value
               self%nek_opts%ptol = atol_p*self%tol_ratio_pr         ! set pressure tolerance based on current value
               call set_nek_opts(self%nek_opts, transpose = .true., stamp_log=.true.)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  
                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst .and. vec_in%has_rst_fields()) then
                     call vec_in%get_rst(vec_rst, istep)
                     call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
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
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, irst)
               end do
               ! Reset iteration count and time.
               istep = itmp
               time  = rtmp

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

      ! Reset tolerances and endtime.
               self%nek_opts%vtol = atol_v
               self%nek_opts%ptol = atol_p
               self%nek_opts%endtime = time
               call set_nek_opts(self%nek_opts, transpose= .true.)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module,'jac_exptA_rmatvec')
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module,'jac_exptA_rmatvec')
         end select
         end procedure jac_exptA_rmatvec
      end submodule
