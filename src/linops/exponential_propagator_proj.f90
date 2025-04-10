      submodule(neklab_linops) exponential_propagator_proj
         implicit none
      contains
         module procedure init_exptA_proj
      ! internal
         integer :: idir, nelx, nely, nelz, i
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_information("Set self%baseflow -> vx, vy, vz, pr, t", this_module, "init_exptA_proj")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         
      ! Define and initialize planar average
         idir = 1
         nelx = 10
         nely = 12
         nelz = 1
         call gtpp_gs_setup(self%hndl, nelx, nely, nelz, idir)
      ! Define projection basis
         do i = 1, lv
            self%cv(i) = cos(self%alpha*xm1(i, 1, 1, 1))
            self%sv(i) = sin(self%alpha*xm1(i, 1, 1, 1))
         end do
         end procedure
      
         module procedure exptA_proj_matvec
         integer :: nrst, itmp
         real(dp) :: rtmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               call setup_linear_solver(transpose     = .false.,
     &                                  silent        = .true.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Project out unwanted wavenumbers
               call self%proj()
      
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
      
      ! Project out unwanted wavenumbers
               call self%proj()
      
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Compute restart fields.
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, 'exptA_proj_matvec')
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst

                  call nek_advance()

                  call self%proj()
                  
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module,'exptA_proj_matvec')
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module,'exptA_proj_matvec')
         end select
         end procedure
      
         module procedure exptA_proj_rmatvec
         integer :: nrst, itmp
         real(dp) :: rtmp
         type(nek_dvector) :: vec_rst
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               
               nrst = abs(param(27)) - 1
               call setup_linear_solver(transpose     = .true.,
     &                                  silent        = .true.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Project out unwanted wavenumbers
               call self%proj()
      
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
      
      ! Project out unwanted wavenumbers
               call self%proj()

      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Compute restart fields.
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, 'exptA_rmatvec_proj')
               ! We don't need to reset the end time but we do it to get a clean logfile
               itmp = nsteps
               rtmp = time
               call setup_linear_solver(endtime = time + nrst*dt)
               nsteps = itmp
               do istep = nsteps + 1, nsteps + nrst

                  call nek_advance()

                  call self%proj()

                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst)
               end do
               ! Reset iteration count and time
               istep = itmp
               time  = rtmp

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module,'exptA_proj_rmatvec')
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module,'exptA_proj_rmatvec')
         end select
         end procedure

         module procedure proj_alpha
            ! internal
            integer :: ntot
            real(dp), dimension(lx1,ly1,lz1,lelv) :: ip
            real(dp), dimension(lx1,ly1,lz1,lelv) :: utmp
            real(dp), dimension(lx1,ly1,lz1,lelv) :: avg
            
            ntot = lx1*ly1*lz1*nelv
            ! x
            call rzero(utmp, ntot)
            call rzero(ip, ntot); call admcol3(ip,vxp,self%cv,2.0_dp,ntot) ! compute integrand
            call planar_avg(avg,ip,self%hndl)                              ! streamwise integral
            call add2col2(utmp,self%cv,avg,ntot)                           ! reconstruct mode
            call rzero(ip, ntot); call admcol3(ip,vxp,self%sv,2.0_dp,ntot)
            call planar_avg(avg,ip,self%hndl)
            call add2col2(utmp,self%sv,avg,ntot)
            call copy(vxp, utmp, ntot)
            ! y
            call rzero(utmp, ntot)
            call rzero(ip, ntot); call admcol3(ip,vyp,self%cv,2.0_dp,ntot)
            call planar_avg(avg,ip,self%hndl)
            call add2col2(utmp,self%cv,avg,ntot)
            call rzero(ip, ntot); call admcol3(ip,vyp,self%sv,2.0_dp,ntot)
            call planar_avg(avg,ip,self%hndl)
            call add2col2(utmp,self%sv,avg,ntot)
            call copy(vyp, utmp, ntot)
            ! z
            if (if3d) then
               call rzero(utmp, ntot)
               call rzero(ip, ntot); call admcol3(ip,vzp,self%cv,2.0_dp,ntot)
               call planar_avg(avg,ip,self%hndl)
               call add2col2(utmp,self%cv,avg,ntot)
               call rzero(ip, ntot); call admcol3(ip,vzp,self%sv,2.0_dp,ntot)
               call planar_avg(avg,ip,self%hndl)
               call add2col2(utmp,self%sv,avg,ntot)
               call copy(vzp, utmp, ntot)
            endif
      
         end procedure proj_alpha
      end submodule
