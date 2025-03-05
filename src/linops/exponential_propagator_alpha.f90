      submodule(neklab_linops) exponential_propagator_alpha
         implicit none
      contains
         module procedure init_exptA_alpha
         ! internal
         integer :: idir, nelx, nely, nelz, i
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call logger%log_message("Set self%baseflow -> vx, vy, vz, pr, t", module=this_module, procedure="init_exptA")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow=.false., endtime=self%tau, recompute_dt=.true., cfl_limit=0.5_dp)
         idir = 1
         nelx = 10
         nely = 12
         nelz = 1
         call gtpp_gs_setup(self%hndl, nelx, nely, nelz, idir)
         do i = 1, lv
            self%cv(i) = cos(self%alpha*xm1(i, 1, 1, 1))
            self%sv(i) = sin(self%alpha*xm1(i, 1, 1, 1))
         end do
         end procedure
      
         module procedure exptA_matvec_alpha
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               call setup_linear_solver(transpose=.false., silent=.false.)
      ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Keep only chosen modes
               call self%proj()
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Keep only chosen modes
               call self%proj()
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call stop_error('Output must be a nek_dvector', module=this_module, procedure='exptA_matvec')
            end select
         class default
            call stop_error('Input must be a nek_dvector', module=this_module, procedure='exptA_matvec')
         end select
         end procedure
      
         module procedure exptA_rmatvec_alpha
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               call setup_linear_solver(transpose=.true., silent=.true.)
      ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Keep only chosen modes
               call self%proj()
               ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Keep only chosen modes
               call self%proj()
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call stop_error('Output must be a nek_dvector', module=this_module, procedure='exptA_rmatvec')
            end select
         class default
            call stop_error('Input must be a nek_dvector', module=this_module, procedure='exptA_rmatvec')
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
