      module neklab_nek_forcing
         use LightKrylov, only: dp
         use LightKrylov_Logger
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_forcing'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      
         public :: get_neklab_forcing, set_neklab_forcing
         public :: zero_neklab_forcing, zero_neklab_forcing_ipert
         public :: neklab_forcing
      
         real(dp), private, dimension(lv, lpert + 1) :: neklab_ffx = 0.0_dp
         real(dp), private, dimension(lv, lpert + 1) :: neklab_ffy = 0.0_dp
         real(dp), private, dimension(lv, lpert + 1) :: neklab_ffz = 0.0_dp
      
      contains
      
         subroutine zero_neklab_forcing()
      ! internal
            integer :: ipert
            character(len=128) :: msg
            write (msg, *) 'Setting all neklab_forcing vectors to zero.'
            call logger%log_debug(msg, module=this_module, procedure='zero_neklab_forcing')
            neklab_ffx = 0.0_dp
            neklab_ffy = 0.0_dp
            neklab_ffz = 0.0_dp
            return
         end subroutine zero_neklab_forcing
      
         subroutine get_neklab_forcing(fx, fy, fz, ipert)
            real(dp), dimension(lv), intent(out) :: fx
            real(dp), dimension(lv), intent(out) :: fy
            real(dp), dimension(lv), intent(out) :: fz
            integer, intent(in) :: ipert
      ! internal
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call logger%log_message(msg, module=this_module, procedure='get_neklab_forcing')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               write (msg, '(A,I0)') 'Retrieving value of the neklab forcing. ipert = ', ipert
               call logger%log_debug(msg, module=this_module, procedure='get_neklab_forcing')
            end if
            fx = neklab_ffx(:, ipert + 1)
            fy = neklab_ffy(:, ipert + 1)
            fz = neklab_ffz(:, ipert + 1)
            return
         end subroutine get_neklab_forcing
      
         subroutine set_neklab_forcing(fx, fy, fz, ipert)
            real(dp), dimension(lv), intent(in) :: fx
            real(dp), dimension(lv), intent(in) :: fy
            real(dp), dimension(lv), intent(in) :: fz
            integer, intent(in) :: ipert
      ! internal
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call logger%log_message(msg, module=this_module, procedure='set_neklab_forcing')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               write (msg, '(A,I0)') 'Setting value of the neklab forcing. ipert = ', ipert
               call logger%log_debug(msg, module=this_module, procedure='get_neklab_forcing')
            end if
            neklab_ffx(:, ipert + 1) = fx
            neklab_ffy(:, ipert + 1) = fy
            neklab_ffz(:, ipert + 1) = fz
            return
         end subroutine set_neklab_forcing
      
         subroutine zero_neklab_forcing_ipert(ipert)
            integer, intent(in) :: ipert
      ! internal
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call logger%log_message(msg, module=this_module, procedure='zero_neklab_forcing_ipert')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               write (msg, '(A,I0)') 'Setting value of the neklab forcing to zero. ipert = ', ipert
               call logger%log_debug(msg, module=this_module, procedure='zero_neklab_forcing_ipert')
            end if
            neklab_ffx(:, ipert + 1) = 0.0_dp
            neklab_ffy(:, ipert + 1) = 0.0_dp
            neklab_ffz(:, ipert + 1) = 0.0_dp
            return
         end subroutine zero_neklab_forcing_ipert
      
         subroutine neklab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg, ipert)
            real(dp), intent(inout) :: ffx
            real(dp), intent(inout) :: ffy
            real(dp), intent(inout) :: ffz
            integer, intent(in) :: ix
            integer, intent(in) :: iy
            integer, intent(in) :: iz
            integer, intent(in) :: ieg
            integer, intent(in) :: ipert
      ! internal
            integer :: e, ijke
            e = gllel(ieg)
            ijke = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
      ! note: ipert/jp == 0 : forcing for nonlinear solver
      ! note: ipert/jp >  0 : forcing for linear solver
            ffx = ffx + neklab_ffx(ijke, ipert + 1)
            ffy = ffy + neklab_ffy(ijke, ipert + 1)
            ffz = ffz + neklab_ffz(ijke, ipert + 1)
            return
         end subroutine neklab_forcing
      end module neklab_nek_forcing
