module probe_mod
   !
   ! Module for artificial sensors to compare time series
   ! of physical data with experiments.
   ! Probes are located in a radially symmetrical way
   ! on a radial surface given by r_probe (in terms of r_cmb),
   ! theta_probe in degrees between 0 and 90 and n_phi_probes
   ! denoting the number of probes in phi.
   ! Probes will be located at 'n_phi_probes' points
   ! at two equatorially symmetric latitudes - theta_probe and
   ! (180 - theta_probe) on r = r_probe.
   !
   ! version 1.0: Works only for v_phi, for now. Will be extended for other data later.

   use precision_mod
   use parallel_mod, only: rank
   use truncation, only: n_r_max, n_phi_max, n_theta_max
   use radial_data, only: nRstart, nRstop
   use radial_functions, only: r_cmb, orho1, or1, or2, r
   use num_param, only: vScale
   use horizontal_data, only: O_sin_theta, theta_ord, n_theta_ord2cal
   use output_data, only: tag
   use constants, only: pi, one
   use logic, only: l_save_out

   implicit none

   private

   real(cp), public :: r_probe,theta_probe  !probe locations, r_probe in terms of r_cmb and theta in degrees
   integer, public  :: n_phi_probes         !number of probes in phi - symmetrically distributed
   integer  :: n_theta_usr, rad_usr, rad_rank
   integer  :: n_probeVp, n_probeBr, n_probeBt
   character(len=72) :: probe_fileVp, probe_fileBr, probe_fileBt

   public   :: initialize_probes, finalize_probes, probe_out

contains

   subroutine initialize_probes

      real(cp) :: deg2rad

      deg2rad = pi/180.0_cp
      probe_fileVp = "probeVp."//tag
      probe_fileBr = "probeBr."//tag
      probe_fileBt = "probeBt."//tag

      r_probe = r_probe * r_cmb

      theta_probe = mod(abs(theta_probe),180.0_cp)  ! Make sure theta is positive and between 0 and 90
      if (theta_probe > 90.0_cp) theta_probe = 180.0_cp - theta_probe

      rad_usr = minloc(abs(r_probe - r),1)

      if ((nRstart <= rad_usr) .and. (rad_usr <= nRstop)) then
         if ( .not. l_save_out ) then
            open(newunit=n_probeVp, file=probe_fileVp, status='new')
            open(newunit=n_probeBr, file=probe_fileBr, status='new')
            open(newunit=n_probeBt, file=probe_fileBt, status='new')
         end if
         rad_rank = rank
      end if

      n_theta_usr = minloc(abs(theta_probe*deg2rad - theta_ord),1)
      n_theta_usr = n_theta_ord2cal(n_theta_usr)

   end subroutine initialize_probes
!-------------------------------------------------------------------------------
   subroutine finalize_probes

      if ( rank == rad_rank .and. (.not. l_save_out) ) then
         close(n_probeVp)
         close(n_probeBr)
         close(n_probeBt)
      end if

   end subroutine finalize_probes
!-------------------------------------------------------------------------------
   subroutine probe_out(time,n_r,vp,br,bt)

      real(cp), intent(in) :: time ! Time
      integer,  intent(in) :: n_r ! radial grod point no.
      real(cp), intent(in) :: vp(:,:), br(:,:), bt(:,:)

      !-- Local variables:
      integer  :: n_theta       ! counter for colatitude
      integer  :: probe_phi_step
      integer  :: n_theta_probe
      logical  :: theta_found
      real(cp) :: fac,fac_r
      character(len=10) :: fmtstr !format string

      if ( n_r /= rad_usr ) return

      theta_found = .false.

      do n_theta=1,n_theta_max,2
         if( n_theta == n_theta_usr) then
            theta_found = .true.
            n_theta_probe = n_theta
            exit
         end if
      end do

      if (.not. theta_found) return

      probe_phi_step = n_phi_max/n_phi_probes

      write(fmtstr,'(i3)') 2*n_phi_probes       ! 2*n_phi_probes columns for data

      if ( rank == rad_rank ) then
         if ( l_save_out ) then
            open(newunit=n_probeVp,file=probe_fileVp,status='unknown', &
            &    position='append')
            open(newunit=n_probeBr,file=probe_fileBr,status='unknown', &
            &    position='append')
            open(newunit=n_probeBt,file=probe_fileBt,status='unknown', &
            &    position='append')
         end if

         !-- Vp
         fac_r=or1(n_r)*vScale*orho1(n_r)
         fac=fac_r*O_sin_theta(n_theta)
         write(n_probeVp,'(ES20.12,'//trim(fmtstr)//'ES16.8)')     &
         & time,fac*vp(n_theta_probe,1:n_phi_max:probe_phi_step),  &
         &      fac*vp(n_theta_probe+1,1:n_phi_max:probe_phi_step)

         !-- Br
         fac=or2(n_r)
         write(n_probeBr,'(ES20.12,'//trim(fmtstr)//'ES16.8)')     &
         & time,fac*br(n_theta_probe,1:n_phi_max:probe_phi_step),  &
         &      fac*br(n_theta_probe+1,1:n_phi_max:probe_phi_step)

         !-- Btheta
         fac=or1(n_r)*O_sin_theta(n_theta)
         write(n_probeBt,'(ES20.12,'//trim(fmtstr)//'ES16.8)')     &
         & time,fac*bt(n_theta_probe,1:n_phi_max:probe_phi_step),  &
         &      fac*bt(n_theta_probe+1,1:n_phi_max:probe_phi_step)

         if ( l_save_out ) then
            close(n_probeVp)
            close(n_probeBr)
            close(n_probeBt)
         end if
      end if

   end subroutine probe_out
!-------------------------------------------------------------------------------
end module probe_mod
