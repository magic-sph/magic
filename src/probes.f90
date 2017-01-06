module probe_mod

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

   use parallel_mod, only: rank
   use precision_mod
   use truncation, only: n_r_max, n_phi_max, nrp
   use radial_functions, only: r_cmb, orho1, or1, or2, r, r_icb
   use num_param, only: vScale
   use blocking, only: nThetaBs, sizeThetaB, nfs
   use horizontal_data, only: O_sin_theta, theta 
   use output_data, only: tag
   use constants, only: pi

   implicit none

   private

   real(cp) :: r_probe,theta_probe              !probe locations, r_probe in terms of r_cmb and theta in degrees
   integer  :: n_phi_probes,n_theta_usr,rad_usr !number of probes in phi - symmetrically distributed
   character(len=72) :: probe_file
    
   public   :: initialize_probes, probe_out 
   public   :: r_probe,theta_probe,n_phi_probes
 
contains

   subroutine initialize_probes
      
      integer           :: out_file_handle
      integer           :: rIter,thIter
      real(cp)          :: deg2rad

      deg2rad = pi/180.0_cp
      probe_file = "probeVp."//tag
      
      if(rank==0) then
         open(newunit=out_file_handle, file=probe_file, status='new')
         close(out_file_handle)
      end if
      
      r_probe = r_probe + r_icb
      
      theta_probe = mod(abs(theta_probe),180.0_cp)            !Make sure theta is positive and between 0 and 90
      if(theta_probe > 90.0_cp) theta_probe = 180 - theta_probe

      rad_usr = minloc(abs(r_probe - r),1)
      n_theta_usr = minloc(abs(theta_probe*deg2rad - theta),1)

   end subroutine initialize_probes

   subroutine probe_out(time,n_r,vp, n_theta_start,n_theta_block_size)
      
      real(cp), intent(in) :: time
      integer,  intent(in) :: n_r                    ! radial grod point no.
      integer,  intent(in) :: n_theta_start          ! start theta no.
      integer,  intent(in) :: n_theta_block_size     ! size of theta block
      real(cp), intent(in) :: vp(nrp,*)

      !-- Local variables:
      integer  :: n_theta       ! counter for colatitude
      integer  :: n_theta_cal   ! position of block colat in all colats
      integer  :: probe_phi_step
      integer  :: n_theta_probe
      integer  :: out_file_handle 
      logical  :: theta_found
      real(cp) :: fac,fac_r
      character(len=10) :: fmtstr !format string

      
      if ( n_r /= rad_usr ) return

      theta_found = .false.

      do n_theta=1,n_theta_block_size,2
         n_theta_cal=n_theta_start+n_theta-1
         if( n_theta_cal == n_theta_usr) then
            theta_found = .true.
            n_theta_probe = n_theta
            exit
         end if
      end do   

      if (.not. theta_found) return
      
      probe_phi_step = n_phi_max/n_phi_probes
      
      fac_r=or1(n_r)*vScale*orho1(n_r)
      fac=fac_r*O_sin_theta(n_theta_cal)
      
      write(fmtstr,'(i3)') 2*n_phi_probes       ! 2*n_phi_probes columns for data

      open(newunit=out_file_handle,file=probe_file,status='unknown', &
      &    position='append')

      write(out_file_handle,'(ES20.12,'//trim(fmtstr)//'ES16.8)')      &
      & time,fac*vp(1:n_phi_max:probe_phi_step,n_theta_probe), &
      &      fac*vp(1:n_phi_max:probe_phi_step,n_theta_probe+1)

      close(out_file_handle)

   end subroutine probe_out

end module probe_mod
