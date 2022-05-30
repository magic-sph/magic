module outPar_mod
   !
   ! This module is used to compute several time-averaged radial profiles:
   ! fluxes, boundary layers, etc.
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use grid_blocking, only: radlatlon2spat
   use communications, only: gather_from_Rloc
   use truncation, only: n_r_max, n_r_maxMag, l_max, lm_max, l_maxMag, &
       &                 n_theta_max, n_phi_max
   use logic, only: l_viscBcCalc, l_anel, l_fluxProfs, l_mag_nl, &
       &            l_perpPar, l_save_out, l_temperature_diff,   &
       &            l_anelastic_liquid, l_mag
   use horizontal_data, only: gauss, osn2, gauss, O_sin_theta_E2, cosTheta, sn2
   use fields, only: s_Rloc, ds_Rloc, p_Rloc, dp_Rloc
   use physical_parameters, only: ek, prmag, OhmLossFac, ViscHeatFac, &
       &                          opr, kbots, ktops, ThExpNb, ekScaled
   use num_param, only: eScale
   use constants, only: pi, mass, osq4pi, sq4pi, half, two, four, third, one
   use radial_functions, only: r, or2, sigma, rho0, kappa, temp0, &
       &                       rscheme_oc, orho1, dLalpha0,       &
       &                       dLtemp0, beta, alpha0, or1, orho2, &
       &                       visc
   use radial_data, only: n_r_icb, nRstart, nRstop, nRstartMag, nRstopMag, &
       &                  n_r_cmb
   use num_param, only: tScale
   use output_data, only: tag
   use useful, only: cc2real, round_off
   use mean_sd, only: mean_sd_type
   use integration, only: rInt_R

   implicit none

   private
   
   type(mean_sd_type) :: fcond, fconv, fkin, fvisc, fres,  fpoyn
   type(mean_sd_type) :: Eperp, Epar, Eperpaxi, Eparaxi
   type(mean_sd_type) :: uh, duh, gradT2, entropy
   type(mean_sd_type) :: dlV, dlVc, Rm, Rol, uRol, dlPolpeak
   real(cp), allocatable :: fconvASr(:), fkinASr(:), fviscASr(:), gradT2ASr(:)
   real(cp), allocatable :: fpoynASr(:), fresASr(:), uhASr(:), duhASr(:)
   real(cp), allocatable :: EperpASr(:), EparASr(:), EperpaxiASr(:), EparaxiASr(:)

   integer :: n_perpPar_file, n_calls
   character(len=72) :: perpPar_file

   public initialize_outPar_mod, finalize_outPar_mod, outPar, outPerpPar, &
   &      get_fluxes, get_nlBlayers, get_perpPar

contains

   subroutine initialize_outPar_mod()
      !
      ! Memory allocation and file openings of several outputs (perpPar, fluxes,
      ! bLayers). Mostly time-averaged radial outputs.
      !

      n_calls = 0 
      call dlV%initialize(1,n_r_max)
      call dlVc%initialize(1,n_r_max)
      call Rm%initialize(1,n_r_max)
      call Rol%initialize(1,n_r_max)
      call uRol%initialize(1,n_r_max)
      call dlPolPeak%initialize(1,n_r_max)

      perpPar_file='perpPar.'//tag

      if ( l_viscBcCalc ) then
         call entropy%initialize(1,n_r_max)
         call uh%initialize(1,n_r_max)
         call duh%initialize(1,n_r_max)
         call gradT2%initialize(1,n_r_max)
         allocate( uhASr(nRstart:nRstop), duhASr(nRstart:nRstop) )
         allocate( gradT2ASr(nRstart:nRstop) )
         uhASr(:)    =0.0_cp
         duhASr(:)   =0.0_cp
         gradT2ASr(:)=0.0_cp
         bytes_allocated=bytes_allocated+3*(nRstop-nRstart+1)*SIZEOF_DEF_REAL
      end if

      if ( l_fluxProfs ) then
         call fcond%initialize(1,n_r_max)
         call fconv%initialize(1,n_r_max)
         call fkin%initialize(1,n_r_max)
         call fvisc%initialize(1,n_r_max)
         call fres%initialize(1,n_r_max)
         call fpoyn%initialize(1,n_r_max)
         allocate( fconvASr(nRstart:nRstop), fkinASr(nRstart:nRstop) )
         allocate( fviscASr(nRstart:nRstop) )
         fconvASr(:)=0.0_cp
         fkinASr(:) =0.0_cp
         fviscASr(:)=0.0_cp
         bytes_allocated=bytes_allocated+3*(nRstop-nRstart+1)*SIZEOF_DEF_REAL
         if ( l_mag ) then
            allocate( fresASr(nRstart:nRstop), fpoynASr(nRstartMag:nRstopMag) )
            fresASr(:) =0.0_cp
            fpoynASr(:)=0.0_cp
            bytes_allocated=bytes_allocated+2*(nRstop-nRstart+1)*SIZEOF_DEF_REAL
         end if
      end if

      if ( l_perpPar ) then
         call Eperp%initialize(1,n_r_max)
         call Epar%initialize(1,n_r_max)
         call Eperpaxi%initialize(1,n_r_max)
         call Eparaxi%initialize(1,n_r_max)
         allocate( EperpASr(nRstart:nRstop), EparASr(nRstart:nRstop) )
         allocate( EperpaxiASr(nRstart:nRstop), EparaxiASr(nRstart:nRstop) )
         EperpASr(:)   =0.0_cp
         EparASr(:)    =0.0_cp
         EperpaxiASr(:)=0.0_cp
         EparaxiASr(:) =0.0_cp
         bytes_allocated=bytes_allocated+4*(nRstop-nRstart+1)*SIZEOF_DEF_REAL

         if ( rank == 0 .and. (.not. l_save_out) ) then
            open(newunit=n_perpPar_file, file=perpPar_file, status='new')
         end if
      end if

   end subroutine initialize_outPar_mod
!-----------------------------------------------------------------------
   subroutine finalize_outPar_mod()
      !
      ! Closing and memory deallocation of outPar related outputs: fluxesR.TAG,
      ! perpar.TAG, bLayersR.TAG, ...
      !

      call dlV%finalize()
      call dlVc%finalize()
      call Rm%finalize()
      call Rol%finalize()
      call uRol%finalize()
      call dlPolPeak%finalize()

      if ( l_viscBcCalc ) then
         deallocate( uhASr, duhASr, gradT2ASr )
         call entropy%finalize()
         call uh%finalize()
         call duh%finalize()
         call gradT2%finalize()
      end if

      if ( l_fluxProfs ) then
         deallocate( fkinASr, fconvASr, fviscASr, fpoynASr, fresASr )
         call fcond%finalize()
         call fconv%finalize()
         call fkin%finalize()
         call fvisc%finalize()
         call fpoyn%finalize()
         call fres%finalize()
      end if

      if ( l_perpPar ) then
         deallocate( EperpaxiASr, EparaxiASr, EperpASr, EparASr)
         call Eperp%finalize()
         call Epar%finalize()
         call Eperpaxi%finalize()
         call Eparaxi%finalize()
         if ( rank == 0 .and. (.not. l_save_out) ) close(n_perpPar_file)
      end if

   end subroutine finalize_outPar_mod
!-----------------------------------------------------------------------
   subroutine outPar(timePassed, timeNorm, l_stop_time, ekinR, RolRu2,  &
              &      dlVR, dlVRc, dlPolPeakR, RmR )

      !--- Input of variables
      real(cp), intent(in) :: timePassed,timeNorm
      logical,  intent(in) :: l_stop_time
      real(cp), intent(in) :: RolRu2(n_r_max),dlPolPeakR(n_r_max)
      real(cp), intent(in) :: dlVR(n_r_max),dlVRc(n_r_max)
      real(cp), intent(in) :: ekinR(n_r_max)     ! kinetic energy w radius

      !--- Output of variables
      real(cp), intent(out):: RmR(n_r_max)

      !-- Local variables
      integer :: nR,fileHandle
      real(cp) :: ReR(n_r_max), RoR(n_r_max), RolR(n_r_max)
      character(len=76) :: filename
      real(cp) :: sR(nRstart:nRstop),fcR(nRstart:nRstop)
      real(cp) :: duhR_global(n_r_max), uhR_global(n_r_max)
      real(cp) :: gradT2R_global(n_r_max), sR_global(n_r_max)
      real(cp) :: fkinR_global(n_r_max), fcR_global(n_r_max)
      real(cp) :: fconvR_global(n_r_max), fviscR_global(n_r_max)
      real(cp) :: fresR_global(n_r_maxMag), fpoynR_global(n_r_maxMag)

      n_calls = n_calls+1

      if ( l_viscBcCalc ) then
         sR(:) = real(s_Rloc(1,:))

         duhASr(:)=half*duhASr(:) ! Normalisation for the theta integration
         uhASr(:) =half* uhASr(:) ! Normalisation for the theta integration
         gradT2ASr(:)=half*gradT2ASr(:) ! Normalisation for the theta integration

         call gather_from_RLoc(duhASR, duhR_global, 0)
         call gather_from_RLoc(uhASR, uhR_global, 0)
         call gather_from_RLoc(gradT2ASR, gradT2R_global, 0)
         call gather_from_RLoc(sR, sR_global, 0)
      end if

      if ( l_fluxProfs ) then
         if ( l_anelastic_liquid ) then
            if ( l_temperature_diff ) then
               do nR=nRstart,nRstop
                  fcR(nR)=-real(ds_Rloc(1,nR))*kappa(nR)*rho0(nR)* &
                  &        r(nR)*r(nR)*sq4pi
               end do
            else
               do nR=nRstart,nRstop
                  fcR(nR)=-kappa(nR)*r(nR)*r(nR)*                 &
                  &       sq4pi*( rho0(nR)*(real(ds_Rloc(1,nR))-  &
                  &       dLtemp0(nR)*real(s_Rloc(1,nR)))-ThExpNb*&
                  &       ViscHeatFac*alpha0(nR)*temp0(nR)*(      &
                  &       real(dp_Rloc(1,nR))+(dLalpha0(nR)-      &
                  &       beta(nR))*real(p_Rloc(1,nR))) )
               end do
            end if
         else
            if  ( l_temperature_diff ) then
               do nR=nRstart,nRstop
                  fcR(nR)=-sq4pi*r(nR)*r(nR)*kappa(nR)*rho0(nR)*temp0(nR)*&
                  &        (dLtemp0(nR)*real(s_Rloc(1,nR)) +              &
                  &                     real(ds_Rloc(1,nR))+              &
                  &        ViscHeatFac*ThExpNb*alpha0(nR)*                &
                  &        orho1(nR)*((dLalpha0(nR)+dLtemp0(nR)-beta(nR))*&
                  &                     real(p_Rloc(1,nR))+               &
                  &                     real(dp_Rloc(1,nR))))
               end do
            else
               do nR=nRstart,nRstop
                  fcR(nR)=-real(ds_Rloc(1,nR))*kappa(nR)*rho0(nR)* &
                  &        temp0(nR)*r(nR)*r(nR)*sq4pi
               end do
            end if
         end if

         call gather_from_Rloc(fkinASr, fkinR_global, 0)
         call gather_from_Rloc(fconvASr, fconvR_global, 0)
         call gather_from_Rloc(fviscASr, fviscR_global, 0)
         call gather_from_Rloc(fcR, fcR_global, 0)
         if ( l_mag_nl ) then
            call gather_from_Rloc(fpoynASr, fpoynR_global, 0)
            call gather_from_Rloc(fresASr, fresR_global, 0)
         end if
      end if

      if ( rank == 0 ) then
         do nR=1,n_r_max
            ! Re must be independant of the timescale
            ReR(nR)=sqrt(two*ekinR(nR)*or2(nR)/(4*pi*mass)/eScale)
            RoR(nR)=ReR(nR)*ekScaled
            if ( dlVR(nR) /= 0.0_cp ) then
               RolR(nR)=RoR(nR)/dlVR(nR)
            else
               RolR(nR)=RoR(nR)
            end if
            if ( l_mag_nl ) then
               RmR(nR)=ReR(nR)*prmag*sigma(nR)*r(nR)*r(nR)
            else
               RmR(nR)=ReR(nR)*r(nR)*r(nR)
            end if
         end do

         call dlV%compute(dlVR, n_calls, timePassed, timeNorm)
         call dlVc%compute(dlVRc, n_calls, timePassed, timeNorm)
         call Rol%compute(RolR, n_calls, timePassed, timeNorm)
         if ( l_anel ) then
            call uRol%compute(RolRu2, n_calls, timePassed, timeNorm)
         else
            call uRol%compute(RolR, n_calls, timePassed, timeNorm)
         end if
         RmR(:)=RmR(:)*sqrt(mass*orho1(:))*or2(:)
         call Rm%compute(RmR, n_calls, timePassed, timeNorm)
         call dlPolPeak%compute(dlPolPeakR, n_calls, timePassed, timeNorm)

         if ( l_viscBcCalc ) then
            call entropy%compute(osq4pi*sR_global, n_calls, timePassed, timeNorm)
            call uh%compute(uhR_global, n_calls, timePassed, timeNorm)
            call duh%compute(duhR_global, n_calls, timePassed, timeNorm)
            call gradT2%compute(gradT2R_global, n_calls, timePassed, timeNorm)
         end if

         if ( l_fluxProfs ) then
            call fcond%compute(opr*fcR_global, n_calls, timePassed, timeNorm)
            call fconv%compute(fconvR_global, n_calls, timePassed, timeNorm)
            call fkin%compute(ViscHeatFac*fkinR_global, n_calls, timePassed, &
                 &            timeNorm)
            call fvisc%compute(ViscHeatFac*fviscR_global, n_calls, timePassed, &
                 &             timeNorm)
            if ( l_mag_nl ) then
               call fres%compute(OhmLossFac*fresR_global, n_calls, timePassed, &
                    &            timeNorm)
               call fpoyn%compute(prmag*OhmLossFac*fpoynR_global, n_calls, &
                    &             timePassed, timeNorm)
            end if
         end if

         if ( l_stop_time ) then

            call dlV%finalize_SD(timeNorm)
            call dlVc%finalize_SD(timeNorm)
            call Rol%finalize_SD(timeNorm)
            call uRol%finalize_SD(timeNorm)
            call Rm%finalize_SD(timeNorm)
            call dlPolPeak%finalize_SD(timeNorm)

            if ( l_viscBcCalc ) then
               call entropy%finalize_SD(timeNorm)
               call uh%finalize_SD(timeNorm)
               call duh%finalize_SD(timeNorm)
               call gradT2%finalize_SD(timeNorm)
            end if

            if ( l_fluxProfs ) then
               call fcond%finalize_SD(timeNorm)
               call fconv%finalize_SD(timeNorm)
               call fkin%finalize_SD(timeNorm)
               call fvisc%finalize_SD(timeNorm)
               if ( l_mag_nl ) then
                  call fpoyn%finalize_SD(timeNorm)
                  call fres%finalize_SD(timeNorm)
               end if
            end if

            !----- Output into paR.TAG file:
            filename='parR.'//tag
            open(newunit=fileHandle, file=filename, status='unknown')
            do nR=1,n_r_max
               write(fileHandle,'(ES20.10,6ES15.7,6ES13.5)')               &
               &     r(nR),round_off(Rm%mean(nR),maxval(Rm%mean)),         &
               &     round_off(Rol%mean(nR),maxval(Rol%mean)),             &
               &     round_off(uRol%mean(nR),maxval(uRol%mean)),           &
               &     round_off(dlV%mean(nR),maxval(dlV%mean)),             &
               &     round_off(dlVc%mean(nR),maxval(dlVc%mean)),           &
               &     round_off(dlPolPeak%mean(nR),maxval(dlPolPeak%mean)), &
               &     round_off(Rm%SD(nR),maxval(Rm%SD)),                   &
               &     round_off(Rol%SD(nR),maxval(Rol%SD)),                 &
               &     round_off(uRol%SD(nR),maxval(uRol%SD)),               &
               &     round_off(dlV%SD(nR),maxval(dlV%SD)),                 &
               &     round_off(dlVc%SD(nR),maxval(dlVc%SD)),               &
               &     round_off(dlPolPeak%SD(nR),maxval(dlPolPeak%SD))
            end do
            close(fileHandle)

            if ( l_viscBcCalc ) then
               filename='bLayersR.'//tag
               open(newunit=fileHandle, file=filename, status='unknown')
               do nR=1,n_r_max
                  write(fileHandle,'(ES20.10,4ES15.7,4ES13.4)')                 &
                  &     r(nR),round_off(entropy%mean(nR),maxval(entropy%mean)), &
                  &     round_off(uh%mean(nR),maxval(uh%mean)),                 &
                  &     round_off(duh%mean(nR),maxval(duh%mean)),               &
                  &     round_off(gradT2%mean(nR),maxval(gradT2%mean)),         &
                  &     round_off(entropy%SD(nR),maxval(entropy%SD)),           &
                  &     round_off(uh%SD(nR),maxval(uh%SD)),                     &
                  &     round_off(duh%SD(nR),maxval(duh%SD)),                   &
                  &     round_off(gradT2%SD(nR),maxval(gradT2%SD))
               end do
               close(fileHandle)
            end if

            if ( l_fluxProfs ) then
               filename='fluxesR.'//tag
               open(newunit=fileHandle, file=filename, status='unknown')
               do nR=1,n_r_max
                  write(fileHandle,'(ES20.10,7ES15.7,7ES13.5)')             &
                  &     r(nR),round_off(fcond%mean(nR),maxval(fcond%mean)), &
                  &     round_off(fconv%mean(nR),maxval(fconv%mean)),       &
                  &     round_off(fkin%mean(nR),maxval(fkin%mean)),         &
                  &     round_off(fvisc%mean(nR),maxval(fvisc%mean)),       &
                  &     round_off(fpoyn%mean(nR),maxval(fpoyn%mean)),       &
                  &     round_off(fres%mean(nR),maxval(fres%mean)),         &
                  &     round_off(fcond%SD(nR),maxval(fcond%SD)),           &
                  &     round_off(fconv%SD(nR),maxval(fconv%SD)),           &
                  &     round_off(fkin%SD(nR),maxval(fkin%SD)),             &
                  &     round_off(fvisc%SD(nR),maxval(fvisc%SD)),           &
                  &     round_off(fpoyn%SD(nR),maxval(fpoyn%SD)),           &
                  &     round_off(fres%SD(nR),maxval(fres%SD))
               end do
               close(fileHandle)
            end if

         end if ! l_stop_time ?

      end if ! rank0

   end subroutine outPar
!----------------------------------------------------------------------------
   subroutine outPerpPar(time,timePassed,timeNorm,l_stop_time)
      !
      ! This subroutine handles the writing the time series perpar.tag which
      ! stores kinetic energy content perpendicular and parallel to rotation
      ! axis.
      !

      !--- Input of variables
      real(cp), intent(in) :: time,timePassed,timeNorm
      logical,  intent(in) :: l_stop_time

      !--- Local variables
      integer :: nR,fileHandle
      character(len=76) :: filename

      real(cp) :: EperpR_global(n_r_max), EparR_global(n_r_max)
      real(cp) :: EperpaxiR_global(n_r_max), EparaxiR_global(n_r_max)
      real(cp) :: EperpT,EparT,EperpaxT,EparaxT

      EperpASr(:)   =half*EperpASr(:)    ! Normalisation for the theta integration
      EparASr(:)    =half*EparASr(:)     ! Normalisation for the theta integration
      EperpaxiASr(:)=half*EperpaxiASr(:) ! Normalisation for the theta integration
      EparaxiASr(:) =half*EparaxiASr(:)  ! Normalisation for the theta integration

      call gather_from_Rloc(EperpASr, EperpR_global, 0)
      call gather_from_Rloc(EparASr, EparR_global, 0)
      call gather_from_Rloc(EperpaxiASr, EperpaxiR_global, 0)
      call gather_from_Rloc(EparaxiASr, EparaxiR_global, 0)

      if ( rank == 0 ) then
         EperpT  =four*pi*rInt_R(EperpR_global*r*r,r,rscheme_oc)
         EparT   =four*pi*rInt_R(EparR_global*r*r,r,rscheme_oc)
         EperpaxT=four*pi*rInt_R(EperpaxiR_global*r*r,r,rscheme_oc)
         EparaxT =four*pi*rInt_R(EparaxiR_global*r*r,r,rscheme_oc)

         !-- Output
         if ( l_save_out ) then
            open(newunit=n_perpPar_file, file=perpPar_file, &
            &    status='unknown', position='append')
         end if
         write(n_perpPar_file,'(1P,ES20.12,4ES16.8)')  time*tScale,     & ! 1
         &                                       EperpT,EparT, EperpaxT,EparaxT
         if ( l_save_out ) close(n_perpPar_file)

         call Eperp%compute(EperpR_global, n_calls, timePassed, timeNorm)
         call Epar%compute(EparR_global, n_calls, timePassed, timeNorm)
         call Eperpaxi%compute(EperpaxiR_global, n_calls, timePassed, timeNorm)
         call Eparaxi%compute(EparaxiR_global, n_calls, timePassed, timeNorm)

         if ( l_stop_time ) then
            call Eperp%finalize_SD(timeNorm)
            call Epar%finalize_SD(timeNorm)
            call Eperpaxi%finalize_SD(timeNorm)
            call Eparaxi%finalize_SD(timeNorm)

            filename='perpParR.'//tag
            open(newunit=fileHandle, file=filename, status='unknown')
            do nR=1,n_r_max
               write(fileHandle,'(ES20.10,4ES15.7,4ES13.5)')               &
               &     r(nR),round_off(Eperp%mean(nR),maxval(Eperp%mean)),   &
               &     round_off(Epar%mean(nR),maxval(Epar%mean)),           &
               &     round_off(Eperpaxi%mean(nR),maxval(Eperpaxi%mean)),   &
               &     round_off(Eparaxi%mean(nR),maxval(Eparaxi%mean)),     &
               &     round_off(Eperp%SD(nR),maxval(Eperp%SD)),             &
               &     round_off(Epar%SD(nR),maxval(Epar%SD)),               &
               &     round_off(Eperpaxi%SD(nR),maxval(Eperpaxi%SD)),       &
               &     round_off(Eparaxi%SD(nR),maxval(Eperpaxi%SD))
            end do
            close(fileHandle)
         end if
      end if

   end subroutine outPerpPar
!----------------------------------------------------------------------------
   subroutine get_fluxes(vr,vt,vp,dvrdr,dvtdr,dvpdr,dvrdt,dvrdp,sr,pr,br,bt, &
              &          bp,cbt,cbp,nR)
      !
      !   This routine computes the various contribution to heat fluxes:
      !
      !     * Convective flux: :math:`F_c= \rho T (u_r s)`
      !     * Kinetic flux: :math:`F_k = 1/2\,\rho u_r (u_r^2+u_\theta^2+u_\phi^2)`
      !     * Viscous flux: :math:`F_= -(u \cdot S )_r`)
      !
      !   If the run is magnetic, then this routine also computes:
      !
      !     * Poynting flux
      !     * Resistive flux
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(*),vt(*),vp(*)
      real(cp), intent(in) :: dvrdr(*),dvtdr(*),dvpdr(*)
      real(cp), intent(in) :: dvrdt(*),dvrdp(*)
      real(cp), intent(in) :: sr(*),pr(*)
      real(cp), intent(in) :: br(*),bt(*),bp(*)
      real(cp), intent(in) :: cbt(*),cbp(*)

      !-- Local variables:
      real(cp) :: fkinAS,fconvAS,fviscAS,fresAS,fpoynAS
      real(cp) :: fkin,fconv,phiNorm,fvisc,fpoyn,fres
      integer :: nTheta,nThetaNHS,nPhi,nelem

      phiNorm=two*pi/real(n_phi_max,cp)

      fkinAS =0.0_cp
      fconvAS=0.0_cp
      fviscAS=0.0_cp
      fvisc  =0.0_cp
      !$omp parallel do default(shared)                                  &
      !$omp& private(nTheta, nPhi, nelem, nThetaNHS, fconv, fkin, fvisc) &
      !$omp& reduction(+:fkinAS,fconvAS,fviscAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nelem = radlatlon2spat(nTheta,nPhi,nR)
            nThetaNHS=(nTheta+1)/2
            if ( l_anelastic_liquid ) then
               fconv=vr(nelem)*sr(nelem)
            else
               fconv=temp0(nr)*vr(nelem)*sr(nelem)     +            &
               &          ViscHeatFac*ThExpNb*alpha0(nr)*temp0(nr)* &
               &          orho1(nr)*vr(nelem)*pr(nelem)
            end if

            fkin=half*or2(nR)*orho2(nR)*(osn2(nThetaNHS)*(       &
            &    vt(nelem)*vt(nelem) + vp(nelem)*vp(nelem) )+    &
            &    or2(nR)*vr(nelem)*vr(nelem) )*vr(nelem)

            if ( nR/=n_r_icb .and. nR/=n_r_cmb ) then
               fvisc=-two*visc(nR)*orho1(nR)*vr(nelem)*or2(nR)* (     &
               &                             dvrdr(nelem)             &
               & -(two*or1(nR)+two*third*beta(nR))*vr(nelem) )-       &
               &                       visc(nR)*orho1(nR)*vt(nelem)*  &
               &                      osn2(nThetaNHS)* (              &
               &                       or2(nR)*dvrdt(nelem)           &
               &                              +dvtdr(nelem)           &
               &       -(two*or1(nR)+beta(nR))*vt(nelem) )  -         &
               &       visc(nR)*orho1(nR)*vp(nelem)*                  &
               &                         osn2(nThetaNHS)* (           &
               &                       or2(nR)*dvrdp(nelem)           &
               &                              +dvpdr(nelem)           &
               &       -(two*or1(nR)+beta(nR))*vp(nelem) )
            end if

            fkinAS = fkinAS+phiNorm*gauss(nThetaNHS)*fkin
            fconvAS=fconvAS+phiNorm*gauss(nThetaNHS)*fconv
            fviscAS=fviscAS+phiNorm*gauss(nThetaNHS)*fvisc
         end do
      end do
      !$omp end parallel do

      fkinASr(nR) =fkinAS
      fconvASr(nR)=fconvAS
      fviscASr(nR)=fviscAS

      if ( l_mag_nl) then
         fresAS =0.0_cp
         fpoynAS=0.0_cp
         !$omp parallel do default(shared)                      &
         !$omp& private(nTheta, nPhi, nThetaNHS, fres, fpoyn)   &
         !$omp& reduction(+:fresAS,fpoynAS)
         do nPhi=1,n_phi_max
            do nTheta=1,n_theta_max
               nelem = radlatlon2spat(nTheta,nPhi,nR)
               nThetaNHS=(nTheta+1)/2
               fres =osn2(nThetaNHS)*(cbt(nelem)*bp(nelem)  - &
               &                      cbp(nelem)*bt(nelem) )

               fpoyn=-orho1(nR)*or2(nR)*osn2(nThetaNHS)*(   &
               &           vp(nelem)*br(nelem)*bp(nelem)  - &
               &           vr(nelem)*bp(nelem)*bp(nelem)  - &
               &           vr(nelem)*bt(nelem)*bt(nelem)  + &
               &           vt(nelem)*br(nelem)*bt(nelem) )

               fresAS = fresAS+phiNorm*gauss(nThetaNHS)*fres
               fpoynAS=fpoynAS+phiNorm*gauss(nThetaNHS)*fpoyn
            end do
         end do
         !$omp end parallel do

         fpoynASr(nR)=fpoynAS
         fresASr(nR) =fresAS
      end if

   end subroutine get_fluxes
!----------------------------------------------------------------------------
   subroutine get_nlBLayers(vt,vp,dvtdr,dvpdr,dsdr,dsdt,dsdp,nR)
      !
      !   This subroutine calculates the axisymmetric contributions of:
      !
      !     * the horizontal velocity :math:`u_h = \sqrt{u_\theta^2+u_\phi^2}`
      !     * its radial derivative :math:`|\partial u_h/\partial r|`
      !     * The thermal dissipation rate :math:`(\nabla T)^2`
      !
      !   This subroutine is used when one wants to evaluate viscous and thermal
      !   dissipation layers
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vt(*),vp(*)
      real(cp), intent(in) :: dvtdr(*),dvpdr(*)
      real(cp), intent(in) :: dsdr(*),dsdt(*),dsdp(*)

      !-- Local variables:
      real(cp):: uhAS,duhAS,gradsAS,uh,duh,phiNorm,grads
      integer :: nTheta,nPhi,nThetaNHS,nelem

      phiNorm=one/real(n_phi_max,cp)
      uhAS   =0.0_cp
      duhAS  =0.0_cp
      gradsAS=0.0_cp

      !--- Horizontal velocity uh and duh/dr + (grad T)**2
      !$omp parallel do default(shared)                              &
      !$omp& private(nTheta, nThetaNHS, nPhi, uh, duh, grads, nelem) &
      !$omp& reduction(+:uhAS,duhAS,gradsAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nelem = radlatlon2spat(nTheta,nPhi,nR)
            nThetaNHS=(nTheta+1)/2
            uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(     &
            &             vt(nelem)*vt(nelem)+vp(nelem)*vp(nelem)  )
            duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(                     &
            &   dvtdr(nelem)*vt(nelem)-(or1(nR)+beta(nR))*vt(nelem)*vt(nelem)+ &
            &   dvpdr(nelem)*vp(nelem)-(or1(nR)+beta(nR))*vp(nelem)*vp(nelem) )

            grads = dsdr(nelem)*dsdr(nelem)+or2(nR)*O_sin_theta_E2(nTheta)*(  &
            &       dsdt(nelem)*dsdt(nelem)+dsdp(nelem)*dsdp(nelem) )

            uhAS=uhAS+phiNorm*gauss(nThetaNHS)*sqrt(uh)
            if (uh /= 0.0_cp) then
               duhAS=duhAS+phiNorm*gauss(nThetaNHS)*abs(duh)/sqrt(uh)
            end if
            gradsAS=gradsAS+phiNorm*gauss(nThetaNHS)*grads
         end do
      end do
      !$omp end parallel do

      uhASr(nR)    =uhAS
      duhASr(nR)   =duhAS
      gradT2ASr(nR)=gradsAS

   end subroutine get_nlBLayers
!----------------------------------------------------------------------------
   subroutine get_perpPar(vr,vt,vp,nR)
      !
      !   This subroutine calculates the energies parallel and perpendicular
      !   to the rotation axis
      !
      !     * :math:`E_\perp = 0.5 (v_s^2+v_\phi^2)` with
      !       :math:`v_s= v_r\sin\theta+v_\theta\cos\theta`
      !     * :math:`E_\parallel  = 0.5v_z^2` with
      !       :math:`v_z= v_r\cos\theta-v_\theta*\sin\theta`
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(*),vt(*),vp(*)

      !-- Local variables:
      real(cp) :: vras(n_theta_max),vtas(n_theta_max),vpas(n_theta_max),phiNorm
      real(cp) :: Eperp,Epar,Eperpaxi,Eparaxi
      real(cp) :: EperpAS,EparAS,EperpaxiAS,EparaxiAS
      integer :: nTheta,nPhi,nThetaNHS,nelem

      phiNorm=one/real(n_phi_max,cp)
      EperpAS   =0.0_cp
      EparAS    =0.0_cp
      EperpaxiAS=0.0_cp
      EparaxiAS =0.0_cp
      vras(:)   =0.0_cp
      vtas(:)   =0.0_cp
      vpas(:)   =0.0_cp

      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nelem = radlatlon2spat(nTheta,nPhi,nR)
            vras(nTheta)=vras(nTheta)+vr(nelem)
            vtas(nTheta)=vtas(nTheta)+vt(nelem)
            vpas(nTheta)=vpas(nTheta)+vp(nelem)
         end do
      end do
      vras(:)=vras(:)*phiNorm
      vtas(:)=vtas(:)*phiNorm
      vpas(:)=vpas(:)*phiNorm

      !$omp parallel do default(shared)                 &
      !$omp& private(nTheta,nPhi,nThetaNHS,nelem)       &
      !$omp& private(Eperp, Epar, Eperpaxi, Eparaxi)    &
      !$omp& reduction(+:EparAS,EperpAS,EparaxiAS,EperpaxiAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nelem = radlatlon2spat(nTheta,nPhi,nR)
            nThetaNHS=(nTheta+1)/2

            Eperp=half*or2(nR)*orho2(nR)*(                             &
            &       or2(nR)*sn2(nThetaNHS)*      vr(nelem)*vr(nelem) + &
            &       (osn2(nThetaNHS)-one)*       vt(nelem)*vt(nelem) + &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nelem)*vt(nelem) + &
            &       osn2(nThetaNHS)*             vp(nelem)*vp(nelem) )

            Epar =half*or2(nR)*orho2(nR)*(                             &
            &       or2(nR)*(one-sn2(nThetaNHS))*vr(nelem)*vr(nelem) + &
            &                                    vt(nelem)*vt(nelem) - &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nelem)*vt(nelem) )

            Eperpaxi=half*or2(nR)*orho2(nR)*(                                  &
            &         or2(nR)*sn2(nThetaNHS)*      vras(nTheta)*vras(nTheta) + &
            &         (osn2(nThetaNHS)-one)*       vtas(nTheta)*vtas(nTheta) + &
            &         two*or1(nR)*cosTheta(nTheta)*vras(nTheta)*vtas(nTheta) + &
            &         osn2(nThetaNHS)*             vpas(nTheta)*vpas(nTheta) )

            Eparaxi =half*or2(nR)*orho2(nR)*(                                  &
            &         or2(nR)*(one-sn2(nThetaNHS))*vras(nTheta)*vras(nTheta) + &
            &                                      vtas(nTheta)*vtas(nTheta) - &
            &         two*or1(nR)*cosTheta(nTheta)*vras(nTheta)*vtas(nTheta) )

            EperpAS   =   EperpAS+phiNorm*gauss(nThetaNHS)*Eperp
            EparAS    =    EparAS+phiNorm*gauss(nThetaNHS)*Epar
            EperpaxiAS=EperpaxiAS+phiNorm*gauss(nThetaNHS)*Eperpaxi
            EparaxiAS = EparaxiAS+phiNorm*gauss(nThetaNHS)*Eparaxi
         end do
      end do
      !$omp end parallel do

      EperpASr(nR)   =EperpAS
      EparASr(nR)    =EparAS
      EperpaxiASr(nR)=EperpaxiAS
      EparaxiASr(nR) =EparaxiAS

   end subroutine get_perpPar
!----------------------------------------------------------------------------
end module outPar_mod
