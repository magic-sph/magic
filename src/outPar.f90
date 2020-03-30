module outPar_mod
   !
   ! This module is used to compute several time-averaged radial profiles:
   ! fluxes, boundary layers, etc.
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use communications, only: gather_from_Rloc
   use truncation, only: n_r_max, n_r_maxMag, l_max, lm_max, l_maxMag,&
       &                 n_r_icb, nRstart, nRstop, nRstartMag,        &
       &                 nRstopMag
   use blocking, only: nfs, nThetaBs, sizeThetaB, lm2m
   use logic, only: l_viscBcCalc, l_anel, l_fluxProfs, l_mag_nl, &
       &            l_perpPar, l_save_out, l_temperature_diff,   &
       &            l_anelastic_liquid
   use horizontal_data, only: gauss
   use fields, only: s_Rloc, ds_Rloc, p_Rloc, dp_Rloc
   use physical_parameters, only: ek, prmag, OhmLossFac, ViscHeatFac, &
       &                          opr, kbots, ktops, ThExpNb, ekScaled
   use num_param, only: eScale
   use constants, only: pi, mass, osq4pi, sq4pi, half, two, four
   use radial_functions, only: r, or2, sigma, rho0, kappa, temp0, &
       &                       rscheme_oc, orho1, dLalpha0,       &
       &                       dLtemp0, beta, alpha0
   use num_param, only: tScale
   use output_data, only: tag
   use useful, only: cc2real, round_off
   use mean_sd, only: mean_sd_type
   use integration, only: rInt_R
#ifdef WITH_SHTNS
   use shtns, only: axi_to_spat
#else
   use legendre_spec_to_grid, only: lmAS2pt
#endif

   implicit none

   private
   
   type(mean_sd_type) :: fcond, fconv, fkin, fvisc, fres,  fpoyn
   type(mean_sd_type) :: Eperp, Epar, Eperpaxi, Eparaxi
   type(mean_sd_type) :: uh, duh, gradT2, entropy
   type(mean_sd_type) :: dlV, dlVc, Rm, Rol, uRol, dlPolpeak

   integer :: n_perpPar_file, n_calls
   character(len=72) :: perpPar_file

   public initialize_outPar_mod, finalize_outPar_mod, outPar, outPerpPar

contains

   subroutine initialize_outPar_mod
      !
      ! Memory allocation and file openings
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
      end if

      if ( l_fluxProfs ) then
         call fcond%initialize(1,n_r_max)
         call fconv%initialize(1,n_r_max)
         call fkin%initialize(1,n_r_max)
         call fvisc%initialize(1,n_r_max)
         call fres%initialize(1,n_r_max)
         call fpoyn%initialize(1,n_r_max)
      end if

      if ( l_perpPar ) then
         call Eperp%initialize(1,n_r_max)
         call Epar%initialize(1,n_r_max)
         call Eperpaxi%initialize(1,n_r_max)
         call Eparaxi%initialize(1,n_r_max)

         if ( l_master_rank .and. (.not. l_save_out) ) then
            open(newunit=n_perpPar_file, file=perpPar_file, status='new')
         end if
      end if

   end subroutine initialize_outPar_mod
!-----------------------------------------------------------------------
   subroutine finalize_outPar_mod

      call dlV%finalize()
      call dlVc%finalize()
      call Rm%finalize()
      call Rol%finalize()
      call uRol%finalize()
      call dlPolPeak%finalize()

      if ( l_viscBcCalc ) then
         call entropy%finalize()
         call uh%finalize()
         call duh%finalize()
         call gradT2%finalize()
      end if

      if ( l_fluxProfs ) then
         call fcond%finalize()
         call fconv%finalize()
         call fkin%finalize()
         call fvisc%finalize()
         call fpoyn%finalize()
         call fres%finalize()
      end if

      if ( l_perpPar ) then
         call Eperp%finalize()
         call Epar%finalize()
         call Eperpaxi%finalize()
         call Eparaxi%finalize()
         if ( l_master_rank .and. (.not. l_save_out) ) close(n_perpPar_file)
      end if

   end subroutine finalize_outPar_mod
!-----------------------------------------------------------------------
   subroutine outPar(timePassed, timeNorm, l_stop_time, ekinR, RolRu2,  &
              &      dlVR, dlVRc, dlPolPeakR, uhLMr, duhLMr, gradsLMr,  &
              &      fconvLMr, fkinLMr, fviscLMr, fpoynLMr, fresLMr, RmR)

      !--- Input of variables
      real(cp), intent(in) :: timePassed,timeNorm
      logical,  intent(in) :: l_stop_time
      real(cp), intent(in) :: RolRu2(n_r_max),dlPolPeakR(n_r_max)
      real(cp), intent(in) :: dlVR(n_r_max),dlVRc(n_r_max)
      real(cp), intent(in) :: ekinR(n_r_max)     ! kinetic energy w radius
      real(cp), intent(in) :: uhLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: duhLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: gradsLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: fkinLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: fconvLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: fviscLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp), intent(in) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)

      !--- Output of variables
      real(cp), intent(out):: RmR(n_r_max)

      !-- Local variables
      integer :: nR,n
      real(cp) :: ReR(n_r_max), RoR(n_r_max), RolR(n_r_max)
      character(len=76) :: filename
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS
      real(cp) :: duhR(nRstart:nRstop), uhR(nRstart:nRstop)
      real(cp) :: gradT2R(nRstart:nRstop), sR(nRstart:nRstop)
      real(cp) :: fkinR(nRstart:nRstop), fcR(nRstart:nRstop)
      real(cp) :: fconvR(nRstart:nRstop), fviscR(nRstart:nRstop)
      real(cp) :: fresR(nRstartMag:nRstopMag),fpoynR(nRstartMag:nRstopMag)
      real(cp) :: duhR_global(n_r_max), uhR_global(n_r_max)
      real(cp) :: gradT2R_global(n_r_max), sR_global(n_r_max)
      real(cp) :: fkinR_global(n_r_max), fcR_global(n_r_max)
      real(cp) :: fconvR_global(n_r_max), fviscR_global(n_r_max)
      real(cp) :: fresR_global(n_r_maxMag), fpoynR_global(n_r_maxMag)
      real(cp) :: duhTh(nfs), uhTh(nfs), gradT2Th(nfs)
      real(cp) :: fkinTh(nfs), fconvTh(nfs), fviscTh(nfs), fresTh(nfs), fpoynTh(nfs)

      integer :: fileHandle

      n_calls = n_calls+1


      if ( l_viscBcCalc ) then
         sR(:) = real(s_Rloc(1,:))

         do nR=nRstart,nRstop
            uhR(nR)    =0.0_cp
            gradT2R(nR)=0.0_cp
            duhR(nR)   =0.0_cp
#ifdef WITH_SHTNS
            call axi_to_spat(duhLMr(:,nR),duhTh)
            call axi_to_spat(uhLMr(:,nR),uhTh)
            call axi_to_spat(gradsLMr(:,nR),gradT2Th)
#endif
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
#ifndef WITH_SHTNS
               call lmAS2pt(duhLMr(:,nR),duhTh,nThetaStart,sizeThetaB)
               call lmAS2pt(uhLMr(:,nR),uhTh,nThetaStart,sizeThetaB)
               call lmAS2pt(gradsLMr(:,nR),gradT2Th,nThetaStart,sizeThetaB)
#endif
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
                  duhR(nR)=duhR(nR)+gauss(nThetaNHS)*duhTh(nThetaBlock)
                  uhR(nR) =uhR(nR) +gauss(nThetaNHS)* uhTh(nThetaBlock)
                  gradT2R(nR)=gradT2R(nR)+gauss(nThetaNHS)*gradT2Th(nThetaBlock)
               end do
            end do
         end do
         duhR=half*duhR ! Normalisation for the theta integration
         uhR =half* uhR ! Normalisation for the theta integration
         gradT2R =half*gradT2R ! Normalisation for the theta integration

         call gather_from_RLoc(duhR, duhR_global, 0)
         call gather_from_RLoc(uhR, uhR_global, 0)
         call gather_from_RLoc(gradT2R, gradT2R_global, 0)
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
         do nR=nRstart,nRstop
            fkinR(nR) =0.0_cp
            fconvR(nR)=0.0_cp
            fviscR(nR)=0.0_cp
#ifdef WITH_SHTNS
            call axi_to_spat(fkinLMr(:,nR),fkinTh)
            call axi_to_spat(fconvLMr(:,nR),fconvTh)
            call axi_to_spat(fviscLMr(:,nR),fviscTh)
#endif
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
#ifndef WITH_SHTNS
               call lmAS2pt(fkinLMr(:,nR),fkinTh,nThetaStart,sizeThetaB)
               call lmAS2pt(fconvLMr(:,nR),fconvTh,nThetaStart,sizeThetaB)
               call lmAS2pt(fviscLMr(:,nR),fviscTh,nThetaStart,sizeThetaB)
#endif
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
                  fkinR(nR) =fkinR(nR) +gauss(nThetaNHS)* fkinTh(nThetaBlock)
                  fconvR(nR)=fconvR(nR)+gauss(nThetaNHS)*fconvTh(nThetaBlock)
                  fviscR(nR)=fviscR(nR)+gauss(nThetaNHS)*fviscTh(nThetaBlock)
               end do
            end do
         end do

         if ( l_mag_nl ) then
            do nR=nRstart,nRstop
               fresR(nR) =0.0_cp
               fpoynR(nR)=0.0_cp
#ifdef WITH_SHTNS
               call axi_to_spat(fpoynLMr(:,nR),fpoynTh)
               call axi_to_spat(fresLMr(:,nR),fresTh)
#endif
               do n=1,nThetaBs ! Loop over theta blocks
                  nTheta=(n-1)*sizeThetaB
                  nThetaStart=nTheta+1
#ifndef WITH_SHTNS
                  call lmAS2pt(fpoynLMr(:,nR),fpoynTh,nThetaStart,sizeThetaB)
                  call lmAS2pt(fresLMr(:,nR),fresTh,nThetaStart,sizeThetaB)
#endif
                  do nThetaBlock=1,sizeThetaB
                     nTheta=nTheta+1
                     nThetaNHS=(nTheta+1)/2
                     fpoynR(nR)=fpoynR(nR)+gauss(nThetaNHS)*fpoynTh(nThetaBlock)
                     fresR(nR) =fresR(nR) +gauss(nThetaNHS)*fresTh(nThetaBlock)
                  end do
               end do
            end do
         end if

         call gather_from_Rloc(fkinR, fkinR_global, 0)
         call gather_from_Rloc(fconvR, fconvR_global, 0)
         call gather_from_Rloc(fviscR, fviscR_global, 0)
         call gather_from_Rloc(fcR, fcR_global, 0)
         if ( l_mag_nl ) then
            call gather_from_Rloc(fpoynR, fpoynR_global, 0)
            call gather_from_Rloc(fresR, fresR_global, 0)
         end if

      end if


      if ( coord_r == 0 ) then
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
            if (l_master_rank) then
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
            end if

         end if ! l_stop_time ?

      end if ! rank0

   end subroutine outPar
!----------------------------------------------------------------------------
   subroutine outPerpPar(time,timePassed,timeNorm,l_stop_time, &
              &          EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)


      !--- Input of variables
      real(cp), intent(in) :: time,timePassed,timeNorm
      logical,  intent(in) :: l_stop_time
      real(cp), intent(in) :: EparLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: EperpLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: EparaxiLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: EperpaxiLMr(l_max+1,nRstart:nRstop)

      !--- Local variables
      integer :: nR,n,nTheta,nThetaStart,nThetaBlock,nThetaNHS
      character(len=76) :: filename

      real(cp) ::EperpaxiR(nRstart:nRstop), EparaxiR(nRstart:nRstop)
      real(cp) :: EperpR(nRstart:nRstop), EparR(nRstart:nRstop)
      real(cp) :: EperpR_global(n_r_max), EparR_global(n_r_max)
      real(cp) :: EperpaxiR_global(n_r_max), EparaxiR_global(n_r_max)
      real(cp) :: EperpTh(nfs), EparTh(nfs), EperpaxiTh(nfs), EparaxiTh(nfs)
      real(cp) :: EperpT,EparT,EperpaxT,EparaxT

      integer :: fileHandle

      do nR=nRstart,nRstop
         EperpR(nR)   =0.0_cp
         EparR(nR)    =0.0_cp
         EparaxiR(nR) =0.0_cp
         EperpaxiR(nR)=0.0_cp
#ifdef WITH_SHTNS
         call axi_to_spat(EperpLMr(:,nR),EperpTh)
         call axi_to_spat(EparLMr(:,nR),EparTh)
         call axi_to_spat(EperpaxiLMr(:,nR),EperpaxiTh)
         call axi_to_spat(EparaxiLMr(:,nR),EparaxiTh)
#endif
         do n=1,nThetaBs ! Loop over theta blocks
            nTheta=(n-1)*sizeThetaB
            nThetaStart=nTheta+1
#ifndef WITH_SHTNS
            call lmAS2pt(EperpLMr(:,nR),EperpTh,nThetaStart,sizeThetaB)
            call lmAS2pt(EparLMr(:,nR),EparTh,nThetaStart,sizeThetaB)
            call lmAS2pt(EperpaxiLMr(:,nR),EperpaxiTh,nThetaStart,sizeThetaB)
            call lmAS2pt(EparaxiLMr(:,nR),EparaxiTh,nThetaStart,sizeThetaB)
#endif
            do nThetaBlock=1,sizeThetaB
               nTheta=nTheta+1
               nThetaNHS=(nTheta+1)/2
               EperpR(nR)=EperpR(nR)+gauss(nThetaNHS)*EperpTh(nThetaBlock)
               EparR(nR) =EparR(nR) +gauss(nThetaNHS)* EparTh(nThetaBlock)
               EperpaxiR(nR)=EperpaxiR(nR)+gauss(nThetaNHS)*EperpaxiTh(nThetaBlock)
               EparaxiR(nR)=EparaxiR(nR)+gauss(nThetaNHS)*EparaxiTh(nThetaBlock)
            end do
         end do
      end do
      EperpR   =half*EperpR    ! Normalisation for the theta integration
      EparR    =half*EparR     ! Normalisation for the theta integration
      EperpaxiR=half*EperpaxiR ! Normalisation for the theta integration
      EparaxiR =half*EparaxiR  ! Normalisation for the theta integration

      call gather_from_Rloc(EperpR, EperpR_global, 0)
      call gather_from_Rloc(EparR, EparR_global, 0)
      call gather_from_Rloc(EperpaxiR, EperpaxiR_global, 0)
      call gather_from_Rloc(EparaxiR, EparaxiR_global, 0)

      if ( coord_r == 0 ) then
         EperpT  =four*pi*rInt_R(EperpR_global*r*r,r,rscheme_oc)
         EparT   =four*pi*rInt_R(EparR_global*r*r,r,rscheme_oc)
         EperpaxT=four*pi*rInt_R(EperpaxiR_global*r*r,r,rscheme_oc)
         EparaxT =four*pi*rInt_R(EparaxiR_global*r*r,r,rscheme_oc)

         !-- Output
         if (l_master_rank) then
            if ( l_save_out ) then
               open(newunit=n_perpPar_file, file=perpPar_file, &
               &    status='unknown', position='append')
            end if
            write(n_perpPar_file,'(1P,ES20.12,4ES16.8)')  time*tScale,     & ! 1
            &                                       EperpT,EparT, EperpaxT,EparaxT
            if ( l_save_out ) close(n_perpPar_file)
         end if

         call Eperp%compute(EperpR_global, n_calls, timePassed, timeNorm)
         call Epar%compute(EparR_global, n_calls, timePassed, timeNorm)
         call Eperpaxi%compute(EperpaxiR_global, n_calls, timePassed, timeNorm)
         call Eparaxi%compute(EparaxiR_global, n_calls, timePassed, timeNorm)

         if ( l_stop_time ) then
            call Eperp%finalize_SD(timeNorm)
            call Epar%finalize_SD(timeNorm)
            call Eperpaxi%finalize_SD(timeNorm)
            call Eparaxi%finalize_SD(timeNorm)

            if (l_master_rank) then
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
      end if

   end subroutine outPerpPar
!----------------------------------------------------------------------------
end module outPar_mod
