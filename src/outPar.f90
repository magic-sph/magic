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
   use useful, only: cc2real
   use integration, only: rInt_R
#ifdef WITH_SHTNS
   use shtns, only: axi_to_spat
#else
   use legendre_spec_to_grid, only: lmAS2pt
#endif

   implicit none

   private

   real(cp), allocatable :: dlVMeanR(:),dlVcMeanR(:)
   real(cp), allocatable :: dlVu2MeanR(:),dlVu2cMeanR(:)
   real(cp), allocatable :: RolMeanR(:),RolMeanRu2(:),RmMeanR(:)
   real(cp), allocatable :: sMeanR(:),Svar(:),Mvar(:)
   real(cp), allocatable :: uhMeanR(:),duhMeanR(:)
   real(cp), allocatable :: gradT2MeanR(:)
   real(cp), allocatable :: fcondMeanR(:),fconvMeanR(:),fkinMeanR(:)
   real(cp), allocatable :: fviscMeanR(:)
   real(cp), allocatable :: fresMeanR(:), fpoynMeanR(:)
   real(cp), allocatable :: EperpMeanR(:),EparMeanR(:)
   real(cp), allocatable :: EperpaxiMeanR(:),EparaxiMeanR(:)
   integer :: n_perpPar_file
   character(len=72) :: perpPar_file

   public initialize_outPar_mod, finalize_outPar_mod, outPar, outPerpPar

contains

   subroutine initialize_outPar_mod

      allocate( dlVMeanR(n_r_max),dlVcMeanR(n_r_max) )
      allocate( dlVu2MeanR(n_r_max),dlVu2cMeanR(n_r_max) )
      allocate( RolMeanR(n_r_max),RolMeanRu2(n_r_max),RmMeanR(n_r_max) )
      bytes_allocated = bytes_allocated+7*n_r_max*SIZEOF_DEF_REAL

      dlVMeanR(:)     =0.0_cp
      dlVcMeanR(:)    =0.0_cp
      dlVu2MeanR(:)   =0.0_cp
      dlVu2cMeanR(:)  =0.0_cp
      RolMeanR(:)     =0.0_cp
      RolMeanRu2(:)   =0.0_cp
      RmMeanR(:)      =0.0_cp

      perpPar_file='perpPar.'//tag

      if ( l_viscBcCalc ) then
         allocate( sMeanR(n_r_max),Svar(nRstart:nRstop),Mvar(nRstart:nRstop) )
         allocate( uhMeanR(n_r_max),duhMeanR(n_r_max),gradT2MeanR(n_r_max) )
         bytes_allocated = bytes_allocated+6*n_r_max*SIZEOF_DEF_REAL
         Svar(:)         =0.0_cp
         Mvar(:)         =0.0_cp
         sMeanR(:)       =0.0_cp
         uhMeanR(:)      =0.0_cp
         duhMeanR(:)     =0.0_cp
         gradT2MeanR(:)  =0.0_cp
      end if

      if ( l_fluxProfs ) then
         allocate( fcondMeanR(n_r_max),fconvMeanR(n_r_max),fkinMeanR(n_r_max) )
         allocate( fviscMeanR(n_r_max) )
         bytes_allocated = bytes_allocated+4*n_r_max*SIZEOF_DEF_REAL
         fcondMeanR(:)   =0.0_cp
         fconvMeanR(:)   =0.0_cp
         fkinMeanR(:)    =0.0_cp
         fviscMeanR(:)   =0.0_cp
         allocate( fresMeanR(n_r_max),fpoynMeanR(n_r_max) )
         bytes_allocated = bytes_allocated+2*n_r_max*SIZEOF_DEF_REAL
         fresMeanR(:)    =0.0_cp
         fpoynMeanR(:)   =0.0_cp
      end if

      if ( l_perpPar ) then
         allocate( EperpMeanR(n_r_max),EparMeanR(n_r_max) )
         allocate( EperpaxiMeanR(n_r_max),EparaxiMeanR(n_r_max) )
         bytes_allocated = bytes_allocated+4*n_r_max*SIZEOF_DEF_REAL
         EperpMeanR(:)   =0.0_cp
         EparMeanR(:)    =0.0_cp
         EperpaxiMeanR(:)=0.0_cp
         EparaxiMeanR(:) =0.0_cp

         if ( l_master_rank .and. (.not. l_save_out) ) then
            open(newunit=n_perpPar_file, file=perpPar_file, status='new')
         end if
      end if

   end subroutine initialize_outPar_mod
!-----------------------------------------------------------------------
   subroutine finalize_outPar_mod

      deallocate( dlVMeanR, dlVcMeanR )
      deallocate( dlVu2MeanR, dlVu2cMeanR )
      deallocate( RolMeanR, RolMeanRu2,RmMeanR )

      if ( l_viscBcCalc ) then
         deallocate( sMeanR,Svar,Mvar )
         deallocate( uhMeanR,duhMeanR,gradT2MeanR )
      end if

      if ( l_fluxProfs ) then
         deallocate( fcondMeanR,fconvMeanR,fkinMeanR )
         deallocate( fviscMeanR )
      end if

      if ( l_perpPar ) then
         deallocate( EperpMeanR, EparMeanR )
         deallocate( EperpaxiMeanR, EparaxiMeanR )
         if ( l_master_rank .and. (.not. l_save_out) ) close(n_perpPar_file)
      end if

   end subroutine finalize_outPar_mod
!-----------------------------------------------------------------------
   subroutine outPar(timePassed,timeNorm,nLogs,l_stop_time,       &
              &      ekinR,RolRu2,dlVR,dlVRc,dlVRu2,dlVRu2c,      &
              &      uhLMr,duhLMr,gradsLMr,fconvLMr,fkinLMr,      &
              &      fviscLMr,fpoynLMr,fresLMr,RmR)

      !--- Input of variables
      real(cp), intent(in) :: timePassed,timeNorm
      logical,  intent(in) :: l_stop_time
      integer,  intent(in) :: nLogs
      real(cp), intent(in) :: RolRu2(n_r_max),dlVRu2(n_r_max),dlVRu2c(n_r_max)
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
      integer :: nR,n,m,lm
      real(cp) :: ReR(n_r_max), RoR(n_r_max), RolR(n_r_max)
      character(len=76) :: filename
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS
      real(cp) :: duhR(nRstart:nRstop), uhR(nRstart:nRstop)
      real(cp) :: gradT2R(nRstart:nRstop), sR(nRstart:nRstop), sR2(nRstart:nRstop)
      real(cp) :: fkinR(nRstart:nRstop), fcR(nRstart:nRstop)
      real(cp) :: fconvR(nRstart:nRstop), fviscR(nRstart:nRstop)
      real(cp) :: fresR(nRstartMag:nRstopMag),fpoynR(nRstartMag:nRstopMag)
      real(cp) :: duhR_global(n_r_max), uhR_global(n_r_max)
      real(cp) :: gradT2R_global(n_r_max), sR_global(n_r_max)
      real(cp) :: Svar_global(n_r_max)
      real(cp) :: fkinR_global(n_r_max), fcR_global(n_r_max)
      real(cp) :: fconvR_global(n_r_max), fviscR_global(n_r_max)
      real(cp) :: fresR_global(n_r_maxMag), fpoynR_global(n_r_maxMag)
      real(cp) :: duh(nfs), uh(nfs), gradT2(nfs)
      real(cp) :: fkin(nfs), fconv(nfs), fvisc(nfs), fres(nfs), fpoyn(nfs)

      integer :: fileHandle


      if ( l_viscBcCalc ) then
         do nR=nRstart,nRstop
            sR(nR) = real(s_Rloc(1,nR))
            ! calculate entropy/temperature variance:
            sR2(nR)=0.0_cp
            do lm=1,lm_max
              m=lm2m(lm)
              sR2(nR)=sR2(nR)+cc2real(s_Rloc(lm,nR),m)
            end do
            if ( nLogs  <=  1) then
               Mvar(nR)=sR(nR)
               Svar(nR)=sR2(nR)-sR(nR)**2
            else
               Mvar(nR)=Mvar(nR)+(sR(nR)-Mvar(nR))/nLogs
               Svar(nR)=Svar(nR)+(sR2(nR)-Mvar(nR)**2)
            end if
         end do

         do nR=nRstart,nRstop
            uhR(nR)    =0.0_cp
            gradT2R(nR)=0.0_cp
            duhR(nR)   =0.0_cp
#ifdef WITH_SHTNS
            call axi_to_spat(duhLMr(:,nR),duh)
            call axi_to_spat(uhLMr(:,nR),uh)
            call axi_to_spat(gradsLMr(:,nR),gradT2)
#endif
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
#ifndef WITH_SHTNS
               call lmAS2pt(duhLMr(:,nR),duh,nThetaStart,sizeThetaB)
               call lmAS2pt(uhLMr(:,nR),uh,nThetaStart,sizeThetaB)
               call lmAS2pt(gradsLMr(:,nR),gradT2,nThetaStart,sizeThetaB)
#endif
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
                  duhR(nR)=duhR(nR)+gauss(nThetaNHS)*duh(nThetaBlock)
                  uhR(nR) =uhR(nR) +gauss(nThetaNHS)* uh(nThetaBlock)
                  gradT2R(nR)=gradT2R(nR)+gauss(nThetaNHS)*gradT2(nThetaBlock)
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
         call gather_from_RLoc(Svar, Svar_global, 0)

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
            call axi_to_spat(fkinLMr(:,nR),fkin)
            call axi_to_spat(fconvLMr(:,nR),fconv)
            call axi_to_spat(fviscLMr(:,nR),fvisc)
#endif
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
#ifndef WITH_SHTNS
               call lmAS2pt(fkinLMr(:,nR),fkin,nThetaStart,sizeThetaB)
               call lmAS2pt(fconvLMr(:,nR),fconv,nThetaStart,sizeThetaB)
               call lmAS2pt(fviscLMr(:,nR),fvisc,nThetaStart,sizeThetaB)
#endif
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
                  fkinR(nR) =fkinR(nR) +gauss(nThetaNHS)* fkin(nThetaBlock)
                  fconvR(nR)=fconvR(nR)+gauss(nThetaNHS)*fconv(nThetaBlock)
                  fviscR(nR)=fviscR(nR)+gauss(nThetaNHS)*fvisc(nThetaBlock)
               end do
            end do
         end do

         if ( l_mag_nl ) then
            do nR=nRstart,nRstop
               fresR(nR) =0.0_cp
               fpoynR(nR)=0.0_cp
#ifdef WITH_SHTNS
               call axi_to_spat(fpoynLMr(:,nR),fpoyn)
               call axi_to_spat(fresLMr(:,nR),fres)
#endif
               do n=1,nThetaBs ! Loop over theta blocks
                  nTheta=(n-1)*sizeThetaB
                  nThetaStart=nTheta+1
#ifndef WITH_SHTNS
                  call lmAS2pt(fpoynLMr(:,nR),fpoyn,nThetaStart,sizeThetaB)
                  call lmAS2pt(fresLMr(:,nR),fres,nThetaStart,sizeThetaB)
#endif
                  do nThetaBlock=1,sizeThetaB
                     nTheta=nTheta+1
                     nThetaNHS=(nTheta+1)/2
                     fpoynR(nR)=fpoynR(nR)+gauss(nThetaNHS)*fpoyn(nThetaBlock)
                     fresR(nR) =fresR(nR) +gauss(nThetaNHS)*fres(nThetaBlock)
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

         dlVMeanR   =dlVMeanR   +timePassed*dlVR
         dlVcMeanR  =dlVcMeanR  +timePassed*dlVRc
         dlVu2MeanR =dlVu2MeanR +timePassed*dlVRu2
         dlVu2cMeanR=dlVu2cMeanR+timePassed*dlVRu2c
         RolMeanR   =RolMeanR   +timePassed*RolR
         RolMeanRu2 =RolMeanRu2 +timePassed*RolRu2
         RmMeanR    =RmMeanR    +timePassed*RmR*sqrt(mass/rho0)*or2
         !write(*,"(A,ES20.12)") "dlVcMeanR(n_r_icb) = ",dlVcMeanR(n_r_icb)
         ! this is to get u2 value for RmR(r) to plot in parR.tag
         ! and also remove r**2, so it has to be volume-averaged
         ! like RolR
         if ( l_viscBcCalc ) then
            sMeanR     =sMeanR    +timePassed*sR_global
            uhMeanR    =uhMeanR   +timePassed*uhR_global
            duhMeanR   =duhMeanR  +timePassed*duhR_global
            gradT2MeanR=gradT2MeanR+timePassed*gradT2R_global
         end if

         if ( l_fluxProfs ) then
            fkinMeanR =fkinMeanR  +timePassed*fkinR_global
            fcondMeanR=fcondMeanR +timePassed*fcR_global
            fconvMeanR=fconvMeanR +timePassed*fconvR_global
            fviscMeanR=fviscMeanR +timePassed*fviscR_global
            if ( l_mag_nl ) then
               fresMeanR =fresMeanR +timePassed*fresR_global
               fpoynMeanR=fpoynMeanR+timePassed*fpoynR_global
            end if
         end if

         if ( l_stop_time ) then

            dlVMeanR   =dlVMeanR/timeNorm
            dlVcMeanR  =dlVcMeanR/timeNorm
            RolMeanR   =RolMeanR/timeNorm
            if ( l_anel ) then
               dlVu2MeanR =dlVu2MeanR/timeNorm
               dlVu2cMeanR=dlVu2cMeanR/timeNorm
               RolMeanRu2 =RolMeanRu2/timeNorm
            else
               dlVu2MeanR =dlVMeanR
               dlVu2cMeanR=dlVcMeanR
               RolMeanRu2 =RolMeanR
            end if
            RmMeanR    =RmMeanR/timeNorm

            if ( l_viscBcCalc ) then
               sMeanR     =sMeanR/timeNorm
               Svar_global=Svar_global/(nLogs)
               if ( ktops == 1 ) Svar_global(1) = 0.0_cp
               if ( kbots == 1 ) Svar_global(n_r_max) = 0.0_cp
               duhMeanR   =duhMeanR/timeNorm
               uhMeanR    =uhMeanR/timeNorm
               gradT2MeanR=gradT2MeanR/timeNorm
            end if

            if ( l_fluxProfs ) then
               fkinMeanR =ViscHeatFac*fkinMeanR/timeNorm
               fcondMeanR=opr*fcondMeanR/timeNorm
               fconvMeanR=fconvMeanR/timeNorm
               fviscMeanR=ViscHeatFac*fviscMeanR/timeNorm
               if ( l_mag_nl ) then
                  fresMeanR =OhmLossFac*fresMeanR/timeNorm
                  fpoynMeanR=prmag*OhmLossFac*fpoynMeanR/timeNorm
               end if
            end if

            !----- Output into paR.TAG file:
            filename='parR.'//tag
            open(newunit=fileHandle, file=filename, status='unknown')
            do nR=1,n_r_max
               write(fileHandle,'(ES20.10,8ES15.7)') &
               &              r(nR),                 &! 1) radius
               &              RmMeanR(nR),           &! 2) magnetic Reynolds number
               &              RolMeanR(nR),          &! 3) local Rossby number
               &              RolMeanRu2(nR),        &! 4) u squared local Rossby number
               &              dlVMeanR(nR),          &! 5) local length scale
               &              dlVcMeanR(nR),         &! 6) conv. local length scale
               &              dlVu2MeanR(nR),        &! 7) u squared local length scale
               &              dlVu2cMeanR(nR)         ! 8) u squared conv. local length scale
            end do
            close(fileHandle)

            if ( l_viscBcCalc ) then
               filename='bLayersR.'//tag
               open(newunit=fileHandle, file=filename, status='unknown')
               do nR=1,n_r_max
                  write(fileHandle,'(ES20.10,6ES15.7)')    &
                  &           r(nR),                       &! 1) radius
                  &           sMeanR(nR)*osq4pi,           &! 2) entropy
                  &           Svar_global(nR)/(four*pi),   &! 3) entropy variance
                  &           uhMeanR(nR),                 &! 4) uh
                  &           duhMeanR(nR),                &! 5) duh/dr
                  &           gradT2MeanR(nR)               ! 6) (grad T)**2
               end do
               close(fileHandle)
            end if

            if ( l_fluxProfs ) then
               filename='fluxesR.'//tag
               open(newunit=fileHandle, file=filename, status='unknown')
               do nR=1,n_r_max
                  write(fileHandle,'(ES20.10,7ES15.7)')  &
                  &           r(nR),                     &! 1) radius
                  &           fcondMeanR(nR),            &! 2) Fcond
                  &           fconvMeanR(nR),            &! 3) Fconv
                  &           fkinMeanR(nR),             &! 4) Fkin
                  &           fviscMeanR(nR),            &! 5) Fvisc
                  &           fpoynMeanR(nR),            &! 6) Fpoyn
                  &           fresMeanR(nR)               ! 7) Fres
               end do
               close(fileHandle)
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
      real(cp) :: Eperp(nfs), Epar(nfs), Eperpaxi(nfs), Eparaxi(nfs)
      real(cp) :: EperpT,EparT,EperpaxT,EparaxT

      integer :: fileHandle

      do nR=nRstart,nRstop
         EperpR(nR)   =0.0_cp
         EparR(nR)    =0.0_cp
         EparaxiR(nR) =0.0_cp
         EperpaxiR(nR)=0.0_cp
#ifdef WITH_SHTNS
         call axi_to_spat(EperpLMr(:,nR),Eperp)
         call axi_to_spat(EparLMr(:,nR),Epar)
         call axi_to_spat(EperpaxiLMr(:,nR),Eperpaxi)
         call axi_to_spat(EparaxiLMr(:,nR),Eparaxi)
#endif
         do n=1,nThetaBs ! Loop over theta blocks
            nTheta=(n-1)*sizeThetaB
            nThetaStart=nTheta+1
#ifndef WITH_SHTNS
            call lmAS2pt(EperpLMr(:,nR),Eperp,nThetaStart,sizeThetaB)
            call lmAS2pt(EparLMr(:,nR),Epar,nThetaStart,sizeThetaB)
            call lmAS2pt(EperpaxiLMr(:,nR),Eperpaxi,nThetaStart,sizeThetaB)
            call lmAS2pt(EparaxiLMr(:,nR),Eparaxi,nThetaStart,sizeThetaB)
#endif
            do nThetaBlock=1,sizeThetaB
               nTheta=nTheta+1
               nThetaNHS=(nTheta+1)/2
               EperpR(nR)=EperpR(nR)+gauss(nThetaNHS)*Eperp(nThetaBlock)
               EparR(nR) =EparR(nR) +gauss(nThetaNHS)* Epar(nThetaBlock)
               EperpaxiR(nR)=EperpaxiR(nR)+gauss(nThetaNHS)*Eperpaxi(nThetaBlock)
               EparaxiR(nR)=EparaxiR(nR)+gauss(nThetaNHS)*Eparaxi(nThetaBlock)
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
         if ( l_save_out ) then
            open(newunit=n_perpPar_file, file=perpPar_file, &
            &    status='unknown', position='append')
         end if
         write(n_perpPar_file,'(1P,ES20.12,4ES16.8)')  time*tScale,     & ! 1
         &                                       EperpT,EparT, EperpaxT,EparaxT
         if ( l_save_out ) close(n_perpPar_file)

         EperpMeanR    =EperpMeanR     +timePassed*EperpR_global
         EparMeanR     =EparMeanR      +timePassed*EparR_global
         EperpaxiMeanR =EperpaxiMeanR  +timePassed*EperpaxiR_global
         EparaxiMeanR  =EparaxiMeanR   +timePassed*EparaxiR_global
         if ( l_stop_time ) then
             EperpMeanR     =EperpMeanR/timeNorm
             EparMeanR      =EparMeanR/timeNorm
             EperpaxiMeanR  =EperpaxiMeanR/timeNorm
             EparaxiMeanR   =EparaxiMeanR/timeNorm
             filename='perpParR.'//tag
             open(newunit=fileHandle, file=filename, status='unknown')
             do nR=1,n_r_max
                write(fileHandle,'(ES20.10,4ES15.7)')&
                &              r(nR),                &! 1) radius
                &              EperpMeanR(nR),       &! 2) E perpendicular
                &              EparMeanR(nR),        &! 3) E parallel
                &              EperpaxiMeanR(nR),    &! 4) E perp (axisymetric)
                &              EparaxiMeanR(nR)       ! 5) E parallel (axisymetric)
             end do
             close(fileHandle)
         end if
      end if

   end subroutine outPerpPar
!----------------------------------------------------------------------------
end module outPar_mod
