module outMisc_mod
   !
   ! This module contains several subroutines that can compute and store
   ! various informations: helicity, heat transfer.
   !

   use parallel_mod
   use precision_mod
   use communications, only: gather_from_Rloc
   use truncation, only: l_max, n_r_max, lm_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: r_icb, rscheme_oc, kappa,         &
       &                       r_cmb,temp0, r, rho0, dLtemp0,    &
       &                       dLalpha0, beta, orho1, alpha0,    &
       &                       otemp1, ogrun, rscheme_oc
   use physical_parameters, only: ViscHeatFac, ThExpNb, oek, GrunNb
   use num_param, only: lScale
   use blocking, only: nThetaBs, nfs, sizeThetaB, llm, ulm, lo_map
   use mean_sd, only: mean_sd_type
   use horizontal_data, only: gauss
   use logic, only: l_save_out, l_anelastic_liquid, l_heat, l_hel, &
       &            l_temperature_diff, l_chemical_conv, l_mag_hel,&
       &            l_gw
   use output_data, only: tag
   use constants, only: pi, vol_oc, osq4pi, sq4pi, one, two, four, ci, &
        &               y20_norm, y21_norm, y22_norm
   use start_fields, only: topcond, botcond, deltacond, topxicond, botxicond, &
       &                   deltaxicond
   use useful, only: cc2real, round_off
   use integration, only: rInt_R
#ifdef WITH_SHTNS
   use shtns, only: axi_to_spat
#else
   use legendre_spec_to_grid, only: lmAS2pt
#endif

   implicit none

   private

   type(mean_sd_type) :: TMeanR, SMeanR, PMeanR, XiMeanR,RhoMeanR
   integer :: n_heat_file, n_helicity_file, n_magHel_file, n_calls
   integer :: n_gw_S_file, n_gw_P_file
   character(len=72) :: heat_file, helicity_file, magHel_file, &
        &               gw_S_file, gw_P_file
   public :: outHelicity, outHeat, initialize_outMisc_mod, finalize_outMisc_mod
   public :: outMagneticHelicity, outGWentropy, outGWpressure

contains

   subroutine initialize_outMisc_mod

      if (l_heat .or. l_chemical_conv) then
         call TMeanR%initialize(1,n_r_max)
         call SMeanR%initialize(1,n_r_max)
         call PMeanR%initialize(1,n_r_max)
         call XiMeanR%initialize(1,n_r_max)
         call RhoMeanR%initialize(1,n_r_max)
      endif
      n_calls = 0

      helicity_file='helicity.'//tag
      heat_file    ='heat.'//tag
      magHel_file  ='helicity_mag.'//tag
      gw_P_file    ='gwPressure.'//tag
      gw_S_file    ='gwEntropy.'//tag
      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_hel ) then
            open(newunit=n_helicity_file, file=helicity_file, status='new')
         end if
         if ( l_heat .or. l_chemical_conv ) then
            open(newunit=n_heat_file, file=heat_file, status='new')
         end if
         if (l_mag_hel) then
            open(newunit=n_magHel_file, file=magHel_file, status='new')
         end if
         if ( l_gw ) then
            open(newunit=n_gw_P_file, file=gw_P_file, &
                 status='new',form='unformatted')
            if ( l_heat .or. l_chemical_conv ) then
               open(newunit=n_gw_S_file, file=gw_S_file, &
                    status='new',form='unformatted')
            end if
         endif
      end if
    end subroutine initialize_outMisc_mod
!---------------------------------------------------------------------------
   subroutine finalize_outMisc_mod

      if ( l_heat .or. l_chemical_conv ) then
         !deallocate( TMeanR, SMeanR, PMeanR, XiMeanR )
         call TMeanR%finalize()
         call SMeanR%finalize()
         call PMeanR%finalize()
         call XiMeanR%finalize()
         call RhoMeanR%finalize()
      end if

      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_hel ) close(n_helicity_file)
         if ( l_mag_hel ) close(n_magHel_file)
         if ( l_heat .or. l_chemical_conv ) close(n_heat_file)
         if (l_gw) then
            close(n_gw_P_file)
            if ( l_heat .or. l_chemical_conv ) close(n_gw_S_file)
         end if
      end if
   end subroutine finalize_outMisc_mod
!---------------------------------------------------------------------------
   subroutine outHelicity(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr)
      !
      ! This subroutine is used to store informations about kinetic
      ! helicity
      !

      !-- Input of variables:
      real(cp), intent(in) :: timeScaled
      real(cp), intent(in) :: HelLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: Helna2LMr(l_max+1,nRstart:nRstop)

      !-- Local stuff:
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS,n
      real(cp) :: HelNr(nRstart:nRstop), HelSr(nRstart:nRstop)
      real(cp) :: HelnaNr(nRstart:nRstop), HelnaSr(nRstart:nRstop)
      real(cp) :: Hel2Nr(nRstart:nRstop), Hel2Sr(nRstart:nRstop)
      real(cp) :: Helna2Nr(nRstart:nRstop), Helna2Sr(nRstart:nRstop)
      real(cp) :: HelEAr(nRstart:nRstop)
      real(cp) :: HelNr_global(n_r_max), HelSr_global(n_r_max)
      real(cp) :: HelnaNr_global(n_r_max), HelnaSr_global(n_r_max)
      real(cp) :: Helna2Nr_global(n_r_max), Helna2Sr_global(n_r_max)
      real(cp) :: Hel2Nr_global(n_r_max), Hel2Sr_global(n_r_max)
      real(cp) :: HelEAr_global(n_r_max)
      real(cp) :: Hel(nfs), Hel2(nfs), Helna(nfs), Helna2(nfs), r2
      real(cp) :: HelN,HelS
      real(cp) :: HelnaN,HelnaS
      real(cp) :: HelnaRMSN,HelnaRMSS
      real(cp) :: HelRMSN,HelRMSS,HelEA,HelRMS,HelnaRMS

      integer :: n_r

      !------ Integration of Helicity, on input the Helicity is
      !       already axisymmetric !
      do n_r=nRstart,nRstop
         r2=r(n_r)*r(n_r)
         HelNr(n_r)   =0.0_cp
         HelSr(n_r)   =0.0_cp
         HelnaNr(n_r) =0.0_cp
         HelnaSr(n_r) =0.0_cp
         HelEAr(n_r)  =0.0_cp
         Hel2Nr(n_r)  =0.0_cp
         Hel2Sr(n_r)  =0.0_cp
         Helna2Nr(n_r)=0.0_cp
         Helna2Sr(n_r)=0.0_cp

#ifdef WITH_SHTNS
         call axi_to_spat(HelLMr(:,n_r), Hel)
         call axi_to_spat(Hel2LMr(:,n_r), Hel2)
         call axi_to_spat(HelnaLMr(:,n_r), Helna)
         call axi_to_spat(Helna2LMr(:,n_r), Helna2)
#endif

         do n=1,nThetaBs ! Loop over theta blocks
            nTheta=(n-1)*sizeThetaB
            nThetaStart=nTheta+1
#ifndef WITH_SHTNS
            call lmAS2pt(HelLMr(:,n_r),Hel,nThetaStart,sizeThetaB)
            call lmAS2pt(Hel2LMr(:,n_r),Hel2,nThetaStart,sizeThetaB)
            call lmAS2pt(HelnaLMr(:,n_r),Helna,nThetaStart,sizeThetaB)
            call lmAS2pt(Helna2LMr(:,n_r),Helna2,nThetaStart,sizeThetaB)
#endif
            do nThetaBlock=1,sizeThetaB
               nTheta=nTheta+1
               nThetaNHS=(nTheta+1)/2

               !------ Integration over theta:
               if ( mod(nTheta,2) == 1 ) then ! NHS
                  Hel2Nr(n_r)=Hel2Nr(n_r)+gauss(nThetaNHS)*r2*Hel2(nThetaBlock)
                  Helna2Nr(n_r)=Helna2Nr(n_r)+gauss(nThetaNHS)*r2*Helna2(nThetaBlock)
                  HelEAr(n_r)=HelEAr(n_r)+gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                  HelNr(n_r) =HelNr(n_r)+gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                  HelnaNr(n_r) =HelnaNr(n_r)+gauss(nThetaNHS)*r2*Helna(nThetaBlock)
               else
                  Hel2Sr(n_r)=Hel2Sr(n_r)+gauss(nThetaNHS)*r2*Hel2(nThetaBlock)
                  Helna2Sr(n_r)=Helna2Sr(n_r)+gauss(nThetaNHS)*r2*Helna2(nThetaBlock)
                  HelEAr(n_r)=HelEAr(n_r)-gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                  HelSr(n_r) =HelSr(n_r)+gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                  HelnaSr(n_r)=HelnaSr(n_r)+gauss(nThetaNHS)*r2*Helna(nThetaBlock)
               end if
            end do
         end do

      end do

      ! Now we have to gather the results on rank 0 for
      ! the arrays: Hel2Nr,Helna2Nr,HelEAr,HelNr,HelnaNr
      ! Hel2Sr,Helna2Sr,HelSr,HelnaSr

      call gather_from_Rloc(Hel2Nr, Hel2Nr_global, 0)
      call gather_from_Rloc(Helna2Nr, Helna2Nr_global, 0)
      call gather_from_Rloc(HelEAr, HelEAr_global, 0)
      call gather_from_Rloc(HelNr, HelNr_global, 0)
      call gather_from_Rloc(HelnaNr, HelnaNr_global, 0)
      call gather_from_Rloc(HelSr, HelSr_global, 0)
      call gather_from_Rloc(Helna2Sr, Helna2Sr_global, 0)
      call gather_from_Rloc(Hel2Sr, Hel2Sr_global, 0)
      call gather_from_Rloc(HelnaSr, HelnaSr_global, 0)

      if ( rank == 0 ) then
         !------ Integration over r without the boundaries and normalization:
         HelN  =rInt_R(HelNr_global,r,rscheme_oc)
         HelS  =rInt_R(HelSr_global,r,rscheme_oc)
         HelnaN=rInt_R(HelnaNr_global,r,rscheme_oc)
         HelnaS=rInt_R(HelnaSr_global,r,rscheme_oc)
         HelEA =rInt_R(HelEAr_global,r,rscheme_oc)
         HelRMSN=rInt_R(Hel2Nr_global,r,rscheme_oc)
         HelRMSS=rInt_R(Hel2Sr_global,r,rscheme_oc)
         HelnaRMSN=rInt_R(Helna2Nr_global,r,rscheme_oc)
         HelnaRMSS=rInt_R(Helna2Sr_global,r,rscheme_oc)

         HelN  =two*pi*HelN/(vol_oc/2) ! Note integrated over half spheres only !
         HelS  =two*pi*HelS/(vol_oc/2) ! Factor 2*pi is from phi integration
         HelnaN=two*pi*HelnaN/(vol_oc/2) ! Note integrated over half spheres only !
         HelnaS=two*pi*HelnaS/(vol_oc/2) ! Factor 2*pi is from phi integration
         HelEA =two*pi*HelEA/vol_oc
         HelRMSN=sqrt(two*pi*HelRMSN/(vol_oc/2))
         HelRMSS=sqrt(two*pi*HelRMSS/(vol_oc/2))
         HelnaRMSN=sqrt(two*pi*HelnaRMSN/(vol_oc/2))
         HelnaRMSS=sqrt(two*pi*HelnaRMSS/(vol_oc/2))
         HelRMS=HelRMSN+HelRMSS
         HelnaRMS=HelnaRMSN+HelnaRMSS

         if ( HelnaRMS /= 0 ) then
            HelnaN =HelnaN/HelnaRMSN
            HelnaS =HelnaS/HelnaRMSS
         else
            HelnaN =0.0_cp
            HelnaS =0.0_cp
         end if
         if ( HelRMS /= 0 ) then
            HelN =HelN/HelRMSN
            HelS =HelS/HelRMSS
            HelEA=HelEA/HelRMS
         else
            HelN =0.0_cp
            HelS =0.0_cp
            HelEA=0.0_cp
         end if

         if ( l_save_out ) then
            open(newunit=n_helicity_file, file=helicity_file,   &
            &    status='unknown', position='append')
         end if

         write(n_helicity_file,'(1P,ES20.12,8ES16.8)')   &
         &     timeScaled,HelN, HelS, HelRMSN, HelRMSS,  &
         &     HelnaN, HelnaS, HelnaRMSN, HelnaRMSS

         if ( l_save_out ) close(n_helicity_file)

      end if

    end subroutine outHelicity
!---------------------------------------------------------------------------
    subroutine outGWentropy(timeScaled,s)
      !-- Input of variables:
      real(cp),    intent(in) :: timeScaled
      !-- Input of scalar fields:
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)

      !-- Local stuff:
      integer  :: nR, ierr
      integer  :: lm20,lm21,lm22, lm, l, m
      real(cp) :: r4, prefactor
      !-- entropy contributions to density fluctuations
      real(cp) :: Qc_entropy_20_r(n_r_max)
      real(cp) :: Qc_entropy_21_r(n_r_max)
      real(cp) :: Qs_entropy_21_r(n_r_max)
      real(cp) :: Qc_entropy_22_r(n_r_max)
      real(cp) :: Qs_entropy_22_r(n_r_max)

      real(cp) :: Qc_entropy_20_r_global(n_r_max)
      real(cp) :: Qc_entropy_21_r_global(n_r_max)
      real(cp) :: Qs_entropy_21_r_global(n_r_max)
      real(cp) :: Qc_entropy_22_r_global(n_r_max)
      real(cp) :: Qs_entropy_22_r_global(n_r_max)

      real(cp) :: Qc_entropy_20
      real(cp) :: Qc_entropy_21
      real(cp) :: Qs_entropy_21
      real(cp) :: Qc_entropy_22
      real(cp) :: Qs_entropy_22

      !-- phi derivative
      real(cp) :: dPhiQc_entropy_21_r(n_r_max)
      real(cp) :: dPhiQs_entropy_21_r(n_r_max)
      real(cp) :: dPhiQc_entropy_22_r(n_r_max)
      real(cp) :: dPhiQs_entropy_22_r(n_r_max)
      real(cp) :: dPhiQc_entropy_21_r_global(n_r_max)
      real(cp) :: dPhiQs_entropy_21_r_global(n_r_max)
      real(cp) :: dPhiQc_entropy_22_r_global(n_r_max)
      real(cp) :: dPhiQs_entropy_22_r_global(n_r_max)
      real(cp) :: dPhiQc_entropy_21
      real(cp) :: dPhiQs_entropy_21
      real(cp) :: dPhiQc_entropy_22
      real(cp) :: dPhiQs_entropy_22
      !-- 2nd phi derivative
      real(cp) :: ddPhiQc_entropy_21_r(n_r_max)
      real(cp) :: ddPhiQs_entropy_21_r(n_r_max)
      real(cp) :: ddPhiQc_entropy_22_r(n_r_max)
      real(cp) :: ddPhiQs_entropy_22_r(n_r_max)
      real(cp) :: ddPhiQc_entropy_21_r_global(n_r_max)
      real(cp) :: ddPhiQs_entropy_21_r_global(n_r_max)
      real(cp) :: ddPhiQc_entropy_22_r_global(n_r_max)
      real(cp) :: ddPhiQs_entropy_22_r_global(n_r_max)
      real(cp) :: ddPhiQc_entropy_21
      real(cp) :: ddPhiQs_entropy_21
      real(cp) :: ddPhiQc_entropy_22
      real(cp) :: ddPhiQs_entropy_22

      !-- quadrupole indexes
      lm20 = lo_map%lm2(2,0)
      lm21 = lo_map%lm2(2,1)
      lm22 = lo_map%lm2(2,2)

      !-- radial loop
      do nR=1,n_r_max
         r4 = r(nR)**4
         Qc_entropy_20_r(nR)=0.0
         Qc_entropy_21_r(nR)=0.0
         Qs_entropy_21_r(nR)=0.0
         Qc_entropy_22_r(nR)=0.0
         Qs_entropy_22_r(nR)=0.0
         !-- dphi
         dPhiQc_entropy_21_r(nR)=0.0
         dPhiQs_entropy_21_r(nR)=0.0
         dPhiQc_entropy_22_r(nR)=0.0
         dPhiQs_entropy_22_r(nR)=0.0
         !-- ddphi
         ddPhiQc_entropy_21_r(nR)=0.0
         ddPhiQs_entropy_21_r(nR)=0.0
         ddPhiQc_entropy_22_r(nR)=0.0
         ddPhiQs_entropy_22_r(nR)=0.0

         if ( l_anelastic_liquid .or. l_TP_form ) then
            ! rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &
            !      &               real(s(1,n_r))+ViscHeatFac*(ThExpNb*     &
            !      &               alpha0(n_r)*temp0(n_r)+ogrun(n_r))*      &
            !      &               real(p(1,n_r)) )
            prefactor = -ThExpNb*alpha0(nR)*rho0(nR)
         else
            ! rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &
            !      &               temp0(n_r)*real(s(1,n_r))+ViscHeatFac*   &
            !      &               ogrun(n_r)*real(p(1,n_r)) )
            prefactor = -ThExpNb*alpha0(nR)*rho0(nR)*temp0(nR)
         end if

         ! Qc_entropy_20_r(nR) = y20_norm*r2*prefactor* real(s(lm20,nR))
         ! Qc_entropy_21_r(nR) = y21_norm*r2*prefactor* real(s(lm21,nR))
         ! Qc_entropy_22_r(nR) = y22_norm*r2*prefactor* real(s(lm22,nR))
         ! Qs_entropy_21_r(nR) = y21_norm*r2*prefactor*aimag(s(lm21,nR))
         ! Qs_entropy_22_r(nR) = y22_norm*r2*prefactor*aimag(s(lm22,nR))
         ! !-- dphi
         ! dPhiQc_entropy_21_r(nR)  = y21_norm*r2*prefactor* real(  ci*s(lm21,nR))
         ! dPhiQc_entropy_22_r(nR)  = y22_norm*r2*prefactor* real(2*ci*s(lm22,nR))
         ! dPhiQs_entropy_21_r(nR)  = y21_norm*r2*prefactor*aimag(  ci*s(lm21,nR))
         ! dPhiQs_entropy_22_r(nR)  = y22_norm*r2*prefactor*aimag(2*ci*s(lm22,nR))
         ! !-- ddphi
         ! ddPhiQc_entropy_21_r(nR) = y21_norm*r2*prefactor* real(  -s(lm21,nR))
         ! ddPhiQc_entropy_22_r(nR) = y22_norm*r2*prefactor* real(-4*s(lm22,nR))
         ! ddPhiQs_entropy_21_r(nR) = y21_norm*r2*prefactor*aimag(  -s(lm21,nR))
         ! ddPhiQs_entropy_22_r(nR) = y22_norm*r2*prefactor*aimag(-4*s(lm22,nR))
         do lm=max(2,llm), ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if (l == 2) then
               if (m == 0) then
                  Qc_entropy_20_r(nR) = y20_norm*r4*prefactor* real(s(lm,nR))
               else if (m == 1) then
                  Qc_entropy_21_r(nR) = y21_norm*r4*prefactor* real(s(lm,nR))
                  Qs_entropy_21_r(nR) = y21_norm*r4*prefactor*aimag(s(lm,nR))

                  dPhiQc_entropy_21_r(nR)  = y21_norm*r4*prefactor* real(  ci*s(lm,nR))
                  dPhiQs_entropy_21_r(nR)  = y21_norm*r4*prefactor*aimag(  ci*s(lm,nR))

                  ddPhiQc_entropy_21_r(nR) = y21_norm*r4*prefactor* real(  -s(lm,nR))
                  ddPhiQs_entropy_21_r(nR) = y21_norm*r4*prefactor*aimag(  -s(lm,nR))
               else if (m == 2) then
                  Qc_entropy_22_r(nR) = y22_norm*r4*prefactor* real(s(lm,nR))
                  Qs_entropy_22_r(nR) = y22_norm*r4*prefactor*aimag(s(lm,nR))

                  dPhiQc_entropy_22_r(nR)  = y22_norm*r4*prefactor* real(2*ci*s(lm,nR))
                  dPhiQs_entropy_22_r(nR)  = y22_norm*r4*prefactor*aimag(2*ci*s(lm,nR))

                  ddPhiQc_entropy_22_r(nR) = y22_norm*r4*prefactor* real(-4*s(lm,nR))
                  ddPhiQs_entropy_22_r(nR) = y22_norm*r4*prefactor*aimag(-4*s(lm,nR))
               end if
            end if
         end do
      end do

      ! reduce over the ranks
#ifdef WITH_MPI
      call MPI_Reduce(Qc_entropy_20_r, Qc_entropy_20_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qc_entropy_21_r, Qc_entropy_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qc_entropy_22_r, Qc_entropy_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qs_entropy_21_r, Qs_entropy_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qs_entropy_22_r, Qs_entropy_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_Reduce(dPhiQc_entropy_21_r, dPhiQc_entropy_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dPhiQc_entropy_22_r, dPhiQc_entropy_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dPhiQs_entropy_21_r, dPhiQs_entropy_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dPhiQs_entropy_22_r, dPhiQs_entropy_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_Reduce(ddPhiQc_entropy_21_r, ddPhiQc_entropy_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ddPhiQc_entropy_22_r, ddPhiQc_entropy_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ddPhiQs_entropy_21_r, ddPhiQs_entropy_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ddPhiQs_entropy_22_r, ddPhiQs_entropy_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

#else
      Qc_entropy_20_r_global = Qc_entropy_20_r
      Qc_entropy_21_r_global = Qc_entropy_21_r
      Qc_entropy_22_r_global = Qc_entropy_22_r
      Qs_entropy_21_r_global = Qs_entropy_21_r
      Qs_entropy_22_r_global = Qs_entropy_22_r

      dPhiQc_entropy_21_r_global = dPhiQc_entropy_21_r
      dPhiQc_entropy_22_r_global = dPhiQc_entropy_22_r
      dPhiQs_entropy_21_r_global = dPhiQs_entropy_21_r
      dPhiQs_entropy_22_r_global = dPhiQs_entropy_22_r

      ddPhiQc_entropy_21_r_global = ddPhiQc_entropy_21_r
      ddPhiQc_entropy_22_r_global = ddPhiQc_entropy_22_r
      ddPhiQs_entropy_21_r_global = ddPhiQs_entropy_21_r
      ddPhiQs_entropy_22_r_global = ddPhiQs_entropy_22_r
#endif

      if ( rank == 0 ) then
         !-- Radial Integrals:
         Qc_entropy_20 = rInt_R(Qc_entropy_20_r_global,r,rscheme_oc)
         Qc_entropy_21 = rInt_R(Qc_entropy_21_r_global,r,rscheme_oc)
         Qc_entropy_22 = rInt_R(Qc_entropy_22_r_global,r,rscheme_oc)
         Qs_entropy_21 = rInt_R(Qs_entropy_21_r_global,r,rscheme_oc)
         Qs_entropy_22 = rInt_R(Qs_entropy_22_r_global,r,rscheme_oc)

         dPhiQc_entropy_21 = 2*oek * rInt_R(dPhiQc_entropy_21_r_global,r,rscheme_oc)
         dPhiQc_entropy_22 = 2*oek * rInt_R(dPhiQc_entropy_22_r_global,r,rscheme_oc)
         dPhiQs_entropy_21 = 2*oek * rInt_R(dPhiQs_entropy_21_r_global,r,rscheme_oc)
         dPhiQs_entropy_22 = 2*oek * rInt_R(dPhiQs_entropy_22_r_global,r,rscheme_oc)

         ddPhiQc_entropy_21 = oek**2 * rInt_R(ddPhiQc_entropy_21_r_global,r,rscheme_oc)
         ddPhiQc_entropy_22 = oek**2 * rInt_R(ddPhiQc_entropy_22_r_global,r,rscheme_oc)
         ddPhiQs_entropy_21 = oek**2 * rInt_R(ddPhiQs_entropy_21_r_global,r,rscheme_oc)
         ddPhiQs_entropy_22 = oek**2 * rInt_R(ddPhiQs_entropy_22_r_global,r,rscheme_oc)

         !-- Write outputs
         if ( l_save_out ) then
            open(newunit=n_gw_S_file, file=gw_S_file,      &
                 &    status='unknown', position='append', &
                 &    form='unformatted')
         end if
         write(n_gw_S_file)  timeScaled,                 & ! 1
              &      Qc_entropy_20,                      & ! 2
              &      Qc_entropy_21,      Qs_entropy_21,  & ! 3,4
              &      Qc_entropy_22,      Qs_entropy_22,  & ! 5,6
              &  dPhiQc_entropy_21,  dPhiQs_entropy_21,  & ! 7,8
              &  dPhiQc_entropy_22,  dPhiQs_entropy_22,  & ! 9,10
              & ddPhiQc_entropy_21, ddPhiQs_entropy_21,  & ! 11,12
              & ddPhiQc_entropy_22, ddPhiQs_entropy_22     ! 13,14

         if ( l_save_out ) close(n_gw_S_file)
      end if

    end subroutine outGWentropy

    subroutine outGWpressure(timeScaled,p)
      !
      ! This subroutine is used to compute the coefficient
      ! that appear in the quadrupole formula describing
      ! the gravitational wave signal due to density fluctuations
      !

      !-- Input of variables:
      real(cp),    intent(in) :: timeScaled
      !-- Input of scalar fields:
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)

      !-- Local stuff:
      integer  :: nR, ierr
      integer  :: lm20,lm21,lm22, lm, l, m
      real(cp) :: r4, prefactor
      !-- pressure contributions to density fluctuations
      real(cp) :: Qc_pressure_20_r(n_r_max)
      real(cp) :: Qc_pressure_21_r(n_r_max)
      real(cp) :: Qs_pressure_21_r(n_r_max)
      real(cp) :: Qc_pressure_22_r(n_r_max)
      real(cp) :: Qs_pressure_22_r(n_r_max)

      real(cp) :: Qc_pressure_20_r_global(n_r_max)
      real(cp) :: Qc_pressure_21_r_global(n_r_max)
      real(cp) :: Qs_pressure_21_r_global(n_r_max)
      real(cp) :: Qc_pressure_22_r_global(n_r_max)
      real(cp) :: Qs_pressure_22_r_global(n_r_max)

      real(cp) :: Qc_pressure_20
      real(cp) :: Qc_pressure_21
      real(cp) :: Qs_pressure_21
      real(cp) :: Qc_pressure_22
      real(cp) :: Qs_pressure_22

      !-- phi derivative
      real(cp) :: dPhiQc_pressure_21_r(n_r_max)
      real(cp) :: dPhiQs_pressure_21_r(n_r_max)
      real(cp) :: dPhiQc_pressure_22_r(n_r_max)
      real(cp) :: dPhiQs_pressure_22_r(n_r_max)
      real(cp) :: dPhiQc_pressure_21_r_global(n_r_max)
      real(cp) :: dPhiQs_pressure_21_r_global(n_r_max)
      real(cp) :: dPhiQc_pressure_22_r_global(n_r_max)
      real(cp) :: dPhiQs_pressure_22_r_global(n_r_max)
      real(cp) :: dPhiQc_pressure_21
      real(cp) :: dPhiQs_pressure_21
      real(cp) :: dPhiQc_pressure_22
      real(cp) :: dPhiQs_pressure_22
      !-- 2nd phi derivative
      real(cp) :: ddPhiQc_pressure_21_r(n_r_max)
      real(cp) :: ddPhiQs_pressure_21_r(n_r_max)
      real(cp) :: ddPhiQc_pressure_22_r(n_r_max)
      real(cp) :: ddPhiQs_pressure_22_r(n_r_max)
      real(cp) :: ddPhiQc_pressure_21_r_global(n_r_max)
      real(cp) :: ddPhiQs_pressure_21_r_global(n_r_max)
      real(cp) :: ddPhiQc_pressure_22_r_global(n_r_max)
      real(cp) :: ddPhiQs_pressure_22_r_global(n_r_max)
      real(cp) :: ddPhiQc_pressure_21
      real(cp) :: ddPhiQs_pressure_21
      real(cp) :: ddPhiQc_pressure_22
      real(cp) :: ddPhiQs_pressure_22

      !-- quadrupole indexes
      lm20 = lo_map%lm2(2,0)
      lm21 = lo_map%lm2(2,1)
      lm22 = lo_map%lm2(2,2)

      !-- radial loop
      do nR=1,n_r_max
         r4 = r(nR)**4
         Qc_pressure_20_r(nR)=0.0
         Qc_pressure_21_r(nR)=0.0
         Qs_pressure_21_r(nR)=0.0
         Qc_pressure_22_r(nR)=0.0
         Qs_pressure_22_r(nR)=0.0
         !-- dphi
         dPhiQc_pressure_21_r(nR)=0.0
         dPhiQs_pressure_21_r(nR)=0.0
         dPhiQc_pressure_22_r(nR)=0.0
         dPhiQs_pressure_22_r(nR)=0.0
         !-- ddphi
         ddPhiQc_pressure_21_r(nR)=0.0
         ddPhiQs_pressure_21_r(nR)=0.0
         ddPhiQc_pressure_22_r(nR)=0.0
         ddPhiQs_pressure_22_r(nR)=0.0

         if ( l_anelastic_liquid .or. l_TP_form ) then
            ! rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &
            !      &               real(s(1,n_r))+ViscHeatFac*(ThExpNb*     &
            !      &               alpha0(n_r)*temp0(n_r)+ogrun(n_r))*      &
            !      &               real(p(1,n_r)) )
            prefactor = ThExpNb*alpha0(nR)*ViscHeatFac*(ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))
         else
            ! rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &
            !      &               temp0(n_r)*real(s(1,n_r))+ViscHeatFac*   &
            !      &               ogrun(n_r)*real(p(1,n_r)) )
            !prefactor = osq4pi*ThExpNb*alpha0(nR)*ViscHeatFac*ogrun(nR)
            !-- Rem: ogrun normalized at the outer radius in radial
            ! but it is then rescaled in preCalculation
            prefactor = ThExpNb*ViscHeatFac*alpha0(nR)*ogrun(nR)
            !prefactor = osq4pi*ViscHeatFac*alpha0(nR)*beta(nR)/dLtemp0(nR) ! equivalent formula if adiabatic background
         end if

         ! Qc_pressure_20_r(nR) = y20_norm*r2*prefactor* real(p(lm20,nR))
         ! Qc_pressure_21_r(nR) = y21_norm*r2*prefactor* real(p(lm21,nR))
         ! Qc_pressure_22_r(nR) = y22_norm*r2*prefactor* real(p(lm22,nR))
         ! Qs_pressure_21_r(nR) = y21_norm*r2*prefactor*aimag(p(lm21,nR))
         ! Qs_pressure_22_r(nR) = y22_norm*r2*prefactor*aimag(p(lm22,nR))
         ! !-- dphi
         ! dPhiQc_pressure_21_r(nR)  = y21_norm*r2*prefactor* real(  ci*p(lm21,nR))
         ! dPhiQc_pressure_22_r(nR)  = y22_norm*r2*prefactor* real(2*ci*p(lm22,nR))
         ! dPhiQs_pressure_21_r(nR)  = y21_norm*r2*prefactor*aimag(  ci*p(lm21,nR))
         ! dPhiQs_pressure_22_r(nR)  = y22_norm*r2*prefactor*aimag(2*ci*p(lm22,nR))
         ! !-- ddphi
         ! ddPhiQc_pressure_21_r(nR) = y21_norm*r2*prefactor* real(  -p(lm21,nR))
         ! ddPhiQc_pressure_22_r(nR) = y22_norm*r2*prefactor* real(-4*p(lm22,nR))
         ! ddPhiQs_pressure_21_r(nR) = y21_norm*r2*prefactor*aimag(  -p(lm21,nR))
         ! ddPhiQs_pressure_22_r(nR) = y22_norm*r2*prefactor*aimag(-4*p(lm22,nR))

         do lm=max(2,llm),ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if (l == 2) then
               if (m == 0) then
                  Qc_pressure_20_r(nR) = y20_norm*r4*prefactor* real(p(lm,nR))
               else if (m == 1) then
                  Qc_pressure_21_r(nR) = y21_norm*r4*prefactor* real(p(lm,nR))
                  Qs_pressure_21_r(nR) = y21_norm*r4*prefactor*aimag(p(lm,nR))

                  dPhiQc_pressure_21_r(nR)  = y21_norm*r4*prefactor* real(  ci*p(lm,nR))
                  dPhiQs_pressure_21_r(nR)  = y21_norm*r4*prefactor*aimag(  ci*p(lm,nR))

                  ddPhiQc_pressure_21_r(nR) = y21_norm*r4*prefactor* real(  -p(lm,nR))
                  ddPhiQs_pressure_21_r(nR) = y21_norm*r4*prefactor*aimag(  -p(lm,nR))
               else if (m == 2) then
                  Qc_pressure_22_r(nR) = y22_norm*r4*prefactor* real(p(lm,nR))
                  Qs_pressure_22_r(nR) = y22_norm*r4*prefactor*aimag(p(lm,nR))

                  dPhiQc_pressure_22_r(nR)  = y22_norm*r4*prefactor* real(2*ci*p(lm,nR))
                  dPhiQs_pressure_22_r(nR)  = y22_norm*r4*prefactor*aimag(2*ci*p(lm,nR))

                  ddPhiQc_pressure_22_r(nR) = y22_norm*r4*prefactor* real(-4*p(lm,nR))
                  ddPhiQs_pressure_22_r(nR) = y22_norm*r4*prefactor*aimag(-4*p(lm,nR))
               end if
            end if
         end do
      end do

      ! reduce over the ranks
#ifdef WITH_MPI
      call MPI_Reduce(Qc_pressure_20_r, Qc_pressure_20_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qc_pressure_21_r, Qc_pressure_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qc_pressure_22_r, Qc_pressure_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qs_pressure_21_r, Qs_pressure_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Qs_pressure_22_r, Qs_pressure_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_Reduce(dPhiQc_pressure_21_r, dPhiQc_pressure_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dPhiQc_pressure_22_r, dPhiQc_pressure_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dPhiQs_pressure_21_r, dPhiQs_pressure_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dPhiQs_pressure_22_r, dPhiQs_pressure_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_Reduce(ddPhiQc_pressure_21_r, ddPhiQc_pressure_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ddPhiQc_pressure_22_r, ddPhiQc_pressure_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ddPhiQs_pressure_21_r, ddPhiQs_pressure_21_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ddPhiQs_pressure_22_r, ddPhiQs_pressure_22_r_global, n_r_max,MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      Qc_pressure_20_r_global = Qc_pressure_20_r
      Qc_pressure_21_r_global = Qc_pressure_21_r
      Qc_pressure_22_r_global = Qc_pressure_22_r
      Qs_pressure_21_r_global = Qs_pressure_21_r
      Qs_pressure_22_r_global = Qs_pressure_22_r

      dPhiQc_pressure_21_r_global = dPhiQc_pressure_21_r
      dPhiQc_pressure_22_r_global = dPhiQc_pressure_22_r
      dPhiQs_pressure_21_r_global = dPhiQs_pressure_21_r
      dPhiQs_pressure_22_r_global = dPhiQs_pressure_22_r

      ddPhiQc_pressure_21_r_global = ddPhiQc_pressure_21_r
      ddPhiQc_pressure_22_r_global = ddPhiQc_pressure_22_r
      ddPhiQs_pressure_21_r_global = ddPhiQs_pressure_21_r
      ddPhiQs_pressure_22_r_global = ddPhiQs_pressure_22_r
#endif

      if ( rank == 0 ) then
         !-- Radial Integrals:
         Qc_pressure_20 = rInt_R(Qc_pressure_20_r_global,r,rscheme_oc)
         Qc_pressure_21 = rInt_R(Qc_pressure_21_r_global,r,rscheme_oc)
         Qc_pressure_22 = rInt_R(Qc_pressure_22_r_global,r,rscheme_oc)
         Qs_pressure_21 = rInt_R(Qs_pressure_21_r_global,r,rscheme_oc)
         Qs_pressure_22 = rInt_R(Qs_pressure_22_r_global,r,rscheme_oc)

         dPhiQc_pressure_21 = 2*oek * rInt_R(dPhiQc_pressure_21_r_global,r,rscheme_oc)
         dPhiQc_pressure_22 = 2*oek * rInt_R(dPhiQc_pressure_22_r_global,r,rscheme_oc)
         dPhiQs_pressure_21 = 2*oek * rInt_R(dPhiQs_pressure_21_r_global,r,rscheme_oc)
         dPhiQs_pressure_22 = 2*oek * rInt_R(dPhiQs_pressure_22_r_global,r,rscheme_oc)

         ddPhiQc_pressure_21 = oek**2 * rInt_R(ddPhiQc_pressure_21_r_global,r,rscheme_oc)
         ddPhiQc_pressure_22 = oek**2 * rInt_R(ddPhiQc_pressure_22_r_global,r,rscheme_oc)
         ddPhiQs_pressure_21 = oek**2 * rInt_R(ddPhiQs_pressure_21_r_global,r,rscheme_oc)
         ddPhiQs_pressure_22 = oek**2 * rInt_R(ddPhiQs_pressure_22_r_global,r,rscheme_oc)

         !-- Write outputs
         if ( l_save_out ) then
            open(newunit=n_gw_P_file, file=gw_P_file,      &
                 &    status='unknown', position='append', &
                 &    form='unformatted')
         end if
         write(n_gw_P_file)  timeScaled,                   & ! 1
              &      Qc_pressure_20,                       & ! 2
              &      Qc_pressure_21,      Qs_pressure_21,  & ! 3,4
              &      Qc_pressure_22,      Qs_pressure_22,  & ! 5,6
              &  dPhiQc_pressure_21,  dPhiQs_pressure_21,  & ! 7,8
              &  dPhiQc_pressure_22,  dPhiQs_pressure_22,  & ! 9,10
              & ddPhiQc_pressure_21, ddPhiQs_pressure_21,  & ! 11,12
              & ddPhiQc_pressure_22, ddPhiQs_pressure_22     ! 13,14

         if ( l_save_out ) close(n_gw_P_file)
      end if
    end subroutine outGWpressure
!---------------------------------------------------------------------------
   subroutine outMagneticHelicity(timeScaled,magHelLMr)

     !-- Input of variables:
     real(cp), intent(in) :: timeScaled
     real(cp), intent(in) :: magHelLMr(l_max+1,nRstart:nRstop)

     !-- Local variables
     integer  :: n_r
     integer  :: nTheta,nThetaStart,nThetaBlock,nThetaNHS,n
     real(cp) :: magHelR(nRstart:nRstop)
     real(cp) :: magHelR_global(n_r_max)
     real(cp) :: maghel(nfs)
     real(cp) :: mag_Hel

     do n_r=nRstart,nRstop
         magHelR(n_r)=0.0_cp
#ifdef WITH_SHTNS
         call axi_to_spat(magHelLMr(:,n_r), maghel)
#endif

         do n=1,nThetaBs ! Loop over theta blocks
            nTheta=(n-1)*sizeThetaB
            nThetaStart=nTheta+1
#ifndef WITH_SHTNS
            call lmAS2pt(magHelLMr(:,n_r),maghel,nThetaStart,sizeThetaB)
#endif
            do nThetaBlock=1,sizeThetaB
               nTheta=nTheta+1
               nThetaNHS=(nTheta+1)/2

               !-- Integration over theta
               !-- r2 factor taken into account in mhel calculation
               magHelR(n_r)=magHelR(n_r)+gauss(nThetaNHS)*maghel(nThetaBlock)
            end do
         end do
     end do

     !-- Gather on rank 0
     call gather_from_Rloc(magHelR, magHelR_global, 0)

      if ( rank == 0 ) then
         !------ Integration over r:
         mag_Hel = rInt_R(magHelR_global,r,rscheme_oc)

         if ( l_save_out ) then
            open(newunit=n_magHel_file, file=magHel_file,   &
                 &    status='unknown', position='append')
         end if

         write(n_magHel_file,'(1P,ES20.12,1ES16.8)')   &
              &     timeScaled, mag_Hel

         if ( l_save_out ) close(n_magHel_file)
      end if

   end subroutine outMagneticHelicity
!---------------------------------------------------------------------------
   subroutine outHeat(time,timePassed,timeNorm,l_stop_time,s,ds,p,dp,xi,dxi)
      !
      ! This subroutine is used to store informations about heat transfer
      ! (i.e. Nusselt number, temperature, entropy, ...)
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time
      real(cp),    intent(in) :: timePassed
      real(cp),    intent(in) :: timeNorm
      logical,     intent(in) :: l_stop_time

      !-- Input of scalar fields:
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dp(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dxi(llm:ulm,n_r_max)

      !-- Local stuff:
      real(cp) :: tmp(n_r_max)
      real(cp) :: topnuss,botnuss,deltanuss
      real(cp) :: topsherwood,botsherwood,deltasherwood
      real(cp) :: toptemp,bottemp
      real(cp) :: topxi,botxi
      real(cp) :: toppres,botpres,mass
      real(cp) :: topentropy, botentropy
      real(cp) :: topflux,botflux
      character(len=76) :: filename
      integer :: n_r, filehandle

      if ( rank == 0 ) then
         n_calls = n_calls + 1
         if ( l_anelastic_liquid ) then
            if ( l_heat ) then
               call TMeanR%compute(osq4pi*real(s(1,:)),n_calls,timePassed,timeNorm)
               tmp(:)   = otemp1(:)*TMeanR%mean(:)-ViscHeatFac*ThExpNb* &
               &          alpha0(n_r)*orho1(n_r)*PmeanR%mean(n_r)
               call SMeanR%compute(tmp(:),n_calls,timePassed,timeNorm)
            endif
            if ( l_chemical_conv ) then
               call XiMeanR%compute(osq4pi*real(xi(1,:)),n_calls,timePassed,timeNorm)
            endif
            call PMeanR%compute(osq4pi*real(p(1,:)),n_calls,timePassed,timeNorm)
            tmp(:) = osq4pi*ThExpNb*alpha0(:)*( -rho0(:)*    &
               &     real(s(1,:))+ViscHeatFac*(ThExpNb*      &
               &     alpha0(:)*temp0(:)+ogrun(:))*           &
               &     real(p(1,:)) )
            call RhoMeanR%compute(tmp(:),n_calls,timePassed,timeNorm)
         else
            if ( l_heat ) then
               call SMeanR%compute(osq4pi*real(s(1,:)),n_calls,timePassed,timeNorm)
               tmp(:) = temp0(:)*SMeanR%mean(:)+ViscHeatFac*ThExpNb* &
               &        alpha0(:)*temp0(:)*orho1(:)*PMeanR%mean(:)
               call TMeanR%compute(tmp(:), n_calls, timePassed, timeNorm)
            endif
            if ( l_chemical_conv ) then
               call XiMeanR%compute(osq4pi*real(xi(1,:)),n_calls,timePassed,timeNorm)
            endif
            call PMeanR%compute(osq4pi*real(p(1,:)),n_calls,timePassed,timeNorm)
            tmp(:) = osq4pi*ThExpNb*alpha0(:)*( -rho0(:)* &
               &     temp0(:)*real(s(1,:))+ViscHeatFac*   &
               &     ogrun(:)*real(p(1,:)) )
            call RhoMeanR%compute(tmp(:),n_calls,timePassed,timeNorm)
         end if

         !-- Evaluate nusselt numbers (boundary heat flux density):
         toppres=osq4pi*real(p(1,n_r_cmb))
         botpres=osq4pi*real(p(1,n_r_icb))
         if ( topcond /= 0.0_cp ) then

            if ( l_anelastic_liquid ) then

               bottemp=osq4pi*real(s(1,n_r_icb))
               toptemp=osq4pi*real(s(1,n_r_cmb))

               botentropy=otemp1(n_r_icb)*bottemp-ViscHeatFac*ThExpNb*   &
               &          orho1(n_r_icb)*alpha0(n_r_icb)*botpres
               topentropy=otemp1(n_r_cmb)*toptemp-ViscHeatFac*ThExpNb*   &
               &          orho1(n_r_cmb)*alpha0(n_r_cmb)*toppres

               if ( l_temperature_diff ) then

                  botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
                  topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
                  botflux=-rho0(n_r_max)*real(ds(1,n_r_max))*osq4pi &
                  &        *r_icb**2*four*pi*kappa(n_r_max)
                  topflux=-rho0(1)*real(ds(1,1))*osq4pi &
                  &        *r_cmb**2*four*pi*kappa(1)

                  deltanuss = deltacond/(bottemp-toptemp)

               else

                  botnuss=-osq4pi/botcond*(otemp1(n_r_icb)*( -dLtemp0(n_r_icb)* &
                  &        real(s(1,n_r_icb)) + real(ds(1,n_r_icb))) -          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_icb)*orho1(n_r_icb)*( &
                  &         ( dLalpha0(n_r_icb)-beta(n_r_icb) )*                &
                  &        real(p(1,n_r_icb)) + real(dp(1,n_r_icb)) ) ) / lScale
                  topnuss=-osq4pi/topcond*(otemp1(n_r_cmb)*( -dLtemp0(n_r_cmb)* &
                  &        real(s(1,n_r_cmb)) + real(ds(1,n_r_cmb))) -          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_cmb)*orho1(n_r_cmb)*( &
                  &         ( dLalpha0(n_r_cmb)-beta(n_r_cmb) )*                &
                  &        real(p(1,n_r_cmb)) + real(dp(1,n_r_cmb)) ) ) / lScale

                  botflux=four*pi*r_icb**2*kappa(n_r_icb)*rho0(n_r_icb) *      &
                  &       botnuss*botcond*lScale*temp0(n_r_icb)
                  topflux=four*pi*r_cmb**2*kappa(n_r_cmb)*rho0(n_r_cmb) *      &
                  &       topnuss*topcond*lScale*temp0(n_r_cmb)

                  deltanuss = deltacond/(botentropy-topentropy)

               end if

            else ! s corresponds to entropy

               botentropy=osq4pi*real(s(1,n_r_icb))
               topentropy=osq4pi*real(s(1,n_r_cmb))

               bottemp   =temp0(n_r_icb)*botentropy+ViscHeatFac*ThExpNb*   &
               &          orho1(n_r_icb)*temp0(n_r_icb)*alpha0(n_r_icb)*   &
               &          botpres
               toptemp   =temp0(n_r_cmb)*topentropy+ViscHeatFac*ThExpNb*   &
               &          orho1(n_r_cmb)*temp0(n_r_cmb)*alpha0(n_r_cmb)*   &
               &          toppres

               if ( l_temperature_diff ) then

                  botnuss=-osq4pi/botcond*temp0(n_r_icb)*( dLtemp0(n_r_icb)*   &
                  &        real(s(1,n_r_icb)) + real(ds(1,n_r_icb)) +          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_icb)*orho1(n_r_icb)*(&
                  &     ( dLalpha0(n_r_icb)+dLtemp0(n_r_icb)-beta(n_r_icb) )*  &
                  &        real(p(1,n_r_icb)) + real(dp(1,n_r_icb)) ) ) / lScale
                  topnuss=-osq4pi/topcond*temp0(n_r_cmb)*( dLtemp0(n_r_cmb)*   &
                  &        real(s(1,n_r_cmb)) + real(ds(1,n_r_cmb)) +          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_cmb)*orho1(n_r_cmb)*(&
                  &     ( dLalpha0(n_r_cmb)+dLtemp0(n_r_cmb)-beta(n_r_cmb) )*  &
                  &        real(p(1,n_r_cmb)) + real(dp(1,n_r_cmb)) ) ) / lScale

                  botflux=four*pi*r_icb**2*kappa(n_r_icb)*rho0(n_r_icb) *      &
                  &       botnuss*botcond*lScale
                  topflux=four*pi*r_cmb**2*kappa(n_r_cmb)*rho0(n_r_cmb) *      &
                  &       topnuss*topcond*lScale

                  deltanuss = deltacond/(bottemp-toptemp)

               else

                  botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
                  topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
                  botflux=-rho0(n_r_max)*temp0(n_r_max)*real(ds(1,n_r_max))* &
                  &        r_icb**2*sq4pi*kappa(n_r_max)/lScale
                  topflux=-rho0(1)*temp0(1)*real(ds(1,1))/lScale*r_cmb**2* &
                  &        sq4pi*kappa(1)
                  if ( botentropy /= topentropy ) then
                     deltanuss = deltacond/(botentropy-topentropy)
                  else
                     deltanuss = one
                  end if

               end if

            end if
         else
            botnuss   =one
            topnuss   =one
            botflux   =0.0_cp
            topflux   =0.0_cp
            bottemp   =0.0_cp
            toptemp   =0.0_cp
            botentropy=0.0_cp
            topentropy=0.0_cp
            deltanuss =one
         end if

         if ( l_chemical_conv ) then
            if ( topxicond/=0.0_cp ) then
               topxi=osq4pi*real(xi(1,n_r_cmb))
               botxi=osq4pi*real(xi(1,n_r_icb))
               botsherwood=-osq4pi/botxicond*real(dxi(1,n_r_icb))/lScale
               topsherwood=-osq4pi/topxicond*real(dxi(1,n_r_cmb))/lScale
               deltasherwood = deltaxicond/(botxi-topxi)
            else
               topxi=0.0_cp
               botxi=0.0_cp
               botsherwood=one
               topsherwood=one
               deltasherwood=one
            end if
         else
            topxi=0.0_cp
            botxi=0.0_cp
            botsherwood=one
            topsherwood=one
            deltasherwood=one
         end if

         tmp(:)=tmp(:)*r(:)*r(:)
         mass=four*pi*rInt_R(tmp,r,rscheme_oc)

         if ( l_save_out ) then
            open(newunit=n_heat_file, file=heat_file, status='unknown', &
            &    position='append')
         end if

         !-- avoid too small number in output
         if ( abs(toppres) <= 1e-11_cp ) toppres=0.0_cp

         if ( abs(mass) <= 1e-11_cp ) mass=0.0_cp

         write(n_heat_file,'(1P,ES20.12,16ES16.8)')          &
         &     time, botnuss, topnuss, deltanuss,            &
         &     bottemp, toptemp, botentropy, topentropy,     &
         &     botflux, topflux, toppres, mass, topsherwood, &
         &     botsherwood, deltasherwood, botxi, topxi

         if ( l_save_out ) close(n_heat_file)

         if ( l_stop_time ) then
            call SMeanR%finalize_SD(timeNorm)
            call TMeanR%finalize_SD(timeNorm)
            call PMeanR%finalize_SD(timeNorm)
            call XiMeanR%finalize_SD(timeNorm)
            call RhoMeanR%finalize_SD(timeNorm)

            filename='heatR.'//tag
            open(newunit=filehandle, file=filename, status='unknown')
            do n_r=1,n_r_max
               write(filehandle, '(ES20.10,5ES15.7,5ES13.5)' )                &
               &      r(n_r),round_off(SMeanR%mean(n_r),maxval(SMeanR%mean)), &
               &      round_off(TMeanR%mean(n_r),maxval(TMeanR%mean)),        &
               &      round_off(PMeanR%mean(n_r),maxval(PMeanR%mean)),        &
               &      round_off(RhoMeanR%mean(n_r),maxval(RhoMeanR%mean)),    &
               &      round_off(XiMeanR%mean(n_r),maxval(XiMeanR%mean)),      &
               &      round_off(SMeanR%SD(n_r),maxval(SMeanR%SD)),            &
               &      round_off(TMeanR%SD(n_r),maxval(TMeanR%SD)),            &
               &      round_off(PMeanR%SD(n_r),maxval(PMeanR%SD)),            &
               &      round_off(RhoMeanR%SD(n_r),maxval(RhoMeanR%SD)),        &
               &      round_off(XiMeanR%SD(n_r),maxval(XiMeanR%SD))
            end do

            close(filehandle)
         end if

      end if ! rank == 0

   end subroutine outHeat
!---------------------------------------------------------------------------
end module outMisc_mod
