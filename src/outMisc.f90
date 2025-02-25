module outMisc_mod
   !
   ! This module contains several subroutines that can compute and store
   ! various informations: helicity (helicity.TAG), heat transfer (heat.TAG),
   ! phase field (phase.TAG) and North/South hemisphericity of energies (hemi.TAG)
   !

   use parallel_mod
   use precision_mod

!
   use mem_alloc, only: bytes_allocated
   use communications, only: gather_from_Rloc, gather_from_lo_to_rank0
   use truncation, only: l_max, n_r_max, nlat_padded, n_theta_max, n_r_maxMag, &
       &                 n_phi_max, lm_max, m_min, m_max, minc
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop, nRstartMag, nRstopMag
   use radial_functions, only: r_icb, rscheme_oc, kappa, r_cmb,temp0, r, rho0, &
       &                       dLtemp0, dLalpha0, beta, orho1, alpha0, otemp1, &
       &                       ogrun, rscheme_oc, or2, orho2, or4
   use physical_parameters, only: ViscHeatFac, ThExpNb, opr, stef, LFfac, oek
   use num_param, only: lScale, eScale, vScale
   use blocking, only: llm, ulm, lo_map, lm2
   use radial_der, only: get_dr
   use mean_sd, only: mean_sd_type
   use horizontal_data, only: gauss, theta_ord, n_theta_cal2ord,  &
       &                      O_sin_theta_E2
   use logic, only: l_save_out, l_anelastic_liquid, l_heat, l_hel, l_hemi, &
       &            l_temperature_diff, l_chemical_conv, l_phase_field,    &
       &            l_mag, l_onset, l_gw
   use output_data, only: tag
   use constants, only: pi, vol_oc, osq4pi, sq4pi, one, two, four, half, zero, ci
   use start_fields, only: topcond, botcond, deltacond, topxicond, botxicond, &
       &                   deltaxicond
   use useful, only: cc2real, round_off, abortRun
   use integration, only: rInt_R
   use sht, only: axi_to_spat

   implicit none

   private

   type(mean_sd_type) :: TMeanR, SMeanR, PMeanR, XiMeanR, RhoMeanR, PhiMeanR
   integer :: n_heat_file, n_helicity_file, n_calls, n_phase_file, n_gw_S_file, n_gw_P_file
   integer :: n_rmelt_file, n_hemi_file, n_growth_sym_file, n_growth_asym_file
   integer :: n_drift_sym_file, n_drift_asym_file
   character(len=72) :: heat_file, helicity_file, phase_file, rmelt_file, gw_S_file, gw_P_file
   character(len=72) :: hemi_file, sym_file, asym_file, drift_sym_file
   character(len=72) :: drift_asym_file
   real(cp) :: TPhiOld, Tphi
   real(cp), allocatable :: ekinSr(:), ekinLr(:), volSr(:)
   real(cp), allocatable :: hemi_ekin_r(:,:), hemi_vrabs_r(:,:)
   real(cp), allocatable :: hemi_emag_r(:,:), hemi_brabs_r(:,:)
   real(cp), allocatable :: HelASr(:,:), Hel2ASr(:,:)
   real(cp), allocatable :: HelnaASr(:,:), Helna2ASr(:,:)
   real(cp), allocatable :: HelEAASr(:)
   complex(cp), allocatable :: coeff_old(:)

   public :: outHelicity, outHeat, initialize_outMisc_mod, finalize_outMisc_mod, &
   &         outPhase, outHemi, get_ekin_solid_liquid, get_helicity, get_hemi,   &
   &         get_onset, outGWentropy, outGWpressure

contains

   subroutine initialize_outMisc_mod()
      !
      ! This subroutine handles the opening of some output diagnostic files that
      ! have to do with heat transfer, helicity, phase field, hemisphericity or GWs
      !

      if (l_heat .or. l_chemical_conv) then
         call TMeanR%initialize(1,n_r_max)
         call SMeanR%initialize(1,n_r_max)
         call PMeanR%initialize(1,n_r_max)
         call XiMeanR%initialize(1,n_r_max)
         call RhoMeanR%initialize(1,n_r_max)
      endif
      if ( l_phase_field ) call PhiMeanR%initialize(1,n_r_max)
      n_calls = 0

      TPhiOld = 0.0_cp
      TPhi = 0.0_cp

      if ( l_hel ) then
         helicity_file='helicity.'//tag
         allocate( HelASr(nRstart:nRstop,2), Hel2ASr(nRstart:nRstop,2) )
         allocate( HelnaASr(nRstart:nRstop,2), Helna2ASr(nRstart:nRstop,2) )
         allocate( HelEAASr(nRstart:nRstop) )
         HelASr(:,:)   =0.0_cp
         Hel2ASr(:,:)  =0.0_cp
         HelnaASr(:,:) =0.0_cp
         Helna2ASr(:,:)=0.0_cp
         HelEAASr(:)   =0.0_cp
         bytes_allocated=bytes_allocated+9*(nRstop-nRstart+1)*SIZEOF_DEF_REAL
      end if

      if ( l_hemi ) then
         hemi_file    ='hemi.'//tag
         allocate( hemi_ekin_r(nRstart:nRstop,2), hemi_vrabs_r(nRstart:nRstop,2) )
         hemi_ekin_r(:,:) =0.0_cp
         hemi_vrabs_r(:,:)=0.0_cp
         bytes_allocated=bytes_allocated+(nRstop-nRstart+1)*4*SIZEOF_DEF_REAL

         if ( l_mag ) then
            allocate( hemi_emag_r(nRstart:nRstop,2), hemi_brabs_r(nRstart:nRstop,2) )
            hemi_emag_r(:,:) =0.0_cp
            hemi_brabs_r(:,:)=0.0_cp
            bytes_allocated=bytes_allocated+(nRstop-nRstart+1)*4*SIZEOF_DEF_REAL
         end if
      end if

      heat_file    ='heat.'//tag
      gw_P_file    ='gwPressure.'//tag
      gw_S_file    ='gwEntropy.'//tag
      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_hel ) open(newunit=n_helicity_file, file=helicity_file, status='new')
         if ( l_hemi ) open(newunit=n_hemi_file, file=hemi_file, status='new')
         if ( l_heat .or. l_chemical_conv ) then
            open(newunit=n_heat_file, file=heat_file, status='new')
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

      if ( l_phase_field ) then
         phase_file='phase.'//tag
         rmelt_file='rmelt.'//tag
         allocate( ekinSr(nRstart:nRstop), ekinLr(nRstart:nRstop), volSr(nRstart:nRstop) )
         ekinSr(:)=0.0_cp
         ekinLr(:)=0.0_cp
         volSr(:) =0.0_cp
         bytes_allocated=bytes_allocated+3*(nRstop-nRstart+1)*SIZEOF_DEF_REAL
         if ( rank == 0 .and. (.not. l_save_out) ) then
            open(newunit=n_phase_file, file=phase_file, status='new')
            open(newunit=n_rmelt_file, file=rmelt_file, status='new', form='unformatted')
         end if
      end if

      if ( l_onset ) then
         allocate(coeff_old(llm:ulm) )
         bytes_allocated=bytes_allocated+(ulm-llm+1)*SIZEOF_DEF_COMPLEX
         coeff_old(:)=zero
         sym_file       ='growth_sym.'//tag
         asym_file      ='growth_asym.'//tag
         drift_sym_file ='drift_sym.'//tag
         drift_asym_file='drift_asym.'//tag
         if ( rank == 0 .and. (.not. l_save_out)) then
            open(newunit=n_growth_sym_file, file=sym_file, status='new')
            open(newunit=n_growth_asym_file, file=asym_file, status='new')
            open(newunit=n_drift_sym_file, file=drift_sym_file, status='new')
            open(newunit=n_drift_asym_file, file=drift_asym_file, status='new')
         end if
      end if

   end subroutine initialize_outMisc_mod
!----------------------------------------------------------------------------------
   subroutine finalize_outMisc_mod()
      !
      ! This subroutine handles the closing of the time series of
      ! heat.TAG, hel.TAG, hemi.TAG, phase.TAG and GW*.TAG
      !

      if ( l_heat .or. l_chemical_conv ) then
         call TMeanR%finalize()
         call SMeanR%finalize()
         call PMeanR%finalize()
         call XiMeanR%finalize()
         call RhoMeanR%finalize()
      end if
      if ( l_onset ) then
         deallocate( coeff_old )
         if ( rank == 0 .and. (.not. l_save_out) ) then
            close(n_growth_sym_file)
            close(n_growth_asym_file)
            close(n_drift_sym_file)
            close(n_drift_asym_file)
         end if
      end if
      if ( l_hemi ) then 
         deallocate( hemi_ekin_r, hemi_vrabs_r )
         if ( l_mag ) deallocate( hemi_emag_r, hemi_brabs_r )
      end if
      if ( l_hel ) deallocate( HelASr, Hel2ASr, HelnaASr, Helna2ASr, HelEAASr )
      
      if ( l_phase_field ) then
         deallocate( ekinSr, ekinLr, volSr )
         call PhiMeanR%finalize()
      end if

      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_hel ) close(n_helicity_file)
         if ( l_hemi ) close(n_hemi_file)
         if ( l_heat .or. l_chemical_conv ) close(n_heat_file)
         if (l_gw) then
            close(n_gw_P_file)
            if ( l_heat .or. l_chemical_conv ) close(n_gw_S_file)
         end if
         if ( l_phase_field ) then
            close(n_phase_file)
            close(n_rmelt_file)
         end if
      end if
   end subroutine finalize_outMisc_mod
!----------------------------------------------------------------------------------
   subroutine outHemi(timeScaled)
      !
      ! This function handles the writing of outputs related to hemisphericity of 
      ! the kinetic and magnetic energy between Northern and Southern hemispheres.
      ! This is based on Wieland Dietrich's work (see Dietrich & Wicht, 2013).
      ! Outputs are stored in the time series hemi.TAG
      !

      !-- Input of variables:
      real(cp), intent(in) :: timeScaled

      !-- Local variables:
      real(cp) :: hemi_ekin_r_N(n_r_max), hemi_ekin_r_S(n_r_max)
      real(cp) :: hemi_vrabs_r_N(n_r_max), hemi_vrabs_r_S(n_r_max)
      real(cp) :: hemi_emag_r_N(n_r_maxMag), hemi_emag_r_S(n_r_maxMag)
      real(cp) :: hemi_brabs_r_N(n_r_maxMag), hemi_brabs_r_S(n_r_maxMag)
      real(cp) :: ekinN, ekinS, vrabsN, vrabsS, hemi_ekin, hemi_vr
      real(cp) :: emagN, emagS, brabsN, brabsS, hemi_emag, hemi_br, hemi_cmb

      call gather_from_Rloc(hemi_ekin_r(:,1), hemi_ekin_r_N, 0)
      call gather_from_Rloc(hemi_ekin_r(:,2), hemi_ekin_r_S, 0)
      call gather_from_Rloc(hemi_vrabs_r(:,1), hemi_vrabs_r_N, 0)
      call gather_from_Rloc(hemi_vrabs_r(:,2), hemi_vrabs_r_S, 0)
      if ( l_mag ) then
         call gather_from_Rloc(hemi_emag_r(:,1), hemi_emag_r_N, 0)
         call gather_from_Rloc(hemi_emag_r(:,2), hemi_emag_r_S, 0)
         call gather_from_Rloc(hemi_brabs_r(:,1), hemi_brabs_r_N, 0)
         call gather_from_Rloc(hemi_brabs_r(:,2), hemi_brabs_r_S, 0)
      end if

      if ( rank == 0 ) then
         !------ Integration over r
         ekinN =eScale*rInt_R(hemi_ekin_r_N,r,rscheme_oc)
         ekinS =eScale*rInt_R(hemi_ekin_r_S,r,rscheme_oc)
         vrabsN=vScale*rInt_R(hemi_vrabs_r_N,r,rscheme_oc)
         vrabsS=vScale*rInt_R(hemi_vrabs_r_S,r,rscheme_oc)
         if ( l_mag ) then
            emagN =LFfac*eScale*rInt_R(hemi_emag_r_N,r,rscheme_oc)
            emagS =LFfac*eScale*rInt_R(hemi_emag_r_S,r,rscheme_oc)
            brabsN=rInt_R(hemi_brabs_r_N,r,rscheme_oc)
            brabsS=rInt_R(hemi_brabs_r_S,r,rscheme_oc)
            if ( emagN+emagS > 0.0_cp ) then
               hemi_emag=abs(emagN-emagS)/(emagN+emagS)
               hemi_br  =abs(brabsN-brabsS)/(brabsN+brabsS)
               hemi_cmb =abs(hemi_brabs_r_N(n_r_cmb)-hemi_brabs_r_S(n_r_cmb)) / &
               &         (hemi_brabs_r_N(n_r_cmb)+hemi_brabs_r_S(n_r_cmb))
            else
               hemi_emag=0.0_cp
               hemi_br  =0.0_cp
               hemi_cmb =0.0_cp
            end if
         else
            emagN    =0.0_cp
            emagS    =0.0_cp
            hemi_emag=0.0_cp
            hemi_br  =0.0_cp
            hemi_cmb =0.0_cp
         end if

         if ( ekinN+ekinS > 0.0_cp ) then
            hemi_ekin=abs(ekinN-ekinS)/(ekinN+ekinS)
            hemi_vr  =abs(vrabsN-vrabsS)/(vrabsN+vrabsS)
         else
            hemi_ekin=0.0_cp
            hemi_vr  =0.0_cp
         end if

         if ( l_save_out ) then
            open(newunit=n_hemi_file, file=hemi_file,   &
            &    status='unknown', position='append')
         end if

         write(n_hemi_file,'(1P,ES20.12,7ES16.8)')                &
         &     timeScaled, round_off(hemi_vr,one),                &
         &     round_off(hemi_ekin,one), round_off(hemi_br,one),  &
         &     round_off(hemi_emag,one), round_off(hemi_cmb,one), &
         &     ekinN+ekinS, emagN+emagS

         if ( l_save_out ) close(n_hemi_file)

      end if

   end subroutine outHemi
!----------------------------------------------------------------------------------
   subroutine outHelicity(timeScaled)
      !
      ! This subroutine is used to store informations about kinetic
      ! helicity. Outputs are stored in the time series helicity.TAG
      !

      !-- Input of variables:
      real(cp), intent(in) :: timeScaled

      !-- Local stuff:
      real(cp) :: HelNr_global(n_r_max), HelSr_global(n_r_max)
      real(cp) :: HelnaNr_global(n_r_max), HelnaSr_global(n_r_max)
      real(cp) :: Helna2Nr_global(n_r_max), Helna2Sr_global(n_r_max)
      real(cp) :: Hel2Nr_global(n_r_max), Hel2Sr_global(n_r_max)
      real(cp) :: HelEAr_global(n_r_max)
      real(cp) :: HelN,HelS,HelnaN,HelnaS,HelnaRMSN,HelnaRMSS
      real(cp) :: HelRMSN,HelRMSS,HelEA,HelRMS,HelnaRMS

      ! Now we have to gather the results on rank 0 for
      ! the arrays: Hel2Nr,Helna2Nr,HelEAr,HelNr,HelnaNr
      ! Hel2Sr,Helna2Sr,HelSr,HelnaSr

      call gather_from_Rloc(Hel2Asr(:,1), Hel2Nr_global, 0)
      call gather_from_Rloc(Helna2ASr(:,1), Helna2Nr_global, 0)
      call gather_from_Rloc(HelEAASr, HelEAr_global, 0)
      call gather_from_Rloc(HelASr(:,1), HelNr_global, 0)
      call gather_from_Rloc(HelnaASr(:,1), HelnaNr_global, 0)
      call gather_from_Rloc(HelASr(:,2), HelSr_global, 0)
      call gather_from_Rloc(Helna2ASr(:,2), Helna2Sr_global, 0)
      call gather_from_Rloc(Hel2ASr(:,2), Hel2Sr_global, 0)
      call gather_from_Rloc(HelnaASr(:,2), HelnaSr_global, 0)

      if ( rank == 0 ) then
         HelN  =rInt_R(HelNr_global*r*r,r,rscheme_oc)
         HelS  =rInt_R(HelSr_global*r*r,r,rscheme_oc)
         HelnaN=rInt_R(HelnaNr_global*r*r,r,rscheme_oc)
         HelnaS=rInt_R(HelnaSr_global*r*r,r,rscheme_oc)
         HelEA =rInt_R(HelEAr_global*r*r,r,rscheme_oc)
         HelRMSN=rInt_R(Hel2Nr_global*r*r,r,rscheme_oc)
         HelRMSS=rInt_R(Hel2Sr_global*r*r,r,rscheme_oc)
         HelnaRMSN=rInt_R(Helna2Nr_global*r*r,r,rscheme_oc)
         HelnaRMSS=rInt_R(Helna2Sr_global*r*r,r,rscheme_oc)

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

         if ( l_anelastic_liquid ) then
            ! rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &                 
            !      &               real(s(1,n_r))+ViscHeatFac*(ThExpNb*     &            
            !      &               alpha0(n_r)*temp0(n_r)+ogrun(n_r))*      &            
            !      &               real(p(1,n_r)) )             
            call abortRun('This setup is not ready')
         else
            prefactor = -alpha0(nR)*rho0(nR)*temp0(nR)/ViscHeatFac
         end if

         do lm=max(2,llm), ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if (l == 2) then
               if (m == 0) then
                  Qc_entropy_20_r(nR) = r4*prefactor* real(s(lm,nR))
               else if (m == 1) then
                  Qc_entropy_21_r(nR) = r4*prefactor* real(s(lm,nR))
                  Qs_entropy_21_r(nR) = r4*prefactor*aimag(s(lm,nR))

                  dPhiQc_entropy_21_r(nR)  = r4*prefactor* real(  ci*s(lm,nR))
                  dPhiQs_entropy_21_r(nR)  = r4*prefactor*aimag(  ci*s(lm,nR))

                  ddPhiQc_entropy_21_r(nR) = r4*prefactor* real(  -s(lm,nR))
                  ddPhiQs_entropy_21_r(nR) = r4*prefactor*aimag(  -s(lm,nR))
               else if (m == 2) then
                  Qc_entropy_22_r(nR) = r4*prefactor* real(s(lm,nR))
                  Qs_entropy_22_r(nR) = r4*prefactor*aimag(s(lm,nR))

                  dPhiQc_entropy_22_r(nR)  = r4*prefactor* real(2*ci*s(lm,nR))
                  dPhiQs_entropy_22_r(nR)  = r4*prefactor*aimag(2*ci*s(lm,nR))

                  ddPhiQc_entropy_22_r(nR) = r4*prefactor* real(-4*s(lm,nR))
                  ddPhiQs_entropy_22_r(nR) = r4*prefactor*aimag(-4*s(lm,nR))
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

         if ( l_anelastic_liquid ) then
            call abortRun('This setup is not ready')
         else
            !-- Rem: ogrun normalized at the outer radius in radial                      
            ! but it is then rescaled in preCalculation         
            prefactor = alpha0(nR)*ogrun(nR)
         end if

         do lm=max(2,llm),ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if (l == 2) then
               if (m == 0) then
                  Qc_pressure_20_r(nR) = r4*prefactor* real(p(lm,nR))
               else if (m == 1) then
                  Qc_pressure_21_r(nR) = r4*prefactor* real(p(lm,nR))
                  Qs_pressure_21_r(nR) = r4*prefactor*aimag(p(lm,nR))

                  dPhiQc_pressure_21_r(nR)  = r4*prefactor* real(  ci*p(lm,nR))
                  dPhiQs_pressure_21_r(nR)  = r4*prefactor*aimag(  ci*p(lm,nR))

                  ddPhiQc_pressure_21_r(nR) = r4*prefactor* real(  -p(lm,nR))
                  ddPhiQs_pressure_21_r(nR) = r4*prefactor*aimag(  -p(lm,nR))
               else if (m == 2) then
                  Qc_pressure_22_r(nR) = r4*prefactor* real(p(lm,nR))
                  Qs_pressure_22_r(nR) = r4*prefactor*aimag(p(lm,nR))

                  dPhiQc_pressure_22_r(nR)  = r4*prefactor* real(2*ci*p(lm,nR))
                  dPhiQs_pressure_22_r(nR)  = r4*prefactor*aimag(2*ci*p(lm,nR))

                  ddPhiQc_pressure_22_r(nR) = r4*prefactor* real(-4*p(lm,nR))
                  ddPhiQs_pressure_22_r(nR) = r4*prefactor*aimag(-4*p(lm,nR))
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
   subroutine outHeat(time,timePassed,timeNorm,l_stop_time,s,ds,p,xi,dxi)
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
      complex(cp), intent(in) :: s(llm:ulm,n_r_max) ! Entropy/temperature
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max) ! Radial derivative of entropy/temp
      complex(cp), intent(in) :: p(llm:ulm,n_r_max) ! Pressure
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max) ! Chemical composition
      complex(cp), intent(in) :: dxi(llm:ulm,n_r_max) ! Radial derivative of xi

      !-- Local stuff:
      real(cp) :: tmp(n_r_max), dp(n_r_max)
      real(cp) :: topnuss,botnuss,deltanuss
      real(cp) :: topsherwood,botsherwood,deltasherwood
      real(cp) :: toptemp,bottemp,topxi,botxi,topflux,botflux
      real(cp) :: toppres,botpres,mass,topentropy,botentropy
      character(len=76) :: filename
      integer :: n_r, filehandle

      if ( rank == 0 ) then
         tmp(:) = real(p(1,:))
         call get_dr(tmp, dp, n_r_max, rscheme_oc)
         n_calls = n_calls + 1
         if ( l_anelastic_liquid ) then
            if ( l_heat ) then
               call TMeanR%compute(osq4pi*real(s(1,:)),n_calls,timePassed,timeNorm)
               tmp(:)   = otemp1(:)*TMeanR%mean(:)-ViscHeatFac*ThExpNb* &
               &          alpha0(:)*orho1(:)*PmeanR%mean(:)
               call SMeanR%compute(tmp(:),n_calls,timePassed,timeNorm)
            endif
            if ( l_chemical_conv ) then
               call XiMeanR%compute(osq4pi*real(xi(1,:)),n_calls,timePassed,timeNorm)
            endif
            call PMeanR%compute(osq4pi*real(p(1,:)),n_calls,timePassed,timeNorm)
            tmp(:) = osq4pi*ThExpNb*alpha0(:)*( -rho0(:)*real(s(1,:)) +  &
            &        ViscHeatFac*(ThExpNb*alpha0(:)*temp0(:)+ogrun(:))*  &
            &        real(p(1,:)) )
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
            tmp(:) = osq4pi*ThExpNb*alpha0(:)*( -rho0(:)*temp0(:)*real(s(1,:)) + &
            &        ViscHeatFac*ogrun(:)*real(p(1,:)) )
            call RhoMeanR%compute(tmp(:),n_calls,timePassed,timeNorm)
         end if

         !-- Evaluate nusselt numbers (boundary heat flux density):
         toppres=osq4pi*real(p(1,n_r_cmb))
         botpres=osq4pi*real(p(1,n_r_icb))
         if ( l_anelastic_liquid ) then

            bottemp=osq4pi*real(s(1,n_r_icb))
            toptemp=osq4pi*real(s(1,n_r_cmb))

            botentropy=otemp1(n_r_icb)*bottemp-ViscHeatFac*ThExpNb*   &
            &          orho1(n_r_icb)*alpha0(n_r_icb)*botpres
            topentropy=otemp1(n_r_cmb)*toptemp-ViscHeatFac*ThExpNb*   &
            &          orho1(n_r_cmb)*alpha0(n_r_cmb)*toppres

            if ( l_temperature_diff ) then

               if ( abs(botcond) >= 1e-10_cp ) then
                  botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
               else
                  botnuss=one
               end if
               if ( abs(topcond) >= 1e-10_cp ) then
                  topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
               else
                  topnuss=one
               end if
               botflux=-rho0(n_r_max)*real(ds(1,n_r_max))*osq4pi &
               &        *r_icb**2*four*pi*kappa(n_r_max)
               topflux=-rho0(1)*real(ds(1,1))*osq4pi &
               &        *r_cmb**2*four*pi*kappa(1)

               if ( bottemp /= toptemp ) then
                  deltanuss=deltacond/(bottemp-toptemp)
               else
                  deltanuss=one
               end if

            else

               if ( abs(botcond) >= 1e-10_cp ) then
                  botnuss=-osq4pi/botcond*(otemp1(n_r_icb)*( -dLtemp0(n_r_icb)* &
                  &        real(s(1,n_r_icb)) + real(ds(1,n_r_icb))) -          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_icb)*orho1(n_r_icb)*( &
                  &         ( dLalpha0(n_r_icb)-beta(n_r_icb) )*                &
                  &        real(p(1,n_r_icb)) + dp(n_r_icb) ) ) / lScale
               else
                  botnuss=one
               end if
               if ( abs(topcond) >= 1e-10_cp ) then
                  topnuss=-osq4pi/topcond*(otemp1(n_r_cmb)*( -dLtemp0(n_r_cmb)* &
                  &        real(s(1,n_r_cmb)) + real(ds(1,n_r_cmb))) -          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_cmb)*orho1(n_r_cmb)*( &
                  &         ( dLalpha0(n_r_cmb)-beta(n_r_cmb) )*                &
                  &        real(p(1,n_r_cmb)) + dp(n_r_cmb) ) ) / lScale
               else
                  topnuss=one
               end if

               botflux=four*pi*r_icb**2*kappa(n_r_icb)*rho0(n_r_icb) *      &
               &       botnuss*botcond*lScale*temp0(n_r_icb)
               topflux=four*pi*r_cmb**2*kappa(n_r_cmb)*rho0(n_r_cmb) *      &
               &       topnuss*topcond*lScale*temp0(n_r_cmb)

               if ( botentropy /= topentropy ) then
                  deltanuss=deltacond/(botentropy-topentropy)
               else
                  deltanuss=one
               end if

            end if

         else ! s corresponds to entropy

            botentropy=osq4pi*real(s(1,n_r_icb))
            topentropy=osq4pi*real(s(1,n_r_cmb))

            bottemp=temp0(n_r_icb)*botentropy+ViscHeatFac*ThExpNb*   &
            &       orho1(n_r_icb)*temp0(n_r_icb)*alpha0(n_r_icb)*   &
            &       botpres
            toptemp=temp0(n_r_cmb)*topentropy+ViscHeatFac*ThExpNb*   &
            &       orho1(n_r_cmb)*temp0(n_r_cmb)*alpha0(n_r_cmb)*   &
            &       toppres

            if ( l_temperature_diff ) then

               if ( abs(botcond) >= 1e-10_cp ) then
                  botnuss=-osq4pi/botcond*temp0(n_r_icb)*( dLtemp0(n_r_icb)*   &
                  &        real(s(1,n_r_icb)) + real(ds(1,n_r_icb)) +          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_icb)*orho1(n_r_icb)*(&
                  &     ( dLalpha0(n_r_icb)+dLtemp0(n_r_icb)-beta(n_r_icb) )*  &
                  &        real(p(1,n_r_icb)) + dp(n_r_icb) ) ) / lScale
               else
                  botnuss=one
               end if
               if ( abs(topcond) >= 1e-10_cp ) then
                  topnuss=-osq4pi/topcond*temp0(n_r_cmb)*( dLtemp0(n_r_cmb)*   &
                  &        real(s(1,n_r_cmb)) + real(ds(1,n_r_cmb)) +          &
                  &        ViscHeatFac*ThExpNb*alpha0(n_r_cmb)*orho1(n_r_cmb)*(&
                  &     ( dLalpha0(n_r_cmb)+dLtemp0(n_r_cmb)-beta(n_r_cmb) )*  &
                  &        real(p(1,n_r_cmb)) + dp(n_r_cmb) ) ) / lScale
               else
                  topnuss=one
               end if

               botflux=four*pi*r_icb**2*kappa(n_r_icb)*rho0(n_r_icb) *      &
               &       botnuss*botcond*lScale
               topflux=four*pi*r_cmb**2*kappa(n_r_cmb)*rho0(n_r_cmb) *      &
               &       topnuss*topcond*lScale

               if ( bottemp /= toptemp ) then
                  deltanuss=deltacond/(bottemp-toptemp)
               else
                  deltanuss=one
               end if

            else

               if ( abs(botcond) >= 1e-10_cp ) then
                  botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
               else
                  botnuss=one
               end if
               if ( abs(topcond) >= 1e-10_cp ) then
                  topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
               else
                  topnuss=one
               end if
               botflux=-rho0(n_r_max)*temp0(n_r_max)*real(ds(1,n_r_max))* &
               &        r_icb**2*sq4pi*kappa(n_r_max)/lScale
               topflux=-rho0(1)*temp0(1)*real(ds(1,1))/lScale*r_cmb**2* &
               &        sq4pi*kappa(1)
               if ( botentropy /= topentropy ) then
                  deltanuss=deltacond/(botentropy-topentropy)
               else
                  deltanuss=one
               end if

            end if

         end if

         if ( l_chemical_conv ) then
            topxi=osq4pi*real(xi(1,n_r_cmb))
            botxi=osq4pi*real(xi(1,n_r_icb))
            if ( abs(botxicond) >= 1e-10_cp ) then
               botsherwood=-osq4pi/botxicond*real(dxi(1,n_r_icb))/lScale
            else
               botsherwood=one
            end if
            if ( abs(topxicond) >= 1e-10_cp ) then
               topsherwood=-osq4pi/topxicond*real(dxi(1,n_r_cmb))/lScale
            else
               topsherwood=one
            end if
            if ( botxi /= topxi ) then
               deltasherwood=deltaxicond/(botxi-topxi)
            else
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
!----------------------------------------------------------------------------------
   subroutine outPhase(time, timePassed, timeNorm, l_stop_time, nLogs, s, ds, phi)
      !
      ! This subroutine handles the writing of time series related with phase
      ! field: phase.TAG
      !

      !-- Input variables
      real(cp),    intent(in) :: time                   ! Time
      real(cp),    intent(in) :: timePassed             ! Time passed since last call
      real(cp),    intent(in) :: timeNorm
      logical,     intent(in) :: l_stop_time            ! Last iteration
      integer,     intent(in) :: nLogs                  ! Number of log outputs
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)     ! Entropy/Temperature
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)    ! Radial der. of Entropy/Temperature
      complex(cp), intent(in) :: phi(llm:ulm,n_r_max)   ! Phase field

      !-- Local variables
      character(len=72) :: filename
      integer :: lm00, n_r, n_r_phase, filehandle, l, m, lm, n_t, n_t_ord
      complex(cp) :: phi_axi(l_max+1,n_r_max), phi_axi_loc(l_max+1,n_r_max)
      real(cp) :: phi_axi_g(n_r_max,nlat_padded), rmelt(n_theta_max)
      real(cp) :: ekinSr_global(n_r_max), ekinLr_global(n_r_max), volSr_global(n_r_max)
      real(cp) :: tmp(n_r_max), phi_theta(nlat_padded)
      real(cp) :: ekinL, ekinS, fcmb, ficb, dtTPhi, slope, intersect, volS
      real(cp) :: rphase, tphase

      !-- MPI gather on rank=0
      call gather_from_Rloc(ekinSr,ekinSr_global,0)
      call gather_from_Rloc(ekinLr,ekinLr_global,0)
      call gather_from_Rloc(volSr,volSr_global,0)

      !-- Re-arange m=0 modes and communicate them to rank=0
      do n_r=1,n_r_max
         phi_axi_loc(:,n_r)=zero
         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            m = lo_map%lm2m(lm)
            if ( m == 0 ) phi_axi_loc(l+1,n_r)=phi(lm,n_r)
         end do

#ifdef WITH_MPI
         call MPI_Reduce(phi_axi_loc(:,n_r), phi_axi(:,n_r), l_max+1, &
              &          MPI_DEF_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
#else
         phi_axi(:,n_r)=phi_axi_loc(:,n_r)
#endif
      end do

      if ( rank == 0 ) then

         !-- Get axisymmetric phase field on the grid
         do n_r=1,n_r_max
            call axi_to_spat(phi_axi(:,n_r), phi_theta)
            phi_axi_g(n_r,:)=phi_theta(:)
         end do

         !-- Now get the melting points for each colatitude and compute a linear
         !-- interpolation to get rmelt(theta)
         do n_t=1,n_theta_max
            n_r_phase=2
            do n_r=2,n_r_max
               if ( phi_axi_g(n_r,n_t) < half .and. phi_axi_g(n_r-1,n_t) > half ) then
                  n_r_phase=n_r
                  exit
               end if
            end do
            n_t_ord=n_theta_cal2ord(n_t)
            if ( n_r_phase /= 2 ) then
               slope=(phi_axi_g(n_r_phase,n_t)-phi_axi_g(n_r_phase-1,n_t)) / &
               &     (r(n_r_phase)-r(n_r_phase-1))
               intersect=phi_axi_g(n_r_phase,n_t)-slope*r(n_r_phase)
               rmelt(n_t_ord)=(half-intersect)/slope
            else
               rmelt(n_t_ord)=r_cmb
            end if
         end do

         !-- Now save the melting line into a binary file
         if ( l_save_out ) then
            open(newunit=n_rmelt_file, file=rmelt_file, status='unknown', &
            &    position='append', form='unformatted')
         end if
         !-- Write header when first called
         if ( nLogs == 1 ) then
            write(n_rmelt_file) n_theta_max
            write(n_rmelt_file) real(theta_ord(:),outp)
         end if
         write(n_rmelt_file) real(time,outp), real(rmelt(:),outp)
         if ( l_save_out ) close(n_rmelt_file)

         lm00 = lo_map%lm2(0,0) ! l=m=0
         !-- Mean phase field
         call PhiMeanR%compute(osq4pi*real(phi(lm00,:)),n_calls,timePassed,timeNorm)

         !-- Integration of kinetic energy over radius
         ekinL=eScale*rInt_R(ekinLr_global,r,rscheme_oc)
         ekinS=eScale*rInt_R(ekinSr_global,r,rscheme_oc)

         !-- Get the volume of the solid phase
         volS=eScale*rInt_R(volSr_global,r,rscheme_oc)

         !-- Fluxes
         fcmb = -opr * real(ds(lm00,n_r_cmb))*osq4pi*four*pi*r_cmb**2*kappa(n_r_cmb)
         ficb = -opr * real(ds(lm00,n_r_icb))*osq4pi*four*pi*r_icb**2*kappa(n_r_icb)

         !-- Integration of T-St*Phi
         tmp(:)=osq4pi*(real(s(lm00,:))-stef*real(phi(lm00,:)))*r(:)*r(:)
         TPhiOld=TPhi
         TPhi   =four*pi*rInt_R(tmp,r,rscheme_oc)
         dtTPhi =(TPhi-TPhiOld)/timePassed

         !-- Determine the radial level where \phi=0.5
         tmp(:)=osq4pi*real(phi(lm00,:)) ! Reuse tmp as work array
         n_r_phase=2
         do n_r=2,n_r_max
            if ( tmp(n_r) < half .and. tmp(n_r-1) > half ) then
               n_r_phase=n_r
            end if
         end do

         !-- Linear interpolation of melting point
         if ( n_r_phase /= 2 ) then
            slope=(tmp(n_r_phase)-tmp(n_r_phase-1))/(r(n_r_phase)-r(n_r_phase-1))
            intersect=tmp(n_r_phase)-slope*r(n_r_phase)
            rphase=(half-intersect)/slope
            tmp(:)=osq4pi*real(s(lm00,:)) ! Reuse tmp as work array
            slope=(tmp(n_r_phase)-tmp(n_r_phase-1))/(r(n_r_phase)-r(n_r_phase-1))
            intersect=tmp(n_r_phase)-slope*r(n_r_phase)
            tphase=slope*rphase+intersect
         else
            rphase=r_cmb
            tphase=osq4pi*real(s(lm00,n_r_cmb))
         end if

         if ( nLogs > 1 ) then
            if ( l_save_out ) then
               open(newunit=n_phase_file, file=phase_file, status='unknown', &
               &    position='append')
            end if

            write(n_phase_file,'(1P,ES20.12,8ES16.8)')   &
            &     time, rphase, tphase, volS, ekinS, ekinL, fcmb, ficb, dtTPhi

            if ( l_save_out ) close(n_phase_file)
         end if

         if ( l_stop_time ) then
            call PhiMeanR%finalize_SD(timeNorm)
            filename='phiR.'//tag
            open(newunit=filehandle, file=filename, status='unknown')
            do n_r=1,n_r_max
               write(filehandle, '(ES20.10,ES15.7,ES13.5)' )                     &
               &     r(n_r),round_off(PhiMeanR%mean(n_r),maxval(PhiMeanR%mean)), &
               &     round_off(PhiMeanR%SD(n_r),maxval(PhiMeanR%SD))
            end do

            close(filehandle)
         end if

      end if

   end subroutine outPhase
!----------------------------------------------------------------------------------
   subroutine get_hemi(vr,vt,vp,nR,field)
      !
      !   This subroutine is used to compute kinetic or magnetic energy
      !   in Northern or Southern hemipshere.
      !

      !-- Input of variables
      integer,          intent(in) :: nR ! radial level
      real(cp),         intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      character(len=1), intent(in) :: field

      !-- Output variables:

      !-- Local variables:
      real(cp) :: enAS(2) ! energy in North/South hemi at radius nR
      real(cp) :: vrabsAS(2)! abs(vr or Br) in North/South hemi at radius nR
      real(cp) :: en, vrabs, phiNorm, fac
      integer :: nTheta, nPhi, nTh

      enAS(:)   =0.0_cp
      vrabsAS(:)=0.0_cp
      phiNorm=two*pi/real(n_phi_max,cp)
      if ( field == 'V' ) then
         fac = orho1(nR)
      else if ( field == 'B' ) then
         fac = one
      end if
      !--- Helicity:
      !$omp parallel do default(shared)   &
      !$omp& private(nTheta,vrabs,en,nTh) &
      !$omp& reduction(+:enAS,vrabsAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nTh=n_theta_cal2ord(nTheta)
            vrabs=fac*abs(vr(nTheta,nPhi))
            en   =half*fac*(                                                &
            &                     or2(nR)*vr(nTheta,nPhi)*vr(nTheta,nPhi) + &
            &      O_sin_theta_E2(nTheta)*vt(nTheta,nPhi)*vt(nTheta,nPhi) + &
            &      O_sin_theta_E2(nTheta)*vp(nTheta,nPhi)*vp(nTheta,nPhi) )

            if ( nTh <= n_theta_max/2 ) then ! Northern Hemisphere
               enAS(1)   =enAS(1) +phiNorm*gauss(nTheta)*en
               vrabsAS(1)=vrabsAS(1) +phiNorm*gauss(nTheta)*vrabs
            else
               enAS(2)   =enAS(2) +phiNorm*gauss(nTheta)*en
               vrabsAS(2)=vrabsAS(2) +phiNorm*gauss(nTheta)*vrabs
            end if
         end do
      end do
      !$omp end parallel do

      if ( field == 'V' ) then
         hemi_ekin_r(nR,:) =enAS(:)
         hemi_vrabs_r(nR,:)=vrabsAS(:)
      else if ( field == 'B' ) then
         hemi_emag_r(nR,:) =enAS(:)
         hemi_brabs_r(nR,:)=vrabsAS(:)
      end if

   end subroutine get_hemi
!----------------------------------------------------------------------------------
   subroutine get_helicity(vr,vt,vp,cvr,dvrdt,dvrdp,dvtdr,dvpdr,nR)
      !
      ! This subroutine calculates axisymmetric and non-axisymmetric contributions to
      ! kinetic helicity and squared helicity.
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: cvr(:,:),dvrdt(:,:),dvrdp(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvpdr(:,:)

      !-- Local variables:
      integer :: nTheta,nPhi,nTh
      real(cp) :: Helna,Hel,phiNorm
      real(cp) :: HelAS(2), Hel2AS(2), HelnaAS(2), Helna2AS(2), HelEAAS
      real(cp) :: vrna,vtna,vpna,cvrna,dvrdtna,dvrdpna,dvtdrna,dvpdrna
      real(cp) :: vras(n_theta_max),vtas(n_theta_max),vpas(n_theta_max)
      real(cp) :: cvras(n_theta_max),dvrdtas(n_theta_max),dvrdpas(n_theta_max)
      real(cp) :: dvtdras(n_theta_max),dvpdras(n_theta_max)

      !-- Remark: 2pi not used the normalization below
      !-- this is why we have a 2pi factor after radial integration
      !-- in the subroutine outHelicity()
      phiNorm=one/real(n_phi_max,cp)
      HelAS(:)   =0.0_cp
      Hel2AS(:)  =0.0_cp
      HelnaAS(:) =0.0_cp
      Helna2AS(:)=0.0_cp
      HelEAAS    =0.0_cp

      vras(:)   =0.0_cp
      cvras(:)  =0.0_cp
      vtas(:)   =0.0_cp
      vpas(:)   =0.0_cp
      dvrdpas(:)=0.0_cp
      dvpdras(:)=0.0_cp
      dvtdras(:)=0.0_cp
      dvrdtas(:)=0.0_cp
      do nPhi=1,n_phi_max
         vras(:)   =vras(:)   +   vr(1:n_theta_max,nPhi)
         cvras(:)  =cvras(:)  +  cvr(1:n_theta_max,nPhi)
         vtas(:)   =vtas(:)   +   vt(1:n_theta_max,nPhi)
         vpas(:)   =vpas(:)   +   vp(1:n_theta_max,nPhi)
         dvrdpas(:)=dvrdpas(:)+dvrdp(1:n_theta_max,nPhi)
         dvpdras(:)=dvpdras(:)+dvpdr(1:n_theta_max,nPhi)
         dvtdras(:)=dvtdras(:)+dvtdr(1:n_theta_max,nPhi)
         dvrdtas(:)=dvrdtas(:)+dvrdt(1:n_theta_max,nPhi)
      end do
      vras(:)   =vras(:)   *phiNorm
      cvras(:)  =cvras(:)  *phiNorm
      vtas(:)   =vtas(:)   *phiNorm
      vpas(:)   =vpas(:)   *phiNorm
      dvrdpas(:)=dvrdpas(:)*phiNorm
      dvpdras(:)=dvpdras(:)*phiNorm
      dvtdras(:)=dvtdras(:)*phiNorm
      dvrdtas(:)=dvrdtas(:)*phiNorm

      !--- Helicity:
      !$omp parallel do default(shared)                     &
      !$omp& private(nTheta, nPhi, nTh, Hel, Helna)         &
      !$omp& private(vrna, cvrna, vtna, vpna)               &
      !$omp& private(dvrdpna, dvpdrna, dvtdrna, dvrdtna)    &
      !$omp& reduction(+:HelAS,Hel2AS,HelnaAS,Helna2AS,HelEAAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nTh=n_theta_cal2ord(nTheta)
            vrna   =   vr(nTheta,nPhi)-vras(nTheta)
            cvrna  =  cvr(nTheta,nPhi)-cvras(nTheta)
            vtna   =   vt(nTheta,nPhi)-vtas(nTheta)
            vpna   =   vp(nTheta,nPhi)-vpas(nTheta)
            dvrdpna=dvrdp(nTheta,nPhi)-dvrdpas(nTheta)
            dvpdrna=dvpdr(nTheta,nPhi)-beta(nR)*vp(nTheta,nPhi) &
            &       -dvpdras(nTheta)+beta(nR)*vpas(nTheta)
            dvtdrna=dvtdr(nTheta,nPhi)-beta(nR)*vt(nTheta,nPhi) &
            &       -dvtdras(nTheta)+beta(nR)*vtas(nTheta)
            dvrdtna=dvrdt(nTheta,nPhi)-dvrdtas(nTheta)
            Hel=or4(nR)*orho2(nR)*vr(nTheta,nPhi)*cvr(nTheta,nPhi) +  &
            &             or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                                       vt(nTheta,nPhi) * &
            &                          ( or2(nR)*dvrdp(nTheta,nPhi) - &
            &                                    dvpdr(nTheta,nPhi) + &
            &                         beta(nR)*   vp(nTheta,nPhi) ) + &
            &                                       vp(nTheta,nPhi) * &
            &                          (         dvtdr(nTheta,nPhi) - &
            &                           beta(nR)*   vt(nTheta,nPhi) - &
            &                            or2(nR)*dvrdt(nTheta,nPhi) ) )
            Helna=                      or4(nR)*orho2(nR)*vrna*cvrna + &
            &              or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                       vtna*( or2(nR)*dvrdpna-dvpdrna ) + &
            &                       vpna*( dvtdrna-or2(nR)*dvrdtna ) )

            if ( nTh <= n_theta_max/2 ) then ! Northern Hemisphere
               HelAS(1)   =HelAS(1) +phiNorm*gauss(nTheta)*Hel
               Hel2AS(1)  =Hel2AS(1)+phiNorm*gauss(nTheta)*Hel*Hel
               HelnaAS(1) =HelnaAS(1) +phiNorm*gauss(nTheta)*Helna
               Helna2AS(1)=Helna2AS(1)+phiNorm*gauss(nTheta)*Helna*Helna
               HelEAAS    =HelEAAS +phiNorm*gauss(nTheta)*Hel
            else
               HelAS(2)   =HelAS(2) +phiNorm*gauss(nTheta)*Hel
               Hel2AS(2)  =Hel2AS(2)+phiNorm*gauss(nTheta)*Hel*Hel
               HelnaAS(2) =HelnaAS(2) +phiNorm*gauss(nTheta)*Helna
               Helna2AS(2)=Helna2AS(2)+phiNorm*gauss(nTheta)*Helna*Helna
               HelEAAS    =HelEAAS -phiNorm*gauss(nTheta)*Hel
            end if
         end do
      end do
      !$omp end parallel do

      HelASr(nR,:)   =HelAS(:)
      Hel2ASr(nR,:)  =Hel2AS(:)
      HelnaASr(nR,:) =HelnaAS(:)
      Helna2ASr(nR,:)=Helna2AS(:)
      HelEAASr(nR)   =HelEAAS

   end subroutine get_helicity
!----------------------------------------------------------------------------------
   subroutine get_ekin_solid_liquid(vr,vt,vp,phi,nR)
      !
      ! This subroutine computes the kinetic energy content in the solid
      ! and in the liquid phase when phase field is employed.
      !

      !-- Input variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:),phi(:,:)

      !-- Output variables:

      !-- Local variables:
      real(cp) :: phiNorm, ekin
      real(cp) :: ekinS ! Kinetic energy in the solid phase
      real(cp) :: ekinL ! Kinetic energy in the liquid phase
      real(cp) :: volS  ! volume of the solid
      integer :: nTheta,nPhi

      phiNorm=two*pi/real(n_phi_max,cp)
      ekinL=0.0_cp
      ekinS=0.0_cp
      volS =0.0_cp

      !$omp parallel do default(shared) &
      !$omp& private(nTheta,nPhi,ekin)  &
      !$omp& reduction(+:ekinS,ekinL,volS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            ekin = half*orho1(nR)*(                                         &
            &          or2(nR)*           vr(nTheta,nPhi)*vr(nTheta,nPhi) + &
            &      O_sin_theta_E2(nTheta)*vt(nTheta,nPhi)*vt(nTheta,nPhi) + &
            &      O_sin_theta_E2(nTheta)*vp(nTheta,nPhi)*vp(nTheta,nPhi) )

            if ( phi(nTheta,nPhi) >= half ) then
               ekinS=ekinS+phiNorm*gauss(nTheta)*ekin
               volS =volS +phiNorm*gauss(nTheta)*r(nR)*r(nR)
            else
               ekinL=ekinL+phiNorm*gauss(nTheta)*ekin
            end if
         end do
      end do
      !$omp end parallel do

      ekinSr(nR)=ekinS
      ekinLr(nR)=ekinL
      volSr(nR) =volS

   end subroutine get_ekin_solid_liquid
!----------------------------------------------------------------------------------
   subroutine get_onset(time, w, dt, l_log, nLogs)
      !
      ! This subroutine is used to estimate the growth rate and the drifting
      ! frequencies of a series of m modes for both equatorially-symmetric
      ! and equatorially-antisymmetric modes. This is used to compute the critical
      ! Rayleigh number. This makes uses of the radial integration of the poloidal
      ! potential at different (l,m) tuples.
      !

      real(cp),    intent(in) :: time  ! Time
      real(cp),    intent(in) :: dt    ! Timestep size
      logical,     intent(in) :: l_log ! Do we need to store outputs
      integer,     intent(in) :: nLogs ! Do not write at the first time step
      complex(cp), intent(in) :: w(llm:ulm, n_r_max) ! Poloidal potential

      !-- Local variables
      character(len=3) :: ad
      complex(cp) :: coeff, tau(llm:ulm)
      real(cp) :: tmpr(n_r_max), tmpi(n_r_max)
      real(cp) :: facR, facI
      complex(cp), allocatable :: tau_glob(:)
      integer :: lm, l, m

      do lm=llm,ulm
         tmpr(:) = real(w(lm,:)) * r * r
         tmpi(:) = aimag(w(lm,:)) * r * r
         facR = rInt_R(tmpr, r, rscheme_oc)
         facI = rInt_R(tmpi, r, rscheme_oc)
         coeff = cmplx(facR, facI, kind=cp)
         if ( abs(coeff) > 0.0_cp .and. abs(coeff_old(lm)) > 0.0_cp ) then
            tau(lm) = cmplx((abs(coeff) - abs(coeff_old(lm)))/abs(coeff), &
            &               aimag((coeff - coeff_old(lm))/coeff), kind=cp)/dt
         else
            tau(lm) = zero
         end if
         coeff_old(lm) = coeff
      end do

      if ( l_log ) then
         if ( rank == 0 ) then
            allocate( tau_glob(lm_max) )
         else
            allocate( tau_glob(1) )
         end if
         call gather_from_lo_to_rank0(tau, tau_glob)
      end if

      if ( rank == 0 .and. l_log .and. nLogs > 1 ) then
         if ( l_save_out ) then
            open(newunit=n_growth_sym_file, file=sym_file, status='unknown', &
            &    position='append')
            open(newunit=n_growth_asym_file, file=asym_file, status='unknown',   &
            &    position='append')
            open(newunit=n_drift_sym_file, file=drift_sym_file, status='unknown',   &
            &    position='append')
            open(newunit=n_drift_asym_file, file=drift_asym_file, status='unknown', &
            &    position='append')
         end if
         ad='no'
         do m=m_min,m_max,minc
            if ( m == m_max .and. (.not. l_save_out)) ad='yes'
            l =max(m,1)
            lm=lm2(l,m)
            if ( m == m_min ) then
               write(n_growth_sym_file, '(es16.8,es14.6)', advance=ad) &
               &     time, real(tau_glob(lm))
               write(n_drift_sym_file, '(es16.8,es14.6)', advance=ad) &
               &     time, aimag(tau_glob(lm))
            else
               write(n_growth_sym_file, '(es14.6)', advance=ad) real(tau_glob(lm))
               write(n_drift_sym_file, '(es14.6)', advance=ad) aimag(tau_glob(lm))
            end if
            lm=lm2(l+1,m)
            if ( m == m_min ) then
               write(n_growth_asym_file, '(es16.8,es14.6)', advance=ad) &
               &     time, real(tau_glob(lm))
               write(n_drift_asym_file, '(es16.8,es14.6)', advance=ad) &
               &     time, aimag(tau_glob(lm))
            else
               write(n_growth_asym_file, '(es14.6)', advance=ad) real(tau_glob(lm))
               write(n_drift_asym_file, '(es14.6)', advance=ad) aimag(tau_glob(lm))
            end if

         end do
         if ( l_save_out ) then
            close(n_growth_sym_file)
            close(n_growth_asym_file)
            close(n_drift_sym_file)
            close(n_drift_asym_file)
         end if
      end if

      if ( l_log ) deallocate(tau_glob)

   end subroutine get_onset
!----------------------------------------------------------------------------------
end module outMisc_mod
