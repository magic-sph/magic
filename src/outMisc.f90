module outMisc_mod
   !
   ! This module contains several subroutines that can compute and store
   ! various informations: helicity (helicity.TAG), heat transfer (heat.TAG),
   ! phase field (phase.TAG) and North/South hemisphericity of energies (hemi.TAG)
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use communications, only: gather_from_Rloc, gather_from_lo_to_rank0
   use truncation, only: n_r_max, n_theta_max, n_r_maxMag, n_phi_max, lm_max, &
       &                 m_min, m_max, minc
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop, radial_balance
   use radial_functions, only: r_icb, rscheme_oc, kappa, r_cmb,temp0, r, rho0, &
       &                       dLtemp0, dLalpha0, beta, orho1, alpha0, otemp1, &
       &                       ogrun, rscheme_oc, or2, orho2, or4
   use physical_parameters, only: ViscHeatFac, ThExpNb, opr, stef, LFfac
   use num_param, only: lScale, eScale, vScale
   use blocking, only: llm, ulm, lo_map, lm2
   use radial_der, only: get_dr
   use mean_sd, only: mean_sd_type
   use movie_data, only: n_movie_fields, n_movie_file, n_movie_field_type, &
       &                 n_movies, n_movie_field_start, frames
   use horizontal_data, only: gauss, theta_ord, n_theta_cal2ord, O_sin_theta_E2
   use logic, only: l_save_out, l_anelastic_liquid, l_heat, l_hel, l_hemi, &
       &            l_temperature_diff, l_chemical_conv, l_phase_field,    &
       &            l_mag, l_onset, l_dtphaseMovie
   use output_data, only: tag
   use constants, only: pi, vol_oc, osq4pi, sq4pi, one, two, four, half, zero
   use start_fields, only: topcond, botcond, deltacond, topxicond, botxicond, &
       &                   deltaxicond
   use useful, only: round_off, lagrange_interp
   use integration, only: rInt_R

   implicit none

   private

   type(mean_sd_type) :: TMeanR, SMeanR, PMeanR, XiMeanR, RhoMeanR, PhiMeanR
   integer :: n_heat_file, n_helicity_file, n_calls, n_phase_file
   integer :: n_rmelt_file, n_hemi_file, n_growth_sym_file, n_growth_asym_file
   integer :: n_drift_sym_file, n_drift_asym_file
   character(len=72) :: heat_file, helicity_file, phase_file, rmelt_file
   character(len=72) :: hemi_file, sym_file, asym_file, drift_sym_file
   character(len=72) :: drift_asym_file
   real(cp) :: TPhiOld, Tphi
   real(cp), allocatable :: ekinSr(:), ekinLr(:), volSr(:)
   real(cp), allocatable :: hemi_ekin_r(:,:), hemi_vrabs_r(:,:)
   real(cp), allocatable :: hemi_emag_r(:,:), hemi_brabs_r(:,:)
   real(cp), allocatable :: HelASr(:,:), Hel2ASr(:,:)
   real(cp), allocatable :: HelnaASr(:,:), Helna2ASr(:,:)
   real(cp), allocatable :: HelEAASr(:)
   real(cp), allocatable :: temp_Rloc(:,:,:), phase_Rloc(:,:,:), dtemp_Rloc(:,:,:)
   real(cp), allocatable :: temp_Ploc(:,:,:), phase_Ploc(:,:,:), dtemp_Ploc(:,:,:)
   complex(cp), allocatable :: coeff_old(:)
   integer :: nPstart ! Starting nPhi index when MPI distributed
   integer :: nPstop  ! Stoping nPhi index when MPI distributed
   type(load), allocatable :: phi_balance(:) ! phi-distributed balance
   real(cp), allocatable :: rmelt_loc(:,:) ! Melting radius (theta,phi)
   real(cp), allocatable :: dt_rmelt_loc(:,:) ! Temp. gradient at melting radius (theta,phi)

   public :: outHelicity, outHeat, initialize_outMisc_mod, finalize_outMisc_mod, &
   &         outPhase, outHemi, get_ekin_solid_liquid, get_helicity, get_hemi,   &
   &         get_onset, calc_melt_frame

contains

   subroutine initialize_outMisc_mod()
      !
      ! This subroutine handles the opening of some output diagnostic files that
      ! have to do with heat transfer, helicity, phase field or hemisphericity
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
      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_hel ) open(newunit=n_helicity_file, file=helicity_file, status='new')
         if ( l_hemi ) open(newunit=n_hemi_file, file=hemi_file, status='new')
         if ( l_heat .or. l_chemical_conv ) then
            open(newunit=n_heat_file, file=heat_file, status='new')
         end if
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

         !-- Distribute over the ranks
         allocate(phi_balance(0:n_procs-1))
         call getBlocks(phi_balance, n_phi_max, n_procs)
         nPstart = phi_balance(rank)%nStart
         nPstop = phi_balance(rank)%nStop

         allocate( phase_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
         allocate( temp_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
         allocate( phase_Ploc(n_theta_max,nPstart:nPstop,n_r_max) )
         allocate( temp_Ploc(n_theta_max,nPstart:nPstop,n_r_max) )
         bytes_allocated=bytes_allocated+4*(nRstop-nRstart+1)*n_phi_max* &
         &               n_theta_max*SIZEOF_DEF_REAL
         phase_Rloc(:,:,:)=0.0_cp
         temp_Rloc(:,:,:) =0.0_cp
         phase_Ploc(:,:,:)=0.0_cp
         temp_Ploc(:,:,:) =0.0_cp
         allocate( rmelt_loc(n_theta_max,nPstart:nPstop) )
         bytes_allocated=bytes_allocated+(nPstop-nPstart+1)*n_theta_max* &
         &               SIZEOF_DEF_REAL
         rmelt_loc(:,:)   =0.0_cp

         if ( l_dtphaseMovie ) then
            allocate( dtemp_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
            allocate( dtemp_Ploc(n_theta_max,nPstart:nPstop,n_r_max) )
            bytes_allocated=bytes_allocated+2*(nRstop-nRstart+1)*n_phi_max* &
            &               n_theta_max*SIZEOF_DEF_REAL
            dtemp_Rloc(:,:,:)=0.0_cp
            dtemp_Ploc(:,:,:)=0.0_cp
            allocate( dt_rmelt_loc(n_theta_max,nPstart:nPstop) )
            bytes_allocated=bytes_allocated+(nPstop-nPstart+1)*n_theta_max* &
            &               SIZEOF_DEF_REAL
            dt_rmelt_loc(:,:)=0.0_cp
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
      ! heat.TAG, hel.TAG, hemi.TAG and phase.TAG
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
         if ( l_dtphaseMovie ) deallocate( dtemp_Rloc, dtemp_Ploc, dt_rmelt_loc)
         deallocate( temp_Rloc, phase_Rloc, temp_Ploc, phase_Ploc, phi_balance )
         deallocate( ekinSr, ekinLr, volSr, rmelt_loc )
         call PhiMeanR%finalize()
      end if

      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_hel ) close(n_helicity_file)
         if ( l_hemi ) close(n_hemi_file)
         if ( l_heat .or. l_chemical_conv ) close(n_heat_file)
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
!----------------------------------------------------------------------------------
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
      integer :: lm00, n_r, filehandle, n_t, n_p, n_t_ord
      real(cp) :: ekinSr_global(n_r_max), ekinLr_global(n_r_max), volSr_global(n_r_max)
      real(cp) :: phase(n_r_max), temp(n_r_max)
      real(cp) :: ekinL, ekinS, fcmb, ficb, dtTPhi, volS, phase_max, phase_min
      real(cp) :: rphase, tphase, norm, rmelt_mean_loc, tmelt_mean_loc
      real(cp) :: rmelt_mean, tmelt_mean, rmelt_max, rmelt_min
      real(cp) :: phase_max_loc, phase_min_loc
      real(cp) :: rmelt_min_loc, rmelt_max_loc, tmelt_loc
      real(cp) :: rmelt_axi(n_theta_max), rmelt_axi_loc(n_theta_max)

      !-- MPI gather on rank=0
      call gather_from_Rloc(ekinSr,ekinSr_global,0)
      call gather_from_Rloc(ekinLr,ekinLr_global,0)
      call gather_from_Rloc(volSr,volSr_global,0)

      phase_max_loc = maxval(phase_Rloc)
      phase_min_loc = minval(phase_RLoc)
#ifdef WITH_MPI
      call MPI_Reduce(phase_max_loc, phase_max, 1, MPI_DEF_REAL, MPI_MAX, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(phase_min_loc, phase_min, 1, MPI_DEF_REAL, MPI_MIN, &
           &          0, MPI_COMM_WORLD, ierr)
#else
      phase_max=phase_max_loc
      phase_min=phase_min_loc
#endif

      !-- MPI transpose
      call transp_R2Phi(temp_Rloc, temp_Ploc)
      if ( l_dtphaseMovie ) call transp_R2Phi(dtemp_Rloc, dtemp_Ploc)
      call transp_R2Phi(phase_Rloc, phase_Ploc)

      rmelt_axi_loc(:)=0.0_cp
      rmelt_mean_loc=0.0_cp
      tmelt_mean_loc=0.0_cp
      norm=half/n_phi_max
      do n_p=nPstart,nPstop
         do n_t=1,n_theta_max
            n_t_ord=n_theta_cal2ord(n_t)
            if ( l_dtphaseMovie ) then
               call get_rmelt_tmelt(phase_Ploc(n_t,n_p,:), temp_Ploc(n_t,n_p,:), &
                    &               rmelt_loc(n_t,n_p), tmelt_loc,               &
                    &               dtemp_Ploc(n_t,n_p,:), dt_rmelt_loc(n_t,n_p))
            else
               call get_rmelt_tmelt(phase_Ploc(n_t,n_p,:), temp_Ploc(n_t,n_p,:), &
                    &               rmelt_loc(n_t,n_p), tmelt_loc)
            end if
            rmelt_axi_loc(n_t_ord)=rmelt_axi_loc(n_t_ord)+two*norm*rmelt_loc(n_t,n_p)
            rmelt_mean_loc=rmelt_mean_loc+gauss(n_t)*norm*rmelt_loc(n_t,n_p)
            tmelt_mean_loc=tmelt_mean_loc+gauss(n_t)*norm*tmelt_loc
         end do
      end do
      rmelt_max_loc=maxval(rmelt_loc)
      rmelt_min_loc=minval(rmelt_loc)

#ifdef WITH_MPI
      call MPI_Reduce(rmelt_axi_loc, rmelt_axi, n_theta_max, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(rmelt_mean_loc, rmelt_mean, 1, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(tmelt_mean_loc, tmelt_mean, 1, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(rmelt_min_loc, rmelt_min, 1, MPI_DEF_REAL, MPI_MIN, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(rmelt_max_loc, rmelt_max, 1, MPI_DEF_REAL, MPI_MAX, &
           &          0, MPI_COMM_WORLD, ierr)
#else
      rmelt_axi(:)=rmelt_axi_loc(:)
      rmelt_mean=rmelt_mean_loc
      tmelt_mean=tmelt_mean_loc
      rmelt_max =rmelt_max_loc
      rmelt_min =rmelt_min_loc
#endif

      if ( rank == 0 ) then

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
         write(n_rmelt_file) real(time,outp), real(rmelt_axi(:),outp)
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
         phase(:)=osq4pi*(real(s(lm00,:))-stef*real(phi(lm00,:)))*r(:)*r(:)
         TPhiOld=TPhi
         TPhi   =four*pi*rInt_R(phase,r,rscheme_oc)
         dtTPhi =(TPhi-TPhiOld)/timePassed

         !-- Determine the radial level where \phi=0.5
         phase(:)=osq4pi*real(phi(lm00,:))
         temp(:) =osq4pi*real(s(lm00,:))
         call get_rmelt_tmelt(phase, temp, rphase, tphase)

         if ( nLogs > 1 ) then
            if ( l_save_out ) then
               open(newunit=n_phase_file, file=phase_file, status='unknown', &
               &    position='append')
            end if

            write(n_phase_file,'(1P,ES20.12,11ES16.8,3ES13.5)')    &
            &     time, rphase, tphase, rmelt_mean, tmelt_mean,    &
            &     rmelt_min, rmelt_max, volS, ekinS, ekinL, fcmb,  &
            &     ficb, dtTPhi, phase_min, phase_max

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
   subroutine get_ekin_solid_liquid(s,ds,vr,vt,vp,phi,nR)
      !
      ! This subroutine computes the kinetic energy content in the solid
      ! and in the liquid phase when phase field is employed.
      !

      !-- Input variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: s(:,:),ds(:,:),vr(:,:),vt(:,:),vp(:,:),phi(:,:)

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

            phase_Rloc(nTheta,nPhi,nR)=phi(nTheta,nPhi)
            temp_Rloc(nTheta,nPhi,nR) =s(nTheta,nPhi)
            if ( l_dtphaseMovie ) dtemp_Rloc(nTheta,nPhi,nR)=ds(nTheta,nPhi)
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
   subroutine get_rmelt_tmelt(phase, temp, rphase, tphase, dtemp, dtphase)
      !
      ! This subroutine determines the melting point by approximating it by
      ! the radius where phi=0.5. It returns the radius and the temperature.
      ! It computes a 4th order Lagrangian interpolation between the four closest
      ! radii.
      !

      !-- Input variables
      real(cp), intent(in) :: phase(:) ! Phase field
      real(cp), intent(in) :: temp(:)  ! Temperature
      real(cp), optional, intent(in) :: dtemp(:) ! Temperature gradient

      !-- Output variables
      real(cp), intent(out) :: rphase ! Radius of the melting point
      real(cp), intent(out) :: tphase ! Temperature of the melting point
      real(cp), optional, intent(out) :: dtphase ! Temperature gradient at rm

      !-- Local variables
      integer :: n_r, n_r_phase, n_r_start, n_r_stop
      real(cp) :: x(4), y(4)

      !-- Determine the radial level where \phi=0.5
      n_r_phase=1
      do n_r=2,n_r_max
         if ( phase(n_r) < half .and. phase(n_r-1) > half ) then
            n_r_phase=n_r
         end if
      end do

      !-- 4th order Lagrangian interpolation of melting point
      if ( n_r_phase == n_r_cmb ) then
         rphase=r_cmb
         tphase=temp(n_r_cmb)
      else
         if ( n_r_phase == n_r_cmb+1 ) then
            n_r_start=n_r_phase-1
            n_r_stop =n_r_phase+2
         else if ( n_r_phase == n_r_icb ) then
            n_r_start=n_r_phase-3
            n_r_stop =n_r_phase
         else
            n_r_start=n_r_phase-2
            n_r_stop =n_r_phase+1
         end if
         x(:)=phase(n_r_start:n_r_stop)
         y(:)=r(n_r_start:n_r_stop)
         rphase=lagrange_interp(x,half,y)
         x(:)=r(n_r_start:n_r_stop)
         y(:)=temp(n_r_start:n_r_stop)
         tphase=lagrange_interp(x,rphase,y)
         if ( present(dtemp) ) then
            y(:)=dtemp(n_r_start:n_r_stop)
            dtphase=lagrange_interp(x,rphase,y)
         end if
      end if

   end subroutine get_rmelt_tmelt
!------------------------------------------------------------------------------------
   subroutine calc_melt_frame()
      !
      ! This subroutine handles the computation and the writing of the movie
      ! which contains rmelt(phi,theta) or dt(phi,theta,r=rmelt) as a function of time.
      !

      !-- Local variables
      integer :: n_movie, n_field, n_fields, n_out, n_field_type, n_store_last
      integer :: n_theta, n_phi, n_o
      real(cp) :: data(n_phi_max,n_theta_max)

      do n_movie=1,n_movies
         n_fields=n_movie_fields(n_movie)
         n_out   =n_movie_file(n_movie)

         do n_field=1,n_fields
            n_field_type=n_movie_field_type(n_field,n_movie)
            n_store_last=n_movie_field_start(n_field,n_movie)-1

            !-- MPI gather
            if ( n_field_type == 117 ) then ! melt radius
               call gather_Ploc(rmelt_loc, data)
            else if ( n_field_type == 118 ) then ! radial derivative of temp. along interface
               call gather_Ploc(dt_rmelt_loc, data)
            else
               cycle
            end if

            !-- Store outputs in frames(*)
            if ( rank == 0 ) then
               do n_phi=1,n_phi_max
                  do n_theta=1,n_theta_max
                     n_o=n_store_last+(n_theta-1)*n_phi_max
                     frames(n_phi+n_o)=real(data(n_phi,n_theta),kind=outp)
                  end do
               end do
            end if
         end do
      end do

   end subroutine calc_melt_frame
!------------------------------------------------------------------------------------
   subroutine transp_R2Phi(arr_Rloc, arr_Ploc)
      !
      ! This subroutine is used to compute a MPI transpose between a R-distributed
      ! array and a Phi-distributed array
      !

      !-- Input array
      real(cp), intent(in) :: arr_Rloc(n_theta_max,n_phi_max,nRstart:nRstop)

      !-- Output array
      real(cp), intent(out) :: arr_Ploc(n_theta_max,nPstart:nPstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_t, n_p
#ifdef WITH_MPI
      integer, allocatable :: rcounts(:), scounts(:), rdisp(:), sdisp(:)
      real(cp), allocatable :: sbuff(:), rbuff(:)
      integer :: p, ii, my_phi_counts

      !-- Set displacements vectors and buffer sizes
      allocate( rcounts(0:n_procs-1), scounts(0:n_procs-1) )
      allocate( rdisp(0:n_procs-1), sdisp(0:n_procs-1) )
      do p=0,n_procs-1
         my_phi_counts=phi_balance(p)%n_per_rank
         scounts(p)=nR_per_rank*my_phi_counts*n_theta_max
         rcounts(p)=radial_balance(p)%n_per_rank*(nPStop-nPStart+1)*n_theta_max
      end do

      rdisp(0)=0
      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do
      allocate( sbuff(sum(scounts)), rbuff(sum(rcounts)) )
      sbuff(:)=0.0_cp
      rbuff(:)=0.0_cp

      !-- Prepare buffer
      do p=0,n_procs-1
         ii=sdisp(p)+1
         do n_r=nRstart,nRstop
            do n_p=phi_balance(p)%nStart,phi_balance(p)%nStop
               do n_t=1,n_theta_max
                  sbuff(ii)=arr_Rloc(n_t,n_p,n_r)
                  ii=ii+1
               end do
            end do
         end do
      end do

      !-- All to all
      call MPI_Alltoallv(sbuff, scounts, sdisp, MPI_DEF_REAL, &
           &             rbuff, rcounts, rdisp, MPI_DEF_REAL, &
           &             MPI_COMM_WORLD, ierr)

      !-- Reassemble array
      do p=0,n_procs-1
         ii=rdisp(p)+1
         do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
            do n_p=nPstart,nPstop
               do n_t=1,n_theta_max
                  arr_Ploc(n_t,n_p,n_r)=rbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end do

      !-- Clear memory from temporary arrays
      deallocate( rcounts, scounts, rdisp, sdisp, rbuff, sbuff )
#else
      arr_Ploc(:,:,:)=arr_Rloc(:,:,:)
#endif

   end subroutine transp_R2Phi
!------------------------------------------------------------------------------------
   subroutine gather_Ploc(arr_Ploc, arr_full)
      !
      ! This subroutine gathers and transpose a phi-distributed array on rank0.
      !

      !-- Input variable
      real(cp), intent(in) :: arr_Ploc(n_theta_max,nPstart:nPstop)

      !-- Output variable
      real(cp), intent(out) :: arr_full(n_phi_max,n_theta_max)

      !-- Local variables
      real(cp) :: tmp(n_theta_max,n_phi_max)
      integer :: n_t, n_p, n_t_ord
#ifdef WITH_MPI
      integer :: rcounts(0:n_procs-1), rdisp(0:n_procs-1), scount
      integer :: p

      scount = (nPstop-nPstart+1)*n_theta_max
      do p=0,n_procs-1
         rcounts(p)=phi_balance(p)%n_per_rank*n_theta_max
      end do
      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      call MPI_GatherV(arr_Ploc, scount, MPI_DEF_REAL, tmp, rcounts, rdisp, &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         do n_p=1,n_phi_max
            do n_t=1,n_theta_max
               n_t_ord=n_theta_cal2ord(n_t)
               arr_full(n_p,n_t_ord)=tmp(n_t,n_p)
            end do
         end do
      end if

#else
      do n_p=1,n_phi_max
         do n_t=1,n_theta_max
            n_t_ord=n_theta_cal2ord(n_t)
            arr_full(n_p,n_t_ord)=arr_Ploc(n_t,n_p)
         end do
      end do
#endif

   end subroutine gather_Ploc
!------------------------------------------------------------------------------------
end module outMisc_mod
