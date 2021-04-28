module RMS
   !
   ! This module contains the calculation of the RMS force balance and induction
   ! terms.
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use blocking, only: st_map, lo_map, lm2, lm2m, llm, ulm, llmMag, ulmMag, &
       &               lm2lmA, lm2lmP, lm2l, lm2lmS
   use finite_differences, only: type_fd
   use chebyshev, only: type_cheb_odd
   use radial_scheme, only: type_rscheme
   use truncation, only: n_r_max, n_cheb_max, n_r_maxMag, lm_max, lm_maxMag, &
       &                 l_max, n_phi_max, n_theta_max, minc, n_r_max_dtB,   &
       &                 lm_max_dtB, fd_ratio, fd_stretch, lmP_max, nlat_padded
   use physical_parameters, only: ra, ek, pr, prmag, radratio, CorFac, n_r_LCR, &
       &                          BuoFac, ChemFac, ThExpNb, ViscHeatFac
   use radial_data, only: nRstop, nRstart, radial_balance, nRstartMag, nRstopMag
   use radial_functions, only: rscheme_oc, r, r_cmb, r_icb, or1, or2, or3, or4, &
       &                       rho0, rgrav, beta, dLvisc, dbeta, ogrun, alpha0, &
       &                       temp0, visc, l_R
   use logic, only: l_save_out, l_heat, l_chemical_conv, l_conv_nl, l_mag_LF,    &
       &            l_conv, l_corr, l_mag, l_finite_diff, l_newmap, l_2D_RMS,    &
       &            l_parallel_solve, l_mag_par_solve, l_adv_curl, l_double_curl,&
       &            l_anelastic_liquid, l_mag_nl
   use num_param, only: tScale, alph1, alph2
   use horizontal_data, only: phi, theta_ord, cosTheta, sinTheta, O_sin_theta_E2,  &
       &                      cosn_theta_E2, O_sin_theta, dTheta2A, dPhi, dTheta2S,&
       &                      dLh, hdif_V
   use constants, only: zero, one, half, four, third, vol_oc, pi, two, three
   use integration, only: rInt_R
   use radial_der, only: get_dr, get_dr_Rloc
   use output_data, only: rDea, rCut, tag, runid
   use cosine_transform_odd
   use RMS_helpers, only: hInt2dPol, get_PolTorRms, hInt2dPolLM, hIntRms
   use dtB_mod, only: PdifLM_LMloc, TdifLM_LMloc, PstrLM_LMloc, PadvLM_LMloc, &
       &              TadvLM_LMloc, TstrLM_LMloc, TomeLM_LMloc
   use useful, only: abortRun
   use mean_sd, only: mean_sd_type, mean_sd_2D_type
   use time_schemes, only: type_tscheme
   use sht, only: spat_to_sphertor, spat_to_qst, scal_to_SH, scal_to_grad_spat

   implicit none

   private

   class(type_rscheme), pointer :: rscheme_RMS
   integer, public :: n_r_maxC     ! Number of radial points
   integer, public :: n_cheb_maxC  ! Number of Chebyshevs
   integer, public :: nCut         ! Number of points for the cut-off

   real(cp), allocatable :: rC(:)        ! Cut-off radii
   real(cp), public, allocatable :: dr_facC(:)

   real(cp), public, allocatable :: dtBPol2hInt(:,:), dtBTor2hInt(:,:)
   complex(cp), public, allocatable :: dtBPolLMr(:,:)

   real(cp), public, allocatable :: DifPol2hInt(:,:), DifTor2hInt(:,:)
   complex(cp), public, allocatable :: DifPolLMr(:,:)

   real(cp), allocatable :: Adv2hInt(:,:), Cor2hInt(:,:), LF2hInt(:,:), Buo_temp2hInt(:,:)
   real(cp), allocatable :: Buo_xi2hInt(:,:), Pre2hInt(:,:), Geo2hInt(:,:)
   real(cp), allocatable :: Mag2hInt(:,:), Arc2hInt(:,:), ArcMag2hInt(:,:), CIA2hInt(:,:)
   real(cp), allocatable :: CLF2hInt(:,:), PLF2hInt(:,:), Iner2hInt(:,:)

   !-- Arrays on the physical grid
   real(cp), allocatable :: Advt2(:,:), Advp2(:,:), dpdtc(:,:), dpdpc(:,:)
   real(cp), allocatable :: CFt2(:,:), CFp2(:,:), LFt2(:,:), LFp2(:,:)
   real(cp), allocatable :: dpkindrc(:,:), dtVr(:,:), dtVt(:,:), dtVp(:,:)
   real(cp), allocatable :: vr_old(:,:,:), vt_old(:,:,:), vp_old(:,:,:)

   complex(cp), allocatable :: dtVrLM(:), dtVtLM(:), dtVpLM(:)
   complex(cp), allocatable :: dpkindrLM(:), Advt2LM(:), Advp2LM(:)
   complex(cp), allocatable :: PFt2LM(:), PFp2LM(:), LFt2LM(:), LFp2LM(:)
   complex(cp), allocatable :: CFt2LM(:), CFp2LM(:), LFrLM(:)

   !-- Time-averaged spectra
   type(mean_sd_type) :: InerRmsL, CorRmsL, LFRmsL, AdvRmsL
   type(mean_sd_type) :: DifRmsL, Buo_tempRmsL, Buo_xiRmsL,PreRmsL, GeoRmsL
   type(mean_sd_type) :: MagRmsL, ArcRmsL, ArcMagRmsL
   type(mean_sd_type) :: CIARmsL, CLFRmsL, PLFRmsL

   type(mean_sd_2D_type) :: CorRmsLnR, AdvRmsLnR, LFRmsLnR
   type(mean_sd_2D_type) :: Buo_tempRmsLnR, Buo_xiRmsLnR
   type(mean_sd_2D_type) :: PreRmsLnR, DifRmsLnR, InerRmsLnR, GeoRmsLnR
   type(mean_sd_2D_type) :: MagRmsLnR, ArcRmsLnR, ArcMagRmsLnR
   type(mean_sd_2D_type) :: CIARmsLnR, CLFRmsLnR, PLFRmsLnR

   integer :: n_dtvrms_file, n_dtbrms_file
   character(len=72) :: dtvrms_file, dtbrms_file

   public :: dtVrms, dtBrms, initialize_RMS, zeroRms, finalize_RMS, get_nl_RMS, &
   &         transform_to_lm_RMS, compute_lm_forces, transform_to_grid_RMS

contains

   subroutine initialize_RMS
      !
      ! Memory allocation of arrays used in the computation of r.m.s. force balance
      !

      if ( l_mag_par_solve ) then
         allocate( dtBPol2hInt(lm_maxMag,nRstartMag:nRstopMag) )
         allocate( dtBTor2hInt(lm_maxMag,nRstartMag:nRstopMag) )
         allocate( dtBPolLMr(lm_maxMag,nRstartMag:nRstopMag) )
         bytes_allocated = bytes_allocated+2*lm_maxMag*(nRstopMag-nRstartMag+1)*&
         &                 SIZEOF_DEF_REAL+lm_maxMag*(nRstopMag-nRstartMag+1)*  &
         &                 SIZEOF_DEF_COMPLEX
      else
         allocate( dtBPol2hInt(llmMag:ulmMag,n_r_maxMag) )
         allocate( dtBTor2hInt(llmMag:ulmMag,n_r_maxMag) )
         allocate( dtBPolLMr(llmMag:ulmMag,n_r_maxMag) )
         bytes_allocated = bytes_allocated+2*(ulmMag-llmMag+1)*n_r_maxMag*   &
         &                 SIZEOF_DEF_REAL+(llmMag-ulmMag+1)*n_r_maxMag*     &
         &                 SIZEOF_DEF_COMPLEX
      end if

      allocate( DifPol2hInt(0:l_max,n_r_max), DifTor2hInt(0:l_max,n_r_max) )
      if ( l_parallel_solve ) then
         allocate( DifPolLMr(lm_max,nRstart:nRstop) )
      else
         allocate( DifPolLMr(llm:ulm,n_r_max) )
      end if
      bytes_allocated = bytes_allocated+                      &
      &                 2*(l_max+1)*n_r_max*SIZEOF_DEF_REAL+  &
      &                 (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

      allocate( Adv2hInt(0:l_max,n_r_max), Cor2hInt(0:l_max,n_r_max) )
      allocate( LF2hInt(0:l_max,n_r_max), Iner2hInt(0:l_max,n_r_max) )
      allocate( Pre2hInt(0:l_max,n_r_max), Geo2hInt(0:l_max,n_r_max) )
      allocate( Mag2hInt(0:l_max,n_r_max), Arc2hInt(0:l_max,n_r_max) )
      allocate( ArcMag2hInt(0:l_max,n_r_max), CIA2hInt(0:l_max,n_r_max) )
      allocate( CLF2hInt(0:l_max,n_r_max), PLF2hInt(0:l_max,n_r_max) )
      bytes_allocated = bytes_allocated+12*(l_max+1)*n_r_max*SIZEOF_DEF_REAL

      if ( l_heat ) then
         allocate( Buo_temp2hInt(0:l_max,n_r_max) )
         bytes_allocated = bytes_allocated+(l_max+1)*n_r_max*SIZEOF_DEF_REAL
      end if
      if ( l_chemical_conv ) then
         allocate( Buo_xi2hInt(0:l_max,n_r_max) )
         bytes_allocated = bytes_allocated+(l_max+1)*n_r_max*SIZEOF_DEF_REAL
      end if

      !-- RMS Calculations
      allocate( Advt2(nlat_padded,n_phi_max), Advp2(nlat_padded,n_phi_max) )
      allocate( dtVr(nlat_padded,n_phi_max), dtVt(nlat_padded,n_phi_max) )
      allocate( dtVp(nlat_padded,n_phi_max) )
      allocate( LFt2(nlat_padded,n_phi_max), LFp2(nlat_padded,n_phi_max) )
      allocate( CFt2(nlat_padded,n_phi_max), CFp2(nlat_padded,n_phi_max) )
      allocate( dpdtc(nlat_padded,n_phi_max), dpdpc(nlat_padded,n_phi_max) )
      bytes_allocated=bytes_allocated + 11*n_phi_max*nlat_padded*SIZEOF_DEF_REAL

      allocate( vt_old(nlat_padded,n_phi_max,nRstart:nRstop) )
      allocate( vp_old(nlat_padded,n_phi_max,nRstart:nRstop) )
      allocate( vr_old(nlat_padded,n_phi_max,nRstart:nRstop) )
      bytes_allocated=bytes_allocated + 3*n_phi_max*nlat_padded*(nRstop-nRstart+1)*&
      &               SIZEOF_DEF_REAL

      dtVr(:,:)=0.0_cp
      dtVt(:,:)=0.0_cp
      dtVp(:,:)=0.0_cp
      vt_old(:,:,:) =0.0_cp
      vr_old(:,:,:) =0.0_cp
      vp_old(:,:,:) =0.0_cp

      if ( l_adv_curl ) then
         allocate ( dpkindrc(nlat_padded,n_phi_max) )
         bytes_allocated=bytes_allocated + n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      end if

      allocate( dtVrLM(lmP_max), dtVtLM(lmP_max), dtVpLM(lmP_max), LFrLM(lmP_max) )
      allocate( Advt2LM(lmP_max), Advp2LM(lmP_max), LFt2LM(lmP_max), LFp2LM(lmP_max) )
      allocate( CFt2LM(lmP_max), CFp2LM(lmP_max), PFt2LM(lmP_max), PFp2LM(lmP_max) )
      bytes_allocated = bytes_allocated + 12*lmP_max*SIZEOF_DEF_COMPLEX
      if ( l_adv_curl ) then
         allocate( dpkindrLM(lmP_max) )
         bytes_allocated = bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
      end if

      call InerRmsL%initialize(0,l_max)
      call CorRmsL%initialize(0,l_max)
      call LFRmsL%initialize(0,l_max)
      call AdvRmsL%initialize(0,l_max)
      call DifRmsL%initialize(0,l_max)
      call Buo_tempRmsL%initialize(0,l_max)
      call Buo_xiRmsL%initialize(0,l_max)
      call PreRmsL%initialize(0,l_max)
      call GeoRmsL%initialize(0,l_max)
      call MagRmsL%initialize(0,l_max)
      call ArcRmsL%initialize(0,l_max)
      call ArcMagRmsL%initialize(0,l_max)
      call CIARmsL%initialize(0,l_max)
      call CLFRmsL%initialize(0,l_max)
      call PLFRmsL%initialize(0,l_max)

      if ( l_2D_RMS ) then
         call CorRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call AdvRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call LFRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call Buo_tempRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call Buo_xiRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call PreRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call DifRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call InerRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call GeoRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call MagRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call ArcRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call ArcMagRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call CIARmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call CLFRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call PLFRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
      end if

      if ( .not. l_finite_diff ) then
         allocate ( type_cheb_odd :: rscheme_RMS )
      else
         allocate ( type_fd :: rscheme_RMS )
      end if

      !--- Initialize new cut-back grid:
      call init_rNB(r,rCut,rDea,rC,n_r_maxC,n_cheb_maxC,nCut,rscheme_RMS)

      dtvrms_file='dtVrms.'//tag
      dtbrms_file='dtBrms.'//tag
      if ( rank == 0 .and. (.not. l_save_out) ) then
         open(newunit=n_dtvrms_file, file=dtvrms_file, status='new')
         if ( l_mag ) then
            open(newunit=n_dtbrms_file, file=dtbrms_file, status='new')
         end if
      end if

      call zeroRms() ! zero's the fields

   end subroutine initialize_RMS
!----------------------------------------------------------------------------
   subroutine finalize_RMS
      !
      ! Deallocate arrays used for r.m.s. force balance computation
      !

      deallocate( rC, dtBPol2hInt, dtBTor2hInt, dtBPolLMr )
      deallocate( DifPol2hInt, DifTor2hInt, DifPolLMr )
      deallocate( Adv2hInt, Cor2hInt, LF2hInt, Iner2hInt, Pre2hInt )
      deallocate( Geo2hInt, Mag2hInt, Arc2hInt, CIA2hInt )
      deallocate( CLF2hInt, PLF2hInt, ArcMag2hInt )

      if ( l_chemical_conv )  deallocate( Buo_xi2hInt )
      if ( l_heat )  deallocate( Buo_temp2hInt )

      deallocate ( Advt2, Advp2, LFt2, LFp2, CFt2, CFp2, dpdtc, dpdpc )
      deallocate ( dtVr, dtVt, dtVp, vr_old, vt_old, vp_old )
      deallocate( Advt2LM, Advp2LM, LFt2LM, LFp2LM, CFt2LM, CFp2LM, PFt2LM, PFp2LM )
      deallocate( dtVrLM, dtVtLM, dtVpLM, LFrLM )
      if ( l_adv_curl ) deallocate ( dpkindrLM, dpkindrc )

      call InerRmsL%finalize()
      call CorRmsL%finalize()
      call LFRmsL%finalize()
      call AdvRmsL%finalize()
      call DifRmsL%finalize()
      call Buo_tempRmsL%finalize()
      call Buo_xiRmsL%finalize()
      call PreRmsL%finalize()
      call GeoRmsL%finalize()
      call MagRmsL%finalize()
      call ArcMagRmsL%finalize()
      call ArcRmsL%finalize()
      call CIARmsL%finalize()
      call CLFRmsL%finalize()
      call PLFRmsL%finalize()

      if ( l_2D_RMS ) then
          call CorRmsLnR%finalize()
          call AdvRmsLnR%finalize()
          call LFRmsLnR%finalize()
          call Buo_tempRmsLnR%finalize()
          call Buo_xiRmsLnR%finalize()
          call PreRmsLnR%finalize()
          call DifRmsLnR%finalize()
          call InerRmsLnR%finalize()
          call GeoRmsLnR%finalize()
          call MagRmsLnR%finalize()
          call ArcMagRmsLnR%finalize()
          call ArcRmsLnR%finalize()
          call CIARmsLnR%finalize()
          call CLFRmsLnR%finalize()
          call PLFRmsLnR%finalize()
      end if

      call rscheme_RMS%finalize()

      if ( rank == 0 .and. (.not. l_save_out) ) then
         close(n_dtvrms_file)
         if ( l_mag ) close(n_dtbrms_file)
      end if

   end subroutine finalize_RMS
!----------------------------------------------------------------------------
   subroutine zeroRms
      !
      !  Zeros integrals that are set in get_td, update_z,
      !  update_wp, update_b, dtVrms and dtBrms
      !

      DifPol2hInt(:,:)=0.0_cp
      DifTor2hInt(:,:)=0.0_cp
      dtBPol2hInt(:,:)=0.0_cp
      dtBTor2hInt(:,:)=0.0_cp

      Adv2hInt(:,:)   =0.0_cp
      Iner2hInt(:,:)  =0.0_cp
      Cor2hInt(:,:)   =0.0_cp
      LF2hInt(:,:)    =0.0_cp
      Pre2hInt(:,:)   =0.0_cp
      Geo2hInt(:,:)   =0.0_cp
      Mag2hInt(:,:)   =0.0_cp
      Arc2hInt(:,:)   =0.0_cp
      ArcMag2hInt(:,:)=0.0_cp
      CIA2hInt(:,:)   =0.0_cp
      CLF2hInt(:,:)   =0.0_cp
      PLF2hInt(:,:)   =0.0_cp
      if ( l_chemical_conv ) Buo_xi2hInt(:,:)=0.0_cp
      if ( l_heat ) Buo_temp2hInt(:,:)=0.0_cp

      Advt2LM(:)=zero
      Advp2LM(:)=zero
      LFp2LM(:) =zero
      LFt2LM(:) =zero
      CFt2LM(:) =zero
      CFp2LM(:) =zero
      PFt2LM(:) =zero
      PFp2LM(:) =zero
      dtVtLM(:) =zero
      dtVpLM(:) =zero
      dtVrLM(:) =zero
      if ( l_adv_curl ) dpkindrLM(:)=zero

      DifPolLMr(:,:)=zero
      dtBPolLMr(:,:)=zero

   end subroutine zeroRms
!----------------------------------------------------------------------------
   subroutine init_rNB(r,rCut,rDea,r2,n_r_max2,n_cheb_max2,nS,rscheme_RMS)
      !
      ! Prepares the usage of a cut back radial grid where nS points
      ! on both boundaries are discarded.
      ! The aim actually is to discard boundary effects, but just
      ! not considering the boundary grid points does not work when
      ! you also want radial derivatives and integrals. For these
      ! we use the Chebychev transform which needs are particular number
      ! of grid points so that the fast cosine transform can be
      ! applied. Therefor more than just 2 points have to be
      ! thrown away, which may make sense anyway.
      !

      !--- Input variables:
      real(cp),              intent(in) :: r(n_r_max)
      real(cp),              intent(in) :: rCut,rDea

      !--- Output variables:
      integer,               intent(out) :: nS,n_r_max2,n_cheb_max2
      real(cp), allocatable, intent(out) :: r2(:)
      class(type_rscheme),   intent(out) :: rscheme_RMS

      ! Local stuff
      real(cp) :: r2C(n_r_max)
      real(cp) :: dr2(n_r_max)
      real(cp) :: r_icb2, r_cmb2
      real(cp) :: ratio1, ratio2
      integer :: nRs(63), n_in_2

      logical :: lStop
      integer :: n,nR

      !--- New radial grid:

      !--- Find number of points to be cut away at either side:
      lStop=.true.
      do nS=1,(n_r_max-1)/2
         if ( r(1)-r(nS) >= rCut ) then
            lStop=.false.
            exit
         end if
      end do
      if ( lStop ) then
         call abortRun('No nS found in init_rNB!')
      end if
      nS=nS-1
      n_r_max2=n_r_max-2*nS

      if ( .not. l_finite_diff ) then
         ! Allowed number of radial grid points:
         nRs = [25, 33, 37, 41, 49, 61, 65, 73, 81, 97, 101, 109, 121,  &
         &      129, 145, 161, 181, 193, 201, 217, 241, 257, 289, 301,  &
         &      321, 325, 361, 385, 401, 433, 481, 501, 513, 541, 577,  &
         &      601, 641, 649, 721, 769, 801, 865, 901, 961, 973, 1001, &
         &      1025, 1081, 1153, 1201, 1281, 1297, 1441, 1501, 1537,   &
         &      1601, 1621, 1729, 1801, 1921, 1945, 2001, 2049]
         lStop=.true.
         do n=size(nRs),1,-1
            if ( nRs(n) <= n_r_max2 ) then
               lStop=.false.
               exit
            end if
         end do
         if ( lStop ) then
            call abortRun('No n_r_max2 found in init_rNB!')
         end if

         n_r_max2=nRs(n)
         nS=(n_r_max-n_r_max2)/2
         n_cheb_max2=min(int((one-rDea)*n_r_max2),n_cheb_max)

         allocate( r2(n_r_max2) )
         bytes_allocated = bytes_allocated+n_r_max2*SIZEOF_DEF_REAL

         do nR=1,n_r_max2
            r2(nR)=r(nR+nS)
         end do
         r_icb2=r2(n_r_max2)
         r_cmb2=r2(1)

         if ( l_newmap ) then
            n_in_2 = 1
            ratio1 = alph1
            ratio2 = alph2
         else
            n_in_2 = 0
            ratio1 = 0.0_cp
            ratio2 = 0.0_cp
         end if
         call rscheme_RMS%initialize(n_r_max2, n_cheb_max2, n_in_2)
         call rscheme_RMS%get_grid(n_r_max2, r_icb2, r_cmb2, ratio1, ratio2, r2C)

         if ( n_r_max2 == n_r_max ) then
            rscheme_RMS%drx(:)=rscheme_oc%drx(:)
         else
            do nR=1,n_r_max2
               rscheme_RMS%drx(nR)=one
            end do
            call get_dr(r2,dr2,n_r_max2,rscheme_RMS)
            do nR=1,n_r_max2
               rscheme_RMS%drx(nR)=one/dr2(nR)
            end do
         end if

      else ! finite differences

         allocate( r2(n_r_max2) )
         bytes_allocated = bytes_allocated+n_r_max2*SIZEOF_DEF_REAL

         do nR=1,n_r_max2
            r2(nR)=r(nR+nS)
         end do
         r_icb2=r2(n_r_max2)
         r_cmb2=r2(1)

         call rscheme_RMS%initialize(n_r_max, rscheme_oc%order, &
              &                      rscheme_oc%order_boundary)
         ratio1 = fd_stretch
         ratio2 = fd_ratio
         call rscheme_RMS%get_grid(n_r_max, r_icb, r_cmb, ratio1, ratio2, r2C)
         call rscheme_oc%get_der_mat(n_r_max)

      end if

   end subroutine init_rNB
!----------------------------------------------------------------------------
   subroutine get_nl_RMS(nR,vr,vt,vp,dvrdr,dvrdt,dvrdp,dvtdr,dvtdp,dvpdr,dvpdp, &
              &          cvr,Advt,Advp,LFt,LFp,tscheme,lRmsCalc)
      !
      ! This subroutine computes the r.m.s. force balance terms which need to
      ! be computed on the grid
      !

      !-- Input variables
      real(cp),            intent(in) :: vr(:,:), vt(:,:), vp(:,:), cvr(:,:)
      real(cp),            intent(in) :: dvrdr(:,:), dvrdt(:,:), dvrdp(:,:)
      real(cp),            intent(in) :: dvtdp(:,:), dvpdp(:,:)
      real(cp),            intent(in) :: dvtdr(:,:), dvpdr(:,:)
      real(cp),            intent(in) :: Advt(:,:),Advp(:,:),LFt(:,:),LFp(:,:)
      class(type_tscheme), intent(in) :: tscheme ! time scheme
      integer,             intent(in) :: nR      ! radial level
      logical,             intent(in) :: lRmsCalc

      !-- Local variables
      real(cp) ::  O_dt
      integer :: nPhi, nPhStart, nPhStop

      !$omp parallel default(shared) private(nPhStart,nPhStop,nPhi)
      nPhStart=1; nPhStop=n_phi_max
      call get_openmp_blocks(nPhStart,nPhStop)

      do nPhi=nPhStart,nPhStop

         if ( lRmsCalc ) then
            dpdtc(:,nPhi)=dpdtc(:,nPhi)*or1(nR)
            dpdpc(:,nPhi)=dpdpc(:,nPhi)*or1(nR)
            CFt2(:,nPhi)=-two*CorFac*cosTheta(:)*vp(:,nPhi)*or1(nR)
            CFp2(:,nPhi)= two*CorFac*sinTheta(:)* (or1(nR)*cosTheta(:)*&
            &                  O_sin_theta(:)*vt(:,nPhi)+or2(nR)*sinTheta(:)*vr(:,nPhi) )
            if ( l_conv_nl ) then
               Advt2(:,nPhi)=r(nR)*Advt(:,nPhi)
               Advp2(:,nPhi)=r(nR)*Advp(:,nPhi)
            end if
            if ( l_mag_LF .and. nR > n_r_LCR ) then
               LFt2(:,nPhi)=r(nR)*LFt(:,nPhi)
               LFp2(:,nPhi)=r(nR)*LFp(:,nPhi)
            end if

            if ( l_adv_curl ) then
               dpdtc(:,nPhi)=dpdtc(:,nPhi)-or3(nR)*( or2(nR)*  &
               &             vr(:,nPhi)*dvrdt(:,nPhi) -        &
               &             vt(:,nPhi)*(dvrdr(:,nPhi)+        &
               &             dvpdp(:,nPhi)+cosn_theta_E2(:) *  &
               &             vt(:,nPhi))+ vp(:,nPhi)*(         &
               &             cvr(:,nPhi)+dvtdp(:,nPhi)-        &
               &             cosn_theta_E2(:)*vp(:,nPhi)) )
               dpdpc(:,nPhi)=dpdpc(:,nPhi)- or3(nR)*( or2(nR)*  &
               &             vr(:,nPhi)*dvrdp(:,nPhi) +         &
               &             vt(:,nPhi)*dvtdp(:,nPhi) +         &
               &             vp(:,nPhi)*dvpdp(:,nPhi) )
               if ( l_conv_nl ) then
                  Advt2(:,nPhi)=Advt2(:,nPhi)-or3(nR)*( or2(nR)*  &
                  &             vr(:,nPhi)*dvrdt(:,nPhi) -        &
                  &             vt(:,nPhi)*(dvrdr(:,nPhi)+        &
                  &             dvpdp(:,nPhi)+cosn_theta_E2(:) *  &
                  &             vt(:,nPhi))+vp(:,nPhi)*(          &
                  &             cvr(:,nPhi)+dvtdp(:,nPhi)-        &
                  &             cosn_theta_E2(:)*vp(:,nPhi)) )
                  Advp2(:,nPhi)=Advp2(:,nPhi)-or3(nR)*( or2(nR)* &
                  &             vr(:,nPhi)*dvrdp(:,nPhi) +       &
                  &             vt(:,nPhi)*dvtdp(:,nPhi) +       &
                  &             vp(:,nPhi)*dvpdp(:,nPhi) )
               end if

               !- dpkin/dr = 1/2 d (u^2) / dr = ur*dur/dr+ut*dut/dr+up*dup/dr
               dpkindrc(:,nPhi)=or4(nR)*vr(:,nPhi)*(dvrdr(:,nPhi)-           &
               &                two*or1(nR)*vr(:,nPhi))+or2(nR)*             &
               &                O_sin_theta_E2(:)*(         vt(:,nPhi)*(     &
               &                      dvtdr(:,nPhi)-or1(nR)*vt(:,nPhi) ) +   &
               &                vp(:,nPhi)*(dvpdr(:,nPhi)-or1(nR)*vp(:,nPhi) ) )
            end if
         end if

         if ( tscheme%istage == 1 ) then
            O_dt = 1.0_cp/tscheme%dt(1)
            dtVr(:,nPhi)=O_dt*or2(nR)*(vr(:,nPhi)-vr_old(:,nPhi,nR))
            dtVt(:,nPhi)=O_dt*or1(nR)*(vt(:,nPhi)-vt_old(:,nPhi,nR))
            dtVp(:,nPhi)=O_dt*or1(nR)*(vp(:,nPhi)-vp_old(:,nPhi,nR))

            vr_old(:,nPhi,nR)=vr(:,nPhi)
            vt_old(:,nPhi,nR)=vt(:,nPhi)
            vp_old(:,nPhi,nR)=vp(:,nPhi)
         end if

      end do
      !$omp end parallel

   end subroutine get_nl_RMS
!----------------------------------------------------------------------------
   subroutine transform_to_grid_RMS(nR, p_Rloc)
      !
      ! This subroutine is used to transform the arrays used in r.m.s. force
      ! calculations from the spectral space to the physical grid.
      !

      !-- Input variables
      integer,     intent(in) :: nR ! radial level
      complex(cp), intent(inout) :: p_Rloc(lm_max,nRstart:nRstop) ! pressure in LM space

      call scal_to_grad_spat(p_Rloc(:,nR), dpdtc, dpdpc, l_R(nR))

   end subroutine transform_to_grid_RMS
!----------------------------------------------------------------------------
   subroutine transform_to_lm_RMS(nR, LFr)
      !
      ! This subroutine is used to transform the arrays used in r.m.s. force
      ! calculations from the physical grid to spectral space.
      !

      !-- Input variables
      integer,  intent(in) :: nR ! radial level
      real(cp), intent(inout) :: LFr(:,:) ! radial component of the Lorentz force

      Advt2LM(:)=zero
      Advp2LM(:)=zero
      LFrLM(:)  =zero
      LFp2LM(:) =zero
      LFt2LM(:) =zero
      CFt2LM(:) =zero
      CFp2LM(:) =zero
      PFt2LM(:) =zero
      PFp2LM(:) =zero
      dtVtLM(:) =zero
      dtVpLM(:) =zero
      dtVrLM(:) =zero
      if ( l_adv_curl ) dpkindrLM(:)=zero

      if ( l_mag_LF .and. nR>n_r_LCR ) call scal_to_SH(LFr, LFrLM, l_R(nR))
      call spat_to_sphertor(dpdtc, dpdpc, PFt2LM, PFp2LM, l_R(nR))
      call spat_to_sphertor(CFt2, CFp2, CFt2LM, CFp2LM, l_R(nR))
      call spat_to_qst(dtVr, dtVt, dtVp, dtVrLM, dtVtLM, dtVpLM, l_R(nR))
      if ( l_conv_nl ) call spat_to_sphertor(Advt2, Advp2, Advt2LM, Advp2LM, l_R(nR))
      !-- Kinetic pressure : 1/2 d u^2 / dr
      if ( l_adv_curl ) call scal_to_SH(dpkindrc, dpkindrLM, l_R(nR))
      if ( l_mag_nl .and. nR>n_r_LCR ) call spat_to_sphertor(LFt2, LFp2,  &
                                            &                LFt2LM, LFp2LM, l_R(nR))

   end subroutine transform_to_lm_RMS
!----------------------------------------------------------------------------
   subroutine compute_lm_forces(nR, w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, s_Rloc, &
              &                 xi_Rloc, p_Rloc, dp_Rloc, AdvrLM)
      !
      ! This subroutine finalizes the computation of the r.m.s. spectra once
      ! the quantities are back in spectral space.
      !

      !-- Input variables
      integer,     intent(in) :: nR
      complex(cp), intent(in) :: w_Rloc(:), dw_Rloc(:), z_Rloc(:)
      complex(cp), intent(in) :: s_Rloc(:), p_Rloc(:), dp_Rloc(:)
      complex(cp), intent(in) :: xi_Rloc(:), ddw_Rloc(:)
      complex(cp), intent(in) :: AdvrLM(:)

      !-- Local variables
      complex(cp) :: dpdr(lm_max), Buo_temp(lm_max), Buo_xi(lm_max)
      complex(cp) :: LFPol(lm_max), AdvPol(lm_max), CorPol(lm_max)
      complex(cp) :: Geo(lm_max),CLF(lm_max),PLF(lm_max)
      complex(cp) :: ArcMag(lm_max),Mag(lm_max),CIA(lm_max),Arc(lm_max)
      complex(cp) :: CorPol_loc, AdvPol_loc
      integer :: lm, lmA, lmP, lmS, l, m

      !-- l=m=0 spherically-symmetric contributions
      lm=1
      lmA=lm2lmA(lm)
      lmP=lm2lmP(lm)

      if ( l_conv_nl ) then
         AdvPol(lm)=or2(nR)*AdvrLM(lmP)
         if ( l_adv_curl ) AdvPol(lm)=AdvPol(lm)-dpkindrLM(lmP)
      else
         AdvPol(lm)=zero
      end if

      if ( l_corr ) then
         CorPol(lm)=two*CorFac*or1(nR)*dTheta2A(lm)*z_Rloc(lmA)
      else
         CorPol(lm)=zero
      end if

      if (l_heat) then
         Buo_temp(lm)=BuoFac*rgrav(nR)*rho0(nR)*s_Rloc(lm)
      else
         Buo_temp(lm)=0.0_cp
      end if

      if (l_chemical_conv) then
         Buo_xi(lm)=ChemFac*rgrav(nR)*rho0(nR)*xi_Rloc(lm)
      else
         Buo_xi(lm)=0.0_cp
      end if

      if ( l_mag_LF .and. nR>n_r_LCR ) then
         LFPol(lm) =or2(nR)*LFrLM(lmP)
         AdvPol(lm)=AdvPol(lm)-LFPol(lm)
      end if

      if ( l_double_curl ) then
         !-- Recalculate the pressure gradient based on the poloidal
         !-- equation equilibrium
         dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)+beta(nR)*p_Rloc(lm)+AdvPol(lm)+CorPol(lm)
      else
         dpdr(lm)=dp_Rloc(lm)
      end if

      !-- Loop over the other (l,m) modes
      !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
      !$omp private(AdvPol_loc,CorPol_loc)
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmS =lm2lmS(lm)
         lmA =lm2lmA(lm)
         lmP =lm2lmP(lm)

         if ( l_anelastic_liquid ) then
            if (l_heat) then
               Buo_temp(lm) =BuoFac*alpha0(nR)*rgrav(nR)*(      &
               &       rho0(nR)*s_Rloc(lm)-ViscHeatFac*         &
               &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
               &     p_Rloc(lm) )
            else
               Buo_temp(lm)=zero
            end if

            if (l_chemical_conv) then
               Buo_xi(lm)=ChemFac*alpha0(nR)*rgrav(nR)*rho0(nR)*xi_Rloc(lm)
            else
               Buo_xi(lm)=zero
            end if
         else
            if (l_heat) then
               Buo_temp(lm)=BuoFac*rho0(nR)*rgrav(nR)*s_Rloc(lm)
            else
               Buo_temp(lm)=zero
            end if
            if (l_chemical_conv) then
               Buo_xi(lm) =ChemFac*rho0(nR)*rgrav(nR)*xi_Rloc(lm)
            else
               Buo_xi(lm)=zero
            end if
         end if

         !-- We need to compute the Coriolis term once again
         if ( l_corr ) then
            if ( l < l_max .and. l > m ) then
               CorPol_loc =two*CorFac*or1(nR) * (  &
               &           dPhi(lm)*dw_Rloc(lm) +  & ! phi-deriv of dw/dr
               &       dTheta2A(lm)*z_Rloc(lmA) -  & ! sin(theta) dtheta z
               &       dTheta2S(lm)*z_Rloc(lmS) )
            else if ( l == l_max ) then
               CorPol_loc=zero
            else if ( l == m ) then
               CorPol_loc = two*CorFac*or1(nR) * (  &
               &            dPhi(lm)*dw_Rloc(lm)  + &
               &        dTheta2A(lm)*z_Rloc(lmA) )
            end if
         else
            CorPol_loc=zero
         end if

         ! We also need to recompute AdvPol_loc here
         if ( l_conv_nl ) then
            AdvPol_loc=or2(nR)*AdvrLM(lmP)
            if ( l_adv_curl ) AdvPol_loc=AdvPol_loc-dpkindrLM(lmP)
         else
            AdvPol_loc=zero
         end if

         if ( l_double_curl ) then
            !-- Recalculate the pressure gradient based on the poloidal
            !-- equation equilibrium
            dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)-dtVrLM(lmP)+              &
            &        dLh(lm)*or2(nR)*hdif_V(l)*visc(nR)*(              &
            &                                        ddw_Rloc(lm)+     &
            &         (two*dLvisc(nR)-third*beta(nR))*dw_Rloc(lm)-     &
            &         ( dLh(lm)*or2(nR)+four*third*( dbeta(nR)+        &
            &         dLvisc(nR)*beta(nR)+(three*dLvisc(nR)+beta(nR))* &
            &         or1(nR)) )*                      w_Rloc(lm))+    &
            &        beta(nR)*p_Rloc(lm)+AdvPol_loc+CorPol_loc
         else
            dpdr(lm)=dp_Rloc(lm)
         end if

         if ( l_mag_LF .and. nR>n_r_LCR ) then
            LFPol(lm) =or2(nR)*LFrLM(lmP)
            AdvPol(lm)=AdvPol_loc-LFPol(lm)
         else
            AdvPol(lm)=AdvPol_loc
         end if
         CorPol(lm)=CorPol_loc
      end do

      !-- Now compute R.M.S spectra
      if ( l_conv_nl ) then
         call hIntRms(AdvPol,nR,1,lm_max,0,Adv2hInt(:,nR),st_map, .false.)
         call hIntRms(Advt2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map,.true.)
         call hIntRms(Advp2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map,.true.)
         do lm=1,lm_max
            lmP=lm2lmP(lm)
            !-- Use Geo as work array
            Geo(lm)=AdvPol(lm)-dtVrLM(lmP)
         end do
         call hIntRms(Geo,nR,1,lm_max,0,Iner2hInt(:,nR),st_map, .false.)
         call hIntRms(Advt2LM-dtVtLM,nR,1,lmP_max,1,Iner2hInt(:,nR),st_map,.true.)
         call hIntRms(Advp2LM-dtVpLM,nR,1,lmP_max,1,Iner2hInt(:,nR),st_map,.true.)
      end if

      if ( l_anelastic_liquid ) then
         call hIntRms(dpdr,nR,1,lm_max,0,Pre2hInt(:,nR),st_map,.false.)
      else
         ! rho* grad(p/rho) = grad(p) - beta*p
         !-- Geo is used to store the pressure Gradient
         if ( l_adv_curl ) then
            do lm=1,lm_max
               lmP=lm2lmP(lm)
               !-- Use Geo as work array
               Geo(lm)=dpdr(lm)-beta(nR)*p_Rloc(lm)-dpkindrLM(lmP)
            end do
         else
            do lm=1,lm_max
               !-- Use Geo as work array
               Geo(lm)=dpdr(lm)-beta(nR)*p_Rloc(lm)
            end do
         end if
         call hIntRms(Geo,nR,1,lm_max,0,Pre2hInt(:,nR),st_map,.false.)
      end if
      call hIntRms(PFt2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),st_map,.true.)
      call hIntRms(PFp2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),st_map,.true.)

      if ( l_heat ) then
         call hIntRms(Buo_temp,nR,1,lm_max,0,Buo_temp2hInt(:,nR),st_map,.false.)
      end if
      if ( l_chemical_conv ) then
         call hIntRms(Buo_xi,nR,1,lm_max,0,Buo_xi2hInt(:,nR),st_map,.false.)
      end if
      if ( l_corr ) then
         call hIntRms(CorPol,nR,1,lm_max,0,Cor2hInt(:,nR),st_map,.false.)
         call hIntRms(CFt2LM,nR,1,lmP_max,1,Cor2hInt(:,nR),st_map,.true.)
         call hIntRms(CFp2LM,nR,1,lmP_max,1,Cor2hInt(:,nR),st_map,.true.)
      end if
      if ( l_mag_LF .and. nR>n_r_LCR ) then
         call hIntRms(LFPol,nR,1,lm_max,0,LF2hInt(:,nR),st_map,.false.)
         call hIntRms(LFt2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map,.true.)
         call hIntRms(LFp2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map,.true.)
      end if

      do lm=1,lm_max
         lmP=lm2lmP(lm)
         Geo(lm)=CorPol(lm)-dpdr(lm)+beta(nR)*p_Rloc(lm)
         PLF(lm)=LFPol(lm)-dpdr(lm)+beta(nR)*p_Rloc(lm)
         if ( l_adv_curl ) then
            Geo(lm)=Geo(lm)+dpkindrLM(lmP)
            PLF(lm)=PLF(lm)+dpkindrLM(lmP)
         end if
         CLF(lm)=CorPol(lm)+LFPol(lm)
         Mag(lm)=Geo(lm)+LFPol(lm)
         Arc(lm)=Geo(lm)+Buo_temp(lm)+Buo_xi(lm)
         ArcMag(lm)=Mag(lm)+Buo_temp(lm)+Buo_xi(lm)
         CIA(lm)=ArcMag(lm)+AdvPol(lm)-dtVrLM(lmP)
         !CIA(lm)=CorPol(lm)+Buo_temp(lm)+Buo_xi(lm)+AdvPol(lm)
      end do
      call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.false.)
      call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.false.)
      call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.false.)
      call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.false.)
      call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.false.)
      call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.false.)
      call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.false.)

      do lm=1,lm_max
         lmP=lm2lmP(lm)
         Geo(lm)=-CFt2LM(lmP)-PFt2LM(lmP)
         CLF(lm)=-CFt2LM(lmP)+LFt2LM(lmP)
         PLF(lm)=LFt2LM(lmP)-PFt2LM(lmP)
         Mag(lm)=Geo(lm)+LFt2LM(lmP)
         Arc(lm)=Geo(lm)
         ArcMag(lm)=Mag(lm)
         CIA(lm)=ArcMag(lm)+Advt2LM(lmP)-dtVtLM(lmP)
         !CIA(lm)=-CFt2LM(lmP)+Advt2LM(lmP)
      end do
      call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.true.)
      call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.true.)
      call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.true.)
      call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.true.)
      call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.true.)
      call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.true.)
      call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.true.)

      do lm=1,lm_max
         lmP=lm2lmP(lm)
         Geo(lm)=-CFp2LM(lmP)-PFp2LM(lmP)
         CLF(lm)=-CFp2LM(lmP)+LFp2LM(lmP)
         PLF(lm)=LFp2LM(lmP)-PFp2LM(lmP)
         Mag(lm)=Geo(lm)+LFp2LM(lmP)
         Arc(lm)=Geo(lm)
         ArcMag(lm)=Mag(lm)
         CIA(lm)=ArcMag(lm)+Advp2LM(lmP)-dtVpLM(lmP)
         !CIA(lm)=-CFp2LM(lmP)+Advp2LM(lmP)
      end do
      call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.true.)
      call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.true.)
      call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.true.)
      call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.true.)
      call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.true.)
      call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.true.)
      call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.true.)

   end subroutine compute_lm_forces
!----------------------------------------------------------------------------
   subroutine get_force(Force2hInt,ForceRms,ForceRmsL,ForceRmsLnR,      &
              &         volC,nRMS_sets,timePassed,timeNorm,l_stop_time, &
              &         ForcePol2hInt,ForceTor2hInt)
      !
      ! This subroutine is used to compute the contributions of the
      ! forces in the Navier-Stokes equation
      !

      !-- Input variables
      real(cp), intent(in) :: Force2hInt(0:l_max,n_r_max)
      real(cp), optional, intent(in) :: ForcePol2hInt(0:l_max,n_r_max)
      real(cp), optional, intent(in) :: ForceTor2hInt(0:l_max,n_r_max)
      real(cp), intent(in) :: timePassed
      real(cp), intent(in) :: timeNorm
      real(cp), intent(in) :: volC
      integer,  intent(in) :: nRMS_sets
      logical,  intent(in) :: l_stop_time

      !-- Output variables
      real(cp), intent(out) :: ForceRms
      type(mean_sd_type), intent(inout) :: ForceRmsL
      type(mean_sd_2D_type), intent(inout) :: ForceRmsLnR

      !-- Local variables
      real(cp) :: ForceRms_L,ForceRms_LnR(n_r_max)
      real(cp) :: Rms(n_r_max), tmp(0:l_max)
      integer :: l,nR,nRC

      nRC=nCut+1

      ForceRms=0.0_cp
      do l=0,l_max
         if ( present(ForcePol2hInt) .and. present(ForceTor2hInt) ) then
            do nR=1,n_r_max
               Rms(nR)=0.0_cp
               Rms(nR)=Rms(nR)+ForcePol2hInt(l,nR)+ForceTor2hInt(l,nR)
            end do
         else
            do nR=1,n_r_max
               Rms(nR)=Force2hInt(l,nR)
            end do
         end if

         if ( l_2D_RMS ) then
            ForceRms_LnR=sqrt(Rms/vol_oc)
            call ForceRmsLnR%compute(ForceRms_LnR, nRMS_sets,   &
                 &                   timePassed, timeNorm, l)
         end if

         ForceRms_L=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
         ForceRms  =ForceRms+ForceRms_L
         tmp(l)    =sqrt(ForceRms_L/volC)
      end do
      call ForceRmsL%compute(tmp, nRMS_sets, timePassed, timeNorm)

      ForceRms=sqrt(ForceRms/volC)

      if ( l_stop_time ) then
         call ForceRmsL%finalize_SD(timeNorm)
      end if

   end subroutine get_force
!-----------------------------------------------------------------------------
   subroutine dtVrms(time,nRMS_sets,timePassed,timeNorm,l_stop_time)
      !
      ! This routine calculates and stores the different contributions
      ! of the forces entering the Navier-Stokes equation.
      !

      !-- Input variables:
      real(cp), intent(in) :: time
      real(cp), intent(in) :: timePassed
      real(cp), intent(in) :: timeNorm
      logical,  intent(in) :: l_stop_time
      integer,  intent(inout) :: nRMS_sets

      !-- Output:
      real(cp) :: InerRms    =0.0_cp
      real(cp) :: CorRms     =0.0_cp
      real(cp) :: AdvRms     =0.0_cp
      real(cp) :: LFRms      =0.0_cp
      real(cp) :: DifRms     =0.0_cp
      real(cp) :: Buo_tempRms=0.0_cp
      real(cp) :: Buo_xiRms  =0.0_cp
      real(cp) :: PreRms     =0.0_cp
      real(cp) :: GeoRms     =0.0_cp
      real(cp) :: MagRms     =0.0_cp
      real(cp) :: ArcRms     =0.0_cp
      real(cp) :: ArcMagRms  =0.0_cp
      real(cp) :: CIARms     =0.0_cp
      real(cp) :: CLFRms     =0.0_cp
      real(cp) :: PLFRms     =0.0_cp

      !-- Local variables:
      integer :: nR,nRC,l,fileHandle,version
      real(cp) :: volC
      real(cp) :: Dif2hInt(n_r_max)

      complex(cp) :: workA(llm:ulm,n_r_max), work_Rloc(lm_max,nRstart:nRstop)
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(cp) :: global_sum(l_max+1,n_r_max)
      integer :: irank,sendcount
      character(len=80) :: fileName

      nRC=nCut+1

      !-- Diffusion
      DifRms=0.0_cp
      if ( rscheme_RMS%version == 'cheb' ) then
         call get_dr(DifPolLMr(llm:ulm,:),workA(llm:ulm,:),ulm-llm+1,1, &
              &      ulm-llm+1,n_r_max,rscheme_oc,nocopy=.true.)
      else
         if ( l_parallel_solve ) then
            call get_dr_Rloc(DifPolLMr,work_Rloc,lm_max,nRstart,nRstop,n_r_max,&
                 &           rscheme_oc)
         else
            call get_dr(DifPolLMr(llm:ulm,:),workA(llm:ulm,:),ulm-llm+1,1, &
                 &      ulm-llm+1,n_r_max,rscheme_oc)
         end if
      end if

      if ( l_parallel_solve ) then
         do nR=nRstart,nRstop
            call hInt2dPol(work_Rloc(:,nR),1,lm_max,DifPol2hInt(:,nR),st_map)
         end do
      else
         do nR=1,n_r_max
            call hInt2dPol(workA(llm:ulm,nR),llm,ulm,DifPol2hInt(:,nR),lo_map)
         end do
      end if

#ifdef WITH_MPI
      ! The following fields are only 1D and R distributed.
      sendcount = nR_per_rank*(l_max+1)
      do irank=0,n_procs-1
         recvcounts(irank) = radial_balance(irank)%n_per_rank*(l_max+1)
      end do
      displs(0)=0
      do irank=1,n_procs-1
         displs(irank) = displs(irank-1)+recvcounts(irank-1)
      end do

      if ( l_parallel_solve ) then
         call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,         &
              &              DifPol2hInt,recvcounts,displs,MPI_DEF_REAL,  &
              &              MPI_COMM_WORLD,ierr)
      else
         call MPI_Reduce(DifPol2hInt(:,:),global_sum,n_r_max*(l_max+1), &
              &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         if ( rank == 0 ) DifPol2hInt(:,:)=global_sum
      end if
#endif

      ! First gather all needed arrays on rank 0
      ! some more arrays to gather for the dtVrms routine
      ! we need some more fields for the dtBrms routine
#ifdef WITH_MPI
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Cor2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Adv2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Iner2hInt,recvcounts,displs,MPI_DEF_REAL,  &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              LF2hInt,recvcounts,displs,MPI_DEF_REAL,    &
           &              MPI_COMM_WORLD,ierr)

      if ( l_heat ) then
         call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
              &              Buo_temp2hInt,recvcounts,displs,           &
              &              MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      end if

      if ( l_chemical_conv ) then
         call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
              &              Buo_xi2hInt,recvcounts,displs,MPI_DEF_REAL,&
              &              MPI_COMM_WORLD,ierr)
      end if

      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Pre2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Geo2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Mag2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Arc2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              ArcMag2hInt,recvcounts,displs,MPI_DEF_REAL,&
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              CIA2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              CLF2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              PLF2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      if ( l_parallel_solve ) then
         call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,         &
              &              DifTor2hInt,recvcounts,displs,MPI_DEF_REAL,  &
              &              MPI_COMM_WORLD,ierr)
      else
         call MPI_Reduce(DifTor2hInt(:,:),global_sum,n_r_max*(l_max+1), &
              &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         if ( rank == 0 ) DifTor2hInt(:,:)=global_sum
      end if
#endif

      if ( rank == 0 ) then

         nRMS_sets=nRMS_sets+1
         volC=four*third*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)

         !-- Coriolis force
         if ( l_corr ) then
            call get_force(Cor2hInt,CorRms,CorRmsL,CorRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Advection
         if ( l_conv_nl ) then
            call get_force(Adv2hInt,AdvRms,AdvRmsL,AdvRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Lorentz force
         if ( l_mag_LF ) then
            call get_force(LF2hInt,LFRms,LFRmsL,LFRmsLnR,volC,        &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Thermal Buoyancy
         if ( l_heat ) then
            call get_force(Buo_temp2hInt,Buo_tempRms,Buo_tempRmsL,   &
                 &         Buo_tempRmsLnR,volC,nRMS_sets,timePassed, &
                 &         timeNorm,l_stop_time)
         end if

         !-- Chemical Buoyancy
         if ( l_chemical_conv ) then
            call get_force(Buo_xi2hInt,Buo_xiRms,Buo_xiRmsL,Buo_xiRmsLnR, &
                 &         volC,nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Pressure gradient
         call get_force(Pre2hInt,PreRms,PreRmsL,PreRmsLnR,volC,    &
              &         nRMS_sets,timePassed,timeNorm,l_stop_time)

         !-- Geostrophic balance
         if ( l_corr ) then
            call get_force(Geo2hInt,GeoRms,GeoRmsL,GeoRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Magnetostrophic balance
         if ( l_corr .and. l_mag_LF ) then
            call get_force(Mag2hInt,MagRms,MagRmsL,MagRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Coriolis/Lorentz balance:
         if ( l_corr .and. l_mag_LF ) then
            call get_force(CLF2hInt,CLFRms,CLFRmsL,CLFRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Pressure/Lorentz balance:
         if ( l_mag_LF ) then
            call get_force(PLF2hInt,PLFRms,PLFRmsL,PLFRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Buoyancy/Pressure/Coriolis/Lorentz balance:
         if (l_corr .and. l_mag_LF) then
            call get_force(ArcMag2hInt,ArcMagRms,ArcMagRmsL,ArcMagRmsLnR,  &
                 &         volC,nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Buoyancy/Pressure/Coriolis balance:
         if ( l_corr ) then
            call get_force(Arc2hInt,ArcRms,ArcRmsL,ArcRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Coriolis/Inertia/Archimedian balance:
         if (l_corr) then
            call get_force(CIA2hInt,CIARms,CIARmsL,CIARmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Advection+du/dt
         if ( l_conv_nl ) then
            call get_force(Iner2hInt,InerRms,InerRmsL,InerRmsLnR,volC, &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Viscosity
         call get_force(Dif2hInt,DifRms,DifRmsL,DifRmsLnR,volC,    &
              &         nRMS_sets,timePassed,timeNorm,l_stop_time, &
              &         DifPol2hInt,DifTor2hInt)


         !----- Output:
         if ( l_save_out ) then
            open(newunit=n_dtvrms_file, file=dtvrms_file, &
            &    form='formatted', status='unknown', position='append')
         end if
         write(n_dtvrms_file,'(1P,ES20.12,8ES16.8,7ES14.6)')          &
         &     time*tScale, InerRms, CorRms, LFRms, AdvRms, DifRms,   &
         &     Buo_tempRms, Buo_xiRms, PreRms, GeoRms/(CorRms+PreRms),&
         &     MagRms/(CorRms+PreRms+LFRms),                          &
         &     ArcRms/(CorRms+PreRms+Buo_tempRms+Buo_xiRms),          &
         &     ArcMagRms/(CorRms+PreRms+LFRms+Buo_tempRms+Buo_xiRms), &
         &     CLFRms/(CorRms+LFRms), PLFRms/(PreRms+LFRms),          &
         &     CIARms/(CorRms+PreRms+Buo_tempRms+Buo_xiRms+InerRms+LFRms)
         if ( l_save_out) then
            close(n_dtvrms_file)
         end if

         if ( l_stop_time ) then
            fileName='dtVrms_spec.'//tag
            open(newunit=fileHandle,file=fileName,form='formatted',status='unknown')
            do l=0,l_max
               write(fileHandle,'(1P,I4,30ES16.8)') l+1,                   &
               &     InerRmsL%mean(l),CorRmsL%mean(l),LFRmsL%mean(l),      &
               &     AdvRmsL%mean(l),DifRmsL%mean(l),Buo_tempRmsL%mean(l), &
               &     Buo_xiRmsL%mean(l), PreRmsL%mean(l),                  &
               &     GeoRmsL%mean(l),MagRmsL%mean(l),                      &
               &     ArcRmsL%mean(l),ArcMagRmsL%mean(l),CLFRmsL%mean(l),   &
               &     PLFRmsL%mean(l),CIARmsL%mean(l),InerRmsL%SD(l),       &
               &     CorRmsL%SD(l),LFRmsL%SD(l),AdvRmsL%SD(l),             &
               &     DifRmsL%SD(l),Buo_tempRmsL%SD(l), Buo_xiRmsL%SD(l),   &
               &     PreRmsL%SD(l), GeoRmsL%SD(l), MagRmsL%SD(l),          &
               &     ArcRmsL%SD(l),  ArcMagRmsL%SD(l),CLFRmsL%SD(l),       &
               &     PLFRmsL%SD(l), CIARmsL%SD(l)
            end do
            close(fileHandle)
         end if

         if ( l_2D_RMS .and. l_stop_time ) then
            version = 1
            fileName='2D_dtVrms_spec.'//tag
            open(newunit=fileHandle,file=fileName,form='unformatted', &
            &    status='unknown')
            write(fileHandle) version
            write(fileHandle) n_r_max, l_max
            write(fileHandle) r
            write(fileHandle) CorRmsLnR%mean(:,:)
            write(fileHandle) AdvRmsLnR%mean(:,:)
            write(fileHandle) LFRmsLnR%mean(:,:)
            write(fileHandle) Buo_tempRmsLnR%mean(:,:)
            write(fileHandle) Buo_xiRmsLnR%mean(:,:)
            write(fileHandle) PreRmsLnR%mean(:,:)
            write(fileHandle) DifRmsLnR%mean(:,:)
            write(fileHandle) InerRmsLnR%mean(:,:)
            write(fileHandle) GeoRmsLnR%mean(:,:)
            write(fileHandle) MagRmsLnR%mean(:,:)
            write(fileHandle) ArcRmsLnR%mean(:,:)
            write(fileHandle) ArcMagRmsLnR%mean(:,:)
            write(fileHandle) CIARmsLnR%mean(:,:)
            write(fileHandle) CLFRmsLnR%mean(:,:)
            write(fileHandle) PLFRmsLnR%mean(:,:)
            close(fileHandle)
         end if

      end if

   end subroutine dtVrms
!----------------------------------------------------------------------------
   subroutine dtBrms(time)

      !-- Input of variables:
      real(cp), intent(in) :: time

      !-- Local
      integer :: nR,l1m0,l1m1,lm,m

      real(cp) :: dtBPolRms,dtBPolAsRms
      real(cp) :: dtBTorRms,dtBTorAsRms
      real(cp) :: DdynRms,DdynAsRms
      real(cp) :: PdynRms,PdynAsRms
      real(cp) :: TdynRms,TdynAsRms
      real(cp) :: dummy1,dummy2,dummy3

      complex(cp) :: PdynLM(llmMag:ulmMag,n_r_max_dtB)
      complex(cp) :: drPdynLM(llmMag:ulmMag,n_r_max_dtB)
      complex(cp) :: TdynLM(llmMag:ulmMag,n_r_max_dtB)
      complex(cp) :: work_LMloc(llmMag:ulmMag,n_r_max_dtB)
      complex(cp) :: work_Rloc(lm_maxMag,nRstartMag:nRstopMag)

      real(cp) :: dtBP(n_r_max),dtBPAs(n_r_max)
      real(cp) :: dtBT(n_r_max),dtBTAs(n_r_max)
      real(cp) :: dtBP_global(n_r_max),dtBPAs_global(n_r_max)
      real(cp) :: dtBT_global(n_r_max),dtBTAs_global(n_r_max)

      real(cp) :: PdifRms, PdifAsRms, TdifRms, TdifAsRms, TomeRms, TomeAsRms

      !--- Stretching
      call get_dr(PstrLM_LMloc(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
           &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,      &
           &      nocopy=.true.)

      !--- Add to the total dynamo term
      do nR=1,n_r_max
         do lm=llmMag,ulmMag
            PdynLM(lm,nR)  =PstrLM_LMloc(lm,nR)
            drPdynLM(lm,nR)=work_LMloc(lm,nR)
         end do
      end do

      !-- Finalize advection
      call get_dr(PadvLM_LMloc(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
           &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,      &
           &      nocopy=.true.)

      !-- Add to total dynamo term:
      do nR=1,n_r_max
         do lm=llmMag,ulmMag
            PdynLM(lm,nR)  =PdynLM(lm,nR)-PadvLM_LMloc(lm,nR)
            drPdynLM(lm,nR)=drPdynLM(lm,nR)-work_LMloc(lm,nR)
            TdynLM(lm,nR)  =TstrLM_LMloc(lm,nR)-TadvLM_LMloc(lm,nR)
         end do
      end do

      !--- Get RMS values of the total dynamo term:
      call get_PolTorRms(PdynLM,drPdynLM,TdynLM,llmMag,ulmMag,PdynRms,TdynRms, &
           &             PdynAsRms,TdynAsRms,lo_map)


      !--- Finalize diffusion:
      call get_dr(PdifLM_LMloc(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
           &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,      &
           &      nocopy=.true.)

      !-- Get RMS values for diffusion
      call get_PolTorRms(PdifLM_LMloc,work_LMloc,TdifLM_LMloc,llmMag,ulmMag, &
           &             PdifRms,TdifRms,PdifAsRms,TdifAsRms,lo_map)

      !--- Get Omega effect rms: total toroidal field changes due to zonal flow
      !    (this is now stretching plus advection, changed May 23 2013):
      !    TomeAsRms is the rms of the more classical Omega effect which
      !    decribes the creation of axisymmetric azimuthal field by zonal flow.
      call get_PolTorRms(PdifLM_LMloc,work_LMloc,TomeLM_LMloc,llmMag,ulmMag, &
           &             dummy1,TomeRms,dummy2,TomeAsRms,lo_map)

      !--- B changes:
      if ( l_mag_par_solve ) then
         call get_dr_Rloc(dtBPolLMr,work_Rloc,lm_maxMag,nRstartMag,nRstopMag, &
              &           n_r_max,rscheme_oc)
         do nR=nRstartMag,nRstopMag
            call hInt2dPolLM(work_Rloc(:,nR),1,lm_max,dtBPol2hInt(:,nR),st_map)
            dtBP(nR)  =0.0_cp
            dtBT(nR)  =0.0_cp
            dtBPAs(nR)=0.0_cp
            dtBTAs(nR)=0.0_cp
            do lm=1,lm_maxMag
               m=st_map%lm2m(lm)
               dtBP(nR)=dtBP(nR)+dtBPol2hInt(lm,nR)
               dtBT(nR)=dtBT(nR)+dtBTor2hInt(lm,nR)
               if ( m == 0 ) then
                  dtBPAs(nR)=dtBPAs(nR)+dtBPol2hInt(lm,nR)
                  dtBTAs(nR)=dtBTAs(nR)+dtBTor2hInt(lm,nR)
               end if
            end do
         end do
      else
         call get_dr(dtBPolLMr(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
              &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,   &
              &      nocopy=.true.)

         do nR=1,n_r_max
            call hInt2dPolLM(work_LMloc(llm:ulm,nR),llm,ulm, &
                 &           dtBPol2hInt(llm:ulm,nR),lo_map)
            dtBP(nR)  =0.0_cp
            dtBT(nR)  =0.0_cp
            dtBPAs(nR)=0.0_cp
            dtBTAs(nR)=0.0_cp
            do lm=llm,ulm
               m=lo_map%lm2m(lm)
               dtBP(nR)=dtBP(nR)+dtBPol2hInt(lm,nR)
               dtBT(nR)=dtBT(nR)+dtBTor2hInt(lm,nR)
               if ( m == 0 ) then
                  dtBPAs(nR)=dtBPAs(nR)+dtBPol2hInt(lm,nR)
                  dtBTAs(nR)=dtBTAs(nR)+dtBTor2hInt(lm,nR)
               end if
            end do
         end do
      end if

#ifdef WITH_MPI
      call MPI_Reduce(dtBP, dtBP_global, n_r_max, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(dtBT, dtBT_global, n_r_max, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(dtBPAs, dtBPAs_global, n_r_max, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(dtBTAs, dtBTAs_global, n_r_max, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
#else
      dtBP_global(:)  =dtBP(:)
      dtBT_global(:)  =dtBT(:)
      dtBPAs_global(:)=dtBPAs(:)
      dtBTAs_global(:)=dtBTAs(:)
#endif

      if ( rank == 0 ) then
         dtBPolRms  =rInt_R(dtBP_global,r,rscheme_oc)
         dtBPolAsRms=rInt_R(dtBPAs_global,r,rscheme_oc)
         dtBTorRms  =rInt_R(dtBT_global,r,rscheme_oc)
         dtBTorAsRms=rInt_R(dtBTAs_global,r,rscheme_oc)

         dtBPolRms  =sqrt(dtBPolRms  /vol_oc)
         dtBPolAsRms=sqrt(dtBPolAsRms/vol_oc)
         dtBTorRms  =sqrt(dtBTorRms  /vol_oc)
         dtBTorAsRms=sqrt(dtBTorAsRms/vol_oc)
      end if

      !-- Get dipole dynamo contribution:
      l1m0=lo_map%lm2(1,0)
      l1m1=lo_map%lm2(1,1)
      do nR=1,n_r_max
         do lm=llmMag,ulmMag
            if ( lm/=l1m0 .and. lm/=l1m1 ) then
               PdynLM(lm,nR)  =zero
               drPdynLM(lm,nR)=zero
            end if
         end do
      end do

      !-- Get dipole dynamo terms:
      call get_PolTorRms(PdynLM,drPdynLM,TdynLM,llm,ulm,DdynRms,dummy1, &
           &             DdynAsRms,dummy3,lo_map)

      if ( rank == 0 ) then
         !-- Output:
         if ( l_save_out) then
            open(newunit=n_dtbrms_file, file=dtbrms_file,  &
            &    form='formatted', status='unknown', position='append')
         end if
         write(n_dtbrms_file,'(1P,ES20.12,10ES16.8)')               &
         &     time*tScale, dtBPolRms, dtBTorRms, PdynRms, TdynRms, &
         &     PdifRms, TdifRms, TomeRms/TdynRms,                   &
         &     TomeAsRms/TdynRms,  DdynRms,DdynAsRms
         if ( l_save_out) close(n_dtbrms_file)

      end if

   end subroutine dtBrms
!----------------------------------------------------------------------------
end module RMS
