module RMS
   !
   ! This module contains the calculation of the RMS force balance and induction
   ! terms.
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use blocking, only: st_map, nThetaBs, nfs, sizeThetaB, lo_map, lm2, &
       &               lm2m
   use finite_differences, only: type_fd
   use chebyshev, only: type_cheb_odd
   use radial_scheme, only: type_rscheme
   use truncation, only: n_r_max, n_cheb_max, n_r_maxMag, lm_max, lm_maxMag, &
       &                 l_max, n_phi_max, n_theta_max, minc, n_r_max_dtB,   &
       &                 lm_max_dtB, fd_ratio, fd_stretch
   use physical_parameters, only: ra, ek, pr, prmag, radratio
   use radial_data, only: nRstop, nRstart, radial_balance
   use radial_functions, only: rscheme_oc, r, r_cmb, r_icb
   use logic, only: l_save_out, l_heat, l_conv_nl, l_mag_LF, l_conv, &
       &            l_corr, l_mag, l_finite_diff, l_newmap, l_2D_RMS
   use num_param, only: tScale, alph1, alph2
   use horizontal_data, only: phi, theta_ord
   use constants, only: zero, one, half, four, third, vol_oc, pi
   use integration, only: rInt_R
   use radial_der, only: get_dr
   use output_data, only: rDea, rCut, tag, runid
   use cosine_transform_odd
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use RMS_helpers, only: hInt2dPol, get_PolTorRms, hInt2dPolLM
   use dtB_mod, only: PdifLM_LMloc, TdifLM_LMloc, PstrLM_LMloc, PadvLM_LMloc, &
       &              TadvLM_LMloc, TstrLM_LMloc, TomeLM_LMloc
   use useful, only: abortRun
   use mean_sd, only: mean_sd_type, mean_sd_2D_type
                                                                  
   implicit none
 
   private

   class(type_rscheme), pointer :: rscheme_RMS
   integer, public :: n_r_maxC     ! Number of radial points
   integer, public :: n_cheb_maxC  ! Number of Chebyshevs
   integer, public :: nCut         ! Number of points for the cut-off

   real(cp), allocatable :: rC(:)        ! Cut-off radii
   real(cp), public, allocatable :: dr_facC(:)   

   real(cp), public, allocatable :: dtBPol2hInt(:,:)
   real(cp), public, allocatable :: dtBTor2hInt(:,:)
   complex(cp), public, allocatable :: dtBPolLMr(:,:)
 
   real(cp), public, allocatable :: dtVPol2hInt(:,:)
   real(cp), public, allocatable :: dtVTor2hInt(:,:)
   complex(cp), public, allocatable :: dtVPolLMr(:,:)

   real(cp), public, allocatable :: DifPol2hInt(:,:)
   real(cp), public, allocatable :: DifTor2hInt(:,:)
   complex(cp), public, allocatable :: DifPolLMr(:,:)
 
   real(cp), public, allocatable :: Adv2hInt(:,:)
   real(cp), public, allocatable :: Cor2hInt(:,:)
   real(cp), public, allocatable :: LF2hInt(:,:)
   real(cp), public, allocatable :: Buo2hInt(:,:)
   real(cp), public, allocatable :: Pre2hInt(:,:)
   real(cp), public, allocatable :: Geo2hInt(:,:)
   real(cp), public, allocatable :: Mag2hInt(:,:)
   real(cp), public, allocatable :: Arc2hInt(:,:)
   real(cp), public, allocatable :: ArcMag2hInt(:,:)
   real(cp), public, allocatable :: CIA2hInt(:,:)
   real(cp), public, allocatable :: CLF2hInt(:,:)
   real(cp), public, allocatable :: PLF2hInt(:,:)

   !-- Time-averaged spectra
   type(mean_sd_type) :: dtVRmsL, CorRmsL, LFRmsL, AdvRmsL
   type(mean_sd_type) :: DifRmsL, BuoRmsL, PreRmsL, GeoRmsL
   type(mean_sd_type) :: MagRmsL, ArcRmsL, ArcMagRmsL
   type(mean_sd_type) :: CIARmsL, CLFRmsL, PLFRmsL

   type(mean_sd_2D_type) :: CorRmsLnR, AdvRmsLnR, LFRmsLnR, BuoRmsLnR
   type(mean_sd_2D_type) :: PreRmsLnR, DifRmsLnR, dtVRmsLnR, GeoRmsLnR
   type(mean_sd_2D_type) :: MagRmsLnR, ArcRmsLnR, ArcMagRmsLnR
   type(mean_sd_2D_type) :: CIARmsLnR, CLFRmsLnR, PLFRmsLnR

   integer :: n_dtvrms_file, n_dtbrms_file
   character(len=72) :: dtvrms_file, dtbrms_file
 
   public :: dtVrms, dtBrms, initialize_RMS, zeroRms, finalize_RMS

contains

   subroutine initialize_RMS
      !
      ! Memory allocation
      !

      allocate( dtBPol2hInt(llmMag:ulmMag,n_r_maxMag) )
      allocate( dtBTor2hInt(llmMag:ulmMag,n_r_maxMag) )
      allocate( dtBPolLMr(llmMag:ulmMag,n_r_maxMag) )
      bytes_allocated = bytes_allocated+ &
      &                 2*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_REAL+&
      &                 (llmMag-ulmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
    
      allocate( dtVPol2hInt(0:l_max,n_r_max) )
      allocate( dtVTor2hInt(0:l_max,n_r_max) )
      allocate( dtVPolLMr(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+ &
      &                 2*(l_max+1)*n_r_max*SIZEOF_DEF_REAL+  &
      &                 (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

      allocate( DifPol2hInt(0:l_max,n_r_max) )
      allocate( DifTor2hInt(0:l_max,n_r_max) )
      allocate( DifPolLMr(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+ &
      &                 2*(l_max+1)*n_r_max*SIZEOF_DEF_REAL+  &
      &                 (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
    
      allocate( Adv2hInt(0:l_max,n_r_max) )
      allocate( Cor2hInt(0:l_max,n_r_max) )
      allocate( LF2hInt(0:l_max,n_r_max) )
      allocate( Buo2hInt(0:l_max,n_r_max) )
      allocate( Pre2hInt(0:l_max,n_r_max) )
      allocate( Geo2hInt(0:l_max,n_r_max) )
      allocate( Mag2hInt(0:l_max,n_r_max) )
      allocate( Arc2hInt(0:l_max,n_r_max) )
      allocate( ArcMag2hInt(0:l_max,n_r_max) )
      allocate( CIA2hInt(0:l_max,n_r_max) )
      allocate( CLF2hInt(0:l_max,n_r_max) )
      allocate( PLF2hInt(0:l_max,n_r_max) )
      bytes_allocated = bytes_allocated+ 12*(l_max+1)*n_r_max*SIZEOF_DEF_REAL

      call dtVRmsL%initialize(0,l_max)
      call CorRmsL%initialize(0,l_max)
      call LFRmsL%initialize(0,l_max)
      call AdvRmsL%initialize(0,l_max)
      call DifRmsL%initialize(0,l_max)
      call BuoRmsL%initialize(0,l_max)
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
         call BuoRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call PreRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call DifRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
         call dtVRmsLnR%initialize(0,l_max,1,n_r_max,.false.)
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

   end subroutine initialize_RMS
!----------------------------------------------------------------------------
   subroutine finalize_RMS

      deallocate( rC )
      deallocate( dtBPol2hInt, dtBTor2hInt, dtBPolLMr )
      deallocate( dtVPol2hInt, dtVTor2hInt, dtVPolLMr )
      deallocate( DifPol2hInt, DifTor2hInt, DifPolLMr )
      deallocate( Adv2hInt, Cor2hInt, LF2hInt )
      deallocate( Buo2hInt, Pre2hInt, Geo2hInt )
      deallocate( Mag2hInt, Arc2hInt, CIA2hInt )
      deallocate( CLF2hInt, PLF2hInt, ArcMag2hInt )

      call dtVRmsL%finalize()
      call CorRmsL%finalize()
      call LFRmsL%finalize()
      call AdvRmsL%finalize()
      call DifRmsL%finalize()
      call BuoRmsL%finalize()
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
          call BuoRmsLnR%finalize()
          call PreRmsLnR%finalize()
          call DifRmsLnR%finalize()
          call dtVRmsLnR%finalize()
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

      !-- Local variables
      integer :: nR,lm,l

      do nR=1,n_r_max
         do l=0,l_max
            dtVPol2hInt(l,nR)=0.0_cp
            dtVTor2hInt(l,nR)=0.0_cp
            DifPol2hInt(l,nR)=0.0_cp
            DifTor2hInt(l,nR)=0.0_cp
         end do
      end do
      do nR=1,n_r_maxMag
         do lm=llmMag,ulmMag
            dtBPol2hInt(lm,nR)=0.0_cp
            dtBTor2hInt(lm,nR)=0.0_cp
         end do
      end do

      do nR=1,n_r_max
         do l=0,l_max
            Adv2hInt(l,nR)   =0.0_cp
            Cor2hInt(l,nR)   =0.0_cp
            LF2hInt(l,nR)    =0.0_cp
            Buo2hInt(l,nR)   =0.0_cp
            Pre2hInt(l,nR)   =0.0_cp
            Geo2hInt(l,nR)   =0.0_cp
            Mag2hInt(l,nR)   =0.0_cp
            Arc2hInt(l,nR)   =0.0_cp
            ArcMag2hInt(l,nR)=0.0_cp
            CIA2hInt(l,nR)   =0.0_cp
            CLF2hInt(l,nR)   =0.0_cp
            PLF2hInt(l,nR)   =0.0_cp
         end do
      end do
      do nR=1,n_r_max
         do lm=llm,ulm
            dtVPolLMr(lm,nR)=zero
            DifPolLMr(lm,nR)=zero
         end do
      end do
      do nR=1,n_r_maxMag
         do lm=llmMag,ulmMag
            dtBPolLMr(lm,nR)=zero
         end do
      end do

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
      real(cp) :: dtV_Rms  =0.0_cp
      real(cp) :: CorRms   =0.0_cp
      real(cp) :: AdvRms   =0.0_cp
      real(cp) :: LFRms    =0.0_cp
      real(cp) :: DifRms   =0.0_cp
      real(cp) :: BuoRms   =0.0_cp
      real(cp) :: PreRms   =0.0_cp
      real(cp) :: GeoRms   =0.0_cp
      real(cp) :: MagRms   =0.0_cp
      real(cp) :: ArcRms   =0.0_cp
      real(cp) :: ArcMagRms=0.0_cp
      real(cp) :: CIARms   =0.0_cp
      real(cp) :: CLFRms   =0.0_cp
      real(cp) :: PLFRms   =0.0_cp

      !-- Local variables:
      integer :: nR,nRC,l,fileHandle
      real(cp) :: volC
      real(cp) :: Dif2hInt(n_r_max),dtV2hInt(n_r_max)

      complex(cp) :: workA(llm:ulm,n_r_max)
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
         call get_dr(DifPolLMr(llm:ulm,:),workA(llm:ulm,:),ulm-llm+1,1, &
              &      ulm-llm+1,n_r_max,rscheme_oc)
      end if

      do nR=1,n_r_max
         call hInt2dPol(workA(llm:ulm,nR),llm,ulm,DifPol2hInt(:,nR), &
              &         lo_map)
      end do
#ifdef WITH_MPI
      call MPI_Reduce(DifPol2hInt(:,:),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPol2hInt(:,:)=global_sum
#endif

      !-- Flow changes
      dtV_Rms=0.0_cp
      if ( rscheme_RMS%version == 'cheb' ) then
         call get_dr(dtVPolLMr(llm:ulm,:),workA(llm:ulm,:),ulm-llm+1,1, &
              &      ulm-llm+1,n_r_max,rscheme_oc,nocopy=.true.)
      else
         call get_dr(dtVPolLMr(llm:ulm,:),workA(llm:ulm,:),ulm-llm+1,1, &
              &      ulm-llm+1,n_r_max,rscheme_oc)
      end if

      do nR=1,n_r_max
         call hInt2dPol( workA(llm:ulm,nR),llm,ulm,dtVPol2hInt(:,nR), &
              &          lo_map)
      end do
#ifdef WITH_MPI
      call MPI_Reduce(dtVPol2hInt(:,:),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPol2hInt(:,:)=global_sum
#endif

      ! First gather all needed arrays on rank 0
      ! some more arrays to gather for the dtVrms routine
      ! we need some more fields for the dtBrms routine
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
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Cor2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Adv2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              LF2hInt,recvcounts,displs,MPI_DEF_REAL,    &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,       &
           &              Buo2hInt,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
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
    
      ! The following fields are LM distributed and have to be gathered:
      ! dtVPolLMr, DifPolLMr
    
      call MPI_Reduce(dtVTor2hInt(:,:),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTor2hInt(:,:)=global_sum
      call MPI_Reduce(DifTor2hInt(:,:),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTor2hInt(:,:)=global_sum
#endif
    
      if ( rank == 0 ) then
    
         !write(*,"(A,ES22.14)") "dtVPol2hInt = ",SUM(dtVPol2hInt)
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
   
         !-- Buoyancy
         if ( l_conv_nl ) then
            call get_force(Buo2hInt,BuoRms,BuoRmsL,BuoRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
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
         if ( l_conv_nl .and. l_corr .and. l_mag_LF ) then
            call get_force(ArcMag2hInt,ArcMagRms,ArcMagRmsL,ArcMagRmsLnR,  &
                 &         volC,nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Buoyancy/Pressure/Coriolis balance:
         if ( l_conv_nl .and. l_corr ) then
            call get_force(Arc2hInt,ArcRms,ArcRmsL,ArcRmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Coriolis/Inertia/Archimedian balance:
         if ( l_conv_nl .and. l_corr ) then
            call get_force(CIA2hInt,CIARms,CIARmsL,CIARmsLnR,volC,    &
                 &         nRMS_sets,timePassed,timeNorm,l_stop_time)
         end if

         !-- Visosity
         call get_force(Dif2hInt,DifRms,DifRmsL,DifRmsLnR,volC,    &
              &         nRMS_sets,timePassed,timeNorm,l_stop_time, &
              &         DifPol2hInt,DifTor2hInt)

         !-- Time derivative
         call get_force(dtV2hInt,dtV_Rms,dtVRmsL,dtVRmsLnR,volC,   &
              &         nRMS_sets,timePassed,timeNorm,l_stop_time, &
              &         dtVPol2hInt,dtVTor2hInt)


         !----- Output:
         if ( l_save_out ) then
            open(newunit=n_dtvrms_file, file=dtvrms_file, &
            &    form='formatted', status='unknown', position='append')
         end if
         write(n_dtvrms_file,'(1P,ES20.12,7ES16.8,7ES14.6)')        &
         &     time*tScale, dtV_Rms, CorRms, LFRms, AdvRms, DifRms, &
         &     BuoRms, PreRms, GeoRms/(CorRms+PreRms),              &
         &     MagRms/(CorRms+PreRms+LFRms),                        &
         &     ArcRms/(CorRms+PreRms+BuoRms),                       &
         &     ArcMagRms/(CorRms+PreRms+LFRms+BuoRms),              &
         &     CLFRms/(CorRms+LFRms), PLFRms/(PreRms+LFRms),        &
         &     CIARms/(CorRms+PreRms+BuoRms+AdvRms+LFRms)
         if ( l_save_out) then
            close(n_dtvrms_file)
         end if

         if ( l_stop_time ) then
            fileName='dtVrms_spec.'//tag
            open(newunit=fileHandle,file=fileName,form='formatted',status='unknown')
            do l=0,l_max
               write(fileHandle,'(1P,I4,28ES16.8)') l+1,                 &
               &     dtVRmsL%mean(l),CorRmsL%mean(l),LFRmsL%mean(l),     &
               &     AdvRmsL%mean(l),DifRmsL%mean(l),BuoRmsL%mean(l),    &
               &     PreRmsL%mean(l),GeoRmsL%mean(l),MagRmsL%mean(l),    &
               &     ArcRmsL%mean(l),ArcMagRmsL%mean(l),CLFRmsL%mean(l), &
               &     PLFRmsL%mean(l),CIARmsL%mean(l),dtVRmsL%SD(l),      &
               &     CorRmsL%SD(l),LFRmsL%SD(l),AdvRmsL%SD(l),           &
               &     DifRmsL%SD(l),BuoRmsL%SD(l),PreRmsL%SD(l),          &
               &     GeoRmsL%SD(l),MagRmsL%SD(l),ArcRmsL%SD(l),          &
               &     ArcMagRmsL%SD(l),CLFRmsL%SD(l),PLFRmsL%SD(l),       &
               &     CIARmsL%SD(l)
            end do
            close(fileHandle)
         end if

         if ( l_2D_RMS ) then
            fileName='2D_dtVrms_spec.'//tag
            open(newunit=fileHandle,file=fileName,form='unformatted', &
            &    status='unknown')
            write(fileHandle) n_r_max, l_max
            write(fileHandle) r
            write(fileHandle) CorRmsLnR%mean(:,:)
            write(fileHandle) AdvRmsLnR%mean(:,:)
            write(fileHandle) LFRmsLnR%mean(:,:)
            write(fileHandle) BuoRmsLnR%mean(:,:)
            write(fileHandle) PreRmsLnR%mean(:,:)
            write(fileHandle) DifRmsLnR%mean(:,:)
            write(fileHandle) dtVRmsLnR%mean(:,:)
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

      real(cp) :: dtBP(n_r_max),dtBPAs(n_r_max)
      real(cp) :: dtBT(n_r_max),dtBTAs(n_r_max)
      real(cp) :: dtBP_global(n_r_max),dtBPAs_global(n_r_max)
      real(cp) :: dtBT_global(n_r_max),dtBTAs_global(n_r_max)

      real(cp) :: PdifRms, PdifAsRms, TdifRms, TdifAsRms, TomeRms, TomeAsRms
    
      !-- For new movie output
      ! character(len=80) :: fileName
      ! integer :: nField,nFields,nFieldSize
      ! integer :: nTheta,nThetaN,nThetaS,nThetaStart
      ! integer :: nPos, fileHandle
      ! real(cp) :: dumm(12),rS
      ! real(cp) :: fOut(n_theta_max*n_r_max)
      ! real(cp) :: outBlock(nfs)
      ! character(len=80) :: version
      ! logical :: lRmsMov
    
    
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


         !-- Output of movie files for axisymmetric toroidal field changes:
         !   Tstr,Tome,Tdyn=Tstr+Tadv,
         ! lRmsMov=.false.
         ! if ( lRmsMov ) then
    
         !    nFieldSize=n_theta_max*n_r_max
         !    nFields=7
         !    fileName='dtTas_mov.'//tag
         !    open(newunit=fileHandle, file=fileName, status='unknown', &
         !    &    form='unformatted')
    
         !    !------ Write header
         !    version='JW_Movie_Version_2'
         !    write(fileHandle) version
         !    dumm(1)=112           ! type of input
         !    dumm(2)=3             ! marker for constant phi plane
         !    dumm(3)=0.0_cp          ! surface constant
         !    dumm(4)=nFields       ! no of fields
         !    write(fileHandle) (real(dumm(n),kind=outp),n=1,4)
    
         !    !------ Define marker for output fields stored in movie field
         !    dumm(1)=101           ! Field marker for AS Br stretching
         !    dumm(2)=102           ! Field marker for AS Br dynamo term
         !    dumm(3)=103           ! Field marker for AS Br diffusion
         !    dumm(4)=104           ! Field marker for AS Bp stretching
         !    dumm(5)=105           ! Field marker for AS Bp dynamo term
         !    dumm(6)=106           ! Field marker for AS Bp omega effect
         !    dumm(7)=107           ! Field marker for AS Bp diffusion
         !    write(fileHandle) (real(dumm(n),kind=outp),n=1,nFields)
    
         !    !------ Now other info about grid and parameters:
         !    write(fileHandle) runid        ! run identifier
         !    dumm( 1)=n_r_max          ! total number of radial points
         !    dumm( 2)=n_r_max          ! no of radial point in outer core
         !    dumm( 3)=n_theta_max      ! no. of theta points
         !    dumm( 4)=n_phi_max        ! no. of phi points
         !    dumm( 5)=minc             ! imposed symmetry
         !    dumm( 6)=ra               ! control parameters
         !    dumm( 7)=ek               ! (for information only)
         !    dumm( 8)=pr               !      -"-
         !    dumm( 9)=prmag            !      -"-
         !    dumm(10)=radratio         ! ratio of inner / outer core
         !    dumm(11)=tScale           ! timescale
         !    write(fileHandle) (real(dumm(n),kind=outp),     n=1,11)
         !    write(fileHandle) (real(r(n)/r_cmb,kind=outp),  n=1,n_r_max)
         !    write(fileHandle) (real(theta_ord(n),kind=outp),n=1,n_theta_max)
         !    write(fileHandle) (real(phi(n),kind=outp),      n=1,n_phi_max)
    
         !    dumm(1)=1    ! time frame number for movie
         !    dumm(2)=0.0_cp ! time
         !    dumm(3)=0.0_cp
         !    dumm(4)=0.0_cp
         !    dumm(5)=0.0_cp
         !    dumm(6)=0.0_cp
         !    dumm(7)=0.0_cp
         !    dumm(8)=0.0_cp
         !    write(fileHandle) (real(dumm(n),kind=outp),n=1,8)
    
         !    !------ Loop over different output field:
         !    do nField=1,nFields
    
         !       !------ Loop over r and theta:
         !       do nR=1,n_r_max ! Loop over radial points
         !          rS=r(nR)
         !          do n=1,nThetaBs ! Loop over theta blocks
         !             nThetaStart=(n-1)*sizeThetaB+1
    
         !             !------ Convert from lm to theta block and store in outBlock:
         !             if ( nField == 1 ) then
         !                call get_RAS(PstrLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 2 ) then
         !                ! Note that PadvLM stores PdynLM=PstrLM+PadvLM at this point!
         !                call get_RAS(PdynLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 3 ) then
         !                ! Note that PdynLM stores PdifLM at this point!
         !                call get_RAS(PdifLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 4 ) then
         !                call get_PASLM(TstrLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 5 ) then
         !                call get_PASLM(TdynLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 6 ) then
         !                call get_PASLM(TomeLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 7 ) then
         !                call get_PASLM(TdifLM(:,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             end if
    
         !             !------ Storage of field in fout for theta block
         !             do nTheta=1,sizeThetaB,2
         !                !-- Convert to correct order in theta grid points 
         !                !-- and store of fOut:
         !                nThetaN=(nThetaStart+nTheta)/2
         !                nPos=(nR-1)*n_theta_max+nThetaN
         !                fOut(nPos)=outBlock(nTheta)
         !                nThetaS=n_theta_max-nThetaN+1
         !                nPos=(nR-1)*n_theta_max+nThetaS
         !                fOut(nPos)=outBlock(nTheta+1)
         !             end do ! Loop over thetas in block
    
         !          end do ! Loop over theta blocks
    
         !       end do ! Loop over R
    
         !       !------ Output of field:
         !       write(fileHandle) (real(fOut(nPos),kind=outp),nPos=1,nFieldSize)
    
         !    end do ! Loop over different fields

         !    close(fileHandle)
    
         ! end if ! output of mov fields ?

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
