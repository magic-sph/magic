module RMS
   !
   ! This module contains the calculation of thr RMS force balance and induction
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
   use radial_data, only: nRstop, nRstart
   use radial_functions, only: rscheme_oc, r, r_cmb, r_icb, alph1, alph2
   use logic, only: l_save_out, l_heat, l_conv_nl, l_mag_LF, l_conv, &
       &            l_corr, l_mag, l_finite_diff, l_newmap
   use num_param, only: tScale
   use horizontal_data, only: phi, theta_ord
   use constants, only: zero, one, half, four, third, vol_oc, pi
   use integration, only: rInt_R
   use chebyshev_polynoms_mod, only: cheb_grid
   use radial_der, only: get_dr, get_drNS
   use output_data, only: rDea, rCut, tag, runid
   use cosine_transform_odd
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use RMS_helpers, only: hInt2dPol, get_PolTorRms, get_PASLM, get_RAS, &
       &                  hInt2dPolLM
   use dtB_mod, only: PdifLM_LMloc, TdifLM_LMloc, PstrLM_LMloc, PadvLM_LMloc, &
       &              TadvLM_LMloc, TstrLM_LMloc, TomeLM_LMloc
   use useful, only: getMSD2
                                                                  
   implicit none
 
   private

   class(type_rscheme), pointer :: rscheme_RMS
   integer, public :: n_r_maxC     ! Number of radial points
   integer, public :: n_cheb_maxC  ! Number of Chebyshevs
   integer, public :: nCut         ! Number of points for the cut-off

   real(cp), allocatable :: rC(:)        ! Cut-off radii
   real(cp), public, allocatable :: dr_facC(:)   

   real(cp), public, allocatable :: dtBPol2hInt(:,:,:)
   real(cp), public, allocatable :: dtBTor2hInt(:,:,:)
   complex(cp), public, allocatable :: dtBPolLMr(:,:)
 
   real(cp), public, allocatable :: dtVPol2hInt(:,:,:)
   real(cp), public, allocatable :: dtVTor2hInt(:,:,:)
   complex(cp), public, allocatable :: dtVPolLMr(:,:)

   real(cp), public, allocatable :: DifPol2hInt(:,:,:)
   real(cp), public, allocatable :: DifTor2hInt(:,:,:)
   complex(cp), public, allocatable :: DifPolLMr(:,:)
 
   real(cp), public, allocatable :: Adv2hInt(:,:)
   real(cp), public, allocatable :: Cor2hInt(:,:)
   real(cp), public, allocatable :: LF2hInt(:,:)
   real(cp), public, allocatable :: Buo2hInt(:,:)
   real(cp), public, allocatable :: Pre2hInt(:,:)
   real(cp), public, allocatable :: Geo2hInt(:,:)
   real(cp), public, allocatable :: Mag2hInt(:,:)
   real(cp), public, allocatable :: Arc2hInt(:,:)
   real(cp), public, allocatable :: CIA2hInt(:,:)
   real(cp), public, allocatable :: CLF2hInt(:,:)
   real(cp), public, allocatable :: PLF2hInt(:,:)

   !-- Time-averaged spectra
   real(cp), allocatable :: dtVRmsL_TA(:), dtVRmsL_SD(:), dtVRmsSD(:)
   real(cp), allocatable :: CorRmsL_TA(:), CorRmsL_SD(:), CorRmsSD(:)
   real(cp), allocatable :: LFRmsL_TA(:), LFRmsL_SD(:), LFRmsSD(:)
   real(cp), allocatable :: AdvRmsL_TA(:), AdvRmsL_SD(:), AdvRmsSD(:)
   real(cp), allocatable :: DifRmsL_TA(:), DifRmsL_SD(:), DifRmsSD(:)
   real(cp), allocatable :: BuoRmsL_TA(:), BuoRmsL_SD(:), BuoRmsSD(:)
   real(cp), allocatable :: PreRmsL_TA(:), PreRmsL_SD(:), PreRmsSD(:)
   real(cp), allocatable :: GeoRmsL_TA(:), GeoRmsL_SD(:), GeoRmsSD(:)
   real(cp), allocatable :: MagRmsL_TA(:), MagRmsL_SD(:), MagRmsSD(:)
   real(cp), allocatable :: ArcRmsL_TA(:), ArcRmsL_SD(:), ArcRmsSD(:)
   real(cp), allocatable :: CIARmsL_TA(:), CIARmsL_SD(:), CIARmsSD(:)
   real(cp), allocatable :: CLFRmsL_TA(:), CLFRmsL_SD(:), CLFRmsSD(:)
   real(cp), allocatable :: PLFRmsL_TA(:), PLFRmsL_SD(:), PLFRmsSD(:)

   real(cp), parameter :: eps = 10.0_cp*epsilon(one)
   integer :: n_dtvrms_file, n_dtbrms_file
   character(len=72) :: dtvrms_file, dtbrms_file
 
   public :: dtVrms, dtBrms, initialize_RMS, zeroRms, finalize_RMS

contains

   subroutine initialize_RMS
      !
      ! Memory allocation
      !

      integer, parameter :: nThreadsMax=1

      allocate( dtBPol2hInt(llmMag:ulmMag,n_r_maxMag,nThreadsMax) )
      allocate( dtBTor2hInt(llmMag:ulmMag,n_r_maxMag,nThreadsMax) )
      allocate( dtBPolLMr(llmMag:ulmMag,n_r_maxMag) )
      bytes_allocated = bytes_allocated+ &
                        2*(ulmMag-llmMag+1)*n_r_maxMag*nThreadsMax*SIZEOF_DEF_REAL+&
                        (llmMag-ulmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
    
      allocate( dtVPol2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( dtVTor2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( dtVPolLMr(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+ &
                        2*(l_max+1)*n_r_max*nThreadsMax*SIZEOF_DEF_REAL+&
                        (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

      allocate( DifPol2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( DifTor2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( DifPolLMr(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+ &
                        2*(l_max+1)*n_r_max*nThreadsMax*SIZEOF_DEF_REAL+&
                        (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
    
      allocate( Adv2hInt(0:l_max,n_r_max) )
      allocate( Cor2hInt(0:l_max,n_r_max) )
      allocate( LF2hInt(0:l_max,n_r_max) )
      allocate( Buo2hInt(0:l_max,n_r_max) )
      allocate( Pre2hInt(0:l_max,n_r_max) )
      allocate( Geo2hInt(0:l_max,n_r_max) )
      allocate( Mag2hInt(0:l_max,n_r_max) )
      allocate( Arc2hInt(0:l_max,n_r_max) )
      allocate( CIA2hInt(0:l_max,n_r_max) )
      allocate( CLF2hInt(0:l_max,n_r_max) )
      allocate( PLF2hInt(0:l_max,n_r_max) )
      bytes_allocated = bytes_allocated+ 11*(l_max+1)*n_r_max*SIZEOF_DEF_REAL

      allocate( dtVRmsL_TA(0:l_max), dtVRmsL_SD(0:l_max), dtVRmsSD(0:l_max) )
      allocate( CorRmsL_TA(0:l_max), CorRmsL_SD(0:l_max), CorRmsSD(0:l_max) )
      allocate( LFRmsL_TA(0:l_max), LFRmsL_SD(0:l_max), LFRmsSD(0:l_max) )
      allocate( AdvRmsL_TA(0:l_max), AdvRmsL_SD(0:l_max), AdvRmsSD(0:l_max) )
      allocate( DifRmsL_TA(0:l_max), DifRmsL_SD(0:l_max), DifRmsSD(0:l_max) )
      allocate( BuoRmsL_TA(0:l_max), BuoRmsL_SD(0:l_max), BuoRmsSD(0:l_max) )
      allocate( PreRmsL_TA(0:l_max), PreRmsL_SD(0:l_max), PreRmsSD(0:l_max) )
      allocate( GeoRmsL_TA(0:l_max), GeoRmsL_SD(0:l_max), GeoRmsSD(0:l_max) )
      allocate( MagRmsL_TA(0:l_max), MagRmsL_SD(0:l_max), MagRmsSD(0:l_max) )
      allocate( ArcRmsL_TA(0:l_max), ArcRmsL_SD(0:l_max), ArcRmsSD(0:l_max) )
      allocate( CIARmsL_TA(0:l_max), CIARmsL_SD(0:l_max), CIARmsSD(0:l_max) )
      allocate( CLFRmsL_TA(0:l_max), CLFRmsL_SD(0:l_max), CLFRmsSD(0:l_max) )
      allocate( PLFRmsL_TA(0:l_max), PLFRmsL_SD(0:l_max), PLFRmsSD(0:l_max) )
      bytes_allocated = bytes_allocated+ 39*(l_max+1)*SIZEOF_DEF_REAL

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
      deallocate( Buo2hInt, Pre2hInt, Geo2hInt)
      deallocate( Mag2hInt, Arc2hInt, CIA2hInt)
      deallocate( CLF2hInt, PLF2hInt)
      deallocate( dtVRmsL_TA, dtVRmsL_SD, dtVRmsSD )
      deallocate( CorRmsL_TA, CorRmsL_SD, CorRmsSD )
      deallocate( LFRmsL_TA, LFRmsL_SD, LFRmsSD )
      deallocate( AdvRmsL_TA, AdvRmsL_SD, AdvRmsSD )
      deallocate( DifRmsL_TA, DifRmsL_SD, DifRmsSD )
      deallocate( BuoRmsL_TA, BuoRmsL_SD, BuoRmsSD )
      deallocate( PreRmsL_TA, PreRmsL_SD, PreRmsSD )
      deallocate( GeoRmsL_TA, GeoRmsL_SD, GeoRmsSD )
      deallocate( MagRmsL_TA, MagRmsL_SD, MagRmsSD )
      deallocate( ArcRmsL_TA, ArcRmsL_SD, ArcRmsSD )
      deallocate( CIARmsL_TA, CIARmsL_SD, CIARmsSD )
      deallocate( CLFRmsL_TA, CLFRmsL_SD, CLFRmsSD )
      deallocate( PLFRmsL_TA, PLFRmsL_SD, PLFRmsSD )

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
      integer :: n,nR,lm,l

      do n=1,1
         do nR=1,n_r_max
            do l=0,l_max
               dtVPol2hInt(l,nR,n)  =0.0_cp
               dtVTor2hInt(l,nR,n)  =0.0_cp
               DifPol2hInt(l,nR,n)  =0.0_cp
               DifTor2hInt(l,nR,n)  =0.0_cp
            end do
         end do
         do nR=1,n_r_maxMag
            do lm=llmMag,ulmMag
               dtBPol2hInt(lm,nR,n)=0.0_cp
               dtBTor2hInt(lm,nR,n)=0.0_cp
            end do
         end do
      end do

      do nR=1,n_r_max
         do l=0,l_max
            Adv2hInt(l,nR)  =0.0_cp
            Cor2hInt(l,nR)  =0.0_cp
            LF2hInt(l,nR)   =0.0_cp
            Buo2hInt(l,nR)  =0.0_cp
            Pre2hInt(l,nR)  =0.0_cp
            Geo2hInt(l,nR)  =0.0_cp
            Mag2hInt(l,nR)  =0.0_cp
            Arc2hInt(l,nR)  =0.0_cp
            CIA2hInt(l,nR)  =0.0_cp
            CLF2hInt(l,nR)  =0.0_cp
            PLF2hInt(l,nR)  =0.0_cp
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
   subroutine init_rNB(r,rCut,rDea,r2,n_r_max2,n_cheb_max2, &
        &              nS,rscheme_RMS)
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
         write(*,*) 'No nS found in init_rNB!'
         stop
      end if
      nS=nS-1
      n_r_max2=n_r_max-2*nS
    
      if ( .not. l_finite_diff ) then
         ! Allowed number of radial grid points:
         nRs = [25, 33, 37, 41, 49, 61, 65, 73, 81, 97, 101, 109, 121,  &
                129, 145, 161, 181, 193, 201, 217, 241, 257, 289, 301,  &
                321, 325, 361, 385, 401, 433, 481, 501, 513, 541, 577,  &
                601, 641, 649, 721, 769, 801, 865, 901, 961, 973, 1001, &
                1025, 1081, 1153, 1201, 1281, 1297, 1441, 1501, 1537,   &
                1601, 1621, 1729, 1801, 1921, 1945, 2001, 2049]
         lStop=.true.
         do n=size(nRs),1,-1
            if ( nRs(n) <= n_r_max2 ) then
               lStop=.false.
               exit
            end if
         end do
         if ( lStop ) then
            write(*,*) 'No n_r_max2 found in init_rNB!'
            stop
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

         do nR=1,n_r_max2
            rscheme_RMS%drx(nR)=one
         end do

         call get_dr(r2,dr2,n_r_max2,rscheme_RMS)

         do nR=1,n_r_max2
            rscheme_RMS%drx(nR)=one/dr2(nR)
         end do

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
!-----------------------------------------------------------------------------
   subroutine dtVrms(time,nRMS_sets,timePassed,timeNorm)
      !
      ! This routine calculates and stores the different contributions
      ! of the forces entering the Navier-Stokes equation.
      !

      !-- Input variables:
      real(cp), intent(in) :: time
      real(cp), intent(in) :: timePassed
      real(cp), intent(in) :: timeNorm
      integer,  intent(inout) :: nRMS_sets
    
      !-- Output:
      real(cp) :: dtV_Rms,dtVRmsL
      real(cp) :: CorRms,CorRmsL
      real(cp) :: AdvRms,AdvRmsL
      real(cp) :: LFRms,LFRmsL
      real(cp) :: DifRms,DifRmsL
      real(cp) :: BuoRms,BuoRmsL
      real(cp) :: PreRms,PreRmsL
      real(cp) :: GeoRms,GeoRmsL
      real(cp) :: MagRms,MagRmsL
      real(cp) :: ArcRms,ArcRmsL
      real(cp) :: CIARms,CIARmsL
      real(cp) :: CLFRms,CLFRmsL
      real(cp) :: PLFRms,PLFRmsL
    
      !-- Local variables:
      integer :: nR,nRC,l,n,fileHandle
      real(cp) :: volC
      real(cp) :: Rms(n_r_max),Dif2hInt(n_r_max),dtV2hInt(n_r_max)
    
      complex(cp) :: workA(llm:ulm,n_r_max)
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(cp) :: global_sum(l_max+1,n_r_max)
      integer :: irank,sendcount
      character(len=80) :: fileName

      nRC=nCut+1

      !-- Diffusion
      DifRms=0.0_cp
      if ( rscheme_RMS%version == 'cheb' ) then
         call get_dr(DifPolLMr(llm:,nRC:),workA(llm:,nRC:), &
              &      ulm-llm+1,1,ulm-llm+1, &
              &      n_r_maxC,rscheme_RMS,nocopy=.true.)
      else
         call get_dr(DifPolLMr(llm:,:),workA(llm:,:), &
              &      ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_RMS)
      end if

      do nR=1,n_r_maxC
         call hInt2dPol( workA(llm:,nR+nCut),llm,ulm,DifPol2hInt(:,nR+nCut,1), &
              &           lo_map )
      end do
#ifdef WITH_MPI
      call MPI_Reduce(DifPol2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPol2hInt(:,:,1)=global_sum
#endif

      !-- Flow changes
      dtV_Rms=0.0_cp
      if ( rscheme_RMS%version == 'cheb' ) then
         call get_dr(dtVPolLMr(llm:,nRC:),workA(llm:,nRC:),ulm-llm+1,1,ulm-llm+1, &
              &      n_r_maxC,rscheme_RMS,nocopy=.true.)
      else
         call get_dr(dtVPolLMr(llm:,:),workA(llm:,:),ulm-llm+1,1,ulm-llm+1, &
              &      n_r_max,rscheme_RMS)
      end if

      do nR=1,n_r_maxC
         call hInt2dPol( workA(llm:,nR+nCut),llm,ulm,dtVPol2hInt(:,nR+nCut,1), &
              &          lo_map)
      end do
#ifdef WITH_MPI
      call MPI_Reduce(dtVPol2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPol2hInt(:,:,1)=global_sum
#endif

      ! First gather all needed arrays on rank 0
      ! some more arrays to gather for the dtVrms routine
      ! we need some more fields for the dtBrms routine
#ifdef WITH_MPI
    
      ! The following fields are only 1D and R distributed.
      sendcount  = (nRstop-nRstart+1)*(l_max+1)
      recvcounts = nR_per_rank*(l_max+1)
      recvcounts(n_procs-1) = nR_on_last_rank*(l_max+1)
      do irank=0,n_procs-1
         displs(irank) = irank*nR_per_rank*(l_max+1)
      end do
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Cor2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Adv2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & LF2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Buo2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Pre2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Geo2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Mag2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Arc2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CIA2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CLF2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & PLF2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
    
      ! The following fields are LM distributed and have to be gathered:
      ! dtVPolLMr, DifPolLMr
    
      call MPI_Reduce(dtVTor2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTor2hInt(:,:,1)=global_sum
      call MPI_Reduce(DifTor2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTor2hInt(:,:,1)=global_sum
#endif
    
      if ( rank == 0 ) then
    
         !write(*,"(A,ES22.14)") "dtVPol2hInt = ",SUM(dtVPol2hInt)
         nRMS_sets=nRMS_sets+1
         volC=four*third*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)
    
         CorRms=0.0_cp
         if ( l_corr ) then
            do l=0,l_max
               !-- Copy each mode on radial array
               do nR=1,n_r_max
                  Rms(nR)=Cor2hInt(l,nR)
               end do
               !-- Integrate in radius
               CorRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               !-- Add for total Rms
               CorRms = CorRms+CorRmsL
               !-- Finish Rms for mode l
               CorRmsL=sqrt(CorRmsL/volC)
               !-- Calculate time average and SD for mode l:
               call getMSD2(CorRmsL_TA(l),CorRmsL_SD(l),CorRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         CorRms=sqrt(CorRms/volC)

         !-- Advection
         AdvRms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Adv2hInt(l,nR)
               end do
               AdvRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               AdvRms =AdvRms+AdvRmsL
               AdvRmsL=sqrt(AdvRmsL/volC)
               call getMSD2(AdvRmsL_TA(l),AdvRmsL_SD(l),AdvRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         AdvRms=sqrt(AdvRms/volC)

         !-- Lorentz force
         LFRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=LF2hInt(l,nR)
               end do
               LFRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               LFRms =LFRms+LFRmsL
               LFRmsL=sqrt(LFRmsL/volC)
               call getMSD2(LFRmsL_TA(l),LFRmsL_SD(l),LFRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         LFRms=sqrt(LFRms/volC)
    
         !-- Buoyancy
         BuoRms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Buo2hInt(l,nR)
               end do
               BuoRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               BuoRms =BuoRms+BuoRmsL
               BuoRmsL=sqrt(BuoRmsL/volC)
               call getMSD2(BuoRmsL_TA(l),BuoRmsL_SD(l),BuoRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         BuoRms=sqrt(BuoRms/volC)
    
         !-- Pressure gradient
         PreRms=0.0_cp
         do l=0,l_max
            do nR=1,n_r_max
               Rms(nR)=Pre2hInt(l,nR)
            end do
            PreRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
            PreRms =PreRms+PreRmsL
            PreRmsL=sqrt(PreRmsL/volC)
            call getMSD2(PreRmsL_TA(l),PreRmsL_SD(l),PreRmsL, &
                         nRMS_sets,timePassed,timeNorm)
         end do
         PreRms=sqrt(PreRms/volC)

         !-- Geostrophic balance
         GeoRms=0.0_cp
         do l=0,l_max
            do nR=1,n_r_max
               Rms(nR)=Geo2hInt(l,nR)
            end do
            GeoRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
            GeoRms =GeoRms+GeoRmsL
            GeoRmsL=sqrt(GeoRmsL/volC)
            call getMSD2(GeoRmsL_TA(l),GeoRmsL_SD(l),GeoRmsL, &
                         nRMS_sets,timePassed,timeNorm)
         end do
         GeoRms=sqrt(GeoRms/volC)

         !-- Magnetostrophic balance
         MagRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Mag2hInt(l,nR)
               end do
               MagRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               MagRms =MagRms+MagRmsL
               MagRmsL=sqrt(MagRmsL/volC)
               call getMSD2(MagRmsL_TA(l),MagRmsL_SD(l),MagRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         MagRms=sqrt(MagRms/volC)

         !-- Coriolis/Lorentz balance:
         CLFRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=CLF2hInt(l,nR)
               end do
               CLFRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               CLFRms =CLFRms+CLFRmsL
               CLFRmsL=sqrt(CLFRmsL/volC)
               call getMSD2(CLFRmsL_TA(l),CLFRmsL_SD(l),CLFRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         CLFRms=sqrt(CLFRms/volC)

         !-- Pressure/Lorentz balance:
         PLFRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=PLF2hInt(l,nR)
               end do
               PLFRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               PLFRms =PLFRms+PLFRmsL
               PLFRmsL=sqrt(PLFRmsL/volC)
               call getMSD2(PLFRmsL_TA(l),PLFRmsL_SD(l),PLFRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         PLFRms=sqrt(PLFRms/volC)

         !-- Archimedian balance:
         ArcRms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Arc2hInt(l,nR)
               end do
               ArcRmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               ArcRms =ArcRms+ArcRmsL
               ArcRmsL=sqrt(ArcRmsL/volC)
               call getMSD2(ArcRmsL_TA(l),ArcRmsL_SD(l),ArcRmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         ArcRms=sqrt(ArcRms/volC)

         !-- Coriolis/Inertia/Archimedia balance:
         CIARms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=CIA2hInt(l,nR)
               end do
               CIARmsL=rInt_R(Rms(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
               CIARms =CIARms+CIARmsL
               CIARmsL=sqrt(CIARmsL/volC)
               call getMSD2(CIARmsL_TA(l),CIARmsL_SD(l),CIARmsL, &
                            nRMS_sets,timePassed,timeNorm)
            end do
         end if
         CIARms=sqrt(CIARms/volC)

         do l=0,l_max
            do nR=1,n_r_maxC
               Dif2hInt(nR+nCut)=0.0_cp
               do n=1,1
                  Dif2hInt(nR+nCut)=Dif2hInt(nR+nCut) +        &
                                    DifPol2hInt(l,nR+nCut,n) + &
                                    DifTor2hInt(l,nR+nCut,n)
               end do
           end do
           DifRmsL=rInt_R(Dif2hInt(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
           DifRms =DifRms+DifRmsL
           DifRmsL=sqrt(DifRmsL/volC)
           call getMSD2(DifRmsL_TA(l),DifRmsL_SD(l),DifRmsL, &
                        nRMS_sets,timePassed,timeNorm)
         end do
         DifRms=sqrt(DifRms/volC)

         do l=0,l_max
            do nR=1,n_r_maxC
               dtV2hInt(nR+nCut)=0.0_cp
               do n=1,1
                  dtV2hInt(nR+nCut)=dtV2hInt(nR+nCut) +         &
                  &                 dtVPol2hInt(l,nR+nCut,n) +  &
                  &                 dtVTor2hInt(l,nR+nCut,n)
               end do
            end do
            dtVRmsL=rInt_R(dtV2hInt(nRC:n_r_max-nRC+1),rC,rscheme_RMS)
            dtV_Rms =dtV_Rms+dtVRmsL
            dtVRmsL=sqrt(dtVRmsL/volC)
            call getMSD2(dtVRmsL_TA(l),dtVRmsL_SD(l),dtVRmsL, &
                 &       nRMS_sets,timePassed,timeNorm)
         end do
         dtV_Rms=sqrt(dtV_Rms/volC)
    
         !----- Output:
         if ( l_save_out) then
            open(newunit=n_dtvrms_file, file=dtvrms_file, &
            &    form='formatted', status='unknown', position='append')
         end if
         write(n_dtvrms_file,'(1P,ES20.12,7ES16.8,6ES14.6)')&
         &    time, dtV_Rms, CorRms, LFRms, AdvRms, DifRms, &
         &    BuoRms, PreRms, GeoRms/(CorRms+PreRms),       &
         &    MagRms/(CorRms+PreRms+LFRms),                 &
         &    ArcRms/(CorRms+PreRms+LFRms+BuoRms),          &
         &    CLFRms/(CorRms+LFRms), PLFRms/(PreRms+LFRms), &
         &    CIARms/(CorRms+PreRms+BuoRms+AdvRms+LFRms)
         if ( l_save_out) then
            close(n_dtvrms_file)
         end if

         !-- RMS time averaged spectra
         do l=0,l_max
            dtVRmsSD(l)=sqrt(dtVRmsL_SD(l)/timeNorm)
            CorRmsSD(l)=sqrt(CorRmsL_SD(l)/timeNorm)
            AdvRmsSD(l)=sqrt(AdvRmsL_SD(l)/timeNorm)
            DifRmsSD(l)=sqrt(DifRmsL_SD(l)/timeNorm)
            BuoRmsSD(l)=sqrt(BuoRmsL_SD(l)/timeNorm)
            PreRmsSD(l)=sqrt(PreRmsL_SD(l)/timeNorm)
            GeoRmsSD(l)=sqrt(GeoRmsL_SD(l)/timeNorm)
            ArcRmsSD(l)=sqrt(ArcRmsL_SD(l)/timeNorm)
            CIARmsSD(l)=sqrt(CIARmsL_SD(l)/timeNorm)
            if ( l_mag_LF ) then
               LFRmsSD(l) =sqrt(LFRmsL_SD(l)/timeNorm)
               MagRmsSD(l)=sqrt(MagRmsL_SD(l)/timeNorm)
               CLFRmsSD(l)=sqrt(CLFRmsL_SD(l)/timeNorm)
               PLFRmsSD(l)=sqrt(PLFRmsL_SD(l)/timeNorm)
            end if 
            dtVRmsL_TA(l)=max(dtVRmsL_TA(l),eps)
            CorRmsL_TA(l)=max(CorRmsL_TA(l),eps)
            LFRmsL_TA(l) =max(LFRmsL_TA(l),eps)
            AdvRmsL_TA(l)=max(AdvRmsL_TA(l),eps)
            DifRmsL_TA(l)=max(DifRmsL_TA(l),eps)
            BuoRmsL_TA(l)=max(BuoRmsL_TA(l),eps)
            PreRmsL_TA(l)=max(PreRmsL_TA(l),eps)
            GeoRmsL_TA(l)=max(GeoRmsL_TA(l),eps)
            MagRmsL_TA(l)=max(MagRmsL_TA(l),eps)
            ArcRmsL_TA(l)=max(ArcRmsL_TA(l),eps)
            CIARmsL_TA(l)=max(CIARmsL_TA(l),eps)
            CLFRmsL_TA(l)=max(CLFRmsL_TA(l),eps)
            PLFRmsL_TA(l)=max(PLFRmsL_TA(l),eps)
         end do

         fileName='dtVrms_spec.'//tag
         open(newunit=fileHandle,file=fileName,form='formatted', &
         &    status='unknown')
         do l=0,l_max
            write(fileHandle,'(1P,I4,26ES16.8)') l+1,                        &
            &     dtVRmsL_TA(l),CorRmsL_TA(l),LFRmsL_TA(l),AdvRmsL_TA(l),    &
            &     DifRmsL_TA(l),BuoRmsL_TA(l),PreRmsL_TA(l),GeoRmsL_TA(l),   &
            &     MagRmsL_TA(l),ArcRmsL_TA(l),CLFRmsL_TA(l),PLFRmsL_TA(l),   &
            &     CIARmsL_TA(l),dtVRmsSD(l),CorRmsSD(l),LFRmsSD(l),          &
            &     AdvRmsSD(l),DifRmsSD(l),BuoRmsSD(l),PreRmsSD(l),           &
            &     GeoRmsSD(l),MagRmsSD(l),ArcRmsSD(l),CLFRmsSD(l),           &
            &     PLFRmsSD(l),CIARmsSD(l)
         end do
         close(fileHandle)
    
      end if

   end subroutine dtVrms
!----------------------------------------------------------------------------
   subroutine dtBrms(time)

      !-- Input of variables:
      real(cp), intent(in) :: time
    
      !-- Local
      integer :: nR,n,l1m0,l1m1,lm,m
    
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
      call get_dr(PstrLM_LMloc(llmMag:,:),work_LMloc(llmMag:,:),ulmMag-llmMag+1, &
           &      1,ulmMag-llmMag+1,n_r_max,rscheme_oc,nocopy=.true.)
    
      !--- Add to the total dynamo term
      do nR=1,n_r_max
         do lm=llmMag,ulmMag
            PdynLM(lm,nR)  =PstrLM_LMloc(lm,nR)
            drPdynLM(lm,nR)=work_LMloc(lm,nR)
         end do
      end do

      !-- Finalize advection
      call get_dr(PadvLM_LMloc(llmMag:,:),work_LMloc(llmMag:,:),ulmMag-llmMag+1, &
           &      1,ulmMag-llmMag+1,n_r_max,rscheme_oc,nocopy=.true.)

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
      call get_dr(PdifLM_LMloc(llmMag:,:),work_LMloc(llmMag:,:),ulmMag-llmMag+1, &
           &      1,ulmMag-llmMag+1,n_r_max,rscheme_oc,nocopy=.true.)

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
      call get_dr(dtBPolLMr(llmMag:,:),work_LMloc(llmMag:,:),ulmMag-llmMag+1, &
           &      1,ulmMag-llmMag+1,n_r_max,rscheme_oc,nocopy=.true.)

      do nR=1,n_r_max
         call hInt2dPolLM(work_LMloc(llm:,nR),llm,ulm,dtBPol2hInt(llm:,nR,1),lo_map)
         dtBP(nR)  =0.0_cp
         dtBT(nR)  =0.0_cp
         dtBPAs(nR)=0.0_cp
         dtBTAs(nR)=0.0_cp
         do n=1,1
            do lm=llm,ulm
               m=lo_map%lm2m(lm)
               dtBP(nR)=dtBP(nR)+dtBPol2hInt(lm,nR,n)
               dtBT(nR)=dtBT(nR)+dtBTor2hInt(lm,nR,n)
               if ( m == 0 ) then
                  dtBPAs(nR)=dtBPAs(nR)+dtBPol2hInt(lm,nR,n)
                  dtBTAs(nR)=dtBTAs(nR)+dtBTor2hInt(lm,nR,n)
               end if
            end do
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
         !                call get_RAS(PstrLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 2 ) then
         !                ! Note that PadvLM stores PdynLM=PstrLM+PadvLM at this point!
         !                call get_RAS(PdynLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 3 ) then
         !                ! Note that PdynLM stores PdifLM at this point!
         !                call get_RAS(PdifLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 4 ) then
         !                call get_PASLM(TstrLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 5 ) then
         !                call get_PASLM(TdynLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 6 ) then
         !                call get_PASLM(TomeLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
         !             else if ( nField == 7 ) then
         !                call get_PASLM(TdifLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
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
         write(n_dtbrms_file,'(1P,ES20.12,10ES16.8)')            &
              time, dtBPolRms, dtBTorRms, PdynRms, TdynRms,      &
              PdifRms, TdifRms, TomeRms/TdynRms,                 &
              TomeAsRms/TdynRms,  DdynRms,DdynAsRms
         if ( l_save_out) close(n_dtbrms_file)

      end if

   end subroutine dtBrms
!----------------------------------------------------------------------------
end module RMS
