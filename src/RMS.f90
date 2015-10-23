module RMS
   !
   ! This module contains the global array used when RMS force balance
   ! is requested
   !

   use precision_mod
   use truncation, only: n_r_max, n_cheb_max, n_r_maxMag, lm_max, lm_maxMag
   use constants, only: zero, one
   use chebyshev_polynoms_mod, only: cheb_grid
   use radial_functions, only: nDi_costf1, nDd_costf1, r
   use radial_der, only: get_dr
   use output_data, only: rDea, rCut
   use cosine_transform_odd
 
   implicit none
 
   private

   type(costf_odd_t), public :: chebt_RMS
   integer, public :: n_r_maxC     ! Number of radial points
   integer, public :: n_cheb_maxC  ! Number of Chebyshevs
   integer, public :: nCut         ! Number of points for the cut-off

   real(cp), allocatable :: rC(:)        ! Cut-off radii
   real(cp), public, allocatable :: dr_facC(:)   

   real(cp), public, allocatable :: dtBPol2hInt(:,:)
   real(cp), public, allocatable :: dtBTor2hInt(:,:)
   real(cp), public, allocatable :: dtBPolAs2hInt(:,:)
   real(cp), public, allocatable :: dtBTorAs2hInt(:,:)
   complex(cp), public, allocatable :: dtBPolLMr(:,:)
 
   real(cp), public, allocatable :: dtVPol2hInt(:,:)
   real(cp), public, allocatable :: dtVTor2hInt(:,:)
   real(cp), public, allocatable :: dtVPolAs2hInt(:,:)
   real(cp), public, allocatable :: dtVTorAs2hInt(:,:)
   complex(cp), public, allocatable :: dtVPolLMr(:,:)
 
   real(cp), public, allocatable :: AdvPol2hInt(:)
   real(cp), public, allocatable :: AdvTor2hInt(:)
   real(cp), public, allocatable :: AdvPolAs2hInt(:)
   real(cp), public, allocatable :: AdvTorAs2hInt(:)
   complex(cp), public, allocatable :: AdvPolLMr(:,:)
 
   real(cp), public, allocatable :: CorPol2hInt(:)
   real(cp), public, allocatable :: CorTor2hInt(:)
   real(cp), public, allocatable :: CorPolAs2hInt(:)
   real(cp), public, allocatable :: CorTorAs2hInt(:)
   complex(cp), public, allocatable :: CorPolLMr(:,:)
 
   real(cp), public, allocatable :: LFPol2hInt(:)
   real(cp), public, allocatable :: LFTor2hInt(:)
   real(cp), public, allocatable :: LFPolAs2hInt(:)
   real(cp), public, allocatable :: LFTorAs2hInt(:)
   complex(cp), public, allocatable :: LFPolLMr(:,:)
 
   real(cp), public, allocatable :: DifPol2hInt(:,:)
   real(cp), public, allocatable :: DifTor2hInt(:,:)
   real(cp), public, allocatable :: DifPolAs2hInt(:,:)
   real(cp), public, allocatable :: DifTorAs2hInt(:,:)
   complex(cp), public, allocatable :: DifPolLMr(:,:)
 
   real(cp), public, allocatable :: Buo2hInt(:)
   real(cp), public, allocatable :: BuoAs2hInt(:)
   complex(cp), public, allocatable :: BuoLMr(:,:)
 
   real(cp), public, allocatable :: Pre2hInt(:)
   real(cp), public, allocatable :: PreAs2hInt(:)
   complex(cp), public, allocatable :: PreLMr(:,:)
 
   real(cp), public, allocatable :: Geo2hInt(:)
   real(cp), public, allocatable :: GeoAs2hInt(:)
   complex(cp), public, allocatable :: GeoLMr(:,:)
 
   real(cp), public, allocatable :: Mag2hInt(:)
   real(cp), public, allocatable :: MagAs2hInt(:)
   complex(cp), public, allocatable :: MagLMr(:,:)
 
   real(cp), public, allocatable :: Arc2hInt(:)
   real(cp), public, allocatable :: ArcAs2hInt(:)
   complex(cp), public, allocatable :: ArcLMr(:,:)
 
   public :: initialize_RMS, zeroRms

contains

   subroutine initialize_RMS
      !
      ! Memory allocation
      !

      integer, parameter :: nThreadsMax=1

      allocate( dtBPol2hInt(n_r_maxMag,nThreadsMax) )
      allocate( dtBTor2hInt(n_r_maxMag,nThreadsMax) )
      allocate( dtBPolAs2hInt(n_r_maxMag,nThreadsMax) )
      allocate( dtBTorAs2hInt(n_r_maxMag,nThreadsMax) )
      allocate( dtBPolLMr(lm_maxMag,n_r_maxMag) )
    
      allocate( dtVPol2hInt(n_r_max,nThreadsMax) )
      allocate( dtVTor2hInt(n_r_max,nThreadsMax) )
      allocate( dtVPolAs2hInt(n_r_max,nThreadsMax) )
      allocate( dtVTorAs2hInt(n_r_max,nThreadsMax) )
      allocate( dtVPolLMr(lm_max,n_r_max) )
    
      allocate( AdvPol2hInt(n_r_max) )
      allocate( AdvTor2hInt(n_r_max) )
      allocate( AdvPolAs2hInt(n_r_max) )
      allocate( AdvTorAs2hInt(n_r_max) )
      allocate( AdvPolLMr(lm_max,n_r_max) )
    
      allocate( CorPol2hInt(n_r_max) )
      allocate( CorTor2hInt(n_r_max) )
      allocate( CorPolAs2hInt(n_r_max) )
      allocate( CorTorAs2hInt(n_r_max) )
      allocate( CorPolLMr(lm_max,n_r_max) )
    
      allocate( LFPol2hInt(n_r_max) )
      allocate( LFTor2hInt(n_r_max) )
      allocate( LFPolAs2hInt(n_r_max) )
      allocate( LFTorAs2hInt(n_r_max) )
      allocate( LFPolLMr(lm_max,n_r_maxMag) )
    
      allocate( DifPol2hInt(n_r_max,nThreadsMax) )
      allocate( DifTor2hInt(n_r_max,nThreadsMax) )
      allocate( DifPolAs2hInt(n_r_max,nThreadsMax) )
      allocate( DifTorAs2hInt(n_r_max,nThreadsMax) )
      allocate( DifPolLMr(lm_max,n_r_max) )
    
      allocate( Buo2hInt(n_r_max) )
      allocate( BuoAs2hInt(n_r_max) )
      allocate( BuoLMr(lm_max,n_r_max) )
    
      allocate( Pre2hInt(n_r_max) )
      allocate( PreAs2hInt(n_r_max) )
      allocate( PreLMr(lm_max,n_r_max) )
    
      allocate( Geo2hInt(n_r_max) )
      allocate( GeoAs2hInt(n_r_max) )
      allocate( GeoLMr(lm_max,n_r_max) )
    
      allocate( Mag2hInt(n_r_max) )
      allocate( MagAs2hInt(n_r_max) )
      allocate( MagLMr(lm_max,n_r_max) )
    
      allocate( Arc2hInt(n_r_max) )
      allocate( ArcAs2hInt(n_r_max) )
      allocate( ArcLMr(lm_max,n_r_max) )

      !--- Initialize new cut-back grid:
      call init_rNB(r,rCut,rDea,rC,n_r_maxC,n_cheb_maxC, &
                    nCut,dr_facC,chebt_RMS,nDi_costf1,nDd_costf1)

   end subroutine initialize_RMS
!----------------------------------------------------------------------------
   subroutine zeroRms
      !
      !  Zeros integrals that are set in get_td, update_z,
      !  update_wp, update_b, dtVrms and dtBrms
      !

      !-- Local variables
      integer :: n,nR,lm

      do n=1,1
         do nR=1,n_r_max
            dtVPol2hInt(nR,n)  =0.0_cp
            dtVTor2hInt(nR,n)  =0.0_cp
            dtVPolAs2hInt(nR,n)=0.0_cp
            dtVTorAs2hInt(nR,n)=0.0_cp
            DifPol2hInt(nR,n)  =0.0_cp
            DifTor2hInt(nR,n)  =0.0_cp
            DifPolAs2hInt(nR,n)=0.0_cp
            DifTorAs2hInt(nR,n)=0.0_cp
         end do
         do nR=1,n_r_maxMag
            dtBPol2hInt(nR,n)  =0.0_cp
            dtBTor2hInt(nR,n)  =0.0_cp
            dtBPolAs2hInt(nR,n)=0.0_cp
            dtBTorAs2hInt(nR,n)=0.0_cp
         end do
      end do

      do nR=1,n_r_max
         AdvPol2hInt(nR)  =0.0_cp
         AdvTor2hInt(nR)  =0.0_cp
         AdvPolAs2hInt(nR)=0.0_cp
         AdvTorAs2hInt(nR)=0.0_cp
         CorPol2hInt(nR)  =0.0_cp
         CorTor2hInt(nR)  =0.0_cp
         CorPolAs2hInt(nR)=0.0_cp
         CorTorAs2hInt(nR)=0.0_cp
         LFPol2hInt(nR)   =0.0_cp
         LFTor2hInt(nR)   =0.0_cp
         LFPolAs2hInt(nR) =0.0_cp
         LFTorAs2hInt(nR) =0.0_cp
         Buo2hInt(nR)     =0.0_cp
         BuoAs2hInt(nR)   =0.0_cp
         Pre2hInt(nR)     =0.0_cp
         PreAs2hInt(nR)   =0.0_cp
         Geo2hInt(nR)     =0.0_cp
         GeoAs2hInt(nR)   =0.0_cp
         Mag2hInt(nR)     =0.0_cp
         MagAs2hInt(nR)   =0.0_cp
         Arc2hInt(nR)     =0.0_cp
         ArcAs2hInt(nR)   =0.0_cp
      end do
      do nR=1,n_r_max
         do lm=1,lm_max
            dtVPolLMr(lm,nR)=zero
            AdvPolLMr(lm,nR)=zero
            CorPolLMr(lm,nR)=zero
            DifPolLMr(lm,nR)=zero
            BuoLMr(lm,nR)   =zero
            PreLMr(lm,nR)   =zero
            GeoLMr(lm,nR)   =zero
            ArcLMr(lm,nR)   =zero
         end do
      end do
      do nR=1,n_r_maxMag
         do lm=1,lm_maxMag
            dtBPolLMr(lm,nR)=zero
            LFPolLMr(lm,nR) =zero
         end do
      end do

   end subroutine zeroRms
!----------------------------------------------------------------------------
   subroutine init_rNB(r,rCut,rDea,r2,n_r_max2,n_cheb_max2, &
        &              nS,dr_fac2,chebt_RMS,nDi_costf1,     &
        &              nDd_costf1)
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
      integer,               intent(in) :: nDi_costf1,nDd_costf1
    
      !--- Output variables:
      integer,               intent(out) :: nS,n_r_max2,n_cheb_max2
      real(cp), allocatable, intent(out) :: r2(:)
      real(cp), allocatable, intent(out) :: dr_fac2(:)
      type(costf_odd_t),     intent(out) :: chebt_RMS
    
      ! Local stuff
      real(cp) :: drx(n_r_max)      ! first derivatives of x(r)
      real(cp) :: r2C(n_r_max)
      real(cp) :: r_cheb2(n_r_max)
      real(cp) :: dr2(n_r_max)
      real(cp) :: w1(n_r_max), w2(n_r_max)
      real(cp) :: r_icb2, r_cmb2, dr_fac
      integer :: nRs(16)
    
      logical :: lStop
      integer :: nR,n
    
      !--- New radial grid:
    
      !--- Find number of points to be cut away at either side:
      lStop=.true.
      do nS=1,(n_r_max-1)/2
         if ( r(1)-r(nS) > rCut ) then
            lStop=.false.
            exit
         end if
      end do
      if ( lStop ) then
         write(*,*) 'No nS found in init_rNB!'
         stop
      end if
      n_r_max2=n_r_max-2*nS
    
      ! Allowed number of radial grid points:
      nRs(1) =25
      nRs(2) =33
      nRs(3) =37
      nRs(4) =41
      nRs(5) =49
      nRs(6) =61
      nRs(7) =65
      nRs(8) =73
      nRs(9) =81
      nRs(10)=97
      nRs(11)=101
      nRs(12)=109
      nRs(13)=121
      nRs(14)=161
      nRs(15)=193
      nRs(16)=257
      lStop=.true.
      do n=14,1,-1
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
      allocate( dr_fac2(n_r_max2) )
    
      do nR=1,n_r_max2
         r2(nR)=r(nR+nS)
      end do
      r_icb2=r2(n_r_max2)
      r_cmb2=r2(1)
      call cheb_grid(r_icb2,r_cmb2,n_r_max2-1,r2C,r_cheb2, &
                     0.0_cp,0.0_cp,0.0_cp,0.0_cp)
      call chebt_RMS%initialize(n_r_max2, nDi_costf1, nDd_costf1)

      dr_fac=one
      do nR=1,n_r_max
         drx(nR)=one
      end do
      call get_dr(r2,dr2,n_r_max2,n_cheb_max2,w1,w2,chebt_RMS,drx)
      do nR=1,n_r_max2
         dr_fac2(nR)=one/dr2(nR)
      end do

   end subroutine init_rNB
!-----------------------------------------------------------------------------
end module RMS
