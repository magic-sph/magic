!$Id$
module RMS

  use truncation, only: n_r_max, n_r_maxMag, lm_max, lm_maxMag

  implicit none

  private

  real(kind=8), public, allocatable :: dtBPol2hInt(:,:)
  real(kind=8), public, allocatable :: dtBTor2hInt(:,:)
  real(kind=8), public, allocatable :: dtBPolAs2hInt(:,:)
  real(kind=8), public, allocatable :: dtBTorAs2hInt(:,:)
  complex(kind=8), public, allocatable :: dtBPolLMr(:,:)

  real(kind=8), public, allocatable :: dtVPol2hInt(:,:)
  real(kind=8), public, allocatable :: dtVTor2hInt(:,:)
  real(kind=8), public, allocatable :: dtVPolAs2hInt(:,:)
  real(kind=8), public, allocatable :: dtVTorAs2hInt(:,:)
  complex(kind=8), public, allocatable :: dtVPolLMr(:,:)

  real(kind=8), public, allocatable :: AdvPol2hInt(:)
  real(kind=8), public, allocatable :: AdvTor2hInt(:)
  real(kind=8), public, allocatable :: AdvPolAs2hInt(:)
  real(kind=8), public, allocatable :: AdvTorAs2hInt(:)
  complex(kind=8), public, allocatable :: AdvPolLMr(:,:)

  real(kind=8), public, allocatable :: CorPol2hInt(:)
  real(kind=8), public, allocatable :: CorTor2hInt(:)
  real(kind=8), public, allocatable :: CorPolAs2hInt(:)
  real(kind=8), public, allocatable :: CorTorAs2hInt(:)
  complex(kind=8), public, allocatable :: CorPolLMr(:,:)

  real(kind=8), public, allocatable :: LFPol2hInt(:)
  real(kind=8), public, allocatable :: LFTor2hInt(:)
  real(kind=8), public, allocatable :: LFPolAs2hInt(:)
  real(kind=8), public, allocatable :: LFTorAs2hInt(:)
  complex(kind=8), public, allocatable :: LFPolLMr(:,:)

  real(kind=8), public, allocatable :: DifPol2hInt(:,:)
  real(kind=8), public, allocatable :: DifTor2hInt(:,:)
  real(kind=8), public, allocatable :: DifPolAs2hInt(:,:)
  real(kind=8), public, allocatable :: DifTorAs2hInt(:,:)
  complex(kind=8), public, allocatable :: DifPolLMr(:,:)

  real(kind=8), public, allocatable :: Buo2hInt(:)
  real(kind=8), public, allocatable :: BuoAs2hInt(:)
  complex(kind=8), public, allocatable :: BuoLMr(:,:)

  real(kind=8), public, allocatable :: Pre2hInt(:)
  real(kind=8), public, allocatable :: PreAs2hInt(:)
  complex(kind=8), public, allocatable :: PreLMr(:,:)

  real(kind=8), public, allocatable :: Geo2hInt(:)
  real(kind=8), public, allocatable :: GeoAs2hInt(:)
  complex(kind=8), public, allocatable :: GeoLMr(:,:)

  real(kind=8), public, allocatable :: Mag2hInt(:)
  real(kind=8), public, allocatable :: MagAs2hInt(:)
  complex(kind=8), public, allocatable :: MagLMr(:,:)

  real(kind=8), public, allocatable :: Arc2hInt(:)
  real(kind=8), public, allocatable :: ArcAs2hInt(:)
  complex(kind=8), public, allocatable :: ArcLMr(:,:)

  public :: initialize_RMS, zeroRms

contains

   subroutine initialize_RMS

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

   end subroutine initialize_RMS
!----------------------------------------------------------------------------
   subroutine zeroRms
      !--------------------------------------------------------------------
      !  Zeros integrals that are set in s_get_td, s_update_z,
      !  s_update_wp, s_update_b, s_dtVrms and s_dtBrms
      !--------------------------------------------------------------------

      !-- Local variables
      integer :: n,nR,lm

      do n=1,1
         do nR=1,n_r_max
            dtVPol2hInt(nR,n)  =0.D0
            dtVTor2hInt(nR,n)  =0.D0
            dtVPolAs2hInt(nR,n)=0.D0
            dtVTorAs2hInt(nR,n)=0.D0
            DifPol2hInt(nR,n)  =0.D0
            DifTor2hInt(nR,n)  =0.D0
            DifPolAs2hInt(nR,n)=0.D0
            DifTorAs2hInt(nR,n)=0.D0
         end do
         do nR=1,n_r_maxMag
            dtBPol2hInt(nR,n)  =0.D0
            dtBTor2hInt(nR,n)  =0.D0
            dtBPolAs2hInt(nR,n)=0.D0
            dtBTorAs2hInt(nR,n)=0.D0
         end do
      end do

      do nR=1,n_r_max
         AdvPol2hInt(nR)  =0.D0
         AdvTor2hInt(nR)  =0.D0
         AdvPolAs2hInt(nR)=0.D0
         AdvTorAs2hInt(nR)=0.D0
         CorPol2hInt(nR)  =0.D0
         CorTor2hInt(nR)  =0.D0
         CorPolAs2hInt(nR)=0.D0
         CorTorAs2hInt(nR)=0.D0
         LFPol2hInt(nR)   =0.D0
         LFTor2hInt(nR)   =0.D0
         LFPolAs2hInt(nR) =0.D0
         LFTorAs2hInt(nR) =0.D0
         Buo2hInt(nR)     =0.D0
         BuoAs2hInt(nR)   =0.D0
         Pre2hInt(nR)     =0.D0
         PreAs2hInt(nR)   =0.D0
         Geo2hInt(nR)     =0.D0
         GeoAs2hInt(nR)   =0.D0
         Mag2hInt(nR)     =0.D0
         MagAs2hInt(nR)   =0.D0
         Arc2hInt(nR)     =0.D0
         ArcAs2hInt(nR)   =0.D0
      end do
      do nR=1,n_r_max
         do lm=1,lm_max
            dtVPolLMr(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
            AdvPolLMr(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
            CorPolLMr(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
            DifPolLMr(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
            BuoLMr(lm,nR)   =cmplx(0.D0,0.D0,kind=kind(0d0))
            PreLMr(lm,nR)   =cmplx(0.D0,0.D0,kind=kind(0d0))
            GeoLMr(lm,nR)   =cmplx(0.D0,0.D0,kind=kind(0d0))
            ArcLMr(lm,nR)   =cmplx(0.D0,0.D0,kind=kind(0d0))
         end do
      end do
      do nR=1,n_r_maxMag
         do lm=1,lm_maxMag
            dtBPolLMr(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
            LFPolLMr(lm,nR) =cmplx(0.D0,0.D0,kind=kind(0d0))
         end do
      end do

   end subroutine zeroRms
!----------------------------------------------------------------------------
end module RMS
