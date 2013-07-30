!$Id$
!------------------------------------------------------------------
!               COMMON block for RMS calculation
!------------------------------------------------------------------

MODULE RMS
  USE truncation
  USE blocking, only: nThreadsMax

  IMPLICIT NONE

  !----  Horizontal integrals of squared forces:

  REAL(kind=8),ALLOCATABLE :: dtBPol2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: dtBTor2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: dtBPolAs2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: dtBTorAs2hInt(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dtBPolLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: dtVPol2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: dtVTor2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: dtVPolAs2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: dtVTorAs2hInt(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dtVPolLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: AdvPol2hInt(:)
  REAL(kind=8),ALLOCATABLE :: AdvTor2hInt(:)
  REAL(kind=8),ALLOCATABLE :: AdvPolAs2hInt(:)
  REAL(kind=8),ALLOCATABLE :: AdvTorAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: AdvPolLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: CorPol2hInt(:)
  REAL(kind=8),ALLOCATABLE :: CorTor2hInt(:)
  REAL(kind=8),ALLOCATABLE :: CorPolAs2hInt(:)
  REAL(kind=8),ALLOCATABLE :: CorTorAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: CorPolLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: LFPol2hInt(:)
  REAL(kind=8),ALLOCATABLE :: LFTor2hInt(:)
  REAL(kind=8),ALLOCATABLE :: LFPolAs2hInt(:)
  REAL(kind=8),ALLOCATABLE :: LFTorAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: LFPolLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: DifPol2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: DifTor2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: DifPolAs2hInt(:,:)
  REAL(kind=8),ALLOCATABLE :: DifTorAs2hInt(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: DifPolLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: Buo2hInt(:)
  REAL(kind=8),ALLOCATABLE :: BuoAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: BuoLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: Pre2hInt(:)
  REAL(kind=8),ALLOCATABLE :: PreAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: PreLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: Geo2hInt(:)
  REAL(kind=8),ALLOCATABLE :: GeoAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: GeoLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: Mag2hInt(:)
  REAL(kind=8),ALLOCATABLE :: MagAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: MagLMr(:,:)

  REAL(kind=8),ALLOCATABLE :: Arc2hInt(:)
  REAL(kind=8),ALLOCATABLE :: ArcAs2hInt(:)
  COMPLEX(kind=8),ALLOCATABLE :: ArcLMr(:,:)

CONTAINS
  SUBROUTINE initialize_RMS

  ALLOCATE( dtBPol2hInt(n_r_maxMag,nThreadsMax) )
  ALLOCATE( dtBTor2hInt(n_r_maxMag,nThreadsMax) )
  ALLOCATE( dtBPolAs2hInt(n_r_maxMag,nThreadsMax) )
  ALLOCATE( dtBTorAs2hInt(n_r_maxMag,nThreadsMax) )
  ALLOCATE( dtBPolLMr(lm_maxMag,n_r_maxMag) )

  ALLOCATE( dtVPol2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( dtVTor2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( dtVPolAs2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( dtVTorAs2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( dtVPolLMr(lm_max,n_r_max) )

  ALLOCATE( AdvPol2hInt(n_r_max) )
  ALLOCATE( AdvTor2hInt(n_r_max) )
  ALLOCATE( AdvPolAs2hInt(n_r_max) )
  ALLOCATE( AdvTorAs2hInt(n_r_max) )
  ALLOCATE( AdvPolLMr(lm_max,n_r_max) )

  ALLOCATE( CorPol2hInt(n_r_max) )
  ALLOCATE( CorTor2hInt(n_r_max) )
  ALLOCATE( CorPolAs2hInt(n_r_max) )
  ALLOCATE( CorTorAs2hInt(n_r_max) )
  ALLOCATE( CorPolLMr(lm_max,n_r_max) )

  ALLOCATE( LFPol2hInt(n_r_max) )
  ALLOCATE( LFTor2hInt(n_r_max) )
  ALLOCATE( LFPolAs2hInt(n_r_max) )
  ALLOCATE( LFTorAs2hInt(n_r_max) )
  ALLOCATE( LFPolLMr(lm_max,n_r_maxMag) )

  ALLOCATE( DifPol2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( DifTor2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( DifPolAs2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( DifTorAs2hInt(n_r_max,nThreadsMax) )
  ALLOCATE( DifPolLMr(lm_max,n_r_max) )

  ALLOCATE( Buo2hInt(n_r_max) )
  ALLOCATE( BuoAs2hInt(n_r_max) )
  ALLOCATE( BuoLMr(lm_max,n_r_max) )

  ALLOCATE( Pre2hInt(n_r_max) )
  ALLOCATE( PreAs2hInt(n_r_max) )
  ALLOCATE( PreLMr(lm_max,n_r_max) )

  ALLOCATE( Geo2hInt(n_r_max) )
  ALLOCATE( GeoAs2hInt(n_r_max) )
  ALLOCATE( GeoLMr(lm_max,n_r_max) )

  ALLOCATE( Mag2hInt(n_r_max) )
  ALLOCATE( MagAs2hInt(n_r_max) )
  ALLOCATE( MagLMr(lm_max,n_r_max) )

  ALLOCATE( Arc2hInt(n_r_max) )
  ALLOCATE( ArcAs2hInt(n_r_max) )
  ALLOCATE( ArcLMr(lm_max,n_r_max) )

  END SUBROUTINE initialize_RMS
END MODULE RMS
