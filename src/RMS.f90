module RMS
   !
   ! This module contains the global array used when RMS force balance
   ! is requested
   !

   use const, only: zero
   use precision_mod
   use truncation, only: n_r_max, n_r_maxMag, lm_max, lm_maxMag
 
   implicit none
 
   private
 
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
end module RMS
