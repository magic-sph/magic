module RMS
   !
   ! This module contains the global array used when RMS force balance
   ! is requested
   !

   use precision_mod
   use truncation, only: n_r_max, n_cheb_max, n_r_maxMag, lm_max, lm_maxMag, &
                         l_max
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
   real(cp), public, allocatable :: CLF2hInt(:,:)
   real(cp), public, allocatable :: PLF2hInt(:,:)
 
   public :: initialize_RMS, zeroRms

contains

   subroutine initialize_RMS
      !
      ! Memory allocation
      !

      integer, parameter :: nThreadsMax=1

      allocate( dtBPol2hInt(lm_maxMag,n_r_maxMag,nThreadsMax) )
      allocate( dtBTor2hInt(lm_maxMag,n_r_maxMag,nThreadsMax) )
      allocate( dtBPolLMr(lm_maxMag,n_r_maxMag) )
    
      allocate( dtVPol2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( dtVTor2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( dtVPolLMr(lm_max,n_r_max) )

      allocate( DifPol2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( DifTor2hInt(0:l_max,n_r_max,nThreadsMax) )
      allocate( DifPolLMr(lm_max,n_r_max) )
    
      allocate( Adv2hInt(0:l_max,n_r_max) )
      allocate( Cor2hInt(0:l_max,n_r_max) )
      allocate( LF2hInt(0:l_max,n_r_max) )
      allocate( Buo2hInt(0:l_max,n_r_max) )
      allocate( Pre2hInt(0:l_max,n_r_max) )
      allocate( Geo2hInt(0:l_max,n_r_max) )
      allocate( Mag2hInt(0:l_max,n_r_max) )
      allocate( Arc2hInt(0:l_max,n_r_max) )
      allocate( CLF2hInt(0:l_max,n_r_max) )
      allocate( PLF2hInt(0:l_max,n_r_max) )

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
            do lm=1,lm_maxMag
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
            CLF2hInt(l,nR)  =0.0_cp
            PLF2hInt(l,nR)  =0.0_cp
         end do
      end do
      do nR=1,n_r_max
         do lm=1,lm_max
            dtVPolLMr(lm,nR)=zero
            DifPolLMr(lm,nR)=zero
         end do
      end do
      do nR=1,n_r_maxMag
         do lm=1,lm_maxMag
            dtBPolLMr(lm,nR)=zero
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
      integer :: nRs(256)
    
      logical :: lStop
      integer :: n,i,nR
    
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
    
      ! Allowed number of radial grid points:
      nRs = [([4*i+1],i=1,256)]
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
