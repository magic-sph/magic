module RMS_helpers
   !
   ! This module contains several useful subroutines required to compute RMS
   ! diagnostics
   !

   use precision_mod
   use truncation, only: l_max, lm_max_dtB, n_r_max, lm_max
   use blocking, only: lm2, st_map
   use radial_functions, only: or2, drx, i_costf_init, d_costf_init, &
                               r
   use horizontal_data, only: osn1, Plm, dPlm, dLh
   use useful, only: cc2real
   use integration, only: rInt_R
   use LMmapping, only: mappings
   use const, only: vol_oc, one
   use chebyshev_polynoms_mod, only: cheb_grid
   use init_costf, only: init_costf1
   use radial_der, only: get_dr

   implicit none

   private

   public :: get_PASLM, get_PolTorRms, hInt2dPol, hInt2Pol, hInt2Tor, &
             init_rNB, get_RAS

contains

   subroutine get_PASLM(Tlm,Bp,rT,nThetaStart,sizeThetaBlock)
      !
      !  Purpose of this subroutine is to calculated the axisymmetric     
      !  phi component Bp of an axisymmetric toroidal field Tlm           
      !  given in spherical harmonic space (l,m=0).                       
      !

      !-- Input variables:
      integer,     intent(in) :: nThetaStart    ! first theta to be treated
      integer,     intent(in) :: sizeThetaBlock ! size of theta block
      real(cp),    intent(in) :: rT             ! radius
      complex(cp), intent(in) :: Tlm(lm_max_dtB) ! field in (l,m)-space for rT

      !-- Output variables:
      real(cp), intent(out) :: Bp(*)

      !-- Local variables:
      integer :: lm,l
      integer :: nTheta,nThetaN
      real(cp) :: fac
      real(cp) :: sign
      real(cp) :: Bp_1,Bp_n,Bp_s

      do nTheta=1,sizeThetaBlock,2 ! loop over thetas in northers HS

         nThetaN=(nThetaStart+nTheta)/2
         fac=osn1(nThetaN)/rT

         sign=-one
         Bp_n=0.0_cp
         Bp_s=0.0_cp
         do l=0,l_max
            lm=lm2(l,0)
            sign=-sign
            Bp_1=-real(Tlm(l+1))*dPlm(lm,nThetaN)
            Bp_n=Bp_n+Bp_1
            Bp_s=Bp_s-sign*Bp_1
         end do  ! Loop over degree
         Bp(nTheta)  =fac*Bp_n
         Bp(nTheta+1)=fac*Bp_s

      end do        ! Loop over colatitudes

   end subroutine get_PASLM
!---------------------------------------------------------------------------
   subroutine get_PolTorRms(Pol,drPol,Tor,PolRms,TorRms,PolAsRms,TorAsRms,map)
      !
      !  calculates integral PolRms=sqrt( Integral (pol^2 dV) )
      !  calculates integral TorRms=sqrt( Integral (tor^2 dV) )
      !  plus axisymmetric parts.
      !  integration in theta,phi by summation of spherical harmonics
      !  integration in r by using Chebycheff integrals
      !  The mapping map gives the mapping lm to l,m for the input
      !  arrays Pol,drPol and Tor
      !  Output: PolRms,TorRms,PolAsRms,TorAsRms
      !
    
      !-- Input variables:
      complex(cp),     intent(in) :: Pol(lm_max,n_r_max)   ! Poloidal field Potential
      complex(cp),     intent(in) :: drPol(lm_max,n_r_max) ! Radial derivative of Pol
      complex(cp),     intent(in) :: Tor(lm_max,n_r_max)   ! Toroidal field Potential
      type(mappings),  intent(in) :: map
    
      !-- Output variables:
      real(cp), intent(out) :: PolRms,PolAsRms
      real(cp), intent(out) :: TorRms,TorAsRms
    
      !-- Local variables:
      real(cp) :: PolRmsTemp,TorRmsTemp
      real(cp) :: PolRms_r(n_r_max)
      real(cp) :: TorRms_r(n_r_max)
      real(cp) :: PolAsRms_r(n_r_max)
      real(cp) :: TorAsRms_r(n_r_max)
    
      integer :: n_r,lm,l,m
      real(cp) :: fac
    
      do n_r=1,n_r_max
    
         PolRms_r(n_r)  =0.0_cp
         TorRms_r(n_r)  =0.0_cp
         PolAsRms_r(n_r)=0.0_cp
         TorAsRms_r(n_r)=0.0_cp
    
         do lm=2,lm_max
            l=map%lm2l(lm)
            m=map%lm2m(lm)
            PolRmsTemp= dLh(st_map%lm2(l,m)) * (                        &
                 dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(Pol(lm,n_r),m) + &
                 cc2real(drPol(lm,n_r),m) )
            TorRmsTemp=   dLh(st_map%lm2(l,m))*cc2real(Tor(lm,n_r),m)
            if ( m == 0 ) then  ! axisymmetric part
               PolAsRms_r(n_r)=PolAsRms_r(n_r) + PolRmsTemp
               TorAsRms_r(n_r)=TorAsRms_r(n_r) + TorRmsTemp
            else
               PolRms_r(n_r)  =PolRms_r(n_r)   + PolRmsTemp
               TorRms_r(n_r)  =TorRms_r(n_r)   + TorRmsTemp
            end if
         end do    ! do loop over lms in block
         PolRms_r(n_r)=PolRms_r(n_r) + PolAsRms_r(n_r)
         TorRms_r(n_r)=TorRms_r(n_r) + TorAsRms_r(n_r)
      end do    ! radial grid points
    
      !-- Radial Integrals:
      PolRms  =rInt_R(PolRms_r,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
      TorRms  =rInt_R(TorRms_r,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
      PolAsRms=rInt_R(PolAsRms_r,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
      TorAsRms=rInt_R(TorAsRms_r,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
      fac=one/vol_oc
      PolRms  =sqrt(fac*PolRms)
      TorRms  =sqrt(fac*TorRms)
      PolAsRms=sqrt(fac*PolAsRms)
      TorAsRms=sqrt(fac*TorAsRms)

   end subroutine get_PolTorRms
!-----------------------------------------------------------------------------
   subroutine hInt2dPol(dPol,lmStart,lmStop,Pol2hInt,PolAs2hInt,map)
    
      !-- Input variables:
      complex(cp),     intent(in) :: dPol(lm_max)   ! Toroidal field Potential
      integer,         intent(in) :: lmStart,lmStop
      type(mappings),  intent(in) :: map
    
      !-- Output variables:
      real(cp), intent(inout) :: Pol2hInt,PolAs2hInt
    
      !-- Local variables:
      real(cp) :: help
      integer :: lm,l,m
    
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=dLh(st_map%lm2(l,m))*cc2real(dPol(lm),m)
         if ( m == 0 ) PolAs2hInt=PolAs2hInt+help
         Pol2hInt=Pol2hInt+help
      end do

   end subroutine hInt2dPol
!-----------------------------------------------------------------------------
   subroutine hInt2Pol(Pol,lb,ub,nR,lmStart,lmStop,PolLMr, &
                       Pol2hInt,PolAs2hInt,map)

      !-- Input variables:
      integer,         intent(in) :: lb,ub
      complex(cp),     intent(in) :: Pol(lb:ub)   ! Poloidal field Potential
      integer,         intent(in) :: nR,lmStart,lmStop
      type(mappings),  intent(in) :: map

      !-- Output variables:
      complex(cp), intent(out) :: PolLMr(lm_max,n_r_max)
      real(cp),    intent(inout) :: Pol2hInt,PolAs2hInt

      !-- Local variables:
      real(cp) :: help,rE2
      integer :: lm,l,m

      rE2=r(nR)*r(nR)
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=rE2*cc2real(Pol(lm),m)
         if ( m == 0 ) PolAs2hInt=PolAs2hInt+help
         Pol2hInt=Pol2hInt+help
         PolLMr(lm,nR)=rE2/dLh(st_map%lm2(l,m))*Pol(lm)
      end do

   end subroutine hInt2Pol
!-----------------------------------------------------------------------------
   subroutine hInt2Tor(Tor,lb,ub,nR,lmStart,lmStop,Tor2hInt,TorAs2hInt,map)

      !-- Input variables:
      integer,         intent(in) :: lb,ub
      complex(cp),     intent(in) :: Tor(lb:ub)   ! Toroidal field Potential
      integer,         intent(in) :: nR
      integer,         intent(in) :: lmStart,lmStop
      type(mappings),  intent(in) :: map
    
      !-- Output variables:
      real(cp),        intent(inout) :: Tor2hInt,TorAs2hInt
    
      !-- Local variables:
      real(cp) :: help,rE4
      integer :: lm,l,m
    
      rE4=r(nR)**4
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=rE4/dLh(st_map%lm2(l,m))*cc2real(Tor(lm),m)
         if ( m == 0 ) TorAs2hInt=TorAs2hInt+help
         Tor2hInt=Tor2hInt+help
      end do
    
   end subroutine hInt2Tor
!-----------------------------------------------------------------------------
   subroutine get_RAS(Blm,Br,rT,nThetaStart,sizeThetaBlock)
      !
      !  Purpose of this subroutine is to calculate the axisymmetric      
      !  radial component Br of an axisymmetric ploidal field Blm         
      !  given in spherical harmonic space (l,m=0).                       
      !

      !-- Input variables:
      integer,         intent(in) :: nThetaStart    ! first theta to be treated
      integer,         intent(in) :: sizeThetaBlock ! last theta
      real(cp),        intent(in) :: rT             ! radius
      complex(cp),     intent(in) :: Blm(lm_max_dtB)! field in (l,m)-space for rT

      !-- Output variables:
      real(cp), intent(out) :: Br(*)

      !-- Local variables:
      integer :: lm,l
      integer :: nTheta,nThetaN
      real(cp) :: fac
      real(cp) :: sign
      real(cp) :: Br_1,Br_n,Br_s

      fac=one/(rT*rT)

      do nTheta=1,sizeThetaBlock,2 ! loop over thetas in northers HS
         nThetaN=(nThetaStart+nTheta)/2

         sign=-one
         Br_n=0.0_cp
         Br_s=0.0_cp
         do l=0,l_max
            lm=lm2(l,0)
            sign=-sign
            Br_1=real(Blm(l+1))*real(l*(l+1),kind=cp)*Plm(lm,nThetaN)
            Br_n=Br_n+Br_1
            Br_s=Br_s+sign*Br_1
         end do  ! Loop over degree
         Br(nTheta)  =fac*Br_n
         Br(nTheta+1)=fac*Br_s

      end do        ! Loop over colatitudes

    end subroutine get_RAS
!-----------------------------------------------------------------------------
   subroutine init_rNB(r,n_r_max,n_cheb_max,rCut,rDea,      &
        &              r2,n_r_max2,n_cheb_max2,             &
        &              nS,dr_fac2,i_costf_init2,nDi_costf1, &
        &              d_costf_init2,nDd_costf1)
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
      real(cp), intent(in) :: r(*),rCut,rDea
      integer,  intent(in) :: n_r_max,n_cheb_max
      integer,  intent(in) :: nDi_costf1,nDd_costf1
    
      !--- Output variables:
      integer,  intent(out) :: nS,n_r_max2,n_cheb_max2
      real(cp), intent(out) :: r2(*),dr_fac2(*)
      integer,  intent(out) :: i_costf_init2(nDi_costf1)   ! info for transform
      real(cp), intent(out) :: d_costf_init2(nDd_costf1)   ! info for tranfor
    
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
    
      do nR=1,n_r_max2
         r2(nR)=r(nR+nS)
      end do
      r_icb2=r2(n_r_max2)
      r_cmb2=r2(1)
      call cheb_grid(r_icb2,r_cmb2,n_r_max2-1,r2C,r_cheb2, &
                     0.0_cp,0.0_cp,0.0_cp,0.0_cp)
      call init_costf1(n_r_max2,i_costf_init2,nDi_costf1, &
                       d_costf_init2,nDd_costf1)
      dr_fac=one
      do nR=1,n_r_max
         drx(nR)=one
      end do
      call get_dr(r2,dr2,n_r_max2,n_cheb_max2, &
                  w1,w2,i_costf_init2,d_costf_init2,drx)
      do nR=1,n_r_max2
         dr_fac2(nR)=one/dr2(nR)
      end do

   end subroutine init_rNB
!-----------------------------------------------------------------------------
end module RMS_helpers
