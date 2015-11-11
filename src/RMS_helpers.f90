module RMS_helpers
   !
   ! This module contains several useful subroutines required to compute RMS
   ! diagnostics
   !

   use precision_mod
   use truncation, only: l_max, lm_max_dtB, n_r_max, lm_max
   use blocking, only: lm2, st_map
   use radial_functions, only: or2, drx, chebt_oc, r
   use horizontal_data, only: osn1, Plm, dPlm, dLh
   use useful, only: cc2real
   use integration, only: rInt_R
   use LMmapping, only: mappings
   use constants, only: vol_oc, one

   implicit none

   private

   public :: get_PASLM, get_PolTorRms, hInt2dPol, hInt2Pol, hInt2Tor, &
             get_RAS, hIntRms, hInt2PolLM, hInt2TorLM

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
      PolRms  =rInt_R(PolRms_r,n_r_max,n_r_max,drx,chebt_oc)
      TorRms  =rInt_R(TorRms_r,n_r_max,n_r_max,drx,chebt_oc)
      PolAsRms=rInt_R(PolAsRms_r,n_r_max,n_r_max,drx,chebt_oc)
      TorAsRms=rInt_R(TorAsRms_r,n_r_max,n_r_max,drx,chebt_oc)
      fac=one/vol_oc
      PolRms  =sqrt(fac*PolRms)
      TorRms  =sqrt(fac*TorRms)
      PolAsRms=sqrt(fac*PolAsRms)
      TorAsRms=sqrt(fac*TorAsRms)

   end subroutine get_PolTorRms
!-----------------------------------------------------------------------------
   subroutine hInt2dPol(dPol,lmStart,lmStop,Pol2hInt,map)
    
      !-- Input variables:
      complex(cp),     intent(in) :: dPol(lm_max)   ! Toroidal field Potential
      integer,         intent(in) :: lmStart,lmStop
      type(mappings),  intent(in) :: map
    
      !-- Output variables:
      real(cp), intent(inout) :: Pol2hInt(0:l_max)
    
      !-- Local variables:
      real(cp) :: help
      integer :: lm,l,m
    
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=dLh(st_map%lm2(l,m))*cc2real(dPol(lm),m)
         Pol2hInt(l)=Pol2hInt(l)+help
      end do

   end subroutine hInt2dPol
!-----------------------------------------------------------------------------
   subroutine hInt2dPolLM(dPol,lmStart,lmStop,Pol2hInt,map)

      !-- Input variables
      complex(cp),     intent(in) :: dPol(lm_max) 
      integer,         intent(in) :: lmStart,lmStop
      type(mappings),  intent(in) :: map

      !-- Output variables
      real(cp), intent(inout) :: Pol2hInt(lm_max)

      !-- Local variables
      real(cp) :: help
      integer :: lm,l,m

      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=dLh(st_map%lm2(l,m))*cc2real(dPol(lm),m)
         Pol2hInt(lm)=Pol2hInt(lm)+help
      end do

   end subroutine hInt2dPolLM
!-----------------------------------------------------------------------------
   subroutine hInt2Pol(Pol,lb,ub,nR,lmStart,lmStop,PolLMr,Pol2hInt,map)

      !-- Input variables:
      integer,         intent(in) :: lb,ub
      complex(cp),     intent(in) :: Pol(lb:ub)   ! Poloidal field Potential
      integer,         intent(in) :: nR,lmStart,lmStop
      type(mappings),  intent(in) :: map

      !-- Output variables:
      complex(cp), intent(out) :: PolLMr(lm_max)
      real(cp),    intent(inout) :: Pol2hInt(0:l_max)

      !-- Local variables:
      real(cp) :: help,rE2
      integer :: lm,l,m

      rE2=r(nR)*r(nR)
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=rE2*cc2real(Pol(lm),m)
         Pol2hInt(l)=Pol2hInt(l)+help
         PolLMr(lm)=rE2/dLh(st_map%lm2(l,m))*Pol(lm)
      end do

   end subroutine hInt2Pol
!-----------------------------------------------------------------------------
   subroutine hInt2PolLM(Pol,lb,ub,nR,lmStart,lmStop,PolLMr,Pol2hInt,map)

      !-- Input variables:
      integer,         intent(in) :: lb,ub
      complex(cp),     intent(in) :: Pol(lb:ub)   ! Poloidal field Potential
      integer,         intent(in) :: nR,lmStart,lmStop
      type(mappings),  intent(in) :: map

      !-- Output variables:
      complex(cp), intent(out) :: PolLMr(lm_max)
      real(cp),    intent(inout) :: Pol2hInt(lm_max)

      !-- Local variables:
      real(cp) :: help,rE2
      integer :: lm,l,m

      rE2=r(nR)*r(nR)
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=rE2*cc2real(Pol(lm),m)
         Pol2hInt(lm)=Pol2hInt(lm)+help
         PolLMr(lm)=rE2/dLh(st_map%lm2(l,m))*Pol(lm)
      end do

   end subroutine hInt2PolLM
!-----------------------------------------------------------------------------
   subroutine hIntRms(f,nR,lmStart,lmStop,lmP,f2hInt,map)

      !-- Input variables
      complex(cp),    intent(in) :: f(*)
      integer,        intent(in) :: nR, lmStart, lmStop, lmP
      type(mappings), intent(in) :: map

      !-- Output variables
      real(cp),       intent(inout) :: f2hInt(0:l_max)

      !-- Local variables
      real(cp) :: help,rE2
      integer :: lm,l,m

      rE2=r(nR)*r(nR)
      do lm=lmStart,lmStop
         if ( lmP == 0 ) then
            l=map%lm2l(lm)
            m=map%lm2m(lm)
         else
            l=map%lmP2l(lm)
            m=map%lmP2m(lm)
         end if

         if ( l <= l_max ) then
            help=rE2*cc2real(f(lm),m)
            f2hInt(l)=f2hInt(l)+help
         end if
      end do

   end subroutine hIntRms
!-----------------------------------------------------------------------------
   subroutine hInt2Tor(Tor,lb,ub,nR,lmStart,lmStop,Tor2hInt,map)

      !-- Input variables:
      integer,         intent(in) :: lb,ub
      complex(cp),     intent(in) :: Tor(lb:ub)   ! Toroidal field Potential
      integer,         intent(in) :: nR
      integer,         intent(in) :: lmStart,lmStop
      type(mappings),  intent(in) :: map
    
      !-- Output variables:
      real(cp),        intent(inout) :: Tor2hInt(0:l_max)
    
      !-- Local variables:
      real(cp) :: help,rE4
      integer :: lm,l,m
    
      rE4=r(nR)**4
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=rE4/dLh(st_map%lm2(l,m))*cc2real(Tor(lm),m)
         Tor2hInt(l)=Tor2hInt(l)+help
      end do
    
   end subroutine hInt2Tor
!-----------------------------------------------------------------------------
   subroutine hInt2TorLM(Tor,lb,ub,nR,lmStart,lmStop,Tor2hInt,map)

      !-- Input variables:
      integer,         intent(in) :: lb,ub
      complex(cp),     intent(in) :: Tor(lb:ub)   ! Toroidal field Potential
      integer,         intent(in) :: nR
      integer,         intent(in) :: lmStart,lmStop
      type(mappings),  intent(in) :: map
    
      !-- Output variables:
      real(cp),        intent(inout) :: Tor2hInt(lm_max)
    
      !-- Local variables:
      real(cp) :: help,rE4
      integer :: lm,l,m
    
      rE4=r(nR)**4
      do lm=lmStart,lmStop
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         help=rE4/dLh(st_map%lm2(l,m))*cc2real(Tor(lm),m)
         Tor2hInt(lm)=Tor2hInt(lm)+help
      end do
    
   end subroutine hInt2TorLM
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
end module RMS_helpers
