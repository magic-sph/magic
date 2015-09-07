!$Id$
module TO_helpers

   use precision_mod
   use truncation, only: l_max
   use blocking, only: lm2
   use horizontal_data, only: dPlm, osn1
   use const, only: one, two, half

   implicit none

   private

   public :: getPAStr, get_PAS, getAStr

contains

   subroutine getPAStr(fZ,flmn,nZmax,nZmaxA,lmMax,lMax, &
                              rMin,rMax,nChebMax,rZ,dPlm,OsinTS)
      !----------------------------------------------------------------------------
      !  Calculates axisymmetric phi component for a given toroidal
      !  potential flmn in spherical harmonic/Cheb space.
      !  This is calculated at radii rZ(nZmax) and matching
      !  colatitutes theta(nZmax) for which dPlm(theta) and OsinTS(theta)
      !  are provided for the south hemisphere.
      !  Points in the northern hemipshere use Plm symmetries.
      !----------------------------------------------------------------------------

      !--- Input variables:
      integer,  intent(in) :: lmMax, lMax
      real(cp), intent(in) :: flmn(lmMax,*)
      integer,  intent(in) :: nZmax, nZmaxA
      real(cp), intent(in) :: rMin, rMax
      integer,  intent(in) :: nChebMax
      real(cp), intent(in) :: rZ(nZmaxA/2+1)
      real(cp), intent(in) :: dPlm(lmMax,nZmaxA/2+1)
      real(cp), intent(in) :: OsinTS(nZmaxA/2+1)

      !--- Output variables:
      real(cp), intent(out) ::    fZ(*)

      !--- Local variables:
      integer :: nCheb
      real(cp) :: cheb(nChebMax)
      integer :: l,nZS,nZN!,nZ
      real(cp) :: x,chebNorm,fac,flmr

      chebNorm=sqrt(two/real(nChebMax-1, kind=cp))
              
      do nZN=1,nZmax
         fZ(nZN)=0.0_cp
      end do

      do nZN=1,nZmax/2 ! Loop over all z-points in northern HS
         nZS=nZmax+1-nZN  ! point in southern HS
         fac=-OsinTS(nZN)/rZ(nZN)

         !--- Map r to cheb intervall [-1,1]:
         !    and calculate the cheb polynomia:
         !    Note: the factor chebNorm is needed
         !    for renormalisation. Its not needed if one used
         !    costf1 for the back transform.
         x=two*(rZ(nZN)-half*(rMin+rMax))/(rMax-rMin)
         cheb(1)=one*chebNorm*fac
         cheb(2)=x*chebNorm*fac
         do nCheb=3,nChebMax
            cheb(nCheb)=two*x*cheb(nCheb-1)-cheb(nCheb-2)
         end do
         cheb(1)       =half*cheb(1)
         cheb(nChebMax)=half*cheb(nChebMax)

         !--- Loop to add all contribution functions:
         do l=0,lMax
            flmr=0.0_cp
            do nCheb=1,nChebMax
               flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
            end do
            fZ(nZN)=fZ(nZN)+flmr*dPlm(l+1,nZN)
            if ( mod(l,2) == 0 ) then ! Odd contribution
               fZ(nZS)=fZ(nZS)-flmr*dPlm(l+1,nZN)
            else
               fZ(nZS)=fZ(nZS)+flmr*dPlm(l+1,nZN)
            end if
         end do

      end do

      if ( mod(nZmax,2) == 1 ) then ! Remaining equatorial point
         nZS=(nZmax-1)/2+1
         fac=-OsinTS(nZS)/rZ(nZS)

         x=two*(rZ(nZS)-half*(rMin+rMax))/(rMax-rMin)
         cheb(1)=one*chebNorm*fac
         cheb(2)=x*chebNorm*fac
         do nCheb=3,nChebMax
            cheb(nCheb)=two*x*cheb(nCheb-1)-cheb(nCheb-2)
         end do
         cheb(1)       =half*cheb(1)
         cheb(nChebMax)=half*cheb(nChebMax)

         do l=0,lMax
            flmr=0.0_cp
            do nCheb=1,nChebMax
               flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
            end do
            fZ(nZS)=fZ(nZS)+flmr*dPlm(l+1,nZS)
         end do

      end if

   end subroutine getPAStr
!---------------------------------------------------------------------------
   subroutine get_PAS(Tlm,Bp,rT,nThetaStart,sizeThetaBlock)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to calculate the axisymmetric      |
      !  |  phi component Bp of an axisymmetric toroidal field Tlm           |
      !  |  given in spherical harmonic space (l,m=0).                       |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables
      integer,  intent(in) :: nThetaStart    ! first theta to be treated
      integer,  intent(in) :: sizeThetaBlock ! size of theta block
      real(cp), intent(in) :: rT             ! radius
      real(cp), intent(in) :: Tlm(*)         ! field in (l,m)-space for rT

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
            Bp_1=-Tlm(l+1)*dPlm(lm,nThetaN)
            Bp_n=Bp_n+Bp_1
            Bp_s=Bp_s-sign*Bp_1
         end do  ! Loop over degree
         Bp(nTheta)  =fac*Bp_n
         Bp(nTheta+1)=fac*Bp_s

      end do        ! Loop over colatitudes

   end subroutine get_PAS
!----------------------------------------------------------------------------
   subroutine getAStr(fZ,flmn,nZmax,nZmaxA,lmMax,lMax,rMin,rMax,nChebMax,rZ,Plm)
      !----------------------------------------------------------------------------
      !  Calculates function value at radii rZ(nZmax) and
      !  colatitudes for which Plm(theta) is given from
      !  the spherical harmonic/Chebychev coefficients of an
      !  axisymmetric function (order=0).
      !----------------------------------------------------------------------------

      !--- Input variables:
      integer,  intent(in) :: lmMax, lMax
      real(cp), intent(in) :: flmn(lmMax,*)
      integer,  intent(in) :: nZmax, nZmaxA
      real(cp), intent(in) :: rMin, rMax
      integer,  intent(in) :: nChebMax
      real(cp), intent(in) :: rZ(nZmaxA/2+1)
      real(cp), intent(in) :: Plm(lmMax,nZmaxA/2+1)

      !--- Output variables :
      real(cp), intent(out) :: fZ(*)


      !--- Local variables:
      integer :: nCheb
      real(cp) :: cheb(nChebMax)
      integer :: l,nZS,nZN
      real(cp) :: x,chebNorm,flmr

      chebNorm=sqrt(two/real(nChebMax-1, kind=cp))

      do nZN=1,nZmax
         fZ(nZN)=0.0_cp
      end do

      do nZN=1,nZmax/2 ! Loop over all z-points in south HS
         nZS=nZmax+1-nZN

         !--- Map r to cheb intervall [-1,1]:
         !    and calculate the cheb polynomia:
         !    Note: the factor chebNorm is needed
         !    for renormalisation. Its not needed if one used
         !    costf1 for the back transform.
         x=two*(rZ(nZN)-half*(rMin+rMax))/(rMax-rMin)
         cheb(1)=one*chebNorm
         cheb(2)=x*chebNorm
         do nCheb=3,nChebMax
            cheb(nCheb)=two*x*cheb(nCheb-1)-cheb(nCheb-2)
         end do
         cheb(1)       =half*cheb(1)
         cheb(nChebMax)=half*cheb(nChebMax)

         !--- Loop to add all contribution functions:
         do l=0,lMax
            flmr=0.0_cp
            do nCheb=1,nChebMax
               flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
            end do
            fZ(nZN)=fZ(nZN)+flmr*Plm(l+1,nZN)
            if ( mod(l,2) == 0 ) then ! Even contribution
               fZ(nZS)=fZ(nZS)+flmr*Plm(l+1,nZN)
            else
               fZ(nZS)=fZ(nZS)-flmr*Plm(l+1,nZN)
            end if
         end do

      end do

      if ( mod(nZmax,2) == 1 ) then ! Remaining equatorial point
         nZS=(nZmax-1)/2+1

         x=two*(rZ(nZS)-half*(rMin+rMax))/(rMax-rMin)
         cheb(1)=one*chebNorm
         cheb(2)=x*chebNorm
         do nCheb=3,nChebMax
            cheb(nCheb)=two*x*cheb(nCheb-1)-cheb(nCheb-2)
         end do
         cheb(1)       =half*cheb(1)
         cheb(nChebMax)=half*cheb(nChebMax)

         do l=0,lMax
            flmr=0.0_cp
            do nCheb=1,nChebMax
               flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
            end do
            fZ(nZS)=fZ(nZS)+flmr*Plm(l+1,nZS)
         end do

      end if

   end subroutine getAStr
!----------------------------------------------------------------------------
end module TO_helpers
