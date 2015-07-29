!$Id$
!****************************************************************************
subroutine getAStr(fZ,flmn,nZmax,nZmaxA,lmMax,lMax, &
                       rMin,rMax,nChebMax,rZ,Plm)
   !----------------------------------------------------------------------------
   !  Calculates function value at radii rZ(nZmax) and
   !  colatitudes for which Plm(theta) is given from
   !  the spherical hamornic/Chebychev coefficients of an
   !  axisymmetric function (order=0).
   !----------------------------------------------------------------------------

   IMPLICIT NONE

   !--- Output variables :
   real(kind=8), intent(out) :: fZ(*)

   !--- Input variables:
   integer,      intent(in) :: lmMax,lMax
   real(kind=8), intent(in) :: flmn(lmMax,*)
   integer,      intent(in) :: nZmax,nZmaxA
   real(kind=8), intent(in) :: rMin,rMax
   integer,      intent(in) :: nChebMax
   real(kind=8), intent(in) :: rZ(nZmaxA/2+1)
   real(kind=8), intent(in) :: Plm(lmMax,nZmaxA/2+1)

   !--- Local variables:
   integer :: nCheb
   real(kind=8) :: cheb(nChebMax)
   integer :: l,nZS,nZN
   real(kind=8) :: x,chebNorm,flmr

   chebNorm=dsqrt(2.D0/dble(nChebMax-1))

   do nZN=1,nZmax
      fZ(nZN)=0.D0
   end do

   do nZN=1,nZmax/2 ! Loop over all z-points in south HS
      nZS=nZmax+1-nZN

      !--- Map r to cheb intervall [-1,1]:
      !    and calculate the cheb polynomia:
      !    Note: the factor chebNorm is needed
      !    for renormalisation. Its not needed if one used
      !    costf1 for the back transform.
      x=2.D0*(rZ(nZN)-0.5D0*(rMin+rMax))/(rMax-rMin)
      cheb(1)=1.D0*chebNorm
      cheb(2)=x*chebNorm
      do nCheb=3,nChebMax
         cheb(nCheb)=2.D0*x*cheb(nCheb-1)-cheb(nCheb-2)
      end do
      cheb(1)       =0.5D0*cheb(1)
      cheb(nChebMax)=0.5D0*cheb(nChebMax)

      !--- Loop to add all contribution functions:
      do l=0,lMax
         flmr=0.D0
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

      x=2.D0*(rZ(nZS)-0.5D0*(rMin+rMax))/(rMax-rMin)
      cheb(1)=1.D0*chebNorm
      cheb(2)=x*chebNorm
      do nCheb=3,nChebMax
         cheb(nCheb)=2.D0*x*cheb(nCheb-1)-cheb(nCheb-2)
      end do
      cheb(1)       =0.5D0*cheb(1)
      cheb(nChebMax)=0.5D0*cheb(nChebMax)

      do l=0,lMax
         flmr=0.D0
         do nCheb=1,nChebMax
            flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
         end do
         fZ(nZS)=fZ(nZS)+flmr*Plm(l+1,nZS)
      end do

   end if

end subroutine getAStr
!----------------------------------------------------------------------
