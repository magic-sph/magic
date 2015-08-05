!$Id$
subroutine lmAS2pt(alm,aij,nThetaStart,nThetaBlockSize)
   !  +-------------+----------------+------------------------------------+
   !  |                                                                   |
   !  | Spherical harmonic transform from alm(l) to aij(theta)            |
   !  | Done within the range [nThetaStart,n_thetaStart+nThetaBlockSize-1]|
   !  | only considering axisymmetric contributions.                      |
   !  | alm contains only m=0 coefficients                                |
   !  |                                                                   |
   !  +-------------------------------------------------------------------+

   use truncation
   use radial_functions
   use blocking
   use horizontal_data

   implicit none

   !-- Input variables:
   integer,      intent(in) :: nThetaStart     ! first theta to be treated
   integer,      intent(in) :: nThetaBlockSize !
         
   real(kind=8), intent(in) :: alm(*)      ! field in (l,m)-space

   !-- Output variables:
   real(kind=8), intent(out) :: aij(*)  ! field in (theta,phi)-space

   !-- Local variables:
   integer :: nTheta        ! last theta treated
   integer :: nThetaNHS     ! counter for theta in one HS
   integer :: nThetaN       ! counter for theta in NHS
   integer :: nThetaS       ! counter for theta in SHS
   integer :: nThetaBlock   ! counter for theta in block
   integer :: l,lm          ! degree/order
   real(kind=8) ::  sign

   nTheta=nThetaStart-1  ! last theta
    
   !-- Zero output field:
   do nThetaBlock=1,nThetaBlockSize
      aij(nThetaBlock)=0.D0
   end do

   !-- Transform l 2 theta:
   nThetaNHS=nTheta/2  ! last theta in one HS

   do nThetaN=1,nThetaBlockSize,2
      nThetaS  =nThetaN+1
      nThetaNHS=nThetaNHS+1

      sign=-1.D0
      do l=0,l_max
         lm=lm2(l,0)
         sign=-sign
         ! Northern hemisphere
         aij(nThetaN)=aij(nThetaN) +      alm(l+1)*Plm(lm,nThetaNHS)
         ! Southern hemisphere
         aij(nThetaS)=aij(nThetaS) + sign*alm(l+1)*Plm(lm,nThetaNHS)
      end do

   end do

end subroutine lmAS2pt
!---------------------------------------------------------------------------------
