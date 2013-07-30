!$Id$
!***********************************************************************
    SUBROUTINE lmAS2pt(alm,aij,nThetaStart,nThetaBlockSize)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  | Spherical harmonic transform from alm(l) to aij(theta)            |
!  | Done within the range [nThetaStart,n_thetaStart+nThetaBlockSize-1]|
!  | only considering axisymmetric contributions.                      |
!  | alm contains only m=0 coefficients                                |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data

    IMPLICIT NONE

    INTEGER :: nThetaStart     ! first theta to be treated
    INTEGER :: nThetaBlockSize !
          
    REAL(kind=8) :: alm(*)      ! field in (l,m)-space

!-- Output:
    REAL(kind=8) :: aij(*)  ! field in (theta,phi)-space

!-- Local:
    INTEGER :: nTheta        ! last theta treated
    INTEGER :: nThetaNHS     ! counter for theta in one HS
    INTEGER :: nThetaN       ! counter for theta in NHS
    INTEGER :: nThetaS       ! counter for theta in SHS
    INTEGER :: nThetaBlock   ! counter for theta in block
    INTEGER :: l,lm          ! degree/order
    REAL(kind=8) ::  sign

!-- End of declaration
!-----------------------------------------------------------------------


    nTheta=nThetaStart-1  ! last theta
     
!-- Zero output field:
    DO nThetaBlock=1,nThetaBlockSize
        aij(nThetaBlock)=0.D0
    END DO

!-- Transform l 2 theta:
    nThetaNHS=nTheta/2  ! last theta in one HS

    DO nThetaN=1,nThetaBlockSize,2
        nThetaS  =nThetaN+1
        nThetaNHS=nThetaNHS+1

        sign=-1.D0
        DO l=0,l_max
            lm=lm2(l,0)
            sign=-sign
            aij(nThetaN)=aij(nThetaN) +   & ! Northern hemisp
                         alm(l+1)*Plm(lm,nThetaNHS)
            aij(nThetaS)=aij(nThetaS) +   & ! Southern hemisp
                         sign*alm(l+1)*Plm(lm,nThetaNHS)
        END DO

    END DO


    RETURN
    end SUBROUTINE lmAS2pt

!---------------------------------------------------------------------------------
