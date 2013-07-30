!$Id$
!***********************************************************************
    SUBROUTINE lm2pt(alm,aij,rT,nThetaStart,lIC,lrComp)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  | Spherical harmonic transform from alm(l,m) to aij(phi,theta)      |
!  | Radial field components are calculated for lrComp=.true.          |
!  | Done within the range [n_theta_min,n_theta_min+n_theta_block-1]   |
!  | Used only for graphic output.                                     |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
#if (FFTLIB==JW)
    USE fft_JW
#elif (FFTLIB==MKL)
    use fft_MKL
#endif

    IMPLICIT NONE

    INTEGER :: nThetaStart ! first theta to be treated
    REAL(kind=8) ::  rT
    LOGICAL :: lIC            ! true for inner core, extra factor !
    LOGICAL :: lrComp         ! true for radial field components
          
    COMPLEX(kind=8) :: alm(*)      ! field in (l,m)-space

!-- output:
    COMPLEX(kind=8) :: aij(ncp,*)  ! field in (theta,phi)-space

!-- local:
    INTEGER :: nThetaN,nThetaS,nThetaHS
    INTEGER :: lm,mca,m,l     ! degree/order
    REAL(kind=8) :: sign
    COMPLEX(kind=8) :: cs1(lm_max) ! help array
    REAL(kind=8) :: O_r_E_2,rRatio

!-- end of declaration
!-----------------------------------------------------------------------

    O_r_E_2=1.d0/(rT*rT)
    rRatio=rT/r_icb
     
!-- Zero output field:
    DO nThetaN=1,sizeThetaB
        DO mca=1,ncp
            aij(mca,nThetaN)=0.D0
        END DO
    END DO
     
!-- Multiplication with l(l+1)/r**2 for radial component:
    cs1(1)=0.D0
    DO lm=2,lm_max
        cs1(lm)=alm(lm)
        IF ( lrComp ) cs1(lm)=cs1(lm)*O_r_E_2*dLh(lm)
        IF ( lIC )    cs1(lm)=rRatio**D_lP1(lm)*cs1(lm)
    END DO

!-- Transform l 2 theta:
    nThetaHS=(nThetaStart-1)/2  ! last theta in one HS
    do nThetaN=1,sizeThetaB,2
        nThetaS =nThetaN+1
        nThetaHS=nThetaHS+1
        lm=0
        mca=0
        DO m=0,m_max,minc
            mca=mca+1
            sign=-1.D0
            DO l=m,l_max
                lm=lm+1
                sign=-sign
                aij(mca,nThetaN)=aij(mca,nThetaN) + & ! Northern hemisp
                                 cs1(lm)*Plm(lm,nThetaHS)
                aij(mca,nThetaS)=aij(mca,nThetaS) + & ! Southern hemisp
                                 sign*cs1(lm)*Plm(lm,nThetaHS)
            END DO
        END DO
    END DO

!-- Transform m 2 phi:
    CALL fft_thetab(aij,1)

    RETURN
    end SUBROUTINE lm2pt

!---------------------------------------------------------------------------------
