!$Id$
!***********************************************************************
    SUBROUTINE legtfAS2(flm1,flm2,ft1,ft2, &
                        lmMax,nThetaStart,sizeThetaB)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
!  |  [grid to spectral] for 2 arrays                                  |
!  |  ancl1/2 (input) to flm1/2 (output)                               |
!  |  One call to this routine does part of the transform              |
!  |  by summation over theta points in on theta block:                |
!  |      n_theta_min,..,n_theta_min+n_theta_block-1                   |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE horizontal_data

    IMPLICIT NONE

!-- input:
! include 'truncation.f'
! include 'c_horizontal.f'

    INTEGER :: lmMax          ! Number of modes to be processed
    INTEGER :: nThetaStart    ! First no of theta on block
    INTEGER :: sizeThetaB     ! Size of theta block

    REAL(kind=8) :: ft1(*),ft2(*)

!-- output: transformed arrays anlc1,anlc2
    REAL(kind=8) :: flm1(*),flm2(*)

!-- local:
    INTEGER :: nThetaN,nThetaS
    INTEGER :: nThetaHS
    INTEGER :: nTheta1,nTheta2
    INTEGER :: nThetaC1,nThetaC2
    INTEGER :: nThetaMin
    INTEGER :: lm1,lm2

    REAL(kind=8) :: f1p(n_theta_max/2),f1m(n_theta_max/2)
    REAL(kind=8) :: f2p(n_theta_max/2),f2m(n_theta_max/2)

!---------------------------------------------------------------------

!-- Prepare arrays of sums and differences:

    nThetaHS=0
    DO nThetaN=1,sizeThetaB,2 ! thetas in NHS
        nThetaS =nThetaN+1         ! thetas in SHS
        nThetaHS=nThetaHS+1        ! thetas in one HS
        f1p(nThetaHS)=ft1(nThetaN)+ft1(nThetaS) ! Symm
        f1m(nThetaHS)=ft1(nThetaN)-ft1(nThetaS) ! ASymm
        f2p(nThetaHS)=ft2(nThetaN)+ft2(nThetaS)
        f2m(nThetaHS)=ft2(nThetaN)-ft2(nThetaS)
    END DO

!-- Start with first two thetas in first theta block:
!   This initalizes the flm1 and flm2
    IF ( nThetaStart == 1 ) THEN
        DO lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm1(lm2)=f1m(1)*wPlm(lm2,1)+f1m(2)*wPlm(lm2,2)
        END DO
        DO lm1=1,lmMax-1,2
            lm2=lm1+1
            flm2(lm1)=f2p(1)*wPlm(lm1,1)+f2p(2)*wPlm(lm1,2)
            flm2(lm2)=f2m(1)*wPlm(lm2,1)+f2m(2)*wPlm(lm2,2)
        END DO
        IF ( lm2 < lmMax ) THEN
            lm1=lmMax
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm2(lm1)=f2p(1)*wPlm(lm1,1)+f2p(2)*wPlm(lm1,2)
        END IF
         
        IF ( sizeThetaB <= 4 ) RETURN
         
    END IF

!-- Following values of n_theta=3,4,5... when n_theta_min=1
!   or all values of n_theta when n_theta_min > 0
     
    nThetaMin=1
    If ( nThetaStart == 1 ) nThetaMin=3 ! 2 already done

!-- Calculating position for all thetas:
!     (nThetaStart-1)/2 are previous blocks,
!     (nThetaMin-1) is what has been done already
    nThetaC1=(nThetaStart-1)/2+nThetaMin-2

!-- Loop over half of the thetas with step 2 unrolling:
    DO nTheta1=nThetaMin,sizeThetaB/2,2
        nTheta2=nTheta1+1
        nThetaC1=nThetaC1+2
        nThetaC2=nThetaC1+1

        DO lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm1(lm2)=flm1(lm2) + f1m(nTheta1)*wPlm(lm2,nThetaC1) + &
                                  f1m(nTheta2)*wPlm(lm2,nThetaC2)
        END DO
        DO lm1=1,lmMax-1,2
            lm2=lm1+1
            flm2(lm1)=flm2(lm1) + f2p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f2p(nTheta2)*wPlm(lm1,nThetaC2)
            flm2(lm2)=flm2(lm2) + f2m(nTheta1)*wPlm(lm2,nThetaC1) + &
                                  f2m(nTheta2)*wPlm(lm2,nThetaC2)
        END DO
        IF ( lm2 < lmMax ) THEN
            lm1=lmMax
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm2(lm1)=flm2(lm1) + f2p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f2p(nTheta2)*wPlm(lm1,nThetaC2)
        END IF
         
    END DO
     

    RETURN
    end SUBROUTINE legtfAS2

!------------------------------------------------------------------------
