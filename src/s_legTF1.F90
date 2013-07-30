!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"

    SUBROUTINE legTF1(nThetaStart,f1LM,f1TM)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
!  |  [grid to spectral] for 2 arrays                                  |
!  |  f1TM (input) to f1LM (output)                                    |
!  |  One call to this routine does part of the transform              |
!  |  by summation over theta points in one theta block:               |
!  |      nThetaStart,..,nThetaStart+n_theta_block-1                   |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE blocking
    USE horizontal_data

    IMPLICIT NONE

!-- input:
! include 'truncation.f'
! include 'c_horizontal.f'
! include 'c_blocking.f'

    INTEGER,INTENT(IN) :: nThetaStart ! First no of theta on block
    COMPLEX(kind=8),INTENT(IN) :: f1TM(ncp,nfs)

!-- output:
    COMPLEX(kind=8),INTENT(OUT) :: f1LM(lmP_max)

!-- local:
    INTEGER :: nThetaN     ! No. of theta in NHS
    INTEGER :: nThetaS     ! No. of theta in SHS
    INTEGER :: nThetaNHS   ! No. of thetas in one HS only
    INTEGER :: nThetaB1    ! No. of theta in block
    INTEGER :: nThetaB2    ! No. of theta in block
    INTEGER :: nTheta1     ! No. of theta (in one HS)
    INTEGER :: nTheta2     ! No. of theta (in one HS)
    INTEGER :: nThetaMin   ! Where to start in block

    INTEGER :: mc          ! counter of spherical order
    INTEGER :: lmS,lm      ! counter of spherical mode

! NTEGER nfs2
! ARAMETER (nfs2=nfs/2)
    COMPLEX(kind=8) :: f1ES(n_m_max,nfs/2),f1ES1,f1ES2
    COMPLEX(kind=8) :: f1EA(n_m_max,nfs/2),f1EA1,f1EA2

!-- end of declaration
!---------------------------------------------------------------------

    PERFON('legTF1')
!-- Unscrambles equatorially symmetric and antisymmetric contributions:
    nThetaNHS=0
    DO nThetaN=1,sizeThetaB,2 ! thetas in NHS
        nThetaS=nThetaN+1      ! thetas in SHS
        nThetaNHS=nThetaNHS+1  ! thetas counted in NHS only
        DO mc=1,n_m_max        ! counts spherical harmonic orders
            f1ES(mc,nThetaNHS)=f1TM(mc,nThetaN)+f1TM(mc,nThetaS)
            f1EA(mc,nThetaNHS)=f1TM(mc,nThetaN)-f1TM(mc,nThetaS)
        END DO
    END DO

!- Start with first two thetas for first theta block:

    IF ( nThetaStart == 1 ) THEN

        DO mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,1)
            f1ES2=f1ES(mc,2)
            f1EA1=f1EA(mc,1)
            f1EA2=f1EA(mc,2)
            DO lm=lStartP(mc),lmS-1,2
                f1LM(lm)  =f1ES1*wPlm(lm,1) +f1ES2*wPlm(lm,2)
                f1LM(lm+1)=f1EA1*wPlm(lm+1,1)+f1EA2*wPlm(lm+1,2)
            END DO
            IF ( lmOddP(mc) ) THEN
                f1LM(lmS) =f1ES1*wPlm(lmS,1)+f1ES2*wPlm(lmS,2)
            END IF
        END DO
         
        IF ( sizeThetaB <= 4 ) goto 500 ! return

    END IF

     
!-- Loop over half of the thetas with step 2 unrolling:
    nThetaMin=1
    IF ( nThetaStart == 1 ) nThetaMin=3
    nTheta1=(nThetaStart-1)/2+nThetaMin-2 ! NHS thetas covered before
    DO nThetaB1=nThetaMin,sizeThetaB/2,2
        nThetaB2=nThetaB1+1
        nTheta1 =nTheta1+2
        nTheta2 =nTheta1+1

        DO mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,nThetaB1)
            f1ES2=f1ES(mc,nThetaB2)
            f1EA1=f1EA(mc,nThetaB1)
            f1EA2=f1EA(mc,nThetaB2)
            DO lm=lStartP(mc),lmS-1,2
                f1LM(lm)  =f1LM(lm)   + f1ES1*wPlm(lm,nTheta1) + &
                                        f1ES2*wPlm(lm,nTheta2)
                f1LM(lm+1)=f1LM(lm+1) + f1EA1*wPlm(lm+1,nTheta1) + &
                                        f1EA2*wPlm(lm+1,nTheta2)
            END DO
            IF ( lmOddP(mc) ) THEN
                f1LM(lmS)=f1LM(lmS) + f1ES1*wPlm(lmS,nTheta1) + &
                                      f1ES2*wPlm(lmS,nTheta2)
            END IF
        END DO

    END DO  !  loop over theta in block

500 CONTINUE
    PERFOFF
    RETURN
    end SUBROUTINE legTF1

!------------------------------------------------------------------------
