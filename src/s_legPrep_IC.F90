!$Id$
!***********************************************************************
    SUBROUTINE legPrep_IC(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc, &
                                  r,r_ICB,lDeriv,lHor,lCondIC, &
                                    dLhw,vhG,vhC,dLhz,cvhG,cvhC)
!***********************************************************************

!    !------------ This is release 3 level 1  --------------!
!    !------------ Created on 12/08/05  by JW. -------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to prepare Legendre transforms     |
!  |  from (r,l,m) space to (r,theta,m) space by calculating           |
!  |  auxiliary arrays dLhw,vhG, ......... which contain               |
!  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
!  |  factors as well as possible the radial dependencies.             |
!  |   lHor=.TRUE. horizontal components required                      |
!  |   lDeriv=.TRUE. field derivatives required  for curl of field     |
!  |  Note: This routine is used for the inner core magnetic field     |
!  |  which has a special radial function ansatz. It can also be       |
!  |  used to prepare the calculation of a field in an insulating      |
!  |  inner core for lCondIC=.FALSE.. For this the w has to be the     |
!  |  outer core poloidal field and nR is the grid point for the ICB.  |
!  |  In any case legTF can be used for the following Legendre         |
!  |  transform and fftJW for the Fourier transform.                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    IMPLICIT NONE

!-- Input:
    INTEGER :: lm_max
    COMPLEX(kind=8) :: w(lm_max)
    COMPLEX(kind=8) :: dw(lm_max)
    COMPLEX(kind=8) :: ddw(lm_max)
    COMPLEX(kind=8) :: z(lm_max)
    COMPLEX(kind=8) :: dz(lm_max)
    REAL(kind=8) :: dLh(lm_max)
    INTEGER :: l_max,minc
    REAL(kind=8) :: r,r_ICB
    LOGICAL :: lHor
    LOGICAL :: lDeriv
    LOGICAL :: lCondIC

!-- Output:
    COMPLEX(kind=8) :: dLhw(lm_max),dLhz(lm_max)
    COMPLEX(kind=8) :: vhG(lm_max),vhC(lm_max)
    COMPLEX(kind=8) :: cvhG(lm_max),cvhC(lm_max)
     
!-- Local:
    INTEGER :: lm,l,m
    COMPLEX(kind=8) :: help1,help2
            
    INTEGER, PARAMETER :: l_maxL=1000
    REAL(kind=8) :: rRatio,rDep(0:l_maxL),rDep2(0:l_maxL)

!-- end of declaration
!--------------------------------------------------------------------------

    IF ( l_max > l_maxL+1 ) THEN
        WRITE(*,*) 'Message from legPrep_IC:'
        WRITE(*,*) 'Increase l_maxL to ',l_max-1
        STOP
    END IF

    rRatio  =r/r_ICB
    rDep(0) =rRatio
    rDep2(0)=1.D0/r_ICB ! rDep2=rDep/r
    DO l=1,l_max
        rDep(l) =rDep(l-1)*rRatio
        rDep2(l)=rDep2(l-1)*rRatio
    END DO

    lm=0
    DO m=0,l_max,minc
        DO l=m,l_max
            lm=lm+1
            dLhw(lm)=rDep(l)*dLh(lm)*w(lm)
            IF ( lHor ) THEN
                IF ( lCondIC ) THEN
                    help1=rDep2(l)*((l+1)*w(lm)+r*dw(lm))
                    help2=rDep(l)*z(lm)
                    vhG(lm)=help1-CMPLX(-AIMAG(help2),REAL(help2),KIND=KIND(0d0) )
                    vhC(lm)=help1+CMPLX(-AIMAG(help2),REAL(help2),KIND=KIND(0d0) )
                ELSE
                    vhG(lm)=rDep2(l)*(l+1)*w(lm) ! Only poloidal
                    vhC(lm)=vhG(lm)
                END IF
            END IF
        END DO
    END DO
    dLhw(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    IF ( lHor ) THEN
        vhG(1) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        vhC(1) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    END IF

    IF ( lDeriv ) THEN
        lm=0
        DO m=0,l_max,minc
            DO l=m,l_max
                lm=lm+1
                IF ( lCondIC ) THEN
                    dLhz(lm)=rDep(l)*dLh(lm)*z(lm)
                    help1=rDep(l)*( (l+1)*z(lm)/r+dz(lm) )
                    help2=rDep(l)*(-2.D0*(l+1)/r*dw(lm)-ddw(lm))
                    cvhG(lm)=help1-CMPLX(-AIMAG(help2),REAL(help2),KIND=KIND(0d0))
                    cvhC(lm)=help1+CMPLX(-AIMAG(help2),REAL(help2),KIND=KIND(0d0))
                ELSE
                    dLhz(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    cvhG(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    cvhC(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
            END DO
        END DO
        dLhz(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        cvhG(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        cvhc(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    END IF


    RETURN
    end SUBROUTINE legPrep_IC

!---------------------------------------------------------------------------
