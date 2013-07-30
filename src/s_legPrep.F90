!$Id$
!***********************************************************************
    SUBROUTINE legPrep(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc, &
                                             r,lDeriv,lHor, &
                                 dLhw,vhG,vhC,dLhz,cvhG,cvhC)
!***********************************************************************

!    !------------ This is release 3 level 1  --------------!
!    !------------ Created on 12/08/05  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to prepare Legendre transforms     |
!  |  from (r,l,m) space to (r,theta,m) space by calculating           |
!  |  auxiliary arrays w, dw, ddw, ....... which contain               |
!  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
!  |  factors as well as possible the radial dependencies.             |
!  |   lHor=.TRUE. horizontal components required                      |
!  |   lDeriv=.TRUE. field derivatives required  for curl of field     |
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
    REAL(kind=8) :: r
    LOGICAL :: lHor
    LOGICAL :: lDeriv

!-- Output:
    COMPLEX(kind=8) :: dLhw(*),dLhz(*)
    COMPLEX(kind=8) :: vhG(*),vhC(*)
    COMPLEX(kind=8) :: cvhG(*),cvhC(*)
     
!-- Local:
    INTEGER :: lm,l,m
    REAL(kind=8) :: Or_e2
    COMPLEX(kind=8) :: help

!-- end of declaration
!--------------------------------------------------------------------------

    lm=0
    DO m=0,l_max,minc
        DO l=m,l_max
            lm=lm+1
            dLhw(lm)=dLh(lm)*w(lm)
            IF ( lHor ) THEN
                vhG(lm) =dw(lm)-CMPLX(-AIMAG(z(lm)),REAL(z(lm)),KIND=KIND(0d0))
                vhC(lm) =dw(lm)+CMPLX(-AIMAG(z(lm)),REAL(z(lm)),KIND=KIND(0d0))
            END IF
        END DO
    END DO
    dLhw(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    IF ( lHor ) THEN
        vhG(1) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        vhC(1) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    END IF

    IF ( lDeriv ) THEN
        Or_e2=1.D0/r**2
        lm=0
        DO m=0,l_max,minc
            DO l=m,l_max
                lm=lm+1
                dLhz(lm)=dLh(lm)*z(lm)
                help=dLh(lm)*Or_e2*w(lm)-ddw(lm)
                cvhG(lm)=dz(lm)-CMPLX(-AIMAG(help),REAL(help),KIND=KIND(0d0))
                cvhC(lm)=dz(lm)+CMPLX(-AIMAG(help),REAL(help),KIND=KIND(0d0))
            END DO
        END DO
        dLhz(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        cvhG(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        cvhc(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    END IF

    RETURN
    end SUBROUTINE legPrep

!------------------------------------------------------------------------
