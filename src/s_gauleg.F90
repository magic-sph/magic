!$Id$
!**********************************************************************
    SUBROUTINE GAULEG(X1,X2,X,W,N)
!**********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!----------------------------------------------------------------------
! Subroutine is based on a NR code.
! Calculates N zeros of legendre polynomial P(l=N) in
! the intervall [X1,X2].
! Zeros are returned in radiants X(i)
! The respective weights for Gauss-integration are given in W(i).
!----------------------------------------------------------------------

    IMPLICIT NONE
     
!-- INPUT:
    REAL(kind=8) :: X1,X2 ! LOWER/UPPER INTERVALL BOUND IN RADIANTS
    INTEGER :: N    ! DESIRED MAXIMUM DEGREE

!-- OUTPUT:
    REAL(kind=8) :: X(*)  ! ZEROS COS(THETA)
    REAL(kind=8) :: W(*)  ! ASSOCIATED GAUSS LEGENDRE WEIGHTS

!-- LOCAL:
    INTEGER :: M,I,J
    REAL(kind=8) :: XM,XL,P1,P2,P3,PP,Z,Z1,PI
    REAL(kind=8), PARAMETER :: EPS=3.D-14

!-- END OF DECLARATION
!-----------------------------------------------------------------------

    PI=4.D0*DATAN(1.D0)
    M=(N+1)/2  ! use symmetry

!-- Map on symmetric intervall:
    XM=0.5D0*(X2+X1)
    XL=0.5D0*(X2-X1)

    DO I=1,M

    !----- Initial guess for zeros:
        Z=DCOS( PI*( (DBLE(I)-0.25D0)/(DBLE(N)+0.5D0)) )

        1 CONTINUE  ! Jump point to refine zero !

    !----- Use recurrence to calulate P(l=n,z=cos(theta))
        P1=1.D0
        P2=0.D0
        DO J=1,N   ! do loop over degree !
            P3=P2
            P2=P1
            P1=( DBLE(2*J-1)*Z*P2-DBLE(J-1)*P3 )/DBLE(J)
        END DO

    !----- Newton method to refine zero: PP is derivative !
        PP=DBLE(N)*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
        IF( DABS( Z-Z1) > EPS ) GO TO 1 ! refine zero !

    !----- Another zero found
        X(I)=DACOS(XM+XL*Z)
        X(N+1-I)=DACOS(XM-XL*Z)
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)

    END DO
     
    RETURN
    END SUBROUTINE GAULEG

!---------------------------------------------------------------
