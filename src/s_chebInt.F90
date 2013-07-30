!$Id$
!**********************************************************************
    SUBROUTINE chebIntInit(zMin,zMax,zNorm,nNorm,nGridPointsMax, &
                           z,nGridPoints,i_costf_init,d_costf_init)
!**********************************************************************
    IMPLICIT NONE

!-- Input:
    REAL(kind=8) ::  zMin,zMax      ! integration intervall !
    REAL(kind=8) ::  zNorm          ! norm intervall length
    INTEGER :: nNorm          ! suggested number of grid points
! for norm length
! will be adjusted to nGridPoints
    INTEGER :: nGridPointsMax ! dimension of z on input

!-- Output:
    REAL(kind=8) ::  z(nGridPointsMax) ! grid points
! dimension at >= nGridPointsMax
    INTEGER :: nGridPoints       ! number of used grid points
    INTEGER :: i_costf_init(2*nGridPointsMax+2)  ! array of constants
! needed for cheb transform
! dimension >= 2*nGridPointsMax+2
    REAL(kind=8) ::  d_costf_init(2*nGridPointsMax+5)  ! array of constants
! needed for cheb transform
! dimension >= 2*nGridPointsMax

!-- Local:
    INTEGER, PARAMETER :: nChebMax=722
    INTEGER :: n
    REAL(kind=8) :: zCheb(nChebMax)
    INTEGER, PARAMETER :: nGridMin=4*6+1

!---------------------------------------------------------------------

!-- Adjust number of z points:
    n=IDINT(DBLE(nNorm-1)/4.D0*((zMax-zMin)/zNorm))
    IF ( n < 2 )  THEN
        n=2
    ELSE IF ( n == 7 )  THEN
        n=8
    ELSE IF ( n == 11 ) THEN
        n=12
    ELSE IF ( n == 13 ) THEN
        n=12
    ELSE IF ( n == 14 ) THEN
        n=15
    ELSE IF ( n == 17 ) THEN
        n=16
    ELSE IF ( n == 19 ) THEN
        n=20
    ELSE IF ( n == 21 .OR. n == 22 .OR. n == 23 ) THEN
        n=24
    !---    After this we use more course steps:
    ELSE IF ( n > 24  .AND. n <= 30 )  THEN
        n=30
    ELSE IF ( n > 30  .AND. n <= 40 )  THEN
        n=40
    ELSE IF ( n > 40  .AND. n <= 50 )  THEN
        n=50
    ELSE IF ( n > 50  .AND. n <= 60 )  THEN
        n=60
    ELSE IF ( n > 60  .AND. n <= 75 )  THEN
        n=75
    ELSE IF ( n > 75  .AND. n <= 90 )  THEN
        n=90
    ELSE IF ( n > 90  .AND. n <= 100 ) THEN
        n=100
    ELSE IF ( n > 100 .AND. n <= 120 ) THEN
        n=120
    ELSE IF ( n > 120 .AND. n <= 135 ) THEN
        n=135
    ELSE IF ( n > 135 .AND. n <= 150 ) THEN
        n=150
    ELSE IF ( n > 135 .AND. n <= 150 ) THEN
        n=150
    ELSE IF ( n > 150 .AND. n <= 180 ) THEN
        n=180
    ELSE IF ( n > 180 .AND. n <= 200 ) THEN
        n=200
    ELSE IF ( n > 200 .AND. n <= 225 ) THEN
        n=225
    ELSE IF ( n > 225 .AND. n <= 250 ) THEN
        n=250
    ELSE IF ( n > 250 ) THEN
    ! Maximum number of grid points set to 1001:
        WRITE(*,*) 'Sorry, no provision for more than 1001 points!'
        WRITE(*,*) 'Sorry, no provision for more than 1001 points!'
        WRITE(*,*) 'Sorry, no provision for more than 1001 points!'
        n=200
    END IF
    nGridPoints=4*n+1
! New minimum number of grid points (see above):
    nGridPoints=MAX(nGridPoints,nGridMin)

!        WRITE(*,*)
!        WRITE(*,*) '! Interval    :',zMin,zMax,zNorm
!        WRITE(*,*) '! No of points:',nGridPoints,nNorm

    IF ( nGridPointsMax < nGridPoints ) THEN
        WRITE(*,*) '! nGridPointsMax too small in chebIntInit!'
        WRITE(*,*) '! Should be at least:',nGridPoints
        STOP
    END IF

    IF ( nChebMax < nGridPoints ) THEN
        WRITE(*,*) '! Increase nChebMax in s_chebIntInit.f!'
        WRITE(*,*) '! Should be at least:',nGridPoints
        STOP
    END IF

!-- Calculate nGridPoints grid points z(*) in interval [zMin,zMax]
!   zCheb are the grid points in the cheb interval [-1,1]
!   These are not really needed for the integration.
    CALL cheb_x_map_e(zMin,zMax,nGridPoints-1,z,zCheb)

!-- Initialize fast cos transform for chebs:
    CALL init_costf1(nGridPoints,i_costf_init,2*nGridPointsMax+2, &
    d_costf_init,2*nGridPointsMax+5)

    RETURN
    end SUBROUTINE chebIntInit
!----------------------------------------------------------------------
    REAL(kind=8) FUNCTION chebInt(f,zMin,zMax,nGridPoints,nGridPointsMax, &
                                  i_costf_init,d_costf_init)
!**********************************************************************
    IMPLICIT NONE

!-- Input:
    REAL(kind=8) ::  f(*)               ! function on grid points
    REAL(kind=8) ::  zMin,zMax          ! integration boundaries
    INTEGER :: nGridPoints        ! No of grid points
    INTEGER :: nGridPointsMax     ! No of max grid points

    INTEGER :: i_costf_init(2*nGridPointsMax+2) ! help array
    REAL(kind=8) :: d_costf_init(2*nGridPointsMax+5) ! help array

!-- Local:
    INTEGER, PARAMETER :: nWorkMax=722  ! dimension for work array
    REAL(kind=8) :: work(nWorkMax)      ! work array
    REAL(kind=8) :: fr(nWorkMax)        ! function in cheb space
    REAL(kind=8) :: chebNorm
    INTEGER :: nCheb              ! counter for chebs
    INTEGER :: nGrid              ! counter for grid points

!----------------------------------------------------------------------

    chebNorm=DSQRT(2.D0/DBLE(nGridPoints-1))

    IF ( nWorkMax < nGridPoints ) THEN
        WRITE(*,*) '! Increase nWorkMax in chebInt!'
        WRITE(*,*) '! Should be at least:',nGridPoints!'
        STOP
    END IF

!-- Copy function:
    DO nGrid=1,nGridPoints
        fr(nGrid)=f(nGrid)
    END DO

!-- Transform to cheb space:
    CALL costf1(fr,1,1,1,work,i_costf_init,d_costf_init)
    fr(1)          =0.5D0*fr(1)
    fr(nGridPoints)=0.5D0*fr(nGridPoints)

!-- Integration:
    chebInt=0.D0
    DO nCheb=1,nGridPoints,2  ! only even chebs contribute
        chebInt=chebInt - (zMax-zMin)/DBLE(nCheb*(nCheb-2))*fr(nCheb)
    END DO

!-- Normalize with intervall:
    chebInt=chebNorm*chebInt/(zMax-zMin)


    RETURN
    END function chebInt
!----------------------------------------------------------------------
    real(kind=8) FUNCTION chebIntD(f,lDeriv,zMin,zMax, &
                           nGridPoints,nGridPointsMax, &
                               i_costf_init,d_costf_init)
!**********************************************************************
    IMPLICIT NONE

!-- Input:
    REAL(kind=8) ::  f(*)               ! function on grid points
! changed on output
    LOGICAL :: lDeriv             !
    REAL(kind=8) ::  zMin,zMax          ! integration boundaries
    INTEGER :: nGridPoints        ! No of grid points
    INTEGER :: nGridPointsMax     ! No of max grid points

    INTEGER :: i_costf_init(2*nGridPointsMax+2) ! help array
    real(kind=8) :: d_costf_init(2*nGridPointsMax+5) ! help array

!-- Output:
!       real(kind=8) f(*) overwritten for z derivative

!-- Local:
    INTEGER, PARAMETER :: nWorkMax=722  ! dimension for work array
    REAL(kind=8) :: work(nWorkMax)      ! work array
    REAL(kind=8) :: chebNorm
    REAL(kind=8) :: drFac               ! transform fac from cheb space
    INTEGER :: nCheb              ! counter for chebs

!----------------------------------------------------------------------

    chebNorm=DSQRT(2.D0/DBLE(nGridPoints-1))

    IF ( nWorkMax < nGridPoints ) THEN
        WRITE(*,*) '! Increase nWorkMax in chebIntD!'
        WRITE(*,*) '! Should be at least:',nGridPoints!'
        STOP
    END IF

!-- Transform to cheb space:
    CALL costf1(f,1,1,1,work,i_costf_init,d_costf_init)

!----- Copy:
    IF ( lDeriv ) THEN
        DO nCheb=1,nGridPoints
            work(nCheb)=f(nCheb)
        END DO
    END IF

!-- Integration:
    f(1)          =0.5D0*f(1)
    f(nGridPoints)=0.5D0*f(nGridPoints)
    chebIntD=0.D0
    DO nCheb=1,nGridPoints,2  ! only even chebs contribute
        chebIntD=chebIntD - (zMax-zMin)/DBLE(nCheb*(nCheb-2))*f(nCheb)
    END DO
!-- Normalize with intervall:
    chebIntD=chebNorm*chebIntD/(zMax-zMin)

!-- Get derivatives:
    IF ( lDeriv ) THEN
        drFac=2.D0/(zMax-zMin)
        CALL get_dcheb(work,f,1,1,1,nGridPointsMax,nGridPoints,drFac)
    !-- Transform back to grid space:
        CALL costf1(f,1,1,1,work,i_costf_init,d_costf_init)
    END IF


    RETURN
    END FUNCTION chebIntD
!----------------------------------------------------------------------
