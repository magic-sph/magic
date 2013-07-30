!$Id$
MODULE omega
  
  IMPLICIT NONE

  CONTAINS
!***********************************************************************
    SUBROUTINE outOmega(z,omega_IC)
!***********************************************************************

!-----------------------------------------------------------------------
!   Output of axisymmetric zonal flow omega(s) into field omega.TAG,
!   where s is the cylindrical radius. This is done for the southern
!   and norther hemispheres at z=+-(r_icb+0.5)
!-----------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE blocking
    USE logic
    USE output_data
    USE const

    IMPLICIT NONE

!-- Input scalar fields:
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    REAL(kind=8) :: omega_IC
            
!-- Local stuff:

    REAL(kind=8) :: dzVpLMr(l_max+1,n_r_max)

    INTEGER :: nR,lm,l ! counter
    INTEGER :: nNS     ! index for NHS and SHS

    INTEGER,PARAMETER :: nSmax=300
    INTEGER :: nS
    REAL(kind=8) ::  sZ,zZ,dsZ
    REAL(kind=8) :: rZ,thetaZ
    REAL(kind=8) :: VpS,omega(2)
    COMPLEX(kind=8) :: workA(lm_max,n_r_max) ! work array

    CHARACTER(len=64) :: fileName

!-- End of declaration
!---------------------------------------------------------------------------

    IF ( lVerbose ) WRITE(*,*) '! Starting outOmega!'

    fileName='omega.'//tag
    OPEN(99,FILE=fileName,STATUS='UNKNOWN')

    dsZ=r_CMB/DBLE(nSmax-1)

!--- Transform to lm-space for all radial grid points:
    DO nR=1,n_r_max
        dzVpLMr(1,nR)=0.D0
        DO l=1,l_max
            lm=lm2(l,0)
            dzVpLMr(l+1,nR)=REAL(z(lm,nR))
        END DO
    END DO

!---- Transform the contributions to cheb space for z-integral:
    CALL costf1(dzVpLMr,l_max+1,1,l_max+1, &
                workA,i_costf_init,d_costf_init)

    sZ=0.D0
    DO nS=1,nSmax
        sZ=sZ+dsZ

        DO nNS=1,2  !  North and south hemisphere !

            IF ( nNS == 1 ) THEN  ! south hemisphere !
                zZ=r_ICB+0.5D0
            ELSE
                zZ=-(r_ICB+0.5D0)
            END IF
            rZ    =DSQRT(zZ*zZ+sZ*sZ)
            thetaZ=DATAN2(sZ,zZ)
            IF ( rZ > r_CMB ) GOTO 100

        !------ Get the function values for (sZ,zCy)
            VpS=lnPAS2tr(dzVpLMr,l_max+1,r_ICB,r_CMB, &
                         l_max,minc,n_r_max,thetaZ,rZ)
            omega(nNS)=VpS/(rZ*DSIN(thetaZ))/omega_IC

        END DO  ! Loop over north and south hemisphere

        WRITE(99,*) sZ,omega(1),omega(2)

    END DO  ! Loop over s
    100 CONTINUE

    CLOSE(99)

    IF ( lVerbose ) WRITE(*,*) '! End of outOmega!'


    RETURN
    end SUBROUTINE outOmega
!-----------------------------------------------------------------------------
    REAL(kind=8) FUNCTION lnPAS2tr(f,lmMax,a,b,lMax,minc,nChebMax,theta,r)
!-----------------------------------------------------------------------------

    USE plms_theta, ONLY: plm_theta

    IMPLICIT NONE

    INTEGER :: lmMax
    REAL(kind=8) :: f(lmMax,*)
    REAL(kind=8) :: a,b
    INTEGER :: lMax
    INTEGER :: minc
    INTEGER :: nChebMax
    REAL(kind=8) :: theta
    REAL(kind=8) :: r
    INTEGER,PARAMETER :: lmMaxLocal=5000
    REAL(kind=8) :: plm(lmMaxLocal),dthetaPlm(lmMaxLocal)
    INTEGER :: nCheb
    INTEGER,PARAMETER :: nChebMaxLocal=192
    REAL(kind=8) :: cheb(nChebMaxLocal)
    REAL(kind=8) :: ftr
    INTEGER :: l
    REAL(kind=8) :: x,chebNorm

!----------------------------------------------------------------------------


    IF ( nChebMax > nChebMaxLocal ) THEN
        WRITE(*,*) '! nChebMaxH too small in lnPAS2tr!'
        WRITE(*,*) '! Should be at least:',nChebMax
        STOP
    END IF
    IF ( lmMax > lmMaxLocal ) THEN
        WRITE(*,*) '! lmMaxLocal small in lnPAS2tr!'
        WRITE(*,*) '! Should be at least:',lmMax
        STOP
    END IF

!--- Map r to cheb intervall [-1,1]:
!    and calculate the cheb polynomia:
!    Note: the factor chebNorm is needed
!    for renormalisation. Its not needed if one used
!    costf1 for the back transform.
    x=2.D0*(r-0.5D0*(a+b))/(b-a)
    chebNorm=DSQRT(2.D0/DBLE(nChebMax-1))
    cheb(1)=1.D0*chebNorm
    cheb(2)=x*chebNorm
    DO nCheb=3,nChebMax
        cheb(nCheb)=2.D0*x*cheb(nCheb-1)-cheb(nCheb-2)
    END DO
    cheb(1)       =0.5D0*cheb(1)
    cheb(nChebMax)=0.5D0*cheb(nChebMax)

!--- Calculate Legendres:
    CALL plm_theta(theta,lMax,0,minc,plm,dthetaPlm,lmMaxLocal,2)

!--- Loop to add all contribution functions:
    ftr=0.D0
    DO l=0,lMax
        IF ( l+1 > lmMax ) THEN
            WRITE(*,*) '! lmMax too small in ln2rt!'
            WRITE(*,*) '! lm',l+1,lMax
            STOP
        END IF
        DO nCheb=1,nChebMax
            ftr=ftr+f(l+1,nCheb)*dthetaPlm(l+1)*cheb(nCheb)
        END DO
    END DO

    lnPAS2tr=-ftr/(r*DSIN(theta))

    END FUNCTION lnPAS2tr
!----------------------------------------------------------------------
END MODULE omega
