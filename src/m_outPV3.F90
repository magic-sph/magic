!$Id$
MODULE outPV3
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data

  IMPLICIT NONE 

  private
  INTEGER,PARAMETER :: nSmaxA=97
  INTEGER,PARAMETER :: nZmaxA=2*nSmaxA
  REAL(kind=8),ALLOCATABLE :: rZ(:,:)
  REAL(kind=8),ALLOCATABLE :: PlmS(:,:,:)
  REAL(kind=8),ALLOCATABLE :: dPlmS(:,:,:)
  REAL(kind=8),ALLOCATABLE :: PlmZ(:,:,:)
  REAL(kind=8),ALLOCATABLE :: dPlmZ(:,:,:)
  REAL(kind=8),ALLOCATABLE :: OsinTS(:,:)
  REAL(kind=8),ALLOCATABLE :: VorOld(:,:,:)

  public :: initialize_outPV3,outPV
  
contains
  SUBROUTINE initialize_outPV3

    ALLOCATE( rZ(nZmaxA/2+1,nSmaxA) )
    ALLOCATE( PlmS(l_max+1,nZmaxA/2+1,nSmaxA) )
    ALLOCATE( dPlmS(l_max+1,nZmaxA/2+1,nSmaxA) )
    ALLOCATE( PlmZ(lm_max,nZmaxA/2+1,nSmaxA) )
    ALLOCATE( dPlmZ(lm_max,nZmaxA/2+1,nSmaxA) )
    ALLOCATE( OsinTS(nZmaxA/2+1,nSmaxA) )
    ALLOCATE( VorOld(nrp,nZmaxA,nSmaxA) )

  END SUBROUTINE initialize_outPV3

  !***********************************************************************
  SUBROUTINE outPV(time,l_stop_time,nPVsets,w,dw,ddw,               &
       &                            z,dz,omega_IC,omega_MA)
    !***********************************************************************

    !-----------------------------------------------------------------------
    !   Output of z-integrated axisymmetric rotation rate Vp/s 
    !   and s derivatives
    !-----------------------------------------------------------------------

    USE plms_theta, ONLY: plm_theta

    !-- Input of variables:
    REAL(kind=8) :: time
    INTEGER :: nPVsets
    REAL(kind=8) :: omega_IC,omega_MA
    LOGICAL :: l_stop_time

    !-- Input of toroidal flow scalar and rotation rates:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: dw(lm_max,n_r_max)
    COMPLEX(kind=8) :: ddw(lm_max,n_r_max)
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    COMPLEX(kind=8) :: dz(lm_max,n_r_max)

    !-- (l,r) Representation of the different contributions
    REAL(kind=8) :: dzVpLMr(l_max+1,n_r_max)

    !--- Work array:
    COMPLEX(kind=8) :: workA(lm_max,n_r_max)

    INTEGER :: lm,l,m ! counter for degree and order

    REAL(kind=8) :: fac

    !--- define Grid
    INTEGER :: nSmax,nS,nSI
    REAL(kind=8) ::  sZ(nSmaxA),dsZ ! cylindrical radius s and s-step

    INTEGER :: nZ,nZmax,nZmaxNS,nZN
    !PARAMETER (nZmaxA=2*nSmaxA)
    INTEGER :: nZC(nSmaxA),nZ2(nZmaxA,nSmaxA),nZS
    REAL(kind=8) :: zZ(nZmaxA),zstep!,zZC
    REAL(kind=8) :: VpAS(nZmaxA),omS(nZmaxA)
    SAVE nZC,nZ2


    !-- For integration along z:

    !-- Plms: Plm,sin
    INTEGER :: nR,nPhi,nC
    REAL(kind=8) :: thetaZ,rZS!,sinT,cosT

    !-- For PV output files: 
    CHARACTER(len=80) :: fileName

    !-- Output of all three field components:
    REAL(kind=8) :: VsS(nrp,nZmaxA)
    REAL(kind=8) :: VpS(nrp,nZmaxA)
    REAL(kind=8) :: VzS(nrp,nZmaxA)
    REAL(kind=8) :: VorS(nrp,nZmaxA)
    REAL(kind=8) :: dpVorS(nrp,nZmaxA)
    REAL(kind=4) :: out1(n_phi_max*nZmaxA)
    REAL(kind=4) :: out2(n_phi_max*nZmaxA)
    REAL(kind=4) :: out3(n_phi_max*nZmaxA)
    REAL(kind=4) :: out4(n_phi_max*nZmaxA)
    REAL(kind=4) :: out5(n_phi_max*nZmaxA)
    REAL(kind=8) :: timeOld
    SAVE timeOld


    !-- This may be deleted later:
    COMPLEX(kind=8) :: wP(lm_max,n_r_max)
    COMPLEX(kind=8) :: dwP(lm_max,n_r_max)
    COMPLEX(kind=8) :: ddwP(lm_max,n_r_max)
    COMPLEX(kind=8) :: zP(lm_max,n_r_max)
    COMPLEX(kind=8) :: dzP(lm_max,n_r_max)


    !-- End of declaration
    !---------------------------------------------------------------------------

    IF ( lVerbose ) WRITE(*,*) '! Starting outPV!'

    nPVsets=nPVsets+1

    !-- Start with calculating advection due to axisymmetric flows:

    nSmax=n_r_max+INT(r_ICB*DBLE(n_r_max))
    nSmax=INT(sDens*nSmax)
    IF ( nSmax.GT.nSmaxA ) THEN
       WRITE(*,*) 'Increase nSmaxA in outPV!'
       WRITE(*,*) 'Should be at least nSmax=',nSmax
       WRITE(*,*) 'But is only=',nSmaxA
       STOP
    END IF
    nZmax=2*nSmax

    IF ( l_stop_time ) THEN
       IF ( l_SRIC  .AND. omega_IC.NE.0 ) THEN
          fac=1.D0/omega_IC
       ELSE
          fac=1.D0
       END IF
       DO nR=1,n_r_max
          DO l=1,l_max
             lm=lm2(l,0)
             dzVpLMr(l+1,nR)=fac*REAL(z(lm,nR))
          END DO
       END DO

       !---- Transform the contributions to cheb space:
       CALL costf1(dzVpLMr,l_max+1,1,l_max+1,                       &
            &         workA,i_costf_init,d_costf_init)
    END IF

    !--- Transforming of field without the backtransform
    !    Thus this must be the last thing done with the 
    !    fields in a run. See s_output.f and s_step_time.f.
    !    NOTE: output is only non-axisymmetric part!
    DO nR=1,n_r_max
       DO lm=1,lm_max
          m=lm2m(lm)
          !           IF ( m.EQ.0 ) THEN
          !             wP(lm,nR)  =0.D0
          !              dwP(lm,nR) =0.D0
          !              ddwP(lm,nR)=0.D0
          !              zP(lm,nR)  =0.D0
          !              dzP(lm,nR) =0.D0
          !          ELSE
          wP(lm,nR)  =w(lm,nR)*dLh(lm)
          dwP(lm,nR) =dw(lm,nR)
          ddwP(lm,nR)=ddw(lm,nR)
          zP(lm,nR)  =z(lm,nR)
          dzP(lm,nR) =dz(lm,nR)
          !          END IF
       END DO
    END DO

    !---- Transform the contributions to cheb space for z-integral:
    CALL costf1(wP,lm_max_real,1,lm_max_real,                       &
         &          workA,i_costf_init,d_costf_init)
    CALL costf1(dwP,lm_max_real,1,lm_max_real,                      &
         &           workA,i_costf_init,d_costf_init)
    CALL costf1(ddwP,lm_max_real,1,lm_max_real,                     &
         &            workA,i_costf_init,d_costf_init)
    CALL costf1(zP,lm_max_real,1,lm_max_real,                       &
         &          workA,i_costf_init,d_costf_init)
    CALL costf1(dzP,lm_max_real,1,lm_max_real,                      &
         &           workA,i_costf_init,d_costf_init)

    dsZ=r_CMB/DBLE(nSmax)  ! Step in s controlled by nSmax
    nSI=0                  ! Inner core position
    DO nS=1,nSmax
       sZ(nS)=(nS-0.5D0)*dsZ
       IF ( sZ(nS).LT.r_ICB .AND. nS.GT.nSI ) nSI=nS
    END DO
    zstep=2*r_CMB/DBLE(nZmax-1)
    DO nZ=1,nZmax
       zZ(nZ)=r_CMB-(nZ-1)*zstep
    END DO

    !--- Open file for output:
    IF ( l_stop_time ) THEN
       fileName='PVZ.'//TAG
       OPEN(95,FILE=fileName,FORM='UNFORMATTED',                    &
            &                             STATUS='UNKNOWN')
       WRITE(95) SNGL(time),FLOAT(nSmax),FLOAT(nZmax),              &
            &                      SNGL(omega_IC),SNGL(omega_ma)
       WRITE(95) (SNGL(sZ(nS)),nS=1,nSmax)
       WRITE(95) (SNGL(zZ(nZ)),nZ=1,nZmax)


       !--- Open file for the three flow components:
       fileName='Vcy.'//TAG
       OPEN(96,FILE=fileName,FORM='UNFORMATTED',                    &
            &                             STATUS='UNKNOWN')
       WRITE(96) SNGL(time),FLOAT(nSmax),FLOAT(nZmax),              &
            &      FLOAT(n_phi_max),SNGL(omega_IC),SNGL(omega_ma),             &
            &      SNGL(radratio),FLOAT(minc)
       WRITE(96) (SNGL(sZ(nS)),nS=1,nSmax)
       WRITE(96) (SNGL(zZ(nZ)),nZ=1,nZmax)
    END IF


    !$OMP DO ORDERED                                                        &
    !$OMP  PRIVATE(nS,nZ,rZS,thetaZ,nZmaxNS,nZS,VpAS,omS,VsS,VzS,VpS,VorS,  &
    !$OMP          dpVorS,out1,out2,out3,out4,out5,nC,nPhi,nZN)

    DO nS=1,nSmax

       !------ Get r,theta,Plm,dPlm for northern hemishere:
       IF ( nPVsets.EQ.1 ) THEN ! DO this only for the first call !
          nZC(nS)=0 ! Points within shell
          DO nZ=1,nZmax
             rZS=DSQRT(zZ(nZ)**2+sZ(nS)**2)
             IF ( rZS.GE.r_ICB .AND. rZS.LE.r_CMB ) THEN
                nZC(nS)=nZC(nS)+1  ! Counts all z within shell
                nZ2(nZ,nS)=nZC(nS) ! No of point within shell
                IF ( zZ(nZ).GT.0 ) THEN ! Onl north hemisphere
                   rZ(nZC(nS),nS)=rZS
                   thetaZ=DATAN2(sZ(nS),zZ(nZ))
                   OsinTS(nZC(nS),nS)=1.D0/DSIN(thetaZ)
                   CALL plm_theta(thetaZ,l_max,0,minc,              &
                        &    PlmS(1,nZC(nS),nS),dPlmS(1,nZC(nS),nS),l_max+1,2)
                   CALL plm_theta(thetaZ,l_max,m_max,minc,          &
                        &        PlmZ(1,nZC(nS),nS),dPlmZ(1,nZC(nS),nS),lm_max,2)
                END IF
             ELSE
                nZ2(nZ,nS)=-1 ! No z found within shell !
             END IF
          END DO
       END IF

       !-- Get azimuthal flow component in the shell
       nZmaxNS=nZC(nS) ! all z points within shell
       IF ( l_stop_time ) THEN
          CALL getPAStr(VpAS,dzVpLMr,nZmaxNS,nZmaxA,l_max+1,        &
               &                          l_max,r_ICB,r_CMB,n_r_max,         &
               &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))

          !-- Copy to array with all z-points
          DO nZ=1,nZmax
             rZS=DSQRT(zZ(nZ)**2+sZ(nS)**2)
             nZS=nZ2(nZ,nS)
             IF ( nZS.GT.0 ) THEN
                omS(nZ)=VpAS(nZS)/sZ(nS)
             ELSE
                IF ( rZS.LE.r_ICB ) THEN
                   omS(nZ)=1.D0
                ELSE
                   omS(nZ)=fac*omega_MA
                END IF
             END IF
          END DO
       END IF

       !-- Get all three components in the shell
       CALL getPVptr(wP,dwP,ddwP,zP,dzP,r_ICB,r_CMB,rZ(1,nS),        &
            &        nZmaxNS,nZmaxA,PlmZ(1,1,nS),dPlmZ(1,1,nS),OsinTS(1,nS),   &
            &                                       VsS,VpS,VzS,VorS,dpVorS)

       !$OMP ORDERED
       IF ( l_stop_time ) THEN
          WRITE(95) (SNGL(omS(nZ)),nZ=1,nZmax)
          WRITE(96) FLOAT(nZmaxNS)
          nC=0
          DO nZ=1,nZmaxNS
             DO nPhi=1,n_phi_max
                nC=nC+1
                out1(nC)=SNGL(VsS(nPhi,nZ)) ! Vs
                out2(nC)=SNGL(VpS(nPhi,nZ)) ! Vphi
                out3(nC)=SNGL(VzS(nPhi,nZ)) ! Vz
                out4(nC)=SNGL(VorS(nPhi,nZ))
                out5(nC)=(SNGL(VorS(nPhi,nZ)-VorOld(nPhi,nZ,nS)))/     &
                     &                    (SNGL(time-timeOld))
             END DO
          END DO
          WRITE(96) (out1(nZ),nZ=1,nC)
          WRITE(96) (out2(nZ),nZ=1,nC)
          WRITE(96) (out3(nZ),nZ=1,nC)
          WRITE(96) (out4(nZ),nZ=1,nC)
          WRITE(96) (out5(nZ),nZ=1,nC)
       ELSE
          timeOld=time
          DO nZ=1,nZmaxNS
             DO nPhi=1,n_phi_max
                VorOld(nPhi,nZ,nS)=VorS(nPhi,nZ)
             END DO
          END DO
       END IF
       !$OMP END ORDERED


    END DO  ! Loop over s 
    !$OMP END DO

    IF ( l_stop_time ) CLOSE (95)
    IF ( l_stop_time ) CLOSE(96)

    RETURN
  END SUBROUTINE outPV

  !---------------------------------------------------------------------------------
END MODULE outPV3
