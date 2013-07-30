!$Id$
MODULE Egeos_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data, only: zDens
  USE const
  USE LMLoop_data,ONLY: llm,ulm,llm_real,ulm_real
  USE communications,ONLY: gather_all_from_lo_to_rank0,gt_OC
  IMPLICIT NONE 

  INTEGER,PARAMETER :: nZmaxL=lGeos*144
  INTEGER,PARAMETER :: nZmaxA=MAX0(2,nZmaxL)  ! Data along z !

  INTEGER,PARAMETER :: nSmaxL=lGeos*73
  INTEGER,PARAMETER :: nSmaxA=MAX0(1,nSmaxL)

  REAL(kind=8),allocatable :: PlmS(:,:,:)  ! This is huge !
  REAL(kind=8),ALLOCATABLE :: dPlmS(:,:,:) ! This is huge !

contains
  SUBROUTINE initialize_Egeos_mod

    allocate( PlmS(lm_maxGeos,nZmaxA/2+1,nSmaxA) )  ! This is huge !
    allocate( dPlmS(lm_maxGeos,nZmaxA/2+1,nSmaxA) ) ! This is huge !

  END SUBROUTINE initialize_Egeos_mod

  !***********************************************************************
  SUBROUTINE getEgeos(time,nGeosSets,w,dw,ddw,z,dz,         &
       &              Egeos,EkNTC,EkSTC,Ekin,               &
       &              dpFlow,dzFlow,CVzOTC,CVorOTC,CHelOTC)
    !***********************************************************************

    !-----------------------------------------------------------------------
    !   Output of axisymmetric zonal flow, its relative strength,
    !   its time variation, and all forces acting on it.
    !   The slowest part in the TO process is the repitions calculation
    !   of Plms by subroutine plm_theta. They are needed in getDVptr 
    !   when I transform on the cylindrical grid. 
    !   The necessary plms could simply be calculated one and then 
    !   be stored for later use! See s_outTOnew.f.
    !-----------------------------------------------------------------------

    USE plms_theta, ONLY: plm_theta

    !-- Input of variables:
    REAL(kind=8),INTENT(IN) :: time
    INTEGER,INTENT(IN) :: nGeosSets 

    !-- Input of scalar fields:
    COMPLEX(kind=8),INTENT(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: dw(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: ddw(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: z(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: dz(llm:ulm,n_r_max)


    !-- Output:
    REAL(kind=8),INTENT(OUT) :: Egeos,EkNTC,EkSTC,Ekin
    REAL(kind=8),INTENT(OUT) :: dpFlow,dzFlow ! RMS length scales
    REAL(kind=8),INTENT(OUT) :: CVzOTC,CVorOTC,CHelOTC

    !-- Local:
    LOGICAL :: lDeriv
    INTEGER :: nSmax,nS,nS_ICB
    REAL(kind=8) :: ofr ! inverse Froude number (anelastic)
    REAL(kind=8) :: zNorm  ! Norm z intervall
    INTEGER :: nNorm  ! No. of grid points for norm intervall
    REAL(kind=8) :: zMin,zMax,help ! integration boundarie, help variable
    LOGICAL :: lAS    ! .TRUE. if axisymmetric (m=0) functions
    REAL(kind=8) :: sZ(nSmaxA),dsZ ! cylindrical radius s and s-step
    INTEGER :: nPhi,nI
    REAL(kind=8) :: phiNorm
    LOGICAL :: lTC

    !-- Local field copies to avoid changes by back and forth transform:
    COMPLEX(kind=8) :: wS(lm_maxGeos,n_r_maxGeos)
    COMPLEX(kind=8) :: dwS(lm_maxGeos,n_r_maxGeos)
    COMPLEX(kind=8) :: ddwS(lm_maxGeos,n_r_maxGeos)
    COMPLEX(kind=8) :: zS(lm_maxGeos,n_r_maxGeos)
    COMPLEX(kind=8) :: dzS(lm_maxGeos,n_r_maxGeos)
    !---- Additional work array
    COMPLEX(kind=8) :: workA(lm_maxGeos,n_r_maxGeos)

    !-- Representation in (phi,z):
    REAL(kind=8) :: VrS(nrpGeos,nZmaxA),VrInt(nZmaxA),VrIntS
    REAL(kind=8) :: VtS(nrpGeos,nZmaxA),VtInt(nZmaxA),VtIntS
    REAL(kind=8) :: VpS(nrpGeos,nZmaxA),VpInt(nZmaxA),VpIntS
    REAL(kind=8) :: VozS(nrpGeos,nZmaxA)
    REAL(kind=8) :: OsinTS(nZmaxA/2+1,nSmaxA),sinT,cosT
    SAVE   OsinTS
    INTEGER :: nInt,nInts   ! index for NHS and SHS integral
    INTEGER :: nZ,nZmax,nZS,nZN
    INTEGER :: lm,nR
    REAL(kind=8) :: EkInt(nZmaxA),EkIntS
    REAL(kind=8) :: EkSTC_s(nSmaxA),EkNTC_s(nSmaxA),EkOTC_s(nSmaxA)
    REAL(kind=8) :: Egeos_s(nSmaxA)
    REAL(kind=8) :: EkOTC
    REAL(kind=8) :: dpEkInt(nZmaxA),dpEkIntS
    REAL(kind=8) :: dzEkInt(nZmaxA),dzEkIntS
    REAL(kind=8) :: dpEk_s(nSmaxA),dzEk_s(nSmaxA)

    !-- For integration along z:
    INTEGER :: i_costf_initZ(2*nZmaxA+2,nSmaxA)
    REAL(kind=8) :: d_costf_initZ(2*nZmaxA+5,nSmaxA)
    REAL(kind=8) :: zZ(nZmaxA,nSmaxA)
    REAL(kind=8) :: rZ(nZmaxA,nSmaxA),thetaZ
    REAL(kind=8) :: rhoZ(nZmaxA,nSmaxA)   ! density (anelastic version)
    REAL(kind=8) :: orhoZ(nZmaxA,nSmaxA)  ! 1/rho   (anelastic version)
    INTEGER :: nZmaxS(nSmaxA)
    SAVE   i_costf_initZ,d_costf_initZ,zZ,rZ,nZmaxS

    REAL(kind=8),EXTERNAL :: chebInt,chebIntD

    LOGICAL lStopRun
    INTEGER :: l,m

    !-- Correlation (new Oct. 4 2007)
    LOGICAL :: lCorrel
    REAL(kind=8) :: VzS,VzN,VorS,VorN,surf,delz
    REAL(kind=8) :: VzSN,VzSS,VzNN,VorSN,VorSS,VorNN,HelZZ,VZZ,VorZZ
    REAL(kind=8) :: CVz_I,CVor_I,CHel_I
    REAL(kind=8) :: CVz_s(nSmaxA),CVor_s(nSmaxA),CHel_S(nSmaxA)

    !-- Movie output
    INTEGER :: nOutFile,n
    CHARACTER(len=66) :: version,movFile
    INTEGER :: nFields,nFieldSize
    REAL(kind=8) :: dumm(40)
    REAL(kind=4) :: CVz(nrpGeos,nZmaxA)
    REAL(kind=4) :: CVor(nrpGeos,nZmaxA)
    REAL(kind=4) :: CHel(nrpGeos,nZmaxA)


    !-- End of declaration
    !---------------------------------------------------------------------------

    IF ( lVerbose ) WRITE(*,*) '! Starting outGeos!'

    CALL gather_all_from_lo_to_rank0(gt_OC,w,wS)
    CALL gather_all_from_lo_to_rank0(gt_OC,dw,dwS)
    CALL gather_all_from_lo_to_rank0(gt_OC,ddw,ddwS)
    CALL gather_all_from_lo_to_rank0(gt_OC,z,zS)
    CALL gather_all_from_lo_to_rank0(gt_OC,dz,dzS)

    IF (rank.EQ.0) THEN
       lCorrel=.TRUE. ! Calculate Vz and Vorz north/south correlation

       phiNorm=2.D0*pi/n_phi_max
       lStopRun=.FALSE.
       lDeriv=.TRUE.

       !---- Get resolution in s and z for z-integral:
       zNorm=1.D0               ! This is r_CMB-r_ICB
       nNorm=INT(zDens*n_r_max) ! Covered with nNorm  points !
       nSmax=n_r_max+INT(r_ICB*DBLE(n_r_max)) 
       nSmax=INT(sDens*nSmax)
       dsZ  =r_CMB/DBLE(nSmax)  ! Step in s controlled by nSmax
       IF ( nSmax.GT.nSmaxA ) THEN
          WRITE(*,*) 'Increase nSmaxA in getGeos!'
          WRITE(*,*) 'Should be at least nSmax=',nSmax
          STOP
       END IF
       lAS=.FALSE.

       DO nS=1,nSmax
          EkSTC_s(nS)=0.D0
          EkNTC_s(nS)=0.D0
          EkOTC_s(nS)=0.D0
          Egeos_s(nS)=0.D0
          dpEk_s(nS) =0.D0
          dzEk_s(nS) =0.D0
          CVz_s(nS)  =0.D0
          CVor_s(nS) =0.D0
          CHel_s(nS) =0.D0
       END DO

       !---- Copy for following costf so that no      
       !     back transform is needed that would change
       !     the field (just a little bit) anyway...
       DO nR=1,n_r_max
          DO lm=1,lm_max
             wS(lm,nR)  =  wS(lm,nR)*dLh(lm)
          END DO
       END DO

       !---- Transform the contributions to cheb space for z-integral:
       !CALL costf1(wS,ulm_real-llm_real+1,1,ulm_real-llm_real+1,    &
       !     &      workA,i_costf_init,d_costf_init)
       !CALL costf1(dwS,ulm_real-llm_real+1,1,ulm_real-llm_real+1,   &
       !     &      workA,i_costf_init,d_costf_init)
       !CALL costf1(ddwS,ulm_real-llm_real+1,1,ulm_real-llm_real+1,  &
       !     &      workA,i_costf_init,d_costf_init)
       !CALL costf1(zS,ulm_real-llm_real+1,1,ulm_real-llm_real+1,    &
       !     &      workA,i_costf_init,d_costf_init)
       !CALL costf1(dzS,ulm_real-llm_real+1,1,ulm_real-llm_real+1,   &
       !     &      workA,i_costf_init,d_costf_init)

       CALL costf1(wS,lm_max_real,1,lm_max_real,    &
            &      workA,i_costf_init,d_costf_init)
       CALL costf1(dwS,lm_max_real,1,lm_max_real,   &
            &      workA,i_costf_init,d_costf_init)
       CALL costf1(ddwS,lm_max_real,1,lm_max_real,  &
            &      workA,i_costf_init,d_costf_init)
       CALL costf1(zS,lm_max_real,1,lm_max_real,    &
            &      workA,i_costf_init,d_costf_init)
       CALL costf1(dzS,lm_max_real,1,lm_max_real,   &
            &      workA,i_costf_init,d_costf_init)

       !---- Contributions are now in fully spectral space!

       !---- Do the z-integral:
       nI=0

       DO nS=1,nSmax
          sZ(nS)=(nS-0.5D0)*dsZ
          IF ( sZ(nS).LT.r_ICB ) THEN
             lTC=.TRUE.
          ELSE
             lTC=.FALSE.
          END IF
          IF ( nS.GT.1 ) THEN
             IF ( sZ(nS-1).LT.r_ICB .AND.                              &
                  &             sZ(ns).GE.r_ICB ) nS_ICB=nS
          END IF

          !------ Get integral boundaries for this s:
          zMin=-DSQRT(r_CMB*r_CMB-sZ(nS)*sZ(nS))
          IF ( lTC ) THEN
             zMax=-DSQRT(r_ICB*r_ICB-sZ(nS)*sZ(nS))
          ELSE
             zMax=-zMin
          ENDIF

          IF ( nGeosSets.EQ.1 ) THEN
             !------ Initialize integration for NHS:
             !       Each processor calculates Cheb transform data 
             !       for HIS nS and the Plms along the Cylinder
             !       chebIntInit returns zZ,nZmaxS,i_costf_initZ and 
             !       d_costfInitZ:
             CALL chebIntInit(zMin,zMax,zNorm,nNorm,                   &
                  &           nZmaxA,zZ(1,nS),nZmaxS(nS),              &
                  &           i_costf_initZ(1,nS),d_costf_initZ(1,nS))
             !------ Calculate and store 1/sin(theta) and Plms,dPlms for 
             !       southern HS:

             IF ( lTC ) THEN 
                nZmax=nZmaxS(nS)  ! nZmax point in each polar region
             ELSE
                nZmax=(nZmaxS(nS)-1)/2+1 ! Odd point number !    
                ! all together nZmaxS(nS) from
                ! south to north including equator
             END IF
             IF ( 2*nZmax.GT.nZmaxA ) THEN ! TC case most critical
                WRITE(*,*) '! nZmaxA too small in getEgeos!'
                WRITE(*,*) '! Should be at least:',2*nZmax
                lStopRun=.TRUE.
             END IF
             DO nZ=1,nZmax
                rZ(nZ,nS)    =DSQRT(zZ(nZ,nS)**2+sZ(nS)**2)
                thetaZ       =DATAN2(sZ(nS),zZ(nZ,nS))
                OsinTS(nZ,nS)=1.D0/DSIN(thetaZ)
                CALL plm_theta(thetaZ,l_max,m_max,minc,                &
                     &            PlmS(1,nZ,nS),dPlmS(1,nZ,nS),lm_max,2)
             END DO
             IF ( l_anel) THEN
                ofr=( DEXP(strat/polind)-1.D0 )/                       &
                     &               ( g0+0.5D0*g1*(1.D0+radratio) +g2/radratio )
                DO nZ=1,nZmax
                   rhoZ(nZ,nS)=DEXP(-ofr*(g0*(rZ(nZ,nS)-r_cmb) +       &
                        &                         g1/(2.d0*r_cmb)*(rZ(nZ,nS)**2-r_cmb**2) -&
                        &                         g2*(r_cmb**2/rZ(nZ,nS)-r_cmb)))
                   orhoZ(nZ,ns)=1.d0/rhoZ(nZ,nS)
                END DO
             ELSE
                DO nZ=1,nZmax
                   rhoZ(nZ,nS)=1.d0
                   orhoZ(nZ,ns)=1.d0
                END DO
             END IF
          END IF

          !--------- Get the flow components for all northern and
          !          southern thetas and all phis:
          IF ( lTC ) THEN
             nZmax=2*nZmaxS(nS) ! north and south points
             ! calculated in one go
          ELSE
             nZmax=nZmaxS(nS)
          END IF

          CALL getDVptr(wS,dwS,ddwS,zS,dzS,r_ICB,r_CMB,rZ(1,nS),       &
               &      nZmax,nZmaxA,PlmS(1,1,nS),dPlmS(1,1,nS),OsinTS(1,nS),       &
               &                           lDeriv,VrS,VtS,VpS,VozS,dpEkInt)

          nZmax=nZmaxS(nS)

          !------- Perform z-integral(s) for all phis:
          DO nPhi=1,n_phi_max

             !------- Two seperate integrals inside TC:
             IF ( lTC ) THEN
                nInts=2 ! separate north and south integral
             ELSE
                nInts=1
             END IF
             DO nInt=1,nInts
                IF ( nInt.EQ.1 ) THEN
                   DO nZ=1,nZmax ! Copy on simpler array
                      VrInt(nZ)=orhoZ(nZ,nS)*VrS(nPhi,nZ)
                      VtInt(nZ)=orhoZ(nZ,nS)*VtS(nPhi,nZ)
                      VpInt(nZ)=orhoZ(nZ,nS)*VpS(nPhi,nZ)
                      EkInt(nZ)=                                       &
                           &               (VrInt(nZ)**2+VtInt(nZ)**2+VpInt(nZ)**2)
                   END DO
                ELSE IF ( nInt.EQ.2 ) THEN
                   help=zMax
                   zMax=-zMin
                   zMin=-help
                   DO nZ=1,nZmax
                      VrInt(nZ)=orhoZ(nZ,nS)*VrS(nPhi,nZ+nZmax)
                      VtInt(nZ)=orhoZ(nZ,nS)*VtS(nPhi,nZ+nZmax)
                      VpInt(nZ)=orhoZ(nZ,nS)*VpS(nPhi,nZ+nZmax)
                      EkInt(nZ)=                                       &
                           &               (VrInt(nZ)**2+VtInt(nZ)**2+VpInt(nZ)**2)
                   END DO
                END IF

                !-------- NOTE: chebIntD replaces VrInt with z-derivative
                !               for lDeriv=.TRUE.
                VrIntS=chebIntD(VrInt,lDeriv,zMin,zMax,nZmax,nZmaxA,   &
                     &                       i_costf_initZ(1,nS),d_costf_initZ(1,nS))
                VtIntS=chebIntD(VtInt,lDeriv,zMin,zMax,nZmax,nZmaxA,   &
                     &                       i_costf_initZ(1,nS),d_costf_initZ(1,nS))
                VpIntS=chebIntD(VpInt,lDeriv,zMin,zMax,nZmax,nZmaxA,   &
                     &                       i_costf_initZ(1,nS),d_costf_initZ(1,nS))
                EkIntS=chebIntD(EkInt,.FALSE.,zMin,zMax,nZmax,nZmaxA,  &
                     &                        i_costf_initZ(1,nS),d_costf_initZ(1,nS))

                !-------- Get volume integral of energies:
                Egeos_s(nS)=Egeos_s(nS) +                              &
                     &                  0.5D0*phiNorm*(zMax-zMin)*sZ(nS)*dsZ *          &
                     &                          (VrIntS**2+VtIntS**2+VpIntS**2)
                IF ( lTC ) THEN
                   IF ( nInt.EQ.1 ) THEN
                      EkSTC_s(nS)=EkSTC_s(nS) +                        &
                           &                 0.5D0*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                   ELSE IF ( nInt.EQ.2 ) THEN
                      EkNTC_s(nS)=EkNTC_s(nS) +                        &
                           &                 0.5D0*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                   END IF
                ELSE
                   EkOTC_s(nS)=EkOTC_s(nS) +                           &
                        &                 0.5D0*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                END IF

                !-------- Note: chebIntD returns the z derivative
                IF ( lDeriv ) THEN
                   DO nZ=1,nZmax
                      dzEkInt(nZ)=VrInt(nZ)*VrInt(nZ) +                &
                           &                 VtInt(nZ)*VtInt(nZ) + VpInt(nZ)*VpInt(nZ)
                   END DO
                   dzEkIntS=chebInt(dzEkInt,zMin,zMax,nZmax,nZmaxA,    &
                        &                   i_costf_initZ(1,nS),d_costf_initZ(1,nS))
                   dzEk_s(nS)=dzEk_s(nS) +                             &
                        &                 0.5D0*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*dzEkIntS
                END IF

             END DO ! Loop over north/south integral


             ! --- All the stuff for North/South correlation only outside TC
             !     and only 0.1 away from boundaries:
             !     Here only integrals over one hemisphere are required. I thus
             !     copy the southern points to the northern hermisphere.
             IF ( lCorrel .AND. .NOT.lTC ) THEN
                VzSN =0.D0
                VzSS =0.D0
                VzNN =0.D0
                VorSN=0.D0
                VorSS=0.D0
                VorNN=0.D0
                HelZZ=0.D0
                VZZ  =0.D0
                VorZZ=0.D0
                DO nZN=1,nZmax/2 ! Dont use equatorial point
                   nZS  =nZmax-nZN+1
                   delz =zZ(nZN,nS)-zZ(nZN+1,nS)
                   sinT =1.D0/OsinTS(nZN,nS)
                   cosT =DSQRT(1.D0-sinT**2)
                   VzN  = cosT*VrS(nPhi,nZN)-sinT*VtS(nPhi,nZN)
                   VzS  =-cosT*VrS(nPhi,nZS)-sinT*VtS(nPhi,nZS)
                   VorN =VozS(nPhi,nZN)
                   VorS =VozS(nPhi,nZS)
                   VzSN =VzSN+delz*VzS*VzN
                   VzNN =VzNN+delz*VzN*VzN
                   VzSS =VzSS+delz*VzS*VzS
                   VorSN=VorSN+delz*VorS*VorN
                   VorNN=VorNN+delz*VorN*VorN
                   VorSS=VorSS+delz*VorS*VorS
                   HelZZ=HelZZ+delz*(VorN*VzN+VorS*VzS)**2
                   VZZ  =VZZ  +delz*(VzN*VzN+VzS*VzS)
                   VorZZ=VorZZ+delz*(VorN*VorN+VorS*VorS)
                END DO

                !------------ Calculate the integrals from zMin+0.1 to zMax-0.1:
                CVz_I = VzSN/DSQRT(VzSS*VzNN)
                CVor_I=VorSN/DSQRT(VorSS*VorNN)
                CHel_I=DSQRT(HelZZ/(VZZ*VorZZ))
                CVz(nPhi,nSmax-nS+1) =SNGL(CVz_I)
                CVor(nPhi,nSmax-nS+1)=SNGL(CVor_I)   
                CHel(nPhi,nSmax-nS+1)=SNGL(CHel_I)   
                CVz_s(nS) =CVz_s(nS) +phiNorm*sZ(nS)*dsZ*CVz_I  
                CVor_s(nS)=CVor_s(nS)+phiNorm*sZ(nS)*dsZ*CVor_I
                CHel_s(nS)=CHel_s(nS)+phiNorm*sZ(nS)*dsZ*CHel_I

             END IF ! lCorrel ?

          END DO  ! Loop over phi

          IF ( lDeriv ) THEN
             !--------- dpEkInt treated differently cause phi intergral has 
             !          already been preformed by getDVptr.f
             IF ( lTC ) THEN
                nInts=2 ! separate north and south integral
             ELSE
                nInts=1
             END IF
             DO nInt=1,nInts
                IF ( nInt.EQ.2 ) THEN
                   help=zMax
                   zMax=-zMin
                   zMin=-help
                   DO nZ=1,nZmax
                      dpEkInt(nZ)=dpEkInt(nZ+nZmax)
                   END DO
                END IF
                dpEkIntS=chebInt(dpEkInt,zMin,zMax,nZmax,nZmaxA,       &
                     &                   i_costf_initZ(1,nS),d_costf_initZ(1,nS))
                dpEk_s(nS)=dpEk_s(nS) +                                &
                     &                     0.5D0*(zMax-zMin)*sZ(nS)*dsZ*dpEkIntS /      &
                     &                     (sZ(nS)**2) ! Convert angle to length
             END DO
          END IF

       END DO ! big loop over S

       !--- Collect:
       EkSTC  =0.D0
       EkNTC  =0.D0
       EkOTC  =0.D0
       Egeos  =0.D0
       CVzOTC =0.D0
       CVorOTC=0.D0
       CHelOTC=0.D0
       DO nS=1,nSmax
          Egeos=Egeos+Egeos_s(nS)
          IF ( sZ(ns).LT.r_ICB ) THEN
             EkSTC=EkSTC+EkSTC_s(nS)
             EkNTC=EkNTC+EkNTC_s(nS)
          ELSE
             EkOTC  =EkOTC  +EkOTC_s(nS)
             IF ( lCorrel ) THEN
                CVzOTC =CVzOTC +CVz_s(nS)
                CVorOTC=CVorOTC+CVor_s(nS)
                CHelOTC=CHelOTC+CHel_s(nS)
             END IF
          END IF
       END DO
       IF ( lCorrel ) THEN
          surf=0.D0
          DO nS=1,nSmax
             IF ( sZ(nS).GE.r_ICB ) surf=surf+sZ(nS)*dsZ
          END DO
          surf   =2.D0*pi*surf
          CVzOTC =CVzOTC/surf
          CVorOTC=CVorOTC/surf
          CHelOTC=CHelOTC/surf
       END IF
       Ekin=EkSTC+EkNTC+EKOTC ! Can be used for testing 

       dpFlow=0.D0
       dzFlow=0.D0
       IF ( lDeriv ) THEN
          DO nS=1,nSmax
             dpFlow=dpFlow+dpEk_s(nS)
             dzFlow=dzFlow+dzEk_s(nS)
          END DO
          dpFlow=DSQRT(Ekin/dpFlow)
          dzFlow=DSQRT(Ekin/dzFlow)
       END IF

       !        WRITE(99,*) 'E:',EkSTC,EkNTC,EkOTC
       !        WRITE(99,*) 'Ekin:',Ekin
       !        WRITE(99,*) 'Egeos:',Egeos

       !--- Write correlation movie:
       IF ( l_corrMov ) THEN

          !--- Determine s used for correl
          n=0
          DO nS=1,nSmax
             IF ( sZ(nS).GE.r_ICB+0.1 .AND.                            &
                  &             sZ(nS).LE.r_CMB-0.1 ) THEN
                n=n+1
             ELSE 
                DO nPhi=1,n_phi_max
                   CVz(nPhi,nSmax-nS+1)=0.E0
                   CVor(nPhi,nSmax-nS+1)=0.E0
                END DO
             END IF
          END DO
          !           ndS=nSmax-n
          !          nSmax=n

          nOutFile=93
          movFile ='CVorz_mov.'//tag
          OPEN(nOutFile,file=movFile,STATUS='UNKNOWN',                 &
               &          FORM='UNFORMATTED',POSITION='APPEND')

          !--- Write header into output file:
          IF ( nGeosSets.EQ.1 ) THEN

             nFields=3
             nFieldSize=(nSmax-nS_ICB+1)*n_phi_max
             version='JW_Movie_Version_2'
             WRITE(nOutFile) version
             dumm(1)=111           ! type of input
             dumm(2)=2             ! marker for constant theta plane
             dumm(3)=90.D0         ! surface constant
             dumm(4)=nFields       ! no of fields
             WRITE(nOutFile) (real(dumm(n),4),n=1,4)

             dumm(1)=92.0          ! Field marker for AS vPhi
             dumm(2)=93.0          ! Field marker for Reynolds Force
             dumm(3)=94.0          ! Field marker for Reynolds Force
             WRITE(nOutFile) (SNGL(dumm(n)),n=1,nFields)

             !------ Now other info about grid and parameters:
             WRITE(nOutFile) runid     ! run identifier
             dumm( 1)=nSmax            ! total number of radial points
             dumm( 2)=nSmax            ! no of radial point in outer core
             dumm( 3)=1                ! no. of theta points
             dumm( 4)=n_phi_max        ! no. of phi points
             dumm( 5)=minc             ! imposed symmetry
             dumm( 6)=ra               ! control parameters
             dumm( 7)=ek               ! (for information only)
             dumm( 8)=pr               !      -"-
             dumm( 9)=prmag            !      -"-
             dumm(10)=radratio         ! ratio of inner / outer core
             dumm(11)=tScale           ! timescale
             WRITE(nOutFile) (SNGL(dumm(n)),   n=1,11)
             WRITE(nOutFile) (SNGL(sZ(nSmax-n+1)/r_CMB),n=1,nSmax)
             WRITE(nOutFile)  SNGL(90.D0)
             WRITE(nOutFile) (SNGL(phi(n)), n=1,n_phi_max)

          END IF ! Write Header ?

          dumm(1)=nGeosSets          ! time frame number for movie
          dumm(2)=time              ! time
          dumm(3)=0.D0
          dumm(4)=0.D0
          dumm(5)=0.D0
          dumm(6)=0.D0
          dumm(7)=0.D0
          dumm(8)=0.D0
          WRITE(nOutFile) (SNGL(dumm(n)),n=1,8)

          WRITE(nOutFile) ((CVz(nPhi,nS) ,nPhi=1,n_phi_max),           &
               &                       nS=1,nSmax)
          WRITE(nOutFile) ((CVor(nPhi,nS),nPhi=1,n_phi_max),           &
               &                       nS=1,nSmax)
          WRITE(nOutFile) ((CHel(nPhi,nS),nPhi=1,n_phi_max),           &
               &                       nS=1,nSmax)

       END IF

    END IF
    IF ( lVerbose ) WRITE(*,*) '! End of getGeos!'


    RETURN
  END SUBROUTINE getEgeos

  !-----------------------------------------------------------------------------
END MODULE Egeos_mod
