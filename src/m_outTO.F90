!$Id$
MODULE outTO_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE torsional_oscillations
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE const

  IMPLICIT NONE 
  
  INTEGER :: lmMaxS
  INTEGER :: nZmaxA,nZmaxL
  INTEGER :: nSmaxA,nSmaxL

  !-- Plms: Plm,sin
  REAL(kind=8),ALLOCATABLE :: PlmS(:,:,:)
  REAL(kind=8),ALLOCATABLE :: dPlmS(:,:,:)
  REAL(kind=8),ALLOCATABLE :: OsinTS(:,:)
  REAL(kind=4),ALLOCATABLE :: VpM(:,:),LFM(:,:),dVpM(:,:),AstrM(:,:)
  REAL(kind=4),ALLOCATABLE :: RstrM(:,:),CorM(:,:),StrM(:,:),CLM(:,:)
  INTEGER,ALLOCATABLE :: i_costf_initZ(:,:),nZmaxS(:)
  REAL(kind=8),ALLOCATABLE :: d_costf_initZ(:,:),zZ(:,:),rZ(:,:)


contains
  SUBROUTINE initialize_outTO_mod
    nZmaxL=lStressMem*722
    nZmaxA=MAX0(2,nZmaxL)
    nSmaxL=lStressMem*625
    nSmaxA=MAX0(3,nSmaxL)
    lmMaxS = l_max+1
    ALLOCATE( PlmS(lmMaxS,nZmaxA/2+1,nSmaxA) )
    ALLOCATE( dPlmS(lmMaxS,nZmaxA/2+1,nSmaxA) )
    ALLOCATE( OsinTS(nZmaxA/2+1,nSmaxA) )
    ALLOCATE( vpM(nZmaxA/2,nSmaxA) )
    ALLOCATE( LFM(nZmaxA/2,nSmaxA) )
    ALLOCATE( dVpM(nZmaxA/2,nSmaxA) )
    ALLOCATE( AstrM(nZmaxA/2,nSmaxA) )
    ALLOCATE( RstrM(nZmaxA/2,nSmaxA) )
    ALLOCATE( CorM(nZmaxA/2,nSmaxA) )
    ALLOCATE( StrM(nZmaxA/2,nSmaxA) )
    ALLOCATE( CLM(nZmaxA/2,nSmaxA) )
    ALLOCATE( i_costf_initZ(2*nZmaxA+2,nSmaxA) )
    ALLOCATE( d_costf_initZ(2*nZmaxA+5,nSmaxA) )
    ALLOCATE( zZ(nZmaxA,nSmaxA) )
    ALLOCATE( rZ(nZmaxA/2+1,nSmaxA) )
    ALLOCATE( nZmaxS(nSmaxA) )

  END SUBROUTINE initialize_outTO_mod

  !***********************************************************************
  SUBROUTINE outTO(time,n_time_step,                                      &
       &                     eKin,eKinTAS,nOutFile,nOutFile2,             &
       &                 TOfileNhs,TOfileShs,movFile,tayFile,             &
       &                       nTOsets,nTOmovSets,nTOrmsSets,             &
       &                             lTOmov,lTOrms,lTOZwrite,             &
       &                                 z,omega_ic,omega_ma)
    !***********************************************************************

    !-----------------------------------------------------------------------
    !   Output of axisymmetric zonal flow, its relative strength,
    !   its time variation, and all forces acting on it.
    !   The slowest part in the TO process is the repetitious calculation
    !   of plms by subroutine plm_theta. They are needed in getAStr and 
    !   getPAStr when I transform on the cylindrical grid. 
    !   The necessary plms could simply be calculated one and then 
    !   be stored for later use! See s_outTOnew.f.
    !-----------------------------------------------------------------------

    !-- Input of variables:
    USE charmanip, ONLY: dble2str
    USE integration, ONLY: rInt_R
    USE plms_theta, ONLY: plm_theta

    REAL(kind=8) :: time
    INTEGER :: n_time_step
    REAL(kind=8) :: eKin,eKinTAS
    INTEGER :: nOutFile,nOutFile2
    CHARACTER(len=*) :: TOfileNhs,TOfileShs,movFile,tayFile
    INTEGER :: nTOsets,nTOmovSets,nTOrmsSets
    LOGICAL :: lTOmov,lTOrms   
    LOGICAL :: lTOZwrite

    !-- Input of stuff calculated in getTO and getTOfinish
    ! include 'c_TO.f'

    !-- Input of toroidal flow scalar and rotation rates:
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    REAL(kind=8) :: omega_ic,omega_ma

    !-- Output field:
    REAL(kind=8) :: fOut(n_theta_maxStr*n_r_maxStr)

    !-- Local:
    !REAL(kind=8) dzStrMA,dzStrIC
    !REAL(kind=8) dzLFMA,dzLFIC
    LOGICAL :: lTC,lStopRun

    !-- (l,r) Representation of the different contributions
    REAL(kind=8) :: dzVpLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: V2LMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: Bs2LMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BszLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BspLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BpzLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BspdLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BpsdLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BzpdLMr(lmMaxS,n_r_maxStr)
    REAL(kind=8) :: BpzdLMr(lmMaxS,n_r_maxStr)

    !---- Work array:
    COMPLEX(kind=8) :: workA(lmMaxS,n_r_maxStr)

    INTEGER :: lm,l ! counter for degree

    INTEGER :: nSmax,nS,nSI
    REAL(kind=8) :: zNorm  ! Norm z intervall
    INTEGER :: nNorm  ! No. of grid points for norm intervall
    REAL(kind=8) :: zMin,zMax!,help ! integration boundarie, help variable
    LOGICAL :: lAS    ! .TRUE. if axisymmetric (m=0) functions
    REAL(kind=8) :: sZ(nSmaxA),dsZ ! cylindrical radius s and s-step
    REAL(kind=8) :: h(nSmaxA),Oh(nSmaxA)
    REAL(kind=8) :: Os2(nSmaxA)
    !REAL(kind=8)  norm

    INTEGER :: nR     ! counter for radial grid point
    INTEGER :: n      ! counter for theta blocks
    !INTEGER :: nNS    ! index for NHS and SHS
    INTEGER :: nOut,nFields ! counter for output fields
    INTEGER :: nTheta ! counter for all thetas
    INTEGER :: nThetaBlock ! counter for thetas in block
    INTEGER :: nThetaOrd ! counter for ordered thetas
    INTEGER :: nThetaNHS
    INTEGER :: nThetaStart

    INTEGER :: nZ,nZmax,nZmaxNS,nZmaxH!,nZP
    REAL(kind=8) :: VpS(nZmaxA)      
    REAL(kind=8) :: dVpS(nZmaxA)      
    REAL(kind=8) :: ddVpS(nZmaxA)      
    REAL(kind=8) :: V2S(nZmaxA)
    REAL(kind=8) :: LFS(nZmaxA)
    REAL(kind=8) :: CorS(nZmaxA)
    REAL(kind=8) :: RstrS(nZmaxA)
    REAL(kind=8) :: AstrS(nZmaxA)
    REAL(kind=8) :: StrS(nZmaxA)
    REAL(kind=8) :: Bs2S(nZmaxA)
    REAL(kind=8) :: BspS(nZmaxA)
    REAL(kind=8) :: BspdS(nZmaxA)
    REAL(kind=8) :: BpsdS(nZmaxA)
    REAL(kind=8) :: TayS(nZmaxA)
    REAL(kind=8) :: TayRS(nZmaxA)
    REAL(kind=8) :: TayVS(nZmaxA)

    REAL(kind=8) :: VpIntN(nSmaxA)  ,VpIntS(nSmaxA)    ! integration results
    REAL(kind=8) :: dVpIntN(nSmaxA) ,dVpIntS(nSmaxA)   ! integration results
    REAL(kind=8) :: ddVpIntN(nSmaxA),ddVpIntS(nSmaxA)  ! integration results
    REAL(kind=8) :: VpRIntN(nSmaxA) ,VpRIntS(nSmaxA)   ! for different s and 
    REAL(kind=8) :: V2IntS(nSmaxA)  ,V2IntN(nSmaxA)
    REAL(kind=8) :: LFIntN(nSmaxA)  ,LFIntS(nSmaxA)   
    REAL(kind=8) :: RstrIntN(nSmaxA),RstrIntS(nSmaxA) 
    REAL(kind=8) :: AstrIntN(nSmaxA),AstrIntS(nSmaxA) 
    REAL(kind=8) :: StrIntN(nSmaxA) ,StrIntS(nSmaxA) 
    REAL(kind=8) :: Bs2IntN(nSmaxA) ,Bs2IntS(nSmaxA)
    REAL(kind=8) :: BspIntN(nSmaxA) ,BspIntS(nSmaxA)
    REAL(kind=8) :: BspdIntN(nSmaxA),BspdIntS(nSmaxA)
    REAL(kind=8) :: BpsdIntN(nSmaxA),BpsdIntS(nSmaxA)
    REAL(kind=8) :: TayIntN(nSmaxA) ,TayIntS(nSmaxA) 
    REAL(kind=8) :: TayRIntN(nSmaxA),TayRIntS(nSmaxA) 
    REAL(kind=8) :: TayVIntN(nSmaxA),TayVIntS(nSmaxA) 
    REAL(kind=8) :: SVpIntN(nSmaxA) ,SVpIntS(nSmaxA)   ! help arrays and values for 
    REAL(kind=8) :: SBs2IntN(nSmaxA),SBs2IntS(nSmaxA)  ! differentiation in s   
    REAL(kind=8) :: SBspIntN(nSmaxA),SBspIntS(nSmaxA)
    REAL(kind=8) :: dSVpIntN, dSVpIntS
    REAL(kind=8) :: d2SVpIntN,d2SVpIntS
    REAL(kind=8) :: dSBspIntN,dSBspIntS
    REAL(kind=8) :: dSBs2IntN,dSBs2IntS
    REAL(kind=8) :: TauN(nSmaxA),TauS(nSmaxA)          ! Taylor integral
    REAL(kind=8) :: TauBN(nSmaxA),TauBS(nSmaxA)       
    REAL(kind=8) :: dTauBN(nSmaxA),dTauBS(nSmaxA)    
    REAL(kind=8) :: dTTauN(nSmaxA),dTTauS(nSmaxA)      ! time change of Tau...
    REAL(kind=8) :: dTTauBN(nSmaxA),dTTauBS(nSmaxA)   

    !-- For integration along z:
    REAL(kind=8) :: zALL(2*nZmaxA)
    REAL(kind=8) :: thetaZ
    REAL(kind=8) :: fac
    REAL(kind=8) :: vSF,fSF,f1,f2

    !-- For boundaries:
    REAL(kind=8) :: rB(2),BspB(2),BspdB(2),BpsdB(2)
    REAL(kind=8) :: Bs2B(2),BszB(2),BpzB(2),BzpdB(2),BpzdB(2)

    !-- For sphere integration:
    REAL(kind=8) :: CorR(n_r_max)
    REAL(kind=8) :: RstrR(n_r_max)
    REAL(kind=8) :: AstrR(n_r_max)
    REAL(kind=8) :: StrR(n_r_max)
    REAL(kind=8) :: VpR(n_r_max),dVpR(n_r_max)
    REAL(kind=8) :: LFR(n_r_max),LFABSR(n_r_max)
    REAL(kind=8) :: TayR(n_r_max),TayRMSR(n_r_max)
    REAL(kind=8) :: TayRMS,TaySRMS
    REAL(kind=8) :: TayRRMS,TayVRMS
    REAL(kind=8) :: VpRMS,VRMS,VgRMS
    REAL(kind=8) :: rS,sS
    REAL(kind=8) :: outBlock(nfs)
    REAL(kind=8) :: timeLast!,tNorm
    REAL(kind=4) :: timeAve,dt
    SAVE timeLast,timeAve

    CHARACTER(len=64) :: version,fileName,fileZ
    INTEGER :: nFieldSize,nPos

    !INTEGER :: nLines
    REAL(kind=8),EXTERNAL :: chebInt

    REAL(kind=8) :: dumm(12)

    !-- Huge arrays for time average ....
    LOGICAL :: l_TOZave

    !-- For TOZ output files:
    INTEGER,SAVE :: nTOZfile!,length
    CHARACTER(len=10) :: string

    !-- End of declaration
    !---------------------------------------------------------------------------

    IF ( lVerbose ) WRITE(*,*) '! Starting outTO!'

    nTOsets=nTOsets+1

    l_TOZave=.TRUE.

    !--- Rescaling for rotation time scale and planetary radius
    !    length scale, for the velocity I use the Rossby number
    !       vSF=ek/r_CMB                   
    !       fSF=ek*ek/(4.D0*pi**2*r_CMB)  
    vSF=1.D0
    fSF=1.D0

    nFieldSize=n_theta_maxStr*n_r_maxStr

    !-- Start with calculating advection due to axisymmetric flows:

    zNorm=1.D0               ! This is r_CMB-r_ICB
    nNorm=INT(zDens*n_r_max) ! Covered with nNorm  points !
    nSmax=n_r_max+INT(r_ICB*DBLE(n_r_max))
    nSmax=INT(sDens*nSmax)
    IF ( nSmax.GT.nSmaxA ) THEN
       WRITE(*,*) 'Increase nSmaxA in ouTO!'
       WRITE(*,*) 'Should be at least nSmax=',nSmax
       STOP
    END IF
    lAS=.TRUE.

!!! TEST
    !       DO nR=1,n_r_max
    !       DO n=1,n_theta_max
    !          V2AS(n,nR)=n+nR
    !       END DO
    !       END DO
!!!

    !--- Transform to lm-space for all radial grid points:

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)                                   &
    !$OMP  PRIVATE(nR,n,nThetaStart)
    DO nR=1,n_r_max
       DO n=1,nThetaBs
          nThetaStart=(n-1)*sizeThetaB+1
          CALL legtfAS(V2LMr(1,nR),V2AS(nThetaStart,nR),            &
               &               l_max+1,nThetaStart,sizeThetaB)
          CALL legtfAS2(Bs2LMr(1,nR),BszLMr(1,nR),                  &
               &               Bs2AS(nThetaStart,nR),BszAS(nThetaStart,nR),       &
               &               l_max+1,nThetaStart,sizeThetaB)
          CALL legtfAS2(BspLMr(1,nR),BpzLMr(1,nR),                  &
               &               BspAS(nThetaStart,nR),BpzAS(nThetaStart,nR),       &
               &               l_max+1,nThetaStart,sizeThetaB)
          CALL legtfAS2(BspdLMr(1,nR),BpsdLMr(1,nR),                &
               &               BspdAS(nThetaStart,nR),BpsdAS(nThetaStart,nR),     &
               &               l_max+1,nThetaStart,sizeThetaB)
          CALL legtfAS2(BzpdLMr(1,nR),BpzdLMr(1,nR),                &
               &               BzpdAS(nThetaStart,nR),BpzdAS(nThetaStart,nR),     &
               &               l_max+1,nThetaStart,sizeThetaB)
       END DO
    END DO ! Loop over radial grid points
    !$OMP END PARALLEL DO

    DO nR=1,n_r_max
       DO l=1,l_max
          lm=lm2(l,0)
          dzVpLMr(l+1,nR)=REAL(z(lm,nR))
       END DO
    END DO

    !---- Transform the contributions to cheb space for z-integral:
    CALL costf1(dzVpLMr,lmMaxS,1,lmMaxS,                            &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(V2LMr,lmMaxS,1,lmMaxS,                              &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzdVpLMr,lmMaxS,1,lmMaxS,                           &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzddVpLMr,lmMaxS,1,lmMaxS,                          &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(Bs2LMr,lmMaxS,1,lmMaxS,                             &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BszLMr,lmMaxS,1,lmMaxS,                             &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BspLMr,lmMaxS,1,lmMaxS,                             &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BpzLMr,lmMaxS,1,lmMaxS,                             &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzRstrLMr,lmMaxS,1,lmMaxS,                          &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzAstrLMr,lmMaxS,1,lmMaxS,                          &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzStrLMr,lmMaxS,1,lmMaxS,                           &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzLFLMr,lmMaxS,1,lmMaxS,                            &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(dzCorLMr,lmMaxS,1,lmMaxS,                           &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BspdLMr,lmMaxS,1,lmMaxS,                            &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BpsdLMr,lmMaxS,1,lmMaxS,                            &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BpzdLMr,lmMaxS,1,lmMaxS,                            &
         &              workA,i_costf_init,d_costf_init)
    CALL costf1(BzpdLMr,lmMaxS,1,lmMaxS,                            &
         &              workA,i_costf_init,d_costf_init)

    dsZ   =r_CMB/DBLE(nSmax)  ! Step in s controlled by nSmax
    nSI   =0
    DO nS=1,nSmax
       sZ(nS)=(nS-0.5D0)*dsZ
       IF ( sZ(nS).LT.r_ICB .AND. nS.GT.nSI ) nSI=nS
    END DO

    IF ( nTOsets.EQ.1 ) nTOZfile=0
    IF ( lTOZwrite ) THEN
       nTOZfile=nTOZfile+1
       CALL dble2str(DBLE(nTOZfile),string)
       fileName='TOZ_'//trim(adjustl(string))//'.'//TAG
       OPEN(95,FILE=fileName,FORM='UNFORMATTED',                    &
            &                             STATUS='UNKNOWN')
       WRITE(95) SNGL(time),FLOAT(nSmax),                           &
            &         SNGL(omega_ic),SNGL(omega_ma)
       WRITE(95) (SNGL(sZ(nS)),nS=1,nSmax)
    END IF
    IF ( nTOsets.GT.1 .AND. l_TOZave ) THEN
       fileName='TOZM.'//TAG
       OPEN(96,FILE=fileName,FORM='UNFORMATTED',                    &
            &                                STATUS='UNKNOWN')
       WRITE(96) FLOAT(nSmax),SNGL(omega_ic),SNGL(omega_ma)
       WRITE(96) (SNGL(sZ(nS)),nS=1,nSmax)
    END IF

    lStopRun=.FALSE.
    !$OMP DO ORDERED                                                        &
    !$OMP  PRIVATE(zMax,zMin,lTC,nZmax,nZ,nZmaxNS,nZmaxH,thetaZ,fileZ,      &
    !$OMP             VpS,dVpS,ddVpS,RstrS,AstrS,StrS,LFS,CorS,TayS,        &
    !$OMP                     TayRS,TayVS,V2S,Bs2S,BspS,BspdS,BpsdS,        &
    !$OMP            rB,BspB,BspdB,BpsdB,Bs2B,BszB,BpzB,BzpdB,BpzdB)

    DO nS=1,nSmax

       !------ Get integral boundaries for this s:
       zMax=DSQRT(r_CMB*r_CMB-sZ(nS)*sZ(nS))
       IF ( sZ(nS).LT.r_ICB ) THEN
          lTC=.TRUE.
          zMin=DSQRT(r_ICB*r_ICB-sZ(nS)*sZ(nS))
       ELSE
          lTC=.FALSE.
          zMin=-zMax
       ENDIF
       h(nS) =zMax-zMin
       Oh(nS)=1.D0/h(nS)

       IF ( nTOsets.EQ.1 ) THEN
          !------ Initialize integration for NHS:
          !       Each processor calculates Cheb transform data
          !       for HIS nS and the Plms along the Cylinder
          !       chebIntInit returns zZ,nZmaxS,i_costf_initZ and
          !       d_costfInitZ:
          !       Note that this returns z in the MAGIC way, 
          !       starting with zMax, ending with zMin
          !       z(1,nS)=zMin, z(nZmax,nS)=zMax
          CALL chebIntInit(zMin,zMax,zNorm,nNorm,                   &
               &                    nZmaxA,zZ(1,nS),nZmaxS(nS),                   &
               &       i_costf_initZ(1,nS),d_costf_initZ(1,nS))

          !--- Points in nothers halfsphere
          IF ( lTC ) THEN
             nZmax=nZmaxS(nS)  ! nZmax point in each polar region
             IF ( 2*nZmax.GT.nZmaxA ) THEN 
                WRITE(*,*) '! nZmaxA too small in outTO!'
                WRITE(*,*) '! Should be at least:',2*nZmax
                lStopRun=.TRUE.
             END IF
          ELSE
             nZmax=(nZmaxS(nS)-1)/2+1 ! Odd point number !
             ! all together nZmaxS(nS) from
             ! south to north including equator
             IF ( nZmaxS(nS).GT.nZmaxA ) THEN 
                WRITE(*,*) '! nZmaxA too small in outTO!'
                WRITE(*,*) '! Should be at least:',nZmaxS(nS)
                lStopRun=.TRUE.
             END IF
          END IF
          IF ( lStopRun ) GOTO 99
          DO nZ=1,nZmax
             rZ(nZ,nS)    =DSQRT(zZ(nZ,nS)**2+sZ(nS)**2)
             thetaZ       =DATAN2(sZ(nS),zZ(nZ,nS))
             OsinTS(nZ,nS)=1.D0/DSIN(thetaZ)
             CALL plm_theta(thetaZ,l_max,0,minc,                    &
                  &            PlmS(1,nZ,nS),dPlmS(1,nZ,nS),lmMaxS,2)
          END DO
       END IF

       !--------- Get the flow components for all northern and
       !          southern thetas and all phis:

       nZmax=nZmaxS(nS)
       IF ( lTC ) THEN
          nZmaxNS=2*nZmax
          nZmaxH =nZmax
          DO nZ=1,nZmax
             zALL(nZ)=zZ(nZ,nS)
             zALL(nZmaxNS-nZ+1)=-zZ(nZ,nS)
          END DO
       ELSE
          nZmaxNS=nZmax
          nZmaxH =(nZmax-1)/2+1
          DO nZ=1,nZmax
             zALL(nZ)=zZ(nZ,nS)
          END DO
       END IF

       CALL getPAStr(VpS,dzVpLMr,nZmaxNS,nZmaxA,lmMaxS,             &
            &                      l_max,r_ICB,r_CMB,n_r_max,             &
            &                 rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(dVpS,dzdVpLMr,nZmaxNS,nZmaxA,lmMaxS,           &
            &                        l_max,r_ICB,r_CMB,n_r_max,           &
            &                   rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(ddVpS,dzddVpLMr,nZmaxNS,nZmaxA,lmMaxS,         &
            &                          l_max,r_ICB,r_CMB,n_r_max,         &
            &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(RstrS,dzRstrLMr,nZmaxNS,nZmaxA,lmMaxS,         &
            &                          l_max,r_ICB,r_CMB,n_r_max,         &
            &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(AstrS,dzAstrLMr,nZmaxNS,nZmaxA,lmMaxS,         &
            &                          l_max,r_ICB,r_CMB,n_r_max,         &
            &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(StrS,dzStrLMr,nZmaxNS,nZmaxA,lmMaxS,           &
            &                        l_max,r_ICB,r_CMB,n_r_max,           &
            &                   rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(LFS,dzLFLMr,nZmaxNS,nZmaxA,lmMaxS,             &
            &                      l_max,r_ICB,r_CMB,n_r_max,             &
            &                 rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       CALL getPAStr(CorS,dzCorLMr,nZmaxNS,nZmaxA,lmMaxS,           &
            &                        l_max,r_ICB,r_CMB,n_r_max,           &
            &                   rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
       DO nZ=1,nZmaxNS
          TayS(nZ) =DABS(LFS(nZ))
          TayRS(nZ)=DABS(RstrS(nZ))
          TayVS(nZ)=DABS(StrS(nZ))
       END DO

       CALL getAStr(V2S,V2LMr,nZmaxNS,nZmaxA,lmMaxS,                &
            &                   l_max,r_ICB,r_CMB,n_r_max,                &
            &                            rZ(1,nS),PlmS(1,1,nS))
       CALL getAStr(Bs2S,Bs2LMr,nZmaxNS,nZmaxA,lmMaxS,              &
            &                     l_max,r_ICB,r_CMB,n_r_max,              &
            &                              rZ(1,nS),PlmS(1,1,nS))
       CALL getAStr(BspS,BspLMr,nZmaxNS,nZmaxA,lmMaxS,              &
            &                     l_max,r_ICB,r_CMB,n_r_max,              &
            &                              rZ(1,nS),PlmS(1,1,nS))
       CALL getAStr(BspdS,BspdLMr,nZmaxNS,nZmaxA,lmMaxS,            &
            &                       l_max,r_ICB,r_CMB,n_r_max,            &
            &                                rZ(1,nS),PlmS(1,1,nS))
       CALL getAStr(BpsdS,BpsdLMr,nZmaxNS,nZmaxA,lmMaxS,            &
            &                       l_max,r_ICB,r_CMB,n_r_max,            &
            &                                rZ(1,nS),PlmS(1,1,nS))

       IF ( l_TOZave ) THEN
          IF ( nTOsets.EQ.1 ) THEN
             timeAve=1.E0
             DO nZ=1,nZmaxNS
                VpM(nZ,nS)  =SNGL(VpS(nZ))
                dVpM(nZ,nS) =SNGL(dVpS(nZ))
                LFM(nZ,nS)  =SNGL(LFfac*LFS(nZ))
                RstrM(nZ,nS)=SNGL(RstrS(nZ))
                AstrM(nZ,nS)=SNGL(AstrS(nZ))
                StrM(nZ,nS) =SNGL(StrS(nZ))
                CorM(nZ,nS) =SNGL(CorS(nZ))
                CLM(nZ,nS)  =SNGL(CorS(nZ)+LFfac*LFS(nZ))
             END DO
          ELSE IF ( nTOsets.EQ.2 ) THEN
             dt=SNGL(time-timeLast)
             timeAve=dt
             DO nZ=1,nZmaxNS
                VpM(nZ,nS)  =dt*(VpM(nZ,nS)  +SNGL(VpS(nZ)))
                dVpM(nZ,nS) =dt*(dVpM(nZ,nS) +SNGL(dVpS(nZ)))
                LFM(nZ,nS)  =dt*(LFM(nZ,nS)  +SNGL(LFfac*LFS(nZ)))
                RstrM(nZ,nS)=dt*(RstrM(nZ,nS)+SNGL(RstrS(nZ)))
                AstrM(nZ,nS)=dt*(AstrM(nZ,nS)+SNGL(AstrS(nZ)))
                StrM(nZ,nS) =dt*(StrM(nZ,nS) +SNGL(StrS(nZ)))
                CorM(nZ,nS) =dt*(CorM(nZ,nS) +SNGL(CorS(nZ)))
                CLM(nZ,nS)  =dt*(CLM(nZ,nS)  +                         &
                     &                            SNGL(CorS(nZ)+LFfac*LFS(nZ)))
             END DO
          ELSE
             dt=SNGL(time-timeLast)
             timeAve=timeAve+dt
             DO nZ=1,nZmaxNS
                VpM(nZ,nS)  =VpM(nZ,nS)  +dt*SNGL(VpS(nZ))
                dVpM(nZ,nS) =dVpM(nZ,nS) +dt*SNGL(dVpS(nZ))
                LFM(nZ,nS)  =LFM(nZ,nS)  +dt*SNGL(LFfac*LFS(nZ))
                RstrM(nZ,nS)=RstrM(nZ,nS)+dt*SNGL(RstrS(nZ))
                AstrM(nZ,nS)=AstrM(nZ,nS)+dt*SNGL(AstrS(nZ))
                StrM(nZ,nS) =StrM(nZ,nS) +dt*SNGL(StrS(nZ))
                CorM(nZ,nS) =CorM(nZ,nS) +dt*SNGL(CorS(nZ))
                CLM(nZ,nS)  =CLM(nZ,nS)  +                             &
                     &                        dt*SNGL(CorS(nZ)+LFfac*LFS(nZ))
             END DO
          END IF
       END IF

       !$OMP ORDERED
       IF ( l_TOZave .AND. nTOsets.GT.1 ) THEN
          WRITE(96) FLOAT(nZmaxNS)
          WRITE(96) (SNGL(zALL(nZ))      ,nZ=1,nZmaxNS),            &
               &                  (VpM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS),            &
               &                  (dVpM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),            &
               &                  (RstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS),            &
               &                  (AstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS),            &
               &                  (LFM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS),            &
               &                  (StrM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),            &
               &                  (CorM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),            &
               &                  (CLM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS)
       END IF
       IF ( lTOZwrite ) THEN
          WRITE(95) FLOAT(nZmaxNS)
          WRITE(95) (SNGL(zALL(nZ)) ,nZ=1,nZmaxNS),                 &
               &                  (SNGL(VpS(nZ))  ,nZ=1,nZmaxNS),                 &
               &                  (SNGL(dVpS(nZ)) ,nZ=1,nZmaxNS),                 &
               &                  (SNGL(RstrS(nZ)),nZ=1,nZmaxNS),                 &
               &                  (SNGL(AstrS(nZ)),nZ=1,nZmaxNS),                 &
               &                  (SNGL(LFfac*LFS(nZ)),nZ=1,nZmaxNS),             &
               &                  (SNGL(StrS(nZ)) ,nZ=1,nZmaxNS),                 &
               &                  (SNGL(CorS(nZ)) ,nZ=1,nZmaxNS)
       END IF
       !$OMP END ORDERED

       !--- Z-integrals:
       VpIntN(nS)  =chebInt(VpS,zMin,zMax,nZmax,nZmaxA,             &
            &             i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       dVpIntN(nS) =chebInt(dVpS,zMin,zMax,nZmax,nZmaxA,            &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       ddVpIntN(nS)=chebInt(ddVpS,zMin,zMax,nZmax,nZmaxA,           &
            &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       LFIntN(nS)  =chebInt(LFS,zMin,zMax,nZmax,nZmaxA,             &
            &             i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       TayIntN(nS) =chebInt(TayS,zMin,zMax,nZmax,nZmaxA,            &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       TayRIntN(nS)=chebInt(TayRS,zMin,zMax,nZmax,nZmaxA,           &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       TayVIntN(nS)=chebInt(TayVS,zMin,zMax,nZmax,nZmaxA,           &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       RstrIntN(nS)=chebInt(RstrS,zMin,zMax,nZmax,nZmaxA,           &
            &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       AstrIntN(nS)=chebInt(AstrS,zMin,zMax,nZmax,nZmaxA,           &
            &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       StrIntN(nS) =chebInt(StrS,zMin,zMax,nZmax,nZmaxA,            &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       V2IntN(nS)  =chebInt(V2S,zMin,zMax,nZmax,nZmaxA,             &
            &             i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       Bs2IntN(nS) =chebInt(Bs2S,zMin,zMax,nZmax,nZmaxA,            &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       BspIntN(nS) =chebInt(BspS,zMin,zMax,nZmax,nZmaxA,            &
            &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       BspdIntN(nS)=chebInt(BspdS,zMin,zMax,nZmax,nZmaxA,           &
            &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
       BpsdIntN(nS)=chebInt(BpsdS,zMin,zMax,nZmax,nZmaxA,           &
            &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))

       IF ( V2IntN(nS).LT.0.D0 ) THEN
          VpRIntN(nS)=1.D0
       ELSE
          VpRIntN(nS)=DABS(VpIntN(nS))/DSQRT(V2IntN(nS))
       END IF
       VpRIntN(nS) =DMIN1(1.D0,VpRIntN(nS))
       TayIntN(nS) =LFIntN(nS)/TayIntN(nS)
       TayRIntN(nS)=RstrIntN(nS)/TayRIntN(nS)
       TayVIntN(nS)=StrIntN(nS)/TayVIntN(nS)

       !--- Z-integration inside northern TC:
       IF ( lTC ) THEN
          VpIntS(nS)  =chebInt(VpS(nZmax+1),zMin,zMax,nZmax,        &
               &           nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          dVpIntS(nS) =chebInt(dVpS(nZmax+1),zMin,zMax,nZmax,       &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          ddVpIntS(nS)=chebInt(ddVpS(nZmax+1),zMin,zMax,nZmax,      &
               &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          LFIntS(nS)  =chebInt(LFS(nZmax+1),zMin,zMax,nZmax,        &
               &           nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          TayIntS(nS) =chebInt(TayS(nZmax+1),zMin,zMax,nZmax,       &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          TayRIntS(nS)=chebInt(TayRS(nZmax+1),zMin,zMax,nZmax,      &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          TayVIntS(nS)=chebInt(TayVS(nZmax+1),zMin,zMax,nZmax,      &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          RstrIntS(nS)=chebInt(RstrS(nZmax+1),zMin,zMax,nZmax,      &
               &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          AstrIntS(nS)=chebInt(AstrS(nZmax+1),zMin,zMax,nZmax,      &
               &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          StrIntS(nS) =chebInt(StrS(nZmax+1),zMin,zMax,nZmax,       &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          V2IntS(nS)  =chebInt(V2S(nZmax+1),zMin,zMax,nZmax,        &
               &           nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          Bs2IntS(nS) =chebInt(Bs2S(nZmax+1),zMin,zMax,nZmax,       &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          BspIntS(nS) =chebInt(BspS(nZmax+1),zMin,zMax,nZmax,       &
               &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          BspdIntS(nS)=chebInt(BspdS(nZmax+1),zMin,zMax,nZmax,      &
               &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          BpsdIntS(nS)=chebInt(BpsdS(nZmax+1),zMin,zMax,nZmax,      &
               &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
          IF ( V2IntS(nS).LT.0.D0 ) THEN
             VpRIntS(nS)=1.D0
          ELSE
             VpRIntS(nS)=DABS(VpIntS(nS))/DSQRT(V2IntS(nS))
          END IF
          VpRIntS(nS) =DMIN1(1.D0,VpRIntS(nS))
          TayIntS(nS) =LFIntS(nS)/TayIntS(nS)
          TayRIntS(nS)=RstrIntS(nS)/TayRIntS(nS)
          TayVIntS(nS)=StrIntS(nS)/TayVIntS(nS)
       ELSE
          VpIntS(nS)  =VpIntN(nS)
          dVpIntS(nS) =dVpIntN(nS)
          ddVpIntS(nS)=ddVpIntN(nS)
          LFIntS(nS)  =LFIntN(nS)
          RstrIntS(nS)=RstrIntN(nS)
          AstrIntS(nS)=AstrIntN(nS)
          StrIntS(nS) =StrIntN(nS)
          V2IntS(nS)  =V2IntN(nS)
          Bs2IntS(nS) =Bs2IntN(nS)
          BspIntS(nS) =BspIntN(nS)
          BspdIntS(nS)=BspdIntN(nS)
          BpsdIntS(nS)=BpsdIntN(nS)
          VpRIntS(nS) =VpRIntN(nS)
          TayIntS(nS) =TayIntN(nS)                 
          TayRIntS(nS)=TayRIntN(nS)                 
          TayVIntS(nS)=TayVIntN(nS)                 
       END IF

       !------ Boundary Values:
       IF ( lTC ) THEN ! inside TC
          rB(1)=r_CMB
          zMax = DSQRT(r_CMB**2-sZ(nS)**2)
          zMin =-Zmax
          BspB(1) =BspS(1)
          BspB(2) =BspS(nZmaxNS)
          BspdB(1)=BspdS(1)
          BspdB(2)=BspdS(nZmaxNS)
          BpsdB(1)=BpsdS(1)
          BpsdB(2)=BpsdS(nZmaxNS)
          Bs2B(1) =Bs2S(1)
          Bs2B(2) =Bs2S(nZmaxNS)
          CALL getAStr(BszB,BszLMr,2,2,lmMaxS,l_max,           &
               &                   r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
          CALL getAStr(BpzB,BpzLMr,2,2,lmMaxS,l_max,           &
               &                   r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
          CALL getAStr(BzpdB,BzpdLMr,2,2,lmMaxS,l_max,         &
               &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
          CALL getAStr(BpzdB,BpzdLMr,2,2,lmMaxS,l_max,         &
               &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))

          TauBS(nS)  =-(BpzB(2)+sZ(nS)/zMin*BspB(2))
          dTauBS(nS) =-(BpzdB(2)+BzpdB(2) +                         &
               &                      sZ(nS)/zMin*(BspdB(2)+BpsdB(2)))
          dTTauBS(nS)=-(BszB(2)+sZ(nS)/zMin*Bs2B(2))
          TauBN(nS)  =  BpzB(1)+sZ(nS)/zMax*BspB(1)
          dTauBN(nS) =  BpzdB(1)+BzpdB(1) +                         &
               &                      sZ(nS)/zMax*(BspdB(1)+BpsdB(1))
          dTTauBN(nS)=  BszB(1)+sZ(nS)/zMax*Bs2B(1)

          rB(1)=r_ICB 
          zMax = DSQRT(r_ICB**2-sZ(nS)**2)
          zMin =-zMax
          BspB(1) =BspS(nZmax)
          BspB(2) =BspS(nZmax+1)
          BspdB(1)=BspdS(nZmax)
          BspdB(2)=BspdS(nZmax+1)
          BpsdB(1)=BpsdS(nZmax)
          BpsdB(2)=BpsdS(nZmax+1)
          Bs2B(1) =Bs2S(nZmax)
          Bs2B(2) =Bs2S(nZmax+1)
          CALL getAStr(BszB,BszLMr,2,2,lmMaxS,l_max,           &
               &               r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
          CALL getAStr(BpzB,BpzLMr,2,2,lmMaxS,l_max,           &
               &               r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
          CALL getAStr(BzpdB,BzpdLMr,2,2,lmMaxS,l_max,         &
               &                 r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
          CALL getAStr(BpzdB,BpzdLMr,2,2,lmMaxS,l_max,         &
               &                 r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))

          TauBS(nS)  =TauBS(nS)  +(BpzB(2)+sZ(nS)/zMin*BspB(2))
          dTauBS(nS) =dTauBS(nS) +(BpzdB(2)+BzpdB(2) +              &
               &                    sZ(nS)/zMin*(BspdB(2)+BpsdB(2)))
          dTTauBS(nS)=dTTauBS(nS)+(BszB(2)+sZ(nS)/zMin*Bs2B(2))
          TauBN(nS)  =TauBN(nS)  -(BpzB(1) +sZ(nS)/zMax*BspB(1))
          dTauBN(nS) =dTauBN(nS) -(BpzdB(1)+BzpdB(1) +              &
               &                    sZ(nS)/zMax*(BspdB(1)+BpsdB(1)))    
          dTTauBN(nS)=dTTauBN(nS)-(BszB(1)+sZ(nS)/zMax*Bs2B(1))

       ELSE ! outside TC

          rB(1)=r_CMB
          zMax = DSQRT(r_CMB**2-sZ(nS)**2)
          zMin =-Zmax
          BspB(1)=BspS(1)
          BspB(2)=BspS(nZmaxNS)
          BspdB(1)=BspdS(1)
          BspdB(2)=BspdS(nZmaxNS)
          BpsdB(1)=BpsdS(1)
          BpsdB(2)=BpsdS(nZmaxNS)
          Bs2B(1) =Bs2S(1)
          Bs2B(2) =Bs2S(nZmaxNS)
          CALL getAStr(BpzB,BpzLMr,2,2,lmMaxS,l_max,           &
               &                   r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
          CALL getAStr(BzpdB,BzpdLMr,2,2,lmMaxS,l_max,         &
               &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
          CALL getAStr(BpzdB,BpzdLMr,2,2,lmMaxS,l_max,         &
               &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))

          TauBS(nS)  = BpzB(1)+sZ(nS)/zMax*BspB(1) -                &
               &                     BpzB(2)-sZ(nS)/zMin*BspB(2)
          dTauBS(nS) = BpzdB(1)+BzpdB(1) +                          &
               &                     sZ(nS)/zMax*(BspdB(1)+BpsdB(1)) -            &
               &                     BpzdB(2)-BzpdB(2) -                          &
               &                     sZ(nS)/zMin*(BspdB(2)+BpsdB(2))
          dTTauBS(nS)= BszB(1)+sZ(nS)/zMax*Bs2B(1) -                &
               &                     BszB(2)-sZ(nS)/zMin*Bs2B(2)
          TauBN(nS)  =TauBS(nS)
          dTauBN(nS) =dTauBS(nS)
          dTTauBN(nS)=dTTauBS(nS)

       END IF

99     CONTINUE

    END DO  ! Loop over s 
    !$OMP END DO
    ! Integration finished

    CLOSE (95)
    IF ( l_TOZave .AND. nTOsets.GT.1 ) CLOSE (96)

    IF ( lStopRun ) STOP

    !--- Integrate Geostrophic azumithal flow energy and Taylor measure:
    VgRMS  =0.D0
    TayRMS =0.D0
    TayRRMS=0.D0
    TayVRMS=0.D0
    !--- Integrals changed to not represent the different volumes
    !    represented by each cylinder on Nov. 2 2007:
    DO nS=1,nSmax
       ! Old form used again for VgRMS for make is comparable with
       ! VpRMS below.
       VgRMS=VgRMS + 2.D0*pi*h(nS)*sZ(nS)*dsZ *                     &
            &                               VpIntN(nS)*VpIntN(nS)
       TayRMS =TayRMS +dsZ*DABS(TayIntN(nS))
       TayRRMS=TayRRMS+dsZ*DABS(TayRIntN(nS))
       TayVRMS=TayVRMS+dsZ*DABS(TayVIntN(nS))
       IF ( nS.LE.nSI ) THEN
          VgRMS=VgRMS + 2.D0*pi*h(nS)*sZ(nS)*dsZ *                  &
               &                         VpIntS(nS)*VpIntS(nS)
          TayRMS =TayRMS +dsZ*DABS(TayIntS(nS))
          TayRRMS=TayRRMS+dsZ*DABS(TayRIntS(nS))
          TayVRMS=TayVRMS+dsZ*DABS(TayVIntS(nS))
       END IF
    END DO

    !--- s-derivatives:
    !------ Create arrays to be differentiated:
    DO nS=1,nSmax
       Os2(nS)=1.D0/(sZ(nS)*sZ(nS))
       SVpIntN(nS) =VpIntN(nS)/sZ(nS)
       SBspIntN(nS)=h(nS)*sZ(nS)*sZ(nS)*BspIntN(nS)
       SBs2IntN(nS)=h(nS)*sZ(nS)**3*Bs2IntN(nS)
       IF ( sZ(nS).LT.r_ICB ) THEN ! inside TC
          SVpIntS(nS) =VpIntS(nS)/sZ(nS)
          SBspIntS(nS)=h(nS)*sZ(nS)*sZ(nS)*BspIntS(nS)
          SBs2IntS(nS)=h(nS)*sZ(nS)**3*Bs2IntS(nS)
       END IF
    END DO
    !------ Create Derivatives:
    f1=1.D0/(12.D0*dsZ)
    f2=f1/dsZ
    DO nS=3,nSmax-2
       dSVpIntN =f1*(     SVpIntN(nS-2)-8.D0*SVpIntN(nS-1) +        &
            &                   8.D0*SVpIntN(nS+1)-     SVpIntN(nS+2) )
       d2SVpIntN=f2*(    -SVpIntN(nS-2)+16.D0*SVpIntN(nS-1) -       &
            &       30.D0*SVpIntN(nS)+16.D0*SVpIntN(nS+1)-SVpIntN(nS+2))
       dSBspIntN=f1*(     SBspIntN(nS-2)-8.D0*SBspIntN(nS-1) +      &
            &                   8.D0*SBspIntN(nS+1)-     SBspIntN(nS+2) )
       dSBs2IntN=f1*(     SBs2IntN(nS-2)-8.D0*SBs2IntN(nS-1) +      &
            &                   8.D0*SBs2IntN(nS+1)-     SBs2IntN(nS+2) )
       TauN(nS)  =Oh(nS)*(Os2(nS)*dSBspIntN+TauBN(nS))
       dTTauN(nS)=sZ(nS)*d2SVpIntN*Bs2IntN(nS) +                    &
            &               Oh(nS)*(Os2(nS)*dSVpIntN*dSBs2IntN +               &
            &                        sZ(nS)*dSVpIntN*dTTauBN(nS) )
    END DO
    TauN(1)        =TauN(3)
    TauN(2)        =TauN(3)
    TauN(nSmax-1)  =TauN(nSmax-2)
    TauN(nSmax)    =TauN(nSmax-2)
    dTTauN(1)      =dTTauN(3)
    dTTauN(2)      =dTTauN(3)
    dTTauN(nSmax-1)=dTTauN(nSmax-2)
    dTTauN(nSmax)  =dTTauN(nSmax-2)

    !------ South hemisphere:
    DO nS=1,nSmax
       IF ( sZ(nS).LT.r_ICB .AND. nS.GT.2 ) THEN
          dSVpIntS =f1*(SVpIntS(nS-2)-8.D0*SVpIntS(nS-1) +          &
               &                          8.D0*SVpIntS(nS+1)-SVpIntS(nS+2) )
          d2SVpIntS=f2*(-SVpIntS(nS-2)+16.D0*SVpIntS(nS-1) -        &
               &          30.D0*SVpIntS(nS)+16.D0*SVpIntS(nS+1)-SVpIntS(nS+2))
          dSBspIntS=f1*(SBspIntS(nS-2)-8.D0*SBspIntS(nS-1) +        &
               &                      8.D0*SBspIntS(nS+1)-SBspIntS(nS+2) )
          dSBs2IntS=f1*(SBs2IntS(nS-2)-8.D0*SBs2IntS(nS-1) +        &
               &                      8.D0*SBs2IntS(nS+1)-SBs2IntS(nS+2) )
          TauS(nS) =Oh(nS)*(Os2(nS)*dSBspIntS+TauBS(nS))
          dTTauS(nS)=sZ(nS)*d2SVpIntS*Bs2IntS(nS) +                 &
               &               Oh(nS)*(Os2(nS)*dSVpIntS*dSBs2IntS +               &
               &                        sZ(nS)*dSVpIntS*dTTauBS(nS) )
       ELSE
          TauS(nS)  =TauN(nS)
          dTTauS(nS)=dTTauN(nS)
       END IF
    END DO

    !--- Output of z-integral:
    OPEN(nOutFile,FILE=TOfileNhs,STATUS='UNKNOWN',                  &
         &       FORM='UNFORMATTED',POSITION='APPEND')
    IF( nTOsets.EQ.1 ) THEN
       IF ( l_save_out ) THEN
          OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN',           &
               &             POSITION='APPEND')
       END IF
       WRITE(n_log_file,*)
       WRITE(n_log_file,*) '! TO: No. of s-values:',INT(nSmax/r_cmb)
       WRITE(n_log_file,*) '! TO: No. of z-values:',nNorm
       WRITE(n_log_file,*)
       IF ( l_save_out ) CLOSE(n_log_file)
       WRITE(nOutFile) FLOAT(nSmax)                     ! 1
       WRITE(nOutFile) (SNGL(sZ(nS)),nS=1,nSmax)        ! 2
    END IF
    WRITE(nOutFile)                                                 &
         &          SNGL(time),                                             &! 3
         &          (SNGL(VpIntN(nS))             ,nS=1,nSmax),             &! 4
         &          (SNGL(dVpIntN(nS))            ,nS=1,nSmax),             &! 5
         &          (SNGL(ddVpIntN(nS))           ,nS=1,nSmax),             &! 6
         &          (SNGL(VpRIntN(nS))            ,nS=1,nSmax),             &! 7
         &          (SNGL(RstrIntN(nS))           ,nS=1,nSmax),             &! 8
         &          (SNGL(AstrIntN(nS))           ,nS=1,nSmax),             &! 9
         &          (SNGL(LFfac*LFIntN(nS))       ,nS=1,nSmax),             &! 10
         &          (SNGL(StrIntN(nS))            ,nS=1,nSmax),             &! 11
         &          (SNGL(TayIntN(nS))            ,nS=1,nSmax),             &! 12
         &          (SNGL(LFfac*TauN(nS))         ,nS=1,nSmax),             &! 13
         &          (SNGL(LFfac*TauBN(nS)/h(nS))  ,nS=1,nSmax),             &! 14
         &          (SNGL(LFfac*BspdIntN(nS))     ,nS=1,nSmax),             &! 15 For first part of dTau
         &          (SNGL(LFfac*BpsdIntN(nS))     ,nS=1,nSmax),             &! 16 For second part of dTau
         &          (SNGL(LFfac*dTauBN(nS)/h(nS)) ,nS=1,nSmax),             &! 17 Boundary contribution !
         &          (SNGL(LFfac*dTTauN(nS))       ,nS=1,nSmax),             &! 18
         &          (SNGL(LFfac*dTTauBN(nS)/h(nS)),nS=1,nSmax),             &! 19
         &          (SNGL(LFfac*Bs2IntN(nS))      ,nS=1,nSmax)  ! 20 New June 7, 2006
    CLOSE(nOutFile)

    OPEN(nOutFile,FILE=TOfileShs,STATUS='UNKNOWN',                  &
         &       FORM='UNFORMATTED',POSITION='APPEND')
    IF( nTOsets.EQ.1 ) THEN
       WRITE(nOutFile) FLOAT(nSmax)
       WRITE(nOutFile) (SNGL(sZ(nS)),nS=1,nSmax)
    END IF
    WRITE(nOutFile)                                                 &
         &          SNGL(time),                                             &
         &          (SNGL(VpIntS(nS))             ,nS=1,nSmax),             &
         &          (SNGL(dVpIntS(nS))            ,nS=1,nSmax),             &
         &          (SNGL(ddVpIntS(nS))           ,nS=1,nSmax),             &
         &          (SNGL(VpRIntS(nS))            ,nS=1,nSmax),             &
         &          (SNGL(RstrIntS(nS))           ,nS=1,nSmax),             &
         &          (SNGL(AstrIntS(nS))           ,nS=1,nSmax),             &
         &          (SNGL(LFfac*LFIntS(nS))       ,nS=1,nSmax),             &
         &          (SNGL(StrIntS(nS))            ,nS=1,nSmax),             &
         &          (SNGL(TayIntS(nS))            ,nS=1,nSmax),             &
         &          (SNGL(LFfac*TauS(nS))         ,nS=1,nSmax),             &
         &          (SNGL(LFfac*TauBS(nS)/h(nS))  ,nS=1,nSmax),             &
         &          (SNGL(LFfac*BspdIntS(nS))     ,nS=1,nSmax),             &! For first part of dTau
         &          (SNGL(LFfac*BpsdIntS(nS))     ,nS=1,nSmax),             &! For second part of dTau
         &          (SNGL(LFfac*dTauBS(nS)/h(nS)) ,nS=1,nSmax),             &! Boundary contribution !
         &          (SNGL(LFfac*dTTauS(nS))       ,nS=1,nSmax),             &
         &          (SNGL(LFfac*dTTauBS(nS)/h(nS)),nS=1,nSmax),             &
         &          (SNGL(LFfac*Bs2IntS(nS))      ,nS=1,nSmax)
    CLOSE(nOutFile)

    !       WRITE(98,*) TauBS(30),TauBN(30),dTauBS(30),dTauBN(30),
    !    &              dTTauBS(30),dTTauBN(30)
    !       WRITE(98,*) 'IntS:',VpIntS(30),RstrIntS(30),TayIntS(30),
    !    &                          TauS(30),Bs2IntS(30)
    !       WRITE(98,*) 'IntN:',VpIntN(30),RstrIntN(30),TayIntN(30),
    !    &                          TauN(30),Bs2IntN(30)


    ! Note: integral changed on Nov. 2 2007 (see above)
    fac=1.D0/sZ(nSmax) ! I integrate (roughly) from s=0 to s=r_cmb.
    VgRMS  =vSF*DSQRT(VgRMS/vol_oc)
    TayRMS =fac*TayRMS
    TayRRMS=fac*TayRRMS
    TayVRMS=fac*TayVRMS

    lTOrms=.TRUE.
    IF ( lTOmov .OR. lTOrms ) THEN

       !--- Output of graphic file:
       !--- Transform back to radial space:
       DO nR=1,n_r_max
          DO l=1,l_max
             lm=lm2(l,0)
             !              dzRstrLMr(lm,nR)=dzRstrLMr(lm,nR) +
             !     +                         dzAstrLMr(lm,nR)
             dzVpLMr(l+1,nR)=REAL(z(lm,nR)) ! instead of transform copy again
          END DO
       END DO

       CALL costf1(dzdVpLMr,lmMaxS,1,lmMaxS,                        &
            &                 workA,i_costf_init,d_costf_init)
       CALL costf1(dzRstrLMr,lmMaxS,1,lmMaxS,                       &
            &                 workA,i_costf_init,d_costf_init)
       CALL costf1(dzAstrLMr,lmMaxS,1,lmMaxS,                       &
            &                 workA,i_costf_init,d_costf_init)
       CALL costf1(dzStrLMr,lmMaxS,1,lmMaxS,                        &
            &                 workA,i_costf_init,d_costf_init)
       CALL costf1(dzLFLMr,lmMaxS,1,lmMaxS,                         &
            &                 workA,i_costf_init,d_costf_init)
       CALL costf1(dzCorLMr,lmMaxS,1,lmMaxS,                        &
            &                 workA,i_costf_init,d_costf_init)

       !--- Open output file
       nFields=7
       IF ( lTOmov ) THEN

          OPEN(nOutFile,file=movFile,STATUS='UNKNOWN',              &
               &             FORM='UNFORMATTED',POSITION='APPEND')

          !--- Write header into output file:
          IF ( nTOmovSets.EQ.0 ) THEN
             version='JW_Movie_Version_2'
             WRITE(nOutFile) version
             dumm(1)=102           ! type of input
             dumm(2)=3             ! marker for constant phi plane
             dumm(3)=0.D0          ! surface constant
             dumm(4)=nFields       ! no of fields 
             WRITE(nOutFile) (real(dumm(n),4),n=1,4)

             dumm(1)=11.0          ! Field marker for AS vPhi 
             dumm(2)=61.0          ! Field marker for Reynolds Force 
             dumm(3)=62.0          ! Field marker for Advective Force  
             dumm(4)=63.0          ! Field marker for Viscous Force 
             dumm(5)=64.0          ! Field marker for Lorentz Force
             dumm(6)=65.0          ! Field marker for Coriolis force
             dumm(7)=66.0          ! Field marker for dtVp
             WRITE(nOutFile) (SNGL(dumm(n)),n=1,nFields)

             !------ Now other info about grid and parameters:
             WRITE(nOutFile) runid        ! run identifier 
             dumm( 1)=n_r_max          ! total number of radial points
             dumm( 2)=n_r_max          ! no of radial point in outer core
             dumm( 3)=n_theta_max      ! no. of theta points
             dumm( 4)=n_phi_max        ! no. of phi points
             dumm( 5)=minc             ! imposed symmetry
             dumm( 6)=ra               ! control parameters
             dumm( 7)=ek               ! (for information only)
             dumm( 8)=pr               !      -"-
             dumm( 9)=prmag            !      -"-
             dumm(10)=radratio         ! ratio of inner / outer core
             dumm(11)=tScale           ! timescale
             WRITE(nOutFile) (SNGL(dumm(n)),     n=1,11)  
             WRITE(nOutFile) (SNGL(r(n)/r_CMB),  n=1,n_r_max)
             WRITE(nOutFile) (SNGL(theta_ord(n)),n=1,n_theta_max)
             WRITE(nOutFile) (SNGL(phi(n)),      n=1,n_phi_max)

          END IF ! Write Header ?

          nTOmovSets=nTOmovSets+1

          dumm(1)=nTOmovSets        ! time frame number for movie
          dumm(2)=time              ! time      
          dumm(3)=0.D0      
          dumm(4)=0.D0     
          dumm(5)=0.D0          
          dumm(6)=0.D0          
          dumm(7)=0.D0         
          dumm(8)=0.D0        
          WRITE(nOutFile) (SNGL(dumm(n)),n=1,8)

       END IF ! Produce movie output ?

       IF ( lTOrms ) THEN

          DO nR=1,n_r_max
             VpR(nR)   =0.D0
             dVpR(nR)  =0.D0
             RstrR(nR) =0.D0
             AstrR(nR) =0.D0
             LFR(nR)   =0.D0
             LFABSR(nR)=0.D0
             StrR(nR)  =0.D0
             CorR(nR)  =0.D0
          END DO

          nTOrmsSets=nTOrmsSets+1
          OPEN(nOutFile2,FILE=tayFile,FORM='UNFORMATTED',           &
               &             STATUS='UNKNOWN',POSITION='APPEND')
          IF ( nTOrmsSets.EQ.1 ) THEN
             WRITE(nOutFile2) FLOAT(n_r_max)
             WRITE(nOutFile2) (SNGL(r(nR)),nR=1,n_r_max)
          END IF

       END IF

       DO nOut=1,nFields ! Loop over four output fields

          DO nR=1,n_r_max ! Loop over radial points
             rS=r(nR)

             DO n=1,nThetaBs ! Loop over theta blocks
                nThetaStart=(n-1)*sizeThetaB+1

                !------ Convert from lm to theta block: 
                IF ( nOut.EQ.1 ) THEN
                   CALL get_PAS(dzVpLMr(1,nR),outBlock,             &
                        &                              rS,nThetaStart,sizeThetaB)
                ELSE IF ( nOut.EQ.2 ) THEN
                   CALL get_PAS(dzRstrLMr(1,nR),outBlock,           &
                        &                              rS,nThetaStart,sizeThetaB)
                ELSE IF ( nOut.EQ.3 ) THEN
                   CALL get_PAS(dzAstrLMr(1,nR),outBlock,           &
                        &                              rS,nThetaStart,sizeThetaB)
                ELSE IF ( nOut.EQ.4 ) THEN
                   CALL get_PAS(dzStrLMr(1,nR),outBlock,            &
                        &                              rS,nThetaStart,sizeThetaB)
                ELSE IF ( nOut.EQ.5 ) THEN
                   CALL get_PAS(dzLFLMr(1,nR),outBlock,             &
                        &                              rS,nThetaStart,sizeThetaB)
                ELSE IF ( nOut.EQ.6 ) THEN
                   CALL get_PAS(dzCorLMr(1,nR),outBlock,            &
                        &                              rS,nThetaStart,sizeThetaB)
                ELSE IF ( nOut.EQ.7 ) THEN
                   CALL get_PAS(dzdVpLMr(1,nR),outBlock,            &
                        &                              rS,nThetaStart,sizeThetaB)
                END IF

                !------ Storage of field in fout for theta block,
                !       integration and Max/Min values
                nTheta=(n-1)*sizeThetaB
                DO nThetaBlock=1,sizeThetaB
                   nTheta=nTheta+1
                   nThetaNHS=(nTheta+1)/2
                   sS       =rS*sinTheta(nTheta)
                   !--------- Convert to correct order in theta grid points:
                   IF ( MOD(nTheta,2).EQ.1 ) THEN
                      nThetaOrd=(nTheta+1)/2
                   ELSE
                      nThetaOrd=n_theta_max-nTheta/2+1
                   END IF
                   nPos=(nR-1)*n_theta_max+nThetaOrd

                   IF ( nOut.EQ.1 ) THEN
                      !--------------- Zonal flow:
                      fOut(nPos)=vSF*outBlock(nThetaBlock)
                      VpR(nR)=VpR(nR) +                             &
                           &                       gauss(nThetaNHS)*fOut(nPos)/sS
                   ELSE IF ( nOut.EQ.2 ) THEN
                      !--------------- Reynolds force:
                      fOut(nPos)=fSF*outBlock(nThetaBlock)
                      RstrR(nR)=RstrR(nR) +                         &
                           &                               gauss(nThetaNHS)*fOut(nPos)/sS
                   ELSE IF ( nOut.EQ.3 ) THEN
                      !--------------- Advective force:
                      fOut(nPos)=fSF*outBlock(nThetaBlock)
                      AstrR(nR)=AstrR(nR) +                         &
                           &                                gauss(nThetaNHS)*fOut(nPos)/sS
                   ELSE IF ( nOut.EQ.4 ) THEN
                      !--------------- Viscous force:
                      fOut(nPos)=fSF*outBlock(nThetaBlock)
                      StrR(nR)=StrR(nR) +                           &
                           &                                gauss(nThetaNHS)*fOut(nPos)/sS           
                   ELSE IF ( nOut.EQ.5 ) THEN
                      !--------------- Lorentz force:
                      fOut(nPos)=fSF*LFfac*outBlock(nThetaBlock)
                      LFR(nR)=LFR(nR) +                             &
                           &                            gauss(nThetaNHS)*fOut(nPos)/sS
                      LFABSR(nR)=LFABSR(nR) +                       &
                           &                            gauss(nThetaNHS)*DABS(fOut(nPos))/sS
                   ELSE IF ( nOut.EQ.6 ) THEN
                      !--------------- Corriolis force:
                      fOut(nPos)=fSF*outBlock(nThetaBlock)
                      CorR(nR)=CorR(nR) +                           &
                           &                                gauss(nThetaNHS)*fOut(nPos)/sS
                   ELSE IF ( nOut.EQ.7 ) THEN
                      !--------------- dtVp:
                      fOut(nPos)=vSF*outBlock(nThetaBlock)
                      dVpR(nR)  =dVpR(nR) +                         &
                           &                                gauss(nThetaNHS)*fOut(nPos)/sS
                   END IF

                END DO ! Loop over thetas in block
             END DO    ! Loop over theta blocks

          END DO ! Loop over R

          !------ Output of stress contributions for one radial grid point:
          !       Write all fields into movie style file
          IF ( lTOmov )                                             &
               &           WRITE(nOutFile) (SNGL(fOut(nPos)),nPos=1,nFieldSize)

       END DO ! Loop over output functions

       IF ( l_save_out ) THEN
          OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN',           &
               &             POSITION='APPEND')
       END IF
       IF ( lTOmov ) THEN 
          CLOSE(nOutFile)
          WRITE(n_log_file,'(1p,/,                                  &
               &           '' ! WRITING TO MOVIE FRAME NO '',i8,                  &
               &           '' AT STEP/TIME :'',i8,d16.6)')                        &
               &              nTOmovSets,n_time_step,time*tScale
       END IF

       !--- Realtive importance of azimuthal and geostrophic flow
       VRMS   =DSQRT(2.D0*eKin/vol_oc)
       VpRMS  =DSQRT(2.D0*eKinTAS/vol_oc)

       IF ( VRMS /= 0d0) THEN
          VpRMS  =VpRMS/VRMS
          VgRMS  =VgRMS/VRMS
       END IF

       DO nR=1,n_r_max
          fac=1.D0/2.D0 ! Cancel factor 2 from theta integral
          VpR(nR)    =fac*VpR(nR)
          dVpR(nR)   =fac*dVpR(nR)
          RstrR(nR)  =fac*RstrR(nR)
          AstrR(nR)  =fac*AstrR(nR)
          StrR(nR)   =fac*StrR(nR)
          LFR(nR)    =fac*LFR(nR)
          CorR(nR)   =fac*CorR(nR)
          IF ( LFABSR(nR) /= 0d0 ) THEN
             TayR(nR)   =DABS(LFR(nR))/(fac*LFABSR(nR))
          ELSE
             TayR(nR)   =0d0
          END IF
          !              TayRMSR(nR)=rS*rS*TayR(nR)
          TayRMSR(nR)=TayR(nR)
       END DO

       !--- Now perform the radial integral: ( not tested )
       TaySRMS=rInt_R(TayRMSR,n_r_max,n_r_max,drx,                  &
            &                      i_costf_init,d_costf_init)
       !--- And finally calculate the mean value, the factor 4*pi comes from
       !    the fact that the surface integral has already been cared for
       !    NOTE: Integral for RMS Taylorisation changed to not respect the 
       !    different volumes represented by each shell on 2 Nov 2007.
       !          fac    =4.D0*pi/vol_oc
       fac    =1.D0/(r_cmb-r_icb) ! =1
       TaySRMS=fac*TaySRMS

       WRITE(nOutFile2)                                             &
            &          SNGL(time),SNGL(VpRMS**2),SNGL(VgRMS**2),               &
            &          SNGL(TayRMS),SNGL(TaySRMS),SNGL(TayRRMS),               &
            &          SNGL(TayVRMS),SNGL(eKin),                               &! 3
            &          (SNGL(VpR(nR))  ,nR=1,n_r_max),                         &! 4
            &          (SNGL(dVpR(nR)) ,nR=1,n_r_max),                         &! 5
            &          (SNGL(RstrR(nR)),nR=1,n_r_max),                         &! 6
            &          (SNGL(AstrR(nR)),nR=1,n_r_max),                         &! 7
            &          (SNGL(LFR(nR))  ,nR=1,n_r_max),                         &! 8
            &          (SNGL(StrR(nR)) ,nR=1,n_r_max),                         &! 9
            &          (SNGL(CorR(nR)) ,nR=1,n_r_max),                         &! 10
            &          (SNGL(TayR(nR)) ,nR=1,n_r_max)  ! 11
       CLOSE(nOutFile2)

    END IF

    timeLast=time
    IF ( lVerbose ) WRITE(*,*) '! End of outTO!'

    lTOZwrite=.FALSE.

    RETURN
  END SUBROUTINE outTO

  !---------------------------------------------------------------------------------
END MODULE outTO_mod
