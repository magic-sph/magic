!$Id$

MODULE radialLoop
    use omp_lib
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE torsional_oscillations
    USE blocking
    USE horizontal_data
    USE logic
    USE dtB_mod
    USE output_data
    USE const
    USE parallel_mod
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
    IMPLICIT NONE

    !---- Nonlinear terms, field components and help arrays
    !     for Legendre transform, field components, and dtB and 
    !     TO output. These are all stored in COMMON BLOCKS
    !     that are thread privat rather than using explicit 
    !     PRIVATE statements in the PARALLEL DO clause.

    !----- Nonlinear terms in lm-space: 
    COMPLEX(kind=8),ALLOCATABLE :: AdvrLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: AdvtLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: AdvpLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: LFrLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: LFtLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: LFpLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: VxBrLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: VxBtLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: VxBpLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: VSrLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: VStLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: VSpLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: ViscHeatLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: OhmLossLM(:)

    !----- Nonlinear terms in phi/theta space: 
    REAL(kind=8),ALLOCATABLE :: Advr(:,:)
    REAL(kind=8),ALLOCATABLE :: Advt(:,:)
    REAL(kind=8),ALLOCATABLE :: Advp(:,:)
    REAL(kind=8),ALLOCATABLE :: LFr(:,:)
    REAL(kind=8),ALLOCATABLE :: LFt(:,:)
    REAL(kind=8),ALLOCATABLE :: LFp(:,:)
    REAL(kind=8),ALLOCATABLE :: VxBr(:,:)
    REAL(kind=8),ALLOCATABLE :: VxBt(:,:)
    REAL(kind=8),ALLOCATABLE :: VxBp(:,:)
    REAL(kind=8),ALLOCATABLE :: VSr(:,:)
    REAL(kind=8),ALLOCATABLE :: VSt(:,:)
    REAL(kind=8),ALLOCATABLE :: VSp(:,:)
    REAL(kind=8),ALLOCATABLE :: ViscHeat(:,:)
    REAL(kind=8),ALLOCATABLE :: OhmLoss(:,:)

    !        COMMON/NLterms/AdvrLM,AdvtLM,AdvpLM,LFrLM,LFtLM,LFpLM,
    !     &                 VxBrLM,VxBtLM,VxBpLM,VSrLM,VStLM,VSpLM,
    !     &                 ViscHeatLM,OhmLossLM,
    !     &                 Advr,Advt,Advp,LFr,LFt,LFp,
    !     &                 VxBr,VxBt,VxBp,VSr,VSt,VSp,
    !     &                 ViscHeat,OhmLoss
    !!$OMP THREADPRIVATE (/NLterms/)
    !$OMP THREADPRIVATE( AdvrLM,AdvtLM,AdvpLM,LFrLM,LFtLM,LFpLM )
    !$OMP THREADPRIVATE( VxBrLM,VxBtLM,VxBpLM,VSrLM,VStLM,VSpLM )
    !$OMP THREADPRIVATE( ViscHeatLM,OhmLossLM )
    !$OMP THREADPRIVATE( Advr,Advt,Advp,LFr,LFt,LFp )
    !$OMP THREADPRIVATE( VxBr,VxBt,VxBp,VSr,VSt,VSp )
    !$OMP THREADPRIVATE( ViscHeat,OhmLoss )

    !----- Fields calculated from these help arrays by legtf:
    COMPLEX(kind=8),ALLOCATABLE :: vrc(:,:),vtc(:,:),vpc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: dvrdrc(:,:),dvtdrc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: dvpdrc(:,:),cvrc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: dvrdtc(:,:),dvrdpc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: dvtdpc(:,:),dvpdpc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: brc(:,:),btc(:,:),bpc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: cbrc(:,:),cbtc(:,:),cbpc(:,:)
    COMPLEX(kind=8),ALLOCATABLE :: sc(:,:),drSc(:,:)

    !        COMMON/components/vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,
    !     &                    cvrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc,
    !     &                  brc,btc,bpc,cbrc,cbtc,cbpc,sc,drSc
    !!$OMP THREADPRIVATE (/components/)
    !$OMP THREADPRIVATE( vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc )
    !$OMP THREADPRIVATE( cvrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc )
    !$OMP THREADPRIVATE( brc,btc,bpc,cbrc,cbtc,cbpc,sc,drSc )

    !----- Help arrays for Legendre transform calculated in legPrepG:
    !      Parallelizatio note: these are the R-distributed versions
    !      of the field scalars.
    COMPLEX(kind=8),ALLOCATABLE :: dLhw(:)
    COMPLEX(kind=8),ALLOCATABLE :: dLhdw(:)
    COMPLEX(kind=8),ALLOCATABLE :: dLhz(:)
    COMPLEX(kind=8),ALLOCATABLE :: dLhb(:)
    COMPLEX(kind=8),ALLOCATABLE :: dLhj(:)
    COMPLEX(kind=8),ALLOCATABLE :: vhG(:)
    COMPLEX(kind=8),ALLOCATABLE :: vhC(:)
    COMPLEX(kind=8),ALLOCATABLE :: dvhdrG(:)
    COMPLEX(kind=8),ALLOCATABLE :: dvhdrC(:)
    COMPLEX(kind=8),ALLOCATABLE :: bhG(:)
    COMPLEX(kind=8),ALLOCATABLE :: bhC(:)
    COMPLEX(kind=8),ALLOCATABLE :: cbhG(:)
    COMPLEX(kind=8),ALLOCATABLE :: cbhC(:)
    !----- R-distributed versions of scalar fields (see c_fields.f):
    COMPLEX(kind=8),ALLOCATABLE :: sR(:),dsR(:)
    COMPLEX(kind=8),ALLOCATABLE :: preR(:),dpR(:)
    REAL(kind=8),ALLOCATABLE :: zAS(:),dzAS(:),ddzAS(:) ! used in TO
    REAL(kind=8) :: omegaIC,omegaMA
    COMPLEX(kind=8),ALLOCATABLE :: bCMB(:)

    !        COMMON/LegHelp/dLhw,dLhdw,dLhz,dLhb,dLhj,vhG,vhC,
    !     &                   dvhdrG,dvhdrC,bhG,bhC,cbhG,cbhC,
    !     &         sR,dsR,preR,dpR,zAS,dzAS,ddzAS,bCMB,omegaIC,omegaMA
!!$OMP THREADPRIVATE(/LegHelp/)
    !$OMP THREADPRIVATE( dLhw,dLhdw,dLhz,dLhb,dLhj,vhG,vhC )
    !$OMP THREADPRIVATE( dvhdrG,dvhdrC,bhG,bhC,cbhG,cbhC )
    !$OMP THREADPRIVATE( sR,dsR,preR,dpR,zAS,dzAS,ddzAS,bCMB,omegaIC,omegaMA )

    !----- Local dtB output stuff:
    COMPLEX(kind=8),ALLOCATABLE :: BtVrLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BpVrLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BrVtLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BrVpLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BtVpLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BpVtLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BtVpCotLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BpVtCotLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BtVpSn2LM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BpVtSn2LM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BrVZLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BtVZLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BtVZcotLM(:)
    COMPLEX(kind=8),ALLOCATABLE :: BtVZsn2LM(:)

    !        COMMON/dtBradialLoop/BtVrLM,BpVrLM,BrVtLM,
    !     &             BrVpLM,BtVpLM,BpVtLM,BtVpCotLM,
    !     &              BpVtCotLM,BtVpSn2LM,BpVtSn2LM,
    !     &          BrVZLM,BtVZLM,BtVZcotLM,BtVZsn2LM
    !!$OMP THREADPRIVATE(/dtBradialLoop/)
    !$OMP THREADPRIVATE( BtVrLM,BpVrLM,BrVtLM )
    !$OMP THREADPRIVATE( BrVpLM,BtVpLM,BpVtLM,BtVpCotLM )
    !$OMP THREADPRIVATE( BpVtCotLM,BtVpSn2LM,BpVtSn2LM )
    !$OMP THREADPRIVATE( BrVZLM,BtVZLM,BtVZcotLM,BtVZsn2LM )

    !----- Local TO output stuff:
    REAL(kind=8),ALLOCATABLE :: dzRstrLM(:),dzAstrLM(:)
    REAL(kind=8),ALLOCATABLE :: dzCorLM(:),dzLFLM(:)

    !        COMMON/TOradialLoop/dzRstrLM,dzAstrLM,dzCorLM,dzLFLM
    !$OMP THREADPRIVATE( dzRstrLM,dzAstrLM,dzCorLM,dzLFLM )


    !----- Saved magnetic field components from last time step:
    !      This is needed for the current TO version. However,
    !      the variables calulated with this don't give any
    !      deep insight. TO should be changes in the future to
    !      eliminate this.
    REAL(kind=8),ALLOCATABLE :: BsLast(:,:,:)
    REAL(kind=8),ALLOCATABLE :: BpLast(:,:,:)
    REAL(kind=8),ALLOCATABLE :: BzLast(:,:,:)

CONTAINS
  SUBROUTINE initialize_radialLoop

    !$OMP PARALLEL
    !----- Nonlinear terms in lm-space: 
    ALLOCATE( AdvrLM(lmP_max) )
    ALLOCATE( AdvtLM(lmP_max) )
    ALLOCATE( AdvpLM(lmP_max) )
    ALLOCATE( LFrLM(lmP_max) )
    ALLOCATE( LFtLM(lmP_max) )
    ALLOCATE( LFpLM(lmP_max) )
    ALLOCATE( VxBrLM(lmP_max) )
    ALLOCATE( VxBtLM(lmP_max) )
    ALLOCATE( VxBpLM(lmP_max) )
    ALLOCATE( VSrLM(lmP_max) )
    ALLOCATE( VStLM(lmP_max) )
    ALLOCATE( VSpLM(lmP_max) )
    ALLOCATE( ViscHeatLM(lmP_max) )
    ALLOCATE( OhmLossLM(lmP_max) )

    !----- Nonlinear terms in phi/theta space: 
    ALLOCATE( Advr(nrp,nfs) )
    ALLOCATE( Advt(nrp,nfs) )
    ALLOCATE( Advp(nrp,nfs) )
    ALLOCATE( LFr(nrp,nfs) )
    ALLOCATE( LFt(nrp,nfs) )
    ALLOCATE( LFp(nrp,nfs) )
    ALLOCATE( VxBr(nrp,nfs) )
    ALLOCATE( VxBt(nrp,nfs) )
    ALLOCATE( VxBp(nrp,nfs) )
    ALLOCATE( VSr(nrp,nfs) )
    ALLOCATE( VSt(nrp,nfs) )
    ALLOCATE( VSp(nrp,nfs) )
    ALLOCATE( ViscHeat(nrp,nfs) )
    ALLOCATE( OhmLoss(nrp,nfs) )

    !----- Fields calculated from these help arrays by legtf:
    ALLOCATE( vrc(ncp,nfs),vtc(ncp,nfs),vpc(ncp,nfs) )
    ALLOCATE( dvrdrc(ncp,nfs),dvtdrc(ncp,nfs) )
    ALLOCATE( dvpdrc(ncp,nfs),cvrc(ncp,nfs) )
    ALLOCATE( dvrdtc(ncp,nfs),dvrdpc(ncp,nfs) )
    ALLOCATE( dvtdpc(ncp,nfs),dvpdpc(ncp,nfs) )
    ALLOCATE( brc(ncp,nfs),btc(ncp,nfs),bpc(ncp,nfs) )
    ALLOCATE( cbrc(ncp,nfs),cbtc(ncp,nfs),cbpc(ncp,nfs) )
    ALLOCATE( sc(ncp,nfs),drSc(ncp,nfs) )

    !----- Help arrays for Legendre transform calculated in legPrepG:
    !      Parallelizatio note: these are the R-distributed versions
    !      of the field scalars.
    ALLOCATE( dLhw(lm_max) )
    ALLOCATE( dLhdw(lm_max) )
    ALLOCATE( dLhz(lm_max) )
    ALLOCATE( dLhb(lm_max) )
    ALLOCATE( dLhj(lm_max) )
    ALLOCATE( vhG(lm_max) )
    ALLOCATE( vhC(lm_max) )
    ALLOCATE( dvhdrG(lm_max) )
    ALLOCATE( dvhdrC(lm_max) )
    ALLOCATE( bhG(lm_max) )
    ALLOCATE( bhC(lm_max) )
    ALLOCATE( cbhG(lm_max) )
    ALLOCATE( cbhC(lm_max) )
    !----- R-distributed versions of scalar fields (see c_fields.f):
    ALLOCATE( sR(lm_max),dsR(lm_max) )
    ALLOCATE( preR(lm_max),dpR(lm_max) )
    ALLOCATE( zAS(l_max+1),dzAS(l_max+1),ddzAS(l_max+1) ) ! used in TO

    ALLOCATE( bCMB(lm_max) )

    !----- Local dtB output stuff:
    ALLOCATE( BtVrLM(lmP_max_dtB) )
    ALLOCATE( BpVrLM(lmP_max_dtB) )
    ALLOCATE( BrVtLM(lmP_max_dtB) )
    ALLOCATE( BrVpLM(lmP_max_dtB) )
    ALLOCATE( BtVpLM(lmP_max_dtB) )
    ALLOCATE( BpVtLM(lmP_max_dtB) )
    ALLOCATE( BtVpCotLM(lmP_max_dtB) )
    ALLOCATE( BpVtCotLM(lmP_max_dtB) )
    ALLOCATE( BtVpSn2LM(lmP_max_dtB) )
    ALLOCATE( BpVtSn2LM(lmP_max_dtB) )
    ALLOCATE( BrVZLM(lmP_max_dtB) )
    ALLOCATE( BtVZLM(lmP_max_dtB) )
    ALLOCATE( BtVZcotLM(lmP_max_dtB) )
    ALLOCATE( BtVZsn2LM(lmP_max_dtB) )

    !----- Local TO output stuff:
    ALLOCATE( dzRstrLM(l_max+2),dzAstrLM(l_max+2) )
    ALLOCATE( dzCorLM(l_max+2),dzLFLM(l_max+2) )
    !$OMP END PARALLEL

    allocate( BsLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr) )
    allocate( BpLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr) )
    allocate( BzLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr) )
    
  END SUBROUTINE initialize_radialLoop

  !***********************************************************************
  SUBROUTINE radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,           &
       &                    lTOCalc,lTONext,lTONext2,lHelCalc,lRmsCalc,   &
       &                  dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM,   &
       &                           lorentz_torque_ic,lorentz_torque_ma,   &
       &                                     br_vt_lm_cmb,br_vp_lm_cmb,   &
       &                                     br_vt_lm_icb,br_vp_lm_icb,   &
       &                                       TstrRLM,TadvRLM,TomeRLM,   &
       &    HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,dtrkc,dthkc)
    !***********************************************************************

    !    !------------ This is release 2 level 10  --------------!
    !    !------------ Created on 2/5/02  by JW. -----------

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine performs the actual time-stepping.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+
    !  |  ruler                                                            |
    !  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
    !--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


    !--- Input of variables:
    LOGICAL :: l_graph,l_cour,l_frame
    LOGICAL :: lTOcalc,lTONext,lTONext2,lHelCalc
    LOGICAL :: lRmsCalc
    REAL(kind=8) :: time,dt,dtLast

    !--- Output:

    !---- Output of explicit time step:
    COMPLEX(kind=8) :: dwdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dzdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dpdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dsdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dbdt(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: djdt(lm_maxMag,n_r_maxMag)
    REAL(kind=8) :: lorentz_torque_ma,lorentz_torque_ic
    !---- Output of contributions to explicit time step that
    !     need a further treatment (radial derivatives required):
    COMPLEX(kind=8) :: dVxBhLM(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: dVSrLM(lm_max,n_r_max) 

    !---- Output of arrays for torsionla oscillation analysis:
    ! include 'c_TO.f'

    !---- Output for axisymmetric helicity:
    REAL(kind=8) :: HelLMr(l_max+1,n_r_max)
    REAL(kind=8) :: Hel2LMr(l_max+1,n_r_max)
    REAL(kind=8) :: HelnaLMr(l_max+1,n_r_max)
    REAL(kind=8) :: Helna2LMr(l_max+1,n_r_max)
    REAL(kind=8) :: uhLMr(l_max+1,n_r_max)
    REAL(kind=8) :: duhLMr(l_max+1,n_r_max)

    !---- Output of nonlinear products for nonlinear
    !     magnetic boundary conditions (needed in s_updateB.f):
    COMPLEX(kind=8) :: br_vt_lm_cmb(lmP_max)    ! product br*vt at CMB
    COMPLEX(kind=8) :: br_vp_lm_cmb(lmP_max)    ! product br*vp at CMB
    COMPLEX(kind=8) :: br_vt_lm_icb(lmP_max)    ! product br*vt at ICB
    COMPLEX(kind=8) :: br_vp_lm_icb(lmP_max)    ! product br*vp at ICB

    !---- Output for Courant criteria:
    REAL(kind=8) :: dtrkc(n_r_max),dthkc(n_r_max)

    !---- dtB output to calculate dynamo action:
    ! include 'c_dtB.f'
    COMPLEX(kind=8) :: TstrRLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: TadvRLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: TomeRLM(lm_max_dtB,n_r_max_dtB)


    !--- Local variables:


    !---- Counter, logicals ...
    INTEGER :: nR,nRC,nRC2,nRC2stop
    INTEGER :: nR_Mag
    INTEGER :: nRstart,nRstop
    INTEGER :: nBc,nPhi,nTh,lm,l
    INTEGER :: nTheta,nThetaB,nThetaLast
    INTEGER :: nThetaStart,nThetaStop
    !INTEGER OMP_GET_THREAD_NUM
    LOGICAL :: lDeriv
    LOGICAL :: lOutBc
    LOGICAL :: lMagNlBc
    LOGICAL :: lGraphHeader    ! Write header into graph file


    !---- End of declaration
    !----------------------------------------------------------------------


    lGraphHeader=l_graph
    IF ( lGraphHeader )                                             &
         &     CALL graphOut(time,nR,ngform,vrc,vtc,vpc,                    &
         &                               brc,btc,bpc,sc,                    &
         &          nThetaStart,sizeThetaB,lGraphHeader)

    IF ( l_cour ) THEN
       dtrkc(n_r_icb)=1.D10
       dthkc(n_r_icb)=1.D10
       dtrkc(n_r_cmb)=1.D10
       dthkc(n_r_cmb)=1.D10
    END IF

    !------ Set nonlinear terms that are possibly needed at the boundaries.
    !       They may be overwritten by get_td later.
    DO lm=1,lm_max
       dVSrLM(lm,n_r_cmb) =zero
       dVSrLM(lm,n_r_icb) =zero
       IF ( l_mag ) THEN
          dVxBhLM(lm,n_r_cmb)=zero
          dVxBhLM(lm,n_r_icb)=zero
       END IF
    END DO

    !------ Having to calculate non-linear boundary terms?
    lMagNlBc=.FALSE.
    IF ( ( l_mag_nl .OR. l_mag_kin ) .AND.                          &
         &       ( ktopv.EQ.1 .OR. l_cond_ma .OR.                           &
         &          ( ktopv.EQ.2 .AND. l_rot_ma ) ) .OR.                    &
         &       ( kbotv.EQ.1 .OR. l_cond_ic .OR.                           &
         &          ( kbotv.EQ.2 .AND. l_rot_ic ) ) )                       &
         &     lMagNlBc=.TRUE.

    !------ When boundary output, Courant criterion, or non-magnetic 
    !       boundary conditions are required I have to calculate 
    !       the fields at the boundaries. This is done in one thread and 
    !       is triggered by lOutBc=.TRUE.
    lOutBc=.FALSE.
    IF ( lTOCalc .OR. lHelCalc .OR. l_frame .OR.                    &
         &       l_cour .OR. l_dtB .OR. lMagNlBc .OR. l_graph               &
         &     ) lOutBc=.TRUE.
    nRstart=n_r_cmb
    nRstop =n_r_icb-1

    !--- Start the big do loop over the radial threads:

    !$OMP PARALLEL DO SCHEDULE(STATIC,sizeRB)                               &
    !$OMP  PRIVATE(nTh,nRC2stop,nRC2,nR,nR_mag,nBc,lDeriv,l,nThetaB,        &
    !$OMP          nThetaLast,nThetaStart,nThetaStop,nTheta,nPhi,lm)

    DO nRC=nRstart,nRstop
#ifdef WITHOMP
       nTh=OMP_GET_THREAD_NUM()+1
#else
       nTh=1
#endif
       IF( nTh.GT.nThreadsRmax ) nThreadsRmax=nTh
       IF ( lVerbose )                                              &
            &           WRITE(*,'(/," ! Starting radial level ",i4)') nR
       IF ( lVerbose )                                              &
            &           WRITE(*,'(" ! using thread no:",i4)') nTh

       IF ( nRC.EQ.n_r_cmb ) THEN 
          IF ( lOutBc ) THEN
             nRC2stop=2 ! Do outer and inner boundary in one thread !
          ELSE
             GOTO 999   ! Nothing needs to be done by thread one !
          END IF
       ELSE
          nRC2stop=1
       END IF

       DO nRC2=1,nRC2stop 

          IF ( nRC.EQ.n_r_cmb ) THEN 
             IF ( nRC2.EQ.1 ) THEN
                nR=n_r_cmb
                nBc=ktopv
             ELSE
                nR=n_r_icb
                nBc=kbotv
             ENDIF
             lDeriv= lTOCalc .OR. lHelCalc .OR. l_frame 
          ELSE 
             nR=nRC
             nBc=0
             lDeriv=.TRUE.
          END IF
          IF ( l_cour ) THEN
             dtrkc(nR)=1.D10
             dthkc(nR)=1.D10
          END IF
          IF ( lTOCalc ) THEN
             !------ Zero lm coeffs for first theta block:
             DO l=0,l_max
                dzRstrLM(l+1)=0.D0
                dzAstrLM(l+1)=0.D0
                dzCorLM(l+1) =0.D0
                dzLFLM(l+1)  =0.D0
             END DO
          END IF

          !----- Prepare legendre transform:
          !      legPrepG collects all the different modes necessary 
          !      to calculate the non-linear terms at a radial grid point nR
          CALL legPrepG(nR,nBc,lDeriv,lRmsCalc,l_frame,             &
               &                            lTOnext,lTOnext2,lTOcalc,             &
               &               dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC,             &
               &                         dLhb,dLhj,bhG,bhC,cbhG,cbhC,             &
               &           sR,dsR,preR,dpR,zAS,dzAS,ddzAS,bCMB,omegaIC,omegaMA)

          !----- Blocking of loops over ic (theta):
          DO nThetaB=1,nThetaBs

             nThetaLast =(nThetaB-1)*sizeThetaB
             nThetaStart=nThetaLast+1
             nThetaStop =nThetaLast+sizeThetaB

             !----- Legendre transform from (r,l,m) to (r,theta,m):
             !      First version with PlmTF needed for first-touch policy  
             IF ( l_mag ) THEN
                CALL legTFG(nBc,lDeriv,nThetaStart,                 &
                     &                        vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc,    &
                     &                                  dvrdtc,dvrdpc,dvtdpc,dvpdpc,    &
                     &                           brc,btc,bpc,cbrc,cbtc,cbpc,sc,drSc,    &
                     &                               dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,    &
                     &                    dvhdrC,dLhb,dLhj,bhG,bhC,cbhG,cbhC,sR,dsR)
             ELSE
                CALL legTFGnomag(nBc,lDeriv,nThetaStart,            &
                     &                        vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc,    &
                     &                          dvrdtc,dvrdpc,dvtdpc,dvpdpc,sc,drSc,    &
                     &                 dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC,sR,dsR)
             END IF

             !------ Fourier transform from (r,theta,m) to (r,theta,phi):
             IF ( l_conv .OR. l_mag_kin ) THEN
                IF ( l_heat ) CALL fft_thetab(sc,1)
                IF ( l_HT ) CALL fft_thetab(drSc,1)
                IF ( nBc.EQ.0 ) THEN
                   CALL fft_thetab(vrc,1)
                   CALL fft_thetab(vtc,1)
                   CALL fft_thetab(vpc,1)
                   IF ( lDeriv ) THEN
                      CALL fft_thetab(dvrdrc,1)
                      CALL fft_thetab(dvtdrc,1)
                      CALL fft_thetab(dvpdrc,1)
                      CALL fft_thetab(cvrc,1)
                      CALL fft_thetab(dvrdtc,1)
                      CALL fft_thetab(dvrdpc,1)
                      CALL fft_thetab(dvtdpc,1)
                      CALL fft_thetab(dvpdpc,1)
                   END IF
                ELSE IF ( nBc.EQ.1 ) THEN ! Stress free
                   DO nTheta=1,nfs
                      DO nPhi=1,ncp
                         vrc(nPhi,nTheta)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                      END DO
                   END DO
                   CALL fft_thetab(vtc,1)
                   CALL fft_thetab(vpc,1)
                   IF ( lDeriv ) THEN
                      DO nTheta=1,nfs
                         DO nPhi=1,ncp
                            dvrdtc(nPhi,nTheta)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                            dvrdpc(nPhi,nTheta)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                         END DO
                      END DO
                      CALL fft_thetab(dvrdrc,1)
                      CALL fft_thetab(dvtdrc,1)
                      CALL fft_thetab(dvpdrc,1)
                      CALL fft_thetab(cvrc,1)
                      CALL fft_thetab(dvtdpc,1)
                      CALL fft_thetab(dvpdpc,1)
                   END IF
                ELSE IF ( nBc.EQ.2 ) THEN 
                   IF ( nR.EQ.n_r_cmb ) THEN
                      CALL v_rigid_boundary(nR,omegaMA,lDeriv,      &
                           &                                                vrc,vtc,vpc,      &
                           &                           cvrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc,      &
                           &                                                nThetaStart)
                   ELSE IF ( nR.EQ.n_r_icb ) THEN
                      CALL v_rigid_boundary(nR,omegaIC,lDeriv,      &
                           &                                                vrc,vtc,vpc,      &
                           &                           cvrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc,      &
                           &                                                nThetaStart)
                   END IF
                   IF ( lDeriv ) THEN
                      CALL fft_thetab(dvrdrc,1)
                      CALL fft_thetab(dvtdrc,1)
                      CALL fft_thetab(dvpdrc,1)
                   END IF
                END IF
             END IF
             IF ( l_mag .OR. l_mag_LF ) THEN
                CALL fft_thetab(brc,1)
                CALL fft_thetab(btc,1)
                CALL fft_thetab(bpc,1)
                IF ( lDeriv ) THEN
                   CALL fft_thetab(cbrc,1)
                   CALL fft_thetab(cbtc,1)
                   CALL fft_thetab(cbpc,1)
                END IF
             END IF

             !--------- Calculation of nonliner products in grid scape:
             IF ( nRC.NE.n_r_cmb .OR. lMagNlBc ) THEN 
                CALL get_nl(vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc,  &
                     &                                  dvrdtc,dvrdpc,dvtdpc,dvpdpc,    &
                     &                                  brc,btc,bpc,cbrc,cbtc,cbpc,sc,  &
                     &                                  Advr,Advt,Advp,LFr,LFt,LFp,     &
                     &                                  VSr,VSt,VSp,VxBr,VxBt,VxBp,     &
                     &                                  ViscHeat,OhmLoss,               &
                     &                                  nR,nBc,nThetaStart)
                IF ( nRC.NE.n_r_cmb .AND.                           &
                     &                   ( l_conv_nl .OR. l_mag_LF ) ) THEN
                   IF ( l_conv_nl .AND. l_mag_LF ) THEN
                      DO nTheta=1,sizeThetaB
                         DO nPhi=1,nrp
                            Advr(nPhi,nTheta)=Advr(nPhi,nTheta) +      &
                                 &                                          LFr(nPhi,nTheta)
                            Advt(nPhi,nTheta)=Advt(nPhi,nTheta) +      &
                                 &                                          LFt(nPhi,nTheta)
                            Advp(nPhi,nTheta)=Advp(nPhi,nTheta) +      &
                                 &                                          LFp(nPhi,nTheta)
                         END DO
                      END DO
                   ELSE IF ( l_mag_LF ) THEN
                      DO nTheta=1,sizeThetaB
                         DO nPhi=1,nrp
                            Advr(nPhi,nTheta)=LFr(nPhi,nTheta)
                            Advt(nPhi,nTheta)=LFt(nPhi,nTheta)
                            Advp(nPhi,nTheta)=LFp(nPhi,nTheta)
                         END DO
                      END DO
                   END IF
                   CALL fft_thetab(Advr,-1)
                   CALL fft_thetab(Advt,-1)
                   CALL fft_thetab(Advp,-1)
                   CALL legTF3(nThetaStart,AdvrLM,AdvtLM,AdvpLM,    &
                        &                                               Advr,Advt,Advp)
                   IF ( lRmsCalc .AND. l_mag_LF ) THEN ! LF treated extra:
                      CALL fft_thetab(LFr,-1)
                      CALL fft_thetab(LFt,-1)
                      CALL fft_thetab(LFp,-1)
                      CALL legTF3(nThetaStart,LFrLM,LFtLM,LFpLM,    &
                           &                                                  LFr,LFt,LFp)
                   END IF
                END IF
                IF ( nRC.NE.n_r_cmb .AND. l_heat ) THEN
                   CALL fft_thetab(VSr,-1)
                   CALL fft_thetab(VSt,-1)
                   CALL fft_thetab(VSp,-1)
                   CALL legTF3(nThetaStart,VSrLM,VStLM,VSpLM,       &
                        &                                               VSr,VSt,VSp)
                   IF (l_anel) THEN ! anelastic stuff 
                      IF (l_mag_nl) THEN
                         CALL fft_thetab(ViscHeat,-1)
                         CALL fft_thetab(OhmLoss,-1)
                         CALL legTF2(nThetaStart,OhmLossLM,         &
                              &                                   ViscHeatLM,OhmLoss,            &
                              &                                   ViscHeat)
                      ELSE
                         CALL fft_thetab(ViscHeat,-1)
                         CALL legTF1(nThetaStart,ViscHeatLM,        &
                              &                                   ViscHeat)
                      END IF
                   END IF
                END IF
                IF ( l_mag_nl ) THEN
                   IF ( nRC.NE.n_r_cmb ) THEN
                      CALL fft_thetab(VxBr,-1)
                      CALL fft_thetab(VxBt,-1)
                      CALL fft_thetab(VxBp,-1)
                      CALL legTF3(nThetaStart,VxBrLM,VxBtLM,VxBpLM, &
                           &                                                  VxBr,VxBt,VxBp)
                   ELSE
                      CALL fft_thetab(VxBt,-1)
                      CALL fft_thetab(VxBp,-1)
                      CALL legTF2(nThetaStart,VxBtLM,VxBpLM,        &
                           &                                                VxBt,VxBp)
                   END IF
                END IF
             ELSE IF ( l_mag ) THEN
                DO lm=1,lmP_max
                   VxBtLM(lm)=0.D0
                   VxBpLM(lm)=0.D0
                END DO
             END IF

             !--------- Calculation of nonlinear products needed for conducting mantle or
             !          conducting inner core if free stress BCs are applied:
             !          input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
             !          ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
             !          and br_vp_lm_cmb in lm-space, respectively the contribution
             !          to these products from the points theta(nThetaStart)-theta(nThetaStop)
             !          These products are used in get_b_nl_bcs.
             IF ( nR.EQ.n_r_cmb .AND. l_b_nl_cmb ) THEN
                CALL get_br_v_bcs(brc,vtc,vpc,omegaMA,              &
                     &                     or2(nR),orho1(nR),nThetaStart,sizeThetaB,    &
                     &                          br_vt_lm_cmb,br_vp_lm_cmb)
             ELSE IF ( nR.EQ.n_r_icb .AND. l_b_nl_icb ) THEN
                CALL get_br_v_bcs(brc,vtc,vpc,omegaIC,              &
                     &                     or2(nR),orho1(nR),nThetaStart,sizeThetaB,    &
                     &                          br_vt_lm_icb,br_vp_lm_icb)
             END IF

             !--------- Calculate Lorentz torque on inner core:
             !          each call adds the contribution of the theta-block to
             !          lorentz_torque_ic
             IF ( nR.EQ.n_r_icb .AND. l_mag_LF .AND.                &
                  &                l_rot_ic .AND. l_cond_ic  )                       &
                  &               CALL get_lorentz_torque(lorentz_torque_ic,         &
                  &                                  nThetaStart,sizeThetaB,         &
                  &                                              brc,bpc,nR)

             !--------- Calculate Lorentz torque on mantle:
             !          note: this calculates a torque of a wrong sign.
             !          sign is reversed at the end of the theta blocking.
             IF ( nR.EQ.n_r_cmb .AND. l_mag_LF .AND.                &
                  &                l_rot_ma .AND. l_cond_ma )                        &
                  &              CALL get_lorentz_torque(lorentz_torque_ma,          &
                  &                                 nThetaStart,sizeThetaB,          &
                  &                                             brc,bpc,nR)

             !--------- Calculate courant condition parameters:
             IF ( l_cour )                                          &
                  &              CALL courant(nR,dtrkc(nR),dthkc(nR),                &
                  &                          vrc,vtc,vpc,brc,btc,bpc,                &
                  &                           nThetaStart,sizeThetaB)

             !--------- Since the fields are given at gridpoints here, this is a good
             !          point for graphical output:
             IF ( l_graph )                                         &
                  &              CALL graphOut(time,nR,ngform,vrc,vtc,vpc,           &
                  &                                        brc,btc,bpc,sc,           &
                  &                   nThetaStart,sizeThetaB,lGraphHeader)

             !--------- Helicity output:
             IF ( lHelCalc )                                        &
                  &              CALL getHelLM(vrc,vtc,vpc,                          &
                  &                  cvrc,dvrdtc,dvrdpc,dvtdrc,dvpdrc,               &
                  &                        HelLMr(1,nR),Hel2LMr(1,nR),               &
                  &                    HelnaLMr(1,nR),Helna2LMr(1,nR),               &
                  &                                    nR,nThetaStart)

             !--------- horizontal velocity :

             IF ( l_viscBcCalc )                                          &
                  &       CALL get_duHorizontal(vtc,vpc,dvtdrc,dvpdrc,    &
                              uhLMr(1,nR),duhLMr(1,nR),nR,nThetaStart)
             !--------- Movie output:
             IF ( l_frame .AND. l_movie_oc .AND.                    &
                  &                l_store_frame )                                   &
                  &              CALL store_movie_frame(nR,vrc,vtc,vpc,              &
                  &                                brc,btc,bpc,sc,drSc,              &
                  &                   dvrdpc,dvpdrc,dvtdrc,dvrdtc,cvrc,              &
                  &              cbrc,cbtc,nThetaStart,sizeThetaB,bCMB)


             !--------- Stuff for special output:
             !--------- Calculation of magnetic field production and advection terms
             !          for graphic output:
             IF ( l_dtB )                                           &
                  &              CALL get_dtBLM(nR,vrc,vtc,vpc,brc,btc,bpc,          &
                  &                              nThetaStart,sizeThetaB,             &
                  &           BtVrLM,BpVrLM,BrVtLM,BrVpLM,BtVpLM,BpVtLM,             &
                  &         BrVZLM,BtVZLM,BtVpCotLM,BpVtCotLM,BtVZcotLM,             &
                  &                       BtVpSn2LM,BpVtSn2LM,BtVZsn2LM)


             !--------- Torsional oscillation terms:
             IF ( ( lTONext .OR. lTONext2 ) .AND. l_mag )                            &
                  &              CALL getTOnext(zAS,brc,btc,bpc,lTONext,             &
                  &            lTONext2,dt,dtLast,nR,nThetaStart,sizeThetaB,         &
                  &                                    BsLast,BpLast,BzLast)

             IF ( lTOCalc )                                         &
                  &              CALL getTO(vrc,vtc,vpc,cvrc,dvpdrc,                 &
                  &                           brc,btc,bpc,cbrc,cbtc,                 &
                  &                            BsLast,BpLast,BzLast,                 &
                  &                dzRstrLM,dzAstrLM,dzCorLM,dzLFLM,                 &
                  &                dtLast,nR,nThetaStart,sizeThetaB)

          END DO ! Loop over theta blocks


          IF ( l_mag .OR. l_mag_LF ) THEN
             nR_Mag=nR
          ELSE
             nR_mag=1
          END IF
          !-- Partial calculation of time derivatives (horizontal parts):
          !   input flm...  is in (l,m) space at radial grid points nR !
          !   Only dVxBh needed for boundaries !
          !   get_td finally calculates the d*dt terms needed for the 
          !   time step performed in s_LMLoop.f . This should be distributed
          !   over the different models that s_LMLoop.f parallizes over. 
          CALL get_td(nR,nBc,lRmsCalc,dVSrLM(1,nR),dVxBhLM(1,nR_Mag),             &
               &                              dwdt(1,nR),dzdt(1,nR),dpdt(1,nR),   &
               &                      dsdt(1,nR),dbdt(1,nR_Mag),djdt(1,nR_Mag),   &
               &                        AdvrLM,AdvtLM,AdvpLM,LFrLM,LFtLM,LFpLM,   &
               &                        VSrLM,VStLM,VSpLM,VxBrLM,VxBtLM,VxBpLM,   &
               &              ViscHeatLM,OhmLossLM,dLhw,dLhdw,dLhz,sR,preR,dpR)

          !-- Finish calculation of TO variables:
          IF ( lTOcalc )                                            &
               &           CALL getTOfinish(nR,dtLast,zAS,dzAS,ddzAS,             &
               &               dzRstrLM,dzAstrLM,dzCorLM,dzLFLM)

          !--- Form partial horizontal derivaties of magnetic production and
          !    advection terms:
          IF ( l_dtB )                                              &
               &           CALL get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,             &
               &        BrVpLM,BtVpLM,BpVtLM,BrVZLM,BtVZLM,BtVpCotLM,             &
               &                       BpVtCotLM,BtVZcotLM,BtVpSn2LM,             &
               &                                 BpVtSn2LM,BtVZsn2LM,             &
               &           TstrRLM(1,nR),TadvRLM(1,nR),TomeRLM(1,nR))

       END DO ! Loop over radial sub level n_r_cmb,n_r_icb

999    CONTINUE ! JUNP POINT OF NOTHING NEEDS TO BE DONE 

    END DO    ! Loop over radial levels 

    !$OMP END PARALLEL DO   ! END OF SMP PARALLEL LOOP OVER RADIAL LEVELS !

    !----- Correct sign of mantel Lorentz torque (see above):
    lorentz_torque_ma=-lorentz_torque_ma


  END SUBROUTINE radialLoopG
END MODULE radialLoop
