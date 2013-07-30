!$Id$
#include "perflib_preproc.cpp"

MODULE radialLoop
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
  USE parallel_mod, ONLY: rank, n_procs
  USE radial_data,ONLY: nRstart,nRstop,n_r_cmb
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
  IMPLICIT NONE

  PRIVATE
  !---- Nonlinear terms, field components and help arrays
  !     for Legendre transform, field components, and dtB and 
  !     TO output. These are all stored in COMMON BLOCKS
  !     that are thread privat rather than using explicit 
  !     PRIVATE statements in the PARALLEL DO clause.

  !----- Nonlinear terms in lm-space: 
  COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: AdvrLM, AdvtLM, AdvpLM
  COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: LFrLM,  LFtLM,  LFpLM
  COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: VxBrLM, VxBtLM, VxBpLM
  COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: VSrLM,  VStLM,  VSpLM
  COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: ViscHeatLM, OhmLossLM
  !$OMP THREADPRIVATE( AdvrLM,AdvtLM,AdvpLM,LFrLM,LFtLM,LFpLM )
  !$OMP THREADPRIVATE( VxBrLM,VxBtLM,VxBpLM,VSrLM,VStLM,VSpLM )
  !$OMP THREADPRIVATE( ViscHeatLM,OhmLossLM )

  !----- Nonlinear terms in phi/theta space: 
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: Advr, Advt, Advp
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: LFr,  LFt,  LFp
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: VxBr, VxBt, VxBp
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: VSr,  VSt,  VSp
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: ViscHeat, OhmLoss
  !$OMP THREADPRIVATE( Advr,Advt,Advp,LFr,LFt,LFp )
  !$OMP THREADPRIVATE( VxBr,VxBt,VxBp,VSr,VSt,VSp )
  !$OMP THREADPRIVATE( ViscHeat,OhmLoss )

  !----- Fields calculated from these help arrays by legtf:
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: vrc, vtc, vpc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvrdrc, dvtdrc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvpdrc, cvrc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvrdtc, dvrdpc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvtdpc, dvpdpc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: brc, btc, bpc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: cbrc, cbtc, cbpc
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: sc, drSc
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

  !$OMP THREADPRIVATE( BtVrLM,BpVrLM,BrVtLM )
  !$OMP THREADPRIVATE( BrVpLM,BtVpLM,BpVtLM,BtVpCotLM )
  !$OMP THREADPRIVATE( BpVtCotLM,BtVpSn2LM,BpVtSn2LM )
  !$OMP THREADPRIVATE( BrVZLM,BtVZLM,BtVZcotLM,BtVZsn2LM )

  !----- Local TO output stuff:
  REAL(kind=8),ALLOCATABLE :: dzRstrLM(:),dzAstrLM(:)
  REAL(kind=8),ALLOCATABLE :: dzCorLM(:),dzLFLM(:)

  !$OMP THREADPRIVATE( dzRstrLM,dzAstrLM,dzCorLM,dzLFLM )


  !----- Saved magnetic field components from last time step:
  !      This is needed for the current TO version. However,
  !      the variables calulated with this don't give any
  !      deep insight. TO should be changes in the future to
  !      eliminate this.
  REAL(kind=8),ALLOCATABLE :: BsLast(:,:,:)
  REAL(kind=8),ALLOCATABLE :: BpLast(:,:,:)
  REAL(kind=8),ALLOCATABLE :: BzLast(:,:,:)

  ! public elements of the module

  PUBLIC :: initialize_radialLoop,radialLoopG

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
    btc=1.0d50
    bpc=1.0d50
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
       &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,dtrkc,dthkc)
    !***********************************************************************

    !    !------------ This is release 2 level 10  --------------!
    !    !------------ Created on 2/5/02  by JW. -----------

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine performs the actual time-stepping.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+

    !--- Input of variables:
    LOGICAL,intent(IN) :: l_graph,l_cour,l_frame
    LOGICAL,intent(IN) :: lTOcalc,lTONext,lTONext2,lHelCalc
    LOGICAL,intent(IN) :: lRmsCalc
    REAL(kind=8),intent(IN) :: time,dt,dtLast

    !---- Output of explicit time step:
    !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
    !     need a further treatment (radial derivatives required):
    COMPLEX(kind=8),INTENT(OUT),DIMENSION(lm_max,nRstart:nRstop) :: dwdt,dzdt,dpdt,dsdt,dVSrLM
    COMPLEX(kind=8),INTENT(OUT),DIMENSION(lm_maxMag,nRstartMag:nRstopMag) :: dbdt,djdt,dVxBhLM
    REAL(kind=8),intent(OUT) :: lorentz_torque_ma,lorentz_torque_ic

    !---- Output for axisymmetric helicity:
    REAL(kind=8),INTENT(OUT),DIMENSION(l_max+1,nRstart:nRstop) :: HelLMr,Hel2LMr,HelnaLMr,Helna2LMr
    REAL(kind=8),INTENT(OUT),DIMENSION(l_max+1,nRstart:nRstop) :: uhLMr,duhLMr

    !---- Output of nonlinear products for nonlinear
    !     magnetic boundary conditions (needed in s_updateB.f):
    COMPLEX(kind=8),intent(OUT) :: br_vt_lm_cmb(lmP_max)    ! product br*vt at CMB
    COMPLEX(kind=8),intent(OUT) :: br_vp_lm_cmb(lmP_max)    ! product br*vp at CMB
    COMPLEX(kind=8),intent(OUT) :: br_vt_lm_icb(lmP_max)    ! product br*vt at ICB
    COMPLEX(kind=8),intent(OUT) :: br_vp_lm_icb(lmP_max)    ! product br*vp at ICB

    !---- Output for Courant criteria:
    REAL(kind=8),INTENT(OUT) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)

    !--- Local variables:


    !---- Counter, logicals ...
    INTEGER :: nR!,nRC,nRC2,nRC2stop
    INTEGER :: nR_Mag
    INTEGER :: nBc,nPhi,nTh,lm,l
    INTEGER :: nTheta,nThetaB,nThetaLast
    INTEGER :: nThetaStart,nThetaStop
    LOGICAL :: lDeriv
    LOGICAL :: lOutBc
    LOGICAL :: lMagNlBc
    LOGICAL :: lGraphHeader    ! Write header into graph file
    logical :: isRadialBoundaryPoint

    !---- End of declaration
    !----------------------------------------------------------------------
    PERFON('rloop')

    lGraphHeader=l_graph
    IF ( lGraphHeader ) THEN
#ifdef WITH_MPI
       CALL graphOut_mpi(time,nR,ngform,vrc,vtc,vpc,brc,btc,bpc,sc,&
            &        nThetaStart,sizeThetaB,lGraphHeader)
#else
       CALL graphOut(time,nR,ngform,vrc,vtc,vpc,brc,btc,bpc,sc,&
            &        nThetaStart,sizeThetaB,lGraphHeader)
#endif
    END IF

    IF ( l_cour ) THEN
       IF (rank.EQ.0) THEN
          dtrkc(n_r_cmb)=1.D10
          dthkc(n_r_cmb)=1.D10
       ELSEIF (rank.EQ.n_procs-1) THEN
          dtrkc(n_r_icb)=1.D10
          dthkc(n_r_icb)=1.D10
       END IF
    END IF

    !------ Set nonlinear terms that are possibly needed at the boundaries.
    !       They may be overwritten by get_td later.
    DO lm=1,lm_max
       IF (rank.EQ.0) THEN
          dVSrLM(lm,n_r_cmb) =zero
          IF ( l_mag ) THEN
             dVxBhLM(lm,n_r_cmb)=zero
          END IF
       ELSEIF (rank.EQ.n_procs-1) then
          dVSrLM(lm,n_r_icb) =zero
          IF ( l_mag ) THEN
             dVxBhLM(lm,n_r_icb)=zero
          END IF
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

    !nRstart=n_r_cmb
    !nRstop =n_r_icb-1

    !--- Start the big do loop over the radial threads:

    !$OMP PARALLEL &
    !$OMP PRIVATE(nTh,nR,nR_Mag,nBc,lDeriv,l,nThetaB,isRadialBoundaryPoint)   &
    !$OMP private(nThetaLast,nThetaStart,nThetaStop,nTheta,nPhi,lm)

#ifdef WITHOMP
    !$OMP MASTER
    nThreadsRmax=omp_get_num_threads()
    !$OMP END MASTER
#else
    nThreadsRmax=1
#endif
    !$OMP DO SCHEDULE(STATIC,sizeRB) 
    DO nR=nRstart,nRstop
#ifdef WITHOMP
       nTh=OMP_GET_THREAD_NUM()+1
#else
       nTh=1
#endif
       !IF( nTh.GT.nThreadsRmax ) nThreadsRmax=nTh
       IF ( lVerbose ) THEN
          WRITE(*,'(/," ! Starting radial level ",i4)') nR
          WRITE(*,'(" ! using thread no:",i4)') nTh
       END IF

       !nR = nRC
       nBc = 0
       lDeriv = .true.
       isRadialBoundaryPoint=(nR.EQ.n_r_cmb).OR.(nR.EQ.n_r_icb)

       IF ( nR.EQ.n_r_cmb ) THEN 
          IF ( lOutBc ) THEN
             !nR  = n_r_cmb
             nBc = ktopv
             lDeriv= lTOCalc .OR. lHelCalc .OR. l_frame 
          ELSE
             CYCLE   ! Nothing needs to be done by thread one !
          END IF
       ELSEif ( nR.eq.n_r_icb ) then
          IF ( lOutBc ) THEN
             !nR = n_r_icb
             nBc = kbotv
             lDeriv= lTOCalc .OR. lHelCalc .OR. l_frame 
          ELSE
             CYCLE
          END IF
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
       PERFON('legPrepG')
       CALL legPrepG(nR,nBc,lDeriv,lRmsCalc,l_frame,             &
            &        lTOnext,lTOnext2,lTOcalc,             &
            &        dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC,             &
            &        dLhb,dLhj,bhG,bhC,cbhG,cbhC,             &
            &        sR,dsR,preR,dpR,zAS,dzAS,ddzAS,bCMB,omegaIC,omegaMA)
       PERFOFF

       !----- Blocking of loops over ic (theta):
       DO nThetaB=1,nThetaBs

          nThetaLast =(nThetaB-1)*sizeThetaB
          nThetaStart=nThetaLast+1
          nThetaStop =nThetaLast+sizeThetaB

          !----- Legendre transform from (r,l,m) to (r,theta,m):
          !      First version with PlmTF needed for first-touch policy  
          IF ( l_mag ) THEN
             PERFON('legTFG')
             CALL legTFG(nBc,lDeriv,nThetaStart,                &
                  &      vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc, &
                  &      dvrdtc,dvrdpc,dvtdpc,dvpdpc,           &
                  &      brc,btc,bpc,cbrc,cbtc,cbpc,sc,drSc,    &
                  &      dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,        &
                  &      dvhdrC,dLhb,dLhj,bhG,bhC,cbhG,cbhC,sR,dsR)
             PERFOFF
          ELSE
             PERFON('legTFGnm')
             CALL legTFGnomag(nBc,lDeriv,nThetaStart,                 &
                  &           vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc,  &
                  &           dvrdtc,dvrdpc,dvtdpc,dvpdpc,sc,drSc,    &
                  &           dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC,sR,dsR)
             PERFOFF
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
                      vrc(nPhi,nTheta)=CMPLX(0.D0,0.D0,kind=kind(vrc))
                   END DO
                END DO
                CALL fft_thetab(vtc,1)
                CALL fft_thetab(vpc,1)
                IF ( lDeriv ) THEN
                   DO nTheta=1,nfs
                      DO nPhi=1,ncp
                         dvrdtc(nPhi,nTheta)=CMPLX(0.D0,0.D0,kind=kind(dvrdtc))
                         dvrdpc(nPhi,nTheta)=CMPLX(0.D0,0.D0,kind=kind(dvrdpc))
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
                   CALL v_rigid_boundary(nR,omegaMA,lDeriv, &
                        &                vrc,vtc,vpc,&
                        &                cvrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc,      &
                        &                nThetaStart)
                ELSE IF ( nR.EQ.n_r_icb ) THEN
                   !write(*,"(I4,A,ES20.13)") nR,", omegaIC = ",omegaIC
                   CALL v_rigid_boundary(nR,omegaIC,lDeriv, &
                        &                vrc,vtc,vpc,      &
                        &                cvrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc,      &
                        &                nThetaStart)
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

          !--------- Calculation of nonlinear products in grid space:
          IF ( (.NOT.isRadialBoundaryPoint) .OR. lMagNlBc ) THEN 
             !write(*,"(I4,A,ES20.13)") nR,", vp = ",sum(real(conjg(vpc)*vpc))
             PERFON('get_nl')
             CALL get_nl(vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc,  &
                  &      dvrdtc,dvrdpc,dvtdpc,dvpdpc,    &
                  &      brc,btc,bpc,cbrc,cbtc,cbpc,sc,  &
                  &      Advr,Advt,Advp,LFr,LFt,LFp,     &
                  &      VSr,VSt,VSp,VxBr,VxBt,VxBp,     &
                  &      ViscHeat,OhmLoss,               &
                  &      nR,nBc,nThetaStart)
             PERFOFF
             IF ( (.not.isRadialBoundaryPoint) .AND. ( l_conv_nl .OR. l_mag_LF ) ) THEN
                IF ( l_conv_nl .AND. l_mag_LF ) THEN
                   DO nTheta=1,sizeThetaB
                      DO nPhi=1,nrp
                         Advr(nPhi,nTheta)=Advr(nPhi,nTheta) + LFr(nPhi,nTheta)
                         Advt(nPhi,nTheta)=Advt(nPhi,nTheta) + LFt(nPhi,nTheta)
                         Advp(nPhi,nTheta)=Advp(nPhi,nTheta) + LFp(nPhi,nTheta)
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
             IF ( (.not.isRadialBoundaryPoint) .AND. l_heat ) THEN
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
                IF ( .not.isRadialBoundaryPoint ) THEN
                   CALL fft_thetab(VxBr,-1)
                   CALL fft_thetab(VxBt,-1)
                   CALL fft_thetab(VxBp,-1)
                   CALL legTF3(nThetaStart,VxBrLM,VxBtLM,VxBpLM, VxBr,VxBt,VxBp)
                ELSE
                   !write(*,"(I4,A,ES20.13)") nR,", VxBt = ",sum(VxBt*VxBt)
                   CALL fft_thetab(VxBt,-1)
                   CALL fft_thetab(VxBp,-1)
                   CALL legTF2(nThetaStart,VxBtLM,VxBpLM, VxBt,VxBp)
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
          PERFON('nl_cmb')
          IF ( nR.EQ.n_r_cmb .AND. l_b_nl_cmb ) THEN
             CALL get_br_v_bcs(brc,vtc,vpc,omegaMA,              &
                  &            or2(nR),orho1(nR),nThetaStart,sizeThetaB,    &
                  &            br_vt_lm_cmb,br_vp_lm_cmb)
          ELSE IF ( nR.EQ.n_r_icb .AND. l_b_nl_icb ) THEN
             CALL get_br_v_bcs(brc,vtc,vpc,omegaIC,              &
                  &            or2(nR),orho1(nR),nThetaStart,sizeThetaB,    &
                  &            br_vt_lm_icb,br_vp_lm_icb)
          END IF
          PERFOFF
          !--------- Calculate Lorentz torque on inner core:
          !          each call adds the contribution of the theta-block to
          !          lorentz_torque_ic
          PERFON('lorentz')
          IF ( nR.EQ.n_r_icb .AND. l_mag_LF .AND. l_rot_ic .AND. l_cond_ic  ) &
               & CALL get_lorentz_torque(lorentz_torque_ic,         &
               &                         nThetaStart,sizeThetaB,         &
               &                         brc,bpc,nR)

          !--------- Calculate Lorentz torque on mantle:
          !          note: this calculates a torque of a wrong sign.
          !          sign is reversed at the end of the theta blocking.
          IF ( nR.EQ.n_r_cmb .AND. l_mag_LF .AND. l_rot_ma .AND. l_cond_ma ) &
               & CALL get_lorentz_torque(lorentz_torque_ma,          &
               &                         nThetaStart,sizeThetaB,          &
               &                         brc,bpc,nR)

          PERFOFF
          !--------- Calculate courant condition parameters:
          IF ( l_cour ) THEN
             !PRINT*,"Calling courant with nR=",nR
             CALL courant(nR,dtrkc(nR),dthkc(nR),   &
                  &       vrc,vtc,vpc,brc,btc,bpc,  &
                  &       nThetaStart,sizeThetaB)
          END IF

          !--------- Since the fields are given at gridpoints here, this is a good
          !          point for graphical output:
          IF ( l_graph ) THEN
             PERFON('graphout')
#ifdef WITH_MPI
             CALL graphOut_mpi(time,nR,ngform,vrc,vtc,vpc, &
                  &        brc,btc,bpc,sc,&
                  &        nThetaStart,sizeThetaB,lGraphHeader)
#else
             CALL graphOut(time,nR,ngform,vrc,vtc,vpc, &
                  &        brc,btc,bpc,sc,&
                  &        nThetaStart,sizeThetaB,lGraphHeader)
#endif
             PERFOFF
          END IF

          !--------- Helicity output:
          IF ( lHelCalc ) THEN
             PERFON('hel_out')
             CALL getHelLM(vrc,vtc,vpc,                          &
                  &                  cvrc,dvrdtc,dvrdpc,dvtdrc,dvpdrc,               &
                  &                        HelLMr(:,nR),Hel2LMr(:,nR),               &
                  &                    HelnaLMr(:,nR),Helna2LMr(:,nR),               &
                  &                                    nR,nThetaStart)
             PERFOFF
          END IF

          !--------- horizontal velocity :

          IF ( l_viscBcCalc ) THEN
             CALL get_duHorizontal(vtc,vpc,dvtdrc,dvpdrc,    &
                  &                uhLMr(1,nR),duhLMr(1,nR),nR,nThetaStart)
          END IF

          !--------- Movie output:
          IF ( l_frame .AND. l_movie_oc .AND. l_store_frame ) THEN
             PERFON('mov_out')
             CALL store_movie_frame(nR,vrc,vtc,vpc,                   &
                  &                 brc,btc,bpc,sc,drSc,              &
                  &                 dvrdpc,dvpdrc,dvtdrc,dvrdtc,cvrc, &
                  &                 cbrc,cbtc,nThetaStart,sizeThetaB,bCMB)
             PERFOFF
          END IF


          !--------- Stuff for special output:
          !--------- Calculation of magnetic field production and advection terms
          !          for graphic output:
          IF ( l_dtB ) THEN
             PERFON('dtBLM')
             CALL get_dtBLM(nR,vrc,vtc,vpc,brc,btc,bpc,          &
                  &                              nThetaStart,sizeThetaB,             &
                  &           BtVrLM,BpVrLM,BrVtLM,BrVpLM,BtVpLM,BpVtLM,             &
                  &         BrVZLM,BtVZLM,BtVpCotLM,BpVtCotLM,BtVZcotLM,             &
                  &                       BtVpSn2LM,BpVtSn2LM,BtVZsn2LM)
             PERFOFF
          END IF


          !--------- Torsional oscillation terms:
          PERFON('TO_terms')
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
          PERFOFF
       END DO ! Loop over theta blocks


       !-- Partial calculation of time derivatives (horizontal parts):
       !   input flm...  is in (l,m) space at radial grid points nR !
       !   Only dVxBh needed for boundaries !
       !   get_td finally calculates the d*dt terms needed for the 
       !   time step performed in s_LMLoop.f . This should be distributed
       !   over the different models that s_LMLoop.f parallelizes over. 
       IF ( l_mag .OR. l_mag_LF ) THEN
          nR_Mag=nR
       ELSE
          nR_Mag=1
       END IF
       !write(*,"(A,I4,2ES20.13)") "before_td: ",nR,sum(real(conjg(VxBtLM)*VxBtLM)),sum(real(conjg(VxBpLM)*VxBpLM))
       CALL get_td(nR,nBc,lRmsCalc,dVSrLM(:,nR),dVxBhLM(:,nR_Mag),   &
            &      dwdt(:,nR),dzdt(:,nR),dpdt(:,nR),   &
            &      dsdt(:,nR),dbdt(:,nR_Mag),djdt(:,nR_Mag),   &
            &      AdvrLM,AdvtLM,AdvpLM,LFrLM,LFtLM,LFpLM,   &
            &      VSrLM,VStLM,VSpLM,VxBrLM,VxBtLM,VxBpLM,   &
            &      ViscHeatLM,OhmLossLM,dLhw,dLhdw,dLhz,sR,preR,dpR)
       !write(*,"(A,I4,ES20.13)") "after_td:  ",nR,sum(real(conjg(dVxBhLM(:,nR_Mag))*dVxBhLM(:,nR_Mag)))
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
            &                                 BpVtSn2LM,BtVZsn2LM)

       !END DO ! Loop over radial sub level n_r_cmb,n_r_icb

    END DO    ! Loop over radial levels 

    !$OMP END DO
    !$OMP END PARALLEL   ! END OF SMP PARALLEL LOOP OVER RADIAL LEVELS !

    !----- Correct sign of mantel Lorentz torque (see above):
    lorentz_torque_ma=-lorentz_torque_ma

    PERFOFF
  END SUBROUTINE radialLoopG
END MODULE radialLoop
