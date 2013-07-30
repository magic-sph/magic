!$Id$
MODULE output_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE fields
  USE kinetic_energy
  USE magnetic_energy
  USE fields_average_mod
  USE output_data
  USE spectrum_average_mod
  USE spectrumC_average_mod
  USE outTO_mod
  USE outPV3
  USE const
  use outPar_mod, only: outPar
  USE power, ONLY: get_power

  IMPLICIT NONE

contains

  !***********************************************************************
  SUBROUTINE output(time,dt,dtNew,n_time_step,l_stop_time,        &
       &              l_Bpot,l_Vpot,l_Tpot,l_log,l_graph,lRmsCalc,        &
       &                                   l_store,l_new_rst_file,        &
       &                    l_spectrum,lTOCalc,lTOframe,lTOZwrite,        &
       &                         l_frame,n_frame,l_cmb,n_cmb_sets,        &
       &                 lorentz_torque_ic,lorentz_torque_ma,dbdt,        &
       &                                  TstrRLM,TadvRLM,TomeRLM,        &
       &           HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr)
    !***********************************************************************

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine controls most of the output.                     |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+

    !--- Input of variables
    USE charmanip, ONLY: dble2str
    USE omega, ONLY: outOmega
    USE integration, ONLY: rInt_R

    IMPLICIT NONE

    REAL(kind=8),INTENT(IN) :: time,dt,dtNew
    INTEGER,INTENT(IN) :: n_time_step
    LOGICAL,INTENT(IN) :: l_stop_time
    LOGICAL,INTENT(IN) :: l_Bpot,l_Vpot,l_Tpot
    LOGICAL,INTENT(IN) :: l_log, l_graph, lRmsCalc, l_store
    LOGICAL,INTENT(IN) :: l_new_rst_file, l_spectrum
    LOGICAL,INTENT(IN) :: lTOCalc,lTOframe,lTOZwrite
    LOGICAL,INTENT(IN) :: l_frame, l_cmb
    INTEGER,INTENT(INOUT) :: n_frame
    INTEGER,intent(IN) :: n_cmb_sets

    !--- Input of Lorentz torques and dbdt calculated in radialLoopG
    !    Parallelization note: Only the contribution at the CMB must be 
    !    collected and is (likely) stored on the processor (#0) that performs 
    !    this routine anyway.
    REAL(kind=8),intent(IN) :: lorentz_torque_ma,lorentz_torque_ic
    COMPLEX(kind=8),intent(IN) :: dbdt(lm_maxMag,n_r_maxMag)

    !--- Input of scales fields via common block in c_fields.f:
    !    Parallelization note: these fields are LM-distributed.
    !    The input fields HelLMr,Hel2LMr,TstrRLM,TadvRLM, and TomeRLM
    !    are R-distributed. More R-distributed fields are hidden 
    !    in c_TO.f, c_RMS.f, and c_dtB.f. 
    !    input fields are R-distributed. This has to be taken into
    !    account when collecting the information from the different
    !    processors!
    !    All the fields contained in c_fields.f are needed on
    !    the processor performing this routine:
    !          w,dw,ddw,z,dz,s,ds,p,b,db,ddb,aj,dj,ddj,
    !          b_ic,db_ic,ddb_ic,aj_ic,dj_ic,omega_ic,omega_ma
    !    omega_ic and omega_ma are likely located on processor #0 
    !    which deals with (l=1,m=0) in s_updateZ.f
    !    Note that many of these only have to be collected when
    !    certain output is required. This is controlled by the 
    !    input logicals.
    ! include 'c_fields.f'

    !--- Input help arrays for magnetic field stretching and advection and
    !    for calculating axisymmetric helicity.
    !    Parallelization note: These fields are R-distribute on input 
    !    and must also be collected on the processor performing this routine.
    COMPLEX(kind=8),intent(IN) :: TstrRLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8),intent(IN) :: TadvRLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8),intent(IN) :: TomeRLM(lm_max_dtB,n_r_max_dtB)
    REAL(kind=8),intent(IN) :: HelLMr(l_max+1,n_r_max)
    REAL(kind=8),intent(IN) :: Hel2LMr(l_max+1,n_r_max)
    REAL(kind=8),intent(IN) :: HelnaLMr(l_max+1,n_r_max)
    REAL(kind=8),intent(IN) :: Helna2LMr(l_max+1,n_r_max)
    REAL(kind=8),intent(IN) :: uhLMr(l_max+1,n_r_max)
    REAL(kind=8),intent(IN) :: duhLMr(l_max+1,n_r_max)

    !--- Local stuff:
    INTEGER :: nBpotSets
    INTEGER :: nVpotSets
    INTEGER :: nTpotSets
    SAVE nVpotSets,nBpotSets,nTpotSets

    !--- Energies:
    REAL(kind=8) :: ekinR(n_r_max)     ! kinetic energy w radius
    REAL(kind=8) :: e_mag,e_mag_ic,e_mag_cmb       
    REAL(kind=8) :: e_mag_p,e_mag_t      
    REAL(kind=8) :: e_mag_p_as,e_mag_t_as   
    REAL(kind=8) :: e_mag_p_ic,e_mag_t_ic   
    REAL(kind=8) :: e_mag_p_as_ic,e_mag_t_as_ic
    REAL(kind=8) :: e_mag_os,e_mag_as_os    
    REAL(kind=8) :: e_kin,e_kin_p,e_kin_t  
    REAL(kind=8) :: e_kin_p_as,e_kin_t_as 
    REAL(kind=8) :: eKinIC,eKinMA        
    REAL(kind=8) :: eTot,eTotOld        
    REAL(kind=8) :: dtE,dtEint         
    REAL(kind=8) :: e_kin_pMean,e_kin_tMean
    REAL(kind=8) :: e_mag_pMean,e_mag_tMean
    SAVE   e_kin_pMean,e_kin_tMean
    SAVE   e_mag_pMean,e_mag_tMean
    INTEGER :: n_e_sets
    SAVE n_e_sets
    INTEGER :: nRMS_sets
    SAVE nRMS_sets
    SAVE eTot,eTotOld,dtEint

    !--- Help arrays:
    COMPLEX(kind=8) :: dbdtCMB(lm_max)        ! SV at CMB !

    !-- Counter for output files/sets:
    INTEGER :: n_dt_cmb_sets
    INTEGER :: n_cmb_setsMov
    INTEGER :: n_v_r_sets(n_coeff_r_max) 
    INTEGER :: n_v_r_mov_sets(n_coeff_r_max) 
    INTEGER :: n_b_r_sets(n_coeff_r_max) 
    INTEGER :: n_b_r_mov_sets(n_coeff_r_max) 
    INTEGER :: n_spec
    INTEGER :: nPVsets
    SAVE n_v_r_sets,n_v_r_mov_sets
    SAVE n_b_r_sets,n_b_r_mov_sets
    SAVE n_spec,nPVsets

    !-- Further counter
    INTEGER :: nR,lm,n

    !--- For TO:
    INTEGER :: nTOsets,nTOmovSets,nTOrmsSets
    CHARACTER(len=64) :: TOfileNhs,TOfileShs,movFile
    CHARACTER(len=66) :: tayFile
    LOGICAL :: lTOrms    
    SAVE nTOsets,nTOmovSets,nTOrmsSets
    INTEGER :: nF1,nF2

    !--- For averaging:
    REAL(kind=8) :: timePassedLog
    REAL(kind=8) :: timeNormLog
    INTEGER :: nLogs  
    SAVE timePassedLog,timeNormLog,nLogs     

    !--- Property parameters:
    REAL(kind=8) :: dlBR(n_r_max),dlBRc(n_r_max),dlVR(n_r_max),dlVRc(n_r_max)
    REAL(kind=8) :: RolRu2(n_r_max),dlVRu2(n_r_max),dlVRu2c(n_r_max)

    REAL(kind=8) :: RmR(n_r_max)
    REAL(kind=8) :: Re,Ro,Rm,El,ElCmb,Rol,Geos,Dip,DipCMB!,ul,um
    REAL(kind=8) :: ReConv,RoConv,e_kin_nas,RolC
    REAL(kind=8) :: elsAnel
    REAL(kind=8) :: dlB,dlBc,dmB
    REAL(kind=8),save :: dlBMean,dmBMean
    REAL(kind=8) :: dlV,dlVc,dmV,dpV,dzV
    REAL(kind=8) :: visDiss,ohmDiss,lvDiss,lbDiss
    REAL(kind=8),save :: lvDissMean,lbDissMean
    REAL(kind=8),save :: RmMean,ElMean,ElCmbMean,RolMean,GeosMean
    REAL(kind=8),save :: DipMean,DipCMBMean
    REAL(kind=8),save :: dlVMean,dlVcMean,dmVMean,dpVMean,dzVMean
    REAL(kind=8) :: ReEquat
    INTEGER :: l,lm0

    LOGICAL :: l_r,l_PVout

    REAL(kind=8) :: timeScaled

    CHARACTER(len=76) :: filename
    CHARACTER(len=96) :: message

    CHARACTER(len=20) :: string

    !--- end of declaration
    !-----------------------------------------------------------------


    l_r= l_r_field .AND. l_cmb
    l_PVout=l_PV .AND. l_log
    timeScaled=tScale*time

    IF ( n_time_step.EQ.1 ) THEN
       DO n=1,n_coeff_r_max
          n_v_r_sets(n)    =0
          n_v_r_mov_sets(n)=0
          n_b_r_sets(n)    =0
          n_b_r_mov_sets(n)=0
       END DO
       n_spec       =0
       n_cmb_setsMov=0
       n_dt_cmb_sets=0
       nTOsets      =0
       nTOmovSets   =0
       nTOrmsSets   =0
       nBpotSets    =0
       nVpotSets    =0
       nTpotSets    =0
       n_e_sets     =0
       nRMS_sets    =0
    END IF

    !--- Get radial derivatives and add dt dtB terms:
    IF ( l_dtB ) CALL get_dtBLMfinish(time,n_time_step,             &
         &                    TstrRLM,TadvRLM,TomeRLM,omega_ic,             &
         &                   b,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic,             &
         &                                  aj_ic,dj_ic,ddj_ic)

    !----- Plot out inner core magnetic field, outer core
    !      field has been written in radialLoop !
    IF ( l_graph .AND. l_mag .AND. n_r_ic_max.gt.0 )                &
         &     CALL graphOut_IC(ngform,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,b)

    !--- Write spectra output that has partially been calculated in LMLoop
    IF ( l_frame .AND. l_rMagSpec .AND. n_time_step.GT.1 ) THEN
       CALL rBrSpec(time,b, b_ic ,'rBrSpecMov',.TRUE.)
       CALL rBpSpec(time,aj,aj_ic,'rBpSpecMov',.TRUE.)
    END IF
    IF ( l_log .AND. l_rMagSpec .AND. n_time_step.GT.1 ) THEN
       CALL rBrSpec(time,b, b_ic ,'rBrSpec',.TRUE.)
       CALL rBpSpec(time,aj,aj_ic,'rBpSpec',.TRUE.)
    END IF

    !--- Movie output and various supplementing things:
    IF ( l_frame ) THEN

       IF ( l_movie_ic .AND. l_store_frame )                        &
            &        CALL store_movie_frame_IC(b,b_ic,db_ic,ddb_ic,            &
            &                                          aj_ic,dj_ic)

       n_frame=n_frame+1
       CALL logWrite(' ')
       WRITE(message,'(1p,A,I8,A,D16.6,I8)') &
            & " ! WRITING MOVIE FRAME NO ",n_frame,&
            & "       at time/step",timeScaled,n_time_step
       CALL logWrite(message)

       !--- Storing the movie frame:
       CALL write_movie_frame(n_frame,timeScaled,                   &
            &             b,db,aj,dj,b_ic,db_ic,aj_ic,dj_ic,                   &
            &                             omega_ic,omega_ma)

       IF ( l_cmb_field )                                           &
            &        CALL write_Bcmb(timeScaled,b(1,n_r_cmb),lm_max,l_max,     &
            &                            l_max_cmb,minc,lm2,n_cmb_setsMov,     &
            &                                   cmbMov_file,n_cmbMov_file)

    END IF ! write movie frame ?


    !--- Store poloidal magnetic coeffs at cmb
    IF ( l_cmb )                                                    &
         &  CALL write_Bcmb(timeScaled,b(1,n_r_cmb),lm_max,l_max,           &
         &                         l_max_cmb,minc,lm2,n_cmb_sets,           &
         &                        cmb_file,n_cmb_file)

    !--- Store SV of poloidal magnetic coeffs at cmb
    IF ( l_dt_cmb_field .AND. l_cmb ) THEN
       nR=8! at CMB dbdt=induction=0, only diffusion !
       DO lm=1,lm_max
          dbdtCMB(lm)=                                              &
               &             dbdt(lm,n_r_cmb)/(dLh(lm)*or2(n_r_cmb)) +            &
               &                                      opm*hdif_B(lm) *            &
               &               (                     ddb(lm,n_r_cmb) -            &
               &                  dLh(lm)*or2(n_r_cmb)*b(lm,n_r_cmb) )
       END DO
       CALL write_Bcmb(timeScaled,dbdtCMB,lm_max,l_max,             &
            &                    l_max_cmb,minc,lm2,n_dt_cmb_sets,             &
            &                dt_cmb_file,n_dt_cmb_file)
    END IF

    !--- Store potential coeffs for velocity fields and magnetic fields
    IF ( l_r ) THEN
       DO n=1,n_coeff_r_max
          nR=n_coeff_r(n)
          CALL write_coeff_r(timeScaled,                            &
               &             w(1,nR),dw(1,nR),ddw(1,nR),z(1,nR),r(nR),            &
               &          lm_max,l_max,l_max_r,minc,lm2,n_v_r_sets(n),            &
               &          v_r_file(n),n_v_r_file(n),l_save_out,.TRUE.)
          IF ( l_mag ) THEN
             CALL write_coeff_r(timeScaled,                         &
                  &               b(1,nR),db(1,nR),ddb(1,nR),aj(1,nR),r(nR),         &
                  &             lm_max,l_max,l_max_r,minc,lm2,n_b_r_sets(n),         &
                  &            b_r_file(n),n_b_r_file(n),l_save_out,.FALSE.)
          END IF
       END DO
    END IF


    !--- Energies and rotation info and a lot of other stuff 
    !    performed for l_log=.TRUE.
    IF ( n_time_step.EQ.1 ) THEN 
       timeNormLog  =0.D0
       timePassedLog=0.D0
       nLogs        =0
       RmMean       =0.D0
       ElMean       =0.D0
       ElCmbMean    =0.D0
       RolMean      =0.D0
       GeosMean     =0.D0
       DipMean      =0.D0
       DipCMBMean   =0.D0
       e_kin_pMean  =0.D0
       e_kin_tMean  =0.D0
       e_mag_pMean  =0.D0
       e_mag_tMean  =0.D0
       dlVMean      =0.D0
       dlVcMean     =0.D0
       dmVMean      =0.D0
       dpVMean      =0.D0
       dzVMean      =0.D0
       dlBMean      =0.D0
       dmBMean      =0.D0
       lvDissmean   =0.D0
       lbDissmean   =0.D0
    END IF
    timePassedLog=timePassedLog+dt

    IF ( l_log ) THEN

       nLogs=nLogs+1
       timeNormLog=timeNormLog+timePassedLog

       !----- Write misc. output:
       CALL outMisc(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,   &
            &                       nLogs,w,dw,ddw,z,dz,s,ds,p,Geos,dpV,dzV)

       !----- Write torques and rotation rates:
       CALL write_rot(time,dt,eKinIC,eKinMA,w,z,dz,b,               &
            &                                 omega_ic,omega_ma,               &
            &               lorentz_torque_ic,lorentz_torque_ma)

       !----- Calculate energies (and write them to file in get_e_mag/get_e_kin):
       n_e_sets=n_e_sets+1
       CALL get_e_kin(time,.TRUE.,l_stop_time,n_e_sets,             &
            &                                              w,dw,z,             &
            &               e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,             &
            &                                      ekinR)
       e_kin=e_kin_p+e_kin_t
       e_kin_nas=e_kin-e_kin_p_as-e_kin_t_as

       CALL get_e_mag(time,.TRUE.,l_stop_time,n_e_sets,             &
            &                            b,db,aj,b_ic,db_ic,aj_ic,             &
            &                          e_mag_p,e_mag_t,e_mag_p_as,             &
            &                    e_mag_t_as,e_mag_p_ic,e_mag_t_ic,             &
            &                         e_mag_p_as_ic,e_mag_t_as_ic,             &
            &           e_mag_os,e_mag_as_os,e_mag_cmb,Dip,DipCMB,             &
            &                                             elsAnel)
       e_mag   =e_mag_p+e_mag_t
       e_mag_ic=e_mag_p_ic+e_mag_t_ic

       IF ( l_average ) THEN
          CALL spectrum_average(nLogs,l_stop_time,                  &
               &           timePassedLog,timeNormLog,w,z,dw,'V')
          CALL spectrumC_average(nLogs,l_stop_time,                 &
               &                  timePassedLog,timeNormLog,s,ds)
          IF ( l_mag )                                              &
               &            CALL spectrum_average(nLogs,l_stop_time,              &
               &              timePassedLog,timeNormLog,b,aj,db,'B')
          CALL fields_average(nLogs,l_stop_time,                    &
               &                    timePassedLog,timeNormLog,                    &
               &                            omega_ic,omega_ma,                    &
               &                        w,z,s,b,aj,b_ic,aj_ic)
       END IF

       IF ( l_power ) THEN
          IF ( nLogs.GT.1 ) THEN
             filename='dtE.'//tag
             OPEN(99,file=filename,status='unknown',                &
                  &                POSITION='APPEND')
             eTotOld=eTot
             eTot   =e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
             dtE    =(eTot-eTotOld)/timePassedLog
             dtEint =dtEint+timePassedLog*(eTot-eTotOld)
             WRITE(99,'(D20.10,3D16.6)') time,dtE,                  &
                  &                    dtEint/timeNormLog,dtE/eTot
             CLOSE(99)
          ELSE
             eTot=e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
             dtEint=0.D0
          END IF
          CALL get_power(time,timePassedLog,timeNormLog,l_stop_time,              &
               &                                    omega_ic,omega_ma,            &
               &                  lorentz_torque_ic,lorentz_torque_ma,            &
               &                             w,ddw,z,dz,s,b,ddb,aj,dj,            &
               &             db_ic,ddb_ic,aj_ic,dj_ic,visDiss,ohmDiss)
       END IF


       !----- Getting the property parameters:
       CALL getDlm(w,dw,z,dlV,dlVR,dmV,dlVc,dlVRc,'V')
       Re=DSQRT(2.D0*e_kin/vol_oc)/DSQRT(mass)
       ReConv=DSQRT(2.D0*e_kin_nas/vol_oc)/DSQRT(mass)

       IF ( l_non_rot ) THEN
          Ro=0.D0
          RoConv=0.D0
       ELSE
          Ro=Re*ek
          RoConv=ReConv*ek
       END IF

       !---- Surface zonal velocity at the equator
       IF ( ktopv==1 ) THEN
          ReEquat=0.d0
          DO l=1,l_max
             lm0=lm2(l,0)
             ReEquat=ReEquat-REAL(z(lm0,n_r_cmb))*dPl0Eq(l+1)*or1(n_r_cmb)
          END DO
       ELSE
          ReEquat=0.d0
       END IF

       IF ( l_mag .OR. l_mag_LF ) THEN 
          CALL getDlm(b,db,aj,dlB,dlBR,dmB,dlBc,dlBRc,'B')
       ELSE
          dlB=0.D0
          dmB=0.D0
       END IF

       IF ( dlV /= 0d0 ) THEN
          Rol=Ro/dlV   ! See Christensen&Aubert 2006, eqn.(27)
          RolC=RoConv/dlVc
       ELSE
          Rol=Ro
          RolC=RoConv
       END IF

       !----- If anelastic additional u**2 outputs
       IF ( l_anel) THEN
          CALL get_u_square(time,w,dw,z,RolRu2,dlVRu2,dlVRu2c)
       END IF

       !----- Radial properties
       CALL outPar(timePassedLog,timeNormLog,n_time_step,l_stop_time,    &
               &       ekinR,RolRu2,dlVR,dlVRc,dlVRu2,dlVRu2c,           &
               &       uhLMr,duhLMr,RmR)

       !--- Rm
       IF ( prmag.NE.0 .AND. nVarCond.GT.0 ) THEN
          Rm=0.d0
          Rm=rInt_R(RmR,n_r_max,n_r_max,drx,                        &
               &                                    i_costf_init,d_costf_init)
          Rm=Rm*3/(r_cmb**3-r_icb**3)
       ELSE IF ( prmag.NE.0 ) THEN
          Rm=Re*prmag
       ELSE
          Rm=Re
       END IF

       !--- El   =2.D0*e_mag/vol_oc/LFfac
       ! Elsasser number is computed from the averaged profile
       IF ( l_mag .OR. l_mag_LF ) THEN
          El   =elsAnel/vol_oc
          ElCmb=2.D0*e_mag_cmb/surf_cmb/LFfac
       ELSE
          El   =0d0
          ElCmb=0d0
       END IF
       IF ( l_power ) THEN
          IF ( visDiss /= 0d0 ) THEN
             lvDiss=DSQRT(e_kin/DABS(visDiss))            ! Viscous diffusion
          ELSE
             lvDiss=0d0
          END IF
          IF ( l_mag .OR. l_mag_LF ) THEN
             IF ( ohmDiss /= 0d0 ) THEN
                lbDiss=DSQRT((e_mag+e_mag_ic)/DABS(ohmDiss)) ! Ohmic diffusion 
             ELSE
                lbDiss=0d0
             END IF
          ELSE
             lbDiss=0.D0
          END IF
       ELSE
          lvDiss=0.D0
          lbDiss=0.D0
       END IF



       !----- Ouput into par file:
       IF ( l_save_out ) THEN
          OPEN(n_par_file,FILE=par_file,STATUS='UNKNOWN',           &
               &             POSITION='APPEND')
       END IF
       WRITE(n_par_file,'(D20.10,18D12.4)')    &
            &                   time,          &! 1) time
            &                     Rm,          &! 2) (magnetic) Reynolds number 
            &                     El,          &! 3) Elsasser number
            &                    Rol,          &! 4) local Rossby number
            &                   Geos,          &! 5) Geostrophy measure
            &                    Dip,          &! 6) Dipolarity
            &                 DipCMB,          &! 7) CMB dipolarity
            &        dlV,dmV,dpV,dzV,          &! 8,9,10,11) flow length scales
            &          lvDiss,lbDiss,          &! 12,13) dissipation length scales
            &                dlB,dmB,          &! 14,15) magnetic length scales
            &                  ElCmb,          &! 16) Elsasser number at CMB
            &                   RolC,          &! 17) Local Rol based on non-as flow
            &                   dlVc,          &! 18) convective flow length scale
            &                 ReEquat           ! 19) CMB flow at the equator
       IF ( l_save_out ) CLOSE(n_par_file)

       !---- Building time mean:
       RmMean     =RmMean     +timePassedLog*Rm
       ElMean     =ElMean     +timePassedLog*El
       ElCmbMean  =ElCmbMean  +timePassedLog*ElCmb
       RolMean    =RolMean    +timePassedLog*Rol
       GeosMean   =GeosMean   +timePassedLog*Geos
       DipMean    =DipMean    +timePassedLog*Dip
       DipCMBMean =DipCMBMean +timePassedLog*DipCMB
       e_kin_pMean=e_kin_pMean+timePassedLog*e_kin_p
       e_kin_tMean=e_kin_tMean+timePassedLog*e_kin_t
       e_mag_pMean=e_mag_pMean+timePassedLog*e_mag_p
       e_mag_tMean=e_mag_tMean+timePassedLog*e_mag_t
       dlVMean    =dlVMean    +timePassedLog*dlV   
       dlVcMean   =dlVcMean   +timePassedLog*dlVc
       !           dlVu2Mean  =dlVu2VMean +timePassedLog*dlVu2   
       !           dlVu2cMean =dlVu2cVMean+timePassedLog*dlVu2c   
       dmVMean    =dmVMean    +timePassedLog*dmV    
       dpVMean    =dpVMean    +timePassedLog*dpV
       dzVMean    =dzVMean    +timePassedLog*dzV
       lvDissMean =lvDissMean +timePassedLog*lvDiss
       lbDissMean =lbDissMean +timePassedLog*lbDiss
       dlBMean    =dlBMean    +timePassedLog*dlB
       dmBMean    =dmBMean    +timePassedLog*dmB

       IF ( l_stop_time ) THEN 

          !--- Time averaged parameters (properties)
          RmMean     =RmMean/timeNormLog
          ElMean     =ElMean/timeNormLog
          ElCmbMean  =ElCmbMean/timeNormLog
          RolMean    =RolMean/timeNormLog
          GeosMean   =GeosMean/timeNormLog 
          DipMean    =DipMean/timeNormLog  
          DipCMBMean =DipCMBMean/timeNormLog  
          e_kin_pMean=e_kin_pMean/timeNormLog
          e_kin_tMean=e_kin_tMean/timeNormLog
          e_mag_pMean=e_mag_pMean/timeNormLog
          e_mag_tMean=e_mag_tMean/timeNormLog 
          dlVMean    =dlVMean/timeNormLog
          dlVcMean   =dlVcMean/timeNormLog
          dmVMean    =dmVMean/timeNormLog
          dpVMean    =dpVMean/timeNormLog
          dzVMean    =dzVMean/timeNormLog
          dlBMean    =dlBMean/timeNormLog
          dmBMean    =dmBMean/timeNormLog
          lvDissMean =lvDissMean/timeNormLog
          lbDissMean =lbDissMean/timeNormLog

          CALL safeOpen(nLF,log_file)

          !--- Write end-energies including energy density:
          !    plus info on movie frames in to STDOUT and log-file
          WRITE(*,'(1p,/,A,/,A,/,A,4D16.6,/,A,4D16.6,/,A,4D16.6)') &
               & " ! Energies at end of time integration:", &
               & " !  (total,poloidal,toroidal,total density)",&
               & " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,&
               & " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc,&
               & " !  IC mag. energies:",e_mag_ic,e_mag_p_ic,e_mag_t_ic,e_mag_ic/vol_ic
          
          WRITE(nLF,'(1p,/,A,/,A,/,A,4D16.6,/,A,4D16.6,/,A,4D16.6)')&
               &          " ! Energies at end of time integration:",&
               &          " !  (total,poloidal,toroidal,total density)",&
               &          " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,&
               &          " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc,&
               &          " !  IC mag. energies:",e_mag_ic,e_mag_p_ic,e_mag_t_ic,e_mag_ic/vol_ic
          
          WRITE(nLF,'(1p,/,A,/,A,/,A,4D16.6,/,A,4D16.6)')&
               &          " ! Time averaged energies :",&
               &          " !  (total,poloidal,toroidal,total density)",&
               &          " !  Kinetic energies:",e_kin_pMean+e_kin_tMean,&
               &            e_kin_pMean,e_kin_tMean,(e_kin_pMean+e_kin_tMean)/vol_oc,&
               &          " !  OC mag. energies:",e_mag_pMean+e_mag_tMean,&
               &            e_mag_pMean,e_mag_tMean,&
               &            (e_mag_pMean+e_mag_tMean)/vol_oc

          WRITE(nLF,'(1p,/,A,7(/,A,D12.4),/,A,4D12.4,/,A,2D12.4,/,A,2D12.4)') &
               &           " ! Time averaged property parameters :",&
               &           " !  Rm (Re)         :",RmMean,&
               &           " !  Elsass          :",ElMean,&
               &           " !  Elsass at CMB   :",ElCmbMean,&
               &           " !  Rol             :",RolMean,&
               &           " !  Geos            :",GeosMean,&
               &           " !  Dip             :",DipMean,   &
               &           " !  DipCMB          :",DipCMBMean,&
               &           " !  l,m,p,z V scales:",dlVMean,dmVMean,dpVMean,dzVmean,&
               &           " !  l,m, B scales   :",dlBMean,dmBMean,&
               &           " !  vis, Ohm scale  :",lvDissMean,lbDissMean

          CALL safeClose(nLF)

       END IF ! l_stop_time ?

       timePassedLog=0.D0

    END IF

    IF ( l_spectrum ) THEN
       n_spec=n_spec+1
       CALL spectrum(time,n_spec,w,dw,z,                            &
            &             b,db,aj,b_ic,db_ic,aj_ic)
       CALL spectrumC(time,n_spec,s,ds)
    END IF

    IF ( l_Bpot )                                                   &
         &     CALL storePot(time,b,aj,b_ic,aj_ic,                          &
         &        nBpotSets,'Bpot.',omega_ma,omega_ic)
    IF ( l_Vpot )                                                   &
         &     CALL storePot(time,w,z,b_ic,aj_ic,                           &
         &        nVpotSets,'Vpot.',omega_ma,omega_ic)
    IF ( l_Tpot )                                                   &
         &     CALL storePot(time,s,z,b_ic,aj_ic,                           &
         &        nVpotSets,'Tpot.',omega_ma,omega_ic)


    IF ( lTOCalc ) THEN
       !------ Output for every log time step:
       IF ( lVerbose ) WRITE(*,*) '! Calling outTO !'
       TOfileNhs='TOnhs.'//tag
       TOfileShs='TOshs.'//tag
       movFile  ='TO_mov.'//tag
       tayFile  ='TaySphere4.'//tag
       nF1      =93
       nF2      =94
       lTOrms   =.TRUE.
       IF ( .NOT.l_log )                                            &
            &        CALL get_e_kin(time,.FALSE.,l_stop_time,0,                &
            &                                           w,dw,z,                &
            &            e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,                &
            &                                   ekinR)
       CALL outTO(time,n_time_step,e_kin,e_kin_t_as,                           &
            &                  nF1,nF2,TOfileNhs,TOfileShs,movFile,tayFile,    &
            &      nTOsets,nTOmovSets,nTOrmsSets,lTOframe,lTOrms,lTOZwrite,    &
            &                                          z,omega_ic,omega_ma)
       !------ Note: time averaging, time differencing done by IDL routine!

       IF ( lVerbose ) WRITE(*,*) '! outTO finished !'

    END IF

    IF ( l_RMS .AND. n_time_step.EQ.1 ) CALL zeroRms
    IF ( lRmsCalc ) THEN
       IF ( lVerbose ) WRITE(*,*) '! Writing RMS output !'
       CALL dtVrms(time,nRMS_sets)
       IF ( l_mag ) CALL dtBrms(time)
       CALL zeroRms
    END IF

    !----- Store current solution
    !      Note: unless l_new_rst_file=.TRUE. .and. .not.l_stop_time
    !            this is written into rst_end.TAG
    IF ( l_store ) THEN
       IF ( l_stop_time .OR. .NOT.l_new_rst_file ) THEN
          rst_file="rst_end."//tag
       ELSE IF ( l_new_rst_file ) THEN
          CALL dble2str(time,string)
          rst_file='rst_t='//trim(string)//'.'//tag
       END IF
       OPEN(n_rst_file,FILE=rst_file,                               &
            &          STATUS='UNKNOWN',FORM='UNFORMATTED')

       !------ Parallelization note:
       !       s_store.f stores the full solution in (lm,r) space.
       !       In addition to the scalar fields the time stepping arrays
       !       in c_dt_fieldsLast.f are also stored and have to be 
       !       collected on the processor executing this routine.
       CALL store(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic)
       CLOSE(n_rst_file)
       WRITE(*,'(/,1P,A,/,A,D20.10,/,A,I15,/,A,A)')&
            " ! Storing restart file:",&
            & "             at time=",time,&
            & "            step no.=",n_time_step,&
            & "           into file=",rst_file
       CALL safeOpen(nLF,log_file)

       WRITE(nLF,'(/,1P,A,/,A,D20.10,/,A,I15,/,A,A)') &
            & " ! Storing restart file:",&
            & "             at time=",time,&
            & "            step no.=",n_time_step,&
            & "           into file=",rst_file
       CALL safeClose(nLF)
    END IF

    IF ( l_SRIC .AND. l_stop_time ) CALL outOmega(z,omega_ic)

    !----- Output of axisymm. rotation rate for potential vorticity analysis:
    !  NOTE: For l_stop_time=.TRUE. outPV transforms the fields without 
    !        transforming them back. This must thus be the very last 
    !        thing done with them. 
    IF ( l_PVout ) CALL outPV(time,l_stop_time,nPVsets,             &
         &                     w,dw,ddw,z,dz,omega_ic,omega_ma)


    RETURN 
  END SUBROUTINE output


  !-----------------------------------------------------------------------
END MODULE output_mod
