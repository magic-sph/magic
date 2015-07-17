!$Id$
#include "intrinsic_sizes.h"
#include "perflib_preproc.cpp"
MODULE output_mod
  USE truncation
  USE radial_functions,ONLY: n_r_cmb,or1,or2,r,drx,i_costf_init,d_costf_init,r_cmb,r_icb
  USE radial_data,ONLY: nRstart,nRstop,nRstartMag,nRstopMag
  USE physical_parameters,ONLY: opm,ek,ktopv,prmag,nVarCond,LFfac
  USE num_param,only: tScale
  USE blocking,ONLY: st_map,lm2,lo_map
  USE horizontal_data,ONLY: dLh,hdif_B,dPl0Eq
  USE logic,ONLY: l_average,l_mag,l_power,l_anel,l_mag_LF,lVerbose,l_dtB, &
       & l_RMS,l_r_field,l_r_fieldT,l_PV,l_SRIC,l_cond_ic,l_rMagSpec,     &
       & l_movie_ic,l_store_frame,l_cmb_field,l_dt_cmb_field,             &
       & l_save_out,l_non_rot,l_perpPar
  USE fields,ONLY: omega_ic,omega_ma,b,db,ddb,aj,dj,ddj, &
       & b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic,           &
       & w,dw,ddw,z,dz,s,ds,p,                           &
       & w_LMloc,dw_LMloc,ddw_LMloc,p_LMloc,             &
       & s_LMloc,ds_LMloc,z_LMloc,dz_LMloc,              &
       & b_LMloc,db_LMloc,ddb_LMloc,                     &
       & aj_LMloc,dj_LMloc,ddj_LMloc,                    &
       & b_ic_LMloc,db_ic_LMloc,ddb_ic_LMloc,            &
       & aj_ic_LMloc,dj_ic_LMloc,ddj_ic_LMloc
  USE fieldsLast,ONLY: dwdtLast,dzdtLast,dpdtLast,dsdtLast,       &
       & dbdtLast,djdtLast,dbdt_icLast,djdt_icLast,               &
       & dwdtLast_LMloc,dzdtLast_lo,dpdtLast_LMloc,dsdtLast_LMloc,&
       & dbdtLast_LMloc,djdtLast_LMloc,dbdt_icLast_LMloc,         &
       & djdt_icLast_LMloc
  USE kinetic_energy,only: get_e_kin
  USE magnetic_energy,only: get_e_mag
  USE fields_average_mod,only: fields_average
  USE spectrum_average_mod,only: spectrum_average
  USE spectrumC_average_mod,only: spectrumC_average
  USE outTO_mod,only: outTO
  USE outPV3, only: outPV
  USE output_data,ONLY: tag,tag_wo_rank,ngform,l_max_cmb,cmbMov_file, &
       & n_cmbMov_file,cmb_file,n_cmb_file,dt_cmb_file,n_dt_cmb_file, &
       & n_coeff_r,l_max_r,n_v_r_file,n_b_r_file,n_t_r_file,          &
       & v_r_file,t_r_file,b_r_file,n_r_array,n_r_step,               &
       & par_file,n_par_file,nLF,log_file,n_coeff_r_max,rst_file,     &
       & n_rst_file
  USE const, ONLY: vol_oc,vol_ic,mass,surf_cmb
  USE parallel_mod
  USE outPar_mod, only: outPar
  USE outPerpPar_mod, only: outPerpPar
  USE power, ONLY: get_power
  USE LMLoop_data,ONLY: lm_per_rank,lm_on_last_rank,llm,ulm,llmMag,ulmMag
  USE communications,ONLY: myAllGather,gather_all_from_lo_to_rank0,   &
       & gt_OC,gt_IC
  USE write_special,only: write_Bcmb, write_coeff_r
  USE getDlm_mod,only: getDlm
  USE movie_data,only: movie_gather_frames_to_rank0
  USE storeCheckPoints

  IMPLICIT NONE

  PRIVATE

  INTEGER :: nBpotSets, nVpotSets, nTpotSets
  !-- Counter for output files/sets:
  INTEGER :: n_dt_cmb_sets, n_cmb_setsMov
  INTEGER,ALLOCATABLE :: n_v_r_sets(:), n_b_r_sets(:), n_T_r_sets(:)
  INTEGER :: n_spec,nPVsets

  INTEGER :: nTOsets,nTOmovSets,nTOrmsSets

  !--- For averaging:
  REAL(kind=8) :: timePassedLog, timeNormLog
  INTEGER :: nLogs  

  REAL(kind=8),SAVE :: dlBMean,dmBMean
  REAL(kind=8),SAVE :: lvDissMean,lbDissMean
  REAL(kind=8),SAVE :: RmMean,ElMean,ElCmbMean,RolMean,GeosMean
  REAL(kind=8),SAVE :: DipMean,DipCMBMean
  REAL(kind=8),SAVE :: dlVMean,dlVcMean,dmVMean,dpVMean,dzVMean

  REAL(kind=8) :: eTot,eTotOld,dtEint
  REAL(kind=8) :: e_kin_pMean, e_kin_tMean
  REAL(kind=8) :: e_mag_pMean, e_mag_tMean
  INTEGER :: n_e_sets, nRMS_sets


  PUBLIC :: output,initialize_output
contains

  SUBROUTINE initialize_output

    integer :: n

    IF ( l_r_field .OR. l_r_fieldT ) THEN
       ALLOCATE ( n_coeff_r(n_coeff_r_max))
       ALLOCATE ( n_v_r_file(n_coeff_r_max), v_r_file(n_coeff_r_max) )
       ALLOCATE ( n_v_r_sets(n_coeff_r_max) ) 
       n_v_r_sets=0

       IF ( l_mag ) THEN
          ALLOCATE ( n_b_r_file(n_coeff_r_max), b_r_file(n_coeff_r_max) )
          ALLOCATE ( n_b_r_sets(n_coeff_r_max) ) 
          n_b_r_sets=0
       END IF

       IF ( l_r_fieldT ) THEN
          ALLOCATE ( n_t_r_file(n_coeff_r_max), t_r_file(n_coeff_r_max) )
          ALLOCATE ( n_t_r_sets(n_coeff_r_max) ) 
          n_T_r_sets=0
       END IF

       IF ( COUNT(n_r_array>0)> 0 ) THEN
          n_coeff_r=n_r_array(1:n_coeff_r_max)
       ELSE
          n_r_step=MAX(n_r_step,1)
          DO n=1,n_coeff_r_max
             n_coeff_r(n)=n*n_r_step  ! used every n_r_step point !
          END DO
       END IF

    END IF

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
    nLogs        =0
    nRMS_sets    =0
    
    timeNormLog  =0.D0
    timePassedLog=0.D0
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
  END SUBROUTINE initialize_output

  !***********************************************************************
  SUBROUTINE output(time,dt,dtNew,n_time_step,l_stop_time,            &
       &            l_Bpot,l_Vpot,l_Tpot,l_log,l_graph,lRmsCalc,      &
       &            l_store,l_new_rst_file,                           &
       &            l_spectrum,lTOCalc,lTOframe,lTOZwrite,            &
       &            l_frame,n_frame,l_cmb,n_cmb_sets,l_r,             &
       &            lorentz_torque_ic,lorentz_torque_ma,dbdt_at_CMB,  &
       &            HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,   &
       &            gradsLMr,fconvLMr,fkinLMr,fviscLMr,fpoynLMr,      &
       &            fresLMr,EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)
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
    LOGICAL,INTENT(IN) :: lTOCalc,lTOframe
    LOGICAL,intent(INOUT) :: lTOZwrite
    LOGICAL,INTENT(IN) :: l_frame, l_cmb, l_r
    INTEGER,INTENT(INOUT) :: n_frame
    INTEGER,intent(INOUT) :: n_cmb_sets

    !--- Input of Lorentz torques and dbdt calculated in radialLoopG
    !    Parallelization note: Only the contribution at the CMB must be 
    !    collected and is (likely) stored on the processor (#0) that performs 
    !    this routine anyway.
    REAL(kind=8),intent(IN) :: lorentz_torque_ma,lorentz_torque_ic
    COMPLEX(kind=8),INTENT(IN),POINTER,dimension(:) :: dbdt_at_CMB

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
    REAL(kind=8),intent(IN) :: HelLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: Hel2LMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: HelnaLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: Helna2LMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: uhLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: gradsLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: duhLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: fconvLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: fkinLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: fviscLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
    REAL(kind=8),intent(IN) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)
    REAL(kind=8),intent(IN) :: EperpLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: EparLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: EperpaxiLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),intent(IN) :: EparaxiLMr(l_max+1,nRstart:nRstop)

    !--- Local stuff:
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
    REAL(kind=8) :: dtE

    !--- Help arrays:
    COMPLEX(kind=8) :: dbdtCMB(lm_max)        ! SV at CMB !

    INTEGER :: nR,lm,n

    !--- For TO:
    CHARACTER(len=64) :: TOfileNhs,TOfileShs,movFile
    CHARACTER(len=66) :: tayFile
    LOGICAL :: lTOrms    
    INTEGER :: nF1,nF2

    !--- Property parameters:
    REAL(kind=8) :: dlBR(n_r_max),dlBRc(n_r_max),dlVR(n_r_max),dlVRc(n_r_max)
    REAL(kind=8) :: RolRu2(n_r_max),dlVRu2(n_r_max),dlVRu2c(n_r_max)
    REAL(kind=8) :: RmR(n_r_max)
    REAL(kind=8) :: Re,Ro,Rm,El,ElCmb,Rol,Geos,Dip,DipCMB
    REAL(kind=8) :: ReConv,RoConv,e_kin_nas,RolC
    REAL(kind=8) :: elsAnel
    REAL(kind=8) :: dlB,dlBc,dmB
    REAL(kind=8) :: dlV,dlVc,dmV,dpV,dzV
    REAL(kind=8) :: visDiss,ohmDiss,lvDiss,lbDiss
    INTEGER :: l,lm0
    REAL(kind=8) :: ReEquat

    LOGICAL :: l_PVout

    REAL(kind=8) :: timeScaled

    CHARACTER(len=76) :: filename
    CHARACTER(len=96) :: message

    CHARACTER(len=20) :: string
    logical :: DEBUG_OUTPUT=.false.

    INTEGER :: length

    !--- end of declaration
    !-----------------------------------------------------------------

    timeScaled=tScale*time
    timePassedLog=timePassedLog+dt

    ! We start with the computation of the energies
    ! in parallel.
    IF (l_log) THEN
       nLogs=nLogs+1
       timeNormLog=timeNormLog+timePassedLog

       !----- Write torques and rotation rates:
       PERFON('out_rot')
       CALL write_rot( time,dt,eKinIC,eKinMA,w_LMloc,z_LMloc,dz_LMloc,b_LMloc,  &
            &          omega_ic,omega_ma,               &
            &          lorentz_torque_ic,lorentz_torque_ma)
       PERFOFF
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  write_rot  on rank ",rank

       PERFON('out_ekin')
       n_e_sets=n_e_sets+1
       CALL get_e_kin(time,.TRUE.,l_stop_time,n_e_sets,     &
            &         w_LMloc,dw_LMloc,z_LMloc,                &
            &         e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,&
            &         ekinR)
       e_kin=e_kin_p+e_kin_t
       !WRITE(*,"(A,3(I4,F20.17))") "e_kin, e_kin_p_as,e_kin_t_as = ",&
       !     &EXPONENT(e_kin),FRACTION(e_kin),&
       !     &EXPONENT(e_kin_p_as),FRACTION(e_kin_p_as),EXPONENT(e_kin_t_as),FRACTION(e_kin_t_as)
       e_kin_nas=e_kin-e_kin_p_as-e_kin_t_as
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  e_kin  on rank ",rank

       CALL get_e_mag(time,.TRUE.,l_stop_time,n_e_sets,             &
            &         b_LMloc,db_LMloc,aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc,  &
            &         e_mag_p,e_mag_t,e_mag_p_as,  &
            &         e_mag_t_as,e_mag_p_ic,e_mag_t_ic,  &
            &         e_mag_p_as_ic,e_mag_t_as_ic,  &
            &         e_mag_os,e_mag_as_os,e_mag_cmb,Dip,DipCMB,  &
            &         elsAnel)
       e_mag   =e_mag_p+e_mag_t
       e_mag_ic=e_mag_p_ic+e_mag_t_ic
       PERFOFF
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  e_mag  on rank ",rank

       IF (l_average) THEN
          PERFON('out_aver')
          CALL spectrum_average(nLogs,l_stop_time,                  &
               &                timePassedLog,timeNormLog,w_LMloc,z_LMloc,dw_LMloc,'V')
          CALL spectrumC_average(nLogs,l_stop_time,                 &
               &                 timePassedLog,timeNormLog,s_LMloc,ds_LMloc)

          IF ( l_mag ) THEN
             CALL spectrum_average(nLogs,l_stop_time, &
                  &                timePassedLog,timeNormLog,b_LMloc,aj_LMloc,db_LMloc,'B')
          END IF

          CALL fields_average(nLogs,l_stop_time,        &
               &              timePassedLog,timeNormLog,&
               &              omega_ic,omega_ma,        &
               &              w_LMloc,z_LMloc,s_LMloc,b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc)
          PERFOFF
          IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  averages  on rank ",rank
       END IF

       IF ( l_power ) THEN

          PERFON('out_pwr')
          IF (rank == 0) THEN
             IF ( nLogs > 1 ) THEN
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
                eTot   =e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
                dtEint=0.D0
             END IF
             !WRITE(*,"(A,7ES22.14)") "eTot = ",eTot,e_kin,e_mag,e_mag_ic,e_mag_os,eKinIC,eKinMA
          END IF
          CALL get_power( time,timePassedLog,timeNormLog,l_stop_time,      &
               &          omega_ic,omega_ma,                   &
               &          lorentz_torque_ic,lorentz_torque_ma, &
               &          w_LMloc,ddw_LMloc,z_LMloc,dz_LMloc,s_LMloc,b_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,&
               &          db_ic_LMloc,ddb_ic_LMloc,aj_ic_LMloc,dj_ic_LMloc,visDiss,ohmDiss)
          PERFOFF
          IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  power  on rank ",rank
       END IF

       !----- If anelastic additional u**2 outputs
       IF ( l_anel) THEN
          CALL get_u_square(time,w_LMloc,dw_LMloc,z_LMloc,RolRu2,dlVRu2,dlVRu2c)
          IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  u_square  on rank ",rank
       ELSE
          dlVRu2  = 0.0D0
          dlVRu2c = 0.0D0
       END IF

       IF ( l_perpPar ) THEN
          CALL outPerpPar(time,timePassedLog,timeNormLog,l_stop_time, &
                         &EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)
       END IF
       !----- Radial properties
       !WRITE(*,"(A,4ES20.12)") "before getDlm, w(n_r_icb,n_r_cmb): ",&
       !     & w_LMloc(n_r_icb),w_LMloc(n_r_cmb)
       !WRITE(*,"(A,4ES20.12)") "before getDlm, dw(n_r_icb,n_r_cmb): ",&
       !     & dw_LMloc(n_r_icb),dw_LMloc(n_r_cmb)
       !WRITE(*,"(A,4ES20.12)") "before getDlm, z(n_r_icb,n_r_cmb): ",&
       !     & z_LMloc(n_r_icb),z_LMloc(n_r_cmb)
       CALL getDlm(w_LMloc,dw_LMloc,z_LMloc,dlV,dlVR,dmV,dlVc,dlVRc,'V')
       !WRITE(*,"(A,ES20.12)") "dlVr,dlVrc(n_r_icb) = ",dlVr(n_r_icb),dlVrc(n_r_icb)
       !WRITE(*,"(A,ES20.12)") "dlVr,dlVrc(n_r_cmb) = ",dlVr(n_r_cmb),dlVrc(n_r_cmb)
       CALL outPar(timePassedLog,timeNormLog,nLogs,l_stop_time,    &
            &      ekinR,RolRu2,dlVR,dlVRc,dlVRu2,dlVRu2c,               &
            &      uhLMr,duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,      &
            &      fpoynLMr,fresLMr,RmR)
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  outPar  on rank ",rank

       !----- Write misc. output:
       CALL outMisc( timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,   &
            &        nLogs,w_LMloc,dw_LMloc,ddw_LMloc,z_LMloc,dz_LMloc,    &
            &        s_LMloc,ds_LMloc,p_LMloc,Geos,dpV,dzV)
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  outMisc  on rank ",rank

       IF ( l_mag .OR. l_mag_LF ) THEN 
          CALL getDlm(b_LMloc,db_LMloc,aj_LMloc,dlB,dlBR,dmB,dlBc,dlBRc,'B')
       ELSE
          dlB=0.D0
          dmB=0.D0
       END IF
    END IF

    IF ( l_spectrum ) THEN
       n_spec=n_spec+1
       CALL spectrum(time,n_spec,w_LMloc,dw_LMloc,z_LMloc,                            &
            &        b_LMloc,db_LMloc,aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc)
       CALL spectrumC(time,n_spec,s_LMloc,ds_LMloc)
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  spectrum  on rank ",rank
    END IF

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
       IF ( .NOT.l_log ) THEN
          CALL get_e_kin(time,.FALSE.,l_stop_time,0,           &
               &         w_LMloc,dw_LMloc,z_LMloc,                &
               &         e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,&
               &         ekinR)
       END IF
       CALL outTO(time,n_time_step,e_kin,e_kin_t_as,                      &
            &     nF1,nF2,TOfileNhs,TOfileShs,movFile,tayFile,            &
            &     nTOsets,nTOmovSets,nTOrmsSets,lTOframe,lTOrms,lTOZwrite,&
            &     z_LMloc,omega_ic,omega_ma)
       !------ Note: time averaging, time differencing done by IDL routine!

       IF ( lVerbose ) WRITE(*,*) '! outTO finished !'
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  TO  on rank ",rank
    END IF

    !--- Get radial derivatives and add dt dtB terms:
    IF ( l_dtB ) THEN
       CALL get_dtBLMfinish(time,n_time_step,                    &
            &               omega_ic,    &
            &               b_LMloc,ddb_LMloc,                   &
            &               aj_LMloc,dj_LMloc,ddj_LMloc,         &
            &               b_ic_LMloc,db_ic_LMloc,ddb_ic_LMloc, &
            &               aj_ic_LMloc,dj_ic_LMloc,ddj_ic_LMloc)
    END IF


    IF ( l_RMS .AND. n_time_step == 1 ) CALL zeroRms
    IF ( lRmsCalc ) THEN
       IF ( lVerbose ) WRITE(*,*) '! Writing RMS output !'
       CALL dtVrms(time,nRMS_sets)
       IF ( l_mag ) CALL dtBrms(time)
       !CALL zeroRms
       IF (DEBUG_OUTPUT) WRITE(*,"(A,I6)") "Written  dtV/Brms  on rank ",rank
    END IF

    !
    ! Parallel writing of the restart file (possible only when HDF5 is used)
    !
#ifdef WITH_HDF5
    if ( l_store ) then

       if ( l_stop_time .or. .not.l_new_rst_file ) then
          rst_file='h5_rst_end.'//tag_wo_rank
       else if ( l_new_rst_file ) then
          call dble2str(time,string)
          rst_file='h5_rst_t='//trim(string)//'.'//tag_wo_rank
       end if
       call storeHdf5_parallel(time,dt,dtNew,w_LMloc,z_LMloc,p_LMloc,s_LMloc, &
                               b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc,       &
                               dwdtLast_LMloc,dzdtLast_lo,dpdtLast_LMloc,     &
                               dsdtLast_LMloc,dbdtLast_LMloc,djdtLast_LMloc,  &
                               dbdt_icLast_LMloc,djdt_icLast_LMloc)

       if ( rank == 0 ) then
          write(*,'(/,1P,A,/,A,D20.10,/,A,I15,/,A,A)')&
               & " ! Storing restart file:",&
               & "             at time=",time,&
               & "            step no.=",n_time_step,&
               & "           into file=",rst_file
          call safeOpen(nLF,log_file)
          
          write(nLF,'(/,1P,A,/,A,D20.10,/,A,I15,/,A,A)') &
               & " ! Storing restart file:",&

               & "             at time=",time,&
               & "            step no.=",n_time_step,&
               & "           into file=",rst_file
          call safeClose(nLF)
       end if
    end if
#endif


    ! ===================================================
    !      GATHERING for output
    ! ===================================================
    ! We have all fields in LMloc space. Thus we gather the whole fields on rank 0.

    l_PVout=l_PV .AND. l_log

    !IF (l_log.OR.l_frame.OR.l_graph.OR.l_cmb.OR.l_r.OR.l_Bpot.OR.l_Vpot&

#ifdef WITH_HDF5
    IF (l_frame.OR.l_graph.OR.l_r.OR.l_Bpot.OR.l_Vpot.OR.l_Tpot &
         .OR.(l_SRIC.AND.l_stop_time).OR.l_PVout .or.l_rMagSpec) THEN
#else
    IF (l_frame.OR.l_graph.OR.l_r.OR.l_Bpot.OR.l_Vpot&
         & .OR.l_Tpot.OR.l_store.OR.(l_SRIC.AND.l_stop_time).OR.l_PVout&
         & .or.l_rMagSpec) THEN
#endif
#if 0
       WRITE(*,"(13(A,L1))") "l_log=",l_log,&
            & ", l_frame=",l_frame,&
            & ", l_graph=",l_graph,&
            & ", l_cmb=",l_cmb,&
            & ", l_r=",l_r,&
            & ", l_Bpot=",l_Bpot,&
            & ", l_Vpot=",l_Vpot,&
            & ", l_Tpot=",l_Tpot,&
            & ", l_store=",l_store,&
            & ", l_SRIC=",l_SRIC,&
            & ", l_stop_time=",l_stop_time,&
            & ", l_PVout=",l_PVout,&
            & ", l_rMagSpec=",l_rMagSpec
#endif
       PERFON('out_comm')
       CALL gather_all_from_lo_to_rank0(gt_OC,w_LMloc,w)
       call gather_all_from_lo_to_rank0(gt_OC,dw_LMloc,dw)
       call gather_all_from_lo_to_rank0(gt_OC,ddw_LMloc,ddw)

       CALL gather_all_from_lo_to_rank0(gt_OC,p_LMloc,p)

       call gather_all_from_lo_to_rank0(gt_OC,s_LMloc,s)

       call gather_all_from_lo_to_rank0(gt_OC,z_LMloc,z)
       call gather_all_from_lo_to_rank0(gt_OC,dz_LMloc,dz)

       IF (l_mag) THEN
          CALL gather_all_from_lo_to_rank0(gt_OC,b_LMloc,b)
          CALL gather_all_from_lo_to_rank0(gt_OC,db_LMloc,db)
          CALL gather_all_from_lo_to_rank0(gt_OC,ddb_LMloc,ddb)
          
          CALL gather_all_from_lo_to_rank0(gt_OC,aj_LMloc,aj)
          CALL gather_all_from_lo_to_rank0(gt_OC,dj_LMloc,dj)
          CALL gather_all_from_lo_to_rank0(gt_OC,ddj_LMloc,ddj)
       END IF

       IF (l_cond_ic) THEN
          CALL gather_all_from_lo_to_rank0(gt_IC,b_ic_LMloc,b_ic)
          CALL gather_all_from_lo_to_rank0(gt_IC,db_ic_LMloc,db_ic)
          CALL gather_all_from_lo_to_rank0(gt_IC,ddb_ic_LMloc,ddb_ic)
          
          CALL gather_all_from_lo_to_rank0(gt_IC,aj_ic_LMloc,aj_ic)
          CALL gather_all_from_lo_to_rank0(gt_IC,dj_ic_LMloc,dj_ic)
          CALL gather_all_from_lo_to_rank0(gt_IC,ddj_ic_LMloc,ddj_ic)
       END IF

       ! for writing a restart file, we also need the d?dtLast arrays, which first have to
       ! be gathered on rank 0

#ifndef WITH_HDF5
       IF (l_store) THEN
          CALL gather_all_from_lo_to_rank0(gt_OC,dwdtLast_LMloc,dwdtLast)
          CALL gather_all_from_lo_to_rank0(gt_OC,dpdtLast_LMloc,dpdtLast)
          CALL gather_all_from_lo_to_rank0(gt_OC,dsdtLast_LMloc,dsdtLast)
          CALL gather_all_from_lo_to_rank0(gt_OC,dzdtLast_lo,dzdtLast)
          
          IF (l_mag) THEN
             CALL gather_all_from_lo_to_rank0(gt_OC,dbdtLast_LMloc,dbdtLast)
             CALL gather_all_from_lo_to_rank0(gt_OC,djdtLast_LMloc,djdtLast)
          END IF

          IF (l_cond_ic) THEN
             CALL gather_all_from_lo_to_rank0(gt_IC,dbdt_icLast_LMloc,dbdt_icLast)
             CALL gather_all_from_lo_to_rank0(gt_IC,djdt_icLast_LMloc,djdt_icLast)
          END IF
       END IF
#endif

       PERFOFF

       IF (DEBUG_OUTPUT) THEN
          IF (rank == 0) THEN
             WRITE(*,"(A,8ES22.14)") "output: w,z,p,s = ",SUM( w ),SUM( z ),SUM( p ),SUM( s )
          END IF
       END IF
    ELSE if (l_cmb) then
       ! just gather B_cmb on rank 0 for the B_cmb output
       IF (l_mag) THEN
          !WRITE(*,"(A)") "Gathering only b to rank 0."
          CALL gather_all_from_lo_to_rank0(gt_OC,b_LMloc,b)
       END IF
    END IF

    IF (l_frame) THEN
       ! The frames array for the movies is distributed over the ranks
       ! and has to be gathered on rank 0 for output.

       ! Each movie uses some consecutive frames in the frames array. They
       ! start at n_movie_field_start(1,n_movie) 
       ! up to    n_movie_field_stop(1+n_fields_oc+n_fields,n_movie) (n_fields_ic>0
       ! or       n_movie_field_stop(1+n_fields,n_movie)             (n_fields_ic=0)

       call movie_gather_frames_to_rank0
    END IF

    ! =======================================================================
    ! ======= compute output on rank 0 ==============
    ! =======================================================================
    IF (rank == 0) THEN
       PERFON('out_out')

       !----- Plot out inner core magnetic field, outer core
       !      field has been written in radialLoop !
       IF ( l_graph .AND. l_mag .AND. n_r_ic_max > 0 )                &
            &     CALL graphOut_IC(ngform,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,b)

       !--- Write spectra output that has partially been calculated in LMLoop
       IF ( l_rMagSpec .AND. n_time_step > 1 ) THEN
          IF ( l_frame ) THEN
             CALL rBrSpec(time,b, b_ic ,'rBrSpecMov',.TRUE.,st_map)
             CALL rBpSpec(time,aj,aj_ic,'rBpSpecMov',.TRUE.,st_map)
          END IF
          IF ( l_log ) THEN
             CALL rBrSpec(time,b, b_ic ,'rBrSpec',.TRUE.,st_map)
             CALL rBpSpec(time,aj,aj_ic,'rBpSpec',.TRUE.,st_map)
          END IF
       END IF

       !--- Movie output and various supplementary things:
       IF ( l_frame ) THEN
          PERFON('out_fram')
          IF ( l_movie_ic .AND. l_store_frame ) THEN
             !WRITE(*,"(A)") "Calling store_movie_frame_IC from output."
             CALL store_movie_frame_IC(b,b_ic,db_ic,ddb_ic,aj_ic,dj_ic)
          END IF

          n_frame=n_frame+1
          CALL logWrite(' ')
          WRITE(message,'(1p,A,I8,A,D16.6,I8)') &
               & " ! WRITING MOVIE FRAME NO ",n_frame,&
               & "       at time/step",timeScaled,n_time_step
          CALL logWrite(message)

          !--- Storing the movie frame:
          CALL write_movie_frame(n_frame,timeScaled,                &
               &                 b,db,aj,dj,b_ic,db_ic,aj_ic,dj_ic, &
               &                 omega_ic,omega_ma)

          IF ( l_cmb_field ) THEN
             CALL write_Bcmb(timeScaled,b(1,n_r_cmb),1,lm_max,l_max, &
                  &          l_max_cmb,minc,lm2,n_cmb_setsMov,     &
                  &          cmbMov_file,n_cmbMov_file)
          END IF
          PERFOFF
       END IF ! write movie frame ?

       !--- Store poloidal magnetic coeffs at cmb
       IF ( l_cmb ) THEN
          PERFON('out_cmb')
          CALL write_Bcmb(timeScaled,b(1,n_r_cmb),1,lm_max,l_max,           &
               &          l_max_cmb,minc,lm2,n_cmb_sets,           &
               &          cmb_file,n_cmb_file)
          
          !--- Store SV of poloidal magnetic coeffs at cmb
          IF ( l_dt_cmb_field .AND. ASSOCIATED(dbdt_at_CMB) ) THEN
             !nR=8! at CMB dbdt=induction=0, only diffusion !
             DO lm=1,lm_max
                dbdtCMB(lm)= dbdt_at_CMB(lm)/(dLh(lm)*or2(n_r_cmb))       &
                     &       + opm*hdif_B(lm) * &
                     &       ( ddb(lm,n_r_cmb) - dLh(lm)*or2(n_r_cmb)*b(lm,n_r_cmb) )
             END DO
             CALL write_Bcmb(timeScaled,dbdtCMB,1,lm_max,l_max,             &
                  &          l_max_cmb,minc,lm2,n_dt_cmb_sets,             &
                  &          dt_cmb_file,n_dt_cmb_file)
          END IF
          PERFOFF
       END IF

       !--- Store potential coeffs for velocity fields and magnetic fields
       IF ( l_r ) THEN
          PERFON('out_r')
          DO n=1,n_coeff_r_max
             nR=n_coeff_r(n)
             CALL write_coeff_r(timeScaled,                            &
                  &             w(1,nR),dw(1,nR),ddw(1,nR),z(1,nR),r(nR),&
                  &             1,lm_max,l_max,l_max_r,minc,lm2,n_v_r_sets(n),&
                  &             v_r_file(n),n_v_r_file(n),1)
             IF ( l_mag ) &
                CALL write_coeff_r(timeScaled,                         &
                     &             b(1,nR),db(1,nR),ddb(1,nR),aj(1,nR),r(nR),&
                     &             1,lm_max,l_max,l_max_r,minc,lm2,n_b_r_sets(n),&
                     &             b_r_file(n),n_b_r_file(n),2)
             IF ( l_r_fieldT ) &
                CALL write_coeff_r(timeScaled, &
                                   s(1,nR),db(1,nR),ddb(1,nR),aj(1,nR),r(nR), &
                                   1,lm_max,l_max,l_max_r,minc,lm2,n_T_r_sets(n), &
                                   T_r_file(n),n_t_r_file(n),3)
          END DO
          PERFOFF
       END IF

       IF ( l_log ) THEN
          !--- Energies and rotation info and a lot of other stuff 
          !    performed for l_log=.TRUE.

          !----- Getting the property parameters:
          Re     = SQRT(2.D0*e_kin/vol_oc)/SQRT(mass)
          ReConv = SQRT(2.D0*e_kin_nas/vol_oc)/SQRT(mass)

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

          IF ( dlV /= 0d0 ) THEN
             Rol=Ro/dlV   ! See Christensen&Aubert 2006, eqn.(27)
          ELSE
             Rol=Ro
          END IF
          IF ( dlVc /= 0d0 ) THEN
             RolC=RoConv/dlVc
          ELSE
             RolC=RoConv
          END IF
          !WRITE(*,"(A,3ES20.12)") "dlVc,RoConv,RolC = ",dlVc,RoConv,RolC

          IF ( prmag.NE.0 .AND. nVarCond > 0 ) THEN
             Rm=0.d0
             Rm=rInt_R(RmR,n_r_max,n_r_max,drx, &
                  &    i_costf_init,d_costf_init)
             Rm=Rm*3/(r_cmb**3-r_icb**3)
          elseif ( prmag.NE.0 ) THEN
             Rm=Re*prmag
          ELSE
             Rm=Re
          END IF
          !El   =2.D0*e_mag/vol_oc/LFfac
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
                   lbDiss=SQRT((e_mag+e_mag_ic)/ABS(ohmDiss)) ! Ohmic diffusion 
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

       END IF ! l_log

       IF ( l_Bpot )                                       &
            &     CALL storePot(time,b,aj,b_ic,aj_ic,      &
            &        nBpotSets,'Bpot.',omega_ma,omega_ic)
       IF ( l_Vpot )                                       &
            &     CALL storePot(time,w,z,b_ic,aj_ic,       &
            &        nVpotSets,'Vpot.',omega_ma,omega_ic)
       IF ( l_Tpot )                                       &
            &     CALL storePot(time,s,z,b_ic,aj_ic,       &
            &        nVpotSets,'Tpot.',omega_ma,omega_ic)
       
       !----- Store current solution
       !      Note: unless l_new_rst_file=.TRUE. .and. .not.l_stop_time
       !            this is written into rst_end.TAG

#ifndef WITH_HDF5
       IF ( l_store ) THEN
!#ifdef WITH_HDF5
!          if ( l_stop_time .or. .not.l_new_rst_file ) then
!             rst_file='ser_h5_rst_end.'//tag_wo_rank
!          else if ( l_new_rst_file ) then
!             CALL dble2str(time,string)
!             rst_file='h5_rst_t='//trim(string)//'.'//tag_wo_rank
!          end if
!          CALL storeHdf5_serial(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic, &
!                                  dwdtLast,dzdtLast,dpdtLast,          &
!                                  dsdtLast,dbdtLast,djdtLast,          &
!                                  dbdt_icLast,djdt_icLast)
!#else
          PERFON('out_rst')
          IF ( l_stop_time .OR. .NOT.l_new_rst_file ) THEN
             rst_file="rst_end."//tag_wo_rank
          ELSE IF ( l_new_rst_file ) THEN
             CALL dble2str(time,string)
             rst_file='rst_t='//TRIM(string)//'.'//tag_wo_rank
          END IF

          OPEN(n_rst_file, FILE=rst_file, STATUS='UNKNOWN', FORM='UNFORMATTED')
          CALL store(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic,dwdtLast,dzdtLast, &
                     dpdtLast,dsdtLast,dbdtLast,djdtLast,dbdt_icLast,djdt_icLast)
          CLOSE(n_rst_file)
!#endif

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
          PERFOFF
       END IF
#endif
       
       IF ( l_SRIC .AND. l_stop_time ) CALL outOmega(z,omega_ic)
       
       !----- Output of axisymm. rotation rate for potential vorticity analysis:
       !  NOTE: For l_stop_time=.TRUE. outPV transforms the fields without 
       !        transforming them back. This must thus be the very last 
       !        thing done with them. 
       IF ( l_PVout ) CALL outPV(time,l_stop_time,nPVsets,             &
            &                     w,dw,ddw,z,dz,omega_ic,omega_ma)
       
       PERFOFF
    END IF

    IF ( l_log ) THEN
       timePassedLog=0.0D0
    END IF

    IF ( lRmsCalc ) THEN
       CALL zeroRms
    END IF
    
    RETURN 
  END SUBROUTINE output


  !-----------------------------------------------------------------------
END MODULE output_mod
