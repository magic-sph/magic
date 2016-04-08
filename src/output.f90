#include "perflib_preproc.cpp"
module output_mod

   use precision_mod
   use truncation, only: n_r_max, n_r_ic_max, minc, l_max, l_maxMag, &
       &                 n_r_maxMag, lm_max
   use parallel_mod, only: rank
   use radial_functions, only: or1, or2, r, drx, chebt_oc, r_cmb, r_icb
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag, n_r_cmb
   use physical_parameters, only: opm,ek,ktopv,prmag,nVarCond,LFfac
   use num_param, only: tScale
   use blocking, only: st_map, lm2, lo_map
   use horizontal_data, only: dLh,hdif_B,dPl0Eq
   use logic, only: l_average, l_mag, l_power, l_anel, l_mag_LF, lVerbose, &
       &            l_dtB, l_RMS, l_r_field, l_r_fieldT, l_PV, l_SRIC,     &
       &            l_cond_ic,l_rMagSpec, l_movie_ic, l_store_frame,       &
       &            l_cmb_field, l_dt_cmb_field, l_save_out, l_non_rot,    &
       &            l_perpPar, l_energy_modes, l_heat, l_hel, l_par
   use fields, only: omega_ic, omega_ma, b, db, aj, dj, b_ic,              &
       &             db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic, w, z,            &
       &             s, p, w_LMloc, dw_LMloc, ddw_LMloc, p_LMloc,          &
       &             s_LMloc, ds_LMloc, z_LMloc, dz_LMloc, b_LMloc,        &
       &             db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc, ddj_LMloc,   &
       &             b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc,   &
       &             dj_ic_LMloc, ddj_ic_LMloc, dp_LMloc
   use fieldsLast, only: dwdtLast, dzdtLast, dpdtLast, dsdtLast, dbdtLast,  &
       &                 djdtLast, dbdt_icLast, djdt_icLast, dwdtLast_LMloc,&
       &                 dzdtLast_lo, dpdtLast_LMloc, dsdtLast_LMloc,       &
       &                 dbdtLast_LMloc, djdtLast_LMloc, dbdt_icLast_LMloc, &
       &                 djdt_icLast_LMloc
   use kinetic_energy, only: get_e_kin, get_u_square
   use magnetic_energy, only: get_e_mag
   use fields_average_mod, only: fields_average
   use spectra, only: spectrum_average, spectrum, spectrum_temp, &
       &              spectrum_temp_average, get_amplitude
   use outTO_mod, only: outTO
   use outPV3, only: outPV
   use output_data, only: tag, l_max_cmb,                           &
       &                  cmbMov_file, n_cmbMov_file, cmb_file,     &
       &                  n_cmb_file, dt_cmb_file, n_dt_cmb_file,   & 
       &                  n_coeff_r, l_max_r, n_v_r_file,           &
       &                  n_b_r_file, n_t_r_file, v_r_file,         &
       &                  t_r_file, b_r_file, n_r_array, n_r_step,  &
       &                  par_file, n_par_file, nLF, log_file,      &
       &                  n_coeff_r_max, rst_file, n_rst_file
   use constants, only: vol_oc, vol_ic, mass, surf_cmb, two
   use outMisc_mod, only: outHelicity, outHeat
   use Egeos_mod, only: getEgeos
   use outRot, only: write_rot
   use charmanip, only: dble2str
   use omega, only: outOmega
   use integration, only: rInt_R
   use outPar_mod, only: outPar, outPerpPar
   use graphOut_mod, only: graphOut_IC
   use power, only: get_power
   use LMLoop_data, only: lm_per_rank, lm_on_last_rank, llm, ulm, llmMag, &
       &                  ulmMag
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC
   use out_coeff, only: write_Bcmb, write_coeff_r
   use getDlm_mod, only: getDlm
   use movie_data, only: movie_gather_frames_to_rank0
   use dtB_mod, only: get_dtBLMfinish
   use out_movie, only: write_movie_frame
   use out_movie_IC, only: store_movie_frame_IC
   use RMS, only: zeroRms, dtVrms, dtBrms
   use store_pot_mod, only: storePot
   use useful, only: safeOpen, safeClose, logWrite
   use radial_spectra  ! rBrSpec, rBpSpec
   use storeCheckPoints

   implicit none
 
   private
 
   integer :: nBpotSets, nVpotSets, nTpotSets
   !-- Counter for output files/sets:
   integer :: n_dt_cmb_sets, n_cmb_setsMov
   integer, allocatable :: n_v_r_sets(:), n_b_r_sets(:), n_T_r_sets(:)
   integer :: n_spec,nPVsets
 
   integer :: nTOsets,nTOmovSets,nTOrmsSets
 
   !--- For averaging:
   real(cp) :: timePassedLog, timeNormLog
   integer :: nLogs  
 
   real(cp), save :: dlBMean,dmBMean
   real(cp), save :: lvDissMean,lbDissMean
   real(cp), save :: RmMean,ElMean,ElCmbMean,RolMean,GeosMean
   real(cp), save :: DipMean,DipCMBMean
   real(cp), save :: dlVMean,dlVcMean,dmVMean,dpVMean,dzVMean
 
   real(cp) :: eTot,eTotOld,dtEint
   real(cp) :: e_kin_pMean, e_kin_tMean
   real(cp) :: e_mag_pMean, e_mag_tMean
   integer :: n_e_sets

   real(cp) :: timePassedRMS, timeNormRMS
   integer :: nRMS_sets

   public :: output, initialize_output

contains

   subroutine initialize_output

      integer :: n

      if ( l_r_field .or. l_r_fieldT ) then
         allocate ( n_coeff_r(n_coeff_r_max))
         allocate ( n_v_r_file(n_coeff_r_max), v_r_file(n_coeff_r_max) )
         allocate ( n_v_r_sets(n_coeff_r_max) ) 
         n_v_r_sets=0

         if ( l_mag ) then
            allocate ( n_b_r_file(n_coeff_r_max), b_r_file(n_coeff_r_max) )
            allocate ( n_b_r_sets(n_coeff_r_max) ) 
            n_b_r_sets=0
         end if

         if ( l_r_fieldT ) then
            allocate ( n_t_r_file(n_coeff_r_max), t_r_file(n_coeff_r_max) )
            allocate ( n_t_r_sets(n_coeff_r_max) ) 
            n_T_r_sets=0
         end if

         if ( count(n_r_array>0)> 0 ) then
            n_coeff_r=n_r_array(1:n_coeff_r_max)
         else
            n_r_step=max(n_r_step,1)
            do n=1,n_coeff_r_max
               n_coeff_r(n)=n*n_r_step  ! used every n_r_step point !
            end do
         end if

      end if

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
      
      timeNormLog  =0.0_cp
      timePassedLog=0.0_cp
      RmMean       =0.0_cp
      ElMean       =0.0_cp
      ElCmbMean    =0.0_cp
      RolMean      =0.0_cp
      GeosMean     =0.0_cp
      DipMean      =0.0_cp
      DipCMBMean   =0.0_cp
      e_kin_pMean  =0.0_cp
      e_kin_tMean  =0.0_cp
      e_mag_pMean  =0.0_cp
      e_mag_tMean  =0.0_cp
      dlVMean      =0.0_cp
      dlVcMean     =0.0_cp
      dmVMean      =0.0_cp
      dpVMean      =0.0_cp
      dzVMean      =0.0_cp
      dlBMean      =0.0_cp
      dmBMean      =0.0_cp
      lvDissmean   =0.0_cp
      lbDissmean   =0.0_cp

   end subroutine initialize_output
!----------------------------------------------------------------------------
   subroutine output(time,dt,dtNew,n_time_step,l_stop_time,               &
        &            l_Bpot,l_Vpot,l_Tpot,l_log,l_graph,lRmsCalc,         &
        &            l_store,l_new_rst_file,                              &
        &            l_spectrum,lTOCalc,lTOframe,lTOZwrite,               &
        &            l_frame,n_frame,l_cmb,n_cmb_sets,l_r,                &
        &            lorentz_torque_ic,lorentz_torque_ma,dbdt_CMB_LMloc,  &
        &            HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,      &
        &            gradsLMr,fconvLMr,fkinLMr,fviscLMr,fpoynLMr,         &
        &            fresLMr,EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)
      !
      !  This subroutine controls most of the output.                     
      !
  
      !--- Input of variables
      real(cp),    intent(in) :: time,dt,dtNew
      integer,     intent(in) :: n_time_step
      logical,     intent(in) :: l_stop_time
      logical,     intent(in) :: l_Bpot,l_Vpot,l_Tpot
      logical,     intent(in) :: l_log, l_graph, lRmsCalc, l_store
      logical,     intent(in) :: l_new_rst_file, l_spectrum
      logical,     intent(in) :: lTOCalc,lTOframe
      logical,     intent(in) :: l_frame, l_cmb, l_r
      logical,     intent(inout) :: lTOZwrite
      integer,     intent(inout) :: n_frame
      integer,     intent(inout) :: n_cmb_sets
  
      !--- Input of Lorentz torques and dbdt calculated in radialLoopG
      !    Parallelization note: Only the contribution at the CMB must be 
      !    collected and is (likely) stored on the processor (#0) that performs 
      !    this routine anyway.
      real(cp),    intent(in) :: lorentz_torque_ma,lorentz_torque_ic
  
      !--- Input of scales fields via common block in fields.f90:
      !    Parallelization note: these fields are LM-distributed.
      !    The input fields HelLMr,Hel2LMr,TstrRLM,TadvRLM, and TomeRLM
      !    are R-distributed. More R-distributed fields are hidden 
      !    in TO.f90, RMS.f90, and dtB.f90. 
      !    input fields are R-distributed. This has to be taken into
      !    account when collecting the information from the different
      !    processors!
      !    All the fields contained in fields.fi90 are needed on
      !    the processor performing this routine:
      !          w,dw,ddw,z,dz,s,p,b,db,aj,dj,ddj,
      !          b_ic,db_ic,ddb_ic,aj_ic,dj_ic,omega_ic,omega_ma
      !    omega_ic and omega_ma are likely located on processor #0 
      !    which deals with (l=1,m=0) in s_updateZ.f
      !    Note that many of these only have to be collected when
      !    certain output is required. This is controlled by the 
      !    input logicals.
  
      !--- Input help arrays for magnetic field stretching and advection and
      !    for calculating axisymmetric helicity.
      !    Parallelization note: These fields are R-distribute on input 
      !    and must also be collected on the processor performing this routine.
      real(cp),    intent(in) :: HelLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: Helna2LMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: uhLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: gradsLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: duhLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: fconvLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: fkinLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: fviscLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp),    intent(in) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp),    intent(in) :: EperpLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: EparLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: EperpaxiLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: EparaxiLMr(l_max+1,nRstart:nRstop)

      complex(cp), intent(in) :: dbdt_CMB_LMloc(llmMag:ulmMag)
  
      !--- Local stuff:
      !--- Energies:
      real(cp) :: ekinR(n_r_max)     ! kinetic energy w radius
      real(cp) :: e_mag,e_mag_ic,e_mag_cmb       
      real(cp) :: e_mag_p,e_mag_t      
      real(cp) :: e_mag_p_as,e_mag_t_as   
      real(cp) :: e_mag_p_ic,e_mag_t_ic   
      real(cp) :: e_mag_p_as_ic,e_mag_t_as_ic
      real(cp) :: e_mag_os,e_mag_as_os    
      real(cp) :: e_kin,e_kin_p,e_kin_t  
      real(cp) :: e_kin_p_as,e_kin_t_as 
      real(cp) :: eKinIC,eKinMA        
      real(cp) :: dtE
  
      integer :: nR,lm,n,m
  
      !--- For TO:
      character(len=64) :: TOfileNhs,TOfileShs,movFile
      character(len=66) :: tayFile
      logical :: lTOrms    
      integer :: nF1,nF2
  
      !--- Property parameters:
      complex(cp) :: dbdtCMB(llmMag:ulmMag)        ! SV at CMB !
      real(cp) :: dlBR(n_r_max),dlBRc(n_r_max),dlVR(n_r_max),dlVRc(n_r_max)
      real(cp) :: RolRu2(n_r_max),dlVRu2(n_r_max),dlVRu2c(n_r_max)
      real(cp) :: RmR(n_r_max)
      real(cp) :: Re,Ro,Rm,El,ElCmb,Rol,Geos,Dip,DipCMB
      real(cp) :: ReConv,RoConv,e_kin_nas,RolC
      real(cp) :: elsAnel
      real(cp) :: dlB,dlBc,dmB
      real(cp) :: dlV,dlVc,dmV,dpV,dzV
      real(cp) :: visDiss,ohmDiss,lvDiss,lbDiss
      integer :: l,lm0
      real(cp) :: ReEquat
  
      logical :: l_PVout
  
      real(cp) :: timeScaled
  
      character(len=76) :: filename
      character(len=96) :: message
  
      character(len=20) :: string
      logical :: DEBUG_OUTPUT=.false.
  
      timeScaled=tScale*time
      timePassedLog=timePassedLog+dt
  
      ! We start with the computation of the energies
      ! in parallel.
      if ( l_log ) then
         nLogs=nLogs+1
         timeNormLog=timeNormLog+timePassedLog
  
         !----- Write torques and rotation rates:
         PERFON('out_rot')
         call write_rot( time,dt,eKinIC,eKinMA,w_LMloc,z_LMloc,dz_LMloc,b_LMloc,  &
              &          omega_ic,omega_ma,lorentz_torque_ic,lorentz_torque_ma)
         PERFOFF
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  write_rot  on rank ",rank
  
         PERFON('out_ekin')
         n_e_sets=n_e_sets+1
         call get_e_kin(time,.true.,l_stop_time,n_e_sets,w_LMloc,    &
              &         dw_LMloc,z_LMloc,e_kin_p,e_kin_t,e_kin_p_as, &
              &         e_kin_t_as,ekinR)
         e_kin=e_kin_p+e_kin_t
         !write(*,"(A,3(I4,F20.17))") "e_kin, e_kin_p_as,e_kin_t_as = ",&
         !     &EXPONENT(e_kin),FRACTION(e_kin),&
         !     &EXPONENT(e_kin_p_as),FRACTION(e_kin_p_as),&
         !     &EXPONENT(e_kin_t_as),FRACTION(e_kin_t_as)
         e_kin_nas=e_kin-e_kin_p_as-e_kin_t_as
         if ( DEBUG_OUTPUT ) write(*,"(A,I6)") "Written  e_kin  on rank ",rank
  
         call get_e_mag(time,.true.,l_stop_time,n_e_sets,b_LMloc,db_LMloc, &
              &         aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc,       &
              &         e_mag_p,e_mag_t,e_mag_p_as,e_mag_t_as,e_mag_p_ic,  &
              &         e_mag_t_ic,e_mag_p_as_ic,e_mag_t_as_ic,            &
              &         e_mag_os,e_mag_as_os,e_mag_cmb,Dip,DipCMB,elsAnel )
         e_mag   =e_mag_p+e_mag_t
         e_mag_ic=e_mag_p_ic+e_mag_t_ic
         PERFOFF
         if ( DEBUG_OUTPUT ) write(*,"(A,I6)") "Written  e_mag  on rank ",rank

         !----- Calculate distribution of energies on all m's
         if ( l_energy_modes ) then  
            PERFON('out_amplitude') 
            call get_amplitude(time,w_LMloc,dw_LMloc,z_LMloc,b_LMloc,&
                 &             db_LMloc,aj_LMloc)
            PERFOFF
            if ( DEBUG_OUTPUT ) &
               & write(*,"(A,I6)") "Written  amplitude  on rank ",rank
         endif
  
         if ( l_average ) then
            PERFON('out_aver')
            call spectrum_average(nLogs,l_stop_time,timePassedLog,  &
                 &                timeNormLog,w_LMloc,z_LMloc,      &
                 &                dw_LMloc,'V')
            call spectrum_temp_average(nLogs,l_stop_time,timePassedLog, &
                 &                     timeNormLog,s_LMloc,ds_LMloc)
  
            if ( l_mag ) then
               call spectrum_average(nLogs,l_stop_time,timePassedLog, &
                    &                timeNormLog,b_LMloc,aj_LMloc,db_LMloc,'B')
            end if
  
            call fields_average(nLogs,l_stop_time,timePassedLog,timeNormLog, &
                 &              omega_ic,omega_ma,w_LMloc,z_LMloc,p_LMloc,   &
                 &              s_LMloc,b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc)
            PERFOFF
            if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  averages  on rank ",rank
         end if
  
         if ( l_power ) then
  
            PERFON('out_pwr')
            if ( rank == 0 ) then
               if ( nLogs > 1 ) then
                  filename='dtE.'//tag
                  open(99,file=filename, status='unknown', position='append')
                  eTotOld=eTot
                  eTot   =e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
                  dtE    =(eTot-eTotOld)/timePassedLog
                  dtEint =dtEint+timePassedLog*(eTot-eTotOld)
                  write(99,'(ES20.10,3ES16.6)') time,dtE,                  &
                       &                    dtEint/timeNormLog,dtE/eTot
                  close(99)
               else
                  eTot   =e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
                  dtEint=0.0_cp
               end if
               !write(*,"(A,7ES22.14)") "eTot = ",eTot,e_kin,e_mag,e_mag_ic, &
               !     &                            e_mag_os,eKinIC,eKinMA
            end if
            call get_power( time,timePassedLog,timeNormLog,l_stop_time,      &
                 &          omega_ic,omega_ma,lorentz_torque_ic,             &
                 &          lorentz_torque_ma,w_LMloc,ddw_LMloc,z_LMloc,     &
                 &          dz_LMloc,s_LMloc,b_LMloc,ddb_LMloc,aj_LMloc,     &
                 &          dj_LMloc,db_ic_LMloc,ddb_ic_LMloc,aj_ic_LMloc,   &
                 &          dj_ic_LMloc,visDiss,ohmDiss)
            PERFOFF
            if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  power  on rank ",rank
         end if
  
         !----- If anelastic additional u**2 outputs
         if ( l_anel ) then
            call get_u_square(time,w_LMloc,dw_LMloc,z_LMloc,RolRu2,dlVRu2,dlVRu2c)
            if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  u_square  on rank ",rank
         else
            RolRu2  = 0.0_cp
            dlVRu2  = 0.0_cp
            dlVRu2c = 0.0_cp
         end if
  
         if ( l_perpPar ) then
            call outPerpPar(time,timePassedLog,timeNormLog,l_stop_time, &
                 &          EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)
         end if
         !----- Radial properties
         !write(*,"(A,4ES20.12)") "before getDlm, w(n_r_icb,n_r_cmb): ",&
         !     & w_LMloc(n_r_icb),w_LMloc(n_r_cmb)
         !write(*,"(A,4ES20.12)") "before getDlm, dw(n_r_icb,n_r_cmb): ",&
         !     & dw_LMloc(n_r_icb),dw_LMloc(n_r_cmb)
         !write(*,"(A,4ES20.12)") "before getDlm, z(n_r_icb,n_r_cmb): ",&
         !     & z_LMloc(n_r_icb),z_LMloc(n_r_cmb)
         call getDlm(w_LMloc,dw_LMloc,z_LMloc,dlV,dlVR,dmV,dlVc,dlVRc,'V')
         !write(*,"(A,ES20.12)") "dlVr,dlVrc(n_r_icb) = ",dlVr(n_r_icb),dlVrc(n_r_icb)
         !write(*,"(A,ES20.12)") "dlVr,dlVrc(n_r_cmb) = ",dlVr(n_r_cmb),dlVrc(n_r_cmb)
         call outPar(timePassedLog,timeNormLog,nLogs,l_stop_time,        &
              &      ekinR,RolRu2,dlVR,dlVRc,dlVRu2,dlVRu2c,             &
              &      uhLMr,duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,    &
              &      fpoynLMr,fresLMr,RmR)
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  outPar  on rank ",rank
  
         if ( l_heat ) then
            call outHeat(timeScaled,timePassedLog,timeNormLog,l_stop_time, &
                 &       s_LMloc,ds_LMloc,p_LMloc,dp_LMloc)
         end if

         if ( l_hel ) then
            call outHelicity(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr)
         end if

         if ( l_par ) then
            call getEgeos(timeScaled,nLogs,w_LMloc,dw_LMloc,ddw_LMloc,    &
                 &        z_LMloc,dz_LMloc,Geos,dpV,dzV)
         end if

         if ( l_mag .or. l_mag_LF ) then 
            call getDlm(b_LMloc,db_LMloc,aj_LMloc,dlB,dlBR,dmB,dlBc,dlBRc,'B')
         else
            dlB=0.0_cp
            dmB=0.0_cp
         end if
      end if
  
      if ( l_spectrum ) then
         n_spec=n_spec+1
         call spectrum(time,n_spec,w_LMloc,dw_LMloc,z_LMloc,b_LMloc,  &
              &        db_LMloc,aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc)
         call spectrum_temp(time,n_spec,s_LMloc,ds_LMloc)
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  spectrum  on rank ",rank
      end if
  
      if ( lTOCalc ) then
         !------ Output for every log time step:
         if ( lVerbose ) write(*,*) '! Calling outTO !'
         TOfileNhs='TOnhs.'//tag
         TOfileShs='TOshs.'//tag
         movFile  ='TO_mov.'//tag
         tayFile  ='TaySphere4.'//tag
         nF1      =93
         nF2      =94
         lTOrms   =.true.
         if ( .not. l_log ) then
            call get_e_kin(time,.false.,l_stop_time,0,w_LMloc,dw_LMloc,  &
                 &         z_LMloc,e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,&
                 &         ekinR)
            e_kin=e_kin_p+e_kin_t
         end if
         call outTO(time,n_time_step,e_kin,e_kin_t_as,                      &
              &     nF1,nF2,TOfileNhs,TOfileShs,movFile,tayFile,            &
              &     nTOsets,nTOmovSets,nTOrmsSets,lTOframe,lTOrms,lTOZwrite,&
              &     z_LMloc,omega_ic,omega_ma)
         !------ Note: time averaging, time differencing done by IDL routine!
  
         if ( lVerbose ) write(*,*) '! outTO finished !'
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  TO  on rank ",rank
      end if
  
      !--- Get radial derivatives and add dt dtB terms:
      if ( l_dtB ) then
         call get_dtBLMfinish(time,n_time_step,omega_ic,b_LMloc,ddb_LMloc, &
              &               aj_LMloc,dj_LMloc,ddj_LMloc,b_ic_LMloc,      &
              &               db_ic_LMloc,ddb_ic_LMloc,aj_ic_LMloc,        &
              &               dj_ic_LMloc,ddj_ic_LMloc)
      end if
  
  
      if ( l_RMS ) then
         if ( n_time_step == 1 ) then
            nRMS_sets    =0
            timeNormRMS  =0.0_cp
            timePassedRMS=0.0_cp
            call zeroRms
         end if
         timePassedRMS=timePassedRMS+dt
         if ( lRmsCalc ) then
            if ( lVerbose ) write(*,*) '! Writing RMS output !'
            timeNormRMS=timeNormRMS+timePassedRMS
            call dtVrms(time,nRMS_sets,timePassedRMS,timeNormRMS)
            if ( l_mag ) call dtBrms(time)
            timePassedRMS=0.0_cp
         end if
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  dtV/Brms  on rank ",rank
      end if
  
      !--- Store poloidal magnetic coeffs at cmb
      if ( l_cmb ) then
         PERFON('out_cmb')
         call write_Bcmb(timeScaled,b_LMloc(:,n_r_cmb),l_max_cmb,n_cmb_sets,   &
              &          cmb_file,n_cmb_file)
         
         !--- Store SV of poloidal magnetic coeffs at cmb
         if ( l_dt_cmb_field ) then
            !nR=8! at CMB dbdt=induction=0, only diffusion !
            !do lm=max(2,llm),ulm
            do lm=max(2,llm),ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               dbdtCMB(lm)= dbdt_CMB_LMloc(lm)/                             &
                    &    (dLh(st_map%lm2(l,m))*or2(n_r_cmb))                       &
                    &    + opm*hdif_B(st_map%lm2(l,m)) * ( ddb_LMloc(lm,n_r_cmb) - &
                    &      dLh(st_map%lm2(l,m))*or2(n_r_cmb)*b_LMloc(lm,n_r_cmb) )
            end do

            call write_Bcmb(timeScaled,dbdtCMB(:),l_max_cmb,n_dt_cmb_sets,  &
                 &          dt_cmb_file,n_dt_cmb_file)
         end if
         PERFOFF
      end if

      if ( l_frame .and. l_cmb_field ) then
         call write_Bcmb(timeScaled,b_LMloc(:,n_r_cmb),l_max_cmb,   &
              &          n_cmb_setsMov,cmbMov_file,n_cmbMov_file)
      end if
  
      !--- Store potential coeffs for velocity fields and magnetic fields
      if ( l_r ) then
         PERFON('out_r')
         do n=1,n_coeff_r_max
            nR=n_coeff_r(n)
            call write_coeff_r(timeScaled,w_LMloc(:,nR),dw_LMloc(:,nR),  &
                 &             ddw_LMloc(:,nR),z_LMloc(:,nR),r(nR),      &
                 &             l_max_r,n_v_r_sets(n),v_r_file(n),        &
                 &             n_v_r_file(n),1)
            if ( l_mag )                                                   &
               call write_coeff_r(timeScaled,b_LMloc(:,nR),db_LMloc(:,nR), &
                    &             ddb_LMloc(:,nR),aj_LMloc(:,nR),r(nR),    &
                    &             l_max_r,n_b_r_sets(n),b_r_file(n),       &
                    &             n_b_r_file(n),2)
            if ( l_r_fieldT )                                              &
               call write_coeff_r(timeScaled,s_LMloc(:,nR),db_LMloc(:,nR), &
                    &             ddb_LMloc(:,nR),aj_LMloc(:,nR),r(nR),    &
                    &             l_max_r,n_T_r_sets(n),T_r_file(n),       &
                    &             n_t_r_file(n),3)
         end do
         PERFOFF
      end if
  
      !
      ! Parallel writing of the restart file (possible only when HDF5 is used)
      !
#ifdef WITH_HDF5
      if ( l_store ) then
  
         if ( l_stop_time .or. .not.l_new_rst_file ) then
            rst_file='h5_rst_end.'//tag
         else if ( l_new_rst_file ) then
            call dble2str(time,string)
            rst_file='h5_rst_t='//trim(string)//'.'//tag
         end if
         call storeHdf5_parallel(time,dt,dtNew,w_LMloc,z_LMloc,p_LMloc,s_LMloc, &
                                 b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc,       &
                                 dwdtLast_LMloc,dzdtLast_lo,dpdtLast_LMloc,     &
                                 dsdtLast_LMloc,dbdtLast_LMloc,djdtLast_LMloc,  &
                                 dbdt_icLast_LMloc,djdt_icLast_LMloc)
  
         if ( rank == 0 ) then
            write(*,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)')&
                 & " ! Storing restart file:",           &
                 & "             at time=",time,         &
                 & "            step no.=",n_time_step,  &
                 & "           into file=",rst_file
            call safeOpen(nLF,log_file)
            
            write(nLF,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)') &
                 & " ! Storing restart file:",              &
                 & "             at time=",time,            &
                 & "            step no.=",n_time_step,     &
                 & "           into file=",rst_file
            call safeClose(nLF)
         end if
      end if
#endif
  
  
      ! ===================================================
      !      GATHERING for output
      ! ===================================================
      ! We have all fields in LMloc space. Thus we gather the whole fields on rank 0.
  
      l_PVout=l_PV .and. l_log
  
      !if (l_log.or.l_frame.or.l_graph.or.or.l_r.or.l_Bpot.or.l_Vpot&
  
#ifdef WITH_HDF5
      if (l_frame.or.l_graph.or.l_Bpot.or.l_Vpot.or.l_Tpot &
           .or.(l_SRIC.and.l_stop_time).or.l_PVout .or.l_rMagSpec) then
#else
      if (l_frame.or.l_graph.or.l_Bpot.or.l_Vpot                   &
           & .or.l_Tpot.or.l_store.or.(l_SRIC.and.l_stop_time).or.l_PVout &
           & .or.l_rMagSpec) then
#endif
#if 0
         write(*,"(13(A,L1))") "l_log=",l_log,     &
              & ", l_frame=",l_frame,              &
              & ", l_graph=",l_graph,              &
              & ", l_Bpot=",l_Bpot,                &
              & ", l_Vpot=",l_Vpot,                &
              & ", l_Tpot=",l_Tpot,                &
              & ", l_store=",l_store,              &
              & ", l_SRIC=",l_SRIC,                &
              & ", l_stop_time=",l_stop_time,      &
              & ", l_PVout=",l_PVout,              &
              & ", l_rMagSpec=",l_rMagSpec
#endif
         PERFON('out_comm')
         call gather_all_from_lo_to_rank0(gt_OC,w_LMloc,w)
         call gather_all_from_lo_to_rank0(gt_OC,p_LMloc,p)
         call gather_all_from_lo_to_rank0(gt_OC,s_LMloc,s)
         call gather_all_from_lo_to_rank0(gt_OC,z_LMloc,z)
  
         if ( l_mag ) then
            call gather_all_from_lo_to_rank0(gt_OC,b_LMloc,b)
            call gather_all_from_lo_to_rank0(gt_OC,db_LMloc,db)
            call gather_all_from_lo_to_rank0(gt_OC,aj_LMloc,aj)
            call gather_all_from_lo_to_rank0(gt_OC,dj_LMloc,dj)
         end if
  
         if ( l_cond_ic ) then
            call gather_all_from_lo_to_rank0(gt_IC,b_ic_LMloc,b_ic)
            call gather_all_from_lo_to_rank0(gt_IC,db_ic_LMloc,db_ic)
            call gather_all_from_lo_to_rank0(gt_IC,ddb_ic_LMloc,ddb_ic)
            
            call gather_all_from_lo_to_rank0(gt_IC,aj_ic_LMloc,aj_ic)
            call gather_all_from_lo_to_rank0(gt_IC,dj_ic_LMloc,dj_ic)
            call gather_all_from_lo_to_rank0(gt_IC,ddj_ic_LMloc,ddj_ic)
         end if
  
         ! for writing a restart file, we also need the d?dtLast arrays, 
         ! which first have to be gathered on rank 0
  
#ifndef WITH_HDF5
         if (l_store) then
            call gather_all_from_lo_to_rank0(gt_OC,dwdtLast_LMloc,dwdtLast)
            call gather_all_from_lo_to_rank0(gt_OC,dpdtLast_LMloc,dpdtLast)
            call gather_all_from_lo_to_rank0(gt_OC,dsdtLast_LMloc,dsdtLast)
            call gather_all_from_lo_to_rank0(gt_OC,dzdtLast_lo,dzdtLast)
            
            if (l_mag) then
               call gather_all_from_lo_to_rank0(gt_OC,dbdtLast_LMloc,dbdtLast)
               call gather_all_from_lo_to_rank0(gt_OC,djdtLast_LMloc,djdtLast)
            end if
  
            if (l_cond_ic) then
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_icLast_LMloc,dbdt_icLast)
               call gather_all_from_lo_to_rank0(gt_IC,djdt_icLast_LMloc,djdt_icLast)
            end if
         end if
#endif
  
         PERFOFF
  
         if (DEBUG_OUTPUT) then
            if ( rank == 0 ) then
               write(*,"(A,8ES22.14)") "output: w,z,p,s = ",SUM( w ), &
                                          SUM( z ),SUM( p ),SUM( s )
            end if
         end if
      end if
  
      if ( l_frame ) then
         ! The frames array for the movies is distributed over the ranks
         ! and has to be gathered on rank 0 for output.
  
         ! Each movie uses some consecutive frames in the frames array. They
         ! start at n_movie_field_start(1,n_movie) 
         ! up to    n_movie_field_stop(1+n_fields_oc+n_fields,n_movie) (n_fields_ic>0
         ! or       n_movie_field_stop(1+n_fields,n_movie)             (n_fields_ic=0)
  
         call movie_gather_frames_to_rank0
      end if
  
      ! =======================================================================
      ! ======= compute output on rank 0 ==============
      ! =======================================================================
      if ( rank == 0 ) then
         PERFON('out_out')
  
         !----- Plot out inner core magnetic field, outer core
         !      field has been written in radialLoop !
         if ( l_graph .and. l_mag .and. n_r_ic_max > 0 )          &
              &     call graphOut_IC(b_ic,db_ic,ddb_ic,aj_ic,dj_ic,b)
  
         !--- Write spectra output that has partially been calculated in LMLoop
         if ( l_rMagSpec .and. n_time_step > 1 ) then
            if ( l_frame ) then
               call rBrSpec(time,b, b_ic ,'rBrSpecMov',.true.,st_map)
               call rBpSpec(time,aj,aj_ic,'rBpSpecMov',.true.,st_map)
            end if
            if ( l_log ) then
               call rBrSpec(time,b, b_ic ,'rBrSpec',.true.,st_map)
               call rBpSpec(time,aj,aj_ic,'rBpSpec',.true.,st_map)
            end if
         end if
  
         !--- Movie output and various supplementary things:
         if ( l_frame ) then
            PERFON('out_fram')
            if ( l_movie_ic .and. l_store_frame ) then
               !write(*,"(A)") "Calling store_movie_frame_IC from output."
               call store_movie_frame_IC(b,b_ic,db_ic,ddb_ic,aj_ic,dj_ic)
            end if
  
            n_frame=n_frame+1
            call logWrite(' ')
            write(message,'(1p,A,I8,A,ES16.6,I8)')            &
                 & " ! WRITING MOVIE FRAME NO ",n_frame,      &
                 & "       at time/step",timeScaled,n_time_step
            call logWrite(message)
  
            !--- Storing the movie frame:
            call write_movie_frame(n_frame,timeScaled,                &
                 &                 b,db,aj,dj,b_ic,db_ic,aj_ic,dj_ic, &
                 &                 omega_ic,omega_ma)
  
            PERFOFF
         end if ! write movie frame ?
  
         if ( l_log ) then
            !--- Energies and rotation info and a lot of other stuff 
            !    performed for l_log=.true.
  
            !----- Getting the property parameters:
            Re     = sqrt(two*e_kin/vol_oc)/sqrt(mass)
            ReConv = sqrt(two*e_kin_nas/vol_oc)/sqrt(mass)
  
            if ( l_non_rot ) then
               Ro=0.0_cp
               RoConv=0.0_cp
            else
               Ro=Re*ek
               RoConv=ReConv*ek
            end if
  
            !---- Surface zonal velocity at the equator
            if ( ktopv==1 ) then
               ReEquat=0.0_cp
               do l=1,l_max
                  lm0=lm2(l,0)
                  ReEquat=ReEquat-real(z(lm0,n_r_cmb))*dPl0Eq(l+1)*or1(n_r_cmb)
               end do
            else
               ReEquat=0.0_cp
            end if
  
            if ( dlV /= 0.0_cp ) then
               Rol=Ro/dlV   ! See Christensen&Aubert 2006, eqn.(27)
            else
               Rol=Ro
            end if
            if ( dlVc /= 0.0_cp ) then
               RolC=RoConv/dlVc
            else
               RolC=RoConv
            end if
            !write(*,"(A,3ES20.12)") "dlVc,RoConv,RolC = ",dlVc,RoConv,RolC
  
            if ( prmag /= 0 .and. nVarCond > 0 ) then
               Rm=0.0_cp
               Rm=rInt_R(RmR,n_r_max,n_r_max,drx,chebt_oc)
               Rm=Rm*3/(r_cmb**3-r_icb**3)
            elseif ( prmag /= 0 ) then
               Rm=Re*prmag
            else
               Rm=Re
            end if
            !El   =two*e_mag/vol_oc/LFfac
            ! Elsasser number is computed from the averaged profile
            if ( l_mag .or. l_mag_LF ) then
               El   =elsAnel/vol_oc
               ElCmb=two*e_mag_cmb/surf_cmb/LFfac
            else
               El   =0.0_cp
               ElCmb=0.0_cp
            end if
            if ( l_power ) then
               if ( visDiss /= 0.0_cp ) then
                  lvDiss=sqrt(e_kin/abs(visDiss))            ! Viscous diffusion
               else
                  lvDiss=0.0_cp
               end if
               if ( l_mag .or. l_mag_LF ) then
                  if ( ohmDiss /= 0.0_cp ) then
                     lbDiss=SQRT((e_mag+e_mag_ic)/ABS(ohmDiss)) ! Ohmic diffusion 
                  else
                     lbDiss=0.0_cp
                  end if
               else
                  lbDiss=0.0_cp
               end if
            else
               lvDiss=0.0_cp
               lbDiss=0.0_cp
            end if
  
            !----- Ouput into par file:
            if ( l_save_out ) then
               open(n_par_file, file=par_file, status='unknown', position='append')
            end if
            write(n_par_file,'(ES20.12,18ES16.8)')  &
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
            if ( l_save_out ) close(n_par_file)
  
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
  
            if ( l_stop_time ) then 
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
  
               call safeOpen(nLF,log_file)
  
               !--- Write end-energies including energy density:
               !    plus info on movie frames in to STDOUT and log-file
               write(*,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6,/,A,4ES16.6)')        &
                 & " ! Energies at end of time integration:",                     &
                 & " !  (total,poloidal,toroidal,total density)",                 &
                 & " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,    &
                 & " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc,    &
                 & " !  IC mag. energies:",e_mag_ic,e_mag_p_ic,e_mag_t_ic,e_mag_ic/vol_ic
  
               write(nLF,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6,/,A,4ES16.6)')     &
                 & " ! Energies at end of time integration:",                    &
                 & " !  (total,poloidal,toroidal,total density)",                &
                 & " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,   &
                 & " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc,   &
                 & " !  IC mag. energies:",e_mag_ic,e_mag_p_ic,e_mag_t_ic,e_mag_ic/vol_ic
  
               write(nLF,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6)')                 &
                 & " ! Time averaged energies :",                                &
                 & " !  (total,poloidal,toroidal,total density)",                &
                 & " !  Kinetic energies:",e_kin_pMean+e_kin_tMean,e_kin_pMean,  &
                 &                         e_kin_tMean,(e_kin_pMean+e_kin_tMean)/&
                 &                         vol_oc,                               &
                 & " !  OC mag. energies:",e_mag_pMean+e_mag_tMean,e_mag_pMean,  &
                 &                         e_mag_tMean,(e_mag_pMean+e_mag_tMean)/&
                 &                         vol_oc
  
               write(nLF,'(1p,/,A,7(/,A,ES12.4),/,A,4ES12.4,/,A,2ES12.4,/,A,2ES12.4)') &
                 & " ! Time averaged property parameters :",                           &
                 & " !  Rm (Re)         :",RmMean,                                     &
                 & " !  Elsass          :",ElMean,                                     &
                 & " !  Elsass at CMB   :",ElCmbMean,                                  &
                 & " !  Rol             :",RolMean,                                    &
                 & " !  Geos            :",GeosMean,                                   &
                 & " !  Dip             :",DipMean,                                    &
                 & " !  DipCMB          :",DipCMBMean,                                 &
                 & " !  l,m,p,z V scales:",dlVMean,dmVMean,dpVMean,dzVmean,            &
                 & " !  l,m, B scales   :",dlBMean,dmBMean,                            &
                 & " !  vis, Ohm scale  :",lvDissMean,lbDissMean
  
               call safeClose(nLF)
  
            end if ! l_stop_time ?
  
         end if ! l_log
  
         if ( l_Bpot )                                       &
              &     call storePot(time,b,aj,b_ic,aj_ic,      &
              &        nBpotSets,'Bpot.',omega_ma,omega_ic)
         if ( l_Vpot )                                       &
              &     call storePot(time,w,z,b_ic,aj_ic,       &
              &        nVpotSets,'Vpot.',omega_ma,omega_ic)
         if ( l_Tpot )                                       &
              &     call storePot(time,s,z,b_ic,aj_ic,       &
              &        nTpotSets,'Tpot.',omega_ma,omega_ic)
         
         !----- Store current solution
         !      Note: unless l_new_rst_file=.true. .and. .not.l_stop_time
         !            this is written into rst_end.TAG
  
#ifndef WITH_HDF5
         if ( l_store ) then
!#ifdef WITH_HDF5
  !          if ( l_stop_time .or. .not.l_new_rst_file ) then
  !             rst_file='ser_h5_rst_end.'//tag
  !          else if ( l_new_rst_file ) then
  !             call dble2str(time,string)
  !             rst_file='h5_rst_t='//trim(string)//'.'//tag
  !          end if
  !          call storeHdf5_serial(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic, &
  !                                  dwdtLast,dzdtLast,dpdtLast,          &
  !                                  dsdtLast,dbdtLast,djdtLast,          &
  !                                  dbdt_icLast,djdt_icLast)
!#else
            PERFON('out_rst')
            if ( l_stop_time .or. .not.l_new_rst_file ) then
               rst_file="rst_end."//tag
            else if ( l_new_rst_file ) then
               call dble2str(time,string)
               rst_file='rst_t='//trim(string)//'.'//tag
            end if
  
            open(n_rst_file, file=rst_file, status='unknown', form='uNformatted')
            call store(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic,dwdtLast,dzdtLast, &
                       dpdtLast,dsdtLast,dbdtLast,djdtLast,dbdt_icLast,djdt_icLast)
            close(n_rst_file)
!#endif
  
            write(*,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)')&
                 & " ! Storing restart file:",           &
                 & "             at time=",time,         &
                 & "            step no.=",n_time_step,  &
                 & "           into file=",rst_file
            call safeOpen(nLF,log_file)
            
            write(nLF,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)') &
                 & " ! Storing restart file:",              &
                 & "             at time=",time,            &
                 & "            step no.=",n_time_step,     &
                 & "           into file=",rst_file
            call safeClose(nLF)
            PERFOFF
         end if
#endif
         
         if ( l_SRIC .and. l_stop_time ) call outOmega(z,omega_ic)
         
         PERFOFF
      end if

      !----- Output of axisymm. rotation rate for potential vorticity analysis:
      !  NOTE: For l_stop_time=.true. outPV transforms the fields without 
      !        transforming them back. This must thus be the very last 
      !        thing done with them. 
      if ( l_PVout ) call outPV(time,l_stop_time,nPVsets,                   &
           &                    w_LMloc,dw_LMloc,ddw_LMloc,z_LMloc,dz_LMloc,&
           &                    omega_ic,omega_ma)
         
  
      if ( l_log ) then
         timePassedLog=0.0_cp
      end if
  
      if ( lRmsCalc ) then
         call zeroRms
      end if
      
   end subroutine output
!----------------------------------------------------------------------------
end module output_mod
