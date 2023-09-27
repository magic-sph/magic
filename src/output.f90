#include "perflib_preproc.cpp"
module output_mod
   !
   ! This module handles the calls to the different output routines.
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use parallel_mod
   use truncation, only: n_r_max, n_r_ic_max, minc, l_max, l_maxMag, &
       &                 n_r_maxMag, lm_max
   use radial_functions, only: or1, or2, r, rscheme_oc, r_cmb, r_icb,  &
       &                       orho1, sigma
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag,    &
       &                  n_r_cmb, n_r_icb
   use physical_parameters, only: opm,ek,ktopv,prmag,nVarCond,LFfac,ekScaled
   use num_param, only: tScale,eScale
   use blocking, only: st_map, lm2, lo_map, llm, ulm, llmMag, ulmMag
   use horizontal_data, only: hdif_B, dPl0Eq
   use logic, only: l_average, l_mag, l_power, l_anel, l_mag_LF, lVerbose,    &
       &            l_dtB, l_RMS, l_r_field, l_r_fieldT, l_r_fieldXi,         &
       &            l_SRIC, l_cond_ic,l_rMagSpec, l_movie_ic, l_store_frame,  &
       &            l_cmb_field, l_dt_cmb_field, l_save_out, l_non_rot,       &
       &            l_perpPar, l_energy_modes, l_heat, l_hel, l_par,          &
       &            l_chemical_conv, l_movie, l_full_sphere, l_spec_avg,      &
       &            l_phase_field, l_hemi
   use fields, only: omega_ic, omega_ma, b_ic,db_ic, ddb_ic, aj_ic, dj_ic,   &
       &             w_LMloc, dw_LMloc, ddw_LMloc, p_LMloc, xi_LMloc,        &
       &             s_LMloc, ds_LMloc, z_LMloc, dz_LMloc, b_LMloc,          &
       &             db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc, ddj_LMloc,     &
       &             b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc,     &
       &             dj_ic_LMloc, ddj_ic_LMloc, xi_LMloc, dp_LMloc,          &
       &             dxi_LMloc,w_Rloc,z_Rloc,p_Rloc,s_Rloc,xi_Rloc,b_Rloc,   &
       &             aj_Rloc, bICB, phi_Rloc, phi_LMloc
   use fieldsLast, only: dwdt, dzdt, dpdt, dsdt, dbdt, djdt, dbdt_ic, dphidt,&
       &                 djdt_ic, dxidt, domega_ic_dt, domega_ma_dt,         &
       &                 lorentz_torque_ma_dt, lorentz_torque_ic_dt
   use kinetic_energy, only: get_e_kin, get_u_square
   use magnetic_energy, only: get_e_mag
   use fields_average_mod, only: fields_average
   use spectra, only: spectrum, get_amplitude
   use outTO_mod, only: outTO
   use output_data, only: tag, l_max_cmb, n_log_file, log_file
   use constants, only: vol_oc, vol_ic, mass, surf_cmb, two, three, zero
   use outMisc_mod, only: outHeat, outHelicity, outHemi, outPhase, get_onset
   use geos, only: outGeos, outOmega
   use outRot, only: write_rot
   use integration, only: rInt_R
   use outPar_mod, only: outPar, outPerpPar
   use graphOut_mod, only: graphOut_IC
   use power, only: get_power
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC,  &
       &                     gather_from_lo_to_rank0
   use out_coeff, only: write_Bcmb, write_coeffs, write_Pot, initialize_coeff, &
       &                finalize_coeff
#ifdef WITH_MPI
   use out_coeff, only: write_Pot_mpi
#endif
   use getDlm_mod, only: getDlm
   use movie_data, only: movie_gather_frames_to_rank0
   use dtB_mod, only: get_dtBLMfinish
   use out_movie, only: write_movie_frame
   use out_movie_IC, only: store_movie_frame_IC
   use RMS, only: zeroRms, dtVrms, dtBrms
   use useful, only:  logWrite
   use radial_spectra  ! rBrSpec, rBpSpec
   use time_schemes, only: type_tscheme
   use storeCheckPoints

   implicit none

   private

   !-- Counter for output files/sets:
   integer :: nPotSets, n_spec
   integer :: n_dt_cmb_sets, n_cmb_setsMov

   !--- For averaging:
   real(cp) :: timePassedLog, timeNormLog
   integer :: nLogs

   real(cp) :: dlBMean,dmBMean
   real(cp) :: lvDissMean,lbDissMean
   real(cp) :: RmMean,ElMean,ElCmbMean,RolMean
   real(cp) :: GeosMean,GeosAMean,GeosZMean,GeosMMean,GeosNAPMean
   real(cp) :: RelA,RelZ,RelM,RelNA
   real(cp) :: DipMean,DipCMBMean
   real(cp) :: dlVMean,dlVcMean,dmVMean,dpVMean,dzVMean

   real(cp) :: eTot,eTotOld,dtEint
   real(cp) :: e_kin_pMean, e_kin_tMean
   real(cp) :: e_mag_pMean, e_mag_tMean
   integer :: n_e_sets

   real(cp) :: timePassedRMS, timeNormRMS
   integer :: nRMS_sets

   integer :: n_dtE_file, n_par_file, n_cmb_file
   integer :: n_cmbMov_file, n_dt_cmb_file
   character(len=72) :: dtE_file, par_file
   character(len=72) :: cmb_file, dt_cmb_file, cmbMov_file

   public :: output, initialize_output, finalize_output

contains

   subroutine initialize_output

      if ( l_r_field .or. l_r_fieldT .or. l_r_fieldXi) call initialize_coeff()

      n_spec       =0
      n_cmb_setsMov=0
      n_dt_cmb_sets=0
      nPotSets     =1
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
      GeosAMean    =0.0_cp
      GeosZMean    =0.0_cp
      GeosMMean    =0.0_cp
      GeosNAPMean  =0.0_cp
      RelA         =0.0_cp
      RelM         =0.0_cp
      RelZ         =0.0_cp
      RelNA        =0.0_cp
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

      par_file='par.'//tag
      if ( l_mag .and. l_cmb_field ) then
         cmb_file   ='B_coeff_cmb.'//tag
         if ( l_movie ) then
            cmbMov_file='B_coeff_cmbMov.'//tag
         end if
      end if

      if ( l_mag .and. l_dt_cmb_field ) then
         dt_cmb_file   ='B_coeff_dt_cmb.'//tag
      end if

      if ( l_power ) then
         dtE_file='dtE.'//tag
      end if

      if ( rank == 0 .and. ( .not. l_save_out ) ) then
         open(newunit=n_par_file, file=par_file, status='new')

         if ( l_mag .and. l_cmb_field ) then
            open(newunit=n_cmb_file, file=cmb_file, &
            &    status='new', form='unformatted')
            if ( l_movie ) then
               open(newunit=n_cmbMov_file, file=cmbMov_file, &
               &    status='new', form='unformatted')
            end if
         end if

         if ( l_mag .and. l_dt_cmb_field ) then
            open(newunit=n_dt_cmb_file, file=dt_cmb_file, &
                 status='new', form='unformatted')
         end if

         if ( l_power ) then
            open(newunit=n_dtE_file, file=dtE_file, status='new')
         end if

      end if

   end subroutine initialize_output
!----------------------------------------------------------------------------
   subroutine finalize_output

      if ( rank == 0 .and. ( .not. l_save_out ) ) then
         if ( l_mag .and. l_cmb_field ) then
            close(n_cmb_file)
            if (l_movie) close(n_cmbMov_file)
         end if
         if ( l_mag .and. l_dt_cmb_field ) then
            close(n_dt_cmb_file)
         end if
         if ( l_power ) close(n_dtE_file)
      end if

      if ( l_r_field .or. l_r_fieldT .or. l_r_fieldXi) call finalize_coeff()

   end subroutine finalize_output
!----------------------------------------------------------------------------
   subroutine output(time,tscheme,n_time_step,l_stop_time,l_pot,l_log,    &
              &      l_graph,lRmsCalc,l_store,l_new_rst_file,lOnsetCalc,  &
              &      l_spectrum,lTOCalc,lTOframe,l_frame,n_frame,l_cmb,   &
              &      n_cmb_sets,l_r,lorentz_torque_ic,lorentz_torque_ma,  &
              &      dbdt_CMB_LMloc)
      !
      !  This subroutine controls most of the output.
      !

      !--- Input of variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      logical,             intent(in) :: l_stop_time
      logical,             intent(in) :: l_pot, lOnsetCalc
      logical,             intent(in) :: l_log, l_graph, lRmsCalc, l_store
      logical,             intent(in) :: l_new_rst_file, l_spectrum
      logical,             intent(in) :: lTOCalc,lTOframe
      logical,             intent(in) :: l_frame, l_cmb, l_r
      integer,             intent(inout) :: n_frame
      integer,             intent(inout) :: n_cmb_sets

      !--- Input of Lorentz torques and dbdt calculated in radialLoopG
      !    Parallelization note: Only the contribution at the CMB must be
      !    collected and is (likely) stored on the processor (#0) that performs
      !    this routine anyway.
      real(cp),            intent(in) :: lorentz_torque_ma,lorentz_torque_ic

      !--- Input of scales fields via common block in fields.f90:
      !    Parallelization note: these fields are LM-distributed.
      !    The input fields HelASr,Hel2ASr,TstrRLM,TadvRLM, and TomeRLM
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
      complex(cp), intent(in) :: dbdt_CMB_LMloc(llmMag:ulmMag)

      !--- Local stuff:
      !--- Energies:
      real(cp) :: ekinR(n_r_max)     ! kinetic energy w radius
      real(cp) :: e_mag,e_mag_ic,e_mag_cmb,e_mag_p,e_mag_t
      real(cp) :: e_mag_p_as,e_mag_t_as,e_mag_p_ic,e_mag_t_ic
      real(cp) :: e_mag_p_as_ic,e_mag_t_as_ic,e_mag_os,e_mag_as_os
      real(cp) :: e_kin,e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as
      real(cp) :: eKinIC,eKinMA,dtE

      integer :: lm,m,l

      !--- Property parameters:
      complex(cp) :: dbdtCMB(llmMag:ulmMag)        ! SV at CMB !
      real(cp) :: dlVR(n_r_max),dlVRc(n_r_max)
      real(cp) :: RolRu2(n_r_max),RmR(n_r_max),dlPolPeakR(n_r_max)
      real(cp) :: Re,Ro,Rm,El,ElCmb,Rol,Geos,GeosA,GeosZ,GeosM,GeosNAP
      real(cp) :: Dip,DipCMB,EC
      real(cp) :: ReConv,RoConv,e_kin_nas,RolC
      real(cp) :: elsAnel,dlVPolPeak,dlBPolPeak
      real(cp) :: dlB,dlBc,dmB
      real(cp) :: dlV,dlVc,dmV,dpV,dzV
      real(cp) :: visDiss,ohmDiss,lvDiss,lbDiss
      real(cp) :: ReEquat,timeScaled,dL
      character(len=96) :: message
      logical :: DEBUG_OUTPUT=.false.

      timeScaled=tScale*time
      timePassedLog=timePassedLog+tscheme%dt(1)

      ! We start with the computation of the energies
      ! in parallel.
      if ( l_log ) then

         nLogs=nLogs+1
         timeNormLog=timeNormLog+timePassedLog

         !----- Write torques and rotation rates:
         PERFON('out_rot')
         call write_rot( time,tscheme%dt(1),eKinIC,eKinMA,w_LMloc,z_LMloc, &
              &          dz_LMloc,b_LMloc,omega_ic,omega_ma,               &
              &          lorentz_torque_ic,lorentz_torque_ma)
         PERFOFF
         if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  write_rot  on rank ",rank

         PERFON('out_ekin')
         n_e_sets=n_e_sets+1
         call get_e_kin(time,.true.,l_stop_time,n_e_sets,w_LMloc,    &
              &         dw_LMloc,z_LMloc,e_kin_p,e_kin_t,e_kin_p_as, &
              &         e_kin_t_as,ekinR)
         e_kin=e_kin_p+e_kin_t
         e_kin_nas=e_kin-e_kin_p_as-e_kin_t_as
         if ( DEBUG_OUTPUT ) write(output_unit,"(A,I6)") "Written  e_kin  on rank ",rank

         call get_e_mag(time,.true.,l_stop_time,n_e_sets,b_LMloc,db_LMloc, &
              &         aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc,       &
              &         e_mag_p,e_mag_t,e_mag_p_as,e_mag_t_as,e_mag_p_ic,  &
              &         e_mag_t_ic,e_mag_p_as_ic,e_mag_t_as_ic,            &
              &         e_mag_os,e_mag_as_os,e_mag_cmb,Dip,DipCMB,elsAnel )
         e_mag   =e_mag_p+e_mag_t
         e_mag_ic=e_mag_p_ic+e_mag_t_ic
         PERFOFF
         if ( DEBUG_OUTPUT ) write(output_unit,"(A,I6)") "Written  e_mag  on rank ",rank

         !----- Calculate distribution of energies on all m's
         if ( l_energy_modes ) then
            PERFON('out_amplitude')
            call get_amplitude(time,w_LMloc,dw_LMloc,z_LMloc,b_LMloc,&
                 &             db_LMloc,aj_LMloc)
            PERFOFF
            if ( DEBUG_OUTPUT ) &
               & write(output_unit,"(A,I6)") "Written  amplitude  on rank ",rank
         endif

         !---- Surface zonal velocity at the equator
         if ( ktopv==1 ) then
            ReEquat=0.0_cp
            do lm=llm,ulm
               l = lo_map%lm2l(lm)
               m = lo_map%lm2m(lm)
               if ( m == 0 ) then
                  ReEquat=ReEquat-real(z_LMloc(lm,n_r_cmb))*dPl0Eq(l+1)*or1(n_r_cmb)
               end if
            end do
#ifdef WITH_MPI
            call MPI_AllReduce(MPI_IN_PLACE, ReEquat, 1, MPI_DEF_REAL, MPI_SUM, &
                 &             MPI_COMM_WORLD, ierr)
#endif
         else
            ReEquat=0.0_cp
         end if

         if ( l_spec_avg ) then
            call spectrum(-1,time,.true.,nLogs,l_stop_time,timePassedLog,      &
                 &        timeNormLog,s_LMloc,ds_LMloc,xi_LMloc,dxi_LMloc,     &
                 &        phi_LMloc,w_LMloc,dw_LMloc,z_LMloc,b_LMloc,db_LMloc, &
                 &        aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc)
         end if

         if ( l_average ) then
            PERFON('out_aver')
            call fields_average(time,tscheme,nLogs,l_stop_time,timePassedLog,  &
                 &              timeNormLog,omega_ic,omega_ma,w_LMloc,z_LMloc, &
                 &              p_LMloc,s_LMloc,xi_LMloc,phi_LMloc,b_LMloc,    &
                 &              aj_LMloc,b_ic_LMloc,aj_ic_LMloc)
            PERFOFF
            if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  averages  on rank ",rank
         end if

         if ( l_power ) then

            PERFON('out_pwr')
            if ( rank == 0 ) then
               if ( nLogs > 1 ) then
                  if ( l_save_out ) then
                     open(newunit=n_dtE_file, file=dtE_file, &
                     &    status='unknown', position='append')
                  end if
                  eTotOld=eTot
                  eTot   =e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
                  dtE    =(eTot-eTotOld)/timePassedLog
                  dtEint =dtEint+timePassedLog*(eTot-eTotOld)
                  write(n_dtE_file,'(ES20.12,3ES16.6)') timeScaled,dtE,   &
                  &     dtEint/timeNormLog,dtE/eTot
                  if ( l_save_out ) close(n_dtE_file)
               else
                  eTot  =e_kin+e_mag+e_mag_ic+e_mag_os+eKinIC+eKinMA
                  dtEint=0.0_cp
               end if
            end if

            call get_power( time,timePassedLog,timeNormLog,l_stop_time,      &
                 &          omega_ic,omega_ma,lorentz_torque_ic,             &
                 &          lorentz_torque_ma,w_LMloc,z_LMloc,               &
                 &          dz_LMloc,s_LMloc,xi_LMloc,                       &
                 &          b_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,db_ic_LMloc, &
                 &          ddb_ic_LMloc,aj_ic_LMloc,dj_ic_LMloc,            &
                 &          visDiss,ohmDiss)
            PERFOFF
            if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  power  on rank ",rank
         end if

         !----- If anelastic additional u**2 outputs
         if ( l_anel ) then
            call get_u_square(time,w_LMloc,dw_LMloc,z_LMloc,RolRu2)
         else
            RolRu2(:)=0.0_cp
         end if

         !-- Get flow lengthscales
         call getDlm(w_LMloc,dw_LMloc,z_LMloc,dlV,dlVR,dmV,dlVc,dlVPolPeak, &
              &      dlVRc,dlPolPeakR,'V')

         !-- Out radial profiles of parameters
         call outPar(s_LMloc, ds_LMloc, xi_LMloc, p_LMloc, dp_LMloc, timePassedLog, &
              &      timeNormLog, l_stop_time, ekinR, RolRu2, dlVR, dlVRc,          &
              &      dlPolPeakR, RmR)

         !-- Perpendicular/parallel
         if ( l_perpPar ) call outPerpPar(time,timePassedLog,timeNormLog,l_stop_time)

         if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  outPar  on rank ",rank

         if ( l_heat .or. l_chemical_conv ) then
            call outHeat(timeScaled,timePassedLog,timeNormLog,l_stop_time, &
                 &       s_LMloc,ds_LMloc,p_LMloc,xi_LMloc,                &
                 &       dxi_LMloc)
         end if

         if ( l_hel ) call outHelicity(timeScaled)

         if ( l_hemi ) call outHemi(timeScaled)

         if ( l_phase_field ) call outPhase(timeScaled,timePassedLog,       &
                                   &        timeNormLog,l_stop_time,nLogs,  &
                                   &        s_LMloc,ds_LMloc,phi_LMloc)

         if ( l_par ) then
            call outGeos(timeScaled,Geos,GeosA,GeosZ,GeosM,GeosNAP,EC)
            dpV=0.0_cp ! To be handled later
            dzV=0.0_cp
         else
            Geos   =0.0_cp
            GeosA  =0.0_cp
            GeosZ  =0.0_cp
            GeosM  =0.0_cp
            GeosNAP=0.0_cp
            dpV    =0.0_cp
            dzV    =0.0_cp
            EC     =0.0_cp ! test kinetic energy
         end if

         if ( l_mag .or. l_mag_LF ) then
            !-- Get magnetic field lengthscales
            call getDlm(b_LMloc,db_LMloc,aj_LMloc,dlB,dlVR,dmB, &
                 &      dlBc,dlBPolPeak,dlVRc,dlPolPeakR,'B')
         else
            dlB=0.0_cp
            dmB=0.0_cp
         end if
      end if

      if ( l_spectrum ) then
         n_spec=n_spec+1
         call spectrum(n_spec,time,.false.,nLogs,l_stop_time,timePassedLog, &
              &        timeNormLog,s_LMloc,ds_LMloc,xi_LMloc,dxi_LMloc,     &
              &        phi_LMloc,w_LMloc,dw_LMloc,z_LMloc,b_LMloc,db_LMloc, &
              &        aj_LMloc,b_ic_LMloc,db_ic_LMloc,aj_ic_LMloc)
         if ( rank == 0 ) then
            write(output_unit,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)')&
            &    " ! Storing spectra:",                            &
            &    "             at time=",timeScaled,               &
            &    "            step no.=",n_time_step

            if ( l_save_out ) then
               open(newunit=n_log_file, file=log_file, status='unknown', &
               &    position='append')
            end if
            write(n_log_file,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)') &
            &    " ! Storing spectra:",                            &
            &    "             at time=",timeScaled,               &
            &    "            step no.=",n_time_step
            if ( l_save_out ) close(n_log_file)
         end if
         if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  spectrum  on rank ",rank
      end if

      if ( lTOCalc ) then
         !------ Output for every log time step:
         if ( lVerbose ) write(output_unit,*) '! Calling outTO !'
         if ( .not. l_log ) then
            call get_e_kin(time,.false.,l_stop_time,0,w_LMloc,dw_LMloc,  &
                 &         z_LMloc,e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,&
                 &         ekinR)
            e_kin=e_kin_p+e_kin_t
         end if
         call outTO(time,n_time_step,e_kin,e_kin_t_as,lTOframe)

         if ( lVerbose ) write(output_unit,*) '! outTO finished !'
         if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  TO  on rank ",rank
      end if

      !--- Get radial derivatives and add dt dtB terms:
      if ( l_dtB ) then
         call get_dtBLMfinish(time,n_time_step,omega_ic,b_LMloc,ddb_LMloc, &
              &               aj_LMloc,dj_LMloc,ddj_LMloc,b_ic_LMloc,      &
              &               db_ic_LMloc,ddb_ic_LMloc,aj_ic_LMloc,        &
              &               dj_ic_LMloc,ddj_ic_LMloc,l_frame)
      end if

      !-- Compute growth rates and drift frequencies of several wavenumbers
      if ( lOnsetCalc .and. (.not. l_stop_time) ) then
         call get_onset(time, w_LMloc, tscheme%dt(1), l_log, nLogs)
      end if

      if ( l_RMS ) then
         if ( n_time_step == 1 ) then
            nRMS_sets    =0
            timeNormRMS  =0.0_cp
            timePassedRMS=0.0_cp
            call zeroRms()
         end if
         timePassedRMS=timePassedRMS+tscheme%dt(1)
         if ( lRmsCalc ) then
            if ( lVerbose ) write(output_unit,*) '! Writing RMS output !'
            timeNormRMS=timeNormRMS+timePassedRMS
            call dtVrms(time,nRMS_sets,timePassedRMS,timeNormRMS,l_stop_time)
            if ( l_mag ) call dtBrms(time)
            timePassedRMS=0.0_cp
         end if
         if (DEBUG_OUTPUT) write(output_unit,"(A,I6)") "Written  dtV/Brms  on rank ",rank
      end if

      !--- Store poloidal magnetic coeffs at cmb
      if ( l_cmb ) then
         PERFON('out_cmb')
         call write_Bcmb(timeScaled,b_LMloc(:,n_r_cmb),l_max_cmb,n_cmb_sets,   &
              &          cmb_file,n_cmb_file)

         !--- Store SV of poloidal magnetic coeffs at cmb
         if ( l_dt_cmb_field ) then
            !nR=8! at CMB dbdt=induction=0, only diffusion !
            do lm=max(2,llm),ulm
               l  =lo_map%lm2l(lm)
               dL =real(l*l+1,cp)
               dbdtCMB(lm)= dbdt_CMB_LMloc(lm)/                                        &
               &         (dL*or2(n_r_cmb)) + opm*hdif_B(l) * ( ddb_LMloc(lm,n_r_cmb) - &
               &          dL*or2(n_r_cmb)*b_LMloc(lm,n_r_cmb) )
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
      if ( l_r ) call write_coeffs(w_LMloc, dw_LMloc, ddw_LMloc, z_LMLoc, b_LMLoc,  &
                      &            db_LMloc, ddb_LMloc, aj_LMloc, s_LMloc,          &
                      &            xi_LMloc, timeScaled)

      if ( l_pot ) then
#ifdef WITH_MPI
         call write_Pot_mpi(time,w_Rloc,z_Rloc,b_ic_LMloc,aj_ic_LMloc, &
              &             nPotSets,'V_lmr.',omega_ma,omega_ic)
         if ( l_heat ) then
           call write_Pot_mpi(time,s_Rloc,z_Rloc,b_ic_LMloc,aj_ic_LMloc, &
                &             nPotSets,'T_lmr.',omega_ma,omega_ic)
         end if
         if ( l_chemical_conv ) then
           call write_Pot_mpi(time,xi_Rloc,z_Rloc,b_ic_LMloc,aj_ic_LMloc, &
                &             nPotSets,'Xi_lmr.',omega_ma,omega_ic)
         end if
         if ( l_mag ) then
            call write_Pot_mpi(time,b_Rloc,aj_Rloc,b_ic_LMloc,aj_ic_LMloc, &
                 &             nPotSets,'B_lmr.',omega_ma,omega_ic)
         end if
#else
         call write_Pot(time,w_LMloc,z_LMloc,b_ic_LMloc,aj_ic_LMloc, &
              &         nPotSets,'V_lmr.',omega_ma,omega_ic)
         if ( l_heat ) then
           call write_Pot(time,s_LMloc,z_LMloc,b_ic_LMloc,aj_ic_LMloc, &
                &         nPotSets,'T_lmr.',omega_ma,omega_ic)
         end if
         if ( l_chemical_conv ) then
           call write_Pot(time,xi_LMloc,z_LMloc,b_ic_LMloc,aj_ic_LMloc, &
                &         nPotSets,'Xi_lmr.',omega_ma,omega_ic)
         end if
         if ( l_mag ) then
            call write_Pot(time,b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc, &
                 &         nPotSets,'B_lmr.',omega_ma,omega_ic)
         end if
#endif
         nPotSets=nPotSets+1
      end if

      !--- Write spectra output that has partially been calculated in LMLoop
      if ( l_rMagSpec .and. n_time_step > 1 ) then
         if ( l_frame ) then
            call rBrSpec(time,b_LMloc, b_ic_LMloc ,'rBrSpecMov',.true.,lo_map)
            call rBpSpec(time,aj_LMloc,aj_ic_LMloc,'rBpSpecMov',.true.,lo_map)
         end if
         if ( l_log ) then
            call rBrSpec(time,b_LMloc, b_ic_LMloc ,'rBrSpec',.true.,lo_map)
            call rBpSpec(time,aj_LMloc,aj_ic_LMloc,'rBpSpec',.true.,lo_map)
         end if
      end if

      !
      ! Writing of the restart file
      !
      if ( l_store ) then
#ifdef WITH_MPI
         call store_mpi(time,tscheme,n_time_step,l_stop_time,l_new_rst_file,  &
              &         .false.,w_Rloc,z_Rloc,p_Rloc,s_Rloc,xi_Rloc,phi_Rloc, &
              &         b_Rloc,aj_Rloc,b_ic_LMloc,aj_ic_LMloc,dwdt,dzdt,dpdt, &
              &         dsdt,dxidt,dphidt,dbdt,djdt,dbdt_ic,djdt_ic,          &
              &         domega_ma_dt,domega_ic_dt,lorentz_torque_ma_dt,       &
              &         lorentz_torque_ic_dt)
#else
         call store(time,tscheme,n_time_step,l_stop_time,l_new_rst_file,.false., &
              &     w_LMloc,z_LMloc,p_LMloc,s_LMloc,xi_LMloc,phi_LMloc,b_LMloc,  &
              &     aj_LMloc,b_ic_LMloc,aj_ic_LMloc,dwdt,dzdt,dpdt,dsdt,dxidt,   &
              &     dphidt,dbdt,djdt,dbdt_ic,djdt_ic,domega_ma_dt,domega_ic_dt,  &
              &     lorentz_torque_ma_dt,lorentz_torque_ic_dt)
#endif
      end if


      ! ===================================================
      !      GATHERING for output
      ! ===================================================
      ! We have all fields in LMloc space. Thus we gather the whole fields on rank 0.

      if ( l_frame .or. (l_graph .and. l_mag .and. n_r_ic_max > 0) ) then
         PERFON('out_comm')

         if ( l_mag ) call gather_from_lo_to_rank0(b_LMloc(:,n_r_icb),bICB)

         if ( l_cond_ic ) then
            call gather_all_from_lo_to_rank0(gt_IC,b_ic_LMloc,b_ic)
            call gather_all_from_lo_to_rank0(gt_IC,db_ic_LMloc,db_ic)
            call gather_all_from_lo_to_rank0(gt_IC,ddb_ic_LMloc,ddb_ic)

            call gather_all_from_lo_to_rank0(gt_IC,aj_ic_LMloc,aj_ic)
            call gather_all_from_lo_to_rank0(gt_IC,dj_ic_LMloc,dj_ic)
         else if ( l_mag ) then ! Set to zero (compat with Leg TF)
            if ( rank == 0 ) then
               db_ic(:,1)=zero
               aj_ic(:,1)=zero
               if ( l_frame ) then
                  ddb_ic(:,1)=zero
                  dj_ic(:,1) =zero
               end if
            end if
         end if

         ! for writing a restart file, we also need the d?dtLast arrays,
         ! which first have to be gathered on rank 0
         PERFOFF

      end if

      !--- Movie output and various supplementary things:
      if ( l_frame ) then
         ! The frames array for the movies is distributed over the ranks
         ! and has to be gathered on rank 0 for output.

         ! Each movie uses some consecutive frames in the frames array. They
         ! start at n_movie_field_start(1,n_movie)
         ! up to    n_movie_field_stop(1+n_fields_oc+n_fields,n_movie) (n_fields_ic>0
         ! or       n_movie_field_stop(1+n_fields,n_movie)             (n_fields_ic=0)

         call movie_gather_frames_to_rank0()

         if ( l_movie_ic .and. l_store_frame .and. rank == 0 ) then
            call store_movie_frame_IC(bICB,b_ic,db_ic,ddb_ic,aj_ic,dj_ic)
         end if

         n_frame=n_frame+1
         call logWrite(' ')
         if ( rank == 0 ) then
            write(message,'(1p,A,I8,A,ES16.6,I8)')            &
            &      " ! WRITING MOVIE FRAME NO ",n_frame,      &
            &      "       at time/step",timeScaled,n_time_step
         end if
         call logWrite(message)

         !--- Storing the movie frame:
         call write_movie_frame(n_frame,timeScaled,b_LMloc,db_LMloc,aj_LMloc, &
              &                 dj_LMloc,b_ic,db_ic,aj_ic,dj_ic,omega_ic,     &
              &                 omega_ma)
      end if

      !----- Plot out inner core magnetic field, outer core
      !      field has been written in radialLoop !
      if ( l_graph .and. l_mag .and. n_r_ic_max > 0 ) then
         call graphOut_IC(b_ic,db_ic,aj_ic,bICB)
      end if

      ! =======================================================================
      ! ======= compute output on rank 0 ==============
      ! =======================================================================
      if ( rank == 0 ) then
         PERFON('out_out')

         if ( l_log ) then
            !--- Energies and rotation info and a lot of other stuff
            !    performed for l_log=.true.

            !----- Getting the property parameters:
            Re     = sqrt(two*e_kin/vol_oc/eScale)/sqrt(mass)
            if ( ( abs(e_kin_nas) <= 10.0_cp * epsilon(mass) ) .or. &
            &     e_kin_nas < 0.0_cp ) e_kin_nas=0.0_cp
            ReConv = sqrt(two*e_kin_nas/vol_oc/eScale)/sqrt(mass)

            if ( l_non_rot ) then
               Ro=0.0_cp
               RoConv=0.0_cp
            else
               Ro=Re*ekScaled
               RoConv=ReConv*ekScaled
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
            !write(output_unit,"(A,3ES20.12)") "dlVc,RoConv,RolC = ",dlVc,RoConv,RolC

            if ( prmag /= 0 .and. nVarCond > 0 ) then
               Rm=0.0_cp
               Rm=rInt_R(RmR,r,rscheme_oc)
               Rm=three*Rm/(r_cmb**3-r_icb**3)
            elseif ( prmag /= 0 ) then
               Rm=Re*prmag
            else
               Rm=Re
            end if
            !El   =two*e_mag/vol_oc/LFfac
            ! Elsasser number is computed from the averaged profile
            if ( l_mag .or. l_mag_LF ) then
               El   =elsAnel/vol_oc
               ElCmb=two*e_mag_cmb/surf_cmb/LFfac*sigma(n_r_cmb)*orho1(n_r_cmb)
               ! Elsasser must not depend of timescale
               ElCmb=ElCmb/eScale
            else
               El   =0.0_cp
               ElCmb=0.0_cp
            end if
            if ( l_power ) then
               if ( visDiss /= 0.0_cp ) then
                  lvDiss=sqrt(two*e_kin/abs(visDiss)) ! Viscous dissipation lengthscale
               else
                  lvDiss=0.0_cp
               end if
               if ( l_mag .or. l_mag_LF ) then
                  if ( ohmDiss /= 0.0_cp ) then
                     lbDiss=sqrt(two*opm*(e_mag+e_mag_ic)/abs(ohmDiss)) ! Ohmic dissipation lengthscale
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
               open(newunit=n_par_file, file=par_file, status='unknown', &
               &    position='append')
            end if
            write(n_par_file,'(ES20.12,19ES16.8)')  &
                 &             timeScaled,          &! 1) time
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
                 &             dlVPolPeak,          &! 19) Peak of the poloidal energy
                 &                ReEquat            ! 20) CMB flow at the equator

            if ( l_save_out ) close(n_par_file)

            !---- Building time mean:
            RmMean     =RmMean     +timePassedLog*Rm
            ElMean     =ElMean     +timePassedLog*El
            ElCmbMean  =ElCmbMean  +timePassedLog*ElCmb
            RolMean    =RolMean    +timePassedLog*Rol
            GeosMean   =GeosMean   +timePassedLog*Geos
            GeosAMean  =GeosAMean  +timePassedLog*GeosA
            GeosZMean  =GeosZMean  +timePassedLog*GeosZ
            GeosMMean  =GeosMMean  +timePassedLog*GeosM
            GeosNAPMean=GeosNAPMean+timePassedLog*GeosNAP
            if ( e_kin > 0.0_cp ) then
               ! Relative axisymmetric kinetic energy:
               RelA       =RelA       +timePassedLog*(e_kin_p_as+e_kin_t_as)/e_kin
               ! Relative zonal kinetic energy:
               RelZ       =RelZ       +timePassedLog*e_kin_t_as/e_kin
               ! Relative meridional kinetic energy:
               RelM       =RelM       +timePassedLog*e_kin_p_as/e_kin
               ! Relative non-axisymmetric kinetic energy:
               RelNA      =RelNA      +timePassedLog*(e_kin-e_kin_p_as-e_kin_t_as)/e_kin
            end if
            DipMean    =DipMean    +timePassedLog*Dip
            DipCMBMean =DipCMBMean +timePassedLog*DipCMB
            e_kin_pMean=e_kin_pMean+timePassedLog*e_kin_p
            e_kin_tMean=e_kin_tMean+timePassedLog*e_kin_t
            e_mag_pMean=e_mag_pMean+timePassedLog*e_mag_p
            e_mag_tMean=e_mag_tMean+timePassedLog*e_mag_t
            dlVMean    =dlVMean    +timePassedLog*dlV
            dlVcMean   =dlVcMean   +timePassedLog*dlVc
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
               GeosAMean  =GeosAMean/timeNormLog
               GeosZMean  =GeosZMean/timeNormLog
               GeosMMean  =GeosMMean/timeNormLog
               GeosNAPMean=GeosNAPMean/timeNormLog
               RelA       =RelA/timeNormLog
               RelZ       =RelZ/timeNormLog
               RelM       =RelM/timeNormLog
               RelNA      =RelNA/timeNormLog
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

               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown', &
                  &    position='append')
               end if

               !--- Write end-energies including energy density:
               !    plus info on movie frames in to STDOUT and log-file
               if ( l_full_sphere ) then
                  write(output_unit,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6)')    &
                  & " ! Energies at end of time integration:",                 &
                  & " !  (total,poloidal,toroidal,total density)",             &
                  & " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,&
                  & " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc
                  write(n_log_file,                                               &
                  &    '(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6)')                    &
                  &    " ! Energies at end of time integration:",                 &
                  &    " !  (total,poloidal,toroidal,total density)",             &
                  &    " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,&
                  &    " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc

               else
                  write(output_unit,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6,/,A,4ES16.6)')  &
                  & " ! Energies at end of time integration:",                 &
                  & " !  (total,poloidal,toroidal,total density)",             &
                  & " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,&
                  & " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc,&
                  & " !  IC mag. energies:",e_mag_ic,e_mag_p_ic,e_mag_t_ic,    &
                  & e_mag_ic/vol_ic

                  write(n_log_file,                                               &
                  &    '(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6,/,A,4ES16.6)')        &
                  &    " ! Energies at end of time integration:",                 &
                  &    " !  (total,poloidal,toroidal,total density)",             &
                  &    " !  Kinetic energies:",e_kin,e_kin_p,e_kin_t,e_kin/vol_oc,&
                  &    " !  OC mag. energies:",e_mag,e_mag_p,e_mag_t,e_mag/vol_oc,&
                  &    " !  IC mag. energies:",e_mag_ic,e_mag_p_ic,e_mag_t_ic,    &
                  &    e_mag_ic/vol_ic
               end if

               write(n_log_file,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6)')        &
               & " ! Time averaged energies :",                                &
               & " !  (total,poloidal,toroidal,total density)",                &
               & " !  Kinetic energies:",e_kin_pMean+e_kin_tMean,e_kin_pMean,  &
               &                         e_kin_tMean,(e_kin_pMean+e_kin_tMean)/&
               &                         vol_oc,                               &
               & " !  OC mag. energies:",e_mag_pMean+e_mag_tMean,e_mag_pMean,  &
               &                         e_mag_tMean,(e_mag_pMean+e_mag_tMean)/&
               &                         vol_oc

               write(n_log_file,                                                &
               & '(1p,/,A,15(/,A,ES12.4),/,A,4ES12.4,/,A,2ES12.4,/,A,2ES12.4)') &
               & " ! Time averaged property parameters :",                      &
               & " !  Rm (Re)          :",RmMean,                               &
               & " !  Elsass           :",ElMean,                               &
               & " !  Elsass at CMB    :",ElCmbMean,                            &
               & " !  Rol              :",RolMean,                              &
               & " !  rel AS  Ekin     :",RelA,                                 &
               & " !  rel Zon Ekin     :",RelZ,                                 &
               & " !  rel Mer Ekin     :",RelM,                                 &
               & " !  rel NA  Ekin     :",RelNA,                                &
               & " !  rel geos Ekin    :",GeosMean,                             &
               & " !  rel geos AS Ekin :",GeosAMean,                            &
               & " !  rel geos Zon Ekin:",GeosZMean,                            &
               & " !  rel geos Mer Ekin:",GeosMMean,                            &
               & " !  rel geos NAP Ekin:",GeosNAPMean,                          &
               & " !  Dip              :",DipMean,                              &
               & " !  DipCMB           :",DipCMBMean,                           &
               & " !  l,m,p,z V scales :",dlVMean,dmVMean,dpVMean,dzVmean,      &
               & " !  l,m, B scales    :",dlBMean,dmBMean,                      &
               & " !  vis, Ohm scale   :",lvDissMean,lbDissMean
               if ( l_par ) then
                  write(n_log_file,*) !' Calculating geostrophic contributions with outEgeos.f90'
                  write(n_log_file,*) '! precision of cyl. transf.  (geos):',abs(EC/e_kin-1)
               end if

               if ( l_save_out ) close(n_log_file)

            end if ! l_stop_time ?

         end if ! l_log

         PERFOFF
      end if

      if ( l_SRIC .and. l_stop_time ) call outOmega(z_LMloc,omega_ic)

      if ( l_log ) timePassedLog=0.0_cp

      if ( lRmsCalc ) call zeroRms()

   end subroutine output
!----------------------------------------------------------------------------
end module output_mod
