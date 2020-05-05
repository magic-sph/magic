#include "perflib_preproc.cpp"
module output_mod

   use precision_mod
   use parallel_mod
   use truncation, only: n_r_max, n_r_ic_max, minc, l_max, l_maxMag, &
       &                 n_r_maxMag, lm_max, nRstart, nRstop,        &
       &                 nRstartMag, nRstopMag, n_r_cmb, n_r_icb,    &
       &                 n_mlo_loc, n_mloMag_loc
   use radial_functions, only: or1, or2, r, rscheme_oc, r_cmb, r_icb,  &
       &                       orho1, sigma
   use physical_parameters, only: opm,ek,ktopv,prmag,nVarCond,LFfac,ekScaled
   use LMmapping, only : map_mlo
   use num_param, only: tScale, eScale
   use horizontal_data, only: dPl0Eq, hdif_B
   use logic, only: l_average, l_mag, l_power, l_anel, l_mag_LF, lVerbose, &
       &            l_dtB, l_RMS, l_r_field, l_r_fieldT, l_SRIC,           &
       &            l_cond_ic,l_rMagSpec, l_movie_ic, l_store_frame,       &
       &            l_cmb_field, l_dt_cmb_field, l_save_out, l_non_rot,    &
       &            l_perpPar, l_energy_modes, l_heat, l_hel, l_par,       &
       &            l_chemical_conv, l_movie, l_full_sphere, l_spec_avg
   use fields, only: omega_ic, omega_ma, b_ic,db_ic, ddb_ic, aj_ic, dj_ic,   &
       &             ddj_ic, w_LMloc, dw_LMloc, ddw_LMloc, p_LMloc, xi_LMloc,&
       &             s_LMloc, ds_LMloc, z_LMloc, dz_LMloc, b_LMloc,          &
       &             db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc, ddj_LMloc,     &
       &             b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc,     &
       &             dj_ic_LMloc, ddj_ic_LMloc, dp_LMloc, xi_LMloc,          &
       &             dxi_LMloc,w_Rloc,z_Rloc,p_Rloc,s_Rloc,xi_Rloc,b_Rloc,   &
       &             aj_Rloc, bICB,  & !@>TODO rmove the LMloc at some point
       &             w_Rdist,z_Rdist,p_Rdist,s_Rdist,xi_Rdist,b_Rdist,  aj_Rdist, &
       &             ddj_ic, w_LMdist, dw_LMdist, ddw_LMdist, p_LMdist, xi_LMdist,&
       &             s_LMdist, ds_LMdist, z_LMdist, dz_LMdist, b_LMdist,          &
       &             db_LMdist, ddb_LMdist, aj_LMdist, dj_LMdist, ddj_LMdist,     &
       &             b_ic_LMdist, db_ic_LMdist, ddb_ic_LMdist, aj_ic_LMdist,     &
       &             dj_ic_LMdist, ddj_ic_LMdist, dp_LMdist, xi_LMdist,          &
       &             dxi_LMdist
   use fieldsLast, only: dwdt, dzdt, dpdt, dsdt, dbdt, djdt, dbdt_ic,  &
       &                 djdt_ic, dxidt, domega_ic_dt, domega_ma_dt,   &
       &                 lorentz_torque_ma_dt, lorentz_torque_ic_dt,   &
       &                 dwdt_dist, dzdt_dist, dpdt_dist, dsdt_dist,   &
       &                 dbdt_dist, djdt_dist, dbdt_ic_dist,  &
       &                 djdt_ic_dist, dxidt_dist
   use kinetic_energy, only: get_e_kin, get_u_square
   use magnetic_energy, only: get_e_mag
   use fields_average_mod, only: fields_average
   use spectra, only: spectrum, spectrum_temp, get_amplitude
   use outTO_mod, only: outTO
   use output_data, only: tag, l_max_cmb, n_coeff_r, l_max_r, n_coeff_r_max,&
       &                  n_r_array, n_r_step,  n_log_file, log_file
   use constants, only: vol_oc, vol_ic, mass, surf_cmb, two, three
   use outMisc_mod, only: outHelicity, outHeat
   use geos_mod, only: getEgeos
   use outRot, only: write_rot
   use omega, only: outOmega
   use integration, only: rInt_R
   use outPar_mod, only: outPar, outPerpPar
   use graphOut_mod, only: graphOut_IC
   use power, only: get_power
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC,  &
       &                     gather_from_lo_to_rank0
   use out_coeff, only: write_Bcmb, write_coeff_r, write_Pot
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
   use mpi_thetap_mod, only: transform_new2old !REMOVEEEEE
   use communications, only: gather_Flm !@> TODO: remove!!


   implicit none

   private

   !-- Counter for output files/sets:
   integer :: nPotSets, n_spec
   integer :: n_dt_cmb_sets, n_cmb_setsMov
   integer, allocatable :: n_v_r_sets(:), n_b_r_sets(:), n_T_r_sets(:)
   integer :: nTOsets,nTOmovSets,nTOrmsSets

   !--- For averaging:
   real(cp) :: timePassedLog, timeNormLog
   integer :: nLogs

   real(cp) :: dlBMean,dmBMean
   real(cp) :: lvDissMean,lbDissMean
   real(cp) :: RmMean,ElMean,ElCmbMean,RolMean
   real(cp) :: GeosMean,GeosAMean,GeosZMean,GeosMMean
   real(cp) :: GeosNAMean,GeosNAPMean
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
   integer, allocatable :: n_v_r_file(:)
   integer, allocatable :: n_t_r_file(:)
   integer, allocatable:: n_b_r_file(:)
   character(len=72) :: dtE_file, par_file
   character(len=72) :: cmb_file, dt_cmb_file, cmbMov_file
   character(len=72), allocatable :: v_r_file(:)
   character(len=72), allocatable :: t_r_file(:)
   character(len=72), allocatable :: b_r_file(:)

   public :: output, initialize_output, finalize_output

contains

   subroutine initialize_output

      integer :: n
      character(len=72) :: string

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

         do n=1,n_coeff_r_max
            write(string,*) n
            v_r_file(n)='V_coeff_r'//trim(adjustl(string))//'.'//tag
            if ( l_mag ) then
               B_r_file(n)='B_coeff_r'//trim(adjustl(string))//'.'//tag
            end if
         end do

         if ( l_r_fieldT ) then
            allocate ( n_t_r_file(n_coeff_r_max), t_r_file(n_coeff_r_max) )
            allocate ( n_t_r_sets(n_coeff_r_max) )
            n_T_r_sets=0

            do n=1,n_coeff_r_max
               write(string,*) n
               t_r_file(n)='T_coeff_r'//trim(adjustl(string))//'.'//tag
            end do
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
      GeosNAMean   =0.0_cp
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

      if ( l_master_rank .and. .not. l_save_out ) then
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

         if ( l_r_field ) then
            do n=1,n_coeff_r_max
               open(newunit=n_v_r_file(n), file=v_r_file(n), &
               &    status='new', form='unformatted')
               if ( l_mag ) then
                  open(newunit=n_b_r_file(n), file=b_r_file(n), &
                  &    status='new', form='unformatted')
               end if
            end do
         end if
         if ( l_r_fieldT ) then
            do n=1,n_coeff_r_max
               open(newunit=n_t_r_file(n), file=t_r_file(n), &
               &    status='new', form='unformatted')
            end do
         end if

      end if

   end subroutine initialize_output
!----------------------------------------------------------------------------
   subroutine finalize_output

      integer :: n

      if ( l_master_rank .and. .not. l_save_out ) then
         if ( l_mag .and. l_cmb_field ) then
            close(n_cmb_file)
            if (l_movie) close(n_cmbMov_file)
         end if
         if ( l_mag .and. l_dt_cmb_field ) then
            close(n_dt_cmb_file)
         end if
         if ( l_r_field ) then
            do n=1,n_coeff_r_max
               close(n_v_r_file(n))
               if ( l_mag ) then
                  close(n_b_r_file(n))
               end if
            end do
         end if

         if ( l_r_fieldT ) then
            do n=1,n_coeff_r_max
               close(n_t_r_file(n))
            end do
         end if

         if ( l_power ) close(n_dtE_file)
      end if

      if ( l_r_field .or. l_r_fieldT ) then
         deallocate ( n_coeff_r, n_v_r_file, v_r_file, n_v_r_sets )

         if ( l_mag ) then
            deallocate ( n_b_r_file, b_r_file, n_b_r_sets )
         end if

         if ( l_r_fieldT ) then
            deallocate ( n_t_r_file, t_r_file, n_t_r_sets )
         end if
      end if

   end subroutine finalize_output
!----------------------------------------------------------------------------
   subroutine output(time,tscheme,n_time_step,l_stop_time,l_pot,l_log,    &
              &      l_graph,lRmsCalc,l_store,l_new_rst_file,             &
              &      l_spectrum,lTOCalc,lTOframe,lTOZwrite,               &
              &      l_frame,n_frame,l_cmb,n_cmb_sets,l_r,                &
              &      lorentz_torque_ic,lorentz_torque_ma,dbdt_CMB_LMdist, &
              &      HelASr,Hel2ASr,HelnaASr,Helna2ASr,HelEAASr,viscASr,  &
              &      uhASr,duhASr,gradsASr,fconvASr,fkinASr,fviscASr,     &
              &      fpoynASr,fresASr,EperpASr,EparASr,EperpaxiASr,EparaxiASr)
      !
      !  This subroutine controls most of the output.
      !

      !--- Input of variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      logical,             intent(in) :: l_stop_time
      logical,             intent(in) :: l_pot
      logical,             intent(in) :: l_log, l_graph, lRmsCalc, l_store
      logical,             intent(in) :: l_new_rst_file, l_spectrum
      logical,             intent(in) :: lTOCalc,lTOframe
      logical,             intent(in) :: l_frame, l_cmb, l_r
      logical,             intent(inout) :: lTOZwrite
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
      real(cp),    intent(in) :: HelASr(2,nRstart:nRstop)
      real(cp),    intent(in) :: Hel2ASr(2,nRstart:nRstop)
      real(cp),    intent(in) :: HelnaASr(2,nRstart:nRstop)
      real(cp),    intent(in) :: Helna2ASr(2,nRstart:nRstop)
      real(cp),    intent(in) :: HelEAASr(nRstart:nRstop)
      real(cp),    intent(in) :: viscASr(nRstart:nRstop)
      real(cp),    intent(inout) :: uhASr(nRstart:nRstop)
      real(cp),    intent(inout) :: gradsASr(nRstart:nRstop)
      real(cp),    intent(inout) :: duhASr(nRstart:nRstop)
      real(cp),    intent(in) :: fconvASr(nRstart:nRstop)
      real(cp),    intent(in) :: fkinASr(nRstart:nRstop)
      real(cp),    intent(in) :: fviscASr(nRstart:nRstop)
      real(cp),    intent(in) :: fpoynASr(nRstartMag:nRstopMag)
      real(cp),    intent(in) :: fresASr(nRstartMag:nRstopMag)
      real(cp),    intent(inout) :: EperpASr(nRstart:nRstop)
      real(cp),    intent(inout) :: EparASr(nRstart:nRstop)
      real(cp),    intent(inout) :: EperpaxiASr(nRstart:nRstop)
      real(cp),    intent(inout) :: EparaxiASr(nRstart:nRstop)

      complex(cp), intent(in) :: dbdt_CMB_LMdist(n_mloMag_loc)

      !--- Local stuff:
      !--- Energies:
      real(cp) :: ekinR(n_r_max)     ! kinetic energy w radius
      real(cp) :: e_mag,e_mag_ic,e_mag_cmb,e_mag_p,e_mag_t
      real(cp) :: e_mag_p_as,e_mag_t_as,e_mag_p_ic,e_mag_t_ic
      real(cp) :: e_mag_p_as_ic,e_mag_t_as_ic,e_mag_os,e_mag_as_os
      real(cp) :: e_kin,e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as
      real(cp) :: eKinIC,eKinMA,dtE

      integer :: nR,lm,n,m

      !--- For TO:
      logical :: lTOrms

      !--- Property parameters:
      complex(cp) :: dbdtCMB(n_mloMag_loc)        ! SV at CMB !
      real(cp) :: volume,EC
      real(cp) :: dlVR(n_r_max),dlVRc(n_r_max)
      real(cp) :: RolRu2(n_r_max),RmR(n_r_max),dlPolPeakR(n_r_max)
      real(cp) :: Re,Ro,Rm,El,ElCmb,Rol,Geos,GeosA,GeosZ,GeosM,GeosNA,GeosNAP
      real(cp) :: Dip,DipCMB,ReConv,RoConv,e_kin_nas,RolC
      real(cp) :: elsAnel,dlVPolPeak,dlBPolPeak
      real(cp) :: dlB,dlBc,dmB,dLh,dlV,dlVc,dmV,dpV,dzV
      real(cp) :: visDiss,ohmDiss,lvDiss,lbDiss
      integer :: l
      real(cp) :: ReEquat
      real(cp) :: timeScaled
      character(len=96) :: message
      logical :: DEBUG_OUTPUT=.false.
      
      timeScaled=tScale*time
      timePassedLog=timePassedLog+tscheme%dt(1)

      !~~~~~~~~~~~~~~~~~~~~~~~ Conversion Loc > Dist ~~~~~~~~~~~~~~~~~~~~~~
      if ( (l_log .and. l_par) .or. l_pot .or. &
      &    l_frame .or. (l_SRIC .and. l_stop_time ) ) then
         call transform_new2old(w_LMdist, w_LMloc, n_r_max)
         call transform_new2old(dw_LMdist, dw_LMloc, n_r_max)
         call transform_new2old(ddw_LMdist, ddw_LMloc, n_r_max)
         call transform_new2old(p_LMdist, p_LMloc, n_r_max)
         call transform_new2old(dp_LMdist, dp_LMloc, n_r_max)
         call transform_new2old(z_LMdist, z_LMloc, n_r_max)
         call transform_new2old(dz_LMdist, dz_LMloc, n_r_max)

         if ( l_heat ) then
            call transform_new2old(s_LMdist, s_LMloc, n_r_max)
            call transform_new2old(ds_LMdist, ds_LMloc, n_r_max)
         end if
         if ( l_chemical_conv ) then
            call transform_new2old(xi_LMdist, xi_LMloc, n_r_max)
            call transform_new2old(dxi_LMdist, dxi_LMloc, n_r_max)
         end if
         if ( l_mag ) then
            call transform_new2old(b_LMdist, b_LMloc, n_r_max)
            call transform_new2old(db_LMdist, db_LMloc, n_r_max)
            call transform_new2old(ddb_LMdist, ddb_LMloc, n_r_max)
            call transform_new2old(aj_LMdist, aj_LMloc, n_r_max)
            call transform_new2old(dj_LMdist, dj_LMloc, n_r_max)
            call transform_new2old(ddj_LMdist, ddj_LMloc, n_r_max)
            if ( l_cond_ic ) then
               call transform_new2old(b_ic_LMdist, b_ic_LMloc, n_r_ic_max)
               call transform_new2old(db_ic_LMdist, db_ic_LMloc, n_r_ic_max)
               call transform_new2old(ddb_ic_LMdist, ddb_ic_LMloc, n_r_ic_max)
               call transform_new2old(aj_ic_LMdist, aj_ic_LMloc, n_r_ic_max)
               call transform_new2old(dj_ic_LMdist, dj_ic_LMloc, n_r_ic_max)
               call transform_new2old(ddj_ic_LMdist, ddj_ic_LMloc, n_r_ic_max)
            end if
         end if
      end if

      if ( l_store .or. l_pot ) then
         do nR=nRstart,nRstop
            call gather_Flm(w_Rdist(:,nR), w_Rloc(:,nR))
            call gather_Flm(z_Rdist(:,nR), z_Rloc(:,nR))
            call gather_Flm(p_Rdist(:,nR), p_Rloc(:,nR))
            if ( l_heat ) call gather_Flm(s_Rdist(:,nR), s_Rloc(:,nR))
            if ( l_chemical_conv ) call gather_Flm(xi_Rdist(:,nR), xi_Rloc(:,nR))
            if ( l_mag ) then
               call gather_Flm(b_Rdist(:,nR), b_Rloc(:,nR))
               call gather_Flm(aj_Rdist(:,nR), aj_Rloc(:,nR))
            end if
         end do
         call dwdt_dist%gather_all(dwdt)
         call dpdt_dist%gather_all(dpdt)
         call dzdt_dist%gather_all(dzdt)
         if ( l_heat ) call dsdt_dist%gather_all(dsdt)
         if ( l_chemical_conv ) call dxidt_dist%gather_all(dxidt)
         if ( l_mag ) then
            call dbdt_dist%gather_all(dbdt)
            call djdt_dist%gather_all(djdt)
            if ( l_cond_ic ) then
               call dbdt_ic_dist%gather_all(dbdt_ic)
               call djdt_ic_dist%gather_all(djdt_ic)
            end if
         end if
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      ! We start with the computation of the energies
      ! in parallel.
      if ( l_log ) then

         nLogs=nLogs+1
         timeNormLog=timeNormLog+timePassedLog

         !----- Write torques and rotation rates:
         PERFON('out_rot')
         call write_rot( time,tscheme%dt(1),eKinIC,eKinMA,w_LMdist,z_LMdist, &
              &          dz_LMdist,b_LMdist,omega_ic,omega_ma,               &
              &          lorentz_torque_ic,lorentz_torque_ma)
         PERFOFF
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  write_rot  on coord_r ",coord_r

         PERFON('out_ekin')
         n_e_sets=n_e_sets+1
         call get_e_kin(time,.true.,l_stop_time,n_e_sets,w_LMdist,     &
              &         dw_LMdist,z_LMdist,e_kin_p,e_kin_t,e_kin_p_as, &
              &         e_kin_t_as,ekinR)
         e_kin=e_kin_p+e_kin_t
         e_kin_nas=e_kin-e_kin_p_as-e_kin_t_as
         if ( DEBUG_OUTPUT ) write(*,"(A,I6)") "Written  e_kin  on coord_r ",coord_r


         call get_e_mag(time,.true.,l_stop_time,n_e_sets,b_LMdist,db_LMdist, &
              &         aj_LMdist,b_ic_LMdist,db_ic_LMdist,aj_ic_LMdist,     &
              &         e_mag_p,e_mag_t,e_mag_p_as,e_mag_t_as,e_mag_p_ic,    &
              &         e_mag_t_ic,e_mag_p_as_ic,e_mag_t_as_ic,              &
              &         e_mag_os,e_mag_as_os,e_mag_cmb,Dip,DipCMB,elsAnel )
         e_mag   =e_mag_p+e_mag_t
         e_mag_ic=e_mag_p_ic+e_mag_t_ic
         PERFOFF
         if ( DEBUG_OUTPUT ) write(*,"(A,I6)") "Written  e_mag  on coord_r ",coord_r

         !----- Calculate distribution of energies on all m's
         if ( l_energy_modes ) then
            call get_amplitude(time,w_LMdist,dw_LMdist,z_LMdist,b_LMdist,&
                 &             db_LMdist,aj_LMdist)
            if ( DEBUG_OUTPUT ) &
               & write(*,"(A,I6)") "Written  amplitude  on coord_r ",coord_r
         endif

         !---- Surface zonal velocity at the equator
         if ( ktopv==1 ) then
            ReEquat=0.0_cp
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               m = map_mlo%i2m(lm)
               if ( m == 0 ) then
                  ReEquat=ReEquat-real(z_LMdist(lm,n_r_cmb))*dPl0Eq(l+1)*or1(n_r_cmb)
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
            call spectrum(n_spec,time,.true.,nLogs,l_stop_time,timePassedLog,    &
                 &        timeNormLog,w_LMdist,dw_LMdist,z_LMdist,b_LMdist,      &
                 &        db_LMdist,aj_LMdist,b_ic_LMdist,db_ic_LMdist,aj_ic_LMdist)
            if ( l_heat ) then
               call spectrum_temp(n_spec,time,.true.,nLogs,l_stop_time,     &
                    &             timePassedLog,timeNormLog,s_LMdist,ds_LMdist)
            end if
         end if

         if ( l_average ) then
            call fields_average(time,tscheme,nLogs,l_stop_time,timePassedLog,    &
                 &              timeNormLog,omega_ic,omega_ma,w_LMdist,z_LMdist, &
                 &              p_LMdist,s_LMdist,xi_LMdist,b_LMdist,aj_LMdist,  &
                 &              b_ic_LMdist,aj_ic_LMdist)
            if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  averages  on coord_r ",coord_r
         end if

         if ( l_power ) then

            PERFON('out_pwr')
            if ( l_master_rank ) then
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

            call get_power( time,timePassedLog,timeNormLog,l_stop_time,          &
                 &          omega_ic,omega_ma,lorentz_torque_ic,                 &
                 &          lorentz_torque_ma,w_LMdist,z_LMdist,                 &
                 &          dz_LMdist,s_LMdist,xi_LMdist,                        &
                 &          b_LMdist,ddb_LMdist,aj_LMdist,dj_LMdist,db_ic_LMdist,&
                 &          ddb_ic_LMdist,aj_ic_LMdist,dj_ic_LMdist,viscASr,     &
                 &          visDiss,ohmDiss)
            PERFOFF
            if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  power  on coord_r ",coord_r
         end if

         !----- If anelastic additional u**2 outputs
         if ( l_anel ) then
            call get_u_square(time,w_LMdist,dw_LMdist,z_LMdist,RolRu2)
         else
            RolRu2(:)=0.0_cp
         end if

         !-- Get flow lengthscales
         call getDlm(w_LMdist,dw_LMdist,z_LMdist,dlV,dlVR,dmV,dlVc,dlVPolPeak, &
              &      dlVRc,dlPolPeakR,'V')

         !-- Out radial profiles of parameters
         call outPar(timePassedLog,timeNormLog,l_stop_time,s_LMdist,ds_LMdist,&
              &      p_LMdist,dp_LMdist,ekinR,RolRu2,dlVR,dlVRc,dlPolPeakR,   &
              &      uhASr,duhASr,gradsASr,fconvASr,fkinASr,fviscASr,fpoynASr,&
              &      fresASr,RmR)

         !-- Perpendicular/parallel
         if ( l_perpPar ) then
            call outPerpPar(time,timePassedLog,timeNormLog,l_stop_time, &
                 &          EperpASr,EparASr,EperpaxiASr,EparaxiASr)
         end if

         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  outPar  on coord_r ",coord_r

         if ( l_heat .or. l_chemical_conv ) then
            call outHeat(timeScaled,timePassedLog,timeNormLog,l_stop_time, &
                 &       s_LMdist,ds_LMdist,p_LMdist,dp_LMdist,xi_LMdist,  &
                 &       dxi_LMdist)
         end if

         if ( l_hel ) then
            call outHelicity(timeScaled,HelASr,Hel2ASr,HelnaASr,Helna2ASr,HelEAASr)
         end if
         
         if ( l_par ) then
            call getEgeos(timeScaled,nLogs,w_LMloc,dw_LMloc,ddw_LMloc,    &
                 &        z_LMloc,dz_LMloc,Geos,GeosA,GeosZ,GeosM,GeosNA, & 
                 &        GeosNAP,dpV,dzV,volume,EC)
         else
            Geos   =0.0_cp
            GeosA  =0.0_cp
            GeosZ  =0.0_cp
            GeosM  =0.0_cp
            GeosNA =0.0_cp
            GeosNAP=0.0_cp
            dpV    =0.0_cp
            dzV    =0.0_cp
            volume =0.0_cp ! test volume
            EC     =0.0_cp ! test kinetic energy
         end if

         if ( l_mag .or. l_mag_LF ) then
            !-- Get magnetic field lengthscales
            call getDlm(b_LMdist,db_LMdist,aj_LMdist,dlB,dlVR,dmB, &
                 &      dlBc,dlBPolPeak,dlVRc,dlPolPeakR,'B')
         else
            dlB=0.0_cp
            dmB=0.0_cp
         end if
      end if

      if ( l_spectrum ) then
         n_spec=n_spec+1
         call spectrum(n_spec,time,.false.,nLogs,l_stop_time,timePassedLog, &
              &        timeNormLog,w_LMdist,dw_LMdist,z_LMdist,b_LMdist,    &
              &        db_LMdist,aj_LMdist,b_ic_LMdist,db_ic_LMdist,aj_ic_LMdist)
         if ( l_heat ) then
            call spectrum_temp(n_spec,time,.false.,nLogs,l_stop_time,     &
                 &             timePassedLog,timeNormLog,s_LMdist,ds_LMdist)
         end if
         if ( l_master_rank ) then
            write(*,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)')&
            &    " ! Storing spectra:",                  &
            &    "             at time=",timeScaled,     &
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
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  spectrum  on coord_r ",coord_r
      end if

      if ( lTOCalc ) then
         !------ Output for every log time step:
         if ( lVerbose ) write(*,*) '! Calling outTO !'
         lTOrms   =.true.
         if ( .not. l_log ) then
            call get_e_kin(time,.false.,l_stop_time,0,w_LMdist,dw_LMdist, &
                 &         z_LMdist,e_kin_p,e_kin_t,e_kin_p_as,e_kin_t_as,&
                 &         ekinR)
            e_kin=e_kin_p+e_kin_t
         end if
         call outTO(time,n_time_step,e_kin,e_kin_t_as,nTOsets,nTOmovSets, &
         &          nTOrmsSets,lTOframe,lTOrms,lTOZwrite,z_LMdist,omega_ic,&
         &          omega_ma)
         !------ Note: time averaging, time differencing done by IDL routine!

         if ( lVerbose ) write(*,*) '! outTO finished !'
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  TO  on coord_r ",coord_r
      end if

      !--- Get radial derivatives and add dt dtB terms:
      if ( l_dtB ) then
         call get_dtBLMfinish(time,n_time_step,omega_ic,b_LMdist,ddb_LMdist, &
              &               aj_LMdist,dj_LMdist,ddj_LMdist,b_ic_LMdist,    &
              &               db_ic_LMdist,ddb_ic_LMdist,aj_ic_LMdist,       &
              &               dj_ic_LMdist,ddj_ic_LMdist,l_frame)
      end if

      if ( l_RMS ) then
         if ( n_time_step == 1 ) then
            nRMS_sets    =0
            timeNormRMS  =0.0_cp
            timePassedRMS=0.0_cp
            call zeroRms
         end if
         timePassedRMS=timePassedRMS+tscheme%dt(1)
         if ( lRmsCalc ) then
            if ( lVerbose ) write(*,*) '! Writing RMS output !'
            timeNormRMS=timeNormRMS+timePassedRMS
            call dtVrms(time,nRMS_sets,timePassedRMS,timeNormRMS,l_stop_time)
            if ( l_mag ) call dtBrms(time)
            timePassedRMS=0.0_cp
         end if
         if (DEBUG_OUTPUT) write(*,"(A,I6)") "Written  dtV/Brms  on coord_r ",coord_r
      end if

      !--- Store poloidal magnetic coeffs at cmb
      if ( l_cmb ) then
         PERFON('out_cmb')
         call write_Bcmb(timeScaled,b_LMdist(:,n_r_cmb),l_max_cmb,n_cmb_sets,   &
              &          cmb_file,n_cmb_file)

         !--- Store SV of poloidal magnetic coeffs at cmb
         if ( l_dt_cmb_field ) then
            do lm=1,n_mloMag_loc
               l=map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               m=map_mlo%i2m(lm)
               dLh = real(l*(l+1),cp)
               dbdtCMB(lm)= dbdt_CMB_LMdist(lm)/(dLh*or2(n_r_cmb))    &
               &         + opm*hdif_B(l) * ( ddb_LMdist(lm,n_r_cmb) - &
               &           dLh*or2(n_r_cmb)*   b_LMdist(lm,n_r_cmb) )
            end do

            call write_Bcmb(timeScaled,dbdtCMB(:),l_max_cmb,n_dt_cmb_sets,  &
                 &          dt_cmb_file,n_dt_cmb_file)
         end if
         PERFOFF
      end if

      if ( l_frame .and. l_cmb_field ) then
         call write_Bcmb(timeScaled,b_LMdist(:,n_r_cmb),l_max_cmb,   &
              &          n_cmb_setsMov,cmbMov_file,n_cmbMov_file)
      end if

      !--- Store potential coeffs for velocity fields and magnetic fields
      if ( l_r ) then
         PERFON('out_r')
         do n=1,n_coeff_r_max
            nR=n_coeff_r(n)
            call write_coeff_r(timeScaled,w_LMdist(:,nR),dw_LMdist(:,nR),  &
                 &             ddw_LMdist(:,nR),z_LMdist(:,nR),r(nR),      &
                 &             l_max_r,n_v_r_sets(n),v_r_file(n),          &
                 &             n_v_r_file(n),1)
            if ( l_mag )                                                     &
               call write_coeff_r(timeScaled,b_LMdist(:,nR),db_LMdist(:,nR), &
                    &             ddb_LMdist(:,nR),aj_LMdist(:,nR),r(nR),    &
                    &             l_max_r,n_b_r_sets(n),b_r_file(n),         &
                    &             n_b_r_file(n),2)
            if ( l_r_fieldT )                                                &
               call write_coeff_r(timeScaled,s_LMdist(:,nR),db_LMdist(:,nR), &
                    &             ddb_LMdist(:,nR),aj_LMdist(:,nR),r(nR),    &
                    &             l_max_r,n_T_r_sets(n),T_r_file(n),         &
                    &             n_t_r_file(n),3)
         end do
         PERFOFF
      end if

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
            call rBrSpec(time,b_LMdist, b_ic_LMdist ,'rBrSpecMov',.true.)
            call rBpSpec(time,aj_LMdist,aj_ic_LMdist,'rBpSpecMov',.true.)
         end if
         if ( l_log ) then
            call rBrSpec(time,b_LMdist, b_ic_LMdist ,'rBrSpec',.true.)
            call rBpSpec(time,aj_LMdist,aj_ic_LMdist,'rBpSpec',.true.)
         end if
      end if

      !
      ! Writing of the restart file
      !
      if ( l_store ) then
#ifdef WITH_MPI
         call store_mpi(time,tscheme,n_time_step,l_stop_time,l_new_rst_file,  &
              &         .false.,w_Rloc,z_Rloc,p_Rloc,s_Rloc,xi_Rloc,b_Rloc,   &
              &         aj_Rloc,b_ic_LMloc,aj_ic_LMloc,dwdt,dzdt,dpdt,dsdt,   &
              &         dxidt,dbdt,djdt,dbdt_ic,djdt_ic,domega_ma_dt,         &
              &         domega_ic_dt,lorentz_torque_ma_dt,lorentz_torque_ic_dt)
#else
         call store(time,tscheme,n_time_step,l_stop_time,l_new_rst_file,.false., &
              &     w_LMloc,z_LMloc,p_LMloc,s_LMloc,xi_LMloc,b_LMloc,aj_LMloc,   &
              &     b_ic_LMloc,aj_ic_LMloc,dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,       &
              &     djdt,dbdt_ic,djdt_ic,domega_ma_dt,domega_ic_dt,              &
              &     lorentz_torque_ma_dt,lorentz_torque_ic_dt)
#endif
      end if


      ! ===================================================
      !      GATHERING for output
      ! ===================================================
      ! We have all fields in LMloc space. Thus we gather the whole fields on coord_r 0.

      if ( l_frame .or. (l_graph .and. l_mag .and. n_r_ic_max > 0) ) then
         PERFON('out_comm')

         if ( l_mag ) call gather_from_lo_to_rank0(b_LMloc(:,n_r_icb),bICB)

         if ( l_cond_ic ) then
            call gather_all_from_lo_to_rank0(gt_IC,b_ic_LMloc,b_ic)
            call gather_all_from_lo_to_rank0(gt_IC,db_ic_LMloc,db_ic)
            call gather_all_from_lo_to_rank0(gt_IC,ddb_ic_LMloc,ddb_ic)

            call gather_all_from_lo_to_rank0(gt_IC,aj_ic_LMloc,aj_ic)
            call gather_all_from_lo_to_rank0(gt_IC,dj_ic_LMloc,dj_ic)
            call gather_all_from_lo_to_rank0(gt_IC,ddj_ic_LMloc,ddj_ic)
         end if

         ! for writing a restart file, we also need the d?dtLast arrays,
         ! which first have to be gathered on coord_r 0
         PERFOFF

      end if

      !--- Movie output and various supplementary things:
      if ( l_frame ) then
         ! The frames array for the movies is distributed over the ranks
         ! and has to be gathered on coord_r 0 for output.

         ! Each movie uses some consecutive frames in the frames array. They
         ! start at n_movie_field_start(1,n_movie)
         ! up to    n_movie_field_stop(1+n_fields_oc+n_fields,n_movie) (n_fields_ic>0
         ! or       n_movie_field_stop(1+n_fields,n_movie)             (n_fields_ic=0)

         call movie_gather_frames_to_rank0()

         if ( l_movie_ic .and. l_store_frame .and. l_master_rank ) then
            call store_movie_frame_IC(bICB,b_ic,db_ic,ddb_ic,aj_ic,dj_ic)
         end if

         n_frame=n_frame+1
         call logWrite(' ')
         if ( l_master_rank ) then 
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

      ! =======================================================================
      ! ======= compute output on coord_r 0 ==============
      ! =======================================================================
      if ( l_master_rank ) then
         PERFON('out_out')

         !----- Plot out inner core magnetic field, outer core
         !      field has been written in radialLoop !
         if ( l_graph .and. l_mag .and. n_r_ic_max > 0 ) then
            call graphOut_IC(b_ic,db_ic,ddb_ic,aj_ic,dj_ic,bICB)
         end if

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
            !write(*,"(A,3ES20.12)") "dlVc,RoConv,RolC = ",dlVc,RoConv,RolC

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
            GeosNAMean =GeosNAMean +timePassedLog*GeosNA
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
               GeosNAMean =GeosNAMean/timeNormLog 
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
                  write(*,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6)')              &
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
                  write(*,'(1p,/,A,/,A,/,A,4ES16.6,/,A,4ES16.6,/,A,4ES16.6)')  &
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
               & '(1p,/,A,16(/,A,ES12.4),/,A,4ES12.4,/,A,2ES12.4,/,A,2ES12.4)') &
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
               & " !  rel geos NA Ekin :",GeosNAMean,                           &
               & " !  rel geos NAP Ekin:",GeosNAPMean,                          &
               & " !  Dip              :",DipMean,                              &
               & " !  DipCMB           :",DipCMBMean,                           &
               & " !  l,m,p,z V scales :",dlVMean,dmVMean,dpVMean,dzVmean,      &
               & " !  l,m, B scales    :",dlBMean,dmBMean,                      &
               & " !  vis, Ohm scale   :",lvDissMean,lbDissMean
               if ( l_par ) then 
                  write(n_log_file,*) !' Calculating geostrophic contributions with outEgeos.f90'
                  write(n_log_file,*) '! precision of z-integration (geos):',abs(volume/vol_oc-1)
                  write(n_log_file,*) '! precision of cyl. transf.  (geos):',abs(EC/e_kin-1)
               end if

               if ( l_save_out ) close(n_log_file)

            end if ! l_stop_time ?

         end if ! l_log



         PERFOFF
      end if

      if ( l_SRIC .and. l_stop_time ) call outOmega(z_LMdist,omega_ic)

      if ( l_log ) then
         timePassedLog=0.0_cp
      end if

      if ( lRmsCalc ) then
         call zeroRms
      end if

   end subroutine output
!----------------------------------------------------------------------------
end module output_mod
