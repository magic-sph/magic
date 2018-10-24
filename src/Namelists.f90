module Namelists
   !
   ! Read and print the input namelist
   !

   use precision_mod
   use constants
   use truncation
   use radial_functions
   use physical_parameters
   use num_param
   use torsional_oscillations
   use init_fields
   use logic
   use output_data
   use parallel_mod
   use special
   use movie_data, only: movie,n_movies, n_movies_max
   use charmanip, only: length_to_blank,capitalize
   use blocking, only: cacheblock_size_in_B
   use probe_mod
   use useful, only: abortRun

   implicit none
 
   private
 
   public :: readNamelists, writeNamelists

contains

   subroutine readNamelists
      !
      !                                                                   
      !  Purpose of this subroutine is to read the input namelists.       
      !  This program also determins logical parameters that are stored   
      !  in logic.f90.                                                    
      !

      !-- Local stuff
      integer :: n
      integer :: nCounts
      real(cp) :: xiCMB(4*n_impXi_max)
      real(cp) :: sCMB(4*n_impS_max),rad ! cmb heat boundary condition
      logical :: log_does_exist, nml_exist
      integer :: length
      integer :: argument_count
      integer :: res
      integer :: inputHandle
      character(len=100) :: input_filename

      !-- Name lists:
      integer :: runHours,runMinutes,runSeconds
      namelist/grid/n_r_max,n_cheb_max,n_phi_tot,n_theta_axi, &
         &  n_r_ic_max,n_cheb_ic_max,minc,nalias,l_axi,       &
         &  fd_order,fd_order_bound,fd_ratio,fd_stretch

      namelist/control/                                     &
         & mode,tag,n_time_steps,                           &
         & n_tScale,n_lScale,alpha,enscale,                 &
         & l_update_v,l_update_b,l_update_s,l_update_xi,    &
         & dtstart,dtMax,courfac,alffac,intfac,n_cour_step, &
         & difnu,difeta,difkap,difchem,ldif,ldifexp,        &
         & l_correct_AMe,l_correct_AMz,tEND,l_non_rot,      &
         & l_newmap,alph1,alph2,l_cour_alf_damp,            &
         & runHours,runMinutes,runSeconds,map_function,     &
         & cacheblock_size_in_B,anelastic_flavour,          &
         & thermo_variable,radial_scheme,polo_flow_eq
      
      namelist/phys_param/                                      &
         & ra,raxi,pr,sc,prmag,ek,epsc0,epscxi0,radratio,       &
         & ktops,kbots,ktopv,kbotv,ktopb,kbotb,kbotxi,ktopxi,   &
         & s_top,s_bot,impS,sCMB,xi_top,xi_bot,impXi,xiCMB,     &
         & nVarCond,con_DecRate,con_RadRatio,con_LambdaMatch,   &
         & con_LambdaOut,con_FuncWidth,ThExpNb,GrunNb,          &
         & strat,polind,DissNb,g0,g1,g2,r_cut_model,thickStrat, &
         & epsS,slopeStrat,rStrat,ampStrat,cmbHflux,r_LCR,      &
         & nVarDiff,nVarVisc,difExp,nVarEps,interior_model,     &
         & nVarEntropyGrad,l_isothermal,ktopp,po,prec_angle,    &
         & po_diff,diff_prec_angle

      namelist/B_external/                                    &
         & rrMP,amp_imp,expo_imp,bmax_imp,n_imp,l_imp,        &
         & l_curr,amp_curr

      namelist/start_field/                                   &
         & l_start_file,start_file,inform,l_reset_t,          &
         & scale_s,scale_xi,scale_b,scale_v,tipdipole,        &
         & init_s1,init_s2,init_v1,init_b1,imagcon,tmagcon,   &
         & amp_s1,amp_s2,amp_v1,amp_b1, init_xi1, init_xi2,   &
         & amp_xi1, amp_xi2

      namelist/output_control/                                &
         & n_graph_step,n_graphs,t_graph,                     &
         & t_graph_start,t_graph_stop,dt_graph,l_graph_time,  &
         & n_stores,n_rst_step,n_rsts,t_rst,                  &
         & t_rst_start,t_rst_stop,dt_rst,                     &
         & n_log_step,n_logs,t_log,t_log_start,t_log_stop,    &
         & dt_log,n_spec_step,n_specs,t_spec,t_spec_start,    &
         & t_spec_stop,dt_spec,n_cmb_step,n_cmbs,t_cmb,       &
         & t_cmb_start,t_cmb_stop,dt_cmb,                     &
         & n_r_field_step,n_r_fields,t_r_field,               &
         & t_r_field_start,t_r_field_stop,dt_r_field,         &
         & n_Bpot_step,n_Bpots,t_Bpot,t_Bpot_start,           &
         & t_Bpot_stop,dt_Bpot,n_Vpot_step,n_Vpots,t_Vpot,    &
         & t_Vpot_start,t_Vpot_stop,dt_Vpot,n_Tpot_step,      &
         & n_Tpots,t_Tpot,t_Tpot_start,t_Tpot_stop,dt_Tpot,   &
         & n_pot_step,n_pots,t_pot,t_pot_start,t_pot_stop,    &
         & dt_pot,runid,movie,n_movie_step,                   &
         & n_movie_frames,t_movie,t_movie_start,t_movie_stop, &
         & dt_movie,n_TO_step,n_TOs,t_TO,t_TO_start,t_TO_stop,&
         & dt_TO,n_TOZ_step,n_TOZs,t_TOZ,t_TOZ_start,         &
         & t_TOZ_stop,dt_TOZ,n_TOmovie_step,n_TOmovie_frames, &
         & t_TOmovie,t_TOmovie_start,t_TOmovie_stop,          &
         & dt_TOmovie,l_movie,l_average,l_save_out,           &
         & l_true_time,l_cmb_field,l_rMagSpec,l_DTrMagSpec,   &
         & l_dt_cmb_field,l_max_cmb,l_r_field,l_r_fieldT,     &
         & n_r_step,l_max_r,n_r_array,l_TO,l_TOmovie,l_hel,   &
         & lVerbose,l_AM,l_power,l_drift,l_storeBpot,         &
         & l_storeVpot,l_storeTpot,l_storePot,sDens,zDens,    &
         & l_RMS,l_par,l_corrMov,rCut,rDea,                   &
         & l_PV,l_iner,l_viscBcCalc,l_fluxProfs,l_perpPar,    &
         & l_PressGraph,l_energy_modes,m_max_modes,l_probe,   &
         & r_probe,theta_probe,n_phi_probes,n_probe_step,     &
         & n_probe_out,t_probe_start,t_probe_stop,dt_probe,   &
         & l_earth_likeness,l_max_comp

      namelist/mantle/conductance_ma,nRotMa,rho_ratio_ma, &
         & omega_ma1,omegaOsz_ma1,tShift_ma1,             &
         & omega_ma2,omegaOsz_ma2,tShift_ma2,             &
         & amp_RiMaAsym,omega_RiMaAsym,m_RiMaAsym,        &
         & amp_RiMaSym,omega_RiMaSym,m_RiMaSym

      namelist/inner_core/sigma_ratio,nRotIc,rho_ratio_ic, &
         & omega_ic1,omegaOsz_ic1,tShift_ic1,              &
         & omega_ic2,omegaOsz_ic2,tShift_ic2,BIC,          &
         & amp_RiIcAsym,omega_RiIcAsym,m_RiIcAsym,         &
         & amp_RiIcSym,omega_RiIcSym,m_RiIcSym


      do n=1,4*n_impS_max
         sCMB(n)=0.0_cp
      end do

      do n=1,4*n_impXi_max
         xiCMB(n)=0.0_cp
      end do

      runHours  =0
      runMinutes=0
      runSeconds=0

      !-- Set default values of control parameters:
      call defaultNamelists


      ! get the filename of the input file as first argument from the command line
      argument_count = command_argument_count()
      if (argument_count == 0) then
         call abortRun('The filename of the input file must be provided as first argument')
      else
         call get_command_argument(1,input_filename)

         inquire(file = input_filename, exist = nml_exist)

         if (.not. nml_exist) then
            call abortRun('! Input namelist file not found!')
         end if

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading control parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading grid parameters!'
         read(inputHandle,nml=grid,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No grid namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading control parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading control parameters!'
         read(inputHandle,nml=control,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No control namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading physical parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading physical parameters!'
         read(inputHandle,nml=phys_param,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No phys_param namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading start information!'
         read(inputHandle,nml=start_field,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No start_field namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading output parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading output information!'
         read(inputHandle,nml=output_control,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No output_control namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading inner-core parameter from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading inner core information!'
         read(inputHandle,nml=inner_core,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No inner_core namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading mantle parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading mantle information!'
         read(inputHandle,nml=mantle,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No mantle namelist found!'
         end if
         close(inputHandle)

         open(newunit=inputHandle,file=trim(input_filename))
         !-- Reading external field parameters for feedback:
         if ( rank == 0 ) write(*,*) '!  Reading B external parameters!'
         read(inputHandle,nml=B_external,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No B_external namelist found!'
         end if
         close(inputHandle)

         !-- Correcting some parameters:
      end if

      !-- Does log-file already exist?
      log_file='log.'//tag
      if ( rank == 0 ) then
         inquire(file=trim(log_file),exist=log_does_exist)
      end if
#ifdef WITH_MPI
      call MPI_Bcast(log_does_exist,1,MPI_logical,0,MPI_COMM_WORLD,ierr)
#endif
      if (log_does_exist) then
         if ( rank == 0 ) then 
            write(*,*)
            write(*,*) '! The log-file exists already !'
            write(*,*) '! I add _BIS to the tag and create new files!'
         end if
         length=length_to_blank(tag)
         tag=tag(1:length)//'_BIS'
      end if

      call capitalize(polo_flow_eq)
      if ( index(polo_flow_eq, 'DC') /= 0 ) then
         l_double_curl = .true.
         l_PressGraph  = .false.
      else
         l_double_curl = .false.
      end if

      call capitalize(radial_scheme)
      if ( index(radial_scheme, 'FD') /= 0 ) then
         l_finite_diff = .true.
      else
         l_finite_diff = .false.
      end if

      if ( l_finite_diff ) then
         l_double_curl=.true.
         l_PressGraph =.false.
         write(*,*) '! Finite differences are used: I use the double-curl form !'
      end if

      n_stores=max(n_stores,n_rsts)

      l_TOmovie=l_TOmovie.and.l_TO

      !-- Determine what has to be calculated depending on mode:
      lMagMem  =1
      ldtBMem  =0
      lStressMem=0
      lMovieMem=0
      l_conv   =.true.
      l_conv_nl=.true.
      l_mag    =.true.
      l_mag_nl =.true.
      l_mag_kin=.false.
      l_mag_LF =.true.
      l_heat   =.true.
      l_heat_nl=.true.
      l_SRIC   =.false.
      l_SRMA   =.false.
      l_Ri     =.false.
      l_AB1    =.false.

      if ( mode == 1 ) then
         !-- Only convection:
         l_mag   =.false.
         l_mag_nl=.false.
         l_mag_LF=.false.
      else if ( mode == 2 ) then
         !-- Kinematic dynamo:
         l_conv=.false.
         l_conv_nl=.false.
         l_mag_kin=.true.
         l_mag_LF =.false.
         l_heat   =.false.
         l_heat_nl=.false.
      else if ( mode == 3 ) then
         !-- Magnetic decay modes:
         l_conv   =.false.
         l_conv_nl=.false.
         l_mag_nl =.false.
         l_mag_LF =.false.
         l_heat   =.false.
         l_heat_nl=.false.
      else if ( mode == 4 ) then
         !-- Magneto convection:
         l_mag    =.false.
         l_mag_nl =.false.
         l_mag_LF =.false.
         l_rot_ic =.false.
         l_rot_ma =.false.
      else if ( mode == 5 ) then
         !----- Onset of convection (linear):
         l_conv_nl=.false.
         l_mag    =.false.
         l_mag_nl =.false.
         l_mag_LF =.false.
         l_rot_ic =.false.
         l_rot_ma =.false.
      else if ( mode == 6 ) then
         !-- Self-consistent dynamo, but no Lorentz Force
         l_mag_LF=.false.
      else if ( mode == 7 ) then
         !-- Super-rotating IC or MA, no convection, no magnetic field
         l_heat   =.false.
         l_mag    =.false.
         l_mag_nl =.false.
         l_mag_LF =.false.
      else if ( mode == 8 ) then
         !-- Super-rotating IC or MA, no convection, plus magnetic field
         l_heat   =.false.
      else if ( mode == 9 ) then
         !-- Super-rotating IC or MA, no convection, magnetic field, no LF
         l_heat   =.false.
         l_mag_LF =.false.
      else if ( mode == 10 ) then
         !-- Super-rotating IC or MA, no convection, no magnetic field, no LF, 
         !-- no advection
         l_heat   =.false.
         l_conv_nl=.false.
         l_mag    =.false.
         l_mag_nl =.false.
         l_mag_LF =.false.
      end if

      if ( mode == 7 .or. mode == 8 .or. mode == 9 .or. mode == 10 ) then
         !kbotv=2
         !ktopv=2
         l_correct_AMz=.false.
         l_heat       =.false.
      end if

      if ( nRotIc > 0 ) then
         l_rot_ic=.true.
      else if ( nRotIc == 0 ) then
         l_rot_ic=.false.
      else if ( nRotIc == -1 ) then
         l_rot_ic=.true.
         l_SRIC  =.true.
      end if

      if ( nRotMa > 1 ) then
         l_rot_ma=.true.
      else if ( nRotMa == 0 ) then
         l_rot_ma=.false.
      else if ( nRotMa == -1 ) then
         l_rot_ma=.true.
         l_SRMA  =.true.
      end if

      !-- Inertial mode forcing at boundaries

      if ( amp_RiIcAsym /= 0.0_cp .or. amp_RiIcSym /= 0.0_cp .or. &
           amp_RiMaAsym /= 0.0_cp .or. amp_RiMaSym /= 0.0_cp) then
         l_Ri   = .true.
      end if

      if ( raxi > 0.0_cp ) then
         l_chemical_conv = .true.
      else
         l_chemical_conv = .false.
      end if

      if ( ra == 0.0_cp ) l_heat=.false.

      if ( ek < 0.0_cp ) l_non_rot= .true. 
      if ( l_non_rot ) then
         l_corr=.false.
         ek=-one ! used as a flag, not used for the calculation
      else
         l_corr=.true.
      end if

      !-- Choose between temperature and entropy (same in the Boussinesq limit)
      call capitalize(thermo_variable)
      if ( index(thermo_variable, 'T') /= 0 ) then
         l_TP_form=.true.
      else if ( index(thermo_variable, 'S') /= 0 .or. &
              & index(thermo_variable, 'ENT') /=0 ) then
         l_TP_form=.false.
      else
         l_TP_form=.false.
      end if

      !-- Choose between entropy diffusion and temperature diffusion 
      call capitalize(anelastic_flavour)
      if ( index(anelastic_flavour, 'LBR') /= 0 .or. &
         & index(anelastic_flavour, 'ENT') /= 0 ) then
         l_temperature_diff = .false.
         l_anelastic_liquid = .false.
         l_single_matrix    = .false.
      else if ( index(anelastic_flavour, 'ALA') /= 0 .or. &
           index(anelastic_flavour, 'LIQ') /= 0 ) then
         l_temperature_diff = .false.
         l_anelastic_liquid = .true.
         l_single_matrix    = .false.
      else if ( index(anelastic_flavour, 'TEMP') /= 0 .or. &
           index(anelastic_flavour, 'TDIFF') /= 0 ) then
         l_temperature_diff = .true.
         l_anelastic_liquid = .false.
         l_single_matrix    = .true.
      else if ( index(anelastic_flavour, 'SINGLEMAT') /= 0 ) then ! Testing purposes
         l_temperature_diff = .false.
         l_anelastic_liquid = .false.
         l_single_matrix    = .true.
      else
         l_temperature_diff = .false.
         l_anelastic_liquid = .false.
         l_single_matrix    = .false.
      end if

      if ( ktops > 2 .or. kbots > 2 .or. l_TP_form ) then
         l_single_matrix    = .true.
      end if

      if ( l_anelastic_liquid .or. l_temperature_diff .or. l_TP_form ) l_anel=.true.

      call capitalize(interior_model)

      if ( strat > 0.0_cp .and. DissNb > 0.0_cp ) then
         call abortRun('! Please give either strat or DissNb in the input Namelist!')
      end if

      if ( strat > 0.0_cp .or. DissNb > 0.0_cp ) l_anel= .true. 

      if ( index(interior_model,'EARTH') /= 0 ) then
         l_anel=.true.
      else if ( index(interior_model, 'JUP') /= 0 ) then
         l_anel=.true.
      else if ( index(interior_model, 'SAT') /= 0 ) then
         l_anel=.true.
      else if ( index(interior_model, 'SUN') /= 0 ) then
         l_anel=.true.
      else if ( index(interior_model, 'GLIESE229B') /= 0 ) then
         l_anel=.true.
      else if ( index(interior_model, 'COROT3B') /= 0 ) then
         l_anel=.true.
      else if ( index(interior_model, 'KOI889B') /= 0 ) then
         l_anel=.true.
      end if

      if ( prmag == 0.0_cp ) then
         l_mag   =.false.
         l_mag_nl=.false.
         l_mag_LF=.false.
      end if
      if ( .not. l_conv ) then
         l_conv_nl=.false.
         l_mag_LF =.false.
      end if

      !-- If Poincaré number is not zero precession in turned on
      if ( po == 0.0_cp ) then
         l_precession = .false.
      else
         l_precession = .true.
      end if

      if ( l_precession ) prec_angle = prec_angle*pi/180.0_cp

      !-- Same as above for differential precession

      if ( po_diff == 0.0_cp ) then
         l_diff_prec = .false.
      else
         l_diff_prec = .true.
         l_rot_ma = .true.
         l_rot_ic = .true.
         l_SRMA = .true.
         l_SRIC = .true.
      end if

      if ( l_diff_prec ) diff_prec_angle = diff_prec_angle*pi/180.0_cp

      !-- New checking of magnetic boundary condition.
      if ( kbotb > 4 ) then
         call abortRun('! Only outer boundary conditions kbotb<=4 implemented!')
      end if
      if ( sigma_ratio == 0.0_cp ) then
         l_cond_ic=.false.
         if ( kbotb == 3 ) then
            call abortRun('! For an insulating  IC with sigma_ratio=0, kbotb=3 is not appropriate')
         end if
      else
         l_cond_ic=.true.      ! tell the code to use a conducting inner core
         if ( kbotb  /=  3 ) then
            call abortRun('! For a conducting IC with sigma_ratio>0, kbotb=3 is appropriate')
         end if
      end if

      if ( ktopb > 4 ) then
         call abortRun('! Only outer boundary conditions ktopb<=4 implemented!')
      end if
      if ( conductance_ma == 0.0_cp ) then
         l_cond_ma=.false.
         if ( ktopb == 3 ) then
            call abortRun('! For an insulating mantle with conductance_ma=0, ktopb=3 is not appropriate')
         end if
      else
         l_cond_ma=.true.      ! tell the code to use a conducting mantle
         if ( ktopb  /=  3 ) then
            call abortRun('! For a conducting mantle with conductance_ma>0, ktopb=3 is appropriate')
         end if
      end if

      if ( .not. l_mag ) then
         prmag    =0.0_cp
         l_mag_nl =.false.
         l_cond_ic=.false.
         lMagMem  =0
      end if

      if ( l_corrMov ) l_par= .true. 

      !--- Stuff for the radial non-linear mapping
      call capitalize(map_function)

      !-- This is the case of a tangent mapping (see Bayliss & Turkel, 1992)
      if ( index(map_function, 'TAN') /= 0 .or. index(map_function, 'BAY') /= 0 ) then
         !     alph1 can be any positive number, above 0
         !     alph2 has to be a value between -1 and 1 (interval in Chebyshev space)
         if ( (alph1 == 0.0_cp) .or. (alph2 < -one) .or. (alph2 > one) ) then
            call abortRun('! Chebyshev mapping parameter is not correct !')
         elseif ( l_newmap .and. (alph1 < 0.0_cp) ) then
            alph1=abs(alph1)
         end if
      !-- This is the case of the Kosloff & Tal-Ezer mapping (1993)
      else if ( index(map_function, 'ARCSIN') /= 0 .or. &
      &         index(map_function, 'KTL') /= 0 ) then
         if ( (alph1 < 0.0_cp) .or. (alph1 >= one) ) then
            call abortRun('! Chebyshev mapping parameter is not correct !')
         end if
      end if

      !--- Stuff for current carrying loop at equator

      if (l_curr .and. amp_curr == 0.0_cp) then
         call abortRun('! For runs with l_curr please provide amp_curr')
      end if

      !--- Stuff for spherical magnetosphere boundary: rrMP=r(magnetosphere)/r_core
      if ( n_imp /= 0 ) imagcon=0
      if ( n_imp == 1 .and. rrMP <= one ) then
         call abortRun('! For runs with n_imp=1 please provide rrMP>1')
      else if ( n_imp >= 2 .and. amp_imp == 0.0_cp ) then
         call abortRun('! For runs with n_imp>=2 please provide amp_imp')
      else if ((n_imp == 2 .or. n_imp ==3) .and. l_imp <= 0) then
         call abortRun('! For runs with n_imp=2 or n_imp=3 l_imp must be >=0')
      end if
      if ( imagcon /= 0 .and. tmagcon == 0 ) tmagcon=1.0e18_cp

      if ( imagcon == -10 ) then
         !----- This is used to test variable conductivity
         !      with an analytical solution, see s_initB.f
         l_conv   =.false.
         l_conv_nl=.false.
         l_mag_kin=.true.
         l_mag_nl =.false.
         l_mag_LF =.false.
         l_heat   =.false.
         l_heat_nl=.false.
         l_cond_ic=.false.
         l_cond_ma=.false.
         kbotv    =2
         ktopv    =2
         kbotb    =1
         ktopb    =1
         nRotMa   =0
         nRotIc   =0
         radratio =half
      end if

      if ( l_rot_ma ) then
         write(*,*)
         write(*,*) '! I ALLOW FOR ROTATING MANTLE.'
         if ( ktopv == 1 .and. .not. l_cond_ma ) then
            write(*,*)
            write(*,*) '! No torques on mantle!'
            write(*,*) '! I dont update mantle rotation omega_ma.'
            l_rot_ma=.false.
         end if
      end if

      if ( .not. l_mag ) then
         l_cmb_field   =.false.
         l_dt_cmb_field=.false.
         l_rMagSpec    =.false.
         l_DTrMagSpec  =.false.
         l_storeBpot   =.false.
      end if

      l_b_nl_icb=.false.
      if ( l_mag_nl .and. kbotv == 1 .and. l_cond_ic ) then
         write(*,*)
         write(*,*) '! Nonlinear magnetic BC required at ICB!'
         l_b_nl_icb=.true.
      end if

      l_b_nl_cmb=.false.
      if ( l_mag_nl .and. ktopv == 1 .and. l_cond_ma ) then
         write(*,*)
         write(*,*) '! Nonlinear magnetic BC required at CMB!'
         l_b_nl_cmb=.true.
      end if

      !-- Special matrix for z(l=1,m=0) which is the solid body rotation:
      l_z10mat=.false.
      if ( ( l_rot_ma .and. ktopv == 2 ) .or. &
           ( l_rot_ic .and. kbotv == 2 )      ) l_z10mat= .true. 

      !-- Check Courant criteria at even time steps:
      if ( mod(n_cour_step,2) /= 0 ) n_cour_step=n_cour_step+1

      !-- Check whether memory has been reserved:
      if ( l_TO ) lStressMem=1
      if ( l_RMS .or. l_DTrMagSpec ) ldtBMem=1
      if ( l_movie .or. l_TOmovie ) lMovieMem=1

      !-- Output of angular moment?
      l_AM=l_AM .or. l_correct_AMe .or. l_correct_AMz
      l_AM=l_AM .or. l_power

      !-- Heat boundary condition
      if ( impS /= 0 ) then
         rad=pi/180
         n_impS=0
         do n=1,4*n_impS_max,4
            if ( sCMB(n) /= 0.0_cp ) then
               n_impS=n_impS+1
               peakS(n_impS) =sCMB(n)
               thetaS(n_impS)=rad*sCMB(n+1)
               phiS(n_impS)  =rad*sCMB(n+2)
               widthS(n_impS)=rad*sCMB(n+3)
            end if
         end do
      end if

      !-- Chemical composition boundary condition
      if ( impXi /= 0 ) then
         rad=pi/180
         n_impXi=0
         do n=1,4*n_impXi_max,4
            if ( xiCMB(n) /= 0.0_cp ) then
               n_impXi=n_impXi+1
               peakXi(n_impXi) =xiCMB(n)
               thetaXi(n_impXi)=rad*xiCMB(n+1)
               phiXi(n_impXi)  =rad*xiCMB(n+2)
               widthXi(n_impXi)=rad*xiCMB(n+3)
            end if
         end do
      end if

      !-- Grenoble stuff:
      lGrenoble=.false.
      if ( BIC /= 0.0_cp .and. l_mag ) then
         lGrenoble=.true.
         write(*,*)
         write(*,*) '! Running the Grenoble case !'
      end if

      ! Setting up truncation is required to set up ldif and l_max_r
      call initialize_truncation

      !-- Coeffs at radial levels:
      if ( l_r_fieldT ) l_r_field=.true.

      if ( l_r_field ) then
         if ( n_r_step == 0 ) n_r_step=2
         if ( l_max_r == 0 .or. l_max_r > l_max ) l_max_r=l_max
         nCounts = count(n_r_array>0)

         if ( nCounts > 0 ) then
             n_coeff_r_max=nCounts
         else
             n_coeff_r_max=5
         end if
      end if

      if ( l_energy_modes ) then
         if ( m_max_modes==0 .or. m_max_modes>l_max ) m_max_modes=l_max
      end if

      !-- ldif determins at which l hyperdiffusivity starts:
      ldif=max(1,ldif)
      ldif=min(l_max,ldif)

      !-- Coeffs at CMB:
      l_max_cmb=min(l_max_cmb,l_max)

      !-- Maximum run time specified?
      runTimeLimit(1)=runHours
      runTimeLimit(2)=runMinutes
      runTimeLimit(3)=runSeconds
      runTimeLimit(4)=0
      l_runTimeLimit =.false.
      if ( runHours+runMinutes+runSeconds > 0 ) l_runTimeLimit= .true. 

   end subroutine readNamelists
!------------------------------------------------------------------------------
   subroutine writeNamelists(n_out)
      !
      !  Purpose of this subroutine is to write the namelist to           
      !  file unit n_out. This file has to be open before calling this    
      !  routine.                                                         
      !

      !-- Input variable:
      integer, intent(in) :: n_out

      !-- Local variables:
      integer :: l,m,n,i
      real(cp) ::  rad
      integer :: length


      !-- Output of name lists:

      write(n_out,*) " "
      write(n_out,*) "&grid"
      write(n_out,'(''  n_r_max         ='',i5,'','')') n_r_max
      write(n_out,'(''  n_cheb_max      ='',i5,'','')') n_cheb_max
      write(n_out,'(''  n_phi_tot       ='',i5,'','')') n_phi_tot
      write(n_out,'(''  n_theta_axi     ='',i5,'','')') n_theta_axi
      write(n_out,'(''  n_r_ic_max      ='',i5,'','')') n_r_ic_max
      write(n_out,'(''  n_cheb_ic_max   ='',i5,'','')') n_cheb_ic_max
      write(n_out,'(''  minc            ='',i5,'','')') minc
      write(n_out,'(''  nalias          ='',i5,'','')') nalias
      write(n_out,'(''  l_axi           ='',l3,'','')') l_axi
      write(n_out,'(''  fd_stretch      ='',ES14.6,'','')') fd_stretch
      write(n_out,'(''  fd_ratio        ='',ES14.6,'','')') fd_ratio
      write(n_out,'(''  fd_order        ='',i5,'','')') fd_order
      write(n_out,'(''  fd_order_bound  ='',i5,'','')') fd_order_bound
      write(n_out,*) "/"

      write(n_out,*) "&control"
      write(n_out,'(''  mode            ='',i3,'','')') mode
      length=length_to_blank(tag)
      write(n_out,*) " tag             = """,tag(1:length),""","
      write(n_out,'(''  n_time_steps    ='',i8,'','')') n_time_steps
      write(n_out,'(''  n_tScale        ='',i3,'','')') n_tScale
      write(n_out,'(''  n_lScale        ='',i3,'','')') n_lScale
      write(n_out,'(''  alpha           ='',ES14.6,'','')')   alpha
      write(n_out,'(''  enscale         ='',ES14.6,'','')') enscale
      write(n_out,'(''  l_update_v      ='',l3,'','')') l_update_v
      write(n_out,'(''  l_update_b      ='',l3,'','')') l_update_b
      write(n_out,'(''  l_update_s      ='',l3,'','')') l_update_s
      write(n_out,'(''  l_update_xi     ='',l3,'','')') l_update_xi
      write(n_out,'(''  l_newmap        ='',l3,'','')') l_newmap
      length=length_to_blank(map_function)
      write(n_out,*) " map_function    = """,map_function(1:length),""","
      write(n_out,'(''  alph1           ='',ES14.6,'','')') alph1
      write(n_out,'(''  alph2           ='',ES14.6,'','')') alph2
      write(n_out,'(''  dtstart         ='',ES14.6,'','')') dtstart*tScale
      write(n_out,'(''  dtMax           ='',ES14.6,'','')') tScale*dtMax
      write(n_out,'(''  courfac         ='',ES14.6,'','')') courfac
      write(n_out,'(''  alffac          ='',ES14.6,'','')')  alffac
      write(n_out,'(''  l_cour_alf_damp ='',l3,'','')') l_cour_alf_damp
      write(n_out,'(''  intfac          ='',ES14.6,'','')')  intfac
      write(n_out,'(''  n_cour_step     ='',i5,'','')') n_cour_step
      write(n_out,'(''  difnu           ='',ES14.6,'','')') difnu
      write(n_out,'(''  difeta          ='',ES14.6,'','')') difeta
      write(n_out,'(''  difkap          ='',ES14.6,'','')') difkap
      write(n_out,'(''  difchem         ='',ES14.6,'','')') difchem
      write(n_out,'(''  ldif            ='',i3,'','')') ldif
      write(n_out,'(''  ldifexp         ='',i3,'','')') ldifexp
      write(n_out,'(''  l_correct_AMe   ='',l3,'','')') l_correct_AMe
      write(n_out,'(''  l_correct_AMz   ='',l3,'','')') l_correct_AMz
      write(n_out,'(''  l_non_rot       ='',l3,'','')') l_non_rot
      write(n_out,'(''  l_runTimeLimit  ='',l3,'','')') l_runTimeLimit
      write(n_out,'(''  runHours        ='',i6,'','')') runTimeLimit(1)
      write(n_out,'(''  runMinutes      ='',i4,'','')') runTimeLimit(2)
      write(n_out,'(''  runSeconds      ='',i4,'','')') runTimeLimit(3)
      write(n_out,'(''  tEND            ='',ES14.6,'','')') tEND
      length=length_to_blank(radial_scheme)
      write(n_out,*) " radial_scheme   = """,radial_scheme(1:length),""","
      length=length_to_blank(polo_flow_eq)
      write(n_out,*) " polo_flow_eq    = """,polo_flow_eq(1:length),""","
      length=length_to_blank(anelastic_flavour)
      write(n_out,*) "anelastic_flavour= """,anelastic_flavour(1:length),""","
      length=length_to_blank(thermo_variable)
      write(n_out,*) " thermo_variable = """,thermo_variable(1:length),""","
      write(n_out,*) "/"

      write(n_out,*) "&phys_param"
      write(n_out,'(''  ra              ='',ES14.6,'','')') ra
      write(n_out,'(''  raxi            ='',ES14.6,'','')') raxi
      write(n_out,'(''  pr              ='',ES14.6,'','')') pr
      write(n_out,'(''  sc              ='',ES14.6,'','')') sc
      write(n_out,'(''  prmag           ='',ES14.6,'','')') prmag
      write(n_out,'(''  ek              ='',ES14.6,'','')') ek
      write(n_out,'(''  po              ='',ES14.6,'','')') po
      write(n_out,'(''  prec_angle      ='',ES14.6,'','')') prec_angle
      write(n_out,'(''  po_diff         ='',ES14.6,'','')') po_diff
      write(n_out,'(''  diff_prec_angle ='',ES14.6,'','')') diff_prec_angle
      write(n_out,'(''  epsc0           ='',ES14.6,'','')') epsc0/sq4pi
      write(n_out,'(''  epscxi0         ='',ES14.6,'','')') epscxi0/sq4pi
      write(n_out,'(''  DissNb          ='',ES14.6,'','')') DissNb
      write(n_out,'(''  strat           ='',ES14.6,'','')') strat
      write(n_out,'(''  polind          ='',ES14.6,'','')') polind
      write(n_out,'(''  ThExpNb         ='',ES14.6,'','')') ThExpNb
      !write(n_out,'(''  GrunNb          ='',ES14.6,'','')') GrunNb
      write(n_out,'(''  epsS            ='',ES14.6,'','')') epsS
      write(n_out,'(''  cmbHflux        ='',ES14.6,'','')') cmbHflux
      write(n_out,'(''  slopeStrat      ='',ES14.6,'','')') slopeStrat
      write(n_out,'(''  rStrat          ='',ES14.6,'','')') rStrat
      write(n_out,'(''  ampStrat        ='',ES14.6,'','')') ampStrat
      write(n_out,'(''  thickStrat      ='',ES14.6,'','')') thickStrat
      write(n_out,'(''  nVarEntropyGrad ='',i3,'','')') nVarEntropyGrad
      write(n_out,'(''  radratio        ='',ES14.6,'','')') radratio
      write(n_out,'(''  l_isothermal    ='',l3,'','')') l_isothermal
      length=length_to_blank(interior_model)
      write(n_out,*) " interior_model  = """,interior_model(1:length),""","
      write(n_out,'(''  g0              ='',ES14.6,'','')') g0
      write(n_out,'(''  g1              ='',ES14.6,'','')') g1
      write(n_out,'(''  g2              ='',ES14.6,'','')') g2
      write(n_out,'(''  ktopv           ='',i3,'','')') ktopv
      write(n_out,'(''  kbotv           ='',i3,'','')') kbotv
      write(n_out,'(''  ktopb           ='',i3,'','')') ktopb
      write(n_out,'(''  kbotb           ='',i3,'','')') kbotb

      !-- Spherically-symmetric pressure
      write(n_out,'(''  ktopp           ='',i3,'','')') ktopp

      !--- Heat boundary condition:
      write(n_out,'(''  ktops           ='',i3,'','')') ktops
      write(n_out,'(''  kbots           ='',i3,'','')') kbots
      write(n_out,'("  Bottom boundary l,m,S:")')
      do m=0,m_max,minc
          do l=m,l_max
              if ( bots(l,m) /= 0.0_cp ) write(n_out,'(1p,4x,2i4,2ES14.6)') &
                   l,m,real(bots(l,m))/sq4pi,aimag(bots(l,m))/sq4pi
          end do
      end do
      write(n_out,'("  Top boundary l,m,S:")')
      do m=0,m_max,minc
          do l=m,l_max
              if ( tops(l,m) /= 0.0_cp ) write(n_out,'(1p,4x,2i4,2ES14.6)') &
                   l,m,real(tops(l,m))/sq4pi,aimag(tops(l,m))/sq4pi
          end do
      end do
      write(n_out,'(''  impS            ='',i3,'','')') impS
      rad=pi/180.0_cp
      do i=1,n_impS
         if ( i == 1 ) then
            write(n_out,'(A)',advance='NO') "  sCMB        ="
         else
            write(n_out,'(A)',advance='NO') "               "
         end if
         write(n_out,'(1p,4(ES14.6,A))') peakS(i)/rad,",", &
              & thetaS(i)/rad,",", phiS(i)/rad,",",        &
              & widthS(i)/rad,","
      end do

      !--- Chemical composition boundary condition:
      if ( l_chemical_conv ) then
         write(n_out,'(''  ktopxi          ='',i3,'','')') ktopxi
         write(n_out,'(''  kbotxi          ='',i3,'','')') kbotxi
         write(n_out,'("  Bottom boundary l,m,Xi:")')
         do m=0,m_max,minc
             do l=m,l_max
                 if ( botxi(l,m) /= 0.0_cp ) write(n_out,'(1p,4x,2i4,2ES14.6)') &
                      l,m,real(botxi(l,m))/sq4pi,aimag(botxi(l,m))/sq4pi
             end do
         end do
         write(n_out,'("  Top boundary l,m,Xi:")')
         do m=0,m_max,minc
             do l=m,l_max
                 if ( topxi(l,m) /= 0.0_cp ) write(n_out,'(1p,4x,2i4,2ES14.6)') &
                      l,m,real(topxi(l,m))/sq4pi,aimag(topxi(l,m))/sq4pi
             end do
         end do
         write(n_out,'(''  impXi           ='',i3,'','')') impXi
         rad=pi/180.0_cp
         do i=1,n_impXi
            if ( i == 1 ) then
               write(n_out,'(A)',advance='NO') "  xiCMB       ="
            else
               write(n_out,'(A)',advance='NO') "               "
            end if
            write(n_out,'(1p,4(ES14.6,A))') peakXi(i)/rad,",", &
                 & thetaXi(i)/rad,",", phiXi(i)/rad,",",        &
                 & widthXi(i)/rad,","
         end do
      end if

      !----- Conductivity variation:
      write(n_out,'(''  nVarCond        ='',i3,'','')') nVarCond
      write(n_out,'(''  con_DecRate     ='',ES14.6,'','')') con_DecRate
      write(n_out,'(''  con_RadRatio    ='',ES14.6,'','')') con_RadRatio
      write(n_out,'(''  con_LambdaMatch ='',ES14.6,'','')') con_LambdaMatch
      write(n_out,'(''  con_LambdaOut   ='',ES14.6,'','')') con_LambdaOut
      write(n_out,'(''  con_FuncWidth   ='',ES14.6,'','')') con_FuncWidth
      write(n_out,'(''  r_LCR           ='',ES14.6,'','')') r_LCR

      !----- Thermal diffusivity variation:
      write(n_out,'(''  difExp          ='',ES14.6,'','')') difExp
      write(n_out,'(''  nVarDiff        ='',i3,'','')') nVarDiff

      !----- Variable kinematic viscosity:
      write(n_out,'(''  nVarVisc        ='',i3,'','')') nVarVisc

      !----- Internal heating form:
      write(n_out,'(''  nVarEps         ='',i3,'','')') nVarEps

      write(n_out,*) "/"

      !----- External field
      write(n_out,*) "&B_external"
      write(n_out,'(''  n_imp           ='',i3,'','')') n_imp
      write(n_out,'(''  l_imp           ='',i3,'','')') l_imp
      write(n_out,'(''  rrMP            ='',ES14.6,'','')') rrMP
      write(n_out,'(''  amp_imp         ='',ES14.6,'','')') amp_imp
      write(n_out,'(''  expo_imp        ='',ES14.6,'','')') expo_imp
      write(n_out,'(''  bmax_imp        ='',ES14.6,'','')') bmax_imp

      write(n_out,'(''  l_curr          ='',l3,'','')') l_curr
      write(n_out,'(''  amp_curr        ='',ES14.6,'','')') amp_curr

      write(n_out,*) "/"


      write(n_out,*) "&start_field"
      write(n_out,'(''  l_start_file    ='',l3,'','')') l_start_file
      length=length_to_blank(start_file)
      write(n_out,*) " start_file      = """,start_file(1:length),""","
      write(n_out,'(''  inform          ='',i3,'','')') inform
      write(n_out,'(''  l_reset_t       ='',l3,'','')') l_reset_t
      write(n_out,'(''  scale_s         ='',ES14.6,'','')') scale_s
      write(n_out,'(''  scale_xi        ='',ES14.6,'','')') scale_xi
      write(n_out,'(''  scale_b         ='',ES14.6,'','')') scale_b
      write(n_out,'(''  scale_v         ='',ES14.6,'','')') scale_v
      write(n_out,'(''  tipdipole       ='',ES14.6,'','')') tipdipole
      write(n_out,'(''  init_s1         ='',i7,'','')') init_s1
      write(n_out,'(''  init_s2         ='',i3,'','')') init_s2
      write(n_out,'(''  init_v1         ='',i3,'','')') init_v1
      write(n_out,'(''  init_b1         ='',i3,'','')') init_b1
      write(n_out,'(''  init_xi1        ='',i7,'','')') init_xi1
      write(n_out,'(''  init_xi2        ='',i3,'','')') init_xi2
      write(n_out,'(''  imagcon         ='',i3,'','')') imagcon
      write(n_out,'(''  amp_s1          ='',ES14.6,'','')') amp_s1
      write(n_out,'(''  amp_s2          ='',ES14.6,'','')') amp_s2
      write(n_out,'(''  amp_v1          ='',ES14.6,'','')') amp_v1
      write(n_out,'(''  amp_b1          ='',ES14.6,'','')') amp_b1
      write(n_out,'(''  amp_xi1         ='',ES14.6,'','')') amp_xi1
      write(n_out,'(''  amp_xi2         ='',ES14.6,'','')') amp_xi2
      write(n_out,*) "/"

      write(n_out,*) "&output_control"
      write(n_out,'(''  n_graph_step    ='',i5,'','')') n_graph_step
      write(n_out,'(''  n_graphs        ='',i5,'','')') n_graphs
      write(n_out,'(''  t_graph_start   ='',ES14.6,'','')') t_graph_start
      write(n_out,'(''  t_graph_stop    ='',ES14.6,'','')') t_graph_stop
      write(n_out,'(''  dt_graph        ='',ES14.6,'','')') dt_graph
      write(n_out,'(''  n_rst_step      ='',i5,'','')') n_rst_step
      write(n_out,'(''  n_rsts          ='',i5,'','')') n_rsts
      write(n_out,'(''  t_rst_start     ='',ES14.6,'','')') t_rst_start
      write(n_out,'(''  t_rst_stop      ='',ES14.6,'','')') t_rst_stop
      write(n_out,'(''  dt_rst          ='',ES14.6,'','')') dt_rst
      write(n_out,'(''  n_stores        ='',i5,'','')') n_stores
      write(n_out,'(''  n_log_step      ='',i5,'','')') n_log_step
      write(n_out,'(''  n_logs          ='',i5,'','')') n_logs
      write(n_out,'(''  t_log_start     ='',ES14.6,'','')') t_log_start
      write(n_out,'(''  t_log_stop      ='',ES14.6,'','')') t_log_stop
      write(n_out,'(''  dt_log          ='',ES14.6,'','')') dt_log
      write(n_out,'(''  n_spec_step     ='',i5,'','')') n_spec_step
      write(n_out,'(''  n_specs         ='',i5,'','')') n_specs
      write(n_out,'(''  t_spec_start    ='',ES14.6,'','')') t_spec_start
      write(n_out,'(''  t_spec_stop     ='',ES14.6,'','')') t_spec_stop
      write(n_out,'(''  dt_spec         ='',ES14.6,'','')') dt_spec
      write(n_out,'(''  n_cmb_step      ='',i5,'','')') n_cmb_step
      write(n_out,'(''  n_cmbs          ='',i5,'','')') n_cmbs
      write(n_out,'(''  t_cmb_start     ='',ES14.6,'','')') t_cmb_start
      write(n_out,'(''  t_cmb_stop      ='',ES14.6,'','')') t_cmb_stop
      write(n_out,'(''  dt_cmb          ='',ES14.6,'','')') dt_cmb
      write(n_out,'(''  n_r_field_step  ='',i5,'','')') n_r_field_step
      write(n_out,'(''  n_r_fields      ='',i5,'','')') n_r_fields
      write(n_out,'(''  t_r_field_start ='',ES14.6,'','')') t_r_field_start
      write(n_out,'(''  t_r_field_stop  ='',ES14.6,'','')') t_r_field_stop
      write(n_out,'(''  dt_r_field      ='',ES14.6,'','')') dt_r_field
      write(n_out,'(''  l_movie         ='',l3,'','')') l_movie
      write(n_out,'(''  n_movie_step    ='',i5,'','')') n_movie_step
      write(n_out,'(''  n_movie_frames  ='',i5,'','')') n_movie_frames
      write(n_out,'(''  t_movie_start   ='',ES14.6,'','')') t_movie_start
      write(n_out,'(''  t_movie_stop    ='',ES14.6,'','')') t_movie_stop
      write(n_out,'(''  dt_movie        ='',ES14.6,'','')') dt_movie
      do n=1,n_movies_max
         length=len(trim(movie(n)))
         if ( length > 0 ) then
            write(n_out,'(''  movie           = '',a,'','')') movie(n)(1:length)
         end if
      end do
      write(n_out,'(''  l_probe         ='',l3,'','')') l_probe
      write(n_out,'(''  n_probe_step    ='',i5,'','')') n_probe_step
      write(n_out,'(''  n_probe_out     ='',i5,'','')') n_probe_out
      write(n_out,'(''  t_probe_start   ='',ES14.6,'','')') t_probe_start
      write(n_out,'(''  t_probe_stop    ='',ES14.6,'','')') t_probe_stop
      write(n_out,'(''  dt_probe        ='',ES14.6,'','')') dt_probe
      write(n_out,'(''  r_probe         ='',ES14.6,'','')') r_probe
      write(n_out,'(''  theta_probe     ='',ES14.6,'','')') theta_probe
      write(n_out,'(''  n_phi_probes    ='',i3,'','')') n_phi_probes
      write(n_out,'(''  l_average       ='',l3,'','')') l_average
      write(n_out,'(''  l_cmb_field     ='',l3,'','')') l_cmb_field
      write(n_out,'(''  l_dt_cmb_field  ='',l3,'','')') l_dt_cmb_field
      write(n_out,'(''  l_save_out      ='',l3,'','')') l_save_out
      write(n_out,'(''  l_true_time     ='',l3,'','')') l_true_time
      write(n_out,'(''  lVerbose        ='',l3,'','')') lVerbose
      write(n_out,'(''  l_rMagSpec      ='',l3,'','')') l_rMagSpec
      write(n_out,'(''  l_DTrMagSpec    ='',l3,'','')') l_DTrMagSpec
      write(n_out,'(''  l_max_cmb       ='',i3,'','')') l_max_cmb
      write(n_out,'(''  l_r_field       ='',l3,'','')') l_r_field
      write(n_out,'(''  l_r_fieldT      ='',l3,'','')') l_r_fieldT
      write(n_out,'(''  l_max_r         ='',i3,'','')') l_max_r
      write(n_out,'(''  n_r_step        ='',i3,'','')') n_r_step
      do n=1,n_coeff_r_max
         write(n_out,'(''  n_coeff_r       ='',i3,'','')') n_coeff_r(n)
      end do
      write(n_out,'(''  l_earth_likeness='',l3,'','')') l_earth_likeness
      write(n_out,'(''  l_max_comp      ='',i3,'','')') l_max_comp
      write(n_out,'(''  l_hel           ='',l3,'','')') l_hel
      write(n_out,'(''  l_AM            ='',l3,'','')') l_AM
      write(n_out,'(''  l_power         ='',l3,'','')') l_power
      write(n_out,'(''  l_viscBcCalc    ='',l3,'','')') l_viscBcCalc
      write(n_out,'(''  l_fluxProfs     ='',l3,'','')') l_fluxProfs
      write(n_out,'(''  l_perpPar       ='',l3,'','')') l_perpPar
      write(n_out,'(''  l_PressGraph    ='',l3,'','')') l_PressGraph
      write(n_out,'(''  l_energy_modes  ='',l3,'','')') l_energy_modes
      write(n_out,'(''  m_max_modes     ='',i3,'','')') m_max_modes
      write(n_out,'(''  l_drift         ='',l3,'','')') l_drift
      write(n_out,'(''  l_iner          ='',l3,'','')') l_iner
      write(n_out,'(''  l_TO            ='',l3,'','')') l_TO
      write(n_out,'(''  l_TOmovie       ='',l3,'','')') l_TOmovie
      write(n_out,'(''  l_PV            ='',l3,'','')') l_PV
      write(n_out,'(''  l_storeBpot     ='',l3,'','')') l_storeBpot
      write(n_out,'(''  l_storeVpot     ='',l3,'','')') l_storeVpot
      write(n_out,'(''  l_RMS           ='',l3,'','')') l_RMS
      write(n_out,'(''  l_par           ='',l3,'','')') l_par
      write(n_out,'(''  l_corrMov       ='',l3,'','')') l_corrMov
      write(n_out,'(''  rCut            ='',ES14.6,'','')') rCut
      write(n_out,'(''  rDea            ='',ES14.6,'','')') rDea
      write(n_out,*) "/"

      write(n_out,*) "&mantle"
      write(n_out,'(''  conductance_ma  ='',ES14.6,'','')') conductance_ma
      write(n_out,'(''  rho_ratio_ma    ='',ES14.6,'','')') rho_ratio_ma
      write(n_out,'(''  nRotMa          ='',i4,'','')') nRotMa
      write(n_out,'(''  omega_ma1       ='',ES14.6,'','')') omega_ma1
      write(n_out,'(''  omegaOsz_ma1    ='',ES14.6,'','')') omegaOsz_ma1
      write(n_out,'(''  tShift_ma1      ='',ES14.6,'','')') tShift_ma1
      write(n_out,'(''  omega_ma2       ='',ES14.6,'','')') omega_ma2
      write(n_out,'(''  omegaOsz_ma2    ='',ES14.6,'','')') omegaOsz_ma2
      write(n_out,'(''  tShift_ma2      ='',ES14.6,'','')') tShift_ma2
      write(n_out,'(''  amp_RiMaAsym    ='',ES14.6,'','')') amp_RiMaAsym
      write(n_out,'(''  omega_RiMaAsym  ='',ES14.6,'','')') omega_RiMaAsym
      write(n_out,'(''  m_RiMaAsym      ='',i4,'','')') m_RiMaAsym
      write(n_out,'(''  amp_RiMaSym     ='',ES14.6,'','')') amp_RiMaSym
      write(n_out,'(''  omega_RiMaSym   ='',ES14.6,'','')') omega_RiMaSym
      write(n_out,'(''  m_RiMaSym       ='',i4,'','')')  m_RiMaSym
      write(n_out,*) "/"

      write(n_out,*) "&inner_core"
      write(n_out,'(''  sigma_ratio     ='',ES14.6,'','')') sigma_ratio
      write(n_out,'(''  rho_ratio_ic    ='',ES14.6,'','')') rho_ratio_ic
      write(n_out,'(''  nRotIc          ='',i4,'','')') nRotIc
      write(n_out,'(''  omega_ic1       ='',ES14.6,'','')') omega_ic1
      write(n_out,'(''  omegaOsz_ic1    ='',ES14.6,'','')') omegaOsz_ic1
      write(n_out,'(''  tShift_ic1      ='',ES14.6,'','')') tShift_ic1
      write(n_out,'(''  omega_ic2       ='',ES14.6,'','')') omega_ic2
      write(n_out,'(''  omegaOsz_ic2    ='',ES14.6,'','')') omegaOsz_ic2
      write(n_out,'(''  tShift_ic2      ='',ES14.6,'','')') tShift_ic2
      write(n_out,'(''  BIC             ='',ES14.6,'','')') BIC
      write(n_out,'(''  amp_RiIcAsym    ='',ES14.6,'','')') amp_RiIcAsym
      write(n_out,'(''  omega_RiIcAsym  ='',ES14.6,'','')') omega_RiIcAsym
      write(n_out,'(''  m_RiIcAsym      ='',i4,'','')') m_RiIcAsym
      write(n_out,'(''  amp_RiIcSym     ='',ES14.6,'','')') amp_RiIcSym
      write(n_out,'(''  omega_RiIcSym   ='',ES14.6,'','')') omega_RiIcSym
      write(n_out,'(''  m_RiIcSym       ='',i4,'','')')  m_RiIcSym
      write(n_out,*) "/"
      write(n_out,*) " "

   end subroutine writeNamelists
!------------------------------------------------------------------------------
   subroutine defaultNamelists
      !
      !  Purpose of this subroutine is to set default parameters          
      !  for the namelists.                                               
      !

      !-- Local variable:
      integer :: n

      !----- Namelist grid
      ! must be of form 4*integer+1
      ! Possible values for n_r_max:
      !  5,9,13,17,21,25,33,37,41,49,
      !  61,65,73,81,97,101,121,129, ...
      n_r_max       =33
      ! max degree-1 of cheb polynomia
      n_cheb_max    =31
      ! number of longitude grid points
      ! Possible values: 
      ! 16,32,48,64,96,128,192,256,288,320,384,
      ! 400,512,576,640,768,864,1024
      n_phi_tot     =192
      n_theta_axi   =0
      ! number of grid points in inner core
      n_r_ic_max    =17
      ! number of chebs in inner core
      n_cheb_ic_max =15
      ! basic wavenumber, longitude symmetry
      minc          =1
      ! controls dealiasing in latitude and 
      ! longitude direction, no aliasing for nalias=20
      !   20 <= nalias <= 30
      nalias        =20
      l_axi         =.false.

      !-- Finite differences
      fd_order      =2
      fd_order_bound=2
      fd_stretch    =0.3_cp
      fd_ratio      =0.1_cp

      !----- Namelist control
      mode          =0            ! self-consistent dynamo !
      tag           ="default"
      n_time_steps  =100
      n_tScale      =0
      n_lScale      =0
      alpha         =half
      enscale       =one
      dtstart       =0.0_cp
      dtMax         =1.0e-4_cp
      courfac       =2.5_cp
      l_cour_alf_damp=.true. ! By default, use Christensen's (GJI, 1999) CFL
      alffac        =one
      intfac        =0.15_cp
      n_cour_step   =10
      anelastic_flavour="None" ! Useless in Boussinesq
      thermo_variable  ="None" 
      polo_flow_eq     ="WP"   ! Choose between 'DC' (double-curl) and 'WP' (Pressure)
      radial_scheme    ="CHEB" ! Choose between 'CHEB' and 'FD'

      cacheblock_size_in_B=4096

      l_update_v    =.true.
      l_update_b    =.true.
      l_update_s    =.true.
      l_update_xi   =.true.
      l_correct_AMe =.false.  ! Correct equatorial AM
      l_correct_AMz =.false.  ! Correct axial AM
      l_non_rot     =.false.  ! No Coriolis force !
      l_anel        =.false.  ! Anelastic stuff !
      l_isothermal  =.false.  ! Isothermal = 0 Grünesein !
      interior_model="None"   ! Name of the interior model

      !---- Run time and number of threads:
      l_runTimeLimit=.false. ! Control of absolute run time
      do n=1,4
         runTimeLimit(n)=0
      end do

      tEND          =0.0_cp    ! numerical time where run should end

      !----- Hyperdiffusion:
      difnu         =0.0_cp
      difeta        =0.0_cp
      difkap        =0.0_cp
      difchem       =0.0_cp
      ldif          =1
      ldifexp       =-1

      !----- Namelist phys_param:
      ra         =0.0_cp
      raxi       =0.0_cp
      ek         =1.0e-3_cp
      pr         =one
      sc         =10.0_cp
      prmag      =5.0_cp
      po         =0.0_cp
      prec_angle =23.5_cp
      po_diff    =0.0_cp
      diff_prec_angle =23.5_cp
      epsc0      =0.0_cp
      epscxi0    =0.0_cp
      radratio   =0.35_cp
      !----- Anelatic stuff
      DissNb     =0.0_cp     ! Dissipation number
      ThExpNb    =one        ! Thermal expansion * temperature
      strat      =0.0_cp     ! Density contrast
      polind     =2.0_cp     ! Polytropic index
      GrunNb     =one/polind ! Gruneisen parameter (by default this is adiabatic)
      r_cut_model=0.98_cp    ! outer radius when using interior model
      !----- Stably stratified layer
      epsS       =0.0_cp
      cmbHflux   =0.0_cp
      slopeStrat =20.0_cp
      rStrat     =1.3_cp
      ampStrat   =10.0_cp
      thickStrat =0.1_cp
      nVarEntropyGrad=0
      !----- Gravity parameters: defaut value g propto r (i.e. g1=1)
      g0         =0.0_cp
      g1         =one
      g2         =0.0_cp        
      !----- Boundary conditions        
      ktops      =1
      kbots      =1
      ktopxi     =1
      kbotxi     =1
      ktopv      =2
      kbotv      =2
      ktopb      =1
      kbotb      =1
      ktopp      =1
      do n=1,4*n_s_bounds
         s_top(n)=0.0_cp
         s_bot(n)=0.0_cp
      end do
      impS=0
      do n=1,n_impS_max
         peakS(n) =0.0_cp
         thetaS(n)=0.0_cp
         phiS(n)  =0.0_cp
         widthS(n)=0.0_cp
      end do
      do n=1,4*n_xi_bounds
         xi_top(n)=0.0_cp
         xi_bot(n)=0.0_cp
      end do
      impXi=0
      do n=1,n_impXi_max
         peakXi(n) =0.0_cp
         thetaXi(n)=0.0_cp
         phiXi(n)  =0.0_cp
         widthXi(n)=0.0_cp
      end do

      !----- Conductivity variation:
      nVarCond       =0
      con_DecRate    =9.0_cp
      con_RadRatio   =0.75_cp
      con_LambdaMatch=0.6_cp
      con_LambdaOut  =0.1
      con_FuncWidth  =0.25_cp
      r_LCR          =two

      !----- Thermal diffusivity variation:
      difExp         =-half
      nVarDiff       =0

      !----- Variable kinematic viscosity:
      nVarVisc       =0

      !----- Internal heating form:
      nVarEps        =0

      !----- Non-linear mapping parameters (Bayliss, 1992):
      l_newmap       =.false.
      alph1          =0.8_cp
      alph2          =0.0_cp
      map_function   ='arcsin' ! By default Kosloff and Tal-Ezer mapping when l_newmap=.true.

      !----- External field
      n_imp          =0    ! No external field
      rrMP           =0.0_cp ! r(Magnetopause)/r_cmb, used for n_imp=1
      amp_imp        =0.0_cp ! amplitude of external field
      expo_imp       =0.0_cp ! oscillation frequency of external field
      bmax_imp       =0.0_cp
      l_imp          =1    ! Default external field is axial dipole

      l_curr         =.false. !No current loop
      amp_curr       =0.0_cp  !Current loop switched off

      !----- Namelist start_field:
      l_start_file  =.false.
      start_file    ="no_start_file"
      inform        =-1   
      runid         ="MAGIC default run"
      l_reset_t     =.false.
      scale_s       =one
      scale_xi      =one
      scale_b       =one
      scale_v       =one
      tipdipole     =0.0_cp
      init_s1       =0
      init_s2       =0
      init_xi1      =0
      init_xi2      =0
      init_b1       =0
      init_v1       =0
      imagcon       =0
      tmagcon       =0.0_cp
      amp_s1        =one
      amp_s2        =0.0_cp
      amp_v1        =0.0_cp
      amp_b1        =one
      amp_xi1       =0.0_cp
      amp_xi2       =0.0_cp

      !----- Namelist output_control:
      l_save_out    =.false.  ! Save output
      l_true_time   =.false.  ! Use exact requested output times
      lVerbose      =.false.  ! Tell me what you are doing
      l_average     =.false.  ! Average various quantities in time

      !----- Restart files:
      n_rst_step    =0
      n_rsts        =1
      t_rst_start   =0.0_cp
      t_rst_stop    =0.0_cp
      dt_rst        =0.0_cp
      n_stores      =0

      !----- Log output:
      n_log_step    =50
      n_logs        =0
      t_log_start   =0.0_cp
      t_log_stop    =0.0_cp
      dt_log        =0.0_cp

      !----- Graphic output:
      n_graph_step  =0
      n_graphs      =1
      t_graph_start =0.0_cp
      t_graph_stop  =0.0_cp
      dt_graph      =0.0_cp
      l_graph_time  =.false.

      !----- Spectrum files:
      n_spec_step   =0
      n_specs       =0
      t_spec_start  =0.0_cp
      t_spec_stop   =0.0_cp
      dt_spec       =0.0_cp

      !----- Output of poloidal magnetic field potential at CMB:
      !      also stored at times of movie frames
      l_cmb_field   =.false.
      l_dt_cmb_field=.false.
      l_max_cmb     =14
      n_cmb_step    =0
      n_cmbs        =0
      t_cmb_start   =0.0_cp
      t_cmb_stop    =0.0_cp
      dt_cmb        =0.0_cp

      !----- Output of magnetic and flow potential af five different radial levels:
      l_r_field     =.false.
      l_r_fieldT    =.false.
      l_max_r       =l_max
      n_r_step      =2
      do n=1,size(n_r_array)
         n_r_array(n)=0
      end do
      n_r_field_step =0
      n_r_fields     =0
      t_r_field_start=0.0_cp
      t_r_field_stop =0.0_cp
      dt_r_field     =0.0_cp

      !----- Compute Earth-likeness (Christensen et al. EPSL 2010)
      l_earth_likeness=.false.
      l_max_comp    =8

      !----- Output of distribution of energies over m's
      l_energy_modes=.false. ! to get emag and ekin for different m
      m_max_modes   =14      ! number of modes

      !----- Movie output:
      l_movie       =.false.
      n_movies      =0
      n_movie_step  =0
      n_movie_frames=0
      t_movie_start =0.0_cp
      t_movie_stop  =0.0_cp
      dt_movie      =0.0_cp
      do n=1,n_movies_max
         movie(n)=' '
      end do

      !----- Output from probes:
      l_probe       =.false.
      n_probe_step  =0
      n_probe_out   =0
      t_probe_start =0.0_cp
      t_probe_stop  =0.0_cp
      dt_probe      =0.0_cp
      n_phi_probes  =0
      r_probe       =0.0_cp
      theta_probe   =0.0_cp

      !----- Output of magnetic potentials:
      l_storeBpot   =.false.
      n_Bpot_step   =0
      n_Bpots       =0
      t_Bpot_start  =0.0_cp
      t_Bpot_stop   =0.0_cp
      dt_Bpot       =0.0_cp

      !----- Output of flow potentials:
      l_storeVpot   =.false.
      n_Vpot_step   =0
      n_Vpots       =0
      t_Vpot_start  =0.0_cp
      t_Vpot_stop   =0.0_cp
      dt_Vpot       =0.0_cp

      !----- Output of T potential:
      l_storeTpot   =.false.
      n_Tpot_step   =0
      n_Tpots       =0
      t_Tpot_start  =0.0_cp
      t_Tpot_stop   =0.0_cp
      dt_Tpot       =0.0_cp

      !----- Output of all potential:
      l_storePot    =.false.
      n_pot_step    =0
      n_pots        =0
      t_pot_start   =0.0_cp
      t_pot_stop    =0.0_cp
      dt_pot        =0.0_cp

      !----- Output TOZ:
      n_TOZ_step    =0
      n_TOZs        =0
      t_TOZ_start   =0.0_cp
      t_TOZ_stop    =0.0_cp
      dt_TOZ        =0.0_cp

      !----- Output TO:
      n_TO_step     =0
      n_TOs         =0
      t_TO_start    =0.0_cp
      t_TO_stop     =0.0_cp
      dt_TO         =0.0_cp

      !----- Times for different output:
      do n=1,n_time_hits
         t_graph(n)  =-one
         t_rst(n)    =-one
         t_log(n)    =-one
         t_spec(n)   =-one
         t_cmb(n)    =-one
         t_r_field(n)=-one
         t_movie(n)  =-one
         t_Vpot(n)   =-one
         t_Bpot(n)   =-one
         t_Tpot(n)   =-one
         t_pot(n)    =-one
         t_TO(n)     =-one
         t_TOZ(n)    =-one
         t_TOmovie(n)=-one
         t_probe     =-one
      end do

      !----- Magnetic spectra for different depths
      !      at times of log output or movie frames:
      l_rMagSpec    =.false.
      l_DTrMagSpec  =.false.

      !----- TO output, output times same as for log output:
      l_TO          =.false. ! TO output in TOnhs.TAG, TOshs.TAG
      l_TOmovie     =.false. ! TO movies 
      sDens         =one     ! relative s-grid point density 
      zDens         =one     ! relative z-grid point density 

      !----- Potential vorticity:
      l_PV          =.false.

      !----- Different output, output times same as for log outout:
      l_hel         =.false. ! Helicity in misc.TAG 
      l_AM          =.false. ! Angular moment in AM.TAG 
      l_power       =.false. ! power budget in power.TAG and dtE.TAG
      l_viscBcCalc  =.false. ! dissipation layer for stress-free BCs
      l_fluxProfs   =.false. ! radial profiles of flux contributions
      l_perpPar     =.false. ! radial profiles and time series of kinetic energy
                             ! perpendicular and parallel to the rotation axi
      l_PressGraph  =.true.  ! store pressure in graphic files
      l_drift       =.false. ! files for calculating drift rates 
      l_iner        =.false. ! files for calculating inertial modes
      l_RMS         =.false. ! RMS force balance and dynamo term 
                             ! balance in dtVrms.TAG and dtBrms.TAG
      l_par         =.false. ! Calculate additional parameters in s_getEgeos.f
      l_corrMov     =.false. ! North/south correlation movie (see s_getEgeos.f)
      rCut          =0.0_cp  ! Thickness of layer to be left out at both
      ! boundaries for RMS calculation.
      ! rCut=0.075 means that 7.5% at the CMB and ICB are disregarded.
      rDea          =0.0_cp  ! Controls dealiazing in  RMS calculation
      ! rDea=0.1 means that highest 10% of cheb modes are set to zero

      !----- Mantle name list:
      conductance_ma=0.0_cp    ! insulation mantle is default
      nRotMa        =0         ! non rotating mantle is default
      rho_ratio_ma  =one       ! same density as outer core
      omega_ma1     =0.0_cp    ! prescribed rotation rate
      omegaOsz_ma1  =0.0_cp    ! oszillation frequency of mantle rotation rate
      tShift_ma1    =0.0_cp    ! time shift
      omega_ma2     =0.0_cp    ! second mantle rotation rate 
      omegaOsz_ma2  =0.0_cp    ! oscillation frequency of second mantle rotation
      tShift_ma2    =0.0_cp    ! time shift for second rotation
      amp_RiMaAsym  =0.0_cp    ! amplitude of Rieutord forcing (eq anti-symm)
      omega_RiMaAsym=0.0_cp    ! frequency of Rieutord forcing (eq anti-symm)
      m_RiMaAsym    =0         ! default forcing -> axisymmetric (eq anti-symm)
      amp_RiMaSym   =0.0_cp    ! amplitude of Rieutord forcing (eq symm)
      omega_RiMaSym =0.0_cp    ! frequency of Rieutord forcing (eq symm)
      m_RiMaSym     =0         ! default forcing -> axisymmetric (eq symm)

      !----- Inner core name list:
      sigma_ratio   =0.0_cp    ! no conducting inner core is default 
      nRotIc        =0         ! non rotating inner core is default
      rho_ratio_ic  =one       ! same density as outer core
      omega_ic1     =0.0_cp    ! prescribed rotation rate, added to first one
      omegaOsz_ic1  =0.0_cp    ! oszillation frequency of IC rotation rate
      tShift_ic1    =0.0_cp    ! time shift
      omega_ic2     =0.0_cp    ! second prescribed rotation rate 
      omegaOsz_ic2  =0.0_cp    ! oszillation frequency of second IC rotation rate
      tShift_ic2    =0.0_cp    ! tims shift for second IC rotation
      BIC           =0.0_cp    ! Imposed dipole field strength at ICB
      amp_RiIcAsym  =0.0_cp    ! amplitude of Rieutord forcing (eq anti-symm)
      omega_RiIcAsym=0.0_cp    ! frequency of Rieutord forcing (eq anti-symm)
      m_RiIcAsym    =0         ! default forcing -> axisymmetric (eq anti-symm)
      amp_RiIcSym  =0.0_cp     ! amplitude of Rieutord forcing (eq symm)
      omega_RiIcSym=0.0_cp     ! frequency of Rieutord forcing (eq symm)
      m_RiIcSym    =0          ! default forcing -> axisymmetric (eq symm)


   end subroutine defaultNamelists
!------------------------------------------------------------------------------
end module Namelists
