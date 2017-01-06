module logic
   !
   ! Module containing the logicals that control the run
   !

   implicit none
 
   logical :: l_update_v     ! Switch off velocity field update
   logical :: l_update_b     ! Switch off magnetic field update
   logical :: l_update_s     ! Switch off entropy update
   logical :: l_update_xi    ! Switch off update of chemical composition
   logical :: l_mag          ! Switch off magnetic terms calculation
   logical :: l_conv         ! Switch off convection
   logical :: l_mag_kin      ! Switch related for kinematic dynamo
   logical :: l_SRIC         ! Switch to rotating IC with prescribed rot. rate
   logical :: l_SRMA         ! Switch to rotating OC with prescribed rot. rate
   logical :: l_heat         ! Switch off heat terms calculation
   logical :: l_heat_nl      ! Switch off non-linear heat terms calculation
   logical :: l_mag_nl       ! Switch off non-linear magnetic terms calculation
   logical :: l_conv_nl      ! Switch off non-linear convection terms
   logical :: l_mag_LF       ! Switch off Lorentz force term
   logical :: l_corr         ! Switch off rotation
   logical :: l_rot_ic       ! Switch off IC rotation
   logical :: l_rot_ma       ! Switch off OC rotation
   logical :: l_z10mat       ! Switch for solid body rotation
   logical :: l_cond_ic      ! Switch for conducting IC
   logical :: l_cond_ma      ! Switch for conducting OC
   logical :: l_time_hits    ! Switch for time for outputs
   logical :: l_average      ! Switch for calculation of time-averages
   logical :: l_energy_modes ! Switch for calculation of distribution of energies over m's
   logical :: l_movie        ! Switch for recording of movie files
   logical :: l_movie_oc     ! Switch for recording of movie files for OC
   logical :: l_movie_ic     ! Switch for recording of movie files for IC
   logical :: l_save_out     ! Switch off outputs
   logical :: l_true_time    ! Switch for times of outputs
   logical :: l_cmb_field    ! Switch for Bcoef files for gauss coefficients
   logical :: l_dt_cmb_field ! Switch for Bcoef files for secular variation of gauss coefs.
   logical :: l_storeBpot    ! Switch for storing magnetic field potentials
   logical :: l_storeVpot    ! Switch for storing velocity field potentials
   logical :: l_storeTpot    ! Switch for storing entropy field potentials
   logical :: l_storePot     ! Switch for storing all field potentials
   logical :: l_r_field      ! Switch for radial coefficients
   logical :: l_r_fieldT     ! Switch for radial T coefficients
   logical :: l_b_nl_cmb     ! Switch for non-linear magnetic field at OC
   logical :: l_b_nl_icb     ! Switch for non-linear magnetic field at IC
   logical :: l_correct_AMe  ! Switch for correction of equatorial angular mom.
   logical :: l_correct_AMz  ! Switch for correction of axial angular momentum
   logical :: l_HTmovie      ! Switch for heat flux movie output
   logical :: l_HT           ! Switch for heat flux movie frame output
   logical :: l_dtBmovie     ! Switch for dtB movie
   logical :: l_dtB          ! Switch to reserve memory for dtB movie
   logical :: l_store_frame  ! Switch for storing movie frames
   logical :: l_non_rot      ! Switch to non-rotating
   logical :: l_rMagSpec     ! Switch for magnetic spectra at different depths at log times
   logical :: l_DTrMagSpec   ! Switch for magnetic spectra at different depths at movie output times
   logical :: l_TO           ! Switch for TO output in TOnhs.TAG, TOshs.TAG
   logical :: l_TOmovie      ! Switch for TO movie output 
   logical :: l_hel          ! Switch for helicity calculation, output in misc.TAG
   logical :: l_anel         ! Switch for anelastic calculation
   logical :: l_isothermal   ! Switch for isothermal calculation
   logical :: l_anelastic_liquid ! Switch for anelastic liquid calculation
   logical :: l_AM           ! Switch for angular momentum calculation
   logical :: l_power        ! Switch for power budget terms calculation
   logical :: l_drift        ! Switch for drift rates calculation
   logical :: l_iner         ! Switch for inertial modes calculation
   logical :: l_runTimeLimit ! Switch for absolute time limit of the run
   logical :: l_RMS          ! Switch for RMS force balances calculation
   logical :: l_par          ! Switch for additional parameters calculation in s_getEgeos.f
   logical :: l_corrMov      ! Switch for North/south correlation movie (see s_getEgeos.f)
   logical :: l_PV           ! Switch for potential vorticity calculation
   logical :: l_newmap       ! Switch for non-linear mapping (see Bayliss and Turkel, 1990)
   logical :: l_viscBcCalc   ! Switch for dissipation layer for stress-free BCs plots
   logical :: l_fluxProfs    ! Switch for calculation of radial profiles of flux contributions
   logical :: l_perpPar      ! Switch for calculation of of kinetic energy perpendicular+parallel to the rotation axis
   logical :: l_LCR          ! Switch for zero electrical conductivity beyond r_LCR
   logical :: lVerbose       ! Switch for detailed information about run progress

   logical :: l_PressGraph   ! Store pressure in graphic files

   logical :: l_single_matrix ! In case entropy, w and P are solved at once implicitely
   logical :: l_temperature_diff ! diffusion of temperature instead of entropy

   logical :: l_chemical_conv ! Switch for chemical convection
   logical :: l_non_adia ! Switch in case the reference state is non-adiabatic

   logical :: l_TP_form ! Use temperature and pressure instead of entropy and pressure
   logical :: l_probe        ! Switch for artifical sensors

end module logic
