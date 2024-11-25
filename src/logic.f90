module logic
   !
   ! Module containing the logicals that control the run
   !

   implicit none

   logical :: l_update_v     ! Switch off velocity field update
   logical :: l_update_b     ! Switch off magnetic field update
   logical :: l_update_s     ! Switch off entropy update
   logical :: l_update_xi    ! Switch off update of chemical composition
   logical :: l_update_phi   ! Switch off update of phase field
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
   logical :: l_average      ! Switch for calculation of time-averages
   logical :: l_spec_avg     ! Switch for calculation of time-averaged spectra
   logical :: l_energy_modes ! Switch for calculation of distribution of energies over m's
   logical :: l_movie        ! Switch for recording of movie files
   logical :: l_movie_oc     ! Switch for recording of movie files for OC
   logical :: l_movie_ic     ! Switch for recording of movie files for IC
   logical :: l_save_out     ! Switch off outputs
   logical :: l_cmb_field    ! Switch for Bcoef files for gauss coefficients
   logical :: l_dt_cmb_field ! Switch for Bcoef files for secular variation of gauss coefs.
   logical :: l_2D_spectra   ! Switch for storing of r-l-spectra
   logical :: l_2D_RMS       ! Switch for storing of time-averaged r-l-spectra of forces
   logical :: l_r_field      ! Switch for radial coefficients
   logical :: l_r_fieldT     ! Switch for radial T coefficients
   logical :: l_r_fieldXi    ! Switch for radial Xi coefficients
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
   logical :: l_newmap       ! Switch for non-linear mapping (see Bayliss and Turkel, 1990)
   logical :: l_viscBcCalc   ! Switch for dissipation layer for stress-free BCs plots
   logical :: l_fluxProfs    ! Switch for calculation of radial profiles of flux contributions
   logical :: l_perpPar      ! Switch for calculation of of kinetic energy perpendicular+parallel to the rotation axis
   logical :: l_LCR          ! Switch for zero electrical conductivity beyond r_LCR
   logical :: lVerbose       ! Switch for detailed information about run progress
   logical :: l_ehd_dep      ! Switch for dilectrophoretic force

   logical :: l_PressGraph   ! Store pressure in graphic files

   logical :: l_single_matrix ! In case entropy, w and P are solved at once implicitely
   logical :: l_temperature_diff ! diffusion of temperature instead of entropy

   logical :: l_chemical_conv ! Switch for chemical convection
   logical :: l_phase_field ! Switch when phase field is used
   logical :: l_non_adia ! Switch in case the reference state is non-adiabatic

   logical :: l_probe        ! Switch for artifical sensors

   logical :: l_finite_diff ! Use finite differences for the radial scheme
   logical :: l_double_curl ! Use the double-curl of the NS equation to get the poloidal equation
   logical :: l_AB1  ! 1st order Adams Bashforth
   logical :: l_cour_alf_damp ! Modified Alfven Courant condition based on Christensen et al., GJI, 1999 (.true. by default)
   logical :: l_earth_likeness ! Compute the Earth-likeness of the CMB field following Christensen et al., EPSL, 2010
   logical :: l_precession ! Use precession
   logical :: l_centrifuge ! Compute centrifugal acceleration
   logical :: l_adv_curl   ! Use \curl{u}\times u for the advection
   logical :: l_full_sphere ! Set to .true. if this is a full sphere calculation
   logical :: l_var_l ! When set to .true., degree varies with radius
   logical :: l_bridge_step ! Used to bridge missing steps when changing the time integrator
   logical :: l_packed_transp ! Pack or don't pack MPI transposes
   logical :: l_parallel_solve ! Use R-distributed parallel solver (work only for F.D.)
   logical :: l_mag_par_solve ! Can be remove once inner core has also been ported
   logical :: l_hemi ! Compute North/South asymmetry of energies
   logical :: l_onset ! A flag to turn MagIC into a linear stability analysis code
   logical :: l_scramble_theta ! A flag to set theta scrambling
   logical :: l_geosMovie ! A flag to trigger the production of geos movies
   logical :: l_phaseMovie ! A flag to trigger the production of a movie for the melting radius
   logical :: l_dtphaseMovie ! A flag to trigger the production of a movie for the temperature gradient at the melting radius

end module logic
