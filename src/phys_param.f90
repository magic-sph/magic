module physical_parameters
   !
   !  Module containing the physical parameters
   !

   use precision_mod

   implicit none
 
   !-- Control parameters for boundary conditions:
   integer :: mode         ! Mode of calculation
   integer :: ktops,kbots  ! Entropy boundary condition
   integer :: ktopv,kbotv  ! Velocity boundary condition
   integer :: ktopb,kbotb  ! Magnetic boundary condition
   !-- Parameters for a localized temperature (entropy) disturbance at CMB
   integer :: impS         ! Heat boundary condition
   integer :: n_impS       ! Heat boundary condition
   integer, parameter :: n_impS_max=20 ! Heat boundary condition
   real(cp) :: thetaS(n_impS_max)
   real(cp) :: phiS(n_impS_max)
   real(cp) :: peakS(n_impS_max)
   real(cp) :: widthS(n_impS_max) 
 
   !-- Dimensionless parameters:
   real(cp) :: radratio       ! aspect ratio
   real(cp) :: ra             ! Rayleigh number
   real(cp) :: ek             ! Ekman number
   real(cp) :: pr             ! Prandtl number
   real(cp) :: prmag          ! magnetic Prandtl number
   real(cp) :: epsc0          ! Internal heat source magnitude
   real(cp) :: epsc           ! Renormalisation of epsc0
   real(cp) :: strat          ! number of density scale heights
   real(cp) :: polind         ! polytropic index
   real(cp) :: ViscHeatFac    ! Dissipation number *Pr/raScaled
   real(cp) :: OhmLossFac     ! Dissipation number ViscHeatFac/(ekScaled*prmag**2)
   real(cp) :: DissNb         ! Dissipation number
   real(cp) :: epsS           ! Deviation from the adiabat
   real(cp) :: cmbHflux       ! stratified Layer
   real(cp) :: slopeStrat     ! stratified Layer
   character(len=72) :: interior_model ! name of the interior model
   real(cp) :: r_cut_model    ! Percentage on the inner part of the interior model to be used
   real(cp) :: g0             ! Set to 1.0 for constant gravity
   real(cp) :: g1             ! Set to 1.0 for linear gravity
   real(cp) :: g2             ! Set to 1.0 for 1/r**2 gravity
   real(cp) :: sigma_ratio    ! Value of IC rotation
   real(cp) :: conductance_ma ! OC conductivity
   real(cp) :: rho_ratio_ic   ! Same density as outer core
   real(cp) :: rho_ratio_ma   ! Same density as outer core
   real(cp) :: opr            ! Inverse of Prandtl number
   real(cp) :: opm            ! Inverse of magnetic Prandtl number
   real(cp) :: CorFac         ! Inverse of ekScaled
   real(cp) :: LFfac          ! Inverse of Pr*Ekman
   real(cp) :: BuoFac         ! Ratio of Rayleigh number over Prandtl number
   real(cp) :: O_sr           ! Inverse of sigma_ratio
   real(cp) :: raScaled       ! Ra*lscale**3
   real(cp) :: ekScaled       ! Ekman*lscale**2
 
   !-- Variable properties:
   integer :: nVarCond        ! Selection of variable conductivity profile
   real(cp) :: difExp         ! Thermal diffusivity variation
   real(cp) :: con_DecRate    ! Slope of electrical conductivity profile (nVarCond=2)
   real(cp) :: con_RadRatio   ! Transition between branches of electrical conductivity profile (nVarCond=1,2)
   real(cp) :: con_LambdaMatch ! Electrical conductivity at con_RadRatio (nVarCond=2)
   real(cp) :: con_LambdaOut  ! nVarCond=1
   real(cp) :: con_FuncWidth  ! nVarCond=1
   real(cp) :: r_LCR          ! Radius beyond which conductivity is zero
   integer :: n_r_LCR         ! Number of radial points where conductivity is zero
   integer :: nVarDiff        ! Selection of variable diffusivity profile
   integer :: nVarVisc        ! Selection of variable viscosity profile
   integer :: nVarEps         ! Selection of internal heating profile

   !-- To avoid circular dependence
   integer :: imagcon         ! Imposed magnetic field for magnetoconvection, at the boundaries
   real(cp) :: tmagcon        ! Time for magnetoconvection calculation

end module physical_parameters
