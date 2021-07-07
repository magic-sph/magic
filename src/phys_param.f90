module physical_parameters
   !
   !  Module containing the physical parameters
   !

   use precision_mod

   implicit none

   !-- Control parameters for boundary conditions:
   integer :: mode         ! Mode of calculation
   integer :: ktops,kbots  ! Entropy boundary condition
   integer :: ktopxi,kbotxi! Boundary conditions for chemical composition
   integer :: ktopv,kbotv  ! Velocity boundary condition
   integer :: ktopb,kbotb  ! Magnetic boundary condition
   integer :: ktopp        ! Boundary condition for spherically-symmetric pressure

   !-- Parameters for a localized temperature (entropy) disturbance at CMB
   integer :: impS         ! Heat boundary condition
   integer :: n_impS       ! Heat boundary condition
   integer, parameter :: n_impS_max=20 ! Heat boundary condition
   real(cp) :: thetaS(n_impS_max)
   real(cp) :: phiS(n_impS_max)
   real(cp) :: peakS(n_impS_max)
   real(cp) :: widthS(n_impS_max)

   !-- Parameters for a localized chemical disturbance at CMB
   integer :: impXi
   integer :: n_impXi
   integer, parameter :: n_impXi_max=20
   real(cp) :: thetaXi(n_impXi_max)
   real(cp) :: phiXi(n_impXi_max)
   real(cp) :: peakXi(n_impXi_max)
   real(cp) :: widthXi(n_impXi_max)

   !-- Dimensionless parameters:
   real(cp) :: radratio       ! aspect ratio
   real(cp) :: ra             ! Rayleigh number
   real(cp) :: raxi           ! Chemical composition-based Rayleigh number
   real(cp) :: sc             ! Schmidt number (i.e. chemical Prandtl number)
   real(cp) :: ek             ! Ekman number
   real(cp) :: pr             ! Prandtl number
   real(cp) :: prmag          ! magnetic Prandtl number
   real(cp) :: epsc0          ! Internal heat source magnitude
   real(cp) :: epscxi0        ! Internal chemical heat source magnitude
   real(cp) :: epsc           ! Renormalisation of epsc0
   real(cp) :: epscxi         ! Renormalisation of epsc0Xi
   real(cp) :: Bn             ! Normalisation of He burning
   real(cp) :: strat          ! number of density scale heights
   real(cp) :: polind         ! polytropic index
   real(cp) :: stef           ! Stefan number
   real(cp) :: ViscHeatFac    ! Prefactor for viscous heating: :math:`Di\,Pr/Ra`
   real(cp) :: OhmLossFac     ! Prefactor for Ohmic heating: :math:`Di\,Pr/(Ra\,E\,Pm^2)`
   real(cp) :: DissNb         ! Dissipation number
   real(cp) :: ThExpNb        ! Thermal expansion * temperature :math:`\alpha_0 T_0`
   real(cp) :: GrunNb         ! Grüneisen paramater :math:`\Gamma=(\gamma-1)/\alpha T`
   real(cp) :: epsS           ! Deviation from the adiabat
   real(cp) :: cmbHflux       ! stratified Layer
   real(cp) :: slopeStrat     ! stratified Layer
   real(cp) :: rStrat         ! stratified Layer
   real(cp) :: ampStrat       ! stratified Layer
   real(cp) :: thickStrat     ! stratified Layer
   integer  :: nVarEntropyGrad! stratified Layer
   character(len=72) :: interior_model ! name of the interior model
   real(cp) :: r_cut_model    ! Percentage on the inner part of the interior model to be used
   real(cp) :: g0             ! Set to 1.0 for constant gravity
   real(cp) :: g1             ! Set to 1.0 for linear gravity
   real(cp) :: g2             ! Set to 1.0 for :math:`1/r^2` gravity
   real(cp) :: sigma_ratio    ! Value of IC rotation
   real(cp) :: conductance_ma ! OC conductivity
   real(cp) :: rho_ratio_ic   ! Same density as outer core
   real(cp) :: rho_ratio_ma   ! Same density as outer core
   real(cp) :: opr            ! Inverse of Prandtl number
   real(cp) :: oek            ! Inverse of the Ekman number
   real(cp) :: osc            ! Inverse of Schmidt number (i.e. chemical Prandtl number)
   real(cp) :: opm            ! Inverse of magnetic Prandtl number
   real(cp) :: CorFac         ! Inverse of ekScaled
   real(cp) :: LFfac          ! Inverse of Pr*Ekman
   real(cp) :: BuoFac         ! Ratio of Rayleigh number over Prandtl number
   real(cp) :: ChemFac        ! Ratio of comp. Rayleigh number over Schmidt number
   real(cp) :: O_sr           ! Inverse of sigma_ratio
   real(cp) :: raScaled       ! :math:`Ra\,l^3`
   real(cp) :: ekScaled       ! :math:`E\,l^2`

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
   real(cp) :: po             ! Poincaré number
   real(cp) :: prec_angle     ! Precession angle
   real(cp) :: dilution_fac   ! Omega^2 d/g_top for centrifugal acceleration, named after Chandrasekhar (1987)
   real(cp) :: epsPhase       ! Cahn number for phase field equatioo
   real(cp) :: phaseDiffFac   ! Diffusion term of phase field
   real(cp) :: penaltyFac     ! Factor that enters the penalty method in the NS equations
   real(cp) :: tmelt          ! Melting temperature

end module physical_parameters
