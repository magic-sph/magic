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
   integer :: impS,n_impS
   integer, parameter :: n_impS_max=20
   real(cp) ::  thetaS(n_impS_max),phiS(n_impS_max)
   real(cp) :: peakS(n_impS_max),widthS(n_impS_max) 
 
   !-- Dimensionless parameters:
   real(cp) :: radratio
   real(cp) :: ra,ek,pr,prmag,epsc0,epsc
   real(cp) :: strat,polind,ViscHeatFac,OhmLossFac
   real(cp) :: DissNb  ! Dissipation number
   real(cp) :: epsS    ! deviation from the adiabat
   real(cp) :: cmbHflux
   real(cp) :: slopeStrat ! stratified Layer
   character(len=72) :: interior_model ! name of the interior model
   real(cp) :: r_cut_model
   real(cp) :: g0,g1,g2
   real(cp) :: sigma_ratio,conductance_ma
   real(cp) :: rho_ratio_ic,rho_ratio_ma
   real(cp) :: opr,opm,CorFac,LFfac,BuoFac,O_sr
   real(cp) :: raScaled,ekScaled
 
   !-- Variable properties:
   integer :: nVarCond
   real(cp) :: difExp
   real(cp) :: con_DecRate,con_RadRatio,con_LambdaMatch
   real(cp) :: con_LambdaOut,con_FuncWidth
   real(cp) :: r_LCR
   integer :: n_r_LCR
   integer :: nVarDiff
   integer :: nVarVisc
   integer :: nVarEps

   !-- To avoid circular dependence
   integer :: imagcon
   real(cp) :: tmagcon

end module physical_parameters
