!$Id$
module physical_parameters
!************************************************************************
!  Module containing the physical parameters
!************************************************************************

   implicit none
 
   !-- Control parameters for boundary conditions:
   integer :: mode         ! Mode of calculation
   integer :: ktops,kbots  ! Entropy boundary condition
   integer :: ktopv,kbotv  ! Velocity boundary condition
   integer :: ktopb,kbotb  ! Magnetic boundary condition
   !-- Parameters for a localized temperature (entropy) disturbance at CMB
   integer :: impS,n_impS
   integer, parameter :: n_impS_max=20
   real(kind=8) ::  thetaS(n_impS_max),phiS(n_impS_max)
   real(kind=8) :: peakS(n_impS_max),widthS(n_impS_max) 
 
   !-- Dimensionless parameters:
   real(kind=8) :: radratio
   real(kind=8) :: ra,ek,pr,prmag,epsc0,epsc
   real(kind=8) :: strat,polind,ViscHeatFac,OhmLossFac
   real(kind=8) :: DissNb,epsS ! Dissipation number, deviation from the adiabat
   real(kind=8) :: cmbHflux
   real(kind=8) :: slopeStrat ! stratified Layer
   character(len=72) :: interior_model ! name of the interior model
   real(kind=8) :: r_cut_model
   real(kind=8) :: g0,g1,g2
   real(kind=8) :: sigma_ratio,conductance_ma
   real(kind=8) :: rho_ratio_ic,rho_ratio_ma
   real(kind=8) :: opr,opm,CorFac,LFfac,BuoFac,O_sr
   real(kind=8) :: raScaled,ekScaled
 
   !-- Variable properties:
   integer :: nVarCond
   real(kind=8) :: difExp
   real(kind=8) :: con_DecRate,con_RadRatio,con_LambdaMatch
   real(kind=8) :: con_LambdaOut,con_FuncWidth
   real(kind=8) :: r_LCR
   integer :: n_r_LCR
   integer :: nVarDiff
   integer :: nVarVisc
   integer :: nVarEps

   !-- To avoid circular dependence
   integer :: imagcon
   real(kind=8) :: tmagcon


end module physical_parameters
