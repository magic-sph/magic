!$Id$
!************************************************************************
!  Module containing the physical parameters
!************************************************************************
 
MODULE physical_parameters

  IMPLICIT NONE

  !-- Control parameters for boundary conditions:
  INTEGER :: mode         ! Mode of calculation
  INTEGER :: ktops,kbots  ! Entropy boundary condition
  INTEGER :: ktopv,kbotv  ! Velocity boundary condition
  INTEGER :: ktopb,kbotb  ! Magnetic boundary condition
  !-- Parameters for a localized temperature (entropy) disturbance at CMB
  INTEGER :: impS,n_impS
  INTEGER,PARAMETER :: n_impS_max=20
  REAL(kind=8) ::  thetaS(n_impS_max),phiS(n_impS_max)
  REAL(kind=8) :: peakS(n_impS_max),widthS(n_impS_max) 

  !-- Dimensionless parameters:
  REAL(kind=8) :: radratio
  REAL(kind=8) :: ra,ek,pr,prmag,epsc0,epsc
  REAL(kind=8) :: strat,polind,PolFac,ViscHeatFac,OhmLossFac
  REAL(kind=8) :: r_cut_model
  REAL(kind=8) :: g0,g1,g2
  REAL(kind=8) :: sigma_ratio,conductance_ma
  REAL(kind=8) :: rho_ratio_ic,rho_ratio_ma
  REAL(kind=8) :: opr,opm,CorFac,LFfac,BuoFac,O_sr
  REAL(kind=8) :: raScaled,ekScaled

  !-- Variable properties:
  INTEGER :: nVarCond
  REAL(kind=8) :: difExp
  REAL(kind=8) :: con_DecRate,con_RadRatio,con_LambdaMatch
  REAL(kind=8) :: con_LambdaOut,con_FuncWidth
  INTEGER :: nVarDiff
  INTEGER :: nVarVisc
  INTEGER :: nVarEps

END MODULE physical_parameters
