!$Id$
!***********************************************************************
SUBROUTINE get_perpPar(vr,vt,vp,EperpLMr,EparLMr, &
                    &  EperpaxiLMr,EparaxiLMr,nR,nThetaStart)
!***********************************************************************

!-----------------------------------------------------------------------

!   Calculates the energies parallel and perpendicular to the rotation axis
!
!       Eperp = 0.5*(vs**2+vphi**2) with vs= vr*sin(theta)+vt*cos(theta)
!       Epar  = 0.5*vz**2           with vz= vr*cos(theta)-vt*sin(theta)
!
!-----------------------------------------------------------------------

  USE truncation
  USE radial_functions, ONLY: orho2,or2,or1
  USE blocking
  USE horizontal_data, ONLY: sn2,osn2,cosTheta

  IMPLICIT NONE

!-- Input of variables
  INTEGER,intent(IN) :: nR
  INTEGER,intent(IN) :: nThetaStart

!-- Input of fields components:
  REAL(kind=8),intent(IN) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)

!-- Output:
  REAL(kind=8),intent(OUT) :: EperpLMr(l_max+1),EparLMr(l_max+1)
  REAL(kind=8),intent(OUT) :: EperpaxiLMr(l_max+1),EparaxiLMr(l_max+1)

!-- Local:
  INTEGER :: nTheta,nThetaB,nThetaNHS
  INTEGER :: nPhi
  REAL(kind=8) :: vras,vtas,vpas,phiNorm
  REAL(kind=8) :: EperpAS(nfs),EparAS(nfs),Eperp,Epar
  REAL(kind=8) :: EperpaxiAS(nfs),EparaxiAS(nfs),Eperpaxi,Eparaxi

!-- End of declaration
!---------------------------------------------------------------------------

  phiNorm=1.D0/DBLE(n_phi_max)

  nTheta=nThetaStart-1
  DO nThetaB=1,sizeThetaB
     nTheta=nTheta+1
     nThetaNHS=(nTheta+1)/2

     EperpAS(nThetaB)   =0.D0
     EparAS(nThetaB)    =0.D0
     EperpaxiAS(nThetaB)=0.D0
     EparaxiAS(nThetaB) =0.D0
     Eperp   =0.D0
     Epar    =0.D0
     Eperpaxi=0.D0
     Eparaxi =0.D0
     vras    =0.D0
     vtas    =0.D0
     vpas    =0.D0

     DO nPhi=1,n_phi_max
        vras=vras+vr(nPhi,nThetaB)
        vtas=vtas+vt(nPhi,nThetaB)
        vpas=vpas+vp(nPhi,nThetaB)
     END DO
     vras=vras*phiNorm
     vtas=vtas*phiNorm
     vpas=vpas*phiNorm

     DO nPhi=1,n_phi_max
        Eperp=0.5D0*or2(nR)*orho2(nR)*(                                             &
                  or2(nR)*sn2(nThetaNHS)*       vr(nPhi,nThetaB)*vr(nPhi,nThetaB) + &
                  (osn2(nThetaNHS)-1.D0)*       vt(nPhi,nThetaB)*vt(nPhi,nThetaB) + &
                  2.D0*or1(nR)*cosTheta(nTheta)*vr(nPhi,nThetaB)*vt(nPhi,nThetaB) + &
                  osn2(nThetaNHS)*              vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )

        Epar =0.5D0*or2(nR)*orho2(nR)*(                                             &
                  or2(nR)*(1.D0-sn2(nThetaNHS))*vr(nPhi,nThetaB)*vr(nPhi,nThetaB) + &
                                                vt(nPhi,nThetaB)*vt(nPhi,nThetaB) - &
                  2.D0*or1(nR)*cosTheta(nTheta)*vr(nPhi,nThetaB)*vt(nPhi,nThetaB) )

        Eperpaxi=0.5D0*or2(nR)*orho2(nR)*(                  &
                  or2(nR)*sn2(nThetaNHS)*       vras*vras + &
                  (osn2(nThetaNHS)-1.D0)*       vtas*vtas + &
                  2.D0*or1(nR)*cosTheta(nTheta)*vras*vtas + &
                  osn2(nThetaNHS)*              vpas*vpas )

        Eparaxi =0.5D0*or2(nR)*orho2(nR)*(                  &
                  or2(nR)*(1.D0-sn2(nThetaNHS))*vras*vras + &
                                                vtas*vtas - &
                  2.D0*or1(nR)*cosTheta(nTheta)*vras*vtas )

        EperpAS(nThetaB)   =   EperpAS(nThetaB)+Eperp
        EparAS(nThetaB)    =    EparAS(nThetaB)+Epar
        EperpaxiAS(nThetaB)=EperpaxiAS(nThetaB)+Eperpaxi
        EparaxiAS(nThetaB) = EparaxiAS(nThetaB)+Eparaxi
     END DO
     EperpAS(nThetaB)   =phiNorm*   EperpAS(nThetaB)
     EparAS(nThetaB)    =phiNorm*    EparAS(nThetaB)
     EperpaxiAS(nThetaB)=phiNorm*EperpaxiAS(nThetaB)
     EparaxiAS(nThetaB) =phiNorm* EparaxiAS(nThetaB)
  END DO

!------ Add contribution from thetas in block:
  CALL legtfAS2(EperpLMr,EparLMr,EperpAS,EparAS,l_max+1,nThetaStart,sizeThetaB)
  CALL legtfAS2(EperpaxiLMr,EparaxiLMr,EperpaxiAS,EparaxiAS,l_max+1,nThetaStart,sizeThetaB)

  RETURN
END SUBROUTINE get_perpPar
!-----------------------------------------------------------------------------
