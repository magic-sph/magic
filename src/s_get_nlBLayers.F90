!$Id$
!***********************************************************************
SUBROUTINE get_nlBLayers(vt,vp,dvtdr,dvpdr,   &
     &                   dsdr,dsdt,dsdp,&
     &                   uhLMr,duhLMr,gradsLMr,nR,nThetaStart)
  !***********************************************************************

  !-----------------------------------------------------------------------

  !   Calculates axisymmetric contributions of the horizontal velocity
  !       uh = sqrt(utheta^2+uphi^2)
  !   its radial derivative
  !       abs(duh/dr)
  !   and the thermal dissipation rate
  !       (grad T)**2
  !
  !   This subroutine is used when one wants to evaluate viscous and thermal
  !   dissipation layers
  !
  !-----------------------------------------------------------------------

  USE truncation
  USE radial_functions, ONLY: orho2,or2,or1,beta
  USE blocking, ONLY: sizeThetaB,nfs
  USE horizontal_data, ONLY: O_sin_theta_E2

  IMPLICIT NONE

  !-- Input of variables
  INTEGER,intent(IN) :: nR
  INTEGER,intent(IN) :: nThetaStart

  !-- Input of fields components:
  REAL(kind=8),intent(IN) :: vt(nrp,nfs),vp(nrp,nfs)
  REAL(kind=8),intent(IN) :: dvtdr(nrp,nfs),dvpdr(nrp,nfs)
  REAL(kind=8),intent(IN) :: dsdr(nrp,nfs),dsdt(nrp,nfs),dsdp(nrp,nfs)

  !-- Output:
  REAL(kind=8),intent(OUT) :: uhLMr(l_max+1)
  REAL(kind=8),intent(OUT) :: duhLMr(l_max+1)
  REAL(kind=8),intent(OUT) :: gradsLMr(l_max+1)

  !-- Local:
  INTEGER :: nTheta,nThetaB
  INTEGER :: nPhi
  REAL(kind=8) :: uhAS(nfs),duhAS(nfs),gradsAS(nfs),uh,duh,phiNorm,grads

  !-- End of declaration
  !---------------------------------------------------------------------------

  phiNorm=1.D0/DBLE(n_phi_max)

  !IF ( nThetaStart == 1 ) THEN
     !------ Zero lm coeffs for first theta block:
  !   uhLMr =0.D0
  !   duhLMr=0.D0
  !   gradsLMr=0.D0
  !END IF

  !--- Horizontal velocity uh and duh/dr + (grad T)**2
  nTheta=nThetaStart-1
  DO nThetaB=1,sizeThetaB
     nTheta=nTheta+1
     uhAS(nThetaB) =0.D0
     duhAS(nThetaB)=0.D0
     gradsAS(nThetaB)=0.D0
     DO nPhi=1,n_phi_max
        !WRITE(*,"(2I4,A,2ES20.12)") nThetaB,nPhi,": vt,vp = ", vt(nPhi,nThetaB),vp(nPhi,nThetaB)
        uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(          &
             vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+  &
             vp(nPhi,nThetaB)*vp(nPhi,nThetaB)  )
        duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(         &
             dvtdr(nPhi,nThetaB)*vt(nPhi,nThetaB)-&
             (or1(nR)+beta(nR))*vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+&
             dvpdr(nPhi,nThetaB)*vp(nPhi,nThetaB)-&
             (or1(nR)+beta(nR))*vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )

        grads =  dsdr(nPhi,nThetaB)*dsdr(nPhi,nThetaB)          &
             +or2(nR)*O_sin_theta_E2(nTheta)*(                &
             dsdt(nPhi,nThetaB)*dsdt(nPhi,nThetaB)    &
             +dsdp(nPhi,nThetaB)*dsdp(nPhi,nThetaB) ) 

        uhAS(nThetaB)=uhAS(nThetaB)+DSQRT(uh)
        IF (uh.NE.0.0d0) THEN
           duhAS(nThetaB)=duhAS(nThetaB)+DABS(duh)/DSQRT(uh)
        END IF
        gradsAS(nThetaB)=gradsAS(nThetaB)+grads
     END DO
     uhAS(nThetaB)=phiNorm*uhAS(nThetaB)
     duhAS(nThetaB)=phiNorm*duhAS(nThetaB)
     gradsAS(nThetaB)=phiNorm*gradsAS(nThetaB)
  END DO

  !------ Add contribution from thetas in block:
  !WRITE(*,"(2I5,A,2ES20.12)") nR,nThetaStart,": uhAS = ",SUM(uhAS),SUM(duhAS)
  CALL legtfAS2(uhLMr,duhLMr,uhAS,duhAS,l_max+1,nThetaStart,sizeThetaB)
  CALL legtfAS(gradsLMr,gradsAS,l_max+1,nThetaStart,sizeThetaB)
  !WRITE(*,"(2I5,A,ES20.12)") nR,nThetaStart,": duhLMr = ",SUM(duhLMr)
  RETURN
end SUBROUTINE get_nlBLayers
!-----------------------------------------------------------------------------
