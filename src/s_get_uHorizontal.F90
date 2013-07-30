!$Id$
!***********************************************************************
    SUBROUTINE get_duHorizontal(vt,vp,dvtdr,dvpdr,uhLMr,duhLMr,nR,nThetaStart)
!***********************************************************************

!-----------------------------------------------------------------------

!   Calculates axisymmetric contributions of the derivative of the
!   horizontal velocity uh = sqrt(utheta^2+uphi^2), i.e.
!
!   abs(duh/dr)
!
!-----------------------------------------------------------------------

    USE truncation
    USE radial_functions, ONLY: orho2,or2,or1,beta
    USE blocking
    USE horizontal_data, ONLY: O_sin_theta_E2

    IMPLICIT NONE

!-- Input of variables
    INTEGER,intent(IN) :: nR
    INTEGER,intent(IN) :: nThetaStart

!-- Input of fields components:
    REAL(kind=8),intent(IN) :: vt(nrp,nfs),vp(nrp,nfs)
    REAL(kind=8),intent(IN) :: dvtdr(nrp,nfs),dvpdr(nrp,nfs)

!-- Output:
    REAL(kind=8),intent(OUT) :: uhLMr(l_max+1)
    REAL(kind=8),intent(OUT) :: duhLMr(l_max+1)

!-- Local:
    INTEGER :: nTheta,nThetaB
    INTEGER :: nPhi
    REAL(kind=8) :: uhAS(nfs),duhAS(nfs),uh,duh,phiNorm

!-- End of declaration
!---------------------------------------------------------------------------

    phiNorm=1.D0/DBLE(n_phi_max)

    IF ( nThetaStart == 1 ) THEN
    !------ Zero lm coeffs for first theta block:
        uhLMr =0.D0
        duhLMr=0.D0
    END IF

!--- Horizontal velocity uh and duh/dr
    nTheta=nThetaStart-1
    DO nThetaB=1,sizeThetaB
        nTheta=nTheta+1
        uhAS(nThetaB) =0.D0
        duhAS(nThetaB)=0.D0
        DO nPhi=1,n_phi_max
            uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(          &
                               vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+  &
                               vp(nPhi,nThetaB)*vp(nPhi,nThetaB)  )
            duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(         &
                              dvtdr(nPhi,nThetaB)*vt(nPhi,nThetaB)-&
              (or1(nR)+beta(nR))*vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+&
                              dvpdr(nPhi,nThetaB)*vp(nPhi,nThetaB)-&
              (or1(nR)+beta(nR))*vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )

            uhAS(nThetaB)=uhAS(nThetaB)+DSQRT(uh)
            duhAS(nThetaB)=duhAS(nThetaB)+DABS(duh)/DSQRT(uh)
        END DO
        uhAS(nThetaB)=phiNorm*uhAS(nThetaB)
        duhAS(nThetaB)=phiNorm*duhAS(nThetaB)
    END DO

!------ Add contribution from thetas in block:
    CALL legtfAS2(uhLMr,duhLMr,uhAS,duhAS,l_max+1,nThetaStart,sizeThetaB)

    RETURN
    end SUBROUTINE get_duHorizontal
!-----------------------------------------------------------------------------
