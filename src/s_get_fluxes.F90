!$Id$
!***********************************************************************
SUBROUTINE get_fluxes(vr,vt,vp,dvrdr,dvtdr,dvpdr,      &
                 &  dvrdt,dvrdp,sr,pr,br,bt,bp,cbt,cbp,&
                 &  fconvLMr,fkinLMr,fviscLMr,fpoynLMr,&
                 &  fresLMR,nR,nThetaStart)
!***********************************************************************

!-----------------------------------------------------------------------

!   Calculates the fluxes
!       Fconv= rho * temp* (u_r s)
!       Fkin = 1/2 * rho * u_r (u_r^2+u_theta^2+u_phi^2)
!       Fvisc= -(u *\tensor(S) )_r
!
!   This subroutine is used when one wants to evaluate viscous and thermal
!   dissipation layers
!
!-----------------------------------------------------------------------

  USE truncation
  USE const, ONLY: pi
  USE physical_parameters, ONLY: ek,ViscHeatFac
  USE radial_functions, ONLY: orho2,or2,temp0,or1,orho1,beta,visc, &
                              n_r_icb,n_r_cmb
  USE blocking
  USE horizontal_data, ONLY: O_sin_theta_E2,osn2

  IMPLICIT NONE

!-- Input of variables
  INTEGER,intent(IN) :: nR
  INTEGER,intent(IN) :: nThetaStart

!-- Input of fields components:
  REAL(kind=8),intent(IN) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
  REAL(kind=8),intent(IN) :: dvrdr(nrp,nfs),dvtdr(nrp,nfs),dvpdr(nrp,nfs)
  REAL(kind=8),intent(IN) :: dvrdt(nrp,nfs),dvrdp(nrp,nfs)
  REAL(kind=8),intent(IN) :: sr(nrp,nfs),pr(nrp,nfs)
  REAL(kind=8),intent(IN) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)
  REAL(kind=8),intent(IN) :: cbt(nrp,nfs),cbp(nrp,nfs)

!-- Output:
  REAL(kind=8),intent(OUT) :: fkinLMr(l_max+1)
  REAL(kind=8),intent(OUT) :: fconvLMr(l_max+1)
  REAL(kind=8),intent(OUT) :: fviscLMr(l_max+1)
  REAL(kind=8),intent(OUT) :: fresLMr(l_maxMag+1),fpoynLMr(l_maxMag+1)

!-- Local:
  INTEGER :: nTheta,nThetaB,nThetaNHS
  INTEGER :: nPhi
  REAL(kind=8) :: fkinAS(nfs),fconvAS(nfs),fkin,fconv,phiNorm
  REAL(kind=8) :: fviscAS(nfs),fvisc
  REAL(kind=8) :: fpoynAS(nfs),fresAS(nfs),fpoyn,fres

!-- End of declaration
!---------------------------------------------------------------------------

  phiNorm=2.D0*pi/DBLE(n_phi_max)


  nTheta=nThetaStart-1
  DO nThetaB=1,sizeThetaB
     nTheta=nTheta+1
     nThetaNHS=(nTheta+1)/2
     fkinAS(nThetaB) =0.D0
     fconvAS(nThetaB)=0.D0
     fviscAS(nThetaB)=0.D0
     fkin=0.D0
     fconv=0.D0
     fvisc=0.D0
     DO nPhi=1,n_phi_max
        fconv=temp0(nr)*vr(nPhi,nThetaB)*sr(nPhi,nThetaB)     +    &
                   ViscHeatFac*orho1(nr)*vr(nPhi,nThetaB)*pr(nPhi,nThetaB)

        fkin=0.5D0*or2(nR)*orho2(nR)*(osn2(nThetaNHS)*(     &
                           vt(nPhi,nThetaB)*vt(nPhi,nThetaB)  +    &
                           vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )+    &
                   or2(nR)*vr(nPhi,nThetaB)*vr(nPhi,nThetaB) )*    &
             vr(nPhi,nThetaB)

        IF ( nR/=n_r_icb .AND. nR/=n_r_cmb ) THEN
           fvisc=-2.D0*visc(nR)*orho1(nR)*vr(nPhi,nThetaB)*or2(nR)* ( &
                                             dvrdr(nPhi,nThetaB)   &
             -(2.D0*or1(nR)+2.D0/3.D0*beta(nR))*vr(nPhi,nThetaB) )-&
                   visc(nR)*orho1(nR)*vt(nPhi,nThetaB)*            &
                                        osn2(nThetaNHS)* (  &
                                   or2(nR)*dvrdt(nPhi,nThetaB)     &
                                          +dvtdr(nPhi,nThetaB)     &
                     -(2.d0*or1(nR)+beta(nR))*vt(nPhi,nThetaB) )  -&
                   visc(nR)*orho1(nR)*vp(nPhi,nThetaB)*    &
                                           osn2(nThetaNHS)* (  &
                                   or2(nR)*dvrdp(nPhi,nThetaB)     &
                                          +dvpdr(nPhi,nThetaB)     &
                     -(2.d0*or1(nR)+beta(nR))*vp(nPhi,nThetaB) ) 
        END IF

        fkinAS(nThetaB) = fkinAS(nThetaB)+fkin
        fconvAS(nThetaB)=fconvAS(nThetaB)+fconv
        fviscAS(nThetaB)=fviscAS(nThetaB)+fvisc
     END DO
     fkinAS(nThetaB) =phiNorm* fkinAS(nThetaB)
     fconvAS(nThetaB)=phiNorm*fconvAS(nThetaB)
     fviscAS(nThetaB)=phiNorm*fviscAS(nThetaB)
  END DO

  IF ( l_mag_nl) THEN
     nTheta=nThetaStart-1
     DO nThetaB=1,sizeThetaB
        nTheta=nTheta+1
        nThetaNHS=(nTheta+1)/2
        fresAS(nThetaB) =0.D0
        fpoynAS(nThetaB)=0.D0
        fres=0.D0
        fpoyn=0.D0
        DO nPhi=1,n_phi_max
            fres =osn2(nThetaNHS)*(                              &
                           cbt(nPhi,nThetaB)*bp(nPhi,nThetaB)  - &
                           cbp(nPhi,nThetaB)*bt(nPhi,nThetaB) )

            fpoyn=-orho1(nR)*or2(nR)*osn2(nThetaNHS)*(                        &
                        vp(nPhi,nThetaB)*br(nPhi,nThetaB)*bp(nPhi,nThetaB)  - &
                        vr(nPhi,nThetaB)*bp(nPhi,nThetaB)*bp(nPhi,nThetaB)  - &
                        vr(nPhi,nThetaB)*bt(nPhi,nThetaB)*bt(nPhi,nThetaB)  + &
                        vt(nPhi,nThetaB)*br(nPhi,nThetaB)*bt(nPhi,nThetaB) )

            fresAS(nThetaB) = fresAS(nThetaB)+fres
            fpoynAS(nThetaB)=fpoynAS(nThetaB)+fpoyn
        END DO
        fresAS(nThetaB) =phiNorm* fresAS(nThetaB)
        fpoynAS(nThetaB)=phiNorm*fpoynAS(nThetaB)
     END DO
     CALL legtfAS2(fresLMr,fpoynLMr,fresAS,fpoynAS,l_max+1,nThetaStart,sizeThetaB)
  END IF

!------ Add contribution from thetas in block:
  CALL legtfAS(fviscLMr,fviscAS,l_max+1,nThetaStart,sizeThetaB)
  CALL legtfAS2(fconvLMr,fkinLMr,fconvAS,fkinAS,l_max+1,nThetaStart,sizeThetaB)

  RETURN
END SUBROUTINE get_fluxes
!-----------------------------------------------------------------------------
