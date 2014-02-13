! Id: s_getHelLM.f 358 2013-02-08 12:32:19Z dannert $
!***********************************************************************
    SUBROUTINE getHelLM(vr,vt,vp,cvr,dvrdt,dvrdp,dvtdr,dvpdr, &
                        HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,nR,nThetaStart)
!***********************************************************************

!-----------------------------------------------------------------------

!   Calculates axisymmetric contributions of helicity HelLMr and
!   helicity**2  Hel2LMr in (l,m=0,r) space.

!-----------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data

    IMPLICIT NONE

!-- Input of variables
    INTEGER,intent(IN) :: nR
    INTEGER,intent(IN) :: nThetaStart

!-- Input of fields components:
    REAL(kind=8),intent(IN) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
    REAL(kind=8),intent(IN) :: cvr(nrp,nfs)
    REAL(kind=8),intent(IN) :: dvrdt(nrp,nfs),dvrdp(nrp,nfs)
    REAL(kind=8),intent(IN) :: dvtdr(nrp,nfs),dvpdr(nrp,nfs)

!-- Output:
    REAL(kind=8),intent(INOUT) :: HelLMr(l_max+1)
    REAL(kind=8),intent(INOUT) :: Hel2LMr(l_max+1)
    REAL(kind=8),intent(INOUT) :: HelnaLMr(l_max+1)
    REAL(kind=8),intent(INOUT) :: Helna2LMr(l_max+1)

!-- Local:
    INTEGER :: nTheta,nThetaB
    INTEGER :: nPhi,l
    REAL(kind=8) :: Helna,HelAS(nfs),Hel2AS(nfs)
    REAL(kind=8) :: Hel,HelnaAS(nfs),Helna2AS(nfs),phiNorm
    REAL(kind=8) :: vras,vtas,vpas,cvras,dvrdtas,dvrdpas,dvtdras,dvpdras
    REAL(kind=8) :: vrna,vtna,vpna,cvrna,dvrdtna,dvrdpna,dvtdrna,dvpdrna


!-- End of declaration
!---------------------------------------------------------------------------

    phiNorm=1.D0/DBLE(n_phi_max)

    IF ( nThetaStart == 1 ) THEN
    !------ Zero lm coeffs for first theta block:
        DO l=1,l_max+1
            HelLMr(l) =0.D0
            Hel2LMr(l)=0.D0
            HelnaLMr(l) =0.D0
            Helna2LMr(l)=0.D0
        END DO
    END IF

!--- Helicity:
    nTheta=nThetaStart-1
    DO nThetaB=1,sizeThetaB
        nTheta=nTheta+1
        HelAS(nThetaB) =0.D0
        Hel2AS(nThetaB)=0.D0
        vras=0.D0
        cvras=0.D0
        vtas=0.D0
        vpas=0.D0
        dvrdpas=0.D0
        dvpdras=0.D0
        dvtdras=0.D0
        dvrdtas=0.D0
        DO nPhi=1,n_phi_max
            vras=vras+vr(nPhi,nThetaB)
            cvras=cvras+cvr(nPhi,nThetaB)
            vtas=vtas+vt(nPhi,nThetaB)
            vpas=vpas+vp(nPhi,nThetaB)
            dvrdpas=dvrdpas+dvrdp(nPhi,nThetaB)
            dvpdras=dvpdras+dvpdr(nPhi,nThetaB)
            dvtdras=dvtdras+dvtdr(nPhi,nThetaB)
            dvrdtas=dvrdtas+dvrdt(nPhi,nThetaB)
        END DO
        vras=vras*phiNorm
        cvras=cvras*phiNorm
        vtas=vtas*phiNorm
        vpas=vpas*phiNorm
        dvrdpas=dvrdpas*phiNorm
        dvpdras=dvpdras*phiNorm
        dvtdras=dvtdras*phiNorm
        dvrdtas=dvrdtas*phiNorm
        DO nPhi=1,n_phi_max
            vrna   =   vr(nPhi,nThetaB)-vras
            cvrna  =  cvr(nPhi,nThetaB)-cvras
            vtna   =   vt(nPhi,nThetaB)-vtas
            vpna   =   vp(nPhi,nThetaB)-vpas
            dvrdpna=dvrdp(nPhi,nThetaB)-dvrdpas
            dvpdrna=dvpdr(nPhi,nThetaB)-beta(nR)*vp(nPhi,nThetaB) &
                    -dvpdras+beta(nR)*vpas
            dvtdrna=dvtdr(nPhi,nThetaB)-beta(nR)*vt(nPhi,nThetaB) &
                    -dvtdras+beta(nR)*vtas
            dvrdtna=dvrdt(nPhi,nThetaB)-dvrdtas
            Hel=or4(nR)*orho2(nR)*vr(nPhi,nThetaB)*cvr(nPhi,nThetaB) + &
                           or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
                                                    vt(nPhi,nThetaB) * &
                                       ( or2(nR)*dvrdp(nPhi,nThetaB) - &
                                                 dvpdr(nPhi,nThetaB) + &
                                      beta(nR)*   vp(nPhi,nThetaB) ) + &
                                                    vp(nPhi,nThetaB) * &
                                       (         dvtdr(nPhi,nThetaB) - &
                                        beta(nR)*   vt(nPhi,nThetaB) - &
                                         or2(nR)*dvrdt(nPhi,nThetaB) ) )
            Helna=                      or4(nR)*orho2(nR)*vrna*cvrna + &
                           or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
                                    vtna*( or2(nR)*dvrdpna-dvpdrna ) + &
                                    vpna*( dvtdrna-or2(nR)*dvrdtna ) )

            HelAS(nThetaB)   =HelAS(nThetaB) +Hel
            Hel2AS(nThetaB)  =Hel2AS(nThetaB)+Hel*Hel
            HelnaAS(nThetaB) =HelAS(nThetaB) +Helna
            Helna2AS(nThetaB)=Hel2AS(nThetaB)+Helna*Helna
        END DO
        HelAS(nThetaB) =phiNorm*HelAS(nThetaB)
        Hel2AS(nThetaB)=phiNorm*Hel2AS(nThetaB)
        HelnaAS(nThetaB) =phiNorm*HelnaAS(nThetaB)
        Helna2AS(nThetaB)=phiNorm*Helna2AS(nThetaB)
    END DO

!------ Add contribution from thetas in block:
    CALL legtfAS2(HelLMr,Hel2LMr,HelAS,Hel2AS, &
                  l_max+1,nThetaStart,sizeThetaB)
    CALL legtfAS2(HelnaLMr,Helna2LMr,HelnaAS,Helna2AS, &
                  l_max+1,nThetaStart,sizeThetaB)

    RETURN
    end SUBROUTINE getHelLM

!-----------------------------------------------------------------------------
