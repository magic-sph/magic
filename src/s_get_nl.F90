!$Id$
!***********************************************************************
  SUBROUTINE get_nl(vr,vt,vp,dvrdr,dvtdr,dvpdr,cvr, &
       &                 dvrdt,dvrdp,dvtdp,dvpdp, &
       &                 br,bt,bp,cbr,cbt,cbp,sr, &
       &                 Advr,Advt,Advp,LFr,LFt,LFp, &
       &                 VSr,VSt,VSp,VxBr,VxBt,VxBp, &
       &                 ViscHeat,OhmLoss, &
       &                 nR,nBc,nThetaStart)
    !***********************************************************************

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    !-----------------------------------------------------------------------

    !  calculates non-linear products in grid-space for radial
    !  level nR and returns them in arrays wnlr1-3, snlr1-3, bnlr1-3

    !  if nBc >0 velocities are zero only the (vxB)
    !  contributions to bnlr2-3 need to be calculated

    !  vr...sr: (input) velocity, magnetic field comp. and derivs, entropy
    !                   on grid points
    !  nR: (input) radial level
    !  i1: (input) range of points in theta for which calculation is done

    !-----------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  IMPLICIT NONE

    !-- Input of variables:
    INTEGER,INTENT(IN) :: nR
    INTEGER,INTENT(IN) :: nBc
    INTEGER,INTENT(IN) :: nThetaStart

    !---- Fields components:
    REAL(kind=8),INTENT(IN) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
    REAL(kind=8),INTENT(IN) :: dvrdr(nrp,nfs),dvtdr(nrp,nfs),dvpdr(nrp,nfs)
    REAL(kind=8),INTENT(IN) :: cvr(nrp,nfs),dvrdt(nrp,nfs),dvrdp(nrp,nfs)
    REAL(kind=8),INTENT(IN) :: dvtdp(nrp,nfs),dvpdp(nrp,nfs)
    REAL(kind=8),INTENT(IN) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)
    REAL(kind=8),INTENT(IN) :: cbr(nrp,nfs),cbt(nrp,nfs),cbp(nrp,nfs)
    REAL(kind=8),INTENT(IN) :: sr(nrp,nfs)

    !-- Output:
    REAL(kind=8),INTENT(OUT) :: Advr(nrp,nfs),Advt(nrp,nfs),Advp(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: LFr(nrp,nfs),LFt(nrp,nfs),LFp(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: VSr(nrp,nfs),VSt(nrp,nfs),VSp(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: VxBr(nrp,nfs),VxBt(nrp,nfs),VxBp(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: ViscHeat(nrp,nfs), OhmLoss(nrp,nfs)

    !-- Local:
    INTEGER :: nTheta
    INTEGER :: nThetaLast,nThetaB,nThetaNHS
    INTEGER :: nPhi
    REAL(kind=8) :: or2sn2,or4sn2,csn2

    !-- end of declaration
    !---------------------------------------------------------------------------


    nThetaLast=nThetaStart-1

    IF ( l_mag_LF .AND. nBc == 0 ) THEN
       !------ Get the Lorentz force:
       nTheta=nThetaLast
       DO nThetaB=1,sizeThetaB

          nTheta   =nTheta+1
          nThetaNHS=(nTheta+1)/2
          or4sn2   =or4(nR)*osn2(nThetaNHS)

          DO nPhi=1,n_phi_max
             !--------- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
             LFr(nPhi,nThetaB)=  LFfac*osn2(nThetaNHS) * ( &
                  cbt(nPhi,nThetaB)*bp(nPhi,nThetaB) - &
                  cbp(nPhi,nThetaB)*bt(nPhi,nThetaB) )
          END DO
          LFr(n_phi_max+1,nThetaB)=0.D0
          LFr(n_phi_max+2,nThetaB)=0.D0

          !--------- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
          DO nPhi=1,n_phi_max
             LFt(nPhi,nThetaB)=           LFfac*or4sn2 * ( &
                  cbp(nPhi,nThetaB)*br(nPhi,nThetaB) - &
                  cbr(nPhi,nThetaB)*bp(nPhi,nThetaB) )
          END DO
          LFt(n_phi_max+1,nThetaB)=0.D0
          LFt(n_phi_max+2,nThetaB)=0.D0

          !--------- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
          DO nPhi=1,n_phi_max
             LFp(nPhi,nThetaB)=           LFfac*or4sn2 * ( &
                  cbr(nPhi,nThetaB)*bt(nPhi,nThetaB) - &
                  cbt(nPhi,nThetaB)*br(nPhi,nThetaB) )
          END DO
          LFp(n_phi_max+1,nThetaB)=0.D0
          LFp(n_phi_max+2,nThetaB)=0.D0

       END DO   ! theta loop
    END IF      ! Lorentz force required ?

    IF ( l_conv_nl .AND. nBc == 0 ) THEN

       !------ Get Advection:
       nTheta=nThetaLast
       DO nThetaB=1,sizeThetaB ! loop over theta points in block
          nTheta   =nTheta+1
          nThetaNHS=(nTheta+1)/2
          or4sn2   =or4(nR)*osn2(nThetaNHS)
          csn2     =cosn2(nThetaNHS)
          IF ( MOD(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

          DO nPhi=1,n_phi_max
             Advr(nPhi,nThetaB)=      -or2(nR)*orho1(nR) * ( &
                  vr(nPhi,nThetaB) * &
                  (       dvrdr(nPhi,nThetaB) - &
                  ( 2.D0*or1(nR)+beta(nR) )*vr(nPhi,nThetaB) ) + &
                  osn2(nThetaNHS) * ( &
                  vt(nPhi,nThetaB) * &
                  (       dvrdt(nPhi,nThetaB) - &
                  r(nR)*      vt(nPhi,nThetaB) ) + &
                  vp(nPhi,nThetaB) * &
                  (       dvrdp(nPhi,nThetaB) - &
                  r(nR)*      vp(nPhi,nThetaB) ) ) )
          END DO
          Advr(n_phi_max+1,nThetaB)=0.D0
          Advr(n_phi_max+2,nThetaB)=0.D0
          DO nPhi=1,n_phi_max
             Advt(nPhi,nThetaB)=         or4sn2*orho1(nR) * ( &
                  -vr(nPhi,nThetaB) * &
                  (   dvtdr(nPhi,nThetaB) - &
                  beta(nR)*vt(nPhi,nThetaB) )   + &
                  vt(nPhi,nThetaB) * &
                  ( csn2*vt(nPhi,nThetaB) + &
                  dvpdp(nPhi,nThetaB) + &
                  dvrdr(nPhi,nThetaB) )   + &
                  vp(nPhi,nThetaB) * &
                  ( csn2*vp(nPhi,nThetaB) - &
                  dvtdp(nPhi,nThetaB) )  )
          END DO
          Advt(n_phi_max+1,nThetaB)=0.D0
          Advt(n_phi_max+2,nThetaB)=0.D0
          DO nPhi=1,n_phi_max
             Advp(nPhi,nThetaB)=         or4sn2*orho1(nR) * ( &
                  -vr(nPhi,nThetaB) * &
                  ( dvpdr(nPhi,nThetaB) - &
                  beta(nR)*vp(nPhi,nThetaB) )   - &
                  vt(nPhi,nThetaB) * &
                  ( dvtdp(nPhi,nThetaB) + &
                  cvr(nPhi,nThetaB) )   - &
                  vp(nPhi,nThetaB) * dvpdp(nPhi,nThetaB) )
          END DO
          Advp(n_phi_max+1,nThetaB)=0.D0
          Advp(n_phi_max+2,nThetaB)=0.D0
       END DO ! theta loop

    END IF  ! Navier-Stokes nonlinear advection term ?

    IF ( l_heat_nl .AND. nBc == 0 ) THEN
       !------ Get V S, the divergence of the is entropy advection:
       nTheta=nThetaLast
       DO nThetaB=1,sizeThetaB
          nTheta   =nTheta+1
          nThetaNHS=(nTheta+1)/2
          or2sn2=or2(nR)*osn2(nThetaNHS)
          DO nPhi=1,n_phi_max     ! calculate v*s components
             VSr(nPhi,nThetaB)= &
                  vr(nPhi,nThetaB)*sr(nPhi,nThetaB)
             VSt(nPhi,nThetaB)= &
                  or2sn2*vt(nPhi,nThetaB)*sr(nPhi,nThetaB)
             VSp(nPhi,nThetaB)= &
                  or2sn2*vp(nPhi,nThetaB)*sr(nPhi,nThetaB)
          END DO
          VSr(n_phi_max+1,nThetaB)=0.D0
          VSr(n_phi_max+2,nThetaB)=0.D0
          VSt(n_phi_max+1,nThetaB)=0.D0
          VSt(n_phi_max+2,nThetaB)=0.D0
          VSp(n_phi_max+1,nThetaB)=0.D0
          VSp(n_phi_max+2,nThetaB)=0.D0
       END DO  ! theta loop
    END IF     ! heat equation required ?

    IF ( l_mag_nl ) THEN

       IF ( nBc == 0 ) THEN

          !------ Get (V x B) , the curl of this is the dynamo term:
          nTheta=nThetaLast
          DO nThetaB=1,sizeThetaB
             nTheta   =nTheta+1
             nThetaNHS=(nTheta+1)/2
             or4sn2=or4(nR)*osn2(nThetaNHS)

             DO nPhi=1,n_phi_max
                VxBr(nPhi,nThetaB)=  orho1(nR)*osn2(nThetaNHS) * ( &
                     vt(nPhi,nThetaB)*bp(nPhi,nThetaB) - &
                     vp(nPhi,nThetaB)*bt(nPhi,nThetaB) )
             END DO
             VxBr(n_phi_max+1,nThetaB)=0.D0
             VxBr(n_phi_max+2,nThetaB)=0.D0

             DO nPhi=1,n_phi_max
                VxBt(nPhi,nThetaB)=  orho1(nR)*or4sn2 * ( &
                     vp(nPhi,nThetaB)*br(nPhi,nThetaB) - &
                     vr(nPhi,nThetaB)*bp(nPhi,nThetaB) )
             END DO
             VxBt(n_phi_max+1,nThetaB)=0.D0
             VxBt(n_phi_max+2,nThetaB)=0.D0

             DO nPhi=1,n_phi_max
                VxBp(nPhi,nThetaB)=   orho1(nR)*or4sn2 * ( &
                     vr(nPhi,nThetaB)*bt(nPhi,nThetaB) - &
                     vt(nPhi,nThetaB)*br(nPhi,nThetaB) )
             END DO
             VxBp(n_phi_max+1,nThetaB)=0.D0
             VxBp(n_phi_max+2,nThetaB)=0.D0
          END DO   ! theta loop

       ELSE IF ( nBc == 1 ) THEN ! stress free boundary

          nTheta=nThetaLast
          DO nThetaB=1,sizeThetaB
             nTheta   =nTheta+1
             nThetaNHS=(nTheta+1)/2
             or4sn2   =or4(nR)*osn2(nThetaNHS)
             DO nPhi=1,n_phi_max
                VxBt(nPhi,nThetaB)= or4sn2 * orho1(nR) * &
                     vp(nPhi,nThetaB)*br(nPhi,nThetaB)
                VxBp(nPhi,nThetaB)= - or4sn2 * orho1(nR) * &
                     vt(nPhi,nThetaB)*br(nPhi,nThetaB)
             END DO
             VxBt(n_phi_max+1,nThetaB)=0.D0
             VxBt(n_phi_max+2,nThetaB)=0.D0
             VxBp(n_phi_max+1,nThetaB)=0.D0
             VxBp(n_phi_max+2,nThetaB)=0.D0
          END DO

       ELSE IF ( nBc == 2 ) THEN  ! rigid boundary :

          !----- Only vp.ne.0 at boundary allowed (rotation of boundaries about z-axis):
          nTheta=nThetaLast
          DO nThetaB=1,sizeThetaB
             nTheta   =nTheta+1
             nThetaNHS=(nTheta+1)/2
             or4sn2   =or4(nR)*osn2(nThetaNHS)
             DO nPhi=1,n_phi_max
                VxBt(nPhi,nThetaB)= or4sn2 * orho1(nR) * &
                     vp(nPhi,nThetaB)*br(nPhi,nThetaB)
                VxBp(nPhi,nThetaB)= 0.D0
             END DO
             VxBt(n_phi_max+1,nThetaB)=0.D0
             VxBt(n_phi_max+2,nThetaB)=0.D0
             VxBp(n_phi_max+1,nThetaB)=0.D0
             VxBp(n_phi_max+2,nThetaB)=0.D0
          END DO

       END IF  ! boundary ?

    END IF ! l_mag_nl ?

    IF ( l_anel .AND. nBc == 0 ) THEN
       !------ Get viscous heating
       nTheta=nThetaLast
       DO nThetaB=1,sizeThetaB ! loop over theta points in block
          nTheta   =nTheta+1
          nThetaNHS=(nTheta+1)/2
          csn2     =cosn2(nThetaNHS)
          IF ( MOD(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

          DO nPhi=1,n_phi_max
             ViscHeat(nPhi,nThetaB)= or4(nR)*                 &
                  orho1(nR)*otemp1(nR)*visc(nR)*( &
                  2.D0*(                     dvrdr(nPhi,nThetaB) - & ! (1)
                  (2.D0*or1(nR)+beta(nR))*vr(nphi,nThetaB) )**2  + &
                  2.D0*( csn2*                  vt(nPhi,nThetaB) + &
                  dvpdp(nphi,nThetaB) + &
                  dvrdr(nPhi,nThetaB) - & ! (2)
                  or1(nR)*               vr(nPhi,nThetaB) )**2  + &
                  2.D0*(                     dvpdp(nphi,nThetaB) + &
                  csn2*                  vt(nPhi,nThetaB) + & ! (3)
                  or1(nR)*               vr(nPhi,nThetaB) )**2  + &
                  ( 2.D0*               dvtdp(nPhi,nThetaB) + &
                  cvr(nPhi,nThetaB) - & ! (6)
                  2.D0*csn2*             vp(nPhi,nThetaB) )**2  + &
                  osn2(nThetaNHS) * ( &
                  ( r(nR)*              dvtdr(nPhi,nThetaB) - &
                  (2.D0+beta(nR)*r(nR))*  vt(nPhi,nThetaB) + & ! (4)
                  or1(nR)*            dvrdt(nPhi,nThetaB) )**2  + &
                  ( r(nR)*              dvpdr(nPhi,nThetaB) - &
                  (2.D0+beta(nR)*r(nR))*  vp(nPhi,nThetaB) + & ! (5)
                  or1(nR)*            dvrdp(nPhi,nThetaB) )**2 )- &
                  2.D0/3.D0*(  beta(nR)*        vr(nPhi,nThetaB) )**2 )
          END DO
          ViscHeat(n_phi_max+1,nThetaB)=0.D0
          ViscHeat(n_phi_max+2,nThetaB)=0.D0
       END DO ! theta loop

       IF (l_mag_nl) THEN
          !------ Get ohmic losses
          nTheta=nThetaLast
          DO nThetaB=1,sizeThetaB ! loop over theta points in block
             nTheta   =nTheta+1
             nThetaNHS=(nTheta+1)/2
             DO nPhi=1,n_phi_max
                OhmLoss(nPhi,nThetaB)=or2(nR)*otemp1(nR)*lambda(nR)* &
                     ( or2(nR)*                cbr(nPhi,nThetaB)**2 + &
                     osn2(nThetaNHS)*        cbt(nPhi,nThetaB)**2 + &
                     osn2(nThetaNHS)*        cbp(nPhi,nThetaB)**2  )
             END DO
             OhmLoss(n_phi_max+1,nThetaB)=0.D0
             OhmLoss(n_phi_max+2,nThetaB)=0.D0
          END DO ! theta loop

       END IF ! if l_mag_nl ?

    END IF  ! Viscous heating and Ohmic losses ?


    RETURN
  END SUBROUTINE get_nl

!-----------------------------------------------------------------------------
