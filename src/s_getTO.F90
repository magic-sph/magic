!$Id$
!***********************************************************************
    subroutine getTO(vr,vt,vp,cvr,dvpdr,br,bt,bp,cbr,cbt, &
                                    BsLast,BpLast,BzLast, &
                        dzRstrLM,dzAstrLM,dzCorLM,dzLFLM, &
                   dtLast,nR,nThetaStart,nThetaBlockSize)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!-----------------------------------------------------------------------
!  This program calculates various axisymmetric linear
!  and nonlinear variables for a radial grid point nR and
!  a theta-block.
!  Input are the fields vr,vt,vp,cvr,dvpdr
!  Output are linear azimuthally averaged  field VpAS (flow phi component),
!  VpAS2 (square of flow phi component), V2AS (V*V),
!  and Coriolis force Cor. These are give in (r,theta)-space.
!  Also in (r,theta)-space are azimuthally averaged correlations of
!  non-axisymmetric flow components and the respective squares:
!  Vsp=Vs*Vp,Vzp,Vsz,VspC,VzpC,VszC. These are used to calulcate
!  the respective correlations and Reynolds stress.
!  In addition three output field are given in (lm,r) space:
!    dzRstrLMr,dzAstrLMr,dzCorLM,dzLFLM.
!  These are used to calculate the total Reynolds stress,
!  advection and viscous stress later. Their calculation
!  retraces the calculations done in the time-stepping part
!  of the code.
!-----------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE torsional_oscillations
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data

    IMPLICIT NONE

!-- Input of variables
    REAL(kind=8) :: dtLast              ! last time step
    INTEGER :: nR                 ! radial grid point
    INTEGER :: nThetaStart        ! theta block
    INTEGER :: nThetaBlockSize

!-- Input of field components:
    REAL(kind=8) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
    REAL(kind=8) :: cvr(nrp,nfs),dvpdr(nrp,nfs)
    REAL(kind=8) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)
    REAL(kind=8) :: cbr(nrp,nfs),cbt(nrp,nfs)
!-- Input of magnetic field components from last time step
!   stored in getTOnext
    REAL(kind=8) :: BsLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
    REAL(kind=8) :: BpLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
    REAL(kind=8) :: BzLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)


!-- Output of arrays needing further treatment in s_getTOfinish.f:
    REAL(kind=8) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
    REAL(kind=8) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

!-- Output via common block to s_outTO.f
! include 'c_TO.f'


!-- Local stuff:
    INTEGER :: nTheta,nThetaBlock
    INTEGER :: nPhi
    REAL(kind=8) :: VrMean,VtMean,VpMean
    REAL(kind=8) :: Vr2Mean,Vt2Mean,Vp2Mean
    REAL(kind=8) :: LFmean
    REAL(kind=8) :: cvrMean,dvpdrMean
    REAL(kind=8) :: VrdVpdrMean,VtcVrMean
    REAL(kind=8) :: Bs2Mean,BszMean
    REAL(kind=8) :: BspMean,BpzMean
    REAL(kind=8) :: BspdMean,BpsdMean
    REAL(kind=8) :: BzpdMean,BpzdMean

    REAL(kind=8) :: Rmean(nfs),Amean(nfs)
    REAL(kind=8) :: dzCorMean(nfs),dzLFmean(nfs)

    REAL(kind=8) :: sinT,Osin,Osin2,cosT,cosOsin2
    REAL(kind=8) :: phiNorm

    REAL(kind=8) :: BsL,BzL,BpL
    REAL(kind=8) :: Bs2F1,Bs2F2,Bs2F3,BspF1,BspF2
    REAL(kind=8) :: BpzF1,BpzF2,BszF1,BszF2,BszF3
    REAL(kind=8) :: BsF1,BsF2,BpF1,BzF1,BzF2

!-- End of declaration
!---------------------------------------------------------------------------

    IF ( lVerbose ) WRITE(*,*) '! Starting getTO!'

    phiNorm=1.D0/DBLE(n_phi_max)

!-- Big loop over thetas in block:
    nTheta=nThetaStart-1
    DO nThetaBlock=1,nThetaBlockSize
        nTheta=nTheta+1
        sinT =sinTheta(nTheta)
        cosT =cosTheta(nTheta)
        Osin =1.D0/sinT
        Osin2=Osin*Osin
        cosOsin2=cosT*Osin2
        Bs2F1=sinT*sinT*or4(nR)
        Bs2F2=cosT*cosT*Osin2*or2(nR)
        Bs2F3=2.D0*cosT*or3(nR)
        BspF1=or3(nR)
        BspF2=cosT*Osin2*or2(nR)
        BpzF1=cosT*Osin*or3(nR)
        BpzF2=or2(nR)*Osin
        BszF1=sinT*cosT*or4(nR)
        BszF2=(2*cosT*cosT-1.D0)*Osin*or3(nR)
        BszF3=cosT*Osin*or2(nR)
        BsF1 =sinT*or2(nR)
        BsF2 =cosT*Osin*or1(nR)
        BpF1 =Osin*or1(nR)
        BzF1 =cosT*or2(nR)
        BzF2 =or1(nR)

    !--- Get zonal means of velocity and derivatives:
        VrMean     =0.D0
        VtMean     =0.D0
        VpMean     =0.D0
        Vr2Mean    =0.D0
        Vt2Mean    =0.D0
        Vp2Mean    =0.D0
        dVpdrMean  =0.D0
        cVrMean    =0.D0
        VrdVpdrMean=0.D0
        VtcVrMean  =0.D0
        LFmean     =0.D0
        Bs2Mean    =0.D0
        BspMean    =0.D0
        BpzMean    =0.D0
        BszMean    =0.D0
        BspdMean   =0.D0
        BpsdMean   =0.D0
        BzpdMean   =0.D0
        BpzdMean   =0.D0
        DO nPhi=1,n_phi_max
            VrMean =VrMean +vr(nPhi,nThetaBlock)
            VtMean =VtMean +vt(nPhi,nThetaBlock)
            VpMean =VpMean +vp(nPhi,nThetaBlock)
            Vr2Mean=Vr2Mean+vr(nPhi,nThetaBlock)* &
                            vr(nPhi,nThetaBlock)
            Vt2Mean=Vt2Mean+vt(nPhi,nThetaBlock)* &
                            vt(nPhi,nThetaBlock)
            Vp2Mean=Vp2Mean+vp(nPhi,nThetaBlock)* &
                            vp(nPhi,nThetaBlock)
            dVpdrMean  =dVpdrMean    + orho1(nR)* & ! dvp/dr
                      ( dvpdr(nPhi,nThetaBlock) - &
                    beta(nR)*vp(nPhi,nThetaBlock) )
            cVrMean    =cVrMean  +cvr(nPhi,nThetaBlock)
            VrdVpdrMean=VrdVpdrMean  + orho1(nR)*           & ! rho * vr * dvp/dr
              vr(nPhi,nThetaBlock)*(dvpdr(nPhi,nThetaBlock) &
                             -beta(nR)*vp(nPhi,nThetaBlock))
            VtcVrMean  =VtcVrMean    + orho1(nR)*  &  ! rho * vt * cvr
                        vt(nPhi,nThetaBlock)*cvr(nPhi,nThetaBlock)
            IF ( l_mag ) THEN
                LFmean=LFmean +                                      &
                    cbr(nPhi,nThetaBlock)*bt(nPhi,nThetaBlock) -     &
                    cbt(nPhi,nThetaBlock)*br(nPhi,nThetaBlock)
                Bs2Mean=Bs2Mean +                                    &
                   Bs2F1*br(nPhi,nThetaBlock)*br(nPhi,nThetaBlock) + &
                   Bs2F2*bt(nPhi,nThetaBlock)*bt(nPhi,nThetaBlock) + &
                   Bs2F3*br(nPhi,nThetaBlock)*bt(nPhi,nThetaBlock)
                BspMean=BspMean +                                    &
                   BspF1*br(nPhi,nThetaBlock)*bp(nPhi,nThetaBlock) + &
                   BspF2*bt(nPhi,nThetaBlock)*bp(nPhi,nThetaBlock)
                BpzMean=BpzMean +                                    &
                   BpzF1*br(nPhi,nThetaBlock)*bp(nPhi,nThetaBlock) - &
                   BpzF2*bt(nPhi,nThetaBlock)*bp(nPhi,nThetaBlock)
                BszMean=BszMean +                                    &
                   BszF1*br(nPhi,nThetaBlock)*br(nPhi,nThetaBlock) + &
                   BszF2*br(nPhi,nThetaBlock)*bt(nPhi,nThetaBlock) - &
                   BszF3*bt(nPhi,nThetaBlock)*bt(nPhi,nThetaBlock)
                BsL=BsF1*br(nPhi,nThetaBlock) + BsF2*bt(nPhi,nThetaBlock)
                BpL=BpF1*bp(nPhi,nThetaBlock)
                BzL=BzF1*br(nPhi,nThetaBlock) - BzF2*bt(nPhi,nThetaBlock)
                BspdMean=BspdMean+BsL*(BpL-BpLast(nPhi,nTheta,nR))
                BpsdMean=BpsdMean+BpL*(BsL-BsLast(nPhi,nTheta,nR))
                BzpdMean=BzpdMean+BzL*(BpL-BpLast(nPhi,nTheta,nR))
                BpzdMean=BpzdMean+BpL*(BzL-BzLast(nPhi,nTheta,nR))
            END IF
         END DO

        Vr2Mean=phiNorm*or4(nR)*Vr2Mean
        Vt2Mean=phiNorm*or2(nR)*Osin2*Vt2Mean
        Vp2Mean=phiNorm*or2(nR)*Osin2*Vp2Mean
        IF ( nR == n_r_CMB ) THEN
            VrMean=0.D0
            Vr2Mean=0.D0
            IF ( ktopv == 2 ) THEN
                VtMean=0.D0
                Vt2Mean=0.D0
                VpMean=0.D0
                Vp2Mean=0.D0
            END IF
        END IF
        IF ( nR == n_r_CMB ) THEN
            VrMean=0.D0
            Vr2Mean=0.D0
            IF ( kbotv == 2 ) THEN
                VtMean=0.D0
                Vt2Mean=0.D0
                VpMean=0.D0
                Vp2Mean=0.D0
            END IF
        END IF
        V2AS(nTheta,nR)=Vr2Mean+Vt2Mean+Vp2Mean
        VpMean =phiNorm*or1(nR)*Osin*VpMean
    !--- This is Coriolis force / r*sin(theta)
        dzCorMean(nThetaBlock)= phiNorm*2.D0*CorFac * &
             (or3(nR)*VrMean+or2(nR)*cosOsin2*VtMean)
        IF ( l_mag ) THEN
        !--- This is Lorentz force/ r*sin(theta)
            dzLFmean(nThetaBlock)=phiNorm*or4(nR)*Osin2*LFmean
            Bs2AS(nTheta,nR) =phiNorm*Bs2Mean
            BspAS(nTheta,nR) =phiNorm*BspMean
            BpzAS(nTheta,nR) =phiNorm*BpzMean
            BszAS(nTheta,nR) =phiNorm*BszMean
            BspdAS(nTheta,nR)=phiNorm*(BspdMean/dtLast)
            BpsdAS(nTheta,nR)=phiNorm*(BpsdMean/dtLast)
            BzpdAS(nTheta,nR)=phiNorm*(BzpdMean/dtLast)
            BpzdAS(nTheta,nR)=phiNorm*(BpzdMean/dtLast)
        END IF

    ! dVpdrMean, VtcVrMean and VrdVpdrMean are already divided by rho
        Rmean(nThetaBlock)=                                  &
                     phiNorm * or4(nR)*Osin2* (VrdVpdrMean - &
                                  phiNorm*VrMean*dVpdrMean + &
                                                 VtcVrMean - &
                          orho1(nR)*phiNorm*VtMean*cVrMean )
        Amean(nThetaBlock)=                                  &
           phiNorm*or4(nR)*Osin2*phiNorm*(VrMean*dVpdrMean + &
                                  orho1(nR)*VtMean*cVrMean )

    END DO ! Loop over thetas in block !

!--- Transform and Add to LM-Space:

!------ Add contribution from thetas in block:
!       note legtfAS2 returns modes l=0 -- l=l_max+1
    CALL legtfAS2(dzRstrLM,dzAstrLM,Rmean,Amean,     &
                  l_max+2,nThetaStart,nThetaBlockSize)
    CALL legtfAS2(dzCorLM,dzLFLM,dzCorMean,dZLFmean, &
                  l_max+2,nThetaStart,nThetaBlockSize)


    IF ( lVerbose ) WRITE(*,*) '! End of getTO!'

    RETURN
    end subroutine getTO
!-----------------------------------------------------------------------------
