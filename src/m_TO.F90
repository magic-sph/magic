!$Id$
module torsional_oscillations
!----------------------------------------------------------------------
!  This module contains information for TO calculation and output
!----------------------------------------------------------------------

   use truncation, only: nrp, n_phi_maxStr, n_r_maxStr, l_max, &
                         n_theta_maxStr
   use radial_data, only: n_r_cmb
   use radial_functions, only: r, or1, or2, or3, or4, beta, orho1, &
                               dbeta
   use physical_parameters, only: CorFac, kbotv, ktopv
   use blocking, only: nfs, lm2
   use horizontal_data, only: sinTheta, cosTheta, hdif_V, dTheta1A, dTheta1S, & 
                              dLh
   use logic, only: lVerbose, l_mag
   use legendre_grid_to_spec, only: legTFAS2

   implicit none

   private
 
   real(kind=8), public, allocatable :: dzStrLMr(:,:)
   real(kind=8), public, allocatable :: dzRstrLMr(:,:)
   real(kind=8), public, allocatable :: dzAstrLMr(:,:)
   real(kind=8), public, allocatable :: dzCorLMr(:,:)
   real(kind=8), public, allocatable :: dzLFLMr(:,:)
   real(kind=8), public, allocatable :: dzdVpLMr(:,:)
   real(kind=8), public, allocatable :: dzddVpLMr(:,:)
   real(kind=8), public, allocatable :: V2AS(:,:)
   real(kind=8), public, allocatable :: Bs2AS(:,:)
   real(kind=8), public, allocatable :: BszAS(:,:)
   real(kind=8), public, allocatable :: BspAS(:,:)
   real(kind=8), public, allocatable :: BpzAS(:,:)
   real(kind=8), public, allocatable :: BspdAS(:,:)
   real(kind=8), public, allocatable :: BpsdAS(:,:)
   real(kind=8), public, allocatable :: BzpdAS(:,:)
   real(kind=8), public, allocatable :: BpzdAS(:,:)
   real(kind=8), public, allocatable :: ddzASL(:,:)

   public :: initialize_TO, getTO, getTOnext, getTOfinish

contains

   subroutine initialize_TO

      allocate( dzStrLMr(l_max+1,n_r_maxStr) )
      allocate( dzRstrLMr(l_max+1,n_r_maxStr) )
      allocate( dzAstrLMr(l_max+1,n_r_maxStr) )
      allocate( dzCorLMr(l_max+1,n_r_maxStr) )
      allocate( dzLFLMr(l_max+1,n_r_maxStr) )
      allocate( dzdVpLMr(l_max+1,n_r_maxStr) )
      allocate( dzddVpLMr(l_max+1,n_r_maxStr) )
      allocate( V2AS(n_theta_maxStr,n_r_maxStr) )
      allocate( Bs2AS(n_theta_maxStr,n_r_maxStr) )
      allocate( BszAS(n_theta_maxStr,n_r_maxStr) )
      allocate( BspAS(n_theta_maxStr,n_r_maxStr) )
      allocate( BpzAS(n_theta_maxStr,n_r_maxStr) )
      allocate( BspdAS(n_theta_maxStr,n_r_maxStr) )
      allocate( BpsdAS(n_theta_maxStr,n_r_maxStr) )
      allocate( BzpdAS(n_theta_maxStr,n_r_maxStr) )
      allocate( BpzdAS(n_theta_maxStr,n_r_maxStr) )
      allocate( ddzASL(l_max+1,n_r_maxStr) ) 

   end subroutine initialize_TO
!-----------------------------------------------------------------------------
   subroutine getTO(vr,vt,vp,cvr,dvpdr,br,bt,bp,cbr,cbt,  &
                                    BsLast,BpLast,BzLast, &
                        dzRstrLM,dzAstrLM,dzCorLM,dzLFLM, &
                   dtLast,nR,nThetaStart,nThetaBlockSize)
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

      !-- Input of variables
      real(kind=8), intent(in) :: dtLast              ! last time step
      integer,      intent(in) :: nR                 ! radial grid point
      integer,      intent(in) :: nThetaStart        ! theta block
      integer,      intent(in) :: nThetaBlockSize
      real(kind=8), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
      real(kind=8), intent(in) :: cvr(nrp,nfs),dvpdr(nrp,nfs)
      real(kind=8), intent(in) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)
      real(kind=8), intent(in) :: cbr(nrp,nfs),cbt(nrp,nfs)
      real(kind=8), intent(in) :: BsLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
      real(kind=8), intent(in) :: BpLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
      real(kind=8), intent(in) :: BzLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)

      !-- Output of arrays needing further treatment in s_getTOfinish.f:
      real(kind=8), intent(out) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
      real(kind=8), intent(out) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

      !-- Local variables:
      integer :: nTheta,nThetaBlock
      integer :: nPhi
      real(kind=8) :: VrMean,VtMean,VpMean
      real(kind=8) :: Vr2Mean,Vt2Mean,Vp2Mean
      real(kind=8) :: LFmean
      real(kind=8) :: cvrMean,dvpdrMean
      real(kind=8) :: VrdVpdrMean,VtcVrMean
      real(kind=8) :: Bs2Mean,BszMean
      real(kind=8) :: BspMean,BpzMean
      real(kind=8) :: BspdMean,BpsdMean
      real(kind=8) :: BzpdMean,BpzdMean
      real(kind=8) :: Rmean(nfs),Amean(nfs)
      real(kind=8) :: dzCorMean(nfs),dzLFmean(nfs)
      real(kind=8) :: sinT,Osin,Osin2,cosT,cosOsin2
      real(kind=8) :: phiNorm

      real(kind=8) :: BsL,BzL,BpL
      real(kind=8) :: Bs2F1,Bs2F2,Bs2F3,BspF1,BspF2
      real(kind=8) :: BpzF1,BpzF2,BszF1,BszF2,BszF3
      real(kind=8) :: BsF1,BsF2,BpF1,BzF1,BzF2

      if ( lVerbose ) write(*,*) '! Starting getTO!'

      phiNorm=1.D0/dble(n_phi_maxStr)

      !-- Big loop over thetas in block:
      nTheta=nThetaStart-1
      do nThetaBlock=1,nThetaBlockSize
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
         do nPhi=1,n_phi_maxStr
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
            if ( l_mag ) then
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
            end if
         end do

         Vr2Mean=phiNorm*or4(nR)*Vr2Mean
         Vt2Mean=phiNorm*or2(nR)*Osin2*Vt2Mean
         Vp2Mean=phiNorm*or2(nR)*Osin2*Vp2Mean
         if ( nR == n_r_CMB ) then
            VrMean=0.D0
            Vr2Mean=0.D0
            if ( ktopv == 2 ) then
               VtMean=0.D0
               Vt2Mean=0.D0
               VpMean=0.D0
               Vp2Mean=0.D0
            end if
         end if
         if ( nR == n_r_CMB ) then
            VrMean=0.D0
            Vr2Mean=0.D0
            if ( kbotv == 2 ) then
               VtMean=0.D0
               Vt2Mean=0.D0
               VpMean=0.D0
               Vp2Mean=0.D0
            end if
         end if
         V2AS(nTheta,nR)=Vr2Mean+Vt2Mean+Vp2Mean
         VpMean =phiNorm*or1(nR)*Osin*VpMean
         !--- This is Coriolis force / r*sin(theta)
         dzCorMean(nThetaBlock)= phiNorm*2.D0*CorFac * &
               (or3(nR)*VrMean+or2(nR)*cosOsin2*VtMean)
         if ( l_mag ) then
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
         end if

         ! dVpdrMean, VtcVrMean and VrdVpdrMean are already divided by rho
         Rmean(nThetaBlock)=                                  &
                      phiNorm * or4(nR)*Osin2* (VrdVpdrMean - &
                                   phiNorm*VrMean*dVpdrMean + &
                                                  VtcVrMean - &
                            orho1(nR)*phiNorm*VtMean*cVrMean )
         Amean(nThetaBlock)=                                  &
            phiNorm*or4(nR)*Osin2*phiNorm*(VrMean*dVpdrMean + &
                                   orho1(nR)*VtMean*cVrMean )

      end do ! Loop over thetas in block !

      !--- Transform and Add to LM-Space:
      !------ Add contribution from thetas in block:
      !       note legtfAS2 returns modes l=0 -- l=l_max+1
      call legTFAS2(dzRstrLM,dzAstrLM,Rmean,Amean,     &
                    l_max+2,nThetaStart,nThetaBlockSize)
      call legTFAS2(dzCorLM,dzLFLM,dzCorMean,dZLFmean, &
                    l_max+2,nThetaStart,nThetaBlockSize)

      if ( lVerbose ) write(*,*) '! End of getTO!'

   end subroutine getTO
!-----------------------------------------------------------------------------
   subroutine getTOnext(zAS,br,bt,bp,lTONext,lTONext2,    &
                dt,dtLast,nR,nThetaStart,nThetaBlockSize, &
                                    BsLast,BpLast,BzLast)
      !-----------------------------------------------------------------------
      !  Preparing TO calculation by storing flow and magnetic field
      !  contribution to build time derivative.
      !-----------------------------------------------------------------------

      !-- Input of variables:
      real(kind=8), intent(in) :: dt,dtLast
      integer,      intent(in) :: nR
      integer,      intent(in) :: nThetaStart
      integer,      intent(in) :: nThetaBlockSize
      logical,      intent(in) :: lTONext,lTONext2
      real(kind=8), intent(in) :: zAS(l_max+1)
      real(kind=8), intent(in) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)

      !-- Output variables:
      real(kind=8), intent(out) :: BsLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
      real(kind=8), intent(out) :: BpLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
      real(kind=8), intent(out) :: BzLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)

      !-- Local variables:
      integer :: l,lm
      integer :: nTheta,nThetaBlock
      integer :: nPhi

      real(kind=8) :: sinT,cosT
      real(kind=8) :: BsF1,BsF2,BpF1,BzF1,BzF2

      if ( lVerbose ) write(*,*) '! Starting getTOnext!',dtLast

      nTheta=nThetaStart-1

      if ( lTONext2 .and. nThetaStart == 1 ) then

         dzddVpLMr(1,nR)=0.D0
         do l=1,l_max
            lm=lm2(l,0)
            dzddVpLMr(l+1,nR)=zAS(l+1)
         end do

      else if ( lTOnext ) then

         do nThetaBlock=1,nThetaBlockSize
            nTheta=nTheta+1
            sinT =sinTheta(nTheta)
            cosT =cosTheta(nTheta)
            BsF1 =sinT*or2(nR)
            BsF2 =cosT/sinT*or1(nR)
            BpF1 =or1(nR)/sinT
            BzF1 =cosT*or2(nR)
            BzF2 =or1(nR)

            !--- Get zonal means of velocity and derivatives:
            do nPhi=1,n_phi_maxStr
               BsLast(nPhi,nTheta,nR)=BsF1*br(nPhi,nThetaBlock) + &
                                      BsF2*bt(nPhi,nThetaBlock)
               BpLast(nPhi,nTheta,nR)=BpF1*bp(nPhi,nThetaBlock)
               BzLast(nPhi,nTheta,nR)=BzF1*br(nPhi,nThetaBlock) - &
                                      BzF2*bt(nPhi,nThetaBlock)
            end do
                            
         end do ! Loop over thetas in block !
                  
         if ( nThetaStart == 1 ) then
            dzdVpLMr(1,nR) =0.D0
            dzddVpLMr(1,nR)=0.D0
            do l=1,l_max
               lm=lm2(l,0)
               dzdVpLMr(l+1,nR) = zAS(l+1)
               dzddVpLMr(l+1,nR)= ( dzddVpLMr(l+1,nR) - &
                                  ((dtLast+dt)/dt)*zAS(l+1) )/dtLast
            end do
         end if

      end if

      if ( lVerbose ) write(*,*) '! End of getTOnext!'

   end subroutine getTOnext
!-----------------------------------------------------------------------------
   subroutine getTOfinish(nR,dtLast,zAS,dzAS,ddzAS,dzRstrLM, &
                          dzAstrLM,dzCorLM,dzLFLM)
      !-----------------------------------------------------------------------
      !  This program was previously part of s_getTO.f.
      !  It has now been seperated to get it out of the theta-block loop.
      !-----------------------------------------------------------------------

      !-- Input of variables:
      integer,      intent(in) :: nR
      real(kind=8), intent(in) :: dtLast
      real(kind=8), intent(in) :: zAS(l_max+1)
      real(kind=8), intent(in) :: dzAS(l_max+1) ! anelastic
      real(kind=8), intent(in) :: ddzAS(l_max+1)
      real(kind=8), intent(in) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
      real(kind=8), intent(in) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

      !-- Local variables:
      integer :: l,lS,lA,lm

      !------ When all thetas are done calculate viscous stress in LM space:
      dzStrLMr(1,nR) =0.D0
      dzRstrLMr(1,nR)=0.D0
      dzAstrLMr(1,nR)=0.D0
      dzCorLMr(1,nR) =0.D0
      dzLFLMr(1,nR)  =0.D0
      dzdVpLMr(1,nR) =0.D0
      dzddVpLMr(1,nR)=0.D0
      do l=1,l_max
         lS=(l-1)+1
         lA=(l+1)+1
         lm=lm2(l,0)
         dzStrLMr(l+1,nR)= hdif_V(lm) * (                      &
                                                  ddzAS(l+1) - &
                                         beta(nR)* dzAS(l+1) - &
            (dLh(lm)*or2(nR)+dbeta(nR)+2.d0*beta(nR)*or1(nR))* &
                                zAS(l+1) )
      !---- -r**2/(l(l+1)) 1/sin(theta) dtheta sin(theta)**2
      !     minus sign to bring stuff on the RHS of NS equation !
         dzRstrLMr(l+1,nR)=-r(nR)*r(nR)/dLh(lm) * ( &
                        dTheta1S(lm)*dzRstrLM(lS) - &
                        dTheta1A(lm)*dzRstrLM(lA) )
         dzAstrLMr(l+1,nR)=-r(nR)*r(nR)/dLh(lm) * ( &
                        dTheta1S(lm)*dzAstrLM(lS) - &
                        dTheta1A(lm)*dzAstrLM(lA) )
         dzCorLMr(l+1,nR) =-r(nR)*r(nR)/dLh(lm) * ( &
                         dTheta1S(lm)*dzCorLM(lS) - &
                         dTheta1A(lm)*dzCorLM(lA) )
         dzLFLMr(l+1,nR)  = r(nR)*r(nR)/dLh(lm) * ( &
                          dTheta1S(lm)*dzLFLM(lS) - &
                          dTheta1A(lm)*dzLFLM(lA) )
         dzdVpLMr(l+1,nR) =(zAS(l+1)-dzdVpLMr(l+1,nR))/dtLast
         dzddVpLMr(l+1,nR)=(zAS(l+1)/dtLast+dzddVpLMr(l+1,nR))/dtLast
      end do

   end subroutine getTOfinish
!-----------------------------------------------------------------------------
end module torsional_oscillations
