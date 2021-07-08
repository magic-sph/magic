module nl_special_calc
   !
   ! This module allows to calculcate several diagnostics that need to be
   ! computed in the physical space (non-linear quantities)
   !

   use precision_mod
   use truncation, only: n_phi_max, l_max, l_maxMag, n_theta_max
   use constants, only: pi, one, two, third, half
   use logic, only: l_mag_nl, l_anelastic_liquid
   use physical_parameters, only: ek, ViscHeatFac, ThExpNb
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: orho1, orho2, or2, or1, beta, temp0, &
       &                       visc, or4, r, alpha0
   use horizontal_data, only: O_sin_theta_E2, cosTheta, sn2, osn2, cosn2, gauss

   implicit none

   private

   public :: get_nlBLayers, get_perpPar, get_fluxes, get_helicity, &
        &    get_visc_heat, get_ekin_solid_liquid

contains

   subroutine get_nlBLayers(vt,vp,dvtdr,dvpdr,dsdr,dsdt,dsdp,&
              &             uhAS,duhAS,gradsAS,nR)
      !
      !   Calculates axisymmetric contributions of:
      !
      !     * the horizontal velocity :math:`u_h = \sqrt{u_\theta^2+u_\phi^2}`
      !     * its radial derivative :math:`|\partial u_h/\partial r|`
      !     * The thermal dissipation rate :math:`(\nabla T)^2`
      !
      !   This subroutine is used when one wants to evaluate viscous and thermal
      !   dissipation layers
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vt(:,:),vp(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dsdr(:,:),dsdt(:,:),dsdp(:,:)

      !-- Output variables:
      real(cp), intent(out) :: uhAS
      real(cp), intent(out) :: duhAS
      real(cp), intent(out) :: gradsAS

      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS
      real(cp) :: uh,duh,phiNorm,grads

      phiNorm=one/real(n_phi_max,cp)
      uhAS   =0.0_cp
      duhAS  =0.0_cp
      gradsAS=0.0_cp

      !--- Horizontal velocity uh and duh/dr + (grad T)**2
      !$omp parallel do default(shared)                       &
      !$omp& private(nTheta, nThetaNHS, nPhi, uh, duh, grads) &
      !$omp& reduction(+:uhAS,duhAS,gradsAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nThetaNHS=(nTheta+1)/2
            uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(     &
            &             vt(nTheta,nPhi)*vt(nTheta,nPhi)+    &
            &             vp(nTheta,nPhi)*vp(nTheta,nPhi)  )
            duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(          &
            &                    dvtdr(nTheta,nPhi)*vt(nTheta,nPhi)-&
            &    (or1(nR)+beta(nR))*vt(nTheta,nPhi)*vt(nTheta,nPhi)+&
            &                    dvpdr(nTheta,nPhi)*vp(nTheta,nPhi)-&
            &    (or1(nR)+beta(nR))*vp(nTheta,nPhi)*vp(nTheta,nPhi) )

            grads =  dsdr(nTheta,nPhi)*dsdr(nTheta,nPhi)          &
            &      +or2(nR)*O_sin_theta_E2(nTheta)*(              &
            &              dsdt(nTheta,nPhi)*dsdt(nTheta,nPhi)    &
            &             +dsdp(nTheta,nPhi)*dsdp(nTheta,nPhi) )

            uhAS=uhAS+phiNorm*gauss(nThetaNHS)*sqrt(uh)
            if (uh /= 0.0_cp) then
               duhAS=duhAS+phiNorm*gauss(nThetaNHS)*abs(duh)/sqrt(uh)
            end if
            gradsAS=gradsAS+phiNorm*gauss(nThetaNHS)*grads
         end do
      end do
      !$omp end parallel do

   end subroutine get_nlBLayers
!------------------------------------------------------------------------------
   subroutine get_perpPar(vr,vt,vp,EperpAS,EparAS,EperpaxiAS,EparaxiAS,nR)
      !
      !   Calculates the energies parallel and perpendicular to the rotation axis
      !
      !     * :math:`E_\perp = 0.5 (v_s^2+v_\phi^2)` with
      !       :math:`v_s= v_r\sin\theta+v_\theta\cos\theta`
      !     * :math:`E_\parallel  = 0.5v_z^2` with
      !       :math:`v_z= v_r\cos\theta-v_\theta*\sin\theta`
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)

      !-- Output variables:
      real(cp), intent(out) :: EperpAS,EparAS,EperpaxiAS,EparaxiAS

      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS
      real(cp) :: vras(n_theta_max),vtas(n_theta_max),vpas(n_theta_max),phiNorm
      real(cp) :: Eperp,Epar,Eperpaxi,Eparaxi

      phiNorm=one/real(n_phi_max,cp)
      EperpAS   =0.0_cp
      EparAS    =0.0_cp
      EperpaxiAS=0.0_cp
      EparaxiAS =0.0_cp
      vras(:)   =0.0_cp
      vtas(:)   =0.0_cp
      vpas(:)   =0.0_cp

      do nPhi=1,n_phi_max
         vras(:)=vras(:)+vr(1:n_theta_max,nPhi)
         vtas(:)=vtas(:)+vt(1:n_theta_max,nPhi)
         vpas(:)=vpas(:)+vp(1:n_theta_max,nPhi)
      end do
      vras(:)=vras(:)*phiNorm
      vtas(:)=vtas(:)*phiNorm
      vpas(:)=vpas(:)*phiNorm

      !$omp parallel do default(shared)                 &
      !$omp& private(nTheta,nPhi,nThetaNHS)             &
      !$omp& private(Eperp, Epar, Eperpaxi, Eparaxi)    &
      !$omp& reduction(+:EparAS,EperpAS,EparaxiAS,EperpaxiAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nThetaNHS=(nTheta+1)/2

            Eperp=half*or2(nR)*orho2(nR)*(                                         &
            &       or2(nR)*sn2(nThetaNHS)*      vr(nTheta,nPhi)*vr(nTheta,nPhi) + &
            &       (osn2(nThetaNHS)-one)*       vt(nTheta,nPhi)*vt(nTheta,nPhi) + &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nTheta,nPhi)*vt(nTheta,nPhi) + &
            &       osn2(nThetaNHS)*             vp(nTheta,nPhi)*vp(nTheta,nPhi) )

            Epar =half*or2(nR)*orho2(nR)*(                                         &
            &       or2(nR)*(one-sn2(nThetaNHS))*vr(nTheta,nPhi)*vr(nTheta,nPhi) + &
            &                                    vt(nTheta,nPhi)*vt(nTheta,nPhi) - &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nTheta,nPhi)*vt(nTheta,nPhi) )

            Eperpaxi=half*or2(nR)*orho2(nR)*(                                  &
            &         or2(nR)*sn2(nThetaNHS)*      vras(nTheta)*vras(nTheta) + &
            &         (osn2(nThetaNHS)-one)*       vtas(nTheta)*vtas(nTheta) + &
            &         two*or1(nR)*cosTheta(nTheta)*vras(nTheta)*vtas(nTheta) + &
            &         osn2(nThetaNHS)*             vpas(nTheta)*vpas(nTheta) )

            Eparaxi =half*or2(nR)*orho2(nR)*(                                  &
            &         or2(nR)*(one-sn2(nThetaNHS))*vras(nTheta)*vras(nTheta) + &
            &                                      vtas(nTheta)*vtas(nTheta) - &
            &         two*or1(nR)*cosTheta(nTheta)*vras(nTheta)*vtas(nTheta) )

            EperpAS   =   EperpAS+phiNorm*gauss(nThetaNHS)*Eperp
            EparAS    =    EparAS+phiNorm*gauss(nThetaNHS)*Epar
            EperpaxiAS=EperpaxiAS+phiNorm*gauss(nThetaNHS)*Eperpaxi
            EparaxiAS = EparaxiAS+phiNorm*gauss(nThetaNHS)*Eparaxi
         end do
      end do
      !$omp end parallel do

   end subroutine get_perpPar
!------------------------------------------------------------------------------
   subroutine get_fluxes(vr,vt,vp,dvrdr,dvtdr,dvpdr,dvrdt,dvrdp,sr,pr,br,bt, &
              &          bp,cbt,cbp,fconvAS,fkinAS,fviscAS,fpoynAS,fresAS,nR)
      !
      !   Calculates the fluxes:
      !
      !     * Convective flux: :math:`F_c= \rho T (u_r s)`
      !     * Kinetic flux: :math:`F_k = 1/2\,\rho u_r (u_r^2+u_\theta^2+u_\phi^2)`
      !     * Viscous flux: :math:`F_= -(u \cdot S )_r`)
      !
      !   If the run is magnetic, then this routine also computes:
      !
      !     * Poynting flux
      !     * resistive flux
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: dvrdr(:,:),dvtdr(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvrdt(:,:),dvrdp(:,:)
      real(cp), intent(in) :: sr(:,:),pr(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: cbt(:,:),cbp(:,:)

      !-- Output variables:
      real(cp), intent(out) :: fkinAS,fconvAS,fviscAS
      real(cp), intent(out) :: fresAS,fpoynAS

      !-- Local variables:
      integer :: nTheta,nThetaNHS,nPhi
      real(cp) :: fkin,fconv,phiNorm,fvisc,fpoyn,fres

      phiNorm=two*pi/real(n_phi_max,cp)

      fkinAS =0.0_cp
      fconvAS=0.0_cp
      fviscAS=0.0_cp
      fvisc  =0.0_cp
      !$omp parallel do default(shared)                             &
      !$omp& private(nTheta, nPhi, nThetaNHS, fconv, fkin, fvisc)   &
      !$omp& reduction(+:fkinAS,fconvAS,fviscAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nThetaNHS=(nTheta+1)/2
            if ( l_anelastic_liquid ) then
               fconv=vr(nTheta,nPhi)*sr(nTheta,nPhi)
            else
               fconv=temp0(nr)*vr(nTheta,nPhi)*sr(nTheta,nPhi)     +    &
               &          ViscHeatFac*ThExpNb*alpha0(nr)*temp0(nr)*     &
               &          orho1(nr)*vr(nTheta,nPhi)*pr(nTheta,nPhi)
            end if

            fkin=half*or2(nR)*orho2(nR)*(osn2(nThetaNHS)*(           &
            &                  vt(nTheta,nPhi)*vt(nTheta,nPhi)  +    &
            &                  vp(nTheta,nPhi)*vp(nTheta,nPhi) )+    &
            &          or2(nR)*vr(nTheta,nPhi)*vr(nTheta,nPhi) )*    &
            &                             vr(nTheta,nPhi)

            if ( nR/=n_r_icb .and. nR/=n_r_cmb ) then
               fvisc=-two*visc(nR)*orho1(nR)*vr(nTheta,nPhi)*or2(nR)* (     &
               &                             dvrdr(nTheta,nPhi)             &
               & -(two*or1(nR)+two*third*beta(nR))*vr(nTheta,nPhi) )-       &
               &                       visc(nR)*orho1(nR)*vt(nTheta,nPhi)*  &
               &                            osn2(nThetaNHS)* (              &
               &                       or2(nR)*dvrdt(nTheta,nPhi)           &
               &                              +dvtdr(nTheta,nPhi)           &
               &       -(two*or1(nR)+beta(nR))*vt(nTheta,nPhi) )  -         &
               &       visc(nR)*orho1(nR)*vp(nTheta,nPhi)*                  &
               &                               osn2(nThetaNHS)* (           &
               &                       or2(nR)*dvrdp(nTheta,nPhi)           &
               &                              +dvpdr(nTheta,nPhi)           &
               &       -(two*or1(nR)+beta(nR))*vp(nTheta,nPhi) )
            end if

            fkinAS = fkinAS+phiNorm*gauss(nThetaNHS)*fkin
            fconvAS=fconvAS+phiNorm*gauss(nThetaNHS)*fconv
            fviscAS=fviscAS+phiNorm*gauss(nThetaNHS)*fvisc
         end do
      end do
      !$omp end parallel do

      if ( l_mag_nl) then
         fresAS =0.0_cp
         fpoynAS=0.0_cp
         !$omp parallel do default(shared)                      &
         !$omp& private(nTheta, nPhi, nThetaNHS, fres, fpoyn)   &
         !$omp& reduction(+:fresAS,fpoynAS)
         do nPhi=1,n_phi_max
            do nTheta=1,n_theta_max
               nThetaNHS=(nTheta+1)/2
               fres =osn2(nThetaNHS)*(                            &
               &              cbt(nTheta,nPhi)*bp(nTheta,nPhi)  - &
               &              cbp(nTheta,nPhi)*bt(nTheta,nPhi) )

               fpoyn=-orho1(nR)*or2(nR)*osn2(nThetaNHS)*(                     &
               &           vp(nTheta,nPhi)*br(nTheta,nPhi)*bp(nTheta,nPhi)  - &
               &           vr(nTheta,nPhi)*bp(nTheta,nPhi)*bp(nTheta,nPhi)  - &
               &           vr(nTheta,nPhi)*bt(nTheta,nPhi)*bt(nTheta,nPhi)  + &
               &           vt(nTheta,nPhi)*br(nTheta,nPhi)*bt(nTheta,nPhi) )

               fresAS = fresAS+phiNorm*gauss(nThetaNHS)*fres
               fpoynAS=fpoynAS+phiNorm*gauss(nThetaNHS)*fpoyn
            end do
         end do
         !$omp end parallel do
      end if

   end subroutine get_fluxes
!------------------------------------------------------------------------------
   subroutine get_helicity(vr,vt,vp,cvr,dvrdt,dvrdp,dvtdr,dvpdr,HelAS, &
              &            Hel2AS,HelnaAS,Helna2AS,HelEAAS,nR)
      !
      !   Calculates axisymmetric contributions of helicity HelLMr and
      !   helicity**2  Hel2LMr in (l,m=0,r) space.
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: cvr(:,:),dvrdt(:,:),dvrdp(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvpdr(:,:)

      !-- Output variables:
      real(cp), intent(out) :: HelAS(2)
      real(cp), intent(out) :: Hel2AS(2)
      real(cp), intent(out) :: HelnaAS(2)
      real(cp), intent(out) :: Helna2AS(2)
      real(cp), intent(out) :: HelEAAS

      !-- Local variables:
      integer :: nTheta,nThetaNHS,nPhi
      real(cp) :: Helna,Hel,phiNorm
      real(cp) :: vrna,vtna,vpna,cvrna,dvrdtna,dvrdpna,dvtdrna,dvpdrna
      real(cp) :: vras(n_theta_max),vtas(n_theta_max),vpas(n_theta_max)
      real(cp) :: cvras(n_theta_max),dvrdtas(n_theta_max),dvrdpas(n_theta_max)
      real(cp) :: dvtdras(n_theta_max),dvpdras(n_theta_max)

      !-- Remark: 2pi not used the normalization below
      !-- this is why we have a 2pi factor after radial integration
      !-- in the subroutine outHelicity()
      phiNorm=one/real(n_phi_max,cp)
      HelAS(:)   =0.0_cp
      Hel2AS(:)  =0.0_cp
      HelnaAS(:) =0.0_cp
      Helna2AS(:)=0.0_cp
      HelEAAS    =0.0_cp

      vras(:)   =0.0_cp
      cvras(:)  =0.0_cp
      vtas(:)   =0.0_cp
      vpas(:)   =0.0_cp
      dvrdpas(:)=0.0_cp
      dvpdras(:)=0.0_cp
      dvtdras(:)=0.0_cp
      dvrdtas(:)=0.0_cp
      do nPhi=1,n_phi_max
         vras(:)   =vras(:)   +   vr(1:n_theta_max,nPhi)
         cvras(:)  =cvras(:)  +  cvr(1:n_theta_max,nPhi)
         vtas(:)   =vtas(:)   +   vt(1:n_theta_max,nPhi)
         vpas(:)   =vpas(:)   +   vp(1:n_theta_max,nPhi)
         dvrdpas(:)=dvrdpas(:)+dvrdp(1:n_theta_max,nPhi)
         dvpdras(:)=dvpdras(:)+dvpdr(1:n_theta_max,nPhi)
         dvtdras(:)=dvtdras(:)+dvtdr(1:n_theta_max,nPhi)
         dvrdtas(:)=dvrdtas(:)+dvrdt(1:n_theta_max,nPhi)
      end do
      vras(:)   =vras(:)   *phiNorm
      cvras(:)  =cvras(:)  *phiNorm
      vtas(:)   =vtas(:)   *phiNorm
      vpas(:)   =vpas(:)   *phiNorm
      dvrdpas(:)=dvrdpas(:)*phiNorm
      dvpdras(:)=dvpdras(:)*phiNorm
      dvtdras(:)=dvtdras(:)*phiNorm
      dvrdtas(:)=dvrdtas(:)*phiNorm

      !--- Helicity:
      !$omp parallel do default(shared)                     &
      !$omp& private(nTheta, nThetaNHS, nPhi, Hel, Helna)   &
      !$omp& private(vrna, cvrna, vtna, vpna)               &
      !$omp& private(dvrdpna, dvpdrna, dvtdrna, dvrdtna)    &
      !$omp& reduction(+:HelAS,Hel2AS,HelnaAS,Helna2AS,HelEAAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nThetaNHS = (nTheta+1)/2
            vrna   =   vr(nTheta,nPhi)-vras(nTheta)
            cvrna  =  cvr(nTheta,nPhi)-cvras(nTheta)
            vtna   =   vt(nTheta,nPhi)-vtas(nTheta)
            vpna   =   vp(nTheta,nPhi)-vpas(nTheta)
            dvrdpna=dvrdp(nTheta,nPhi)-dvrdpas(nTheta)
            dvpdrna=dvpdr(nTheta,nPhi)-beta(nR)*vp(nTheta,nPhi) &
            &       -dvpdras(nTheta)+beta(nR)*vpas(nTheta)
            dvtdrna=dvtdr(nTheta,nPhi)-beta(nR)*vt(nTheta,nPhi) &
            &       -dvtdras(nTheta)+beta(nR)*vtas(nTheta)
            dvrdtna=dvrdt(nTheta,nPhi)-dvrdtas(nTheta)
            Hel=or4(nR)*orho2(nR)*vr(nTheta,nPhi)*cvr(nTheta,nPhi) +  &
            &             or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                                       vt(nTheta,nPhi) * &
            &                          ( or2(nR)*dvrdp(nTheta,nPhi) - &
            &                                    dvpdr(nTheta,nPhi) + &
            &                         beta(nR)*   vp(nTheta,nPhi) ) + &
            &                                       vp(nTheta,nPhi) * &
            &                          (         dvtdr(nTheta,nPhi) - &
            &                           beta(nR)*   vt(nTheta,nPhi) - &
            &                            or2(nR)*dvrdt(nTheta,nPhi) ) )
            Helna=                      or4(nR)*orho2(nR)*vrna*cvrna + &
            &              or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                       vtna*( or2(nR)*dvrdpna-dvpdrna ) + &
            &                       vpna*( dvtdrna-or2(nR)*dvrdtna ) )

            if ( mod(nTheta,2)  == 1 ) then ! Northern Hemisphere
               HelAS(1)   =HelAS(1) +phiNorm*gauss(nThetaNHS)*Hel
               Hel2AS(1)  =Hel2AS(1)+phiNorm*gauss(nThetaNHS)*Hel*Hel
               HelnaAS(1) =HelnaAS(1) +phiNorm*gauss(nThetaNHS)*Helna
               Helna2AS(1)=Helna2AS(1)+phiNorm*gauss(nThetaNHS)*Helna*Helna
               HelEAAS    =HelEAAS +phiNorm*gauss(nThetaNHS)*Hel
            else
               HelAS(2)   =HelAS(2) +phiNorm*gauss(nThetaNHS)*Hel
               Hel2AS(2)  =Hel2AS(2)+phiNorm*gauss(nThetaNHS)*Hel*Hel
               HelnaAS(2) =HelnaAS(2) +phiNorm*gauss(nThetaNHS)*Helna
               Helna2AS(2)=Helna2AS(2)+phiNorm*gauss(nThetaNHS)*Helna*Helna
               HelEAAS    =HelEAAS -phiNorm*gauss(nThetaNHS)*Hel
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine get_helicity
!------------------------------------------------------------------------------
   subroutine get_visc_heat(vr,vt,vp,cvr,dvrdr,dvrdt,dvrdp,dvtdr,&
              &             dvtdp,dvpdr,dvpdp,viscAS,nR)
      !
      !   Calculates axisymmetric contributions of the viscous heating
      !
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:),cvr(:,:)
      real(cp), intent(in) :: dvrdr(:,:),dvrdt(:,:),dvrdp(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvtdp(:,:)
      real(cp), intent(in) :: dvpdr(:,:),dvpdp(:,:)

      !-- Output variables:
      real(cp), intent(out) :: viscAS

      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS
      real(cp) :: vischeat,csn2,phiNorm

      phiNorm=two*pi/real(n_phi_max,cp)
      viscAS=0.0_cp

      !$omp parallel do default(shared)                   &
      !$omp& private(nTheta,nThetaNHS,csn2,nPhi,vischeat) &
      !$omp& reduction(+:viscAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nThetaNHS=(nTheta+1)/2
            csn2     =cosn2(nThetaNHS)
            if ( mod(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

            vischeat=       or2(nR)*orho1(nR)*visc(nR)*(        &
            &     two*(                    dvrdr(nTheta,nPhi) - & ! (1)
            &    (two*or1(nR)+beta(nR))*vr(nTheta,nPhi) )**2  + &
            &     two*( csn2*                 vt(nTheta,nPhi) + &
            &                              dvpdp(nTheta,nPhi) + &
            &                              dvrdr(nTheta,nPhi) - & ! (2)
            &     or1(nR)*              vr(nTheta,nPhi) )**2  + &
            &     two*(                    dvpdp(nTheta,nPhi) + &
            &           csn2*                 vt(nTheta,nPhi) + & ! (3)
            &     or1(nR)*              vr(nTheta,nPhi) )**2  + &
            &          ( two*              dvtdp(nTheta,nPhi) + &
            &                                cvr(nTheta,nPhi) - & ! (6)
            &      two*csn2*            vp(nTheta,nPhi) )**2  + &
            &                               osn2(nThetaNHS) * ( &
            &         ( r(nR)*             dvtdr(nTheta,nPhi) - &
            &          (two+beta(nR)*r(nR))*  vt(nTheta,nPhi) + & ! (4)
            &     or1(nR)*           dvrdt(nTheta,nPhi) )**2  + &
            &         ( r(nR)*             dvpdr(nTheta,nPhi) - &
            &          (two+beta(nR)*r(nR))*  vp(nTheta,nPhi) + & ! (5)
            &    or1(nR)*            dvrdp(nTheta,nPhi) )**2 )- &
            &    two*third*(  beta(nR)*     vr(nTheta,nPhi) )**2 )

            viscAS=viscAS+phiNorm*gauss(nThetaNHS)*viscHeat
         end do
      end do
      !$omp end parallel do

   end subroutine get_visc_heat
!------------------------------------------------------------------------------
   subroutine get_ekin_solid_liquid(vr,vt,vp,phi,ekinS,ekinL,nR)
      !
      ! This subroutine computes the kinetic energy content in the solid
      ! and in the liquid phase when phase field is employed.
      !

      !-- Input variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:),phi(:,:)

      !-- Output variables:
      real(cp), intent(out) :: ekinS ! Kinetic energy in the solid phase
      real(cp), intent(out) :: ekinL ! Kinetic energy in the liquid phase

      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS
      real(cp) :: vischeat,phiNorm,ekin

      phiNorm=two*pi/real(n_phi_max,cp)
      ekinL=0.0_cp
      ekinS=0.0_cp

      !$omp parallel do default(shared)            &
      !$omp& private(nTheta, nThetaNHS, nPhi,ekin) &
      !$omp& reduction(+:ekinS,ekinL)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nThetaNHS=(nTheta+1)/2

            ekin = half*orho1(nR)*(                                      &
            &          or2(nR)*        vr(nTheta,nPhi)*vr(nTheta,nPhi) + &
            &          osn2(nThetaNHS)*vt(nTheta,nPhi)*vt(nTheta,nPhi) + &
            &          osn2(nThetaNHS)*vp(nTheta,nPhi)*vp(nTheta,nPhi) )

            if ( phi(nTheta,nPhi) >= half ) then
               ekinS=ekinS+phiNorm*gauss(nThetaNHS)*ekin
            else
               ekinL=ekinL+phiNorm*gauss(nThetaNHS)*ekin
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine get_ekin_solid_liquid
!------------------------------------------------------------------------------
end module nl_special_calc
