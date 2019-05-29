module nl_special_calc
   !
   ! This module allows to calculcate several diagnostics that need to be
   ! computed in the physical space (non-linear quantities)
   !

   use precision_mod
   use truncation, only: nrp, n_phi_max, l_max, l_maxMag
   use constants, only: pi, one, two, third, half
   use logic, only: l_mag_nl, l_TP_form, l_anelastic_liquid
   use physical_parameters, only: ek, ViscHeatFac, ThExpNb
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: orho1, orho2, or2, or1, beta, temp0, &
       &                       visc, or4, r, alpha0
   use blocking, only: sizeThetaB, nfs
   use horizontal_data, only: O_sin_theta_E2, cosTheta, sn2, osn2, cosn2
#ifdef WITH_SHTNS
   use shtns, only: spat_to_SH_axi
#else
   use legendre_grid_to_spec, only: legTFAS, legTFAS2
#endif

   implicit none

   private

   public :: get_nlBLayers, get_perpPar, get_fluxes, get_helicity, &
        &    get_visc_heat

contains

   subroutine get_nlBLayers(vt,vp,dvtdr,dvpdr,dsdr,dsdt,dsdp,&
              &             uhLMr,duhLMr,gradsLMr,nR,nThetaStart)
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
      integer,  intent(in) :: nThetaStart
      real(cp), intent(in) :: vt(nrp,nfs),vp(nrp,nfs)
      real(cp), intent(in) :: dvtdr(nrp,nfs),dvpdr(nrp,nfs)
      real(cp), intent(in) :: dsdr(nrp,nfs),dsdt(nrp,nfs),dsdp(nrp,nfs)

      !-- Output variables:
      real(cp), intent(out) :: uhLMr(l_max+1)
      real(cp), intent(out) :: duhLMr(l_max+1)
      real(cp), intent(out) :: gradsLMr(l_max+1)

      !-- Local variables:
      integer :: nTheta,nThetaB
      integer :: nPhi
      real(cp) :: uhAS(nfs),duhAS(nfs),gradsAS(nfs),uh,duh,phiNorm,grads

      phiNorm=one/real(n_phi_max,cp)

      !--- Horizontal velocity uh and duh/dr + (grad T)**2
      nTheta=nThetaStart-1
#ifdef WITH_SHTNS
      !$OMP PARALLEL DO default(shared)                     &
      !$OMP& private(nThetaB, nTheta, nPhi)                 &
      !$OMP& private(uh, duh, grads)
#endif
      do nThetaB=1,sizeThetaB
         nTheta=nThetaStart+nThetaB-1
         uhAS(nThetaB) =0.0_cp
         duhAS(nThetaB)=0.0_cp
         gradsAS(nThetaB)=0.0_cp
         do nPhi=1,n_phi_max
            uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(     &
            &             vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+  &
            &             vp(nPhi,nThetaB)*vp(nPhi,nThetaB)  )
            duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(            &
            &                    dvtdr(nPhi,nThetaB)*vt(nPhi,nThetaB)-&
            &    (or1(nR)+beta(nR))*vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+&
            &                    dvpdr(nPhi,nThetaB)*vp(nPhi,nThetaB)-&
            &    (or1(nR)+beta(nR))*vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )

            grads =  dsdr(nPhi,nThetaB)*dsdr(nPhi,nThetaB)          &
            &      +or2(nR)*O_sin_theta_E2(nTheta)*(                &
            &              dsdt(nPhi,nThetaB)*dsdt(nPhi,nThetaB)    &
            &             +dsdp(nPhi,nThetaB)*dsdp(nPhi,nThetaB) )

            uhAS(nThetaB)=uhAS(nThetaB)+sqrt(uh)
            if (uh /= 0.0_cp) then
               duhAS(nThetaB)=duhAS(nThetaB)+abs(duh)/sqrt(uh)
            end if
            gradsAS(nThetaB)=gradsAS(nThetaB)+grads
         end do
         uhAS(nThetaB)=phiNorm*uhAS(nThetaB)
         duhAS(nThetaB)=phiNorm*duhAS(nThetaB)
         gradsAS(nThetaB)=phiNorm*gradsAS(nThetaB)
      end do
#ifdef WITH_SHTNS
      !$OMP END PARALLEL DO
#endif

      !------ Add contribution from thetas in block:
#ifdef WITH_SHTNS
      call spat_to_SH_axi(gradsAS,gradsLMr)
      call spat_to_SH_axi(uhAS,uhLMr)
      call spat_to_SH_axi(duhAS,duhLMr)
#else
      call legTFAS2(uhLMr,duhLMr,uhAS,duhAS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS(gradsLMr,gradsAS,l_max+1,nThetaStart,sizeThetaB)
#endif

   end subroutine get_nlBLayers
!------------------------------------------------------------------------------
   subroutine get_perpPar(vr,vt,vp,EperpLMr,EparLMr,  &
              &           EperpaxiLMr,EparaxiLMr,nR,nThetaStart)
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
      integer,  intent(in) :: nThetaStart
      real(cp), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)

      !-- Output variables:
      real(cp), intent(out) :: EperpLMr(l_max+1),EparLMr(l_max+1)
      real(cp), intent(out) :: EperpaxiLMr(l_max+1),EparaxiLMr(l_max+1)

      !-- Local variables:
      integer :: nTheta,nThetaB,nThetaNHS
      integer :: nPhi
      real(cp) :: vras,vtas,vpas,phiNorm
      real(cp) :: EperpAS(nfs),EparAS(nfs),Eperp,Epar
      real(cp) :: EperpaxiAS(nfs),EparaxiAS(nfs),Eperpaxi,Eparaxi

      phiNorm=one/real(n_phi_max,cp)

      nTheta=nThetaStart-1
#ifdef WITH_SHTNS
      !$OMP PARALLEL DO default(shared)                 &
      !$OMP& private(nThetaB, nTheta, nPhi)             &
      !$OMP& private(Eperp, Epar, Eperpaxi, Eparaxi)    &
      !$OMP& private(vras, vtas, vpas)
#endif
      do nThetaB=1,sizeThetaB
         nTheta=nThetaStart+nThetaB-1
         nThetaNHS=(nTheta+1)/2

         EperpAS(nThetaB)   =0.0_cp
         EparAS(nThetaB)    =0.0_cp
         EperpaxiAS(nThetaB)=0.0_cp
         EparaxiAS(nThetaB) =0.0_cp
         Eperp   =0.0_cp
         Epar    =0.0_cp
         Eperpaxi=0.0_cp
         Eparaxi =0.0_cp
         vras    =0.0_cp
         vtas    =0.0_cp
         vpas    =0.0_cp

         do nPhi=1,n_phi_max
            vras=vras+vr(nPhi,nThetaB)
            vtas=vtas+vt(nPhi,nThetaB)
            vpas=vpas+vp(nPhi,nThetaB)
         end do
         vras=vras*phiNorm
         vtas=vtas*phiNorm
         vpas=vpas*phiNorm

         do nPhi=1,n_phi_max
            Eperp=half*or2(nR)*orho2(nR)*(                                           &
            &       or2(nR)*sn2(nThetaNHS)*      vr(nPhi,nThetaB)*vr(nPhi,nThetaB) + &
            &       (osn2(nThetaNHS)-one)*       vt(nPhi,nThetaB)*vt(nPhi,nThetaB) + &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nPhi,nThetaB)*vt(nPhi,nThetaB) + &
            &       osn2(nThetaNHS)*                vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )

            Epar =half*or2(nR)*orho2(nR)*(                                           &
            &       or2(nR)*(one-sn2(nThetaNHS))*vr(nPhi,nThetaB)*vr(nPhi,nThetaB) + &
            &                                    vt(nPhi,nThetaB)*vt(nPhi,nThetaB) - &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nPhi,nThetaB)*vt(nPhi,nThetaB) )

            Eperpaxi=half*or2(nR)*orho2(nR)*(                  &
            &         or2(nR)*sn2(nThetaNHS)*      vras*vras + &
            &         (osn2(nThetaNHS)-one)*       vtas*vtas + &
            &         two*or1(nR)*cosTheta(nTheta)*vras*vtas + &
            &         osn2(nThetaNHS)*             vpas*vpas )

            Eparaxi =half*or2(nR)*orho2(nR)*(                  &
            &         or2(nR)*(one-sn2(nThetaNHS))*vras*vras + &
            &                                      vtas*vtas - &
            &         two*or1(nR)*cosTheta(nTheta)*vras*vtas )

            EperpAS(nThetaB)   =   EperpAS(nThetaB)+Eperp
            EparAS(nThetaB)    =    EparAS(nThetaB)+Epar
            EperpaxiAS(nThetaB)=EperpaxiAS(nThetaB)+Eperpaxi
            EparaxiAS(nThetaB) = EparaxiAS(nThetaB)+Eparaxi
         end do
         EperpAS(nThetaB)   =phiNorm*   EperpAS(nThetaB)
         EparAS(nThetaB)    =phiNorm*    EparAS(nThetaB)
         EperpaxiAS(nThetaB)=phiNorm*EperpaxiAS(nThetaB)
         EparaxiAS(nThetaB) =phiNorm* EparaxiAS(nThetaB)
      end do
#ifdef WITH_SHTNS
      !$OMP END PARALLEL DO
#endif

      !-- Add contribution from thetas in block:
#ifdef WITH_SHTNS
      call spat_to_SH_axi(EperpAS, EperpLMr)
      call spat_to_SH_axi(EparAS, EparLMr)
      call spat_to_SH_axi(EperpaxiAS, EperpaxiLMr)
      call spat_to_SH_axi(EparaxiAS, EparaxiLMr)
#else
      call legTFAS2(EperpLMr,EparLMr,EperpAS,EparAS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS2(EperpaxiLMr,EparaxiLMr,EperpaxiAS,EparaxiAS,l_max+1, &
           &        nThetaStart,sizeThetaB)
#endif

   end subroutine get_perpPar
!------------------------------------------------------------------------------
   subroutine get_fluxes(vr,vt,vp,dvrdr,dvtdr,dvpdr,        &
              &          dvrdt,dvrdp,sr,pr,br,bt,bp,cbt,cbp,&
              &          fconvLMr,fkinLMr,fviscLMr,fpoynLMr,&
              &          fresLMR,nR,nThetaStart)
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
      integer,  intent(in) :: nThetaStart
      real(cp), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
      real(cp), intent(in) :: dvrdr(nrp,nfs),dvtdr(nrp,nfs),dvpdr(nrp,nfs)
      real(cp), intent(in) :: dvrdt(nrp,nfs),dvrdp(nrp,nfs)
      real(cp), intent(in) :: sr(nrp,nfs),pr(nrp,nfs)
      real(cp), intent(in) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)
      real(cp), intent(in) :: cbt(nrp,nfs),cbp(nrp,nfs)

      !-- Output variables:
      real(cp), intent(out) :: fkinLMr(l_max+1)
      real(cp), intent(out) :: fconvLMr(l_max+1)
      real(cp), intent(out) :: fviscLMr(l_max+1)
      real(cp), intent(out) :: fresLMr(l_maxMag+1),fpoynLMr(l_maxMag+1)

      !-- Local variables:
      integer :: nTheta,nThetaB,nThetaNHS
      integer :: nPhi
      real(cp) :: fkinAS(nfs),fconvAS(nfs),fkin,fconv,phiNorm
      real(cp) :: fviscAS(nfs),fvisc
      real(cp) :: fpoynAS(nfs),fresAS(nfs),fpoyn,fres

      phiNorm=two*pi/real(n_phi_max,cp)

      nTheta=nThetaStart-1
#ifdef WITH_SHTNS
      !$OMP PARALLEL DO default(shared)         &
      !$OMP& private(nThetaB, nTheta, nPhi)     &
      !$OMP& private(fkin, fconv, fvisc)
#endif
      do nThetaB=1,sizeThetaB
         nTheta=nThetaStart+nThetaB-1
         nThetaNHS=(nTheta+1)/2
         fkinAS(nThetaB) =0.0_cp
         fconvAS(nThetaB)=0.0_cp
         fviscAS(nThetaB)=0.0_cp
         fkin=0.0_cp
         fconv=0.0_cp
         fvisc=0.0_cp
         do nPhi=1,n_phi_max
            if ( l_anelastic_liquid .or. l_TP_form ) then
               fconv=vr(nPhi,nThetaB)*sr(nPhi,nThetaB)
            else
               fconv=temp0(nr)*vr(nPhi,nThetaB)*sr(nPhi,nThetaB)     +    &
               &          ViscHeatFac*ThExpNb*alpha0(nr)*temp0(nr)*       &
               &          orho1(nr)*vr(nPhi,nThetaB)*pr(nPhi,nThetaB)
            end if

            fkin=half*or2(nR)*orho2(nR)*(osn2(nThetaNHS)*(             &
            &                  vt(nPhi,nThetaB)*vt(nPhi,nThetaB)  +    &
            &                  vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )+    &
            &          or2(nR)*vr(nPhi,nThetaB)*vr(nPhi,nThetaB) )*    &
            &                             vr(nPhi,nThetaB)

            if ( nR/=n_r_icb .and. nR/=n_r_cmb ) then
               fvisc=-two*visc(nR)*orho1(nR)*vr(nPhi,nThetaB)*or2(nR)* (     &
               &                             dvrdr(nPhi,nThetaB)             &
               & -(two*or1(nR)+two*third*beta(nR))*vr(nPhi,nThetaB) )-       &
               &                       visc(nR)*orho1(nR)*vt(nPhi,nThetaB)*  &
               &                            osn2(nThetaNHS)* (               &
               &                       or2(nR)*dvrdt(nPhi,nThetaB)           &
               &                              +dvtdr(nPhi,nThetaB)           &
               &       -(two*or1(nR)+beta(nR))*vt(nPhi,nThetaB) )  -         &
               &       visc(nR)*orho1(nR)*vp(nPhi,nThetaB)*                  &
               &                               osn2(nThetaNHS)* (            &
               &                       or2(nR)*dvrdp(nPhi,nThetaB)           &
               &                              +dvpdr(nPhi,nThetaB)           &
               &       -(two*or1(nR)+beta(nR))*vp(nPhi,nThetaB) )
            end if

            fkinAS(nThetaB) = fkinAS(nThetaB)+fkin
            fconvAS(nThetaB)=fconvAS(nThetaB)+fconv
            fviscAS(nThetaB)=fviscAS(nThetaB)+fvisc
         end do
         fkinAS(nThetaB) =phiNorm* fkinAS(nThetaB)
         fconvAS(nThetaB)=phiNorm*fconvAS(nThetaB)
         fviscAS(nThetaB)=phiNorm*fviscAS(nThetaB)
      end do
#ifdef WITH_SHTNS
      !$OMP END PARALLEL DO
#endif

      if ( l_mag_nl) then
         nTheta=nThetaStart-1
#ifdef WITH_SHTNS
         !$OMP PARALLEL DO default(shared)         &
         !$OMP& private(nThetaB, nTheta, nPhi)     &
         !$OMP& private(fkin, fconv, fvisc)
#endif
         do nThetaB=1,sizeThetaB
            nTheta=nThetaStart+nThetaB-1
            nThetaNHS=(nTheta+1)/2
            fresAS(nThetaB) =0.0_cp
            fpoynAS(nThetaB)=0.0_cp
            fres=0.0_cp
            fpoyn=0.0_cp
            do nPhi=1,n_phi_max
                fres =osn2(nThetaNHS)*(                              &
                &              cbt(nPhi,nThetaB)*bp(nPhi,nThetaB)  - &
                &              cbp(nPhi,nThetaB)*bt(nPhi,nThetaB) )

                fpoyn=-orho1(nR)*or2(nR)*osn2(nThetaNHS)*(                        &
                &           vp(nPhi,nThetaB)*br(nPhi,nThetaB)*bp(nPhi,nThetaB)  - &
                &           vr(nPhi,nThetaB)*bp(nPhi,nThetaB)*bp(nPhi,nThetaB)  - &
                &           vr(nPhi,nThetaB)*bt(nPhi,nThetaB)*bt(nPhi,nThetaB)  + &
                &           vt(nPhi,nThetaB)*br(nPhi,nThetaB)*bt(nPhi,nThetaB) )

                fresAS(nThetaB) = fresAS(nThetaB)+fres
                fpoynAS(nThetaB)=fpoynAS(nThetaB)+fpoyn
            end do
            fresAS(nThetaB) =phiNorm* fresAS(nThetaB)
            fpoynAS(nThetaB)=phiNorm*fpoynAS(nThetaB)
         end do
#ifdef WITH_SHTNS
         !$OMP END PARALLEL DO
#endif

#ifdef WITH_SHTNS
      call spat_to_SH_axi(fresAS,fresLMr)
      call spat_to_SH_axi(fpoynAS,fpoynLMr)
#else
      call legTFAS2(fresLMr,fpoynLMr,fresAS,fpoynAS,l_max+1,nThetaStart,sizeThetaB)
#endif
      end if

      !-- Add contribution from thetas in block:
#ifdef WITH_SHTNS
      call spat_to_SH_axi(fviscAS,fviscLMr)
      call spat_to_SH_axi(fconvAS,fconvLMr)
      call spat_to_SH_axi(fkinAS,fkinLMr)
#else
      call legTFAS(fviscLMr,fviscAS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS2(fconvLMr,fkinLMr,fconvAS,fkinAS,l_max+1,nThetaStart,sizeThetaB)
#endif

   end subroutine get_fluxes
!------------------------------------------------------------------------------
   subroutine get_helicity(vr,vt,vp,cvr,dvrdt,dvrdp,dvtdr,dvpdr,HelLMr, &
              &            Hel2LMr,HelnaLMr,Helna2LMr,nR,nThetaStart)
      !
      !   Calculates axisymmetric contributions of helicity HelLMr and
      !   helicity**2  Hel2LMr in (l,m=0,r) space.
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      integer,  intent(in) :: nThetaStart
      real(cp), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
      real(cp), intent(in) :: cvr(nrp,nfs)
      real(cp), intent(in) :: dvrdt(nrp,nfs),dvrdp(nrp,nfs)
      real(cp), intent(in) :: dvtdr(nrp,nfs),dvpdr(nrp,nfs)

      !-- Output variables:
      real(cp), intent(out) :: HelLMr(l_max+1)
      real(cp), intent(out) :: Hel2LMr(l_max+1)
      real(cp), intent(out) :: HelnaLMr(l_max+1)
      real(cp), intent(out) :: Helna2LMr(l_max+1)

      !-- Local variables:
      integer :: nTheta,nThetaB
      integer :: nPhi,l
      real(cp) :: Helna,HelAS(nfs),Hel2AS(nfs)
      real(cp) :: Hel,HelnaAS(nfs),Helna2AS(nfs),phiNorm
      real(cp) :: vras,vtas,vpas,cvras,dvrdtas,dvrdpas,dvtdras,dvpdras
      real(cp) :: vrna,vtna,vpna,cvrna,dvrdtna,dvrdpna,dvtdrna,dvpdrna

      phiNorm=one/real(n_phi_max,cp)

      !-- Zero lm coeffs for first theta block:
      if ( nThetaStart == 1 ) then
         do l=1,l_max+1
            HelLMr(l) =0.0_cp
            Hel2LMr(l)=0.0_cp
            HelnaLMr(l) =0.0_cp
            Helna2LMr(l)=0.0_cp
         end do
      end if

      !--- Helicity:
      nTheta=nThetaStart-1
#ifdef WITH_SHTNS
      !$OMP PARALLEL DO default(shared) &
      !$OMP& private(nThetaB, nTheta, nPhi) &
      !$OMP& private(vras, cvras, vtas, vpas) &
      !$OMP& private(dvrdpas, dvpdras, dvtdras, dvrdtas) &
      !$OMP& private(vrna, cvrna, vtna) &
      !$OMP& private(dvrdpna, dvpdrna, dvtdrna, dvrdtna, Hel, Helna) &
      !$OMP& private(Hel2AS, HelAS, dvtdr, vpna, or2, HelnaAS, Helna2AS)
#endif
      do nThetaB=1,sizeThetaB
         nTheta=nThetaStart+nThetaB-1
         HelAS(nThetaB) =0.0_cp
         Hel2AS(nThetaB)=0.0_cp
         vras=0.0_cp
         cvras=0.0_cp
         vtas=0.0_cp
         vpas=0.0_cp
         dvrdpas=0.0_cp
         dvpdras=0.0_cp
         dvtdras=0.0_cp
         dvrdtas=0.0_cp
         do nPhi=1,n_phi_max
            vras=vras+vr(nPhi,nThetaB)
            cvras=cvras+cvr(nPhi,nThetaB)
            vtas=vtas+vt(nPhi,nThetaB)
            vpas=vpas+vp(nPhi,nThetaB)
            dvrdpas=dvrdpas+dvrdp(nPhi,nThetaB)
            dvpdras=dvpdras+dvpdr(nPhi,nThetaB)
            dvtdras=dvtdras+dvtdr(nPhi,nThetaB)
            dvrdtas=dvrdtas+dvrdt(nPhi,nThetaB)
         end do
         vras=vras*phiNorm
         cvras=cvras*phiNorm
         vtas=vtas*phiNorm
         vpas=vpas*phiNorm
         dvrdpas=dvrdpas*phiNorm
         dvpdras=dvpdras*phiNorm
         dvtdras=dvtdras*phiNorm
         dvrdtas=dvrdtas*phiNorm
         do nPhi=1,n_phi_max
            vrna   =   vr(nPhi,nThetaB)-vras
            cvrna  =  cvr(nPhi,nThetaB)-cvras
            vtna   =   vt(nPhi,nThetaB)-vtas
            vpna   =   vp(nPhi,nThetaB)-vpas
            dvrdpna=dvrdp(nPhi,nThetaB)-dvrdpas
            dvpdrna=dvpdr(nPhi,nThetaB)-beta(nR)*vp(nPhi,nThetaB) &
            &       -dvpdras+beta(nR)*vpas
            dvtdrna=dvtdr(nPhi,nThetaB)-beta(nR)*vt(nPhi,nThetaB) &
            &       -dvtdras+beta(nR)*vtas
            dvrdtna=dvrdt(nPhi,nThetaB)-dvrdtas
            Hel=or4(nR)*orho2(nR)*vr(nPhi,nThetaB)*cvr(nPhi,nThetaB) + &
            &              or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                                       vt(nPhi,nThetaB) * &
            &                          ( or2(nR)*dvrdp(nPhi,nThetaB) - &
            &                                    dvpdr(nPhi,nThetaB) + &
            &                         beta(nR)*   vp(nPhi,nThetaB) ) + &
            &                                       vp(nPhi,nThetaB) * &
            &                          (         dvtdr(nPhi,nThetaB) - &
            &                           beta(nR)*   vt(nPhi,nThetaB) - &
            &                            or2(nR)*dvrdt(nPhi,nThetaB) ) )
            Helna=                      or4(nR)*orho2(nR)*vrna*cvrna + &
            &              or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                       vtna*( or2(nR)*dvrdpna-dvpdrna ) + &
            &                       vpna*( dvtdrna-or2(nR)*dvrdtna ) )

            HelAS(nThetaB)   =HelAS(nThetaB) +Hel
            Hel2AS(nThetaB)  =Hel2AS(nThetaB)+Hel*Hel
            HelnaAS(nThetaB) =HelAS(nThetaB) +Helna
            Helna2AS(nThetaB)=Hel2AS(nThetaB)+Helna*Helna
         end do
         HelAS(nThetaB) =phiNorm*HelAS(nThetaB)
         Hel2AS(nThetaB)=phiNorm*Hel2AS(nThetaB)
         HelnaAS(nThetaB) =phiNorm*HelnaAS(nThetaB)
         Helna2AS(nThetaB)=phiNorm*Helna2AS(nThetaB)
      end do
#ifdef WITH_SHTNS
      !$OMP END PARALLEL DO
#endif

      !-- Add contribution from thetas in block:
#ifdef WITH_SHTNS
      call spat_to_SH_axi(HelAS, HelLMr)
      call spat_to_SH_axi(Hel2AS, Hel2LMr)
      call spat_to_SH_axi(HelnaAS, HelnaLMr)
      call spat_to_SH_axi(Helna2AS, Helna2LMr)
#else
      call legTFAS2(HelLMr,Hel2LMr,HelAS,Hel2AS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS2(HelnaLMr,Helna2LMr,HelnaAS,Helna2AS,l_max+1,nThetaStart,&
           &        sizeThetaB)
#endif

   end subroutine get_helicity
!------------------------------------------------------------------------------
   subroutine get_visc_heat(vr,vt,vp,cvr,dvrdr,dvrdt,dvrdp,dvtdr,&
              &             dvtdp,dvpdr,dvpdp,viscLMr,nR,nThetaStart)
      !
      !   Calculates axisymmetric contributions of the viscous heating
      !
      !

      !-- Input of variables
      integer,  intent(in) :: nR
      integer,  intent(in) :: nThetaStart
      real(cp), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs),cvr(nrp,nfs)
      real(cp), intent(in) :: dvrdr(nrp,nfs),dvrdt(nrp,nfs),dvrdp(nrp,nfs)
      real(cp), intent(in) :: dvtdr(nrp,nfs),dvtdp(nrp,nfs)
      real(cp), intent(in) :: dvpdr(nrp,nfs),dvpdp(nrp,nfs)

      !-- Output variables:
      real(cp), intent(out) :: viscLMr(l_max+1)

      !-- Local variables:
      integer :: nTheta,nThetaB,nPhi, nThetaNHS
      real(cp) :: viscAS(nfs),vischeat,csn2, phinorm

      phiNorm=two*pi/real(n_phi_max,cp)

      nTheta=nThetaStart-1
#ifdef WITH_SHTNS
      !$OMP PARALLEL DO default(shared)                     &
      !$OMP& private(nThetaB, nTheta, nPhi)                 &
      !$OMP& private(vischeat)
#endif
      do nThetaB=1,sizeThetaB
         nTheta=nThetaStart+nThetaB-1
         nThetaNHS=(nTheta+1)/2
         csn2     =cosn2(nThetaNHS)
         if ( mod(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

         viscAS(nThetaB)=0.0_cp
         do nPhi=1,n_phi_max
            vischeat=         or2(nR)*orho1(nR)*visc(nR)*(        &
            &     two*(                     dvrdr(nPhi,nThetaB) - & ! (1)
            &     (two*or1(nR)+beta(nR))*vr(nphi,nThetaB) )**2  + &
            &     two*( csn2*                  vt(nPhi,nThetaB) + &
            &                               dvpdp(nphi,nThetaB) + &
            &                               dvrdr(nPhi,nThetaB) - & ! (2)
            &     or1(nR)*               vr(nPhi,nThetaB) )**2  + &
            &     two*(                     dvpdp(nphi,nThetaB) + &
            &           csn2*                  vt(nPhi,nThetaB) + & ! (3)
            &     or1(nR)*               vr(nPhi,nThetaB) )**2  + &
            &          ( two*               dvtdp(nPhi,nThetaB) + &
            &                                 cvr(nPhi,nThetaB) - & ! (6)
            &      two*csn2*             vp(nPhi,nThetaB) )**2  + &
            &                                 osn2(nThetaNHS) * ( &
            &         ( r(nR)*              dvtdr(nPhi,nThetaB) - &
            &           (two+beta(nR)*r(nR))*  vt(nPhi,nThetaB) + & ! (4)
            &     or1(nR)*            dvrdt(nPhi,nThetaB) )**2  + &
            &         ( r(nR)*              dvpdr(nPhi,nThetaB) - &
            &           (two+beta(nR)*r(nR))*  vp(nPhi,nThetaB) + & ! (5)
            &     or1(nR)*            dvrdp(nPhi,nThetaB) )**2 )- &
            &    two*third*(  beta(nR)*        vr(nPhi,nThetaB) )**2 )

            viscAS(nThetaB)=viscAS(nThetaB)+vischeat
         end do
         viscAS(nThetaB)=phiNorm*viscAS(nThetaB)
      end do
#ifdef WITH_SHTNS
      !$OMP END PARALLEL DO
#endif

#ifdef WITH_SHTNS
      call spat_to_SH_axi(viscAS, viscLMr)
#else
      call legTFAS(viscLMr,viscAS,l_max+1,nThetaStart,sizeThetaB)
#endif

   end subroutine get_visc_heat
!------------------------------------------------------------------------------
end module nl_special_calc
