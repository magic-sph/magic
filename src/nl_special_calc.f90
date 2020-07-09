module nl_special_calc
   !
   ! This module allows to calculcate several diagnostics that need to be
   ! computed in the physical space (non-linear quantities)
   !

   use precision_mod
   use truncation, only: nrp, n_phi_max, n_r_icb, n_r_cmb, nThetaStart, nThetaStop
   use constants, only: pi, one, two, third, half
   use logic, only: l_mag_nl, l_anelastic_liquid
   use physical_parameters, only: ek, ViscHeatFac, ThExpNb
   use radial_functions, only: orho1, orho2, or2, or1, beta, temp0, &
       &                       visc, or4, r, alpha0
   use horizontal_data, only: O_sin_theta_E2, cosTheta, sn2, osn2, cosn2, gauss
   use parallel_mod

   implicit none

   private

   public :: get_nlBLayers, get_perpPar, get_fluxes, get_helicity, &
   &         get_visc_heat

contains

   subroutine get_nlBLayers(vt,vp,dvtdr,dvpdr,dsdr,dsdt,dsdp,uhAS,duhAS,gradsAS,nR)
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
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvtdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvpdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dsdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dsdt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dsdp(nrp,nThetaStart:nThetaStop)

      !-- Output variables:
      real(cp), intent(out) :: uhAS
      real(cp), intent(out) :: duhAS
      real(cp), intent(out) :: gradsAS

      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS
      real(cp) :: uh, duh, phiNorm, grads

      phiNorm=one/real(n_phi_max,cp)
      uhAS   =0.0_cp
      duhAS  =0.0_cp
      gradsAS=0.0_cp

      !--- Horizontal velocity uh and duh/dr + (grad T)**2
      !$omp parallel do default(shared)              &
      !$omp& private(nTheta, nPhi, uh, duh, grads)   &
      !$omp& reduction(+:uhAS,duhAS,gradsAS)
      do nTheta=nThetaStart,nThetaStop
         nThetaNHS=(nTheta+1)/2
         do nPhi=1,n_phi_max
            uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(   &
            &             vt(nPhi,nTheta)*vt(nPhi,nTheta)+  &
            &             vp(nPhi,nTheta)*vp(nPhi,nTheta)  )
            duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(            &
            &                    dvtdr(nPhi,nTheta)*vt(nPhi,nTheta)-  &
            &    (or1(nR)+beta(nR))*vt(nPhi,nTheta)*vt(nPhi,nTheta)+  &
            &                    dvpdr(nPhi,nTheta)*vp(nPhi,nTheta)-  &
            &    (or1(nR)+beta(nR))*vp(nPhi,nTheta)*vp(nPhi,nTheta) )

            grads =  dsdr(nPhi,nTheta)*dsdr(nPhi,nTheta)          &
            &      +or2(nR)*O_sin_theta_E2(nTheta)*(              &
            &              dsdt(nPhi,nTheta)*dsdt(nPhi,nTheta)    &
            &             +dsdp(nPhi,nTheta)*dsdp(nPhi,nTheta) )

            uhAS=uhAS+phiNorm*gauss(nThetaNHS)*sqrt(uh)
            if (uh /= 0.0_cp) duhAS=duhAS+phiNorm*gauss(nThetaNHS)*abs(duh)/sqrt(uh)
            gradsAS=gradsAS+phiNorm*gauss(nThetaNHS)*grads
         end do
      end do
      !$omp end parallel do

      if ( n_ranks_theta>1 ) then
         call MPI_AllReduce(MPI_IN_PLACE, uhAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, duhAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, gradsAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
      end if

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
      real(cp), intent(in) :: vr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)

      !-- Output variables:
      real(cp), intent(out) :: EperpAS,EparAS
      real(cp), intent(out) :: EperpaxiAS,EparaxiAS

      !-- Local variables:
      integer :: nTheta,nThetaNHS,nPhi
      real(cp) :: vras,vtas,vpas,phiNorm,Eperp,Epar,Eperpaxi,Eparaxi

      phiNorm=one/real(n_phi_max,cp)
      EperpAS   =0.0_cp
      EparAS    =0.0_cp
      EperpaxiAS=0.0_cp
      EparaxiAS =0.0_cp

      !$omp parallel do default(shared)                             &
      !$omp& private(nTheta, nPhi, Eperp, Epar, Eperpaxi, Eparaxi)  &
      !$omp& private(vras, vtas, vpas, nThetaNHS)                   &
      !$omp& reduction(+:EperpAS,EparAS,EperpaxiAS,EparaxiAS)
      do nTheta=nThetaStart,nThetaStop
         nThetaNHS=(nTheta+1)/2

         Eperp   =0.0_cp
         Epar    =0.0_cp
         Eperpaxi=0.0_cp
         Eparaxi =0.0_cp
         vras    =0.0_cp
         vtas    =0.0_cp
         vpas    =0.0_cp

         do nPhi=1,n_phi_max
            vras=vras+vr(nPhi,nTheta)
            vtas=vtas+vt(nPhi,nTheta)
            vpas=vpas+vp(nPhi,nTheta)
         end do
         vras=vras*phiNorm
         vtas=vtas*phiNorm
         vpas=vpas*phiNorm

         do nPhi=1,n_phi_max
            Eperp=half*or2(nR)*orho2(nR)*(                                         &
            &       or2(nR)*sn2(nThetaNHS)*      vr(nPhi,nTheta)*vr(nPhi,nTheta) + &
            &       (osn2(nThetaNHS)-one)*       vt(nPhi,nTheta)*vt(nPhi,nTheta) + &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nPhi,nTheta)*vt(nPhi,nTheta) + &
            &       osn2(nThetaNHS)*             vp(nPhi,nTheta)*vp(nPhi,nTheta) )

            Epar =half*or2(nR)*orho2(nR)*(                                         &
            &       or2(nR)*(one-sn2(nThetaNHS))*vr(nPhi,nTheta)*vr(nPhi,nTheta) + &
            &                                    vt(nPhi,nTheta)*vt(nPhi,nTheta) - &
            &       two*or1(nR)*cosTheta(nTheta)*vr(nPhi,nTheta)*vt(nPhi,nTheta) )

            Eperpaxi=half*or2(nR)*orho2(nR)*(                  &
            &         or2(nR)*sn2(nThetaNHS)*      vras*vras + &
            &         (osn2(nThetaNHS)-one)*       vtas*vtas + &
            &         two*or1(nR)*cosTheta(nTheta)*vras*vtas + &
            &         osn2(nThetaNHS)*             vpas*vpas )

            Eparaxi =half*or2(nR)*orho2(nR)*(                  &
            &         or2(nR)*(one-sn2(nThetaNHS))*vras*vras + &
            &                                      vtas*vtas - &
            &         two*or1(nR)*cosTheta(nTheta)*vras*vtas )

            EperpAS   =   EperpAS+phiNorm*gauss(nThetaNHS)*Eperp
            EparAS    =    EparAS+phiNorm*gauss(nThetaNHS)*Epar
            EperpaxiAS=EperpaxiAS+phiNorm*gauss(nThetaNHS)*Eperpaxi
            EparaxiAS = EparaxiAS+phiNorm*gauss(nThetaNHS)*Eparaxi
         end do
      end do
      !$omp end parallel do

      if ( n_ranks_theta>1 ) then
         call MPI_AllReduce(MPI_IN_PLACE, EperpAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, EparAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, EperpaxiAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, EparaxiAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
      end if

   end subroutine get_perpPar
!------------------------------------------------------------------------------
   subroutine get_fluxes(vr,vt,vp,dvrdr,dvtdr,dvpdr,dvrdt,dvrdp,sr,pr,br,bt,bp, &
              &          cbt,cbp,fconvAS,fkinAS,fviscAS,fpoynAS,fresAS,nR)
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
      real(cp), intent(in) :: vr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvtdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvpdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: sr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: pr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: br(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cbt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cbp(nrp,nThetaStart:nThetaStop)

      !-- Output variables:
      real(cp), intent(out) :: fkinAS, fconvAS
      real(cp), intent(out) :: fviscAS, fresAS, fpoynAS

      !-- Local variables:
      integer :: nTheta,nThetaNHS,nPhi
      real(cp) :: fkin, fconv, phiNorm, fvisc, fpoyn, fres

      phiNorm=two*pi/real(n_phi_max,cp)
      fkinAS =0.0_cp
      fconvAS=0.0_cp
      fviscAS=0.0_cp
      fpoynAS=0.0_cp
      fresAS =0.0_cp

      !$omp parallel do default(shared)                           &
      !$omp& private(nTheta, nPhi, fkin, fconv, fvisc, nThetaNHS) &
      !$omp& reduction(+:fkinAS,fconvAS,fviscAS)
      do nTheta=nThetaStart,nThetaStop
         nThetaNHS=(nTheta+1)/2
         fkin=0.0_cp
         fconv=0.0_cp
         fvisc=0.0_cp
         do nPhi=1,n_phi_max
            if ( l_anelastic_liquid ) then
               fconv=vr(nPhi,nTheta)*sr(nPhi,nTheta)
            else
               fconv=temp0(nr)*vr(nPhi,nTheta)*sr(nPhi,nTheta)     +    &
               &          ViscHeatFac*ThExpNb*alpha0(nr)*temp0(nr)*     &
               &          orho1(nr)*vr(nPhi,nTheta)*pr(nPhi,nTheta)
            end if

            fkin=half*or2(nR)*orho2(nR)*(osn2(nThetaNHS)*(           &
            &                  vt(nPhi,nTheta)*vt(nPhi,nTheta)  +    &
            &                  vp(nPhi,nTheta)*vp(nPhi,nTheta) )+    &
            &          or2(nR)*vr(nPhi,nTheta)*vr(nPhi,nTheta) )*    &
            &                             vr(nPhi,nTheta)

            if ( nR /= n_r_icb .and. nR /= n_r_cmb ) then
               fvisc=-two*visc(nR)*orho1(nR)*vr(nPhi,nTheta)*or2(nR)* (     &
               &                             dvrdr(nPhi,nTheta)             &
               & -(two*or1(nR)+two*third*beta(nR))*vr(nPhi,nTheta) )-       &
               &                       visc(nR)*orho1(nR)*vt(nPhi,nTheta)*  &
               &                           osn2(nThetaNHS)* (               &
               &                       or2(nR)*dvrdt(nPhi,nTheta)           &
               &                              +dvtdr(nPhi,nTheta)           &
               &       -(two*or1(nR)+beta(nR))*vt(nPhi,nTheta) )  -         &
               &       visc(nR)*orho1(nR)*vp(nPhi,nTheta)*                  &
               &                              osn2(nThetaNHS)* (            &
               &                       or2(nR)*dvrdp(nPhi,nTheta)           &
               &                              +dvpdr(nPhi,nTheta)           &
               &       -(two*or1(nR)+beta(nR))*vp(nPhi,nTheta) )
            end if

            fkinAS = fkinAS+phiNorm*gauss(nThetaNHS)*fkin
            fconvAS=fconvAS+phiNorm*gauss(nThetaNHS)*fconv
            fviscAS=fviscAS+phiNorm*gauss(nThetaNHS)*fvisc
         end do
      end do
      !$omp end parallel do

      if ( n_ranks_theta>1 ) then
         call MPI_AllReduce(MPI_IN_PLACE, fkinAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, fconvAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, fviscAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
      end if

      if ( l_mag_nl) then
         !$omp parallel do default(shared)                           &
         !$omp& private(nTheta, nPhi, fkin, fconv, fvisc, nThetaNHS) &
         !$omp& reduction(+:fpoynAS,fresAS)
         do nTheta=nThetaStart,nThetaStop
            nThetaNHS=(nTheta+1)/2
            fres=0.0_cp
            fpoyn=0.0_cp
            do nPhi=1,n_phi_max
                fres =osn2(nThetaNHS)*(                            &
                &              cbt(nPhi,nTheta)*bp(nPhi,nTheta)  - &
                &              cbp(nPhi,nTheta)*bt(nPhi,nTheta) )

                fpoyn=-orho1(nR)*or2(nR)*osn2(nThetaNHS)*(                     &
                &           vp(nPhi,nTheta)*br(nPhi,nTheta)*bp(nPhi,nTheta)  - &
                &           vr(nPhi,nTheta)*bp(nPhi,nTheta)*bp(nPhi,nTheta)  - &
                &           vr(nPhi,nTheta)*bt(nPhi,nTheta)*bt(nPhi,nTheta)  + &
                &           vt(nPhi,nTheta)*br(nPhi,nTheta)*bt(nPhi,nTheta) )

                fresAS = fresAS+phiNorm*gauss(nThetaNHS)*fres
                fpoynAS=fpoynAS+phiNorm*gauss(nThetaNHS)*fpoyn
            end do
         end do
         !$omp end parallel do
         if ( n_ranks_theta>1 ) then
            call MPI_AllReduce(MPI_IN_PLACE, fresAS, 1, MPI_DEF_REAL, &
                 &             MPI_SUM, comm_theta, ierr)
            call MPI_AllReduce(MPI_IN_PLACE, fpoynAS, 1, MPI_DEF_REAL, &
                 &             MPI_SUM, comm_theta, ierr)
         end if
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
      real(cp), intent(in) :: vr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cvr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvtdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvpdr(nrp,nThetaStart:nThetaStop)

      !-- Output variables:
      real(cp), intent(out) :: HelAS(2)
      real(cp), intent(out) :: Hel2AS(2)
      real(cp), intent(out) :: HelnaAS(2)
      real(cp), intent(out) :: Helna2AS(2)
      real(cp), intent(out) :: HelEAAS

      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS
      real(cp) :: Helna, Hel, phiNorm
      real(cp) :: vras,vtas,vpas,cvras,dvrdtas,dvrdpas,dvtdras,dvpdras
      real(cp) :: vrna,vtna,vpna,cvrna,dvrdtna,dvrdpna,dvtdrna,dvpdrna

      !-- Remark: 2pi not used the normalization below
      !-- this is why we have a 2pi factor after radial integration
      !-- in the subroutine outHelicity()
      phiNorm=one/real(n_phi_max,cp)
      HelAS(:)   =0.0_cp
      Hel2AS(:)  =0.0_cp
      HelnaAS(:) =0.0_cp
      Helna2AS(:)=0.0_cp
      HelEAAS    =0.0_cp

      !--- Helicity:
      !$omp parallel do default(shared)                     &
      !$omp& private(nTheta, nPhi, Hel, Helna)              &
      !$omp& private(vrna, cvrna, vtna, vpna)               &
      !$omp& private(dvrdpna, dvpdrna, dvtdrna, dvrdtna)    &
      !$omp& reduction(+:HelAS,Hel2AS,HelnaAS,Helna2AS,HelEAAS)
      do nTheta=nThetaStart,nThetaStop
         nThetaNHS=(nTheta+1)/2
         vras=0.0_cp
         cvras=0.0_cp
         vtas=0.0_cp
         vpas=0.0_cp
         dvrdpas=0.0_cp
         dvpdras=0.0_cp
         dvtdras=0.0_cp
         dvrdtas=0.0_cp
         do nPhi=1,n_phi_max
            vras=vras+vr(nPhi,nTheta)
            cvras=cvras+cvr(nPhi,nTheta)
            vtas=vtas+vt(nPhi,nTheta)
            vpas=vpas+vp(nPhi,nTheta)
            dvrdpas=dvrdpas+dvrdp(nPhi,nTheta)
            dvpdras=dvpdras+dvpdr(nPhi,nTheta)
            dvtdras=dvtdras+dvtdr(nPhi,nTheta)
            dvrdtas=dvrdtas+dvrdt(nPhi,nTheta)
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
            vrna   =   vr(nPhi,nTheta)-vras
            cvrna  =  cvr(nPhi,nTheta)-cvras
            vtna   =   vt(nPhi,nTheta)-vtas
            vpna   =   vp(nPhi,nTheta)-vpas
            dvrdpna=dvrdp(nPhi,nTheta)-dvrdpas
            dvpdrna=dvpdr(nPhi,nTheta)-beta(nR)*vp(nPhi,nTheta) &
            &       -dvpdras+beta(nR)*vpas
            dvtdrna=dvtdr(nPhi,nTheta)-beta(nR)*vt(nPhi,nTheta) &
            &       -dvtdras+beta(nR)*vtas
            dvrdtna=dvrdt(nPhi,nTheta)-dvrdtas
            Hel=or4(nR)*orho2(nR)*vr(nPhi,nTheta)*cvr(nPhi,nTheta) +  &
            &             or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                                       vt(nPhi,nTheta) * &
            &                          ( or2(nR)*dvrdp(nPhi,nTheta) - &
            &                                    dvpdr(nPhi,nTheta) + &
            &                         beta(nR)*   vp(nPhi,nTheta) ) + &
            &                                       vp(nPhi,nTheta) * &
            &                          (         dvtdr(nPhi,nTheta) - &
            &                           beta(nR)*   vt(nPhi,nTheta) - &
            &                            or2(nR)*dvrdt(nPhi,nTheta) ) )
            Helna=                     or4(nR)*orho2(nR)*vrna*cvrna + &
            &             or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)* ( &
            &                      vtna*( or2(nR)*dvrdpna-dvpdrna ) + &
            &                      vpna*( dvtdrna-or2(nR)*dvrdtna ) )

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

      if ( n_ranks_theta>1 ) then
         call MPI_AllReduce(MPI_IN_PLACE, HelAS, 2, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, Hel2AS, 2, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, HelnaAS, 2, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, Helna2AS, 2, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
         call MPI_AllReduce(MPI_IN_PLACE, HelEAAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
      end if

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
      real(cp), intent(in) :: vr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cvr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvrdp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvtdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvtdp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvpdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvpdp(nrp,nThetaStart:nThetaStop)

      !-- Output variables:
      real(cp), intent(out) :: viscAS

      !-- Local variables:
      integer :: nTheta, nPhi, nThetaNHS
      real(cp) :: vischeat,csn2, phinorm

      phiNorm=two*pi/real(n_phi_max,cp)
      viscAS=0.0_cp

      !$omp parallel do default(shared)                       &
      !$omp& private(nTheta, nPhi, vischeat, csn2, nThetaNHS) &
      !$omp& reduction(+:viscAS)
      do nTheta=nThetaStart, nThetaStop
         nThetaNHS=(nTheta+1)/2
         csn2     =cosn2(nThetaNHS)
         if ( mod(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

         do nPhi=1,n_phi_max
            vischeat=         or2(nR)*orho1(nR)*visc(nR)*(       &
            &     two*(                     dvrdr(nPhi,nTheta) - & ! (1)
            &     (two*or1(nR)+beta(nR))*vr(nphi,nTheta) )**2  + &
            &     two*( csn2*                  vt(nPhi,nTheta) + &
            &                               dvpdp(nphi,nTheta) + &
            &                               dvrdr(nPhi,nTheta) - & ! (2)
            &     or1(nR)*               vr(nPhi,nTheta) )**2  + &
            &     two*(                     dvpdp(nphi,nTheta) + &
            &           csn2*                  vt(nPhi,nTheta) + & ! (3)
            &     or1(nR)*               vr(nPhi,nTheta) )**2  + &
            &          ( two*               dvtdp(nPhi,nTheta) + &
            &                                 cvr(nPhi,nTheta) - & ! (6)
            &      two*csn2*             vp(nPhi,nTheta) )**2  + &
            &                                osn2(nThetaNHS) * ( &
            &         ( r(nR)*              dvtdr(nPhi,nTheta) - &
            &           (two+beta(nR)*r(nR))*  vt(nPhi,nTheta) + & ! (4)
            &     or1(nR)*            dvrdt(nPhi,nTheta) )**2  + &
            &         ( r(nR)*              dvpdr(nPhi,nTheta) - &
            &           (two+beta(nR)*r(nR))*  vp(nPhi,nTheta) + & ! (5)
            &     or1(nR)*            dvrdp(nPhi,nTheta) )**2 )- &
            &    two*third*(  beta(nR)*        vr(nPhi,nTheta) )**2 )

            viscAS=viscAS+phiNorm*gauss(nThetaNHS)*vischeat
         end do
      end do
      !$omp end parallel do

      if ( n_ranks_theta>1 ) then
         call MPI_AllReduce(MPI_IN_PLACE, viscAS, 1, MPI_DEF_REAL, &
              &             MPI_SUM, comm_theta, ierr)
      end if

   end subroutine get_visc_heat
!------------------------------------------------------------------------------
end module nl_special_calc
