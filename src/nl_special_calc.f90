module nl_special_calc
 
   use precision_mod
   use truncation, only: nrp, n_phi_max, l_max, l_maxMag
   use const, only: pi, one, two, third, half
   use logic, only: l_mag_nl
   use physical_parameters, only: ek, ViscHeatFac
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: orho1, orho2, or2, or1, beta, temp0, &
                               visc, or4
   use blocking, only: sizeThetaB, nfs
   use horizontal_data, only: O_sin_theta_E2, cosTheta, sn2, osn2
   use legendre_grid_to_spec, only: legTFAS, legTFAS2
 
   implicit none

   private

   public :: get_nlBLayers, get_perpPar, get_fluxes, get_helicity

contains

   subroutine get_nlBLayers(vt,vp,dvtdr,dvpdr,dsdr,dsdt,dsdp,&
              &             uhLMr,duhLMr,gradsLMr,nR,nThetaStart)
      !
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
      do nThetaB=1,sizeThetaB
         nTheta=nTheta+1
         uhAS(nThetaB) =0.0_cp
         duhAS(nThetaB)=0.0_cp
         gradsAS(nThetaB)=0.0_cp
         do nPhi=1,n_phi_max
            uh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(     &
                          vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+  &
                          vp(nPhi,nThetaB)*vp(nPhi,nThetaB)  )
            duh=or2(nR)*orho2(nR)*O_sin_theta_E2(nTheta)*(            &
                                 dvtdr(nPhi,nThetaB)*vt(nPhi,nThetaB)-&
                 (or1(nR)+beta(nR))*vt(nPhi,nThetaB)*vt(nPhi,nThetaB)+&
                                 dvpdr(nPhi,nThetaB)*vp(nPhi,nThetaB)-&
                 (or1(nR)+beta(nR))*vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )
    
            grads =  dsdr(nPhi,nThetaB)*dsdr(nPhi,nThetaB)          &
                   +or2(nR)*O_sin_theta_E2(nTheta)*(                &
                           dsdt(nPhi,nThetaB)*dsdt(nPhi,nThetaB)    &
                          +dsdp(nPhi,nThetaB)*dsdp(nPhi,nThetaB) ) 
    
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
    
      !------ Add contribution from thetas in block:
      call legTFAS2(uhLMr,duhLMr,uhAS,duhAS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS(gradsLMr,gradsAS,l_max+1,nThetaStart,sizeThetaB)

   end subroutine get_nlBLayers
!------------------------------------------------------------------------------
   subroutine get_perpPar(vr,vt,vp,EperpLMr,EparLMr,  &
              &           EperpaxiLMr,EparaxiLMr,nR,nThetaStart)
      !
      !   Calculates the energies parallel and perpendicular to the rotation axis
      !
      !       Eperp = 0.5*(vs**2+vphi**2) with vs= vr*sin(theta)+vt*cos(theta)
      !       Epar  = 0.5*vz**2           with vz= vr*cos(theta)-vt*sin(theta)
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
      do nThetaB=1,sizeThetaB
         nTheta=nTheta+1
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
                    or2(nR)*sn2(nThetaNHS)*      vr(nPhi,nThetaB)*vr(nPhi,nThetaB) + &
                    (osn2(nThetaNHS)-one)*       vt(nPhi,nThetaB)*vt(nPhi,nThetaB) + &
                    two*or1(nR)*cosTheta(nTheta)*vr(nPhi,nThetaB)*vt(nPhi,nThetaB) + &
                    osn2(nThetaNHS)*                vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )
    
            Epar =half*or2(nR)*orho2(nR)*(                                           &
                    or2(nR)*(one-sn2(nThetaNHS))*vr(nPhi,nThetaB)*vr(nPhi,nThetaB) + &
                                                 vt(nPhi,nThetaB)*vt(nPhi,nThetaB) - &
                    two*or1(nR)*cosTheta(nTheta)*vr(nPhi,nThetaB)*vt(nPhi,nThetaB) )
    
            Eperpaxi=half*or2(nR)*orho2(nR)*(                  &
                      or2(nR)*sn2(nThetaNHS)*      vras*vras + &
                      (osn2(nThetaNHS)-one)*       vtas*vtas + &
                      two*or1(nR)*cosTheta(nTheta)*vras*vtas + &
                      osn2(nThetaNHS)*             vpas*vpas )
    
            Eparaxi =half*or2(nR)*orho2(nR)*(                  &
                      or2(nR)*(one-sn2(nThetaNHS))*vras*vras + &
                                                   vtas*vtas - &
                      two*or1(nR)*cosTheta(nTheta)*vras*vtas )
    
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
    
      !-- Add contribution from thetas in block:
      call legTFAS2(EperpLMr,EparLMr,EperpAS,EparAS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS2(EperpaxiLMr,EparaxiLMr,EperpaxiAS,EparaxiAS,l_max+1, &
                    nThetaStart,sizeThetaB)

   end subroutine get_perpPar
!------------------------------------------------------------------------------
   subroutine get_fluxes(vr,vt,vp,dvrdr,dvtdr,dvpdr,      &
                    &  dvrdt,dvrdp,sr,pr,br,bt,bp,cbt,cbp,&
                    &  fconvLMr,fkinLMr,fviscLMr,fpoynLMr,&
                    &  fresLMR,nR,nThetaStart)
      !
      !   Calculates the fluxes
      !       Fconv= rho * temp* (u_r s)
      !       Fkin = 1/2 * rho * u_r (u_r^2+u_theta^2+u_phi^2)
      !       Fvisc= -(u *\tensor(S) )_r
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
      do nThetaB=1,sizeThetaB
         nTheta=nTheta+1
         nThetaNHS=(nTheta+1)/2
         fkinAS(nThetaB) =0.0_cp
         fconvAS(nThetaB)=0.0_cp
         fviscAS(nThetaB)=0.0_cp
         fkin=0.0_cp
         fconv=0.0_cp
         fvisc=0.0_cp
         do nPhi=1,n_phi_max
            fconv=temp0(nr)*vr(nPhi,nThetaB)*sr(nPhi,nThetaB)     +    &
                       ViscHeatFac*orho1(nr)*vr(nPhi,nThetaB)*pr(nPhi,nThetaB)
    
            fkin=half*or2(nR)*orho2(nR)*(osn2(nThetaNHS)*(             &
                               vt(nPhi,nThetaB)*vt(nPhi,nThetaB)  +    &
                               vp(nPhi,nThetaB)*vp(nPhi,nThetaB) )+    &
                       or2(nR)*vr(nPhi,nThetaB)*vr(nPhi,nThetaB) )*    &
                                          vr(nPhi,nThetaB)
    
            if ( nR/=n_r_icb .and. nR/=n_r_cmb ) then
               fvisc=-two*visc(nR)*orho1(nR)*vr(nPhi,nThetaB)*or2(nR)* (     &
                                             dvrdr(nPhi,nThetaB)             &
                 -(two*or1(nR)+two*third*beta(nR))*vr(nPhi,nThetaB) )-       &
                                       visc(nR)*orho1(nR)*vt(nPhi,nThetaB)*  &
                                            osn2(nThetaNHS)* (               &
                                       or2(nR)*dvrdt(nPhi,nThetaB)           &
                                              +dvtdr(nPhi,nThetaB)           &
                       -(two*or1(nR)+beta(nR))*vt(nPhi,nThetaB) )  -         &
                       visc(nR)*orho1(nR)*vp(nPhi,nThetaB)*                  &
                                               osn2(nThetaNHS)* (            &
                                       or2(nR)*dvrdp(nPhi,nThetaB)           &
                                              +dvpdr(nPhi,nThetaB)           &
                       -(two*or1(nR)+beta(nR))*vp(nPhi,nThetaB) ) 
            end if
    
            fkinAS(nThetaB) = fkinAS(nThetaB)+fkin
            fconvAS(nThetaB)=fconvAS(nThetaB)+fconv
            fviscAS(nThetaB)=fviscAS(nThetaB)+fvisc
         end do
         fkinAS(nThetaB) =phiNorm* fkinAS(nThetaB)
         fconvAS(nThetaB)=phiNorm*fconvAS(nThetaB)
         fviscAS(nThetaB)=phiNorm*fviscAS(nThetaB)
      end do
    
      if ( l_mag_nl) then
         nTheta=nThetaStart-1
         do nThetaB=1,sizeThetaB
            nTheta=nTheta+1
            nThetaNHS=(nTheta+1)/2
            fresAS(nThetaB) =0.0_cp
            fpoynAS(nThetaB)=0.0_cp
            fres=0.0_cp
            fpoyn=0.0_cp
            do nPhi=1,n_phi_max
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
            end do
            fresAS(nThetaB) =phiNorm* fresAS(nThetaB)
            fpoynAS(nThetaB)=phiNorm*fpoynAS(nThetaB)
         end do
         call legTFAS2(fresLMr,fpoynLMr,fresAS,fpoynAS,l_max+1,nThetaStart,sizeThetaB)
      end if
    
      !-- Add contribution from thetas in block:
      call legTFAS(fviscLMr,fviscAS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS2(fconvLMr,fkinLMr,fconvAS,fkinAS,l_max+1,nThetaStart,sizeThetaB)

   end subroutine get_fluxes
!------------------------------------------------------------------------------
   subroutine get_helicity(vr,vt,vp,cvr,dvrdt,dvrdp,dvtdr,dvpdr,HelLMr, &
              &        Hel2LMr,HelnaLMr,Helna2LMr,nR,nThetaStart)
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
      do nThetaB=1,sizeThetaB
         nTheta=nTheta+1
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
         end do
         HelAS(nThetaB) =phiNorm*HelAS(nThetaB)
         Hel2AS(nThetaB)=phiNorm*Hel2AS(nThetaB)
         HelnaAS(nThetaB) =phiNorm*HelnaAS(nThetaB)
         Helna2AS(nThetaB)=phiNorm*Helna2AS(nThetaB)
      end do

      !-- Add contribution from thetas in block:
      call legTFAS2(HelLMr,Hel2LMr,HelAS,Hel2AS,l_max+1,nThetaStart,sizeThetaB)
      call legTFAS2(HelnaLMr,Helna2LMr,HelnaAS,Helna2AS,l_max+1,nThetaStart,sizeThetaB)

   end subroutine get_helicity
!------------------------------------------------------------------------------
end module nl_special_calc
