module power
   !
   ! This module handles the writing of the power budget
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use communications, only: gather_from_Rloc, reduce_radial, send_lm_pair_to_master
   use truncation, only: n_r_ic_maxMag, n_r_max, n_r_ic_max, &
       &                 n_r_maxMag, n_phi_max, n_theta_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: r_cmb, r_icb, r, rscheme_oc, chebt_ic, &
       &                       or2, O_r_ic2, lambda, temp0, or1,      &
       &                       O_r_ic, rgrav, r_ic, dr_fac_ic,        &
       &                       alpha0, orho1, otemp1, beta, visc, kappa,&
       &                       dentropy0,dLkappa, dLtemp0
   use radial_der, only:get_ddr
   use physical_parameters, only: kbotv, ktopv, opm, LFfac, BuoFac, &
        &                         ChemFac, ThExpNb, ViscHeatFac,opr,raScaled, &
        &                         w_orbit_th, amp_tidal, ek
   use num_param, only: tScale, eScale
   use blocking, only: lo_map, st_map, llm, ulm, llmMag, ulmMag
   use horizontal_data, only: dLh, gauss, O_sin_theta_E2, cosn_theta_E2, hdif_S
   use logic, only: l_rot_ic, l_SRIC, l_rot_ma, l_SRMA, l_save_out, &
       &            l_conv, l_cond_ic, l_heat, l_mag,               &
       &            l_chemical_conv, l_anelastic_liquid, l_tidal
   use fields, only: wtidal_LMloc
   use output_data, only: tag
   use mean_sd, only: mean_sd_type
   use useful, only: cc2real, cc22real, round_off
   use integration, only: rInt_R, rIntIC
   use outRot, only: get_viscous_torque
   use constants, only: one, two, half, pi, third, ci, four, three, zero
   use init_fields, only: X, dX

   implicit none

   private

   type(mean_sd_type) :: buo_ave, buo_chem_ave, visc_ave, ohm_ave
   type(mean_sd_type) :: buo_tidal_ave, tidal_ave, visc_tidal_ave, heat_tidal_ave
   real(cp) :: powerDiff, eDiffInt
   real(cp), allocatable :: viscASr(:)
   integer :: n_power_file, n_calls
   character(len=72) :: power_file

   public :: initialize_output_power, get_power, finalize_output_power, &
   &         get_visc_heat

contains

   subroutine initialize_output_power
      !
      ! Memory allocation
      !

      call buo_ave%initialize(1,n_r_max)
      call buo_chem_ave%initialize(1,n_r_max)
      call visc_ave%initialize(1,n_r_max)
      call ohm_ave%initialize(1,n_r_max)
      call buo_tidal_ave%initialize(1,n_r_max)
      call visc_tidal_ave%initialize(1,n_r_max)
      call heat_tidal_ave%initialize(1,n_r_max)
      call tidal_ave%initialize(1,n_r_max)

      allocate( viscASr(nRstart:nRstop) )
      viscASr(:)=0.0_cp
      bytes_allocated=bytes_allocated+(nRstop-nRstart+1)*SIZEOF_DEF_REAL

      n_calls = 0
      powerDiff=0.0_cp
      eDiffInt =0.0_cp

      power_file='power.'//tag
      if ( rank == 0 .and. (.not. l_save_out) ) then
         open(newunit=n_power_file, file=power_file, status='new')
      end if

   end subroutine initialize_output_power
!----------------------------------------------------------------------------
   subroutine finalize_output_power

     call buo_ave%finalize()
      call buo_chem_ave%finalize()
      call visc_ave%finalize()
      call ohm_ave%finalize()
      call buo_tidal_ave%finalize()
      call visc_tidal_ave%finalize()
      call heat_tidal_ave%finalize()
      call tidal_ave%finalize()
     
      if ( rank == 0 .and. (.not. l_save_out) ) then
         close(n_power_file)
      end if

   end subroutine finalize_output_power
!----------------------------------------------------------------------------
   subroutine get_power(time,timePassed,timeNorm,l_stop_time, &
              &         omega_IC,omega_MA,lorentz_torque_IC,  &
              &         lorentz_torque_MA,w,z,dz,s,           &
              &         xi,b,ddb,aj,dj,db_ic,ddb_ic,aj_ic,    &
              &         dj_ic,viscDiss,ohmDiss,dw,ddw)
      !
      !  This subroutine calculates power and dissipation of
      !  the core/mantle system.
      !  Energy input into the outer core is by buoyancy and
      !  possibly viscous accelarations at the boundaries if
      !  the rotation rates of inner core or mantle are prescribed
      !  and kept fixed.
      !  The losses are due to Ohmic and viscous dissipation.
      !  If inner core and mantel are allowed to change their
      !  rotation rates due to viscous forces this power is not
      !  lost from the system and has to be respected.
      !
      !  The output is written into a file  power.TAG.
      !

      !-- Input of variables:
      logical,     intent(in) :: l_stop_time
      real(cp),    intent(in) :: time,timePassed,timeNorm
      real(cp),    intent(in) :: omega_IC,omega_MA
      real(cp),    intent(in) :: lorentz_torque_IC,lorentz_torque_MA
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Output:
      real(cp),    intent(out) :: viscDiss,ohmDiss

      !-- Local:
      complex(cp) :: ds(llm:ulm,n_r_max)
      complex(cp) :: dds(llm:ulm,n_r_max)
      integer :: n_r,lm,l,m,fileHandle
      real(cp) :: r_ratio
      real(cp) :: viscHeatR_global(n_r_max)
      real(cp) :: curlB2,buoy,curlB2_IC,buoy_chem,viscHeat
      real(cp) :: buoy_tidal, tidalpow, visc_tidal, heat_tidal
      complex(cp):: f_R, f_S, f_T
      real(cp) :: curlB2_r(n_r_max),curlB2_r_global(n_r_max)
      real(cp) :: buoy_r(n_r_max),buoy_r_global(n_r_max)
      real(cp) :: buoy_tidal_r(n_r_max), buoy_tidal_r_global(n_r_max)
      real(cp) :: tidalpow_r(n_r_max), tidalpow_r_global(n_r_max)
      real(cp) :: visc_tidal_r(n_r_max), visc_tidal_r_global(n_r_max)
      real(cp) :: heat_tidal_r(n_r_max), heat_tidal_r_global(n_r_max)
      real(cp) :: buoy_chem_r(n_r_max),buoy_chem_r_global(n_r_max)
      real(cp) :: curlB2_rIC(n_r_ic_max),curlB2_rIC_global(n_r_ic_max)
      real(cp) :: viscous_torque_ic,viscous_torque_ma
      real(cp) :: z10ICB,z10CMB,drz10ICB,drz10CMB
      real(cp) :: powerIC,powerMA
      real(cp) :: powerDiffOld,powerDiffT
      complex(cp) :: laplace,Bh, eiwt
      character(len=76) :: fileName

      if (l_tidal) then
         eiwt=cmplx(cos(-w_orbit_th*time),sin(-w_orbit_th*time))
         call get_ddr(s,ds,dds,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc)
      end if

      do n_r=1,n_r_max

         !if ( l_conv ) then
         !   curlU2_r(n_r)=0.0_cp
         !   !do lm=2,lm_max
         !   do lm=llm,ulm
         !      l=lo_map%lm2l(lm)
         !      m=lo_map%lm2m(lm)
         !      laplace=dLh(st_map%lm2(l,m))*or2(n_r)*w(lm,n_r)-ddw(lm,n_r)
         !      curlU2_r(n_r)=     curlU2_r(n_r) + dLh(st_map%lm2(l,m)) * ( &
         !           dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(z(lm,n_r),m)  + &
         !           cc2real(dz(lm,n_r),m) + &
         !           cc2real(laplace,m)    )
         !   end do
         !end if

         if ( l_mag ) then
            curlB2_r(n_r)=0.0_cp
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               laplace=dLh(st_map%lm2(l,m))*or2(n_r)*b(lm,n_r)-ddb(lm,n_r)
               curlB2_r(n_r)=curlB2_r(n_r) +  LFfac*opm*eScale*              &
               &             dLh(st_map%lm2(l,m))*lambda(n_r)*(              &
               &             dLh(st_map%lm2(l,m))*or2(n_r)*                  &
               &             cc2real(aj(lm,n_r),m) + cc2real(dj(lm,n_r),m) + &
               &             cc2real(laplace,m)    )
            end do
         end if

         if ( l_heat ) then
            buoy_r(n_r)=0.0_cp
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               buoy_r(n_r)=buoy_r(n_r) + eScale*                 &
               &           dLh(st_map%lm2(l,m))*BuoFac*          &
               &           rgrav(n_r)*cc22real(w(lm,n_r),s(lm,n_r),m)
            end do
         end if
         if ( l_tidal ) then
            buoy_tidal_r(n_r)=0.0_cp
            visc_tidal_r(n_r)=0.0_cp
            tidalpow_r(n_r)=0.0_cp
            heat_tidal_r(n_r)=0.0_cp
            do lm=llm,ulm !max(2,llm) ?
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if (l_heat) then
                  buoy_tidal_r(n_r)=buoy_tidal_r(n_r) - eScale*     &
                       &           dLh(st_map%lm2(l,m))*BuoFac*          &
                       !&           r(n_r)**2*BuoFac* &
                       &           rgrav(n_r)*cc22real(eiwt*wtidal_LMloc(lm,n_r),s(lm,n_r),m)
                  heat_tidal_r(n_r)=heat_tidal_r(n_r) &
                       & + r(n_r)**2*raScaled*opr*rgrav(n_r)*orho1(n_r)*otemp1(n_r)/dentropy0(n_r)*&
                       &     opr*hdif_S(l)*kappa(n_R)*(cc22real(dds(lm,n_r),s(lm,n_r),m)+&
                       &    (beta(n_r)+dLtemp0(n_r)+two*or1(n_r) +dLkappa(n_r)) &
                       & *cc22real(ds(lm,n_r),s(lm,n_r),m)-l*(l+one)*or2(n_r)*cc2real(s(lm,n_r),m)) !anelastic to test !rho factor
               end if
               visc_tidal_r(n_r)=visc_tidal_r(n_r)+(one/6.0D0)*or2(n_r)**2*l*(l+one)*( &
                    & l*(l+one)*cc2real(6.0D0*w(lm,n_r)+two*r(n_r)*w(lm,n_r)*beta(n_r)- &
                    & three*r(n_r)*dw(lm,n_r),m)+three*(cc2real(r(n_r)**2*beta(n_r)*dw(lm,n_r)&
                    & -l*(l+one)*w(lm,n_r)-r(n_r)*(-two*dw(lm,n_r)+r(n_r)*ddw(lm,n_r)),m)&
                    & +r(n_r)**2*(cc2real(two*z(lm,n_r)+r(n_r)*z(lm,n_r)*beta(n_r)-r(n_r)* &
                    & dz(lm,n_r),m)+(-two+l+l**2)*(cc2real(z(lm,n_r),m)+&
                    & cc2real(dw(lm,n_r),m)))))
               if (l==2 .and. m ==2) then
                  f_R=-four*eiwt*amp_tidal*w_orbit_th*or1(n_r)*6.0D0*X(n_r)
                  f_S=-four*eiwt*amp_tidal*w_orbit_th*(X(n_r)+r(n_r)*dX(n_r))
                  tidalpow_r(n_r)=tidalpow_r(n_r)+(cc22real(f_R,-w(lm,n_r),m)+&
                       & cc22real(f_S,-dw(lm,n_r),m))
               end if
               if(m==2.and.l==3) then
                  f_T=8.0D0*ci*eiwt*amp_tidal*w_orbit_th &
                       & *(2.0D0*X(n_r)-r(n_r)*dX(n_r))/sqrt(7.0D0)
                  tidalpow_r(n_r)=tidalpow_r(n_r)+cc22real(f_T,-z(lm,n_r),m)
               end if
            end do
            tidalpow_r(n_r)=half*tidalpow_r(n_r)
            visc_tidal_r(n_r)=visc_tidal_r(n_r)*ek*two*orho1(n_r)*visc(n_r)
         end if
         if ( l_chemical_conv ) then
            buoy_chem_r(n_r)=0.0_cp
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               buoy_chem_r(n_r)=buoy_chem_r(n_r) + eScale*              &
               &                dLh(st_map%lm2(l,m))*ChemFac*           &
               &                rgrav(n_r)*cc22real(w(lm,n_r),xi(lm,n_r),m)
            end do
         end if

      end do    ! radial grid points

      if ( l_mag ) call reduce_radial(curlB2_r, curlB2_r_global, 0)
      if ( l_heat ) call reduce_radial(buoy_r, buoy_r_global, 0)
      if ( l_chemical_conv ) call reduce_radial(buoy_chem_r, buoy_chem_r_global, 0)
      if ( l_conv ) call gather_from_Rloc(viscASr, viscHeatR_global, 0)
      if ( l_tidal ) then
         call reduce_radial(visc_tidal_r, visc_tidal_r_global, 0)
         call reduce_radial(tidalpow_r, tidalpow_r_global, 0)
         if (l_heat) then
            call reduce_radial(buoy_tidal_r, buoy_tidal_r_global, 0)
            call reduce_radial(heat_tidal_r, heat_tidal_r_global, 0)
         end if
      end if

      curlB2   =0.0_cp
      viscHeat =0.0_cp
      buoy     =0.0_cp
      buoy_chem=0.0_cp
      buoy_tidal=0.0_cp
      visc_tidal=0.0_cp
      tidalpow=0.0_cp
      heat_tidal=0.0_cp
      
      if ( rank == 0 ) then
         n_calls = n_calls+1
         !-- Transform to cheb space:
         if ( l_conv ) then
            !curlU2MeanR=curlU2MeanR+timePassed*curlU2_r_global*eScale
            !curlU2=rInt_R(curlU2_r_global,r,rscheme_oc)
            call visc_ave%compute(viscHeatR_global, n_calls, timePassed, timeNorm)
            viscHeat=eScale*rInt_R(viscHeatR_global,r,rscheme_oc)
         end if
         if ( l_mag )  then
            call ohm_ave%compute(curlB2_r_global, n_calls, timePassed, timeNorm)
            curlB2=rInt_R(curlB2_r_global,r,rscheme_oc)
         end if
         if ( l_heat ) then
            call buo_ave%compute(buoy_r_global, n_calls, timePassed, timeNorm)
            buoy=rInt_R(buoy_r_global,r,rscheme_oc)
         end if
         if ( l_tidal ) then
            tidalpow=0.0_cp
            call tidal_ave%compute(tidalpow_r_global, n_calls, timePassed, timeNorm)
            tidalpow=eScale*rInt_R(tidalpow_r_global,r,rscheme_oc)
            call visc_tidal_ave%compute(buoy_tidal_r_global, n_calls, timePassed, timeNorm)
            visc_tidal=-eScale*rInt_R(visc_tidal_r_global,r,rscheme_oc)
            if (l_heat ) then
               call buo_tidal_ave%compute(buoy_tidal_r_global, n_calls, timePassed, timeNorm)
               buoy_tidal=rInt_R(buoy_tidal_r_global,r,rscheme_oc)
               call heat_tidal_ave%compute(heat_tidal_r_global, n_calls, timePassed, timeNorm)
               heat_tidal=rInt_R(heat_tidal_r_global,r,rscheme_oc)
            end if
         end if

         if ( l_chemical_conv ) then
            call buo_chem_ave%compute(buoy_chem_r_global, n_calls, timePassed, &
                 &                    timeNorm)
            buoy_chem=rInt_R(buoy_chem_r_global,r,rscheme_oc)
         end if
      end if

      !-- Inner core:

      if ( l_cond_ic .and. l_mag ) then

         do n_r=1,n_r_ic_max
            r_ratio=r_ic(n_r)/r_ic(1)
            curlB2_rIC(n_r)=0.0_cp
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               Bh=(l+one)*O_r_ic(n_r)*aj_ic(lm,n_r)+dj_ic(lm,n_r)
               laplace=-ddb_ic(lm,n_r) - two*(l+one)*O_r_ic(n_r)*db_ic(lm,n_r)
               curlB2_rIC(n_r)=curlB2_rIC(n_r) + LFfac*opm*eScale*        &
               &               dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) *  ( &
               &               dLh(st_map%lm2(l,m))*O_r_ic2(n_r)*         &
               &               cc2real(aj_ic(lm,n_r),m) +                 &
               &               cc2real(Bh,m) + cc2real(laplace,m)       )
            end do
         end do    ! radial grid points

         call reduce_radial(curlB2_rIC, curlB2_rIC_global, 0)

         if ( rank == 0 ) then
            curlB2_IC=rIntIC(curlB2_rIC_global,n_r_ic_max,dr_fac_ic,chebt_ic)
         end if
      else
         if ( rank == 0 ) then
            curlB2_IC=0.0_cp
         end if
      end if  ! conducting inner core ?

      !-- Add up and correct:
      !  A correction is necessary when using the integral over
      !  (curl U)**2 to calculate the dippisated energy and inner core
      !  or mantle are allowed to rotate.
      !  The only flow component here is l=1,m=0, i.e. lm=2
      !  If the rotation rates of inner core of mantle are kept
      !  fixed, the viscous dissipation may actually be a power source if
      !  the radial derivatives drz10 are negative (positive) at the
      !  ICB (CMB)! The energy transfer is described by the very
      !  correction terms.
      if ( l_rot_ic ) then
         call send_lm_pair_to_master(z(:,n_r_icb),1,0,z10ICB)
         call send_lm_pair_to_master(dz(:,n_r_icb),1,0,drz10ICB)
      else
         z10ICB  =0.0_cp
         drz10ICB=0.0_cp
      end if
      if ( l_rot_ma ) then
         call send_lm_pair_to_master(z(:,n_r_cmb),1,0,z10CMB)
         call send_lm_pair_to_master(dz(:,n_r_cmb),1,0,drz10CMB)
      else
         z10CMB  =0.0_cp
         drz10CMB=0.0_cp
      end if

      if ( rank== 0 ) then
         if ( l_conv ) then
            viscDiss= -viscHeat
            !if ( l_rot_IC ) viscDiss=viscDiss - two*z10ICB*drz10ICB
            !if ( l_rot_MA ) viscDiss=viscDiss + two*z10CMB*drz10CMB
         else
            viscDiss=0.0_cp
         end if

         !--- If the inner core or mantle rotation rate is allowed to change due
         !    to viscous drag, the power transfer has to be taken into account:

         !-- Calculating viscous torques:
         if ( l_rot_ic .and. kbotv >= 2 ) then
            call get_viscous_torque(viscous_torque_ic,z10ICB,drz10ICB,r_icb, &
                 &                  beta(n_r_max),visc(n_r_max))
         else
            viscous_torque_ic=0.0_cp
         end if
         if ( l_rot_ma .and. ktopv >= 2 ) then
            call get_viscous_torque(viscous_torque_ma,z10CMB,drz10CMB,r_cmb, &
                 &                  beta(1),visc(1))
         else
            viscous_torque_ma=0.0_cp
         end if

         if ( l_rot_IC .and. .not. l_SRIC ) then
            if ( kbotv >= 2 ) then
               call get_viscous_torque(viscous_torque_ic,z10ICB,drz10ICB,r_icb,&
                    &                  beta(n_r_max),visc(n_r_max))
            else
               viscous_torque_ic=0.0_cp
            end if
            powerIC=omega_IC*(viscous_torque_ic+lorentz_torque_ic)
         else
            powerIC=0.0_cp
         end if
         if ( l_rot_MA ) then
            if ( ktopv >= 2 ) then
               call get_viscous_torque(viscous_torque_ma,z10CMB,drz10CMB,r_cmb, &
                    &                  beta(1),visc(1))
            else
               viscous_torque_ma=0.0_cp
            end if
            powerMA=omega_MA*(-viscous_torque_ma+lorentz_torque_ma)
         else
            powerMA=0.0_cp
         end if

         !--- Because the two systems are coupled only the total
         !--- ohmic dissipation is useful:
         ohmDiss=-curlB2-curlB2_IC

         powerDiffOld=powerDiff
         powerDiff   =(buoy+buoy_chem+powerIC+powerMA+viscDiss+ohmDiss)

         if ( abs(powerDiffOld) > 10.0_cp*epsilon(timePassed) ) then
            powerDiffT  =1.5_cp*powerDiff-half*powerDiffOld
            eDiffInt=eDiffInt+timePassed*timePassed*powerDiffT
            if ( l_save_out ) then
               open(newunit=n_power_file, file=power_file, status='unknown', &
               &    position='append')
            end if
            write(n_power_file,'(1P,ES20.12,10ES16.8)')   &
            &     time*tScale, buoy, buoy_chem,           &
            &     -two*z10ICB*drz10ICB,                   &
            &     two*z10CMB*drz10CMB, viscDiss,          &
            &     ohmDiss, powerMA, powerIC, powerDiff,   &
            &     eDiffInt/timeNorm, buoy_tidal, tidalpow,&
            &     visc_tidal, heat_tidal
            if ( l_save_out ) close(n_power_file)
         end if

         if ( l_stop_time ) then
            if ( l_heat ) call buo_ave%finalize_SD(timeNorm)
            if (l_tidal) then
               call buo_tidal_ave%finalize_SD(timeNorm)
               call visc_tidal_ave%finalize_SD(timeNorm)
               call tidal_ave%finalize_SD(timeNorm)
               call heat_tidal_ave%finalize_SD(timeNorm)
            end if
            if ( l_chemical_conv ) call buo_chem_ave%finalize_SD(timeNorm)
            if ( l_mag )  call ohm_ave%finalize_SD(timeNorm)
            if ( l_conv ) call visc_ave%finalize_SD(timeNorm)

            fileName='powerR.'//tag
            open(newunit=fileHandle, file=fileName, status='unknown')
            do n_r=1,n_r_max
               write(fileHandle,'(ES20.10,4ES15.7,4ES13.5)')                     &
               &     r(n_r),round_off(buo_ave%mean(n_r),maxval(buo_ave%mean)),   &
               &     round_off(buo_chem_ave%mean(n_r),maxval(buo_chem_ave%mean)),&
               &     round_off(visc_ave%mean(n_r),maxval(visc_ave%mean)),        &
               &     round_off(ohm_ave%mean(n_r),maxval(ohm_ave%mean)),          &
               &     round_off(buo_ave%SD(n_r),maxval(buo_ave%SD)),              &
               &     round_off(buo_chem_ave%SD(n_r),maxval(buo_chem_ave%SD)),    &
               &     round_off(visc_ave%SD(n_r),maxval(visc_ave%SD)),            &
               &     round_off(ohm_ave%SD(n_r),maxval(ohm_ave%SD)), &
               &     round_off(buo_tidal_ave%mean(n_r),maxval(buo_tidal_ave%mean)), &
               &     round_off(tidal_ave%mean(n_r),maxval(tidal_ave%mean)),&
               &     round_off(visc_tidal_ave%mean(n_r),maxval(visc_tidal_ave%mean)),&
               &     round_off(heat_tidal_ave%mean(n_r),maxval(heat_tidal_ave%mean)),&
               &     round_off(buo_tidal_ave%SD(n_r),maxval(buo_tidal_ave%SD)), &
               &     round_off(tidal_ave%SD(n_r),maxval(tidal_ave%SD)),&
               &     round_off(visc_tidal_ave%SD(n_r),maxval(visc_tidal_ave%SD)),&
               &     round_off(heat_tidal_ave%SD(n_r),maxval(heat_tidal_ave%SD))

            end do
            close(fileHandle)
         end if
      end if

  end subroutine get_power
!----------------------------------------------------------------------------
   subroutine get_visc_heat(vr,vt,vp,cvr,dvrdr,dvrdt,dvrdp,dvtdr,&
              &             dvtdp,dvpdr,dvpdp,nR)
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

      !-- Local variables:
      integer :: nTheta,nPhi
      real(cp) :: vischeat,csn2,phiNorm,viscAS

      phiNorm=two*pi/real(n_phi_max,cp)
      viscAS=0.0_cp

      !$omp parallel do default(shared)         &
      !$omp& private(nTheta,csn2,nPhi,vischeat) &
      !$omp& reduction(+:viscAS)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            csn2     =cosn_theta_E2(nTheta)

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
            &                        O_sin_theta_E2(nTheta) * ( &
            &         ( r(nR)*             dvtdr(nTheta,nPhi) - &
            &          (two+beta(nR)*r(nR))*  vt(nTheta,nPhi) + & ! (4)
            &     or1(nR)*           dvrdt(nTheta,nPhi) )**2  + &
            &         ( r(nR)*             dvpdr(nTheta,nPhi) - &
            &          (two+beta(nR)*r(nR))*  vp(nTheta,nPhi) + & ! (5)
            &    or1(nR)*            dvrdp(nTheta,nPhi) )**2 )- &
            &    two*third*(  beta(nR)*     vr(nTheta,nPhi) )**2 )

            viscAS=viscAS+phiNorm*gauss(nTheta)*viscHeat
         end do
      end do
      !$omp end parallel do

      viscASr(nR)=viscAS

   end subroutine get_visc_heat
!----------------------------------------------------------------------------
end module power
