module preCalculations
   !
   ! This module is used to handle some pre-calculations of constants (moment of inertia,
   ! mass, volumes), determine the timesteps for I/O and fix boundary values for
   ! temperature/entropy and chemical composition
   !

   use iso_fortran_env, only: output_unit
   use constants
   use num_param
   use output_data
   use precision_mod
   use truncation, only: n_r_max, l_max, minc, n_r_ic_max, nalias, &
       &                 n_cheb_ic_max, m_max, n_cheb_max, m_min,  &
       &                 lm_max, n_phi_max, n_theta_max
   use init_fields, only: bots, tops, s_bot, s_top, n_s_bounds,    &
       &                  l_reset_t, topxi, botxi, xi_bot, xi_top, &
       &                  n_xi_bounds, omega_ma1, omegaOsz_ma1,    &
       &                  omega_ic1, omegaOsz_ic1
   use parallel_mod, only: rank, n_procs, rank_with_r_LCR
   use logic, only: l_mag, l_cond_ic, l_non_rot, l_mag_LF, l_newmap,     &
       &            l_anel, l_heat, l_anelastic_liquid,                  &
       &            l_cmb_field, l_save_out, l_TO, l_TOmovie, l_r_field, &
       &            l_movie, l_LCR, l_dt_cmb_field, l_non_adia,          &
       &            l_temperature_diff, l_chemical_conv, l_probe,        &
       &            l_precession, l_finite_diff, l_full_sphere
   use radial_data, only: radial_balance
   use radial_functions, only: rscheme_oc, temp0, r_CMB, ogrun,            &
       &                       r_surface, visc, or2, r, r_ICB, dLtemp0,    &
       &                       beta, rho0, rgrav, dbeta, alpha0,           &
       &                       dentropy0, sigma, lambda, dLkappa, kappa,   &
       &                       dLvisc, dLlambda, divKtemp0, radial,        &
       &                       transportProperties, l_R
   use physical_parameters, only: nVarEps, pr, prmag, ra, rascaled, ek,    &
       &                          ekscaled, opr, opm, o_sr, radratio,      &
       &                          sigma_ratio, CorFac, LFfac, BuoFac,      &
       &                          PolInd, nVarCond, nVarDiff, nVarVisc,    &
       &                          rho_ratio_ic, rho_ratio_ma, epsc, epsc0, &
       &                          ktops, kbots, interior_model, r_LCR,     &
       &                          n_r_LCR, mode, tmagcon, oek, Bn, imagcon,&
       &                          ktopxi, kbotxi, epscxi, epscxi0, sc, osc,&
       &                          ChemFac, raxi, Po, prec_angle
   use horizontal_data, only: horizontal
   use integration, only: rInt_R
   use useful, only: logWrite, abortRun
   use special, only: l_curr, fac_loop, loopRadRatio, amp_curr, Le, n_imp, &
       &              l_imp, l_radial_flow_bc,                             &
       &              ellipticity_cmb, ellip_fac_cmb,                      &
       &              ellipticity_icb, ellip_fac_icb,                      &
       &              tide_fac20, tide_fac22p, tide_fac22n,                &
       &              amp_tide
   use time_schemes, only: type_tscheme

   implicit none

   private

   public :: preCalc, preCalcTimes, writeInfo

contains

   subroutine preCalc(tscheme)
      !
      !  Purpose of this subroutine is to initialize the calc values,
      !  arrays, constants that are used all over the code.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !---- Local variables:
      real(cp) :: help,facIH
      real(cp) :: delmin,sr_top,si_top,sr_bot,si_bot
      real(cp) :: xir_top,xii_top,xir_bot,xii_bot
      real(cp) :: topconduc, botconduc
      integer :: n,n_r,l,m,l_bot,m_bot,l_top,m_top
      integer :: fileHandle, p
      character(len=76) :: fileName
      character(len=80) :: message
      real(cp) :: mom(n_r_max)
      real(cp) :: y20_norm, y22_norm

      !-- Determine scales depending on n_tScale,n_lScale :
      if ( n_tScale == 0 ) then
         !----- Viscous time scale:
         tScale=one
      else if ( n_tScale == 1 ) then
         !----- Magnetic time scale:
         tScale=one/prmag
      else if ( n_tScale == 2 ) then
         !----- Thermal time scale:
         tScale=one/pr
      else if ( n_tScale == 3 ) then
         !----- Rotational time scale:
         tScale=one/ek  ! or ekScaled ? (not defined yet...)
         if ( rank==0 ) then
            write(output_unit,*) 'Warning: rotational timescale, be sure to set dtmax large enough !'
         end if
      end if
      if ( n_lScale == 0 ) then
         !----- Outer Core:
         lScale=one
      else if ( n_lScale == 1 ) then
         !----- Total Core:
         lScale=(one-radratio)
      end if

      !---- Scale according to scdIFf:
      vScale  =lScale/tScale
      pScale  =tScale**2/lScale**2
      eScale  =vScale*vScale/enscale
      raScaled=ra/lScale**3
      ekScaled=ek*lScale**2

      if ( l_cond_ic ) O_sr=one/sigma_ratio

      opr=one/pr
      if ( l_mag ) then
         opm=one/prmag
      else
         opm=0.0_cp
      end if

      if ( l_chemical_conv ) osc=one/sc

      ! Note: CorFac is the factor in front of the Coriolis force. In the scaling
      !       used here (viscous time scale) this is simply the inverse Ekman num:
      !          CorFac=1/Ek
      !       In the non-rotating case there is no Coriolis force and
      !          CorFac=0
      !       LFfac is the factor the Lorentz-force in the Navier-Stokes
      !       equation is multiplied with. This depends on the magnetic
      !       scale chosen. In the standart rotating case the B-scale
      !       is sqrt(\rho \lambda \mu \Omega) where \rho is core density,
      !       \lambda is magnetic diffusivity, \mu is magnetic permeability
      !       and \Omega is rotation rate. In this scaling the square of the
      !       dimensionless magnetic field is identical with the Elasser number:
      !          Els=B^2 / (\rho\lambda\mu\Omega)
      !       The factor to be used in front of the LF is then
      !          LFfac=1/(Ek Pm).
      !       This can also be used to bring kinetic and magnetic energy to
      !       the same scale since Ek=Em/(Ek Pm).
      !       If the system is non rotating we chose the magnetic scale
      !       sqrt(\rho \mu) \nu / L where \nu is the viscous diffusivity
      !       and L=ro-ri is the chosen length scale. Then LFfac=1.
      !       Note that for kinematic dynamos the magnetic field strength
      !       is arbitrary. We nevertheless still use the same LFfac.

      if ( l_non_rot ) then
         CorFac=0.0_cp
         if ( l_mag .or. l_mag_LF ) then
            LFfac=one
         else
            LFfac=0.0_cp
         end if
      else
         CorFac=one/ekScaled
         if ( l_mag .or. l_mag_LF ) then
            LFfac=one/(ekScaled*prmag)
         else
            LFfac=0.0_cp
         end if
      end if

      oek = CorFac

      if ( l_precession ) CorFac = CorFac*(one+po*cos(prec_angle))

      ! Note: BuoFac is the factor used in front of the buoyancy force.
      !       In the scaling used here its
      !          BuoFac=Ra/Pr
      !       where Ra= ( g0 \alpha \delta T L^3)/(\nu \kappa)
      !       with g0 the CMB gravity, \alpha the thermal expansivity,
      !       \delta T the temperature scale, and \kappa the thermal
      !       diffusivity
      if ( l_heat ) then
         BuoFac=raScaled/pr
      else
         BuoFac=0.0_cp
      end if

      if ( l_chemical_conv ) then
         ChemFac=raxi/sc
      else
         ChemFac=0.0_cp
      end if

      dtMax  =dtMax/tScale
      if ( .not. l_non_rot ) dtMax=min(dtMax,tscheme%intfac*ekScaled)
      dtMin  =dtMax/1.0e6_cp

      !-- Calculate radial functions for all threads (chebs,r,.....):

      call radial()

      call transportProperties()

      if ( ( l_anel .or. l_non_adia ) .and. ( rank == 0 ) ) then
         ! Write the equilibrium setup in anel.tag
         fileName='anel.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(11a15)') 'radius', 'temp0', 'rho0', 'beta',       &
         &                          'dbeta', 'grav', 'ds0/dr', 'div(kgradT)', &
         &                          'alpha', 'ogrun', 'dLtemp0'
         do n_r=1,n_r_max
            write(fileHandle,'(11ES16.8)') r(n_r), temp0(n_r), rho0(n_r),    &
            &                             beta(n_r), dbeta(n_r), rgrav(n_r), &
            &                             dentropy0(n_r), divKtemp0(n_r),    &
            &                             alpha0(n_r), ogrun(n_r), dLtemp0(n_r)
         end do
         close(fileHandle)
      end if

      !-- Write radial profiles
      if ( l_mag .and. nVarCond > 0 ) then
         fileName='varCond.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(4a15)') 'radius', 'sigma', 'lambda', 'dLlambda'
         do n_r=n_r_max,1,-1
            write(fileHandle,'(4ES16.8)') r(n_r),sigma(n_r),lambda(n_r), &
            &                             dLlambda(n_r)
         end do
         close(fileHandle)
      end if

      if ( ( l_heat .and. nVarDiff > 0  .or. nVarVisc > 0) .and. ( rank == 0 ) ) then
         fileName='varDiff.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(5a15)') 'radius', 'conductivity', 'kappa', &
         &                          'dLkappa', 'Prandtl'
         do n_r=n_r_max,1,-1
            write(fileHandle,'(5ES16.8)') r(n_r),kappa(n_r)*rho0(n_r), &
            &                             kappa(n_r),dLkappa(n_r),     &
            &                             pr*visc(n_r)/kappa(n_r)
         end do
         close(fileHandle)
      end if

      if ( ( nVarVisc > 0 ) .and. (rank == 0) ) then
         fileName='varVisc.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(7a15)') 'radius', 'dynVisc', 'kinVisc', &
         &                          'dLvisc', 'Ekman', 'Prandtl', 'Pm'
         if ( l_mag ) then
            do n_r=n_r_max,1,-1
               write(fileHandle,'(7ES16.8)') r(n_r),visc(n_r)*rho0(n_r), &
               &                             visc(n_r),dLvisc(n_r),      &
               &                             ek*visc(n_r),               &
               &                             pr*visc(n_r)/kappa(n_r),    &
               &                             prmag*visc(n_r)/lambda(n_r)
            end do
         else
            do n_r=n_r_max,1,-1
               write(fileHandle,'(7ES16.8)') r(n_r),visc(n_r)*rho0(n_r), &
               &                             visc(n_r), dLvisc(n_r),     &
               &                             ek*visc(n_r),               &
               &                             pr*visc(n_r)/kappa(n_r),    &
               &                             prmag
            end do
         end if
         close(fileHandle)
      end if

      l_LCR=.false.
      n_r_LCR=0
      do n_r=1,n_r_max
         if ( r_LCR <= r(n_r)/r_CMB ) then
             l_LCR=.true.
             n_r_LCR=n_r
         end if
      end do
      if ( n_r_LCR == 1 ) then
         l_LCR=.false.
         n_r_LCR=0
      end if

      if ( l_LCR ) then ! Determine which ranks carries the radius with n_r_LCR
         do p=0,n_procs-1
            if ( radial_balance(p)%nStart <= n_r_LCR .and. &
            &    radial_balance(p)%nStop >= n_r_LCR ) then
               rank_with_r_LCR=p
               exit
            end if
         end do
      end if

      !-- Compute some constants:
      vol_ic=four*third*pi*r_icb**3             ! Inner core volume
      vol_oc=four*third*pi*(r_cmb**3-r_icb**3)  ! Outer core volume
      surf_cmb=four*pi*r_cmb**2                ! Surface of CMB

      !-- Initialize everything that has to do with the horizontal representation
      !   on all threads:
      call horizontal()

      !-- Computation of the average density (useful to compute Re and Rm)
      if ( l_anel ) then
         mom(:)=r(:)**2 * rho0(:)
         mass=four*pi/vol_oc*rInt_R(mom,r,rscheme_oc)
      else
         mass=one
      end if

      !-- Calculate auxiliary arrays containing effective Courant grid intervals:
      delxh2(:)=r(:)**2/real(l_R(:)*(l_R(:)+1),kind=cp)
      delxr2(1)      =(r(1)-r(2))**2
      delxr2(n_r_max)=(r(n_r_max-1)-r(n_r_max))**2
      do n_r=2,n_r_max-1
         delmin=min((r(n_r-1)-r(n_r)),(r(n_r)-r(n_r+1)))
         delxr2(n_r)=delmin*delmin
      end do

      !-- Constants used for rotating core or mantle:
      y10_norm=half*sqrt(three/pi)  ! y10=y10_norm * cos(theta)
      y11_norm=half*sqrt(1.5_cp/pi) ! y11=y11_norm * sin(theta)

      !----- Proportionality factor of (l=1,m=0) toroidal velocity potential
      !      and inner core rotation rate:
      c_z10_omega_ic=y10_norm*or2(n_r_max)/rho0(n_r_max)

      !----- Proportionality factor of (l=1,m=0) toroidal velocity potential
      !      and mantle rotation rate:
      c_z10_omega_ma=y10_norm*or2(1)/rho0(1)

      !----- Inner-core normalized moment of inertia:
      c_moi_ic=8.0_cp*pi/15.0_cp*r_icb**5*rho_ratio_ic*rho0(n_r_max)

      !----- Outer-core normalized moment of inertia:
      ! _moi_oc=8.0_cp*pi/15.0_cp*(r_cmb**5-r_icb**5) ! rho=cst
      mom(:)=r(:)**4 * rho0(:)
      c_moi_oc=8.0_cp*third*pi*rInt_R(mom,r,rscheme_oc)

      !----- Mantle normalized moment of inertia:
      c_moi_ma=8.0_cp*pi/15.0_cp*(r_surface**5-r_cmb**5)*rho_ratio_ma*rho0(1)

      !----- IC normalised moment of inertia / r_icb**4 * 3/(8 pi)
      c_dt_z10_ic=0.2_cp*r_icb*rho_ratio_ic*rho0(n_r_max)

      !----- Mantle normalised moment of inertia / r_cmb**4 * 3/(8 pi)
      c_dt_z10_ma=0.2_cp*r_cmb*rho_ratio_ma * ( (r_surface/r_cmb)**5 - one )

      !----- Proportionality factor for ic lorentz_torque as used in
      !      ic torque-equation (z10):
      c_lorentz_ic=0.25_cp*sqrt(three/pi)*or2(n_r_max)

      !----- Proportionality factor for mantle lorentz_torque as used in
      !      mantle torque-equation (z10):
      c_lorentz_ma=0.25_cp*sqrt(three/pi)*or2(1)

      !-- Set thermal boundary conditions for fixed temp. on both boundaries:
      !----- Extract tops and bots
      if ( l_heat ) then

         !-- Renormalisation of heat flow coeffs and heat source.
         !      This has been copied from the Christensen and
         !      means the we use a different normalisation of
         !      input tops,bots and epsc0 than in the code.
         !      Spherical harmonics of the input tops, bots and epsc0
         !      are normalised so that the integral of Y*Y(c.c.)
         !      (c.c. stands for conjugate complex)
         !      over a spherical surface of radius 1 gives 4 Pi.
         !      In the code we use totally normalized spherical harmonic,
         !      i.e. the integral of Y*Y(c.c.) over a spherical surface
         !      of radius 1 is unity
         epsc=epsc0*sq4pi

         do m=0,m_max,minc
            do l=m,l_max
               bots(l,m)=zero
               tops(l,m)=zero
               do n=1,n_s_bounds
                  l_bot =int(s_bot(4*n-3))
                  m_bot =int(s_bot(4*n-2))
                  sr_bot=s_bot(4*n-1)
                  si_bot=s_bot(4*n)
                  l_top =int(s_top(4*n-3))
                  m_top =int(s_top(4*n-2))
                  sr_top=s_top(4*n-1)
                  si_top=s_top(4*n)
                  if ( l_bot == l .and. m_bot == m .and. &
                       cmplx(sr_bot,si_bot,kind=cp) /= zero ) then
                     if ( m == 0 ) si_bot=0.0_cp
                     bots(l,m)=sq4pi*cmplx(sr_bot,si_bot,kind=cp)
                     if ( kbots == 2 ) bots(l,m)=bots(l,m)*lScale
                  end if
                  if ( l_top == l .and. m_top == m .and. &
                       cmplx(sr_top,si_top,kind=cp) /= zero ) then
                     if ( m == 0 ) si_top=0.0_cp
                     tops(l,m)=sq4pi*cmplx(sr_top,si_top,kind=cp)
                     if ( ktops == 2 ) tops(l,m)=tops(l,m)*lScale
                  end if
               end do
            end do
         end do

         if ( nVarEps==0 ) then
            facIH=vol_oc
         else if ( nVarEps==1 ) then
            facIH=mass*vol_oc
         else if ( nVarEps == 2 ) then
            mom(:)=r(:)**2*rho0(:)*temp0(:)
            facIH=four*pi*rInt_R(mom,r,rscheme_oc)
         else if ( nVarEps == 3 ) then
            mom(:)=r(:)**2*rho0(:)**2*temp0(:)**(-3)*exp(-Bn/temp0(:))
            facIH=four*pi*rInt_R(mom,r,rscheme_oc)
         end if

         if ( l_temperature_diff .or. l_anelastic_liquid ) then
            topconduc = rho0(1)*kappa(1)
            botconduc = rho0(n_r_max)*kappa(n_r_max)
         else
            topconduc = rho0(1)*kappa(1)*temp0(1)
            botconduc = rho0(n_r_max)*kappa(n_r_max)*temp0(n_r_max)
         end if

         if ( ktops == 1 .and. kbots == 1 ) then ! Fixed entropy

            tops(0,0)=0.0_cp
            bots(0,0)=sq4pi

         else if ( ktops == 3 .and. kbots == 3 ) then ! Fixed temperature contrast

            tops(0,0)=0.0_cp
            bots(0,0)=sq4pi

         else if ( (ktops==2 .and. kbots==2) .or. (ktops == 4 .and. kbots==4) ) then

            if ( real(bots(0,0)) > 0.0_cp ) then
               write(output_unit,*)
               write(output_unit,*) '! NOTE: you have supplied'
               write(output_unit,*) '! s_bot(l=0,m=0)>0 which '
               write(output_unit,*) '! means there is a heat '
               write(output_unit,*) '! flux into the inner core.'
               write(output_unit,*) '! This is unrealistic!'
               write(output_unit,*) '! Use s_bot(l=0,m=0)<0 !'
               call abortRun('Stop run in preCalc')
            end if

            !--- |epsc0|=1 signifies that the heat production rate is used as
            !    temperature scale! I make sure here that the heat flux through
            !    inner and outer boundary, bots and tops, balance the sources.
            if ( abs(epsc0) == one ) then

               !--- Make sure that all the heat comes out at CMB
               if ( tops(0,0) == 0.0_cp ) then

                  !--- Compensate by flux from ICB:
                  !    all over the core :
                  if ( epsc0 >= 0 ) then
                     write(output_unit,*) '! NOTE: when the flux through the '
                     write(output_unit,*) '! outer boundary is zero, sinks in'
                     write(output_unit,*) '! the outer core need to balance  '
                     write(output_unit,*) '! the flux from the ICB. Thus we  '
                     write(output_unit,*) '! need epsc<0 !                   '
                     call abortRun('Stop run in preCalc')
                  end if
                  bots(0,0)=epsc*pr*facIH/(four*pi*r_icb**2 * botconduc )
                  call logWrite('! CMB heat flux set to balance volume sources!')

               else if ( tops(0,0) /= 0.0_cp .and. bots(0,0) == 0.0_cp ) then

                  !--- Correct tops to balance inner sources/sinks:
                  if ( epsc0 <= 0 .and. bots(0,0) == 0.0_cp  ) then
                     write(output_unit,*) '! NOTE: when the flux through the '
                     write(output_unit,*) '! ICB is zero we need sources in  '
                     write(output_unit,*) '! the outer core which means      '
                     write(output_unit,*) '! epsc0>0.                        '
                     call abortRun('Stop run in preCalc')
                  end if
                  help=real(tops(0,0))
                  if ( abs(real(tops(0,0))) == sq4pi ) &
                  &    call logWrite('! You intend to use the CMB flux as buoy. scale??')
                  tops(0,0)=-facIH*epsc*pr/(four*pi*r_cmb**2) + &
                  &          radratio**2*botconduc*bots(0,0)
                  if ( real(tops(0,0)) /= help ) &
                     call logWrite('!!!! WARNING: CMB heat flux corrected !!!!')

               else if ( tops(0,0) /= 0.0_cp .and. bots(0,0) /= 0.0_cp ) then

                  !--- Correct tops to balance inner sources/sinks:
                  if ( epsc0 <= 0 .and. bots(0,0) == 0.0_cp  ) then
                     write(output_unit,*) '! NOTE: when the flux through the '
                     write(output_unit,*) '! ICB is zero we need sources in  '
                     write(output_unit,*) '! the outer core which means      '
                     write(output_unit,*) '! epsc0>0.                        '
                     call abortRun('Stop run in preCalc')
                  end if
                  help=four*pi*opr*facIH *                   &
                  &    (r_icb**2*real(bots(0,0))*botconduc - &
                  &     r_cmb**2*real(tops(0,0))*topconduc)
                  if ( help /= epsc ) then
                     write(output_unit,*) '! NOTE: when flux BC through the '
                     write(output_unit,*) '! ICB and CMB are used the sources '
                     write(output_unit,*) '! have to balance the total flux.'
                     call abortRun('Stop run in preCalc')
                  end if

               end if

            else if ( epsc0 == 0.0_cp .and. ( tops(0,0) /= 0.0_cp .or. &
            &                               bots(0,0) /= 0.0_cp ) ) then
               !--- Correct epsc0 to balance the difference between
               !    flux through the inner and outer boundary:
               epsc=four*pi/pr/facIH *                    &
               &    (r_icb**2*real(bots(0,0))*botconduc - &
               &     r_cmb**2*real(tops(0,0))*topconduc)
               call logWrite('! Sources introduced to balance surface heat flux!')
               write(message,'(''!      epsc0='',ES16.6)') epsc/sq4pi
               call logWrite(message)
            else if ( epsc0 /= 0.0_cp .and. tops(0,0) /= 0.0_cp .and. &
            &                             bots(0,0) /= 0.0_cp ) then
               help=four*pi/pr/facIH *                    &
               &    (r_icb**2*real(bots(0,0))*botconduc - &
               &     r_cmb**2*real(tops(0,0))*topconduc)
               if ( help /= epsc ) then
                  write(output_unit,*) '! NOTE: when flux BC through the '
                  write(output_unit,*) '! ICB and/or CMB is used the sources '
                  write(output_unit,*) '! have to balance it.'
                  call abortRun('Stop run in preCalc')
               end if
            end if

         end if

         if ( l_non_adia ) then
            bots(0,0)=0.0_cp
            tops(0,0)=0.0_cp
            epsc=0.0_cp
            !epsc=four*pi/pr/epsS/vol_oc *          &
            !     (r_icb**2*dtemp0(n_r_max)*rho0(n_r_max)*kappa(n_r_max) - &
            !      r_cmb**2*dtemp0(1)*rho0(1)*kappa(1))*sq4pi
         end if
         if ( ktops == 1 ) then
            write(message,'(''! Const. entropy at outer boundary S ='',ES16.6)') &
            &     real(tops(0,0))/sq4pi
            call logWrite(message)
         else if ( ktops == 2 ) then
            help=surf_cmb*topconduc*real(tops(0,0))/sq4pi
            write(message,'(''! Const. total outer boundary entropy flux    ='',ES16.6)') &
            &     help
            call logWrite(message)
         else if ( ktops == 3 ) then
            write(message,'(''! Const. temp. at outer boundary S ='',ES16.6)') &
            &     real(tops(0,0))/sq4pi
            call logWrite(message)
         else if ( ktops == 4 ) then
            help=surf_cmb*topconduc*real(tops(0,0))/sq4pi
            write(message,'(''! Const. total outer boundary temp. flux    ='',ES16.6)') &
            &     help
            call logWrite(message)
         end if
         if ( kbots == 1 ) then
            write(message,'(''! Const. entropy at inner boundary S ='',ES16.6)') &
            &     real(bots(0,0))/sq4pi
            call logWrite(message)
         else if ( kbots == 2 ) then
            help=surf_cmb*radratio**2*botconduc*real(bots(0,0))/sq4pi
            write(message, '(''! Const. total inner boundary entropy flux    ='',ES16.6)')&
            &     help
            call logWrite(message)
         else if ( kbots == 3 ) then
            write(message,'(''! Const. temp. at inner boundary S ='',ES16.6)') &
            &     real(bots(0,0))/sq4pi
            call logWrite(message)
         else if ( kbots == 4 ) then
            help=surf_cmb*radratio**2*botconduc*real(bots(0,0))/sq4pi
            write(message, '(''! Const. total inner boundary temp. flux    ='',ES16.6)') &
            &     help
            call logWrite(message)
         end if
         help=facIH*pr*epsc/sq4pi
         write(message,'(''! Total vol. buoy. source ='',ES16.6)') help
         call logWrite(message)

      end if

      !-- Set  boundary conditions for chemical composition
      if ( l_chemical_conv ) then
         epscxi=epscxi0*sq4pi

         do m=0,m_max,minc
            do l=m,l_max
               botxi(l,m)=zero
               topxi(l,m)=zero
               do n=1,n_xi_bounds
                  l_bot =int(xi_bot(4*n-3))
                  m_bot =int(xi_bot(4*n-2))
                  xir_bot=xi_bot(4*n-1)
                  xii_bot=xi_bot(4*n)
                  l_top =int(xi_top(4*n-3))
                  m_top =int(xi_top(4*n-2))
                  xir_top=xi_top(4*n-1)
                  xii_top=xi_top(4*n)
                  if ( l_bot == l .and. m_bot == m .and. &
                       cmplx(xir_bot,xii_bot,kind=cp) /= zero ) then
                     if ( m == 0 ) xii_bot=0.0_cp
                     botxi(l,m)=sq4pi*cmplx(xir_bot,xii_bot,kind=cp)
                     if ( kbotxi == 2 ) botxi(l,m)=botxi(l,m)*lScale
                  end if
                  if ( l_top == l .and. m_top == m .and. &
                       cmplx(xir_top,xii_top,kind=cp) /= zero ) then
                     if ( m == 0 ) xii_top=0.0_cp
                     topxi(l,m)=sq4pi*cmplx(xir_top,xii_top,kind=cp)
                     if ( ktopxi == 2 ) topxi(l,m)=topxi(l,m)*lScale
                  end if
               end do
            end do
         end do

         facIH=vol_oc
         topconduc = rho0(1)
         botconduc = rho0(n_r_max)

         if ( ktopxi == 1 .and. kbotxi == 1 ) then ! Fixed chemical comp

            topxi(0,0)=0.0_cp
            botxi(0,0)=sq4pi

         else if ( (ktopxi==2 .and. kbotxi==2) ) then

            if ( real(botxi(0,0)) > 0.0_cp ) then
               write(output_unit,*)
               write(output_unit,*) '! NOTE: you have supplied'
               write(output_unit,*) '! xi_bot(l=0,m=0)>0 which '
               write(output_unit,*) '! means there is a composition '
               write(output_unit,*) '! flux into the inner core.'
               write(output_unit,*) '! This is unrealistic!'
               write(output_unit,*) '! Use xi_bot(l=0,m=0)<0 !'
               call abortRun('Stop run in preCalc')
            end if

            if ( abs(epscxi0) == one ) then

               if ( topxi(0,0) == 0.0_cp ) then

                  !--- Compensate by flux from ICB:
                  !    all over the core :
                  if ( epscxi0 >= 0 ) then
                     write(output_unit,*) '! NOTE: when the flux through the '
                     write(output_unit,*) '! outer boundary is zero, sinks in'
                     write(output_unit,*) '! the outer core need to balance  '
                     write(output_unit,*) '! the flux from the ICB. Thus we  '
                     write(output_unit,*) '! need epscxi<0 !                   '
                     call abortRun('Stop run in preCalc')
                  end if
                  botxi(0,0)=epscxi*sc*facIH/(four*pi*r_icb**2*botconduc)
                  call logWrite('! CMB heat flux set to balance volume sources!')

               else if ( topxi(0,0) /= 0.0_cp .and. botxi(0,0) == 0.0_cp ) then

                  !--- Correct topxi to balance inner sources/sinks:
                  if ( epscxi0 <= 0 .and. botxi(0,0) == 0.0_cp  ) then
                     write(output_unit,*) '! NOTE: when the flux through the '
                     write(output_unit,*) '! ICB is zero we need sources in  '
                     write(output_unit,*) '! the outer core which means      '
                     write(output_unit,*) '! epscxi0>0.                        '
                     call abortRun('Stop run in preCalc')
                  end if
                  if ( abs(real(topxi(0,0))) == sq4pi ) &
                  &    call logWrite('! You intend to use the CMB flux as buoy. scale??')
                  topxi(0,0)=-facIH*epscxi*sc/(four*pi*r_cmb**2) + &
                  &          radratio**2*botxi(0,0)*botconduc
                  if ( topxi(0,0) /= help ) &
                     call logWrite('!!!! WARNING: CMB composition flux corrected !!!!')

               else if ( topxi(0,0) /= 0.0_cp .and. botxi(0,0) /= 0.0_cp ) then

                  !--- Correct tops to balance inner sources/sinks:
                  if ( epscxi0 <= 0 .and. botxi(0,0) == 0.0_cp  ) then
                     write(output_unit,*) '! NOTE: when the flux through the '
                     write(output_unit,*) '! ICB is zero we need sources in  '
                     write(output_unit,*) '! the outer core which means      '
                     write(output_unit,*) '! epscxi0>0.                      '
                     call abortRun('Stop run in preCalc')
                  end if
                  help=four*pi/sc/facIH *                     &
                  &    (r_icb**2*real(botxi(0,0))*botconduc - &
                  &     r_cmb**2*real(topxi(0,0))*topconduc)
                  if ( help /= epscxi ) then
                     write(output_unit,*) '! NOTE: when flux BC through the '
                     write(output_unit,*) '! ICB and CMB are used the sources '
                     write(output_unit,*) '! have to balance the total flux.'
                     call abortRun('Stop run in preCalc')
                  end if

               end if

            else if ( epscxi0 == 0.0_cp .and. ( topxi(0,0) /= 0.0_cp .or. &
                                            botxi(0,0) /= 0.0_cp ) ) then
               !--- Correct epscxi0 to balance the difference between
               !    flux through the inner and outer boundary:
               epscxi=four*pi/sc/facIH *                     &
               &      (r_icb**2*real(botxi(0,0))*botconduc - &
               &       r_cmb**2*real(topxi(0,0))*topconduc)
               call logWrite('! Sources introduced to balance surface Fickian flux!')
               write(message,'(''!      epscxi0='',ES16.6)') epscxi/sq4pi
               call logWrite(message)
            else if ( epscxi0 /= 0.0_cp .and. topxi(0,0) /= 0.0_cp .and. &
            &                             botxi(0,0) /= 0.0_cp ) then
               help=four*pi/sc/facIH *                     &
               &    (r_icb**2*real(botxi(0,0))*botconduc - &
               &     r_cmb**2*real(topxi(0,0))*topconduc)
               if ( help /= epscxi ) then
                  write(output_unit,*) '! NOTE: when flux BC through the '
                  write(output_unit,*) '! ICB and/or CMB is used the sources '
                  write(output_unit,*) '! have to balance it.'
                  call abortRun('Stop run in preCalc')
               end if
            end if

         end if

         if ( ktopxi == 1 ) then
            write(message,'(''! Constant comp. at CMB T ='',ES16.6)') &
            &     real(topxi(0,0))/sq4pi
            call logWrite(message)
         else if ( ktopxi == 2 ) then
            help=surf_cmb*topconduc*real(topxi(0,0))/sq4pi
            write(message,'(''! Const. total CMB comp. flux    ='',ES16.6)') help
            call logWrite(message)
         end if
         if ( kbotxi == 1 ) then
            write(message,'(''! Constant comp. at ICB T ='',ES16.6)') &
            &     real(botxi(0,0))/sq4pi
            call logWrite(message)
         else if ( kbotxi == 2 ) then
            help=surf_cmb*radratio**2*botconduc*real(botxi(0,0))/sq4pi
            write(message, '(''! Const. total ICB comp. flux    ='',ES16.6)') help
            call logWrite(message)
         end if
         help=facIH*sc*epscxi/sq4pi
         write(message,'(''! Total vol. comp. source ='',ES16.6)') help
         call logWrite(message)

      end if

      !--- Compute fac_loop for current carrying loop
      if ( l_curr ) then
         allocate(fac_loop(l_max))

         do l=1,l_max
            fac_loop(l)=0.0_cp
            if (mod(l,2)/=0) then
               if ( l==1 ) then
                  fac_loop(l)= one
               else
                  fac_loop(l)= -fac_loop(l-2)*loopRadRatio**2*real(l,kind=cp)/ &
                  &            real(l-1,kind=cp)
               end if
            end if
         end do

         if ( l_non_rot ) then
            amp_curr = Le
         else
            amp_curr = Le * sqrt(prmag/ek)
         end if
      end if

      if ( (n_imp == 3 .or. n_imp == 4 .or. n_imp == 7) .and. ( l_imp /= 1 ) ) then
         call abortRun('l_imp /= 1 not implemented for this imposed field setup!')
      end if

      if ( ((imagcon /= 0) .or. l_curr .or. (n_imp > 1)) .and. l_LCR ) then
         call abortRun('LCR not compatible with imposed field!')
      end if

      !-- Define n_s_max for cylindral grid
      n_s_max = n_r_max+int(r_icb*n_r_max)
      n_s_max = int(sDens*n_s_max)

      if (l_radial_flow_bc) then
         if ( ellipticity_cmb /= 0.0_cp ) then
            ellip_fac_cmb=-two*r_cmb**3*ellipticity_cmb*omega_ma1*omegaOsz_ma1
         else
            ellip_fac_cmb=0.0_cp
         end if

         if ( ellipticity_icb /= 0.0_cp ) then
            ellip_fac_icb=-two*r_icb**3*ellipticity_icb*omega_ic1*omegaOsz_ic1
         else
            ellip_fac_icb=0.0_cp
         end if

         if ( amp_tide /= 0.0_cp ) then
            y20_norm = 0.5_cp  * sqrt(5.0_cp/pi)
            y22_norm = 0.25_cp * sqrt(7.5_cp/pi)
            tide_fac20  = amp_tide / y20_norm * r_cmb**2 ! (2,0,1) mode of Ogilvie 2014
            tide_fac22p = half * amp_tide / y22_norm / sqrt(6.0_cp) * r_cmb**2 ! Needs a half factor, (2,2,1) mode
            tide_fac22n = -7.0_cp * tide_fac22p                                ! Half factor carried over, (2,2,3) mode,
                                                                               ! has opposite sign to that of the other two (Polfliet & Smeyers, 1990)
         else
            tide_fac20  = 0.0_cp
            tide_fac22p = 0.0_cp
            tide_fac22n = 0.0_cp
         end if
      end if
   end subroutine preCalc
!-------------------------------------------------------------------------------
   subroutine preCalcTimes(time,n_time_step)
      !
      !  Precalc. after time, time and dt has been read from startfile.
      !

      !-- Output variables
      real(cp), intent(inout) :: time ! Current time
      integer,  intent(inout) :: n_time_step ! Index of time loop

      !----- Set time step:
      if ( l_reset_t ) then
         time=0.0_cp
         n_time_step=0
      end if

      tmagcon=tmagcon+time

      !-- Get output times:
      call get_hit_times(t_graph,t_graph_start,t_graph_stop,dt_graph,  &
           &             n_graphs,n_graph_step,'graph',time,tScale)

      call get_hit_times(t_pot,t_pot_start,t_pot_stop,dt_pot,    &
           &             n_pots,n_pot_step,'pot',time,tScale)

      call get_hit_times(t_rst,t_rst_start,t_rst_stop,dt_rst,    &
           &             n_rsts,n_rst_step,'rst',time,tScale)

      call get_hit_times(t_log,t_log_start,t_log_stop,dt_log,    &
           &             n_logs,n_log_step,'log',time,tScale)

      call get_hit_times(t_spec,t_spec_start,t_spec_stop,dt_spec,  &
           &             n_specs,n_spec_step,'spec',time,tScale)

      if ( l_probe ) then
         call get_hit_times(t_probe,t_probe_start,t_probe_stop,dt_probe, &
              &             n_probe_out,n_probe_step,'probe',time,tScale)
      end if

      if ( l_cmb_field ) then
         call get_hit_times(t_cmb,t_cmb_start,t_cmb_stop,dt_cmb,n_cmbs,&
              &             n_cmb_step,'cmb',time,tScale)
      end if
      l_dt_cmb_field=l_dt_cmb_field .and. l_cmb_field

      if ( l_r_field ) then
         call get_hit_times(t_r_field,t_r_field_start,t_r_field_stop,            &
              &             dt_r_field,n_r_fields,n_r_field_step,'r_field',time, &
              &             tScale)
      end iF

      if ( l_movie ) then
         call get_hit_times(t_movie,t_movie_start,t_movie_stop,dt_movie,&
              &             n_movie_frames,n_movie_step,'movie',time,tScale)
      end if

      if ( l_TO ) then
         call get_hit_times(t_TO,t_TO_start,t_TO_stop,dt_TO,n_TOs, &
              &             n_TO_step,'TO',time,tScale)
      end if

      if ( l_TOmovie ) then
         call get_hit_times(t_TOmovie,t_TOmovie_start,t_TOmovie_stop,            &
              &             dt_TOmovie,n_TOmovie_frames,n_TOmovie_step,'TOmovie',&
              &             time,tScale)
      end if

   end subroutine preCalcTimes
!-------------------------------------------------------------------------------
   subroutine get_hit_times(t,t_start,t_stop,dt,n_tot,n_step,string,time,&
              &             tScale)
      !
      ! This subroutine checks whether any specific times t(*) are given
      ! on input. If so, it returns their number n_r and sets l_t to true.
      ! If not, t(:) may also be defined by giving a time step dt or a
      ! number n_tot of desired output times and ``t_stop>t_start``.
      !

      !-- Input variables:
      real(cp),         intent(in) :: time       ! Time of start file
      real(cp),         intent(in) :: tScale     ! Scale unit for time
      character(len=*), intent(in) :: string
      integer,               intent(inout) :: n_tot   ! Number of output (times) if no times defined
      integer,               intent(inout) :: n_step  ! Ouput step in no. of time steps
      real(cp), allocatable, intent(inout) :: t(:)    ! Times for output
      real(cp),              intent(inout) :: t_start ! Starting time for output
      real(cp),              intent(inout) :: t_stop  ! Stop time for output
      real(cp),              intent(inout) :: dt      ! Time step for output

      !-- Local variables:
      logical :: l_t
      real(cp) :: tmp(size(t))
      integer :: n, n_t

      t_start=t_start/tScale
      t_stop =t_stop/tScale
      dt     =dt/tScale

      if ( n_tot /= 0 .and. n_step /= 0 ) then
         write(output_unit,*)
         write(output_unit,*) '! You have to either provide the total'
         write(output_unit,*) '! number or the step for output:'
         write(output_unit,'(A,2(A,I10))') string, "n_tot = ",n_tot,", n_step = ",n_step
         write(output_unit,*) '! I set the step width to zero!'
         n_step=0
      end if

      !-- Check whether any time is given explicitly:
      l_t=.false.
      n_t=0
      where ( t >= 0.0_cp ) t=t/tScale
      if ( maxval(t(:)) >= 0.0_cp ) l_t=.true.; n_t=count(t(:)>=0.0)

      !-- Properly reallocate at the right size
      if ( l_t ) then
         tmp(:)=t(:)
         deallocate(t)
         allocate(t(n_t))
         t(:)=tmp(1:n_t)
      end if

      !-- Check times should be constructed:
      if ( t_start < time ) t_start=time
      if ( .not. l_t .and. ( dt > 0.0_cp .or. ( n_tot > 0 .and. t_stop > t_start ) ) ) then
         if ( n_tot > 0 .and. dt > 0.0_cp ) then
            n_t=n_tot
            n_tot=0 ! This is to ensure that time array is used later in l_correct_step
         else if ( dt > 0.0_cp ) then ! In that case t(:) is an array with one element
            n_t=1
         else if ( n_tot > 0 ) then
            n_t=n_tot
            dt=(t_stop-t_start)/real(n_t-1,kind=cp)
            n_tot=0 ! This is to ensure that time array is used later in l_correct_step
         end if
         l_t=.true.

         deallocate(t) ! Deallocate to get the proper size
         allocate(t(n_t))
         t(:) = [(t_start+(n-1)*dt, n=1,n_t)]
      end if

      !-- In case only n_step or n_tot is specified such that l_t=.false.
      if ( .not. l_t ) then
         deallocate(t)
         allocate(t(1))
         t(1)=-one
      end if

   end subroutine get_hit_times
!------------------------------------------------------------------------------
   subroutine writeInfo(n_out)
      !
      !  Purpose of this subroutine is to write the namelist to
      !  file unit n_out. This file has to be open before calling this
      !  routine.
      !

      !-- Input variable:
      integer, intent(in) :: n_out ! Output unit

      if ( rank == 0 ) then

         !-- Output of mode:
         write(n_out,*)
         if ( mode == 0 ) then
            write(n_out,*) '! Self consistent dynamo integration.'
         else if ( mode == 1 ) then
            write(n_out,*) '! Convection integration.'
         else if ( mode == 2 ) then
            write(n_out,*) '! Kinematic dynamo integration.'
         else if ( mode == 3 ) then
            write(n_out,*) '! Magnetic decay modes.'
         else if ( mode == 4 ) then
            write(n_out,*) '! Magneto convection.'
         else if ( mode == 5 ) then
            write(n_out,*) '! Linear onset of convection.'
         else if ( mode == 6 ) then
            write(n_out,*) '! Self consistent dynamo integration without LF.'
         else if ( mode == 7 ) then
            write(n_out,*) '! Super-rotating IC, no convection, no dynamo.'
         else if ( mode == 8 ) then
            write(n_out,*) '! Super-rotating IC, no convection, dynamo.'
         else if ( mode == 9 ) then
            write(n_out,*) '! Super-rotating IC, no convection, dynamo, no LF.'
         else if ( mode == 10 ) then
            write(n_out,*) '! Super-rotating IC, no advection, no convection, no dynamo.'
         else if ( mode == 11 ) then
            write(n_out,*) '! Viscous flow, no inertia, no rotation, no dynamo.'
         end if
      end if
      if ( mode > 11 ) then
         call abortRun('Mode > 11 not implemented !')
      end if

      if (rank == 0) then
         !-- Output of name lists:
         write(n_out, '('' ! Normalized OC moment of inertia:'',ES14.6)') c_moi_oc
         write(n_out, '('' ! Normalized IC moment of inertia:'',ES14.6)') c_moi_ic
         write(n_out, '('' ! Normalized MA moment of inertia:'',ES14.6)') c_moi_ma
         write(n_out, '('' ! Normalized IC volume           :'',ES14.6)') vol_ic
         write(n_out, '('' ! Normalized OC volume           :'',ES14.6)') vol_oc
         write(n_out, '('' ! Normalized IC surface          :'',ES14.6)')  &
                       surf_cmb*radratio**2
         write(n_out, '('' ! Normalized OC surface          :'',ES14.6)') surf_cmb
         write(n_out,*)
         write(n_out,*) '! Grid parameters:'
         write(n_out,'(''  n_r_max      ='',i6, &
         &        '' = number of radial grid points'')') n_r_max
         if ( .not. l_finite_diff ) then
            write(n_out,'(''  n_cheb_max   ='',i6)') n_cheb_max
            write(n_out,'(''  max cheb deg.='',i6)') n_cheb_max-1
         end if
         write(n_out,'(''  n_phi_max    ='',i6, &
         &        '' = no of longitude grid points'')') n_phi_max
         write(n_out,'(''  n_theta_max  ='',i6, &
         &        '' = no of latitude grid points'')') n_theta_max
         if ( .not. l_full_sphere ) then
            write(n_out,'(''  n_r_ic_max   ='',i6, &
            &        '' = number of radial grid points in IC'')') n_r_ic_max
            write(n_out,'(''  n_cheb_ic_max='',i6)') n_cheb_ic_max-1
            write(n_out,'(''  max cheb deg ='',i6)') 2*(n_cheb_ic_max-1)
         end if
         write(n_out,'(''  l_max        ='',i6, '' = max degree of Plm'')') l_max
         write(n_out,'(''  m_min        ='',i6, '' = min oder of Plm'')') m_min
         write(n_out,'(''  m_max        ='',i6, '' = max oder of Plm'')') m_max
         if ( lm_max < 1000000 ) then
            write(n_out,'(''  lm_max       ='',i6, '' = no of l/m combinations'')') lm_max
         else
            write(n_out,'(''  lm_max       ='',i8, '' = no of l/m combinations'')') lm_max
         end if
         write(n_out,'(''  minc         ='',i6, '' = longitude symmetry wave no'')') minc
         write(n_out,'(''  nalias       ='',i6, &
         &        '' = spher. harm. deal. factor '')') nalias

      end if

   end subroutine writeInfo
!------------------------------------------------------------------------------
end module preCalculations
