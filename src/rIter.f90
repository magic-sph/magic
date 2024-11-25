module rIter_mod
   !
   ! This module actually handles the loop over the radial levels. It contains
   ! the spherical harmonic transforms and the operations on the arrays in
   ! physical space.
   !

#ifdef WITHOMP
   use omp_lib
#endif
   use precision_mod
   use num_param, only: phy2lm_counter, lm2phy_counter, nl_counter, &
       &                td_counter
   use parallel_mod
   use truncation, only: n_phi_max, lm_max, lm_maxMag
   use logic, only: l_mag, l_conv, l_mag_kin, l_heat, l_ht, l_anel,  &
       &            l_mag_LF, l_conv_nl, l_mag_nl, l_b_nl_cmb,       &
       &            l_b_nl_icb, l_rot_ic, l_cond_ic, l_rot_ma,       &
       &            l_cond_ma, l_dtB, l_store_frame, l_movie_oc,     &
       &            l_TO, l_chemical_conv, l_probe, l_full_sphere,   &
       &            l_precession, l_centrifuge, l_adv_curl,          &
       &            l_double_curl, l_parallel_solve, l_single_matrix,&
       &            l_temperature_diff, l_RMS, l_phase_field,        &
       &            l_onset, l_DTrMagSpec, l_ehd_dep
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop, nRstartMag, &
       &                  nRstopMag
   use radial_functions, only: or2, orho1, l_R
   use constants, only: zero
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use torsional_oscillations, only: prep_TO_axi, getTO, getTOnext, getTOfinish
#ifdef WITH_MPI
   use graphOut_mod, only: graphOut_mpi, graphOut_mpi_header
#else
   use graphOut_mod, only: graphOut, graphOut_header
#endif
   use dtB_mod, only: get_dtBLM, get_dH_dtBLM
   use out_movie, only: store_movie_frame
   use outRot, only: get_lorentz_torque
   use courant_mod, only: courant
   use nonlinear_bcs, only: get_br_v_bcs, v_rigid_boundary, v_center_sphere
   use power, only: get_visc_heat
   use outMisc_mod, only: get_ekin_solid_liquid, get_hemi, get_helicity
   use outPar_mod, only: get_fluxes, get_nlBlayers, get_perpPar
   use geos, only: calcGeos
   use sht
   use fields, only: s_Rloc, ds_Rloc, z_Rloc, dz_Rloc, p_Rloc,    &
       &             b_Rloc, db_Rloc, ddb_Rloc, aj_Rloc,dj_Rloc,  &
       &             w_Rloc, dw_Rloc, ddw_Rloc, xi_Rloc, omega_ic,&
       &             omega_ma, phi_Rloc
   use time_schemes, only: type_tscheme
   use physical_parameters, only: ktops, kbots, n_r_LCR, ktopv, kbotv
   use rIteration, only: rIter_t
   use RMS, only: get_nl_RMS, transform_to_lm_RMS, compute_lm_forces, &
       &          transform_to_grid_RMS
   use probe_mod
   use special, only: ellip_fac_icb, l_radial_flow_bc

   implicit none

   private

   type, public, extends(rIter_t) :: rIter_single_t
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: radialLoop
      procedure :: transform_to_grid_space
      procedure :: transform_to_lm_space
   end type rIter_single_t

contains

   subroutine initialize(this)

      class(rIter_single_t) :: this

      call this%gsa%initialize()
      call this%nl_lm%initialize(lm_max)

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(rIter_single_t) :: this

      call this%gsa%finalize()
      call this%nl_lm%finalize()

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine radialLoop(this,l_graph,l_frame,time,timeStage,tscheme,dtLast,  &
              &          lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,        &
              &          lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,          &
              &          lFluxProfCalc,lPerpParCalc,lGeosCalc,lHemiCalc,      &
              &          lPhaseCalc,l_probe_out,dsdt,dwdt,dzdt,dpdt,dxidt,    &
              &          dphidt,dbdt,djdt,dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,     &
              &          lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb,    &
              &          br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,dtrkc,dthkc)
      !
      ! This subroutine handles the main loop over the radial levels. It calls
      ! the SH transforms, computes the nonlinear terms on the grid and bring back
      ! the quantities in spectral space. This is the most time consuming part
      ! of MagIC.
      !

      class(rIter_single_t) :: this

      !--- Input of variables:
      logical,             intent(in) :: l_graph,l_frame
      logical,             intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,             intent(in) :: lPowerCalc,lHemiCalc
      logical,             intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,             intent(in) :: lRmsCalc,lGeosCalc,lPhaseCalc
      logical,             intent(in) :: l_probe_out
      logical,             intent(in) :: lPressCalc
      logical,             intent(in) :: lPressNext
      real(cp),            intent(in) :: time,timeStage,dtLast
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dphidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVSrLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVXirLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in updateB.f90):
      complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma, lorentz_torque_ic

      !-- Courant citeria:
      real(cp),    intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)

      integer :: nR, nBc
      logical :: lMagNlBc, l_bound, lDeriv

      if ( l_graph ) then
#ifdef WITH_MPI
         call graphOut_mpi_header(time)
#else
         call graphOut_header(time)
#endif
      end if

      if ( rank == 0 ) then
         dtrkc(n_r_cmb)=1.e10_cp
         dthkc(n_r_cmb)=1.e10_cp
      elseif (rank == n_procs-1) then
         dtrkc(n_r_icb)=1.e10_cp
         dthkc(n_r_icb)=1.e10_cp
      end if

      !------ Set nonlinear terms that are possibly needed at the boundaries.
      !       They may be overwritten by get_td later.
      if ( rank == 0 ) then
         if ( l_mag ) dVxBhLM(:,n_r_cmb)=zero
         if ( l_double_curl ) dVxVhLM(:,n_r_cmb)=zero
      else if (rank == n_procs-1) then
         if ( l_mag ) dVxBhLM(:,n_r_icb)=zero
         if ( l_double_curl ) dVxVhLM(:,n_r_icb)=zero
      end if

      lorentz_torque_ma = 0.0_cp
      lorentz_torque_ic = 0.0_cp

      !------ Having to calculate non-linear boundary terms?
      lMagNlBc=.false.
      if ( ( l_mag_nl .or. l_mag_kin ) .and.                          &
           &       ( ktopv == 1 .or. l_cond_ma .or.                   &
           &          ( ktopv == 2 .and. l_rot_ma ) ) .or.            &
           &       ( kbotv == 1 .or. l_cond_ic .or.                   &
           &          ( kbotv == 2 .and. l_rot_ic ) ) )               &
           &     lMagNlBc=.true.


      do nR=nRstart,nRstop
         l_Bound = ( nR == n_r_icb ) .or. ( nR == n_r_cmb )

         nBc = 0
         lDeriv = .true.
         if ( nR == n_r_cmb ) then
            nBc = ktopv
            lDeriv= lTOCalc .or. lHelCalc .or. l_frame .or. lPerpParCalc   &
            &       .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
            &       lPowerCalc .or. lGeosCalc .or. lHemiCalc
         else if ( nR == n_r_icb ) then
            nBc = kbotv
            lDeriv= lTOCalc .or. lHelCalc .or. l_frame  .or. lPerpParCalc  &
            &       .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
            &       lPowerCalc .or. lGeosCalc .or. lHemiCalc
         end if

         if ( l_parallel_solve .or. (l_single_matrix .and. l_temperature_diff) ) then
            ! We will need the nonlinear terms on ricb for the pressure l=m=0
            ! equation
            lDeriv=.true.
            nBc=0
            l_Bound=.false.
         end if

         if ( lRmsCalc ) nBc=0 ! One also needs to compute the boundaries in that case

         dtrkc(nR)=1e10_cp
         dthkc(nR)=1e10_cp

         if ( lTOnext .or. lTOnext2 .or. lTOCalc ) then
            call prep_TO_axi(z_Rloc(:,nR), dz_Rloc(:,nR))
         end if

         call this%nl_lm%set_zero()

         if ( .not. l_onset ) then
            call lm2phy_counter%start_count()
            call this%transform_to_grid_space(nR, nBc, lViscBcCalc, lRmsCalc,       &
                 &                            lPressCalc, lTOCalc, lPowerCalc,      &
                 &                            lFluxProfCalc, lPerpParCalc, lHelCalc,&
                 &                            lGeosCalc, lHemiCalc, l_frame, lDeriv)
            call lm2phy_counter%stop_count(l_increment=.false.)
         end if

         !--------- Calculation of nonlinear products in grid space:
         if ( (.not. l_onset ) .and. ((.not. l_bound) .or. lMagNlBc .or. lRmsCalc) ) then

            call nl_counter%start_count()
            call this%gsa%get_nl(timeStage, nR, nBc, lRmsCalc)
            call nl_counter%stop_count(l_increment=.false.)

            !-- Get nl loop for r.m.s. computation
            if ( l_RMS ) then
               call get_nl_RMS(nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%dvrdrc,&
                    &          this%gsa%dvrdtc,this%gsa%dvrdpc,this%gsa%dvtdrc,          &
                    &          this%gsa%dvtdpc,this%gsa%dvpdrc,this%gsa%dvpdpc,          &
                    &          this%gsa%cvrc,this%gsa%Advt,this%gsa%Advp,this%gsa%LFt,   &
                    &          this%gsa%LFp,tscheme,lRmsCalc)
            end if

            call phy2lm_counter%start_count()
            call this%transform_to_lm_space(nR, lRmsCalc, dVSrLM(:,nR), dVXirLM(:,nR), &
                 &                         dphidt(:,nR))
            call phy2lm_counter%stop_count(l_increment=.false.)
         else if ( l_mag ) then
            this%nl_lm%VxBtLM(:)=zero
            this%nl_lm%VxBpLM(:)=zero
         end if

         !---- Calculation of nonlinear products needed for conducting mantle or
         !     conducting inner core if free stress BCs are applied:
         !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
         !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
         !     and br_vp_lm_cmb in lm-space, respectively the contribution
         !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
         !     These products are used in get_b_nl_bcs.
         if ( nR == n_r_cmb .and. l_b_nl_cmb ) then
            br_vt_lm_cmb(:)=zero
            br_vp_lm_cmb(:)=zero
            call get_br_v_bcs(this%gsa%brc, this%gsa%vtc, this%gsa%vpc,omega_ma,  &
                 &            or2(nR), orho1(nR), br_vt_lm_cmb, br_vp_lm_cmb)
         else if ( nR == n_r_icb .and. l_b_nl_icb ) then
            br_vt_lm_icb(:)=zero
            br_vp_lm_icb(:)=zero
            call get_br_v_bcs(this%gsa%brc, this%gsa%vtc, this%gsa%vpc, omega_ic,  &
                 &            or2(nR), orho1(nR), br_vt_lm_icb, br_vp_lm_icb)
         end if
         !--------- Calculate Lorentz torque on inner core:
         !          each call adds the contribution of the theta-block to
         !          lorentz_torque_ic
         if ( nR == n_r_icb .and. l_mag_LF .and. l_rot_ic .and. l_cond_ic  ) then
            call get_lorentz_torque(lorentz_torque_ic, this%gsa%brc,  &
                 &                  this%gsa%bpc, nR)
         end if

         !--------- Calculate Lorentz torque on mantle:
         !          note: this calculates a torque of a wrong sign.
         !          sign is reversed at the end of the theta blocking.
         if ( nR == n_r_cmb .and. l_mag_LF .and. l_rot_ma .and. l_cond_ma ) then
            call get_lorentz_torque(lorentz_torque_ma, this%gsa%brc, &
                 &                  this%gsa%bpc, nR)
         end if

         !--------- Calculate courant condition parameters:
         if ( (.not. l_full_sphere .or. nR /= n_r_icb) .and. (.not. l_onset) ) then
            call courant(nR, dtrkc(nR), dthkc(nR), this%gsa%vrc,              &
                 &       this%gsa%vtc,this%gsa%vpc,this%gsa%brc,this%gsa%btc, &
                 &       this%gsa%bpc, tscheme%courfac, tscheme%alffac)
         end if

         !--------- Since the fields are given at gridpoints here, this is a good
         !          point for graphical output:
         if ( l_graph ) then
#ifdef WITH_MPI
            call graphOut_mpi(nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc, &
                 &            this%gsa%brc,this%gsa%btc,this%gsa%bpc,    &
                 &            this%gsa%sc,this%gsa%pc,this%gsa%xic,      &
                 &            this%gsa%phic)
#else
            call graphOut(nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,     &
                 &        this%gsa%brc,this%gsa%btc,this%gsa%bpc,        &
                 &        this%gsa%sc,this%gsa%pc,this%gsa%xic,          &
                 &        this%gsa%phic)
#endif
         end if

         if ( l_probe_out ) then
            call probe_out(time, nR, this%gsa%vpc, this%gsa%brc, this%gsa%btc)
         end if

         !--------- Helicity output:
         if ( lHelCalc ) then
            call get_helicity(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,         &
                 &            this%gsa%cvrc,this%gsa%dvrdtc,this%gsa%dvrdpc,  &
                 &            this%gsa%dvtdrc,this%gsa%dvpdrc,nR)
         end if

         !-- North/South hemisphere differences
         if ( lHemiCalc ) then
            call get_hemi(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,nR,'V')
            if ( l_mag ) call get_hemi(this%gsa%brc,this%gsa%btc,this%gsa%bpc,nR,'B')
         end if

         !-- Viscous heating:
         if ( lPowerCalc ) then
            call get_visc_heat(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,          &
                 &             this%gsa%cvrc,this%gsa%dvrdrc,this%gsa%dvrdtc,   &
                 &             this%gsa%dvrdpc,this%gsa%dvtdrc,this%gsa%dvtdpc, &
                 &             this%gsa%dvpdrc,this%gsa%dvpdpc,nR)
         end if

         !-- horizontal velocity :
         if ( lViscBcCalc ) then
            call get_nlBLayers(this%gsa%vtc,this%gsa%vpc,this%gsa%dvtdrc,    &
                 &             this%gsa%dvpdrc,this%gsa%drSc,this%gsa%dsdtc, &
                 &             this%gsa%dsdpc,nR )
         end if

         !-- Radial flux profiles
         if ( lFluxProfCalc ) then
            call get_fluxes(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,            &
                 &          this%gsa%dvrdrc,this%gsa%dvtdrc,this%gsa%dvpdrc,   &
                 &          this%gsa%dvrdtc,this%gsa%dvrdpc,this%gsa%sc,       &
                 &          this%gsa%pc,this%gsa%brc,this%gsa%btc,this%gsa%bpc,&
                 &          this%gsa%cbtc,this%gsa%cbpc,nR)
         end if

         !-- Kinetic energy in the solid and liquid phases
         if ( lPhaseCalc ) then
            call get_ekin_solid_liquid(this%gsa%sc,this%gsa%drsc,this%gsa%vrc, &
                 &                     this%gsa%vtc,this%gsa%vpc,this%gsa%phic,nR)
         end if

         !-- Kinetic energy parallel and perpendicular to rotation axis
         if ( lPerpParCalc ) then
            call get_perpPar(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,nR)
         end if

         !-- Geostrophic/non-geostrophic flow components
         if ( lGeosCalc ) then
            call calcGeos(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc, &
                 &        this%gsa%dvrdpc,this%gsa%dvpdrc,nR)
         end if

         !--------- Movie output:
         if ( l_frame .and. l_movie_oc .and. l_store_frame ) then
            call store_movie_frame(nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,      &
                 &                 this%gsa%brc,this%gsa%btc,this%gsa%bpc,         &
                 &                 this%gsa%sc,this%gsa%drSc,this%gsa%xic,         &
                 &                 this%gsa%phic,this%gsa%dvrdpc,                  &
                 &                 this%gsa%dvpdrc,this%gsa%dvtdrc,this%gsa%dvrdtc,&
                 &                 this%gsa%cvrc,this%gsa%cbrc,this%gsa%cbtc)
         end if

         !--------- Stuff for special output:
         !--------- Calculation of magnetic field production and advection terms
         !          for graphic output:
         if ( l_dtB ) then
            call get_dtBLM(nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,  &
                 &         this%gsa%brc,this%gsa%btc,this%gsa%bpc)
         end if


         !--------- Torsional oscillation terms:
         if ( lTONext .or. lTONext2 ) then
            call getTOnext(this%gsa%brc,this%gsa%btc,this%gsa%bpc,lTONext, &
                 &         lTONext2,tscheme%dt(1),dtLast,nR)
         end if

         if ( lTOCalc ) then
            call getTO(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc,   &
                 &     this%gsa%dvpdrc,this%gsa%brc,this%gsa%btc,this%gsa%bpc, &
                 &     this%gsa%cbrc,this%gsa%cbtc,this%gsa%phic,dtLast,nR)
         end if

         !-- Partial calculation of time derivatives (horizontal parts):
         !   input flm...  is in (l,m) space at radial grid points nR !
         !   Only dVxBh needed for boundaries !
         !   get_td finally calculates the d*dt terms needed for the
         !   time step performed in LMLoop.f90 .
         call td_counter%start_count()
         if ( l_conv ) then
            call this%nl_lm%get_dzdt(nR, nBc, dzdt(:,nR))
            if ( l_double_curl ) then
               call this%nl_lm%get_dwdt_double_curl(nR, nBc, dwdt(:,nR), dVxVhLM(:,nR))
            else
               call this%nl_lm%get_dwdt(nR, nBc, dwdt(:,nR))
            end if
         end if
         if ( (.not. l_double_curl) .or. lPressNext ) then
            call this%nl_lm%get_dpdt(nR, nBc, dpdt(:,nR))
         end if
         if ( l_heat ) call this%nl_lm%get_dsdt(nR, nBc, dsdt(:,nR), dVSrLM(:,nR))
         if ( l_chemical_conv ) then
            call this%nl_lm%get_dxidt(nBc, dxidt(:,nR), dVXirLM(:,nR))
         end if
         if ( l_mag ) then
            call this%nl_lm%get_dbdt(nR, nBc, dbdt(:,nR), djdt(:,nR), dVxBhLM(:,nR))
         end if
         call td_counter%stop_count(l_increment=.false.)

         !-- Finish computation of r.m.s. forces
         if ( lRmsCalc ) then
            call compute_lm_forces(nR, this%nl_lm%AdvrLM)
         end if

         !-- Finish calculation of TO variables:
         if ( lTOcalc ) call getTOfinish(nR, dtLast)

         !--- Form partial horizontal derivaties of magnetic production and
         !    advection terms:
         if ( l_dtB ) call get_dH_dtBLM(nR)

      end do

      !------ Set nonlinear terms that are possibly needed at the boundaries.
      !       They may be overwritten by get_td later.
      if ( rank == 0 ) then
         if ( l_heat ) dVSrLM(:,n_r_cmb) =zero
         if ( l_chemical_conv ) dVXirLM(:,n_r_cmb)=zero
      else if (rank == n_procs-1) then
         if ( l_heat ) dVSrLM(:,n_r_icb) =zero
         if ( l_chemical_conv ) dVXirLM(:,n_r_icb)=zero
      end if

      phy2lm_counter%n_counts=phy2lm_counter%n_counts+1
      lm2phy_counter%n_counts=lm2phy_counter%n_counts+1
      nl_counter%n_counts=nl_counter%n_counts+1
      td_counter%n_counts=td_counter%n_counts+1

      !----- Correct sign of mantle Lorentz torque (see above):
      lorentz_torque_ma=-lorentz_torque_ma

   end subroutine radialLoop
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space(this, nR, nBc, lViscBcCalc, lRmsCalc,      &
              &                       lPressCalc, lTOCalc, lPowerCalc,           &
              &                       lFluxProfCalc, lPerpParCalc, lHelCalc,     &
              &                       lGeosCalc, lHemiCalc, l_frame, lDeriv)
      !
      ! This subroutine actually handles the spherical harmonic transforms from
      ! (\ell,m) space to (\theta,\phi) space.
      !

      class(rIter_single_t) :: this

      !--Input variables
      integer, intent(in) :: nR, nBc
      logical, intent(in) :: lViscBcCalc, lRmsCalc, lPressCalc, lTOCalc, lPowerCalc
      logical, intent(in) :: lFluxProfCalc, lPerpParCalc, lHelCalc, l_frame
      logical, intent(in) :: lDeriv, lGeosCalc, lHemiCalc

      if ( l_conv .or. l_mag_kin ) then
         if ( l_heat ) then
            call scal_to_spat(s_Rloc(:,nR), this%gsa%sc, l_R(nR))
            if ( lViscBcCalc ) then
               call scal_to_grad_spat(s_Rloc(:,nR), this%gsa%dsdtc, this%gsa%dsdpc, &
                    &                 l_R(nR))
               if ( nR == n_r_cmb .and. ktops==1) then
                  this%gsa%dsdtc(:,:)=0.0_cp
                  this%gsa%dsdpc(:,:)=0.0_cp
               end if
               if ( nR == n_r_icb .and. kbots==1) then
                  this%gsa%dsdtc(:,:)=0.0_cp
                  this%gsa%dsdpc(:,:)=0.0_cp
               end if
            end if
         end if

         if ( lRmsCalc ) call transform_to_grid_RMS(nR, p_Rloc)

         !-- Pressure
         if ( lPressCalc ) call scal_to_spat(p_Rloc(:,nR), this%gsa%pc, l_R(nR))

         !-- Composition
         if ( l_chemical_conv ) call scal_to_spat(xi_Rloc(:,nR), this%gsa%xic, l_R(nR))

         !-- Phase field
         if ( l_phase_field ) call scal_to_spat(phi_Rloc(:,nR), this%gsa%phic, l_R(nR))

         if ( l_HT .or. lViscBcCalc ) then
            call scal_to_spat(ds_Rloc(:,nR), this%gsa%drsc, l_R(nR))
         endif
         if ( nBc == 0 ) then ! Bulk points
            !-- pol, sph, tor > ur,ut,up
            call torpol_to_spat(w_Rloc(:,nR), dw_Rloc(:,nR),  z_Rloc(:,nR), &
                 &              this%gsa%vrc, this%gsa%vtc, this%gsa%vpc, l_R(nR))

            !-- Advection is treated as u \times \curl u
            if ( l_adv_curl ) then
               !-- z,dz,w,dd< -> wr,wt,wp
               call torpol_to_curl_spat(or2(nR), w_Rloc(:,nR), ddw_Rloc(:,nR), &
                    &                   z_Rloc(:,nR), dz_Rloc(:,nR),           &
                    &                   this%gsa%cvrc, this%gsa%cvtc,          &
                    &                   this%gsa%cvpc, l_R(nR))

               !-- For some outputs one still need the other terms
               if ( lViscBcCalc .or. lPowerCalc .or. lRmsCalc .or. lFluxProfCalc &
               &    .or. lTOCalc .or. lHelCalc .or. lPerpParCalc .or. lGeosCalc  &
               &    .or. lHemiCalc                                               &
               &    .or. ( l_frame .and. l_movie_oc .and. l_store_frame) ) then
                  call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR),         &
                       &              dz_Rloc(:,nR), this%gsa%dvrdrc,        &
                       &              this%gsa%dvtdrc, this%gsa%dvpdrc, l_R(nR))
                  call pol_to_grad_spat(w_Rloc(:,nR), this%gsa%dvrdtc, &
                       &                this%gsa%dvrdpc, l_R(nR))
                  call torpol_to_dphspat(dw_Rloc(:,nR),  z_Rloc(:,nR), &
                       &                 this%gsa%dvtdpc, this%gsa%dvpdpc, l_R(nR))
               end if

            else ! Advection is treated as u\grad u

               call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR), dz_Rloc(:,nR), &
                 &                 this%gsa%dvrdrc, this%gsa%dvtdrc,             &
                 &                 this%gsa%dvpdrc, l_R(nR))

               call pol_to_curlr_spat(z_Rloc(:,nR), this%gsa%cvrc, l_R(nR))

               call pol_to_grad_spat(w_Rloc(:,nR), this%gsa%dvrdtc, this%gsa%dvrdpc,&
                    &                l_R(nR))
               call torpol_to_dphspat(dw_Rloc(:,nR),  z_Rloc(:,nR), &
                    &                 this%gsa%dvtdpc, this%gsa%dvpdpc, l_R(nR))
            end if

         else if ( nBc == 1 ) then ! Stress free
             ! TODO don't compute vrc as it is set to 0 afterward
            call torpol_to_spat(w_Rloc(:,nR), dw_Rloc(:,nR),  z_Rloc(:,nR), &
                 &              this%gsa%vrc, this%gsa%vtc, this%gsa%vpc, l_R(nR))
            this%gsa%vrc(:,:)=0.0_cp
            if ( lDeriv ) then
               this%gsa%dvrdtc(:,:)=0.0_cp
               this%gsa%dvrdpc(:,:)=0.0_cp
               call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR), dz_Rloc(:,nR), &
                    &              this%gsa%dvrdrc, this%gsa%dvtdrc,             &
                    &              this%gsa%dvpdrc, l_R(nR))
               call pol_to_curlr_spat(z_Rloc(:,nR), this%gsa%cvrc, l_R(nR))
               call torpol_to_dphspat(dw_Rloc(:,nR),  z_Rloc(:,nR), &
                    &                 this%gsa%dvtdpc, this%gsa%dvpdpc, l_R(nR))
            end if
         else if ( nBc == 2 ) then
            if ( nR == n_r_cmb .and. (.not. l_radial_flow_bc) ) then
               call v_rigid_boundary(nR, omega_ma, lDeriv, this%gsa%vrc,      &
                    &              this%gsa%vtc, this%gsa%vpc, this%gsa%cvrc, &
                    &              this%gsa%dvrdtc, this%gsa%dvrdpc,          &
                    &              this%gsa%dvtdpc,this%gsa%dvpdpc)
            else if ( nR == n_r_icb .and. ellip_fac_icb == 0.0_cp ) then
               call v_rigid_boundary(nR, omega_ic, lDeriv, this%gsa%vrc,    &
                    &              this%gsa%vtc, this%gsa%vpc,              &
                    &              this%gsa%cvrc, this%gsa%dvrdtc,          &
                    &              this%gsa%dvrdpc, this%gsa%dvtdpc,        &
                    &              this%gsa%dvpdpc)
            end if
            if ( lDeriv ) then
               call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR), dz_Rloc(:,nR), &
                    &              this%gsa%dvrdrc, this%gsa%dvtdrc,             &
                    &              this%gsa%dvpdrc, l_R(nR))
            end if
         end if

         if ( nR == n_r_icb .and. l_full_sphere ) then
            call v_center_sphere(ddw_Rloc(:,nR), this%gsa%vrc, this%gsa%vtc, &
                 &               this%gsa%vpc)
         end if
      end if

      if ( l_mag .or. l_mag_LF ) then
         call torpol_to_spat(b_Rloc(:,nR), db_Rloc(:,nR),  aj_Rloc(:,nR),    &
              &              this%gsa%brc, this%gsa%btc, this%gsa%bpc, l_R(nR))

         if ( lDeriv ) then
            call torpol_to_curl_spat(or2(nR), b_Rloc(:,nR), ddb_Rloc(:,nR), &
                 &                   aj_Rloc(:,nR), dj_Rloc(:,nR),          &
                 &                   this%gsa%cbrc, this%gsa%cbtc,          &
                 &                   this%gsa%cbpc, l_R(nR))
         end if

         if ( nR == n_r_icb .and. l_full_sphere ) then
            call v_center_sphere(ddb_Rloc(:,nR), this%gsa%brc, this%gsa%btc, &
                 &               this%gsa%bpc)
         end if
      end if

   end subroutine transform_to_grid_space
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space(this, nR, lRmsCalc, dVSrLM, dVXirLM, dphidt)
      !
      ! This subroutine actually handles the spherical harmonic transforms from
      ! (\theta,\phi) space to (\ell,m) space.
      !

      class(rIter_single_t) :: this

      !-- Input variables
      integer, intent(in) :: nR
      logical, intent(in) :: lRmsCalc

      !-- Output variables
      complex(cp), intent(out) :: dVSrLM(lm_max)
      complex(cp), intent(out) :: dVXirLM(lm_max)
      complex(cp), intent(out) :: dphidt(lm_max)

      !-- Local variables
      integer :: nPhi, nPhStart, nPhStop

      if ( l_conv_nl .or. l_mag_LF ) then

         !$omp parallel default(shared) private(nPhStart,nPhStop,nPhi)
         nPhStart=1; nPhStop=n_phi_max
         call get_openmp_blocks(nPhStart,nPhStop)

         do nPhi=nPhStart, nPhStop
            if ( l_conv_nl .and. l_mag_LF ) then
               if ( nR>n_r_LCR ) then
                  this%gsa%Advr(:,nPhi)=this%gsa%Advr(:,nPhi) + this%gsa%LFr(:,nPhi)
                  this%gsa%Advt(:,nPhi)=this%gsa%Advt(:,nPhi) + this%gsa%LFt(:,nPhi)
                  this%gsa%Advp(:,nPhi)=this%gsa%Advp(:,nPhi) + this%gsa%LFp(:,nPhi)
               end if
            else if ( l_mag_LF ) then
               if ( nR > n_r_LCR ) then
                  this%gsa%Advr(:,nPhi) = this%gsa%LFr(:,nPhi)
                  this%gsa%Advt(:,nPhi) = this%gsa%LFt(:,nPhi)
                  this%gsa%Advp(:,nPhi) = this%gsa%LFp(:,nPhi)
               else
                  this%gsa%Advr(:,nPhi)=0.0_cp
                  this%gsa%Advt(:,nPhi)=0.0_cp
                  this%gsa%Advp(:,nPhi)=0.0_cp
               end if
            end if

            if ( l_precession ) then
               this%gsa%Advr(:,nPhi)=this%gsa%Advr(:,nPhi) + this%gsa%PCr(:,nPhi)
               this%gsa%Advt(:,nPhi)=this%gsa%Advt(:,nPhi) + this%gsa%PCt(:,nPhi)
               this%gsa%Advp(:,nPhi)=this%gsa%Advp(:,nPhi) + this%gsa%PCp(:,nPhi)
            end if

            if ( l_centrifuge ) then
               this%gsa%Advr(:, nPhi)=this%gsa%Advr(:,nPhi) + this%gsa%CAr(:,nPhi)
               this%gsa%Advt(:, nPhi)=this%gsa%Advt(:,nPhi) + this%gsa%CAt(:,nPhi)
            end if

            if ( l_ehd_dep ) then
               this%gsa%Advr(:, nPhi)=this%gsa%Advr(:,nPhi) + this%gsa%DEPFr(:,nPhi)
            end if   
         end do
         !$omp end parallel

         call spat_to_qst(this%gsa%Advr, this%gsa%Advt, this%gsa%Advp, &
              &           this%nl_lm%AdvrLM, this%nl_lm%AdvtLM,        &
              &           this%nl_lm%AdvpLM, l_R(nR))
      end if

      if ( l_heat ) then
         call spat_to_qst(this%gsa%VSr, this%gsa%VSt, this%gsa%VSp, &
              &           dVSrLM, this%nl_lm%VStLM, this%nl_lm%VSpLM, l_R(nR))

         if ( l_anel ) call scal_to_SH(this%gsa%heatTerms, this%nl_lm%heatTermsLM, &
                            &          l_R(nR))
      end if
      if ( l_chemical_conv ) then
         call spat_to_qst(this%gsa%VXir, this%gsa%VXit, this%gsa%VXip, &
              &           dVXirLM, this%nl_lm%VXitLM, this%nl_lm%VXipLM, l_R(nR))
      end if
      if ( l_phase_field ) call scal_to_SH(this%gsa%phiTerms, dphidt,l_R(nR))
      if ( l_mag_nl ) then
         if ( nR>n_r_LCR ) then
            call spat_to_qst(this%gsa%VxBr, this%gsa%VxBt, this%gsa%VxBp, &
                 &           this%nl_lm%VxBrLM, this%nl_lm%VxBtLM,        &
                 &           this%nl_lm%VxBpLM, l_R(nR))
         else
            call spat_to_sphertor(this%gsa%VxBt, this%gsa%VxBp, this%nl_lm%VxBtLM, &
                 &                this%nl_lm%VxBpLM, l_R(nR))
         end if
      end if

      if ( lRmsCalc ) call transform_to_lm_RMS(nR, this%gsa%LFr)

   end subroutine transform_to_lm_space
!-------------------------------------------------------------------------------
end module rIter_mod
