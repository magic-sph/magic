#include "perflib_preproc.cpp"
module rIterThetaBlocking_shtns_mod
#ifdef WITHOMP
   use omp_lib
#endif
   use precision_mod
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t
   use num_param, only: phy2lm_counter, lm2phy_counter, nl_counter,  &
       &                td_counter
   use parallel_mod, only: get_openmp_blocks
   use truncation, only: lmP_max, n_theta_max, n_phi_max, n_r_cmb,   &
       &                 n_r_icb, n_lmP_loc, n_theta_loc,            &
       &                 nThetaStart, nThetaStop
   use logic, only: l_mag, l_conv, l_mag_kin, l_heat, l_ht, l_anel,  &
       &            l_mag_LF, l_conv_nl, l_mag_nl, l_b_nl_cmb,       &
       &            l_b_nl_icb, l_rot_ic, l_cond_ic, l_rot_ma,       &
       &            l_cond_ma, l_dtB, l_store_frame, l_movie_oc,     &
       &            l_TO, l_chemical_conv, l_probe,                  &
       &            l_precession, l_centrifuge, l_adv_curl,          &
       &            l_full_sphere
   use radial_functions, only: or2, orho1, l_R
   use constants, only: zero
   use leg_helper_mod, only: leg_helper_t
   use nonlinear_lm_mod, only:nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use TO_arrays_mod, only: TO_arrays_t
   use dtB_arrays_mod, only: dtB_arrays_t
   use torsional_oscillations, only: getTO, getTOnext, getTOfinish
#ifdef WITH_MPI
   use graphOut_mod, only: graphOut_mpi
#else
   use graphOut_mod, only: graphOut
#endif
   use dtB_mod, only: get_dtBLM, get_dH_dtBLM
   use out_movie, only: store_movie_frame
   use outRot, only: get_lorentz_torque
   use courant_mod, only: courant
   use nonlinear_bcs, only: get_br_v_bcs, v_rigid_boundary
   use nl_special_calc
   use shtns
   use horizontal_data
   use fields, only: s_Rloc, ds_Rloc, z_Rloc, dz_Rloc, p_Rloc,   &
       &             b_Rloc, db_Rloc, ddb_Rloc, aj_Rloc,dj_Rloc, &
       &             w_Rloc, dw_Rloc, ddw_Rloc, xi_Rloc
   use time_schemes, only: type_tscheme
   use physical_parameters, only: ktops, kbots, n_r_LCR
   use probe_mod
   use communications  ! DELETEME

   implicit none

   private

   type, public, extends(rIterThetaBlocking_t) :: rIterThetaBlocking_shtns_t
      integer :: nThreads
      type(grid_space_arrays_t) :: gsa
      type(grid_space_arrays_t) :: gsa_dist
      type(TO_arrays_t) :: TO_arrays
      type(dtB_arrays_t) :: dtB_arrays
      type(nonlinear_lm_t) :: nl_lm
      type(nonlinear_lm_t) :: nl_lm_dist
      real(cp) :: lorentz_torque_ic,lorentz_torque_ma
   contains
      procedure :: initialize => initialize_rIterThetaBlocking_shtns
      procedure :: finalize => finalize_rIterThetaBlocking_shtns
      procedure :: do_iteration => do_iteration_ThetaBlocking_shtns
      procedure :: getType => getThisType
      procedure :: transform_to_grid_space_shtns => transform_to_grid_space_shtns
      procedure :: transform_to_lm_space_shtns => transform_to_lm_space_shtns
   end type rIterThetaBlocking_shtns_t

contains

   function getThisType(this)

      class(rIterThetaBlocking_shtns_t) :: this
      character(len=100) :: getThisType
      getThisType="rIterThetaBlocking_shtns_t"

   end function getThisType
!------------------------------------------------------------------------------
   subroutine initialize_rIterThetaBlocking_shtns(this)

      class(rIterThetaBlocking_shtns_t) :: this

      call this%allocate_common_arrays()
      call this%gsa%initialize()
      call this%gsa_dist%initialize_dist()
      if ( l_TO ) call this%TO_arrays%initialize()
      call this%dtB_arrays%initialize()
      call this%nl_lm%initialize(lmP_max)
      call this%nl_lm_dist%initialize(n_lmP_loc)

   end subroutine initialize_rIterThetaBlocking_shtns
!------------------------------------------------------------------------------
   subroutine finalize_rIterThetaBlocking_shtns(this)

      class(rIterThetaBlocking_shtns_t) :: this

      call this%deallocate_common_arrays()
      call this%gsa%finalize()
      call this%gsa_dist%finalize()
      if ( l_TO ) call this%TO_arrays%finalize()
      call this%dtB_arrays%finalize()
      call this%nl_lm%finalize()
      call this%nl_lm_dist%finalize()

   end subroutine finalize_rIterThetaBlocking_shtns
!------------------------------------------------------------------------------
   subroutine do_iteration_ThetaBlocking_shtns(this,nR,nBc,time,timeStage,   &
              &           tscheme,dtLast,dsdt,dwdt,dzdt,dpdt,dxidt,dbdt,djdt,&
              &           dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,                    &
              &           br_vt_lm_cmb,br_vp_lm_cmb,                         &
              &           br_vt_lm_icb,br_vp_lm_icb,                         &
              &           lorentz_torque_ic, lorentz_torque_ma,              &
              &           HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,viscLMr,         &
              &           uhLMr,duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,   &
              &           fpoynLMr,fresLMr,EperpLMr,EparLMr,EperpaxiLMr,     &
              &           EparaxiLMr)

      class(rIterThetaBlocking_shtns_t) :: this
      integer,             intent(in) :: nR,nBc
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time,dtLast,timeStage

      complex(cp), intent(out) :: dwdt(:),dzdt(:),dpdt(:),dsdt(:),dVSrLM(:)
      complex(cp), intent(out) :: dxidt(:),dVXirLM(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:),dVxVhLM(:),dVxBhLM(:)
      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma, lorentz_torque_ic
      real(cp),    intent(out) :: HelLMr(:),Hel2LMr(:),HelnaLMr(:),Helna2LMr(:)
      real(cp),    intent(out) :: viscLMr(:)
      real(cp),    intent(out) :: uhLMr(:), duhLMr(:) ,gradsLMr(:)
      real(cp),    intent(out) :: fconvLMr(:), fkinLMr(:), fviscLMr(:)
      real(cp),    intent(out) :: fpoynLMr(:), fresLMr(:)
      real(cp),    intent(out) :: EperpLMr(:), EparLMr(:), EperpaxiLMr(:), EparaxiLMr(:)

      integer :: lm
      logical :: lGraphHeader=.false.
      logical :: DEBUG_OUTPUT=.false.
      real(cp) :: lorentz_torques_ic
      
      
      ! TODO DELETEMEEEEEEEEEEEE
      complex(cp) :: br_vt_lm_cmb_dist(n_lmP_loc) ! product br*vt at CMB
      complex(cp) :: br_vp_lm_cmb_dist(n_lmP_loc) ! product br*vp at CMB
      complex(cp) :: br_vt_lm_icb_dist(n_lmP_loc) ! product br*vt at ICB
      complex(cp) :: br_vp_lm_icb_dist(n_lmP_loc) ! product br*vp at ICB

      this%nR=nR
      this%nBc=nBc
      this%isRadialBoundaryPoint=(nR == n_r_cmb).or.(nR == n_r_icb)

      this%dtrkc=1.e10_cp
      this%dthkc=1.e10_cp

      if ( this%lTOCalc ) then
         !------ Zero lm coeffs for first theta block:
         call this%TO_arrays%set_zero()
      end if

      call this%leg_helper%legPrepG(this%nR,this%nBc,this%lDeriv,this%lRmsCalc, &
           &                        this%l_frame,this%lTOnext,this%lTOnext2,    &
           &                        this%lTOcalc)

      if (DEBUG_OUTPUT) then
         write(*,"(I3,A,I1,2(A,L1))") this%nR,": nBc = ", &
              & this%nBc,", lDeriv = ",this%lDeriv,", l_mag = ",l_mag
      end if

      this%lorentz_torque_ma = 0.0_cp
      this%lorentz_torque_ic = 0.0_cp
      lorentz_torques_ic = 0.0_cp

      call this%nl_lm%set_zero()

      call lm2phy_counter%start_count()
      call this%transform_to_grid_space_shtns(this%gsa)
      call lm2phy_counter%stop_count(l_increment=.false.)

      call this%gsa%slice_all(this%gsa_dist)
      call this%nl_lm%slice_all(this%nl_lm_dist)
      
      call slice_FlmP_cmplx(br_vt_lm_cmb, br_vt_lm_cmb_dist)
      call slice_FlmP_cmplx(br_vp_lm_cmb, br_vp_lm_cmb_dist)
      call slice_FlmP_cmplx(br_vt_lm_icb, br_vt_lm_icb_dist)
      call slice_FlmP_cmplx(br_vp_lm_icb, br_vp_lm_icb_dist)
      
      !--------- Calculation of nonlinear products in grid space:
      if ( (.not.this%isRadialBoundaryPoint) .or. this%lMagNlBc .or. &
      &     this%lRmsCalc ) then

         call nl_counter%start_count()
         PERFON('get_nl')
!          call this%gsa%get_nl_shtns(timeStage, tscheme, this%nR, this%nBc, &
!               &                     this%lRmsCalc, 1, n_theta_max)
         call this%gsa_dist%get_nl_shtns(timeStage, tscheme, this%nR, this%nBc, &
              &                     this%lRmsCalc)
         PERFOFF
         call nl_counter%stop_count(l_increment=.false.)

         
         call phy2lm_counter%start_count()
         call this%transform_to_lm_space_shtns(this%gsa_dist, this%nl_lm_dist)
         call phy2lm_counter%stop_count(l_increment=.false.)

      else if ( l_mag ) then
         do lm=1,n_lmP_loc
            this%nl_lm_dist%VxBtLM(lm)=0.0_cp
            this%nl_lm_dist%VxBpLM(lm)=0.0_cp
         end do
      end if
      
      !---- Calculation of nonlinear products needed for conducting mantle or
      !     conducting inner core if free stress BCs are applied:
      !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
      !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
      !     and br_vp_lm_cmb in lm-space, respectively the contribution
      !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
      !     These products are used in get_b_nl_bcs.
      if ( this%nR == n_r_cmb .and. l_b_nl_cmb ) then
         br_vt_lm_cmb_dist(:)=zero
         br_vp_lm_cmb_dist(:)=zero
         call get_br_v_bcs(this%gsa_dist%brc,this%gsa_dist%vtc,               &
              &            this%gsa_dist%vpc,this%leg_helper%omegaMA,    &
              &            or2(this%nR),orho1(this%nR),             &
              &            br_vt_lm_cmb_dist,br_vp_lm_cmb_dist)
              
         !TODO DELETE MEEEEEEEEEEEEEEEEEEE
         call gather_FlmP(br_vt_lm_cmb_dist, br_vt_lm_cmb)
         call gather_FlmP(br_vp_lm_cmb_dist, br_vp_lm_cmb)
      else if ( this%nR == n_r_icb .and. l_b_nl_icb ) then
         br_vt_lm_icb_dist(:)=zero
         br_vp_lm_icb_dist(:)=zero
         call get_br_v_bcs(this%gsa_dist%brc,this%gsa_dist%vtc,               &
              &            this%gsa_dist%vpc,this%leg_helper%omegaIC,    &
              &            or2(this%nR),orho1(this%nR),             &
              &            br_vt_lm_icb_dist,br_vp_lm_icb_dist)
         !TODO DELETE MEEEEEEEEEEEEEEEEEEE
         call gather_FlmP(br_vt_lm_icb_dist, br_vt_lm_icb)
         call gather_FlmP(br_vp_lm_icb_dist, br_vp_lm_icb)
      end if
      
      !--------- Calculate Lorentz torque on inner core:
      !          each call adds the contribution of the theta-block to
      !          lorentz_torque_ic
      if ( this%nR == n_r_icb .and. l_mag_LF .and. l_rot_ic .and. l_cond_ic  ) then
         lorentz_torques_ic=0.0_cp
         call get_lorentz_torque(lorentz_torques_ic, this%gsa_dist%brc,this%gsa_dist%bpc,this%nR)
      end if

      !--------- Calculate Lorentz torque on mantle:
      !          note: this calculates a torque of a wrong sign.
      !          sign is reversed at the end of the theta blocking.
      if ( this%nR == n_r_cmb .and. l_mag_LF .and. l_rot_ma .and. l_cond_ma ) then
         call get_lorentz_torque(this%lorentz_torque_ma,this%gsa_dist%brc,this%gsa_dist%bpc,this%nR)
      end if
      

      !--------- Calculate courant condition parameters:
      if ( .not. l_full_sphere .or. this%nR /= n_r_icb ) then
         call courant(this%nR,this%dtrkc,this%dthkc,this%gsa_dist%vrc,          &
              &       this%gsa_dist%vtc,this%gsa_dist%vpc,this%gsa_dist%brc,this%gsa_dist%btc, &
              &       this%gsa_dist%bpc, tscheme%courfac, tscheme%alffac)
      end if
      

      !--------- Since the fields are given at gridpoints here, this is a good
      !          point for graphical output:
      if ( this%l_graph ) then
#ifdef WITH_MPI
            call graphOut_mpi(time,this%nR,this%gsa_dist%vrc,this%gsa_dist%vtc,           &
                 &            this%gsa_dist%vpc,this%gsa_dist%brc,this%gsa_dist%btc,           &
                 &            this%gsa_dist%bpc,this%gsa_dist%sc,this%gsa_dist%pc,this%gsa_dist%xic,&
                 &            lGraphHeader)
#else
            call graphOut(time,this%nR,this%gsa_dist%vrc,this%gsa_dist%vtc,           &
                 &        this%gsa_dist%vpc,this%gsa_dist%brc,this%gsa_dist%btc,           &
                 &        this%gsa_dist%bpc,this%gsa_dist%sc,this%gsa_dist%pc,this%gsa_dist%xic,&
                 &        lGraphHeader)
#endif
      end if

      call this%nl_lm_dist%gather_all(this%nl_lm)
      call this%gsa_dist%gather_all(this%gsa)
      
      if ( this%l_probe_out ) then
         call probe_out(time,this%nR,this%gsa%vpc,this%gsa%brc,this%gsa%btc,1, &
              &         this%sizeThetaB)
      end if

      !--------- Helicity output:
      if ( this%lHelCalc ) then
         HelLMr(:)   =0.0_cp
         Hel2LMr(:)  =0.0_cp
         HelnaLMr(:) =0.0_cp
         Helna2LMr(:)=0.0_cp
         call get_helicity(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,         &
              &            this%gsa%cvrc,this%gsa%dvrdtc,this%gsa%dvrdpc,  &
              &            this%gsa%dvtdrc,this%gsa%dvpdrc,HelLMr,Hel2LMr, &
              &            HelnaLMr,Helna2LMr,this%nR,1 )
      end if

      !-- Viscous heating:
      if ( this%lPowerCalc ) then
         viscLMr(:)=0.0_cp
         call get_visc_heat(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,          &
              &             this%gsa%cvrc,this%gsa%dvrdrc,this%gsa%dvrdtc,   &
              &             this%gsa%dvrdpc,this%gsa%dvtdrc,this%gsa%dvtdpc, &
              &             this%gsa%dvpdrc,this%gsa%dvpdpc,viscLMr,         &
              &             this%nR,1)
      end if

      !-- horizontal velocity :
      if ( this%lViscBcCalc ) then
         gradsLMr(:)=0.0_cp
         uhLMr(:)   =0.0_cp
         duhLMr(:)  =0.0_cp
         call get_nlBLayers(this%gsa%vtc,this%gsa%vpc,this%gsa%dvtdrc,    &
              &             this%gsa%dvpdrc,this%gsa%drSc,this%gsa%dsdtc, &
              &             this%gsa%dsdpc,uhLMr,duhLMr,gradsLMr,nR,1 )
      end if

      !-- Radial flux profiles
      if ( this%lFluxProfCalc ) then
         fconvLMr(:)=0.0_cp
         fkinLMr(:) =0.0_cp
         fviscLMr(:)=0.0_cp
         fpoynLMr(:)=0.0_cp
         fresLMr(:) =0.0_cp
         call get_fluxes(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,            &
              &          this%gsa%dvrdrc,this%gsa%dvtdrc,this%gsa%dvpdrc,   &
              &          this%gsa%dvrdtc,this%gsa%dvrdpc,this%gsa%sc,       &
              &          this%gsa%pc,this%gsa%brc,this%gsa%btc,this%gsa%bpc,&
              &          this%gsa%cbtc,this%gsa%cbpc,fconvLMr,fkinLMr,      &
              &          fviscLMr,fpoynLMr,fresLMr,nR,1 )
      end if

      !-- Kinetic energy parallel and perpendicular to rotation axis
      if ( this%lPerpParCalc ) then
         EperpLMr(:)   =0.0_cp
         EparLMr(:)    =0.0_cp
         EperpaxiLMr(:)=0.0_cp
         EparaxiLMr(:) =0.0_cp
         call get_perpPar(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,EperpLMr, &
              &           EparLMr,EperpaxiLMr,EparaxiLMr,nR,1 )
      end if


      !--------- Movie output:
      if ( this%l_frame .and. l_movie_oc .and. l_store_frame ) then
         call store_movie_frame(this%nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc, &
              &                 this%gsa%brc,this%gsa%btc,this%gsa%bpc,         &
              &                 this%gsa%sc,this%gsa%drSc,this%gsa%dvrdpc,      &
              &                 this%gsa%dvpdrc,this%gsa%dvtdrc,this%gsa%dvrdtc,&
              &                 this%gsa%cvrc,this%gsa%cbrc,this%gsa%cbtc,1,    &
              &                 this%sizeThetaB,this%leg_helper%bCMB)
      end if


      !--------- Stuff for special output:
      !--------- Calculation of magnetic field production and advection terms
      !          for graphic output:
      if ( l_dtB ) then
         call get_dtBLM(this%nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,       &
              &         this%gsa%brc,this%gsa%btc,this%gsa%bpc,               &
              &         1 ,this%sizeThetaB,this%dtB_arrays%BtVrLM,            &
              &         this%dtB_arrays%BpVrLM,this%dtB_arrays%BrVtLM,        &
              &         this%dtB_arrays%BrVpLM,this%dtB_arrays%BtVpLM,        &
              &         this%dtB_arrays%BpVtLM,this%dtB_arrays%BrVZLM,        &
              &         this%dtB_arrays%BtVZLM,this%dtB_arrays%BtVpCotLM,     &
              &         this%dtB_arrays%BpVtCotLM,this%dtB_arrays%BtVZcotLM,  &
              &         this%dtB_arrays%BtVpSn2LM,this%dtB_arrays%BpVtSn2LM,  &
              &         this%dtB_arrays%BtVZsn2LM)
      end if


      !--------- Torsional oscillation terms:
      PERFON('TO_terms')
      if ( ( this%lTONext .or. this%lTONext2 ) .and. l_mag ) then
         call getTOnext(this%leg_helper%zAS,this%gsa%brc,this%gsa%btc,           &
              &         this%gsa%bpc,this%lTONext,this%lTONext2,tscheme%dt(1),   &
              &         dtLast,this%nR,1,this%sizeThetaB,this%BsLast,this%BpLast,&
              &         this%BzLast)
      end if

      if ( this%lTOCalc ) then
         call getTO(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc,   &
              &     this%gsa%dvpdrc,this%gsa%brc,this%gsa%btc,this%gsa%bpc, &
              &     this%gsa%cbrc,this%gsa%cbtc,this%BsLast,this%BpLast,    &
              &     this%BzLast,this%TO_arrays%dzRstrLM,                    &
              &     this%TO_arrays%dzAstrLM,this%TO_arrays%dzCorLM,         &
              &     this%TO_arrays%dzLFLM,dtLast,this%nR,1,this%sizeThetaB)
      end if
      PERFOFF

      lorentz_torque_ic = lorentz_torques_ic
      this%lorentz_torque_ic = lorentz_torques_ic
      lorentz_torque_ma = this%lorentz_torque_ma

      if (DEBUG_OUTPUT) then
         call this%nl_lm%output()
      end if

      !-- Partial calculation of time derivatives (horizontal parts):
      !   input flm...  is in (l,m) space at radial grid points this%nR !
      !   Only dVxBh needed for boundaries !
      !   get_td finally calculates the d*dt terms needed for the
      !   time step performed in s_LMLoop.f . This should be distributed
      !   over the different models that s_LMLoop.f parallelizes over.
      !PERFON('get_td')
      call td_counter%start_count()
      call this%nl_lm%get_td(this%nR, this%nBc, this%lRmsCalc,           &
           &                 this%lPressNext, dVSrLM, dVXirLM,           &
           &                 dVxVhLM, dVxBhLM, dwdt, dzdt, dpdt, dsdt,   &
           &                 dxidt, dbdt, djdt)
      call td_counter%stop_count(l_increment=.false.)

      !PERFOFF
      !-- Finish calculation of TO variables:
      if ( this%lTOcalc ) then
         call getTOfinish(this%nR, dtLast, this%leg_helper%zAS,             &
              &           this%leg_helper%dzAS, this%leg_helper%ddzAS,      &
              &           this%TO_arrays%dzRstrLM, this%TO_arrays%dzAstrLM, &
              &           this%TO_arrays%dzCorLM, this%TO_arrays%dzLFLM)
      end if

      !--- Form partial horizontal derivaties of magnetic production and
      !    advection terms:
      if ( l_dtB ) then
         PERFON('dtBLM')
         call get_dH_dtBLM(this%nR,this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM,&
              &            this%dtB_arrays%BrVtLM,this%dtB_arrays%BrVpLM,        &
              &            this%dtB_arrays%BtVpLM,this%dtB_arrays%BpVtLM,        &
              &            this%dtB_arrays%BrVZLM,this%dtB_arrays%BtVZLM,        &
              &            this%dtB_arrays%BtVpCotLM,this%dtB_arrays%BpVtCotLM,  &
              &            this%dtB_arrays%BtVpSn2LM,this%dtB_arrays%BpVtSn2LM)
         PERFOFF
      end if
    end subroutine do_iteration_ThetaBlocking_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space_shtns(this, gsa)

      class(rIterThetaBlocking_shtns_t) :: this
      type(grid_space_arrays_t) :: gsa

      !-- Local variables
      integer :: nR
      nR = this%nR

      if ( l_conv .or. l_mag_kin ) then
         if ( l_heat ) then
            call scal_to_spat(s_Rloc(:,nR), gsa%sc, l_R(nR))
            if ( this%lViscBcCalc ) then
               call scal_to_grad_spat(s_Rloc(:,nR), gsa%dsdtc, gsa%dsdpc, l_R(nR))
               if (this%nR == n_r_cmb .and. ktops==1) then
                  gsa%dsdtc=0.0_cp
                  gsa%dsdpc=0.0_cp
               end if
               if (this%nR == n_r_icb .and. kbots==1) then
                  gsa%dsdtc=0.0_cp
                  gsa%dsdpc=0.0_cp
               end if
            end if
         end if

         if ( this%lRmsCalc ) call scal_to_grad_spat(p_Rloc(:,nR), gsa%dpdtc, &
                                   &                 gsa%dpdpc, l_R(nR))

         !-- Pressure
         if ( this%lPressCalc ) call scal_to_spat(p_Rloc(:,nR), gsa%pc, l_R(nR))

         !-- Composition
         if ( l_chemical_conv ) call scal_to_spat(xi_Rloc(:,nR), gsa%xic, l_R(nR))

         if ( l_HT .or. this%lViscBcCalc ) then
            call scal_to_spat(ds_Rloc(:,nR), gsa%drsc, l_R(nR))
         endif
         if ( this%nBc == 0 ) then ! Bulk points
            !-- pol, sph, tor > ur,ut,up
            call torpol_to_spat(w_Rloc(:,nR), dw_Rloc(:,nR),  z_Rloc(:,nR), &
                 &              gsa%vrc, gsa%vtc, gsa%vpc, l_R(nR))

            !-- Advection is treated as u \times \curl u
            if ( l_adv_curl ) then
               !-- z,dz,w,dd< -> wr,wt,wp
               call torpol_to_curl_spat(or2(nR), w_Rloc(:,nR), ddw_Rloc(:,nR), &
                    &                   z_Rloc(:,nR), dz_Rloc(:,nR),           &
                    &                   gsa%cvrc, gsa%cvtc, gsa%cvpc, l_R(nR))

               !-- For some outputs one still need the other terms
               if ( this%lViscBcCalc .or. this%lPowerCalc .or. this%lRmsCalc  &
               &    .or. this%lFluxProfCalc .or. this%lTOCalc .or.            &
               &    ( this%l_frame .and. l_movie_oc .and. l_store_frame) ) then

                  call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR),         &
                       &              dz_Rloc(:,nR), gsa%dvrdrc, gsa%dvtdrc, &
                       &              gsa%dvpdrc, l_R(nR))
                  call pol_to_grad_spat(w_Rloc(:,nR),gsa%dvrdtc,gsa%dvrdpc, l_R(nR))
                  call torpol_to_dphspat(dw_Rloc(:,nR),  z_Rloc(:,nR), &
                       &                 gsa%dvtdpc, gsa%dvpdpc, l_R(nR))
               end if

            else ! Advection is treated as u\grad u

               call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR), dz_Rloc(:,nR), &
                 &              gsa%dvrdrc, gsa%dvtdrc, gsa%dvpdrc, l_R(nR))

               call pol_to_curlr_spat(z_Rloc(:,nR), gsa%cvrc, l_R(nR))

               call pol_to_grad_spat(w_Rloc(:,nR), gsa%dvrdtc, gsa%dvrdpc, l_R(nR))
               call torpol_to_dphspat(dw_Rloc(:,nR),  z_Rloc(:,nR), &
                    &                 gsa%dvtdpc, gsa%dvpdpc, l_R(nR))
            end if

         else if ( this%nBc == 1 ) then ! Stress free
             ! TODO don't compute vrc as it is set to 0 afterward
            call torpol_to_spat(w_Rloc(:,nR), dw_Rloc(:,nR),  z_Rloc(:,nR), &
                 &              gsa%vrc, gsa%vtc, gsa%vpc, l_R(nR))
            gsa%vrc = 0.0_cp
            if ( this%lDeriv ) then
               gsa%dvrdtc = 0.0_cp
               gsa%dvrdpc = 0.0_cp
               call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR), dz_Rloc(:,nR), &
                    &              gsa%dvrdrc, gsa%dvtdrc, gsa%dvpdrc, l_R(nR))
               call pol_to_curlr_spat(z_Rloc(:,nR), gsa%cvrc, l_R(nR))
               call torpol_to_dphspat(dw_Rloc(:,nR),  z_Rloc(:,nR), &
                    &                 gsa%dvtdpc, gsa%dvpdpc, l_R(nR))
            end if
         else if ( this%nBc == 2 ) then
            if ( this%nR == n_r_cmb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaMA,this%lDeriv, &
                    &                gsa%vrc,gsa%vtc,gsa%vpc,gsa%cvrc,gsa%dvrdtc, &
                    &                gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,1)
            else if ( this%nR == n_r_icb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaIC,this%lDeriv, &
                    &                gsa%vrc,gsa%vtc,gsa%vpc,gsa%cvrc,gsa%dvrdtc, &
                    &                gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,1)
            end if
            if ( this%lDeriv ) then
               call torpol_to_spat(dw_Rloc(:,nR), ddw_Rloc(:,nR), dz_Rloc(:,nR), &
                    &              gsa%dvrdrc, gsa%dvtdrc, gsa%dvpdrc, l_R(nR))
            end if
         end if
      end if

      if ( l_mag .or. l_mag_LF ) then
         call torpol_to_spat(b_Rloc(:,nR), db_Rloc(:,nR),  aj_Rloc(:,nR),    &
              &              gsa%brc, gsa%btc, gsa%bpc, l_R(nR))

         if ( this%lDeriv ) then
            call torpol_to_curl_spat(or2(nR), b_Rloc(:,nR), ddb_Rloc(:,nR), &
                 &                   aj_Rloc(:,nR), dj_Rloc(:,nR),          &
                 &                   gsa%cbrc, gsa%cbtc, gsa%cbpc, l_R(nR))
         end if
      end if

   end subroutine transform_to_grid_space_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space_shtns(this, gsa_dist, nl_lm_dist)

      class(rIterThetaBlocking_shtns_t) :: this
      type(grid_space_arrays_t) :: gsa_dist
!       type(grid_space_arrays_t) :: gsa
!       type(nonlinear_lm_t) :: nl_lm
      type(nonlinear_lm_t) :: nl_lm_dist

      ! Local variables
      integer :: nTheta, nPhi, nThStart, nThStop

      call shtns_load_cfg(1)

      if ( (.not.this%isRadialBoundaryPoint .or. this%lRmsCalc) &
            .and. ( l_conv_nl .or. l_mag_LF ) ) then

         !$omp parallel default(shared) private(nThStart,nThStop,nTheta,nPhi)
         nThStart=nThetaStart; nThStop=nThetaStop
         call get_openmp_blocks(nThStart,nThStop)

         !PERFON('inner1')
         if ( l_conv_nl .and. l_mag_LF ) then
            if ( this%nR>n_r_LCR ) then
               do nTheta=nThStart, nThStop
                  do nPhi=1, n_phi_max
                     gsa_dist%Advr(nPhi, nTheta)=gsa_dist%Advr(nPhi, nTheta) + gsa_dist%LFr(nPhi, nTheta)
                     gsa_dist%Advt(nPhi, nTheta)=gsa_dist%Advt(nPhi, nTheta) + gsa_dist%LFt(nPhi, nTheta)
                     gsa_dist%Advp(nPhi, nTheta)=gsa_dist%Advp(nPhi, nTheta) + gsa_dist%LFp(nPhi, nTheta)
                  end do
               end do
            end if
         else if ( l_mag_LF ) then
            if ( this%nR > n_r_LCR ) then
               do nTheta=nThStart, nThStop
                  do nPhi=1, n_phi_max
                     gsa_dist%Advr(nPhi, nTheta) = gsa_dist%LFr(nPhi, nTheta)
                     gsa_dist%Advt(nPhi, nTheta) = gsa_dist%LFt(nPhi, nTheta)
                     gsa_dist%Advp(nPhi, nTheta) = gsa_dist%LFp(nPhi, nTheta)
                  end do
               end do
            else
               do nTheta=nThStart, nThStop
                  do nPhi=1, n_phi_max
                     gsa_dist%Advr(nPhi,nTheta)=0.0_cp
                     gsa_dist%Advt(nPhi,nTheta)=0.0_cp
                     gsa_dist%Advp(nPhi,nTheta)=0.0_cp
                  end do
               end do
            end if
         end if

         if ( l_precession ) then
            do nTheta=nThStart, nThStop
               do nPhi=1, n_phi_max
                  gsa_dist%Advr(nPhi, nTheta)=gsa_dist%Advr(nPhi, nTheta) + gsa_dist%PCr(nPhi, nTheta)
                  gsa_dist%Advt(nPhi, nTheta)=gsa_dist%Advt(nPhi, nTheta) + gsa_dist%PCt(nPhi, nTheta)
                  gsa_dist%Advp(nPhi, nTheta)=gsa_dist%Advp(nPhi, nTheta) + gsa_dist%PCp(nPhi, nTheta)
               end do
            end do
         end if

         if ( l_centrifuge ) then
            do nTheta=nThStart, nThStop
               do nPhi=1, n_phi_max
                  gsa_dist%Advr(nPhi, nTheta)=gsa_dist%Advr(nPhi, nTheta) + gsa_dist%CAr(nPhi, nTheta)
                  gsa_dist%Advt(nPhi, nTheta)=gsa_dist%Advt(nPhi, nTheta) + gsa_dist%CAt(nPhi, nTheta)
               end do
            end do
         end if
         !$omp end parallel
         
         call spat_to_SH_dist(gsa_dist%Advr, nl_lm_dist%AdvrLM, l_R(this%nR))
         call spat_to_SH_dist(gsa_dist%Advt, nl_lm_dist%AdvtLM, l_R(this%nR))
         call spat_to_SH_dist(gsa_dist%Advp, nl_lm_dist%AdvpLM, l_R(this%nR))

         if ( this%lRmsCalc .and. l_mag_LF .and. this%nR>n_r_LCR ) then
            ! LF treated extra:
            call spat_to_SH_dist(gsa_dist%LFr, nl_lm_dist%LFrLM, l_R(this%nR))
            call spat_to_SH_dist(gsa_dist%LFt, nl_lm_dist%LFtLM, l_R(this%nR))
            call spat_to_SH_dist(gsa_dist%LFp, nl_lm_dist%LFpLM, l_R(this%nR))
         end if
         !PERFOFF
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_heat ) then
         !PERFON('inner2')
         call spat_to_qst_dist(gsa_dist%VSr, gsa_dist%VSt, gsa_dist%VSp, nl_lm_dist%VSrLM, nl_lm_dist%VStLM, &
              &           nl_lm_dist%VSpLM, l_R(this%nR))
         
         if (l_anel) then ! anelastic stuff
            if ( l_mag_nl .and. this%nR>n_r_LCR ) then
               call spat_to_SH_dist(gsa_dist%ViscHeat, nl_lm_dist%ViscHeatLM, l_R(this%nR))
               call spat_to_SH_dist(gsa_dist%OhmLoss,  nl_lm_dist%OhmLossLM, l_R(this%nR))
            else
               call spat_to_SH_dist(gsa_dist%ViscHeat, nl_lm_dist%ViscHeatLM, l_R(this%nR))
            end if
         end if
         !PERFOFF
      end if
      
      if ( (.not.this%isRadialBoundaryPoint) .and. l_chemical_conv ) then
         call spat_to_qst_dist(gsa_dist%VXir, gsa_dist%VXit, gsa_dist%VXip, nl_lm_dist%VXirLM, &
              &           nl_lm_dist%VXitLM, nl_lm_dist%VXipLM, l_R(this%nR))
      end if
      if ( l_mag_nl ) then
         !PERFON('mag_nl')
         if ( .not.this%isRadialBoundaryPoint .and. this%nR>n_r_LCR ) then
            call spat_to_qst_dist(gsa_dist%VxBr, gsa_dist%VxBt, gsa_dist%VxBp, nl_lm_dist%VxBrLM, &
                 &           nl_lm_dist%VxBtLM, nl_lm_dist%VxBpLM, l_R(this%nR))
         else
            call spat_to_sphertor_dist(gsa_dist%VxBt,gsa_dist%VxBp,nl_lm_dist%VxBtLM,nl_lm_dist%VxBpLM, &
                 &                l_R(this%nR))
         end if
         !PERFOFF
      end if

      if ( this%lRmsCalc ) then
         call spat_to_sphertor_dist(gsa_dist%dpdtc, gsa_dist%dpdpc, nl_lm_dist%PFt2LM, nl_lm_dist%PFp2LM, &
              &                l_R(this%nR))
         call spat_to_sphertor_dist(gsa_dist%CFt2, gsa_dist%CFp2, nl_lm_dist%CFt2LM, nl_lm_dist%CFp2LM, &
              &                l_R(this%nR))
         call spat_to_qst_dist(gsa_dist%dtVr, gsa_dist%dtVt, gsa_dist%dtVp, nl_lm_dist%dtVrLM, &
              &           nl_lm_dist%dtVtLM, nl_lm_dist%dtVpLM, l_R(this%nR))
         if ( l_conv_nl ) then
            call spat_to_sphertor_dist(gsa_dist%Advt2, gsa_dist%Advp2, nl_lm_dist%Advt2LM, &
                 &                nl_lm_dist%Advp2LM, l_R(this%nR))
         end if
         if ( l_adv_curl ) then !-- Kinetic pressure : 1/2 d u^2 / dr
            call spat_to_SH_dist(gsa_dist%dpkindrc, nl_lm_dist%dpkindrLM, l_R(this%nR))
         end if
         if ( l_mag_nl .and. this%nR>n_r_LCR ) then
            call spat_to_sphertor_dist(gsa_dist%LFt2, gsa_dist%LFp2, nl_lm_dist%LFt2LM, &
                 &                nl_lm_dist%LFp2LM, l_R(this%nR))
         end if
      end if
      
      call shtns_load_cfg(0)

   end subroutine transform_to_lm_space_shtns
!-------------------------------------------------------------------------------
end module rIterThetaBlocking_shtns_mod
