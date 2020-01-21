#include "perflib_preproc.cpp"
module rIterThetaBlocking_seq_mod

   use precision_mod
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use TO_arrays_mod, only: TO_arrays_t
   use dtB_arrays_mod, only: dtB_arrays_t
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use num_param, only: lm2phy_counter, phy2lm_counter, nl_counter, td_counter
   use truncation, only: lm_max,lmP_max, nrp, l_max, lmP_max_dtB, &
       &                 n_phi_maxStr, n_theta_maxStr, n_r_maxStr,&
       &                 n_r_cmb, n_r_icb
   use blocking, only: nfs
   use logic, only: l_mag, l_conv, l_mag_kin, l_heat, l_ht, l_anel, l_mag_LF,&
       &            l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic,   &
       &            l_cond_ic, l_rot_ma, l_cond_ma, l_dtB, l_store_frame,    &
       &            l_movie_oc, l_TO, l_probe, l_full_sphere
   use radial_functions, only: or2, orho1
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
   use nonlinear_bcs, only: get_br_v_bcs
   use constants, only: zero
   use time_schemes, only: type_tscheme
   use nl_special_calc
   use probe_mod

   implicit none

   private

   type, public, extends(rIterThetaBlocking_t) :: rIterThetaBlocking_seq_t
      type(grid_space_arrays_t) :: gsa
      type(TO_arrays_t) :: TO_arrays
      type(dtB_arrays_t) :: dtB_arrays
      type(nonlinear_lm_t) :: nl_lm
   contains
      procedure :: initialize => initialize_rIterThetaBlocking_seq
      procedure :: finalize => finalize_rIterThetaBlocking_seq
      procedure :: do_iteration => do_iteration_ThetaBlocking_seq
      procedure :: getType => getThisType
   end type rIterThetaBlocking_seq_t

contains

   function getThisType(this)

      class(rIterThetaBlocking_seq_t) :: this
      character(len=100) :: getThisType

      getThisType="rIterThetaBlocking_seq_t"

   end function getThisType
!------------------------------------------------------------------------------
   subroutine initialize_rIterThetaBlocking_seq(this)

      class(rIterThetaBlocking_seq_t) :: this

      call this%allocate_common_arrays()
      call this%gsa%initialize()
      call this%nl_lm%initialize(lmP_max)
      if ( l_TO ) call this%TO_arrays%initialize()
      call this%dtB_arrays%initialize()

   end subroutine initialize_rIterThetaBlocking_seq
!------------------------------------------------------------------------------
   subroutine finalize_rIterThetaBlocking_seq(this)

      class(rIterThetaBlocking_seq_t) :: this

      call this%deallocate_common_arrays()
      call this%gsa%finalize()
      call this%nl_lm%finalize()
      if ( l_TO ) call this%TO_arrays%finalize()
      call this%dtB_arrays%finalize()

   end subroutine finalize_rIterThetaBlocking_seq
!------------------------------------------------------------------------------
   subroutine do_iteration_ThetaBlocking_seq(this,nR,nBc,time,timeStage, &
              &           tscheme,dtLast,dsdt,dwdt,dzdt,dpdt,dxidt,dbdt, &
              &           djdt,dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,           &
              &           br_vt_lm_cmb,br_vp_lm_cmb,br_vt_lm_icb,        &
              &           br_vp_lm_icb,lorentz_torque_ic,                &
              &           lorentz_torque_ma,HelLMr,Hel2LMr,HelnaLMr,     &
              &           Helna2LMr,viscLMr,uhLMr,duhLMr,gradsLMr,       &
              &           fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr,    &
              &           EperpLMr,EparLMr,EperpaxiLMr,EparaxiLmr)

      class(rIterThetaBlocking_seq_t) :: this

      !-- Input variables
      integer,             intent(in) :: nR,nBc
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time,timeStage,dtLast

      !-- Output variables
      complex(cp), intent(out) :: dwdt(:),dzdt(:),dpdt(:),dsdt(:),dVSrLM(:)
      complex(cp), intent(out) :: dxidt(:), dVXirLM(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:),dVxVhLM(:),dVxBhLM(:)
      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic
      real(cp),    intent(out) :: HelLMr(:),Hel2LMr(:),HelnaLMr(:),Helna2LMr(:)
      real(cp),    intent(out) :: viscLMr(:)
      real(cp),    intent(out) :: uhLMr(:),duhLMr(:),gradsLMr(:)
      real(cp),    intent(out) :: fconvLMr(:),fkinLMr(:),fviscLMr(:)
      real(cp),    intent(out) :: fpoynLMr(:),fresLMr(:)
      real(cp),    intent(out) :: EperpLMr(:),EparLMr(:),EperpaxiLMr(:),EparaxiLMr(:)

      !-- Local variables
      integer :: l,lm,nThetaB,nThetaLast,nThetaStart
      logical :: lGraphHeader=.false.
      logical :: DEBUG_OUTPUT=.false.

      this%nR=nR
      this%nBc=nBc
      this%isRadialBoundaryPoint = (nR == n_r_cmb) .or. (nR == n_r_icb)

      this%dtrkc=1.e10_cp
      this%dthkc=1.e10_cp

      if ( this%lTOCalc ) then
         !------ Zero lm coeffs for first theta block:
         do l=0,l_max
            this%TO_arrays%dzRstrLM(l+1)=0.0_cp
            this%TO_arrays%dzAstrLM(l+1)=0.0_cp
            this%TO_arrays%dzCorLM(l+1) =0.0_cp
            this%TO_arrays%dzLFLM(l+1)  =0.0_cp
         end do
      end if

      !----- Prepare legendre transform:
      !      legPrepG collects all the different modes necessary
      !      to calculate the non-linear terms at a radial grid point nR
      PERFON('legPrepG')
      if ( DEBUG_OUTPUT ) then
         write(*,"(I3,A,I1,2(A,L1))") this%nR,": nBc = ",this%nBc,", lDeriv = ", &
              & this%lDeriv,", l_mag = ",l_mag
      end if

      call this%leg_helper%legPrepG(this%nR,this%nBc,this%lDeriv,this%lRmsCalc, &
           &                        this%l_frame,this%lTOnext,this%lTOnext2,    &
           &                        this%lTOcalc)
      PERFOFF

      lorentz_torque_ma = 0.0_cp
      lorentz_torque_ic = 0.0_cp
      br_vt_lm_cmb=zero
      br_vp_lm_cmb=zero
      br_vt_lm_icb=zero
      br_vp_lm_icb=zero
      HelLMr=0.0_cp
      Hel2LMr=0.0_cp
      HelnaLMr=0.0_cp
      Helna2LMr=0.0_cp
      uhLMr = 0.0_cp
      duhLMr = 0.0_cp
      viscLMr = 0.0_cp
      gradsLMr = 0.0_cp
      fconvLMr=0.0_cp
      fkinLMr=0.0_cp
      fviscLMr=0.0_cp
      fpoynLMr=0.0_cp
      fresLMr=0.0_cp
      EperpLMr=0.0_cp
      EparLMr=0.0_cp
      EperpaxiLMr=0.0_cp
      EparaxiLMr=0.0_cp
      call this%nl_lm%set_zero()

      !if (DEBUG_OUTPUT) then
      !   write(*,"(I3,A,44ES20.12)") this%nR,": legPrepG results = ",&
      !        & SUM(dLhw),SUM(dLhdw),SUM(dLhz),SUM(vhG),SUM(vhC),&
      !        & SUM(dvhdrG),SUM(dvhdrC),&
      !        & SUM(dLhb),SUM(dLhj),SUM(bhG),SUM(bhC),SUM(cbhG),SUM(cbhC),&
      !        & SUM(sR),SUM(dsR),SUM(preR),SUM(dpR),SUM(zAS), SUM(dzAS),&
      !        &SUM(ddzAS),SUM(bCMB),omegaIC,omegaMA

      !   write(*,"(A,I4,A)") "We have ",this%nThetaBs," theta blocks!"
      !end if
      !----- Blocking of loops over ic (theta):
      do nThetaB=1,this%nThetaBs
         nThetaLast =(nThetaB-1) * this%sizeThetaB
         nThetaStart=nThetaLast+1

         call lm2phy_counter%start_count()
         call this%transform_to_grid_space(nThetaStart,this%gsa)
         call lm2phy_counter%stop_count(l_increment=.false.)

         !--------- Calculation of nonlinear products in grid space:
         if ( (.not.this%isRadialBoundaryPoint) .or. this%lMagNlBc .or. &
                this%lRmsCalc ) then
            !write(*,"(I4,A,ES20.13)") this%nR,", vp = ",sum(real(conjg(vpc)*vpc))
            call nl_counter%start_count()
            PERFON('get_nl')
            call this%gsa%get_nl(timeStage, tscheme, this%nR, this%nBc, nThetaStart, &
                 &               this%lRmsCalc)
            PERFOFF
            call nl_counter%stop_count(l_increment=.false.)

            call phy2lm_counter%start_count()
            call this%transform_to_lm_space(nThetaStart,this%gsa,this%nl_lm)
            call phy2lm_counter%stop_count(l_increment=.false.)

         else if ( l_mag ) then
            do lm=1,lmP_max
               this%nl_lm%VxBtLM(lm)=0.0_cp
               this%nl_lm%VxBpLM(lm)=0.0_cp
            end do
         end if

         !---- Calculation of nonlinear products needed for conducting mantle or
         !     conducting inner core if free stress BCs are applied:
         !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
         !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
         !     and br_vp_lm_cmb in lm-space, respectively the contribution
         !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
         !     These products are used in get_b_nl_bcs.
         PERFON('nl_cmb')
         if ( this%nR == n_r_cmb .and. l_b_nl_cmb ) then
            call get_br_v_bcs(this%gsa%brc,this%gsa%vtc,this%gsa%vpc, &
                 &            this%leg_helper%omegaMA,or2(this%nR),   &
                 &            orho1(this%nR),nThetaStart,             &
                 &            this%sizeThetaB,br_vt_lm_cmb,br_vp_lm_cmb)
         else if ( this%nR == n_r_icb .and. l_b_nl_icb ) then
            call get_br_v_bcs(this%gsa%brc,this%gsa%vtc,this%gsa%vpc, &
                 &            this%leg_helper%omegaIC,or2(this%nR),   &
                 &            orho1(this%nR),nThetaStart,             &
                 &            this%sizeThetaB,br_vt_lm_icb,br_vp_lm_icb)
         end if
         PERFOFF
         !--------- Calculate Lorentz torque on inner core:
         !          each call adds the contribution of the theta-block to
         !          lorentz_torque_ic
         PERFON('lorentz')
         if ( this%nR == n_r_icb .and. l_mag_LF .and. l_rot_ic .and. l_cond_ic  ) then
            call get_lorentz_torque(lorentz_torque_ic,nThetaStart,     &
                 &                  this%sizeThetaB,this%gsa%brc,      &
                 &                  this%gsa%bpc,this%nR)
         end if

         !--------- Calculate Lorentz torque on mantle:
         !          note: this calculates a torque of a wrong sign.
         !          sign is reversed at the end of the theta blocking.
         if ( this%nR == n_r_cmb .and. l_mag_LF .and. l_rot_ma .and. l_cond_ma ) then
            call get_lorentz_torque(lorentz_torque_ma,nThetaStart,     &
                 &                  this%sizeThetaB,this%gsa%brc,      &
                 &                  this%gsa%bpc,this%nR)
         end if
         PERFOFF
         !--------- Calculate courant condition parameters:
         if ( .not. l_full_sphere .or. this%nR /= n_r_icb  ) then
            call courant(this%nR,this%dtrkc,this%dthkc,this%gsa%vrc, &
                 &       this%gsa%vtc,this%gsa%vpc,this%gsa%brc,     &
                 &       this%gsa%btc,this%gsa%bpc,nThetaStart,      &
                 &       this%sizeThetaB, tscheme%courfac,           &
                 &       tscheme%alffac)
         end if

         !--------- Since the fields are given at gridpoints here, this is a good
         !          point for graphical output:
         if ( this%l_graph ) then
#ifdef WITH_MPI
            PERFON('graphout')
            call graphOut_mpi(time,this%nR,this%gsa%vrc,this%gsa%vtc,  &
                 &            this%gsa%vpc,this%gsa%brc,this%gsa%btc,  &
                 &            this%gsa%bpc,this%gsa%sc,this%gsa%pc,    &
                 &            this%gsa%xic,nThetaStart,this%sizeThetaB,&
                 &            lGraphHeader)
            PERFOFF
#else
            call graphOut(time,this%nR,this%gsa%vrc,this%gsa%vtc,   &
                 &        this%gsa%vpc,this%gsa%brc,this%gsa%btc,   &
                 &        this%gsa%bpc,this%gsa%sc, this%gsa%pc,    &
                 &        this%gsa%xic,nThetaStart,this%sizeThetaB, &
                 &        lGraphHeader)
#endif
         end if

         if ( this%l_probe_out ) then
            call probe_out(time,this%nR,this%gsa%vpc,this%gsa%brc,this%gsa%btc, &
            &              nThetaStart,this%sizeThetaB)
         end if

         !--------- Helicity output:
         if ( this%lHelCalc ) then
            PERFON('hel_out')
            call get_helicity(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,      &
                 &        this%gsa%cvrc,this%gsa%dvrdtc,this%gsa%dvrdpc,   &
                 &        this%gsa%dvtdrc,this%gsa%dvpdrc,HelLMr,Hel2LMr,  &
                 &        HelnaLMr,Helna2LMr,this%nR,nThetaStart)
            PERFOFF
         end if

         !--------- Viscous heating:
         if ( this%lPowerCalc ) then
            PERFON('hel_out')
            call get_visc_heat(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,     &
                 &        this%gsa%cvrc,this%gsa%dvrdrc,this%gsa%dvrdtc,   &
                 &        this%gsa%dvrdpc,this%gsa%dvtdrc,this%gsa%dvtdpc, &
                 &        this%gsa%dvpdrc,this%gsa%dvpdpc,viscLMr,         &
                 &        this%nR,nThetaStart)
            PERFOFF
         end if

         !--------- horizontal velocity :

         if ( this%lViscBcCalc ) then
            call get_nlBLayers(this%gsa%vtc,this%gsa%vpc,this%gsa%dvtdrc,  &
                    &          this%gsa%dvpdrc,this%gsa%drSc,              &
                    &          this%gsa%dsdtc,this%gsa%dsdpc,uhLMr,duhLMr, &
                    &          gradsLMr,nR,nThetaStart)
         end if

         if ( this%lFluxProfCalc ) then
             call get_fluxes(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%dvrdrc,  &
                    &        this%gsa%dvtdrc,this%gsa%dvpdrc,this%gsa%dvrdtc,         &
                    &        this%gsa%dvrdpc,this%gsa%sc,this%gsa%pc,this%gsa%brc,    &
                    &        this%gsa%btc,this%gsa%bpc,this%gsa%cbtc,this%gsa%cbpc,   &
                    &        fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr,nR,nThetaStart)
         end if

         if ( this%lPerpParCalc ) then
             call get_perpPar(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,  &
                    &        EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr,  &
                    &        nR,nThetaStart)
         end if

         !--------- Movie output:
         if ( this%l_frame .and. l_movie_oc .and. l_store_frame ) then
            PERFON('mov_out')
            call store_movie_frame(this%nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,  &
                 &                 this%gsa%brc,this%gsa%btc,this%gsa%bpc,          &
                 &                 this%gsa%sc,this%gsa%drSc,this%gsa%dvrdpc,       &
                 &                 this%gsa%dvpdrc,this%gsa%dvtdrc,                 &
                 &                 this%gsa%dvrdtc,this%gsa%cvrc,this%gsa%cbrc,     &
                 &                 this%gsa%cbtc,nThetaStart,this%sizeThetaB,       &
                 &                 this%leg_helper%bCMB)
            PERFOFF
         end if


         !--------- Stuff for special output:
         !--------- Calculation of magnetic field production and advection terms
         !          for graphic output:
         if ( l_dtB ) then
            PERFON('dtBLM')
            call get_dtBLM(this%nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,      &
                 &         this%gsa%brc,this%gsa%btc,this%gsa%bpc,              &
                 &         nThetaStart,this%sizeThetaB,                         &
                 &         this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM,       &
                 &         this%dtB_arrays%BrVtLM,this%dtB_arrays%BrVpLM,       &
                 &         this%dtB_arrays%BtVpLM,this%dtB_arrays%BpVtLM,       &
                 &         this%dtB_arrays%BrVZLM,this%dtB_arrays%BtVZLM,       &
                 &         this%dtB_arrays%BtVpCotLM,this%dtB_arrays%BpVtCotLM, &
                 &         this%dtB_arrays%BtVZcotLM,this%dtB_arrays%BtVpSn2LM, &
                 &         this%dtB_arrays%BpVtSn2LM,this%dtB_arrays%BtVZsn2LM)
            PERFOFF
         end if


         !--------- Torsional oscillation terms:
         PERFON('TO_terms')
         if ( ( this%lTONext .or. this%lTONext2 ) .and. l_mag ) then
            call getTOnext(this%leg_helper%zAS,this%gsa%brc,this%gsa%btc,        &
                 &         this%gsa%bpc,this%lTONext,this%lTONext2,tscheme%dt(1),&
                 &         dtLast,this%nR,nThetaStart,this%sizeThetaB,           &
                 &         this%BsLast,this%BpLast,this%BzLast)
         end if

         if ( this%lTOCalc ) then
            call getTO(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc,   &
                 &     this%gsa%dvpdrc,this%gsa%brc,this%gsa%btc,this%gsa%bpc, &
                 &     this%gsa%cbrc,this%gsa%cbtc,this%BsLast,this%BpLast,    &
                 &     this%BzLast,this%TO_arrays%dzRstrLM,                    &
                 &     this%TO_arrays%dzAstrLM,this%TO_arrays%dzCorLM,         &
                 &     this%TO_arrays%dzLFLM,dtLast,this%nR,nThetaStart,       &
                 &     this%sizeThetaB)
         end if
         PERFOFF

      end do ! Loop over theta blocks


      !-- Partial calculation of time derivatives (horizontal parts):
      !   input flm...  is in (l,m) space at radial grid points this%nR !
      !   Only dVxBh needed for boundaries !
      !   get_td finally calculates the d*dt terms needed for the
      !   time step performed in s_LMLoop.f . This should be distributed
      !   over the different models that s_LMLoop.f parallelizes over.
      !write(*,"(A,I4,4ES20.13)") "before_td: ",this%nR,SUM(this%nl_lm%VxBtLM),&
      !     & SUM(this%nl_lm%VxBpLM)
      call td_counter%start_count()
      PERFON('get_td')
      call this%nl_lm%get_td(this%nR,this%nBc,this%lRmsCalc,             &
           &                 this%lPressCalc,dVSrLM,dVXirLM,dVxVhLM,     &
           &                 dVxBhLM,dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt)
      PERFOFF
      call td_counter%stop_count(l_increment=.false.)
      !do lm=1,lm_max
      !   write(*,"(2(I3,A),2ES20.12)") this%nR,": dwdt(",lm,") = ",dwdt(lm)
      !end do
      !write(*,"(A,I4,4ES20.13)") "after_td: dwdt ",this%nR, SUM(dwdt)
      !-- Finish calculation of TO variables:
      if ( this%lTOcalc ) then
         call getTOfinish(this%nR,dtLast,this%leg_helper%zAS,             &
              &           this%leg_helper%dzAS,this%leg_helper%ddzAS,     &
              &           this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM,&
              &           this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM)
      end if

      !--- Form partial horizontal derivaties of magnetic production and
      !    advection terms:
      if ( l_dtB ) then
         call get_dH_dtBLM(this%nR,this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM, &
              &            this%dtB_arrays%BrVtLM,this%dtB_arrays%BrVpLM,         &
              &            this%dtB_arrays%BtVpLM,this%dtB_arrays%BpVtLM,         &
              &            this%dtB_arrays%BrVZLM,this%dtB_arrays%BtVZLM,         &
              &            this%dtB_arrays%BtVpCotLM,this%dtB_arrays%BpVtCotLM,   &
              &            this%dtB_arrays%BtVpSn2LM,this%dtB_arrays%BpVtSn2LM)
      end if

   end subroutine do_iteration_ThetaBlocking_seq
!------------------------------------------------------------------------------
end module rIterThetaBlocking_seq_mod
