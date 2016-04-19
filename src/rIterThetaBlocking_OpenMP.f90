#include "perflib_preproc.cpp"
module rIterThetaBlocking_OpenMP_mod
#ifdef WITHOMP
   use omp_lib
#endif
   use precision_mod
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t

   use truncation, only: lm_max, lmP_max, nrp, l_max, lmP_max_dtB,&
                         n_phi_maxStr, n_theta_maxStr, n_r_maxStr
   use blocking, only: nfs
   use logic, only: l_mag, l_conv, l_mag_kin, l_heat, l_ht, l_anel, l_mag_LF, &
                    l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic,    &
                    l_cond_ic, l_rot_ma, l_cond_ma, l_dtB, l_store_frame,     &
                    l_movie_oc, l_TO
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: or2, orho1
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
   use nonlinear_bcs, only: get_br_v_bcs
   use nl_special_calc

   implicit none

   private

   type, public, extends(rIterThetaBlocking_t) :: rIterThetaBlocking_OpenMP_t
      integer :: nThreads
      type(grid_space_arrays_t), allocatable :: gsa(:)
      type(TO_arrays_t), allocatable :: TO_arrays(:)
      type(dtB_arrays_t), allocatable :: dtB_arrays(:)
      type(nonlinear_lm_t), allocatable :: nl_lm(:)
      real(cp), allocatable :: lorentz_torque_ic(:),lorentz_torque_ma(:)
   contains
      procedure :: initialize => initialize_rIterThetaBlocking_OpenMP
      procedure :: finalize => finalize_rIterThetaBlocking_OpenMP
      procedure :: do_iteration => do_iteration_ThetaBlocking_OpenMP
      procedure :: getType => getThisType
   end type rIterThetaBlocking_OpenMP_t

contains

   function getThisType(this)

      class(rIterThetaBlocking_OpenMP_t) :: this
      character(len=100) :: getThisType
      getThisType="rIterThetaBlocking_OpenMP_t"

   end function getThisType
!------------------------------------------------------------------------------
   subroutine initialize_rIterThetaBlocking_OpenMP(this)

      class(rIterThetaBlocking_OpenMP_t) :: this
      integer :: threadid

#ifdef WITHOMP
      this%nThreads=omp_get_max_threads()
#else
      this%nThreads=1
#endif
      allocate(this%gsa(0:this%nThreads-1))
      allocate(this%nl_lm(0:this%nThreads-1))
      allocate(this%lorentz_torque_ic(0:this%nThreads-1))
      allocate(this%lorentz_torque_ma(0:this%nThreads-1))
      if ( l_TO ) allocate(this%TO_arrays(0:this%nThreads-1))
      allocate(this%dtB_arrays(0:this%nThreads-1))

      call this%allocate_common_arrays()
      !$OMP PARALLEL default(shared) shared(this,lmP_max) private(threadid)
#ifdef WITHOMP
      threadid = omp_get_thread_num()
#else
      threadid = 0
#endif
      call this%gsa(threadid)%initialize()
      if ( l_TO ) call this%TO_arrays(threadid)%initialize()
      call this%dtB_arrays(threadid)%initialize()
      call this%nl_lm(threadid)%initialize(lmP_max)
      !$OMP END PARALLEL
   end subroutine initialize_rIterThetaBlocking_OpenMP
!------------------------------------------------------------------------------
   subroutine finalize_rIterThetaBlocking_OpenMP(this)

      class(rIterThetaBlocking_OpenMP_t) :: this
      integer :: threadid

      call this%deallocate_common_arrays()
      !$OMP PARALLEL default(shared) shared(this) private(threadid)
#ifdef WITHOMP
      threadid = omp_get_thread_num()
#else
      threadid = 0
#endif
      call this%gsa(threadid)%finalize()
      call this%nl_lm(threadid)%finalize()
      !$OMP END PARALLEL
      deallocate(this%gsa)
      if ( l_TO ) deallocate(this%TO_arrays)
      deallocate(this%dtB_arrays)
      deallocate(this%nl_lm)
      deallocate(this%lorentz_torque_ic)
      deallocate(this%lorentz_torque_ma)

   end subroutine finalize_rIterThetaBlocking_OpenMP
!------------------------------------------------------------------------------
   subroutine do_iteration_ThetaBlocking_OpenMP(this,nR,nBc,time,dt,dtLast,&
        &                 dsdt,dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dVxBhLM,     &
        &                 dVSrLM,dVXirLM,br_vt_lm_cmb,br_vp_lm_cmb,        &
        &                 br_vt_lm_icb,br_vp_lm_icb,                       &
        &                 lorentz_torque_ic, lorentz_torque_ma,            &
        &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,  &
        &                 gradsLMr,fconvLMr,fkinLMr,fviscLMr,              &
        &                 fpoynLMr,fresLMr,EperpLMr,EparLMr,EperpaxiLMr,   &
        &                 EparaxiLMr)

      class(rIterThetaBlocking_OpenMP_t) :: this
      integer,  intent(in) :: nR,nBc
      real(cp), intent(in) :: time,dt,dtLast

      complex(cp), intent(out) :: dwdt(:),dzdt(:),dpdt(:),dsdt(:),dVSrLM(:)
      complex(cp), intent(out) :: dxidt(:), dVXirLM(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:),dVxBhLM(:)
      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic
      real(cp),    intent(out) :: HelLMr(:),Hel2LMr(:),HelnaLMr(:),Helna2LMr(:)
      real(cp),    intent(out) :: uhLMr(:),duhLMr(:),gradsLMr(:)
      real(cp),    intent(out) :: fconvLMr(:),fkinLMr(:),fviscLMr(:)
      real(cp),    intent(out) :: fpoynLMr(:),fresLMr(:)
      real(cp),    intent(out) :: EperpLMr(:),EparLMr(:),EperpaxiLMr(:),EparaxiLMr(:)

      integer :: l,lm,nThetaB,nThetaLast,nThetaStart,nThetaStop
      integer :: threadid,iThread
      logical :: lGraphHeader=.false.
      logical :: DEBUG_OUTPUT=.false.
      real(cp) :: lt,y,c,t,lorentz_torques_ic(this%nThetaBs)

      this%nR=nR
      this%nBc=nBc
      this%isRadialBoundaryPoint=(nR == n_r_cmb).or.(nR == n_r_icb)

      if ( this%l_cour ) then
         this%dtrkc=1.e10_cp
         this%dthkc=1.e10_cp
      end if
      !----- Prepare legendre transform:
      !      legPrepG collects all the different modes necessary
      !      to calculate the non-linear terms at a radial grid point nR
      !PERFON('legPrepG')
      if (DEBUG_OUTPUT) then
         write(*,"(I3,A,I1,2(A,L1))") this%nR,": nBc = ", &
              & this%nBc,", lDeriv = ",this%lDeriv,", l_mag = ",l_mag
      end if

      call this%leg_helper%legPrepG(this%nR,this%nBc,this%lDeriv,this%lRmsCalc, &
           &                        this%lPressCalc,this%l_frame,this%lTOnext,  &
           &                        this%lTOnext2,this%lTOcalc)
      !PERFOFF

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
      !$OMP PARALLEL default(shared) &
      !$OMP SHARED(this,l_mag,l_b_nl_cmb,l_b_nl_icb,l_mag_LF,l_rot_ic,l_cond_ic) &
      !$OMP SHARED(l_rot_ma,l_cond_ma,l_movie_oc,l_store_frame,l_dtB) &
      !$OMP SHARED(lmP_max,n_r_cmb,n_r_icb) &
      !$OMP SHARED(or2,orho1,time,dt,dtLast,DEBUG_OUTPUT) &
      !$OMP PRIVATE(threadid,nThetaLast,nThetaStart,nThetaStop,y,t,c,lt) &
      !$OMP shared(br_vt_lm_cmb,br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb) &
      !$OMP SHARED(lorentz_torques_ic) &
      !$OMP shared(HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,gradsLMr) &
      !$OMP shared(fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr) &
      !$OMP shared(EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)
#ifdef WITHOMP
      threadid = omp_get_thread_num()
#else
      threadid = 0
#endif
      this%lorentz_torque_ma(threadid) = 0.0_cp
      this%lorentz_torque_ic(threadid) = 0.0_cp
      if ( this%lTOCalc ) then
         !------ Zero lm coeffs for first theta block:
         call this%TO_arrays(threadid)%set_zero()
      end if
      if ( l_dtB ) then
         call this%dtB_arrays(threadid)%set_zero()
      end if

      c = 0.0_cp
      !$OMP SINGLE
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
      !$OMP END SINGLE
      !$OMP BARRIER
      call this%nl_lm(threadid)%set_zero()
      !$OMP do &
      !$OMP reduction(+:br_vt_lm_cmb,br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb) &
      !$OMP reduction(+:HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,gradsLMr) &
      !$OMP reduction(+:fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr) &
      !$OMP reduction(+:EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)

      do nThetaB=1,this%nThetaBs
         nThetaLast =(nThetaB-1) * this%sizeThetaB
         nThetaStart=nThetaLast+1
         nThetaStop =nThetaLast + this%sizeThetaB
         !write(*,"(I3,A,I4,A,I4)") nThetaB,". theta block from ", &
         !      &                  nThetaStart," to ",nThetaStop

         !PERFON('lm2grid')
         call this%transform_to_grid_space(nThetaStart,nThetaStop,&
              &                            this%gsa(threadid))
         !PERFOFF

         !--------- Calculation of nonlinear products in grid space:
         if ( (.not.this%isRadialBoundaryPoint) .or. this%lMagNlBc .or. &
               this%lRmsCalc ) then

            !if (DEBUG_OUTPUT) then
               !if (this%nR == 2) then
               !   write(*,"(A,I2,A,I2)") &
               !        &  "++++ START gsa(",threadid,") for nThetaB = ",nThetaB
               !   call this%gsa(threadid)%output_nl_input()
               !   write(*,"(A,I2,A,I2)") &
               !        & "---- END   gsa(",threadid,") for nThetaB = ",nThetaB
               !end if
            !end if

            !PERFON('get_nl')
            call this%gsa(threadid)%get_nl(this%nR, this%nBc, nThetaStart, &
                          this%lRmsCalc )
            !PERFOFF

            !if (DEBUG_OUTPUT) then
            !   if (this%nR == 2) then
            !      write(*,"(A,I2,A,I2)") &
            !           & "++++ START gsa(",threadid,") for nThetaB = ",nThetaB
            !      call this%gsa(threadid)%output()
            !      write(*,"(A,I2,A,I2)") &
            !           & "---- END   gsa(",threadid,") for nThetaB = ",nThetaB
            !   end if
            !end if
            !PERFON('grid2lm')
            call this%transform_to_lm_space(nThetaStart,nThetaStop, &
                 &                          this%gsa(threadid),this%nl_lm(threadid))
            !PERFOFF

         else if ( l_mag ) then
            do lm=1,lmP_max
               this%nl_lm(threadid)%VxBtLM(lm)=0.0_cp
               this%nl_lm(threadid)%VxBpLM(lm)=0.0_cp
            end do
         end if

         !---- Calculation of nonlinear products needed for conducting mantle or
         !     conducting inner core if free stress BCs are applied:
         !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
         !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
         !     and br_vp_lm_cmb in lm-space, respectively the contribution
         !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
         !     These products are used in get_b_nl_bcs.
         !PERFON('nl_cmb')
         if ( this%nR == n_r_cmb .and. l_b_nl_cmb ) then
            call get_br_v_bcs(this%gsa(threadid)%brc,this%gsa(threadid)%vtc,     &
                 &            this%gsa(threadid)%vpc,this%leg_helper%omegaMA,    &
                 &            or2(this%nR),orho1(this%nR),nThetaStart,           &
                 &            this%sizeThetaB,br_vt_lm_cmb,br_vp_lm_cmb)
         else if ( this%nR == n_r_icb .and. l_b_nl_icb ) then
            call get_br_v_bcs(this%gsa(threadid)%brc,this%gsa(threadid)%vtc,     &
                 &            this%gsa(threadid)%vpc,this%leg_helper%omegaIC,    &
                 &            or2(this%nR),orho1(this%nR),nThetaStart,           &
                 &            this%sizeThetaB,br_vt_lm_icb,br_vp_lm_icb)
         end if
         !PERFOFF
         !--------- Calculate Lorentz torque on inner core:
         !          each call adds the contribution of the theta-block to
         !          lorentz_torque_ic
         !PERFON('lorentz')
         if ( this%nR == n_r_icb .and. l_mag_LF .and. l_rot_ic .and. l_cond_ic  ) then
            lorentz_torques_ic(nThetaB)=0.0_cp
            call get_lorentz_torque(lorentz_torques_ic(nThetaB),     &
                 &                  nThetaStart,this%sizeThetaB,     &
                 &                  this%gsa(threadid)%brc,          &
                 &                  this%gsa(threadid)%bpc,this%nR)
            !y=lt-c
            !t=this%lorentz_torque_ic(threadid)+y
            !c=(t-this%lorentz_torque_ic(threadid))-y
            !this%lorentz_torque_ic(threadid)=t
         end if

         !--------- Calculate Lorentz torque on mantle:
         !          note: this calculates a torque of a wrong sign.
         !          sign is reversed at the end of the theta blocking.
         if ( this%nR == n_r_cmb .and. l_mag_LF .and. l_rot_ma .and. l_cond_ma ) then
            call get_lorentz_torque(this%lorentz_torque_ma(threadid),   &
                 &                  nThetaStart,this%sizeThetaB,        &
                 &                  this%gsa(threadid)%brc,             &
                 &                  this%gsa(threadid)%bpc,this%nR)
         end if
         !PERFOFF

         !--------- Calculate courant condition parameters:
         if ( this%l_cour ) then
            !PRINT*,"Calling courant with this%nR=",this%nR
            call courant(this%nR,this%dtrkc,this%dthkc,this%gsa(threadid)%vrc, &
                 &       this%gsa(threadid)%vtc,this%gsa(threadid)%vpc,        &
                 &       this%gsa(threadid)%brc,this%gsa(threadid)%btc,        &
                 &       this%gsa(threadid)%bpc,nThetaStart,this%sizeThetaB)
         end if

         !--------- Since the fields are given at gridpoints here, this is a good
         !          point for graphical output:
         if ( this%l_graph ) then
#ifdef WITH_MPI
            PERFON('graphout')
            call graphOut_mpi(time,this%nR,this%gsa(threadid)%vrc,           &
                 &            this%gsa(threadid)%vtc,this%gsa(threadid)%vpc, &
                 &            this%gsa(threadid)%brc,this%gsa(threadid)%btc, &
                 &            this%gsa(threadid)%bpc,this%gsa(threadid)%sc,  &
                 &            this%gsa(threadid)%pc,nThetaStart,             &
                 &            this%sizeThetaB,lGraphHeader)
            PERFOFF
#else
            call graphOut(time,this%nR,this%gsa(threadid)%vrc,           &
                 &        this%gsa(threadid)%vtc,this%gsa(threadid)%vpc, &
                 &        this%gsa(threadid)%brc,this%gsa(threadid)%btc, &
                 &        this%gsa(threadid)%bpc,this%gsa(threadid)%sc,  &
                 &        this%gsa(threadid)%pc,nThetaStart,             &
                 &        this%sizeThetaB,lGraphHeader)
#endif
         end if

         !--------- Helicity output:
         if ( this%lHelCalc ) then
            PERFON('hel_out')
            call get_helicity(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,&
                 &        this%gsa(threadid)%vpc,this%gsa(threadid)%cvrc,   &
                 &        this%gsa(threadid)%dvrdtc,                        &
                 &        this%gsa(threadid)%dvrdpc,                        &
                 &        this%gsa(threadid)%dvtdrc,                        &
                 &        this%gsa(threadid)%dvpdrc,HelLMr,Hel2LMr,         &
                 &        HelnaLMr,Helna2LMr,this%nR,nThetaStart)
            PERFOFF
         end if

         !--------- horizontal velocity :

         if ( this%lViscBcCalc ) then
            !write(*,"(2I3,A,3('(',2ES20.12,')'))") nR,nThetaB,&
            !     &" dsdr,dsdt,dsdp = ",&
            !     & SUM(this%gsa(threadid)%drSc),&
            !     & SUM(this%gsa(threadid)%dsdtc),&
            !     & SUM(this%gsa(threadid)%dsdpc)

            call get_nlBLayers(this%gsa(threadid)%vtc,    &
                 &             this%gsa(threadid)%vpc,    &
                 &             this%gsa(threadid)%dvtdrc, &
                 &             this%gsa(threadid)%dvpdrc, &
                 &             this%gsa(threadid)%drSc,   &
                 &             this%gsa(threadid)%dsdtc,  &
                 &             this%gsa(threadid)%dsdpc,  &
                 &             uhLMr,duhLMr,gradsLMr,nR,nThetaStart)
            !write(*,"(2I3,A,3ES20.12)") nR,nThetaB,&
            !     &" uh,duh,grads = ",&
            !     & SUM(uhLMr(:)),SUM(duhLMr(:)),SUM(gradsLMr(:))
         end if


         if ( this%lFluxProfCalc ) then
             call get_fluxes(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,   &
                    &        this%gsa(threadid)%vpc,this%gsa(threadid)%dvrdrc,&
                    &        this%gsa(threadid)%dvtdrc,                       &
                    &        this%gsa(threadid)%dvpdrc,                       &
                    &        this%gsa(threadid)%dvrdtc,                       &
                    &        this%gsa(threadid)%dvrdpc,this%gsa(threadid)%sc, &
                    &        this%gsa(threadid)%pc,this%gsa(threadid)%brc,    &
                    &        this%gsa(threadid)%btc,this%gsa(threadid)%bpc,   &
                    &        this%gsa(threadid)%cbtc,this%gsa(threadid)%cbpc, &
                    &        fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr,nR,   &
                    &        nThetaStart)
         end if

         if ( this%lPerpParCalc ) then
             call get_perpPar(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc, &
                    &         this%gsa(threadid)%vpc,EperpLMr,EparLMr,       &
                    &         EperpaxiLMr,EparaxiLMr,nR,nThetaStart)
         end if


         !--------- Movie output:
         if ( this%l_frame .and. l_movie_oc .and. l_store_frame ) then
            PERFON('mov_out')
            call store_movie_frame(this%nR,this%gsa(threadid)%vrc,                &
                 &                 this%gsa(threadid)%vtc,this%gsa(threadid)%vpc, &
                 &                 this%gsa(threadid)%brc,this%gsa(threadid)%btc, &
                 &                 this%gsa(threadid)%bpc,this%gsa(threadid)%sc,  &
                 &                 this%gsa(threadid)%drSc,                       &
                 &                 this%gsa(threadid)%dvrdpc,                     &
                 &                 this%gsa(threadid)%dvpdrc,                     &
                 &                 this%gsa(threadid)%dvtdrc,                     &
                 &                 this%gsa(threadid)%dvrdtc,                     &
                 &                 this%gsa(threadid)%cvrc,                       &
                 &                 this%gsa(threadid)%cbrc,                       &
                 &                 this%gsa(threadid)%cbtc,nThetaStart,           &
                 &                 this%sizeThetaB,this%leg_helper%bCMB)
            PERFOFF
         end if


         !--------- Stuff for special output:
         !--------- Calculation of magnetic field production and advection terms
         !          for graphic output:
         if ( l_dtB ) then
            PERFON('dtBLM')
            call get_dtBLM(this%nR,this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,&
                 &         this%gsa(threadid)%vpc,this%gsa(threadid)%brc,        &
                 &         this%gsa(threadid)%btc,this%gsa(threadid)%bpc,        &
                 &         nThetaStart,this%sizeThetaB,                          &
                 &         this%dtB_arrays(threadid)%BtVrLM,                     &
                 &         this%dtB_arrays(threadid)%BpVrLM,                     &
                 &         this%dtB_arrays(threadid)%BrVtLM,                     &
                 &         this%dtB_arrays(threadid)%BrVpLM,                     &
                 &         this%dtB_arrays(threadid)%BtVpLM,                     &
                 &         this%dtB_arrays(threadid)%BpVtLM,                     &
                 &         this%dtB_arrays(threadid)%BrVZLM,                     &
                 &         this%dtB_arrays(threadid)%BtVZLM,                     &
                 &         this%dtB_arrays(threadid)%BtVpCotLM,                  &
                 &         this%dtB_arrays(threadid)%BpVtCotLM,                  &
                 &         this%dtB_arrays(threadid)%BtVZcotLM,                  &
                 &         this%dtB_arrays(threadid)%BtVpSn2LM,                  &
                 &         this%dtB_arrays(threadid)%BpVtSn2LM,                  &
                 &         this%dtB_arrays(threadid)%BtVZsn2LM)
            PERFOFF
         end if


         !--------- Torsional oscillation terms:
         PERFON('TO_terms')
         if ( ( this%lTONext .or. this%lTONext2 ) .and. l_mag ) then
            call getTOnext(this%leg_helper%zAS,this%gsa(threadid)%brc,   &
                 &         this%gsa(threadid)%btc,this%gsa(threadid)%bpc,&
                 &         this%lTONext,this%lTONext2,dt,dtLast,this%nR, &
                 &         nThetaStart,this%sizeThetaB,this%BsLast,      &
                 &         this%BpLast,this%BzLast)
         end if

         if ( this%lTOCalc ) then
            call getTO(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,    &
                 &     this%gsa(threadid)%vpc,this%gsa(threadid)%cvrc,   &
                 &     this%gsa(threadid)%dvpdrc,this%gsa(threadid)%brc, &
                 &     this%gsa(threadid)%btc,this%gsa(threadid)%bpc,    &
                 &     this%gsa(threadid)%cbrc,this%gsa(threadid)%cbtc,  &
                 &     this%BsLast,this%BpLast,this%BzLast,              &
                 &     this%TO_arrays(threadid)%dzRstrLM,                &
                 &     this%TO_arrays(threadid)%dzAstrLM,                &
                 &     this%TO_arrays(threadid)%dzCorLM,                 &
                 &     this%TO_arrays(threadid)%dzLFLM,                  &
                 &     dtLast,this%nR,nThetaStart,this%sizeThetaB)
         end if
         PERFOFF

         !if (this%nR == n_r_icb) then
         !   !$OMP ORDERED
         !   write(*,"(2I3,A,I4,F21.17)") nThetaB,threadid,": lorentz_torque_ic = ",&
         !        & exponent(this%lorentz_torque_ic(threadid)),                     &
         !        & fraction(this%lorentz_torque_ic(threadid))
         !   !$OMP END ORDERED
         !end if
      end do ! Loop over theta blocks
      !$OMP end do

      ! do a reduction over the threads by hand
      ! parallelize over the different arrays
      !PERFON('reduc')
      !$OMP SECTIONS private(iThread)
      !$OMP SECTION
      do iThread=1,this%nThreads-1
         this%nl_lm(0)%AdvrLM=this%nl_lm(0)%AdvrLM + this%nl_lm(iThread)%AdvrLM
         this%nl_lm(0)%AdvtLM=this%nl_lm(0)%AdvtLM + this%nl_lm(iThread)%AdvtLM
         this%nl_lm(0)%AdvpLM=this%nl_lm(0)%AdvpLM + this%nl_lm(iThread)%AdvpLM
      end do

      !$OMP SECTION
      do iThread=1,this%nThreads-1
         this%nl_lm(0)%LFrLM=this%nl_lm(0)%LFrLM + this%nl_lm(iThread)%LFrLM
         this%nl_lm(0)%LFtLM=this%nl_lm(0)%LFtLM + this%nl_lm(iThread)%LFtLM
         this%nl_lm(0)%LFpLM=this%nl_lm(0)%LFpLM + this%nl_lm(iThread)%LFpLM
      end do

      !$OMP SECTION
      do iThread=1,this%nThreads-1
         this%nl_lm(0)%VxBrLM=this%nl_lm(0)%VxBrLM + this%nl_lm(iThread)%VxBrLM
         this%nl_lm(0)%VxBtLM=this%nl_lm(0)%VxBtLM + this%nl_lm(iThread)%VxBtLM
         this%nl_lm(0)%VxBpLM=this%nl_lm(0)%VxBpLM + this%nl_lm(iThread)%VxBpLM
      end do

      !$OMP SECTION
      do iThread=1,this%nThreads-1
         this%nl_lm(0)%VSrLM=this%nl_lm(0)%VSrLM + this%nl_lm(iThread)%VSrLM
         this%nl_lm(0)%VStLM=this%nl_lm(0)%VStLM + this%nl_lm(iThread)%VStLM
         this%nl_lm(0)%VSpLM=this%nl_lm(0)%VSpLM + this%nl_lm(iThread)%VSpLM
      end do

      !$OMP SECTION
      do iThread=1,this%nThreads-1
         this%nl_lm(0)%ViscHeatLM=this%nl_lm(0)%ViscHeatLM +  &
                                  this%nl_lm(iThread)%ViscHeatLM
         this%nl_lm(0)%OhmLossLM=this%nl_lm(0)%OhmLossLM +    &
                                 this%nl_lm(iThread)%OhmLossLM
         this%lorentz_torque_ma(0) = this%lorentz_torque_ma(0) + &
                                     this%lorentz_torque_ma(iThread)
      end do

      !$OMP SECTION
      if ( this%lTOCalc ) then
         do iThread=1,this%nThreads-1
            this%TO_arrays(0)%dzRstrLM=this%TO_arrays(0)%dzRstrLM + &
                                       this%TO_arrays(iThread)%dzRstrLM
            this%TO_arrays(0)%dzAstrLM=this%TO_arrays(0)%dzAstrLM + &
                                       this%TO_arrays(iThread)%dzAstrLM
            this%TO_arrays(0)%dzCorLM=this%TO_arrays(0)%dzCorLM + &
                                       this%TO_arrays(iThread)%dzCorLM
            this%TO_arrays(0)%dzLFLM=this%TO_arrays(0)%dzLFLM + &
                                       this%TO_arrays(iThread)%dzLFLM
         end do
      end if

      !$OMP SECTION
      if ( this%lRmsCalc ) then
         do iThread=1,this%nThreads-1
            this%nl_lm(0)%p1LM=this%nl_lm(0)%p1LM+this%nl_lm(iThread)%p1LM
            this%nl_lm(0)%p2LM=this%nl_lm(0)%p2LM+this%nl_lm(iThread)%p2LM
            this%nl_lm(0)%Advt2LM=this%nl_lm(0)%Advt2LM+this%nl_lm(iThread)%Advt2LM
            this%nl_lm(0)%Advp2LM=this%nl_lm(0)%Advp2LM+this%nl_lm(iThread)%Advp2LM
            this%nl_lm(0)%LFt2LM=this%nl_lm(0)%LFt2LM+this%nl_lm(iThread)%LFt2LM
            this%nl_lm(0)%LFp2LM=this%nl_lm(0)%LFp2LM+this%nl_lm(iThread)%LFp2LM
            this%nl_lm(0)%CFt2LM=this%nl_lm(0)%CFt2LM+this%nl_lm(iThread)%CFt2LM
            this%nl_lm(0)%CFp2LM=this%nl_lm(0)%CFp2LM+this%nl_lm(iThread)%CFp2LM
         end do
      end if 

      !$OMP SECTION
      if ( l_dtB ) then
         do iThread=1,this%nThreads-1
            this%dtB_arrays(0)%BtVrLM = this%dtB_arrays(0)%BtVrLM + &
                                        this%dtB_arrays(iThread)%BtVrLM
            this%dtB_arrays(0)%BpVrLM = this%dtB_arrays(0)%BpVrLM + &
                                        this%dtB_arrays(iThread)%BpVrLM
            this%dtB_arrays(0)%BrVtLM = this%dtB_arrays(0)%BrVtLM + &
                                        this%dtB_arrays(iThread)%BrVtLM
            this%dtB_arrays(0)%BrVpLM = this%dtB_arrays(0)%BrVpLM + &
                                        this%dtB_arrays(iThread)%BrVpLM
            this%dtB_arrays(0)%BtVpLM = this%dtB_arrays(0)%BtVpLM + &
                                        this%dtB_arrays(iThread)%BtVpLM
            this%dtB_arrays(0)%BpVtLM = this%dtB_arrays(0)%BpVtLM + &
                                        this%dtB_arrays(iThread)%BpVtLM
            this%dtB_arrays(0)%BrVZLM = this%dtB_arrays(0)%BrVZLM + &
                                        this%dtB_arrays(iThread)%BrVZLM
            this%dtB_arrays(0)%BtVZLM = this%dtB_arrays(0)%BtVZLM + &
                                        this%dtB_arrays(iThread)%BtVZLM
            this%dtB_arrays(0)%BtVpCotLM = this%dtB_arrays(0)%BtVpCotLM + &
                                           this%dtB_arrays(iThread)%BtVpCotLM
            this%dtB_arrays(0)%BpVtCotLM = this%dtB_arrays(0)%BpVtCotLM + &
                                           this%dtB_arrays(iThread)%BpVtCotLM
            this%dtB_arrays(0)%BtVZcotLM = this%dtB_arrays(0)%BtVZcotLM + &
                                           this%dtB_arrays(iThread)%BtVZcotLM
            this%dtB_arrays(0)%BtVpSn2LM = this%dtB_arrays(0)%BtVpSn2LM + &
                                           this%dtB_arrays(iThread)%BtVpSn2LM
            this%dtB_arrays(0)%BpVtSn2LM = this%dtB_arrays(0)%BpVtSn2LM + &
                                           this%dtB_arrays(iThread)%BpVtSn2LM
            this%dtB_arrays(0)%BtVZsn2LM = this%dtB_arrays(0)%BtVZsn2LM + &
                                           this%dtB_arrays(iThread)%BtVZsn2LM
         end do
      end if

!!$    lorentz_torque_ic=0.0_cp
!!$    c=0.0_cp
!!$    do iThread=0,this%nThreads-1
!!$       y=this%lorentz_torque_ic(iThread)-c
!!$       t=lorentz_torque_ic+y
!!$       c=(t-lorentz_torque_ic)-y
!!$       lorentz_torque_ic=t
!!$    end do

      !$OMP END SECTIONS
      !PERFOFF

      !$OMP END PARALLEL
      lorentz_torque_ic = lorentz_torques_ic(1)
      do nThetaB=2,this%nThetaBs
         lorentz_torque_ic = lorentz_torque_ic + lorentz_torques_ic(nThetaB)
      end do
      !lorentz_torque_ic = this%lorentz_torque_ic(0)
      lorentz_torque_ma = this%lorentz_torque_ma(0)

      !if (this%nR == n_r_icb) then
      !   write(*,"(A,2ES20.12)") "after OMP PARALLEL, lorentz_torque = ",&
      !        & lorentz_torque_ic,lorentz_torque_ma
      !end if

      if (DEBUG_OUTPUT) then
         call this%nl_lm(0)%output()
      end if

      !-- Partial calculation of time derivatives (horizontal parts):
      !   input flm...  is in (l,m) space at radial grid points this%nR !
      !   Only dVxBh needed for boundaries !
      !   get_td finally calculates the d*dt terms needed for the
      !   time step performed in s_LMLoop.f . This should be distributed
      !   over the different models that s_LMLoop.f parallelizes over.
      !write(*,"(A,I4,2ES20.13)") "before_td: ", &
      !     &  this%nR,sum(real(conjg(VxBtLM)*VxBtLM)),sum(real(conjg(VxBpLM)*VxBpLM))
      !PERFON('get_td')
      call this%nl_lm(0)%get_td(this%nR,this%nBc,this%lRmsCalc,         &
           &                    dVSrLM,dVXirLM,dVxBhLM,dwdt,dzdt,dpdt,  &
           &                    dsdt,dxidt,dbdt,djdt,this%leg_helper)

      !PERFOFF
      !write(*,"(A,I4,ES20.13)") "after_td:  ", &
      !     & this%nR,sum(real(conjg(dVxBhLM(:,this%nR_Mag))*dVxBhLM(:,this%nR_Mag)))
      !-- Finish calculation of TO variables:
      if ( this%lTOcalc ) then
         call getTOfinish(this%nR,dtLast,this%leg_helper%zAS,             &
              &           this%leg_helper%dzAS,this%leg_helper%ddzAS,     &
              &           this%TO_arrays(0)%dzRstrLM,this%TO_arrays(0)%dzAstrLM,&
              &           this%TO_arrays(0)%dzCorLM,this%TO_arrays(0)%dzLFLM)
      end if

      !--- Form partial horizontal derivaties of magnetic production and
      !    advection terms:
      if ( l_dtB ) then
         PERFON('dtBLM')
         call get_dH_dtBLM(this%nR,this%dtB_arrays(0)%BtVrLM,                    &
              &            this%dtB_arrays(0)%BpVrLM,                            &
              &            this%dtB_arrays(0)%BrVtLM,this%dtB_arrays(0)%BrVpLM,  &
              &            this%dtB_arrays(0)%BtVpLM,this%dtB_arrays(0)%BpVtLM,  &
              &            this%dtB_arrays(0)%BrVZLM,this%dtB_arrays(0)%BtVZLM,  &
              &            this%dtB_arrays(0)%BtVpCotLM,                         &
              &            this%dtB_arrays(0)%BpVtCotLM,                         &
              &            this%dtB_arrays(0)%BtVpSn2LM,                         &
              &            this%dtB_arrays(0)%BpVtSn2LM)
         PERFOFF
      end if
    end subroutine do_iteration_ThetaBlocking_OpenMP
!-------------------------------------------------------------------------------
end module rIterThetaBlocking_OpenMP_mod
