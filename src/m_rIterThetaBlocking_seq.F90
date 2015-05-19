#include "perflib_preproc.cpp"
MODULE rIterThetaBlocking_seq_mod
  USE rIterThetaBlocking_mod,only: rIterThetaBlocking_t
  USE grid_space_arrays_mod,ONLY: grid_space_arrays_t
  USE nonlinear_lm_mod,only: nonlinear_lm_t

  USE truncation, ONLY: lm_max,lmP_max,nrp,ncp,l_max,lmP_max_dtB,&
       &n_phi_maxStr,n_theta_maxStr,n_r_maxStr
  USE blocking, only: nfs
  USE logic, ONLY: l_mag,l_conv,l_mag_kin,l_heat,l_ht,l_anel,l_mag_LF,&
       & l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic, l_cond_ic, &
       & l_rot_ma, l_cond_ma, l_dtB, l_store_frame, l_movie_oc
  USE radial_data,ONLY: n_r_cmb, n_r_icb
  USE radial_functions, ONLY: or2, orho1
  USE output_data, only: ngform

  IMPLICIT NONE

  TYPE,PUBLIC,EXTENDS(rIterThetaBlocking_t) :: rIterThetaBlocking_seq_t
     type(grid_space_arrays_t) :: gsa
     TYPE(nonlinear_lm_t) :: nl_lm
   CONTAINS
     procedure :: initialize => initialize_rIterThetaBlocking_seq
     procedure :: finalize => finalize_rIterThetaBlocking_seq
     procedure :: do_iteration => do_iteration_ThetaBlocking_seq
     procedure :: getType => getThisType
  END TYPE rIterThetaBlocking_seq_t
CONTAINS
  FUNCTION getThisType(this)
    class(rIterThetaBlocking_seq_t) :: this
    character(len=100) :: getThisType
    getThisType="rIterThetaBlocking_seq_t"
  END FUNCTION getThisType

  SUBROUTINE initialize_rIterThetaBlocking_seq(this)
    class(rIterThetaBlocking_seq_t) :: this

    CALL this%allocate_common_arrays()
    call this%gsa%initialize()
    call this%nl_lm%initialize(lmP_max)
  END SUBROUTINE initialize_rIterThetaBlocking_seq

  SUBROUTINE finalize_rIterThetaBlocking_seq(this)
    class(rIterThetaBlocking_seq_t) :: this

    CALL this%deallocate_common_arrays()
    call this%gsa%finalize()
    call this%nl_lm%finalize()
  END SUBROUTINE finalize_rIterThetaBlocking_seq

  SUBROUTINE do_iteration_ThetaBlocking_seq(this,nR,nBc,time,dt,dtLast,&
       &                 dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM, &
       &                 br_vt_lm_cmb,br_vp_lm_cmb,                    &
       &                 br_vt_lm_icb,br_vp_lm_icb,                    &
       &                 lorentz_torque_ic, lorentz_torque_ma,         &
       &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,      &
       &                 duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,    &
       &                 fpoynLMr,fresLMr,EperpLMr,EparLMr,EperpaxiLMr,&
       &                 EparaxiLmr)
    CLASS(rIterThetaBlocking_seq_t) :: this
    INTEGER, INTENT(IN) :: nR,nBc
    REAL(kind=8),INTENT(IN) :: time,dt,dtLast

    COMPLEX(kind=8),INTENT(OUT),DIMENSION(:) :: dwdt,dzdt,dpdt,dsdt,dVSrLM
    COMPLEX(kind=8),INTENT(OUT),DIMENSION(:) :: dbdt,djdt,dVxBhLM
    !---- Output of nonlinear products for nonlinear
    !     magnetic boundary conditions (needed in s_updateB.f):
    COMPLEX(kind=8),INTENT(OUT) :: br_vt_lm_cmb(:) ! product br*vt at CMB
    COMPLEX(kind=8),INTENT(OUT) :: br_vp_lm_cmb(:) ! product br*vp at CMB
    COMPLEX(kind=8),INTENT(OUT) :: br_vt_lm_icb(:) ! product br*vt at ICB
    COMPLEX(kind=8),INTENT(OUT) :: br_vp_lm_icb(:) ! product br*vp at ICB
    REAL(kind=8),intent(OUT) :: lorentz_torque_ma,lorentz_torque_ic
    REAL(kind=8),INTENT(OUT),DIMENSION(:) :: HelLMr,Hel2LMr,HelnaLMr,Helna2LMr
    REAL(kind=8),INTENT(OUT),DIMENSION(:) :: uhLMr,duhLMr,gradsLMr
    REAL(kind=8),INTENT(OUT),DIMENSION(:) :: fconvLMr,fkinLMr,fviscLMr
    REAL(kind=8),INTENT(OUT),DIMENSION(:) :: fpoynLMr,fresLMr
    REAL(kind=8),INTENT(OUT),DIMENSION(:) :: EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr


    INTEGER :: l,lm,nThetaB,nThetaLast,nThetaStart,nThetaStop

    LOGICAL :: DEBUG_OUTPUT=.FALSE.

    this%nR=nR
    this%nBc=nBc
    this%isRadialBoundaryPoint=(nR.EQ.n_r_cmb).OR.(nR.EQ.n_r_icb)

    IF ( this%l_cour ) THEN
       this%dtrkc=1.D10
       this%dthkc=1.D10
    END IF
    IF ( this%lTOCalc ) THEN
       !------ Zero lm coeffs for first theta block:
       DO l=0,l_max
          this%TO_arrays%dzRstrLM(l+1)=0.D0
          this%TO_arrays%dzAstrLM(l+1)=0.D0
          this%TO_arrays%dzCorLM(l+1) =0.D0
          this%TO_arrays%dzLFLM(l+1)  =0.D0
       END DO
    END IF

    !----- Prepare legendre transform:
    !      legPrepG collects all the different modes necessary 
    !      to calculate the non-linear terms at a radial grid point nR
    PERFON('legPrepG')
    IF (DEBUG_OUTPUT) THEN
       WRITE(*,"(I3,A,I1,2(A,L1))") this%nR,": nBc = ",this%nBc,", lDeriv = ",this%lDeriv,&
            &", l_mag = ",l_mag
    END IF
    CALL legPrepG(this%nR,this%nBc,this%lDeriv,this%lRmsCalc,this%l_frame, &
         &        this%lTOnext,this%lTOnext2,this%lTOcalc,                 &
         &        this%leg_helper)
    PERFOFF

    !IF (DEBUG_OUTPUT) THEN
    !   WRITE(*,"(I3,A,44ES20.12)") this%nR,": legPrepG results = ",&
    !        & SUM(dLhw),SUM(dLhdw),SUM(dLhz),SUM(vhG),SUM(vhC),&
    !        & SUM(dvhdrG),SUM(dvhdrC),&
    !        & SUM(dLhb),SUM(dLhj),SUM(bhG),SUM(bhC),SUM(cbhG),SUM(cbhC),&
    !        & SUM(sR),SUM(dsR),SUM(preR),SUM(dpR),SUM(zAS), SUM(dzAS),&
    !        &SUM(ddzAS),SUM(bCMB),omegaIC,omegaMA

    !   WRITE(*,"(A,I4,A)") "We have ",this%nThetaBs," theta blocks!"
    !END IF
    !----- Blocking of loops over ic (theta):
    DO nThetaB=1,this%nThetaBs
       nThetaLast =(nThetaB-1) * this%sizeThetaB
       nThetaStart=nThetaLast+1
       nThetaStop =nThetaLast + this%sizeThetaB
       !WRITE(*,"(I3,A,I4,A,I4)") nThetaB,". theta block from ",nThetaStart," to ",nThetaStop

       CALL this%transform_to_grid_space(nThetaStart,nThetaStop,this%gsa)

       !--------- Calculation of nonlinear products in grid space:
       IF ( (.NOT.this%isRadialBoundaryPoint) .OR. this%lMagNlBc ) THEN 
          !write(*,"(I4,A,ES20.13)") this%nR,", vp = ",sum(real(conjg(vpc)*vpc))
          PERFON('get_nl')
          CALL get_nl(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,&
               &      this%gsa%dvrdrc,this%gsa%dvtdrc,this%gsa%dvpdrc,this%gsa%cvrc,  &
               &      this%gsa%dvrdtc,this%gsa%dvrdpc,this%gsa%dvtdpc,this%gsa%dvpdpc,    &
               &      this%gsa%brc,this%gsa%btc,this%gsa%bpc,this%gsa%cbrc,&
               &      this%gsa%cbtc,this%gsa%cbpc,this%gsa%sc,  &
               &      this%gsa%Advr,this%gsa%Advt,this%gsa%Advp,&
               &      this%gsa%LFr,this%gsa%LFt,this%gsa%LFp,     &
               &      this%gsa%VSr,this%gsa%VSt,this%gsa%VSp,this%gsa%VxBr,&
               &      this%gsa%VxBt,this%gsa%VxBp,     &
               &      this%gsa%ViscHeat,this%gsa%OhmLoss,               &
               &      this%nR,this%nBc,nThetaStart)
          PERFOFF

          CALL this%transform_to_lm_space(nThetaStart,nThetaStop,this%gsa,this%nl_lm)

       ELSE IF ( l_mag ) THEN
          DO lm=1,lmP_max
             this%nl_lm%VxBtLM(lm)=0.D0
             this%nl_lm%VxBpLM(lm)=0.D0
          END DO
       END IF

       !---- Calculation of nonlinear products needed for conducting mantle or
       !     conducting inner core if free stress BCs are applied:
       !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
       !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
       !     and br_vp_lm_cmb in lm-space, respectively the contribution
       !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
       !     These products are used in get_b_nl_bcs.
       PERFON('nl_cmb')
       IF ( this%nR.EQ.n_r_cmb .AND. l_b_nl_cmb ) THEN
          CALL get_br_v_bcs(this%gsa%brc,this%gsa%vtc,this%gsa%vpc,this%leg_helper%omegaMA,              &
               &            or2(this%nR),orho1(this%nR),nThetaStart,this%sizeThetaB,    &
               &            br_vt_lm_cmb,br_vp_lm_cmb)
       ELSE IF ( this%nR.EQ.n_r_icb .AND. l_b_nl_icb ) THEN
          CALL get_br_v_bcs(this%gsa%brc,this%gsa%vtc,this%gsa%vpc,this%leg_helper%omegaIC,              &
               &            or2(this%nR),orho1(this%nR),nThetaStart,this%sizeThetaB,    &
               &            br_vt_lm_icb,br_vp_lm_icb)
       END IF
       PERFOFF
       !--------- Calculate Lorentz torque on inner core:
       !          each call adds the contribution of the theta-block to
       !          lorentz_torque_ic
       PERFON('lorentz')
       IF ( this%nR.EQ.n_r_icb .AND. l_mag_LF .AND. l_rot_ic .AND. l_cond_ic  ) &
            & CALL get_lorentz_torque(lorentz_torque_ic,         &
            &                         nThetaStart,this%sizeThetaB,         &
            &                         this%gsa%brc,this%gsa%bpc,this%nR)

       !--------- Calculate Lorentz torque on mantle:
       !          note: this calculates a torque of a wrong sign.
       !          sign is reversed at the end of the theta blocking.
       IF ( this%nR.EQ.n_r_cmb .AND. l_mag_LF .AND. l_rot_ma .AND. l_cond_ma ) &
            & CALL get_lorentz_torque(lorentz_torque_ma,          &
            &                         nThetaStart,this%sizeThetaB,          &
            &                         this%gsa%brc,this%gsa%bpc,this%nR)

       PERFOFF
       !--------- Calculate courant condition parameters:
       IF ( this%l_cour ) THEN
          !PRINT*,"Calling courant with this%nR=",this%nR
          CALL courant(this%nR,this%dtrkc,this%dthkc,   &
               &       this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%brc,this%gsa%btc,this%gsa%bpc,  &
               &       nThetaStart,this%sizeThetaB)
       END IF

       !--------- Since the fields are given at gridpoints here, this is a good
       !          point for graphical output:
       IF ( this%l_graph ) THEN
          PERFON('graphout')
          CALL graphOut_mpi(time,this%nR,ngform,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc, &
               &        this%gsa%brc,this%gsa%btc,this%gsa%bpc,this%gsa%sc,&
               &        nThetaStart,this%sizeThetaB,.FALSE.)
          PERFOFF
       END IF

       !--------- Helicity output:
       IF ( this%lHelCalc ) THEN
          PERFON('hel_out')
          CALL getHelLM(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,                        &
               &        this%gsa%cvrc,this%gsa%dvrdtc,this%gsa%dvrdpc,this%gsa%dvtdrc,this%gsa%dvpdrc,   &
               &        HelLMr,Hel2LMr,         &
               &        HelnaLMr,Helna2LMr,     &
               &        this%nR,nThetaStart)
          PERFOFF
       END IF

       !--------- horizontal velocity :

       IF ( this%lViscBcCalc ) THEN
          CALL get_nlBLayers(this%gsa%vtc,this%gsa%vpc,this%gsa%dvtdrc,this%gsa%dvpdrc,    &
               &             this%gsa%drSc,this%gsa%dsdtc,this%gsa%dsdpc,    &
               &             uhLMr,duhLMr,gradsLMr,nR,nThetaStart)
       END IF

       IF ( this%lFluxProfCalc ) THEN
           CALL get_fluxes(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%dvrdrc,  &
                  &        this%gsa%dvtdrc,this%gsa%dvpdrc,this%gsa%dvrdtc,         &
                  &        this%gsa%dvrdpc,this%gsa%sc,this%gsa%pc,this%gsa%brc,    &
                  &        this%gsa%btc,this%gsa%bpc,this%gsa%cbtc,this%gsa%cbpc,   &
                  &        fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr,nR,nThetaStart)
       END IF

       IF ( this%lPerpParCalc ) THEN
           CALL get_perpPar(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,  &
                  &        EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr,nR,nThetaStart)
       END IF

       !--------- Movie output:
       IF ( this%l_frame .AND. l_movie_oc .AND. l_store_frame ) THEN
          PERFON('mov_out')
          CALL store_movie_frame(this%nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,                   &
               &                 this%gsa%brc,this%gsa%btc,this%gsa%bpc,this%gsa%sc,this%gsa%drSc, &
               &                 this%gsa%dvrdpc,this%gsa%dvpdrc,this%gsa%dvtdrc,&
               &                 this%gsa%dvrdtc,this%gsa%cvrc, &
               &                 this%gsa%cbrc,this%gsa%cbtc,nThetaStart,this%sizeThetaB,&
               &                 this%leg_helper%bCMB)
          PERFOFF
       END IF


       !--------- Stuff for special output:
       !--------- Calculation of magnetic field production and advection terms
       !          for graphic output:
       IF ( l_dtB ) THEN
          PERFON('dtBLM')
          CALL get_dtBLM(this%nR,this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,&
               &         this%gsa%brc,this%gsa%btc,this%gsa%bpc,        &
               &         nThetaStart,this%sizeThetaB,                     &
               &         this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM,this%dtB_arrays%BrVtLM,&
               &         this%dtB_arrays%BrVpLM,this%dtB_arrays%BtVpLM,this%dtB_arrays%BpVtLM,  &
               &         this%dtB_arrays%BrVZLM,this%dtB_arrays%BtVZLM,this%dtB_arrays%BtVpCotLM,&
               &         this%dtB_arrays%BpVtCotLM,this%dtB_arrays%BtVZcotLM,&
               &         this%dtB_arrays%BtVpSn2LM,this%dtB_arrays%BpVtSn2LM,this%dtB_arrays%BtVZsn2LM)
          PERFOFF
       END IF


       !--------- Torsional oscillation terms:
       PERFON('TO_terms')
       IF ( ( this%lTONext .OR. this%lTONext2 ) .AND. l_mag ) THEN
          CALL getTOnext(this%leg_helper%zAS,this%gsa%brc,this%gsa%btc,this%gsa%bpc,this%lTONext,     &
               &         this%lTONext2,dt,dtLast,this%nR,nThetaStart,this%sizeThetaB, &
               &         this%BsLast,this%BpLast,this%BzLast)
       END IF

       IF ( this%lTOCalc ) THEN
          CALL getTO(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc,this%gsa%dvpdrc,         &
               &     this%gsa%brc,this%gsa%btc,this%gsa%bpc,this%gsa%cbrc,this%gsa%cbtc,         &
               &     this%BsLast,this%BpLast,this%BzLast,         &
               &     this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM,&
               &     this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM,         &
               &     dtLast,this%nR,nThetaStart,this%sizeThetaB)
       END IF
       PERFOFF

    END DO ! Loop over theta blocks


    !-- Partial calculation of time derivatives (horizontal parts):
    !   input flm...  is in (l,m) space at radial grid points this%nR !
    !   Only dVxBh needed for boundaries !
    !   get_td finally calculates the d*dt terms needed for the 
    !   time step performed in s_LMLoop.f . This should be distributed
    !   over the different models that s_LMLoop.f parallelizes over. 
    !WRITE(*,"(A,I4,4ES20.13)") "before_td: ",this%nR,SUM(this%nl_lm%VxBtLM),&
    !     & SUM(this%nl_lm%VxBpLM)
    PERFON('get_td')
    CALL get_td(this%nR,this%nBc,this%lRmsCalc,dVSrLM,dVxBhLM,   &
         &      dwdt,dzdt,dpdt,   &
         &      dsdt,dbdt,djdt,   &
         &      this%nl_lm, this%leg_helper)
    PERFOFF
    !DO lm=1,lm_max
    !   WRITE(*,"(2(I3,A),2ES20.12)") this%nR,": dwdt(",lm,") = ",dwdt(lm)
    !END DO
    !WRITE(*,"(A,I4,4ES20.13)") "after_td: dwdt ",this%nR, SUM(dwdt)
    !-- Finish calculation of TO variables:
    IF ( this%lTOcalc ) THEN                                   
       CALL getTOfinish(this%nR,dtLast,this%leg_helper%zAS,this%leg_helper%dzAS,this%leg_helper%ddzAS, &
            &           this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM,&
            &           this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM)
    END IF

    !--- Form partial horizontal derivaties of magnetic production and
    !    advection terms:
    IF ( l_dtB ) THEN
       CALL get_dH_dtBLM(this%nR,this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM,&
            &            this%dtB_arrays%BrVtLM,                      &
            &            this%dtB_arrays%BrVpLM,this%dtB_arrays%BtVpLM,&
            &            this%dtB_arrays%BpVtLM,this%dtB_arrays%BrVZLM,&
            &            this%dtB_arrays%BtVZLM,this%dtB_arrays%BtVpCotLM, &
            &            this%dtB_arrays%BpVtCotLM,this%dtB_arrays%BtVZcotLM,&
            &            this%dtB_arrays%BtVpSn2LM,                &
            &            this%dtB_arrays%BpVtSn2LM,this%dtB_arrays%BtVZsn2LM)
    END IF
  END SUBROUTINE do_iteration_ThetaBlocking_seq

END MODULE rIterThetaBlocking_seq_mod
