#include "perflib_preproc.cpp"
MODULE rIterThetaBlocking_OpenMP_mod
#ifdef WITHOMP
  use omp_lib
#endif
  USE rIterThetaBlocking_mod, only: rIterThetaBlocking_t

  USE truncation, ONLY: lm_max,lmP_max,nrp,ncp,l_max,lmP_max_dtB,&
       &n_phi_maxStr,n_theta_maxStr,n_r_maxStr
  USE blocking, only: nfs
  USE logic, ONLY: l_mag,l_conv,l_mag_kin,l_heat,l_ht,l_anel,l_mag_LF,&
       & l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic, l_cond_ic, &
       & l_rot_ma, l_cond_ma, l_viscBcCalc, l_dtB, l_store_frame, l_movie_oc, &
       & l_fluxProfs
  USE radial_data,ONLY: n_r_cmb, n_r_icb
  USE radial_functions, ONLY: or2, orho1
  USE output_data, only: ngform
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
  USE legendre_trafo,only: legTFG
  USE leg_helper_mod, only: leg_helper_t
  USE nonlinear_lm_mod,only:nonlinear_lm_t
  USE grid_space_arrays_mod,only: grid_space_arrays_t
  IMPLICIT NONE


  TYPE,PUBLIC,EXTENDS(rIterThetaBlocking_t) :: rIterThetaBlocking_OpenMP_t
     integer :: nThreads
     TYPE(grid_space_arrays_t),ALLOCATABLE,DIMENSION(:) :: gsa
     TYPE(nonlinear_lm_t),ALLOCATABLE,dimension(:) :: nl_lm
     REAL(kind=8),DIMENSION(:), ALLOCATABLE :: lorentz_torque_ic,lorentz_torque_ma
   CONTAINS
     procedure :: initialize => initialize_rIterThetaBlocking_OpenMP
     procedure :: finalize => finalize_rIterThetaBlocking_OpenMP
     procedure :: do_iteration => do_iteration_ThetaBlocking_OpenMP
     procedure :: getType => getThisType
  END TYPE rIterThetaBlocking_OpenMP_t

CONTAINS
  FUNCTION getThisType(this)
    CLASS(rIterThetaBlocking_OpenMP_t) :: this
    character(len=100) :: getThisType
    getThisType="rIterThetaBlocking_OpenMP_t"
  END FUNCTION getThisType

  SUBROUTINE initialize_rIterThetaBlocking_OpenMP(this)
    class(rIterThetaBlocking_OpenMP_t) :: this
    integer :: threadid

#ifdef WITHOMP
    this%nThreads=omp_get_max_threads()
#else
    this%nThreads=1
#endif
    ALLOCATE(this%gsa(0:this%nThreads-1))
    ALLOCATE(this%nl_lm(0:this%nThreads-1))
    ALLOCATE(this%lorentz_torque_ic(0:this%nThreads-1))
    ALLOCATE(this%lorentz_torque_ma(0:this%nThreads-1))

    CALL this%allocate_common_arrays()
    !$OMP PARALLEL default(shared) shared(this,lmP_max) private(threadid)
#ifdef WITHOMP
    threadid = omp_get_thread_num()
#else
    threadid = 0
#endif
    call this%gsa(threadid)%initialize()
    call this%nl_lm(threadid)%initialize(lmP_max)
    !$OMP END PARALLEL
  END SUBROUTINE initialize_rIterThetaBlocking_OpenMP

  SUBROUTINE finalize_rIterThetaBlocking_OpenMP(this)
    class(rIterThetaBlocking_OpenMP_t) :: this
    integer :: threadid

    CALL this%deallocate_common_arrays()
    !$OMP PARALLEL default(shared) shared(this) private(threadid)
#ifdef WITHOMP
    threadid = omp_get_thread_num()
#else
    threadid = 0
#endif
    call this%gsa(threadid)%finalize()
    call this%nl_lm(threadid)%finalize()
    !$OMP END PARALLEL
    DEALLOCATE(this%gsa)
    DEALLOCATE(this%nl_lm)
    DEALLOCATE(this%lorentz_torque_ic)
    DEALLOCATE(this%lorentz_torque_ma)
  END SUBROUTINE finalize_rIterThetaBlocking_OpenMP

  SUBROUTINE do_iteration_ThetaBlocking_OpenMP(this,nR,nBc,time,dt,dtLast,&
       &                 dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM, &
       &                 br_vt_lm_cmb,br_vp_lm_cmb,   &
       &                 br_vt_lm_icb,br_vp_lm_icb,&
       &                 lorentz_torque_ic, lorentz_torque_ma,&
       &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,&
       &                 gradsLMr,fconvLMr,fkinLMr,fviscLMr,    &
       &                 fpoynLMr,fresLMr)
    CLASS(rIterThetaBlocking_OpenMP_t) :: this
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

    INTEGER :: l,lm,nThetaB,nThetaLast,nThetaStart,nThetaStop
    !INTEGER :: nTheta,nPhi
    !INTEGER :: nR_Mag
    INTEGER :: threadid,iThread
    LOGICAL :: DEBUG_OUTPUT=.false.
    REAL(kind=8) :: lt,y,c,t,lorentz_torques_ic(this%nThetaBs)

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
    !PERFON('legPrepG')
    IF (DEBUG_OUTPUT) THEN
       WRITE(*,"(I3,A,I1,2(A,L1))") this%nR,": nBc = ",this%nBc,", lDeriv = ",this%lDeriv,&
            &", l_mag = ",l_mag
    END IF
    CALL legPrepG(this%nR,this%nBc,this%lDeriv,this%lRmsCalc,this%l_frame, &
         &        this%lTOnext,this%lTOnext2,this%lTOcalc,                 &
         &        this%leg_helper)
    !PERFOFF

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
    !$OMP PARALLEL default(shared) &
    !$OMP SHARED(this,l_mag,l_b_nl_cmb,l_b_nl_icb,l_mag_LF,l_rot_ic,l_cond_ic) &
    !$OMP SHARED(l_rot_ma,l_cond_ma,l_viscBcCalc,l_movie_oc,l_store_frame,l_dtB) &
    !$OMP SHARED(lmP_max,n_r_cmb,n_r_icb) &
    !$OMP SHARED(or2,orho1,time,ngform,dt,dtLast,DEBUG_OUTPUT) &
    !$OMP PRIVATE(threadid,nThetaLast,nThetaStart,nThetaStop,y,t,c,lt) &
    !$OMP shared(br_vt_lm_cmb,br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb) &
    !$OMP SHARED(lorentz_torques_ic) &
    !$OMP shared(HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,gradsLMr) &
    !$OMP shared(fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr)
#ifdef WITHOMP
    threadid = omp_get_thread_num()
#else
    threadid = 0
#endif
    this%lorentz_torque_ma(threadid) = 0.0D0
    this%lorentz_torque_ic(threadid) = 0.0D0
    c = 0.0D0
    !$OMP SINGLE
    br_vt_lm_cmb=CMPLX(0.0,0.0,kind=kind(br_vt_lm_cmb))
    br_vp_lm_cmb=CMPLX(0.0,0.0,kind=kind(br_vp_lm_cmb))
    br_vt_lm_icb=CMPLX(0.0,0.0,kind=kind(br_vt_lm_icb))
    br_vp_lm_icb=CMPLX(0.0,0.0,kind=kind(br_vp_lm_icb))
    HelLMr=0.0D0
    Hel2LMr=0.0D0
    HelnaLMr=0.0D0
    Helna2LMr=0.0D0
    uhLMr = 0.0D0
    duhLMr = 0.0D0
    gradsLMr = 0.0D0
    fconvLMr=0.D0
    fkinLMr=0.D0
    fviscLMr=0.D0
    fpoynLMr=0.D0
    fresLMr=0.D0
    !$OMP END SINGLE
    !$OMP BARRIER
    call this%nl_lm(threadid)%set_zero()
    !$OMP DO &
    !$OMP reduction(+:br_vt_lm_cmb,br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb) &
    !$OMP reduction(+:HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,gradsLMr) &
    !$OMP reduction(+:fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr)

    DO nThetaB=1,this%nThetaBs
       nThetaLast =(nThetaB-1) * this%sizeThetaB
       nThetaStart=nThetaLast+1
       nThetaStop =nThetaLast + this%sizeThetaB
       !WRITE(*,"(I3,A,I4,A,I4)") nThetaB,". theta block from ",nThetaStart," to ",nThetaStop

       !PERFON('lm2grid')
       CALL this%transform_to_grid_space(nThetaStart,nThetaStop,&
            &                            this%gsa(threadid))
       !PERFOFF

       !--------- Calculation of nonlinear products in grid space:
       IF ( (.NOT.this%isRadialBoundaryPoint) .OR. this%lMagNlBc ) THEN 
          
          !IF (DEBUG_OUTPUT) THEN
             !IF (this%nR.EQ.2) THEN
             !   WRITE(*,"(A,I2,A,I2)") "++++ START gsa(",threadid,") for nThetaB = ",nThetaB
             !   CALL this%gsa(threadid)%output_nl_input()
             !   WRITE(*,"(A,I2,A,I2)") "---- END   gsa(",threadid,") for nThetaB = ",nThetaB
             !END IF
          !END IF

          !PERFON('get_nl')
          CALL get_nl(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,this%gsa(threadid)%vpc,&
               &      this%gsa(threadid)%dvrdrc,this%gsa(threadid)%dvtdrc,&
               &      this%gsa(threadid)%dvpdrc,this%gsa(threadid)%cvrc,  &
               &      this%gsa(threadid)%dvrdtc,this%gsa(threadid)%dvrdpc,&
               &      this%gsa(threadid)%dvtdpc,this%gsa(threadid)%dvpdpc,    &
               &      this%gsa(threadid)%brc,this%gsa(threadid)%btc,this%gsa(threadid)%bpc,&
               &      this%gsa(threadid)%cbrc,this%gsa(threadid)%cbtc,&
               &      this%gsa(threadid)%cbpc,this%gsa(threadid)%sc,  &
               &      this%gsa(threadid)%Advr,this%gsa(threadid)%Advt,&
               &      this%gsa(threadid)%Advp,this%gsa(threadid)%LFr,&
               &      this%gsa(threadid)%LFt,this%gsa(threadid)%LFp,     &
               &      this%gsa(threadid)%VSr,this%gsa(threadid)%VSt,&
               &      this%gsa(threadid)%VSp,this%gsa(threadid)%VxBr,&
               &      this%gsa(threadid)%VxBt,this%gsa(threadid)%VxBp,     &
               &      this%gsa(threadid)%ViscHeat,this%gsa(threadid)%OhmLoss,  &
               &      this%nR,this%nBc,nThetaStart)
          !PERFOFF
          !IF (DEBUG_OUTPUT) THEN
          !   IF (this%nR.eq.2) THEN
          !      WRITE(*,"(A,I2,A,I2)") "++++ START gsa(",threadid,") for nThetaB = ",nThetaB
          !      CALL this%gsa(threadid)%output()
          !      WRITE(*,"(A,I2,A,I2)") "---- END   gsa(",threadid,") for nThetaB = ",nThetaB
          !   END IF
          !END IF
          !PERFON('grid2lm')
          CALL this%transform_to_lm_space(nThetaStart,nThetaStop,this%gsa(threadid),this%nl_lm(threadid))
          !PERFOFF

       ELSE IF ( l_mag ) THEN
          DO lm=1,lmP_max
             this%nl_lm(threadid)%VxBtLM(lm)=0.D0
             this%nl_lm(threadid)%VxBpLM(lm)=0.D0
          END DO
       END IF

       !---- Calculation of nonlinear products needed for conducting mantle or
       !     conducting inner core if free stress BCs are applied:
       !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
       !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
       !     and br_vp_lm_cmb in lm-space, respectively the contribution
       !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
       !     These products are used in get_b_nl_bcs.
       !PERFON('nl_cmb')
       IF ( this%nR.EQ.n_r_cmb .AND. l_b_nl_cmb ) THEN
          CALL get_br_v_bcs(this%gsa(threadid)%brc,this%gsa(threadid)%vtc,&
               &            this%gsa(threadid)%vpc,this%leg_helper%omegaMA,              &
               &            or2(this%nR),orho1(this%nR),nThetaStart,this%sizeThetaB,    &
               &            br_vt_lm_cmb,br_vp_lm_cmb)
       ELSE IF ( this%nR.EQ.n_r_icb .AND. l_b_nl_icb ) THEN
          CALL get_br_v_bcs(this%gsa(threadid)%brc,this%gsa(threadid)%vtc,&
               &this%gsa(threadid)%vpc,this%leg_helper%omegaIC,              &
               &            or2(this%nR),orho1(this%nR),nThetaStart,this%sizeThetaB,    &
               &            br_vt_lm_icb,br_vp_lm_icb)
       END IF
       !PERFOFF
       !--------- Calculate Lorentz torque on inner core:
       !          each call adds the contribution of the theta-block to
       !          lorentz_torque_ic
       !PERFON('lorentz')
       IF ( this%nR.EQ.n_r_icb .AND. l_mag_LF .AND. l_rot_ic .AND. l_cond_ic  ) THEN
          lorentz_torques_ic(nThetaB)=0.0D0
          CALL get_lorentz_torque(lorentz_torques_ic(nThetaB),     &
               &                  nThetaStart,this%sizeThetaB,         &
               &                  this%gsa(threadid)%brc,this%gsa(threadid)%bpc,this%nR)
          !y=lt-c
          !t=this%lorentz_torque_ic(threadid)+y
          !c=(t-this%lorentz_torque_ic(threadid))-y
          !this%lorentz_torque_ic(threadid)=t
       END IF

       !--------- Calculate Lorentz torque on mantle:
       !          note: this calculates a torque of a wrong sign.
       !          sign is reversed at the end of the theta blocking.
       IF ( this%nR.EQ.n_r_cmb .AND. l_mag_LF .AND. l_rot_ma .AND. l_cond_ma ) THEN
          CALL get_lorentz_torque(this%lorentz_torque_ma(threadid),          &
               &                  nThetaStart,this%sizeThetaB,          &
               &                  this%gsa(threadid)%brc,this%gsa(threadid)%bpc,this%nR)
       END IF
       !PERFOFF

       !--------- Calculate courant condition parameters:
       IF ( this%l_cour ) THEN
          !PRINT*,"Calling courant with this%nR=",this%nR
          CALL courant(this%nR,this%dtrkc,this%dthkc,   &
               &       this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,this%gsa(threadid)%vpc,&
               &       this%gsa(threadid)%brc,this%gsa(threadid)%btc,this%gsa(threadid)%bpc,  &
               &       nThetaStart,this%sizeThetaB)
       END IF

       !--------- Since the fields are given at gridpoints here, this is a good
       !          point for graphical output:
       IF ( this%l_graph ) THEN
          PERFON('graphout')
          CALL graphOut_mpi(time,this%nR,ngform,this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,&
               &this%gsa(threadid)%vpc, &
               &        this%gsa(threadid)%brc,this%gsa(threadid)%btc,this%gsa(threadid)%bpc,&
               &this%gsa(threadid)%sc,&
               &        nThetaStart,this%sizeThetaB,.FALSE.)
          PERFOFF
       END IF

       !--------- Helicity output:
       IF ( this%lHelCalc ) THEN
          PERFON('hel_out')
          CALL getHelLM(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,this%gsa(threadid)%vpc,                        &
               &        this%gsa(threadid)%cvrc,this%gsa(threadid)%dvrdtc,&
               &this%gsa(threadid)%dvrdpc,this%gsa(threadid)%dvtdrc,this%gsa(threadid)%dvpdrc,   &
               &        HelLMr,Hel2LMr,         &
               &        HelnaLMr,Helna2LMr,     &
               &        this%nR,nThetaStart)
          PERFOFF
       END IF

       !--------- horizontal velocity :

       IF ( l_viscBcCalc ) THEN
          !WRITE(*,"(2I3,A,3('(',2ES20.12,')'))") nR,nThetaB,&
          !     &" dsdr,dsdt,dsdp = ",&
          !     & SUM(this%gsa(threadid)%drSc),&
          !     & SUM(this%gsa(threadid)%dsdtc),&
          !     & SUM(this%gsa(threadid)%dsdpc)

          CALL get_nlBLayers(this%gsa(threadid)%vtc,    &
               &             this%gsa(threadid)%vpc,    &
               &             this%gsa(threadid)%dvtdrc, &
               &             this%gsa(threadid)%dvpdrc, &
               &             this%gsa(threadid)%drSc,   &
               &             this%gsa(threadid)%dsdtc,  &
               &             this%gsa(threadid)%dsdpc,  &
               &             uhLMr,duhLMr,gradsLMr,nR,nThetaStart)
          !WRITE(*,"(2I3,A,3ES20.12)") nR,nThetaB,&
          !     &" uh,duh,grads = ",&
          !     & SUM(uhLMr(:)),SUM(duhLMr(:)),SUM(gradsLMr(:))
       END IF


       IF ( l_fluxProfs ) THEN
           CALL get_fluxes(this%gsa(threadid)%vrc, &
                  &        this%gsa(threadid)%vtc, &
                  &        this%gsa(threadid)%vpc, &
                  &        this%gsa(threadid)%dvrdrc,  &
                  &        this%gsa(threadid)%dvtdrc,  &
                  &        this%gsa(threadid)%dvpdrc,  &
                  &        this%gsa(threadid)%dvrdtc,  &
                  &        this%gsa(threadid)%dvrdpc,  &
                  &        this%gsa(threadid)%sc, &
                  &        this%gsa(threadid)%pc, &
                  &        this%gsa(threadid)%brc, &
                  &        this%gsa(threadid)%btc, &
                  &        this%gsa(threadid)%bpc, &
                  &        this%gsa(threadid)%cbtc,&
                  &        this%gsa(threadid)%cbpc,&
                  &        fconvLMr,fkinLMr,fviscLMr,fpoynLMr,&
                  &        fresLMr,nR,nThetaStart)
       END IF


       !--------- Movie output:
       IF ( this%l_frame .AND. l_movie_oc .AND. l_store_frame ) THEN
          PERFON('mov_out')
          CALL store_movie_frame(this%nR,this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,&
               &this%gsa(threadid)%vpc,                   &
               &                 this%gsa(threadid)%brc,this%gsa(threadid)%btc,&
               &this%gsa(threadid)%bpc,this%gsa(threadid)%sc,this%gsa(threadid)%drSc,              &
               &                 this%gsa(threadid)%dvrdpc,this%gsa(threadid)%dvpdrc,&
               &this%gsa(threadid)%dvtdrc,this%gsa(threadid)%dvrdtc,this%gsa(threadid)%cvrc, &
               &                 this%gsa(threadid)%cbrc,this%gsa(threadid)%cbtc,&
               &nThetaStart,this%sizeThetaB,this%leg_helper%bCMB)
          PERFOFF
       END IF


       !--------- Stuff for special output:
       !--------- Calculation of magnetic field production and advection terms
       !          for graphic output:
       IF ( l_dtB ) THEN
          PERFON('dtBLM')
          CALL get_dtBLM(this%nR,this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,&
               &this%gsa(threadid)%vpc,this%gsa(threadid)%brc,this%gsa(threadid)%btc,&
               &this%gsa(threadid)%bpc,                 &
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
          CALL getTOnext(this%leg_helper%zAS,this%gsa(threadid)%brc,this%gsa(threadid)%btc,&
               &this%gsa(threadid)%bpc,this%lTONext,     &
               &         this%lTONext2,dt,dtLast,this%nR,nThetaStart,this%sizeThetaB, &
               &         this%BsLast,this%BpLast,this%BzLast)
       END IF

       IF ( this%lTOCalc ) THEN
          CALL getTO(this%gsa(threadid)%vrc,this%gsa(threadid)%vtc,this%gsa(threadid)%vpc,&
               &this%gsa(threadid)%cvrc,this%gsa(threadid)%dvpdrc,         &
               &     this%gsa(threadid)%brc,this%gsa(threadid)%btc,this%gsa(threadid)%bpc,&
               &this%gsa(threadid)%cbrc,this%gsa(threadid)%cbtc,         &
               &     this%BsLast,this%BpLast,this%BzLast,         &
               &     this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM,&
               &     this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM,         &
               &     dtLast,this%nR,nThetaStart,this%sizeThetaB)
       END IF
       PERFOFF

       !IF (this%nR.EQ.n_r_icb) THEN
       !   !$OMP ORDERED
       !   WRITE(*,"(2I3,A,I4,F21.17)") nThetaB,threadid,": lorentz_torque_ic = ",&
       !        & EXPONENT(this%lorentz_torque_ic(threadid)),FRACTION(this%lorentz_torque_ic(threadid))
       !   !$OMP END ORDERED
       !END IF
    END DO ! Loop over theta blocks
    !$OMP END DO

    ! do a reduction over the threads by hand
    ! parallelize over the different arrays
    !PERFON('reduc')
    !$OMP SECTIONS private(iThread)
    !$OMP SECTION
    DO iThread=1,this%nThreads-1
       this%nl_lm(0)%AdvrLM=this%nl_lm(0)%AdvrLM + this%nl_lm(iThread)%AdvrLM
       this%nl_lm(0)%AdvtLM=this%nl_lm(0)%AdvtLM + this%nl_lm(iThread)%AdvtLM
       this%nl_lm(0)%AdvpLM=this%nl_lm(0)%AdvpLM + this%nl_lm(iThread)%AdvpLM
    END DO

    !$OMP SECTION
    DO iThread=1,this%nThreads-1
       this%nl_lm(0)%LFrLM=this%nl_lm(0)%LFrLM + this%nl_lm(iThread)%LFrLM
       this%nl_lm(0)%LFtLM=this%nl_lm(0)%LFtLM + this%nl_lm(iThread)%LFtLM
       this%nl_lm(0)%LFpLM=this%nl_lm(0)%LFpLM + this%nl_lm(iThread)%LFpLM
    END DO

    !$OMP SECTION
    DO iThread=1,this%nThreads-1
       this%nl_lm(0)%VxBrLM=this%nl_lm(0)%VxBrLM + this%nl_lm(iThread)%VxBrLM
       this%nl_lm(0)%VxBtLM=this%nl_lm(0)%VxBtLM + this%nl_lm(iThread)%VxBtLM
       this%nl_lm(0)%VxBpLM=this%nl_lm(0)%VxBpLM + this%nl_lm(iThread)%VxBpLM
    END DO

    !$OMP SECTION
    DO iThread=1,this%nThreads-1
       this%nl_lm(0)%VSrLM=this%nl_lm(0)%VSrLM + this%nl_lm(iThread)%VSrLM
       this%nl_lm(0)%VStLM=this%nl_lm(0)%VStLM + this%nl_lm(iThread)%VStLM
       this%nl_lm(0)%VSpLM=this%nl_lm(0)%VSpLM + this%nl_lm(iThread)%VSpLM
    END DO

    !$OMP SECTION
    DO iThread=1,this%nThreads-1
       this%nl_lm(0)%ViscHeatLM=this%nl_lm(0)%ViscHeatLM + this%nl_lm(iThread)%ViscHeatLM
       this%nl_lm(0)%OhmLossLM=this%nl_lm(0)%OhmLossLM + this%nl_lm(iThread)%OhmLossLM
       this%lorentz_torque_ma(0) = this%lorentz_torque_ma(0) + this%lorentz_torque_ma(iThread)
    END DO

!!$    lorentz_torque_ic=0.0D0
!!$    c=0.0D0
!!$    DO iThread=0,this%nThreads-1
!!$       y=this%lorentz_torque_ic(iThread)-c
!!$       t=lorentz_torque_ic+y
!!$       c=(t-lorentz_torque_ic)-y
!!$       lorentz_torque_ic=t
!!$    END DO

    !$OMP END SECTIONS
    !PERFOFF
    !$OMP END PARALLEL
    lorentz_torque_ic = lorentz_torques_ic(1)
    DO nThetaB=2,this%nThetaBs
       lorentz_torque_ic = lorentz_torque_ic + lorentz_torques_ic(nThetaB)
    END DO
    !lorentz_torque_ic = this%lorentz_torque_ic(0)
    lorentz_torque_ma = this%lorentz_torque_ma(0)

    !IF (this%nR.EQ.n_r_icb) THEN
    !   WRITE(*,"(A,2ES20.12)") "after OMP PARALLEL, lorentz_torque = ",&
    !        & lorentz_torque_ic,lorentz_torque_ma
    !END IF

    IF (DEBUG_OUTPUT) THEN
       CALL this%nl_lm(0)%output()
    END IF
    !-- Partial calculation of time derivatives (horizontal parts):
    !   input flm...  is in (l,m) space at radial grid points this%nR !
    !   Only dVxBh needed for boundaries !
    !   get_td finally calculates the d*dt terms needed for the 
    !   time step performed in s_LMLoop.f . This should be distributed
    !   over the different models that s_LMLoop.f parallelizes over. 
    !write(*,"(A,I4,2ES20.13)") "before_td: ",this%nR,sum(real(conjg(VxBtLM)*VxBtLM)),sum(real(conjg(VxBpLM)*VxBpLM))
    !PERFON('get_td')
    CALL get_td(this%nR,this%nBc,this%lRmsCalc,dVSrLM,dVxBhLM,   &
         &      dwdt,dzdt,dpdt,   &
         &      dsdt,dbdt,djdt,   &
         &      this%nl_lm, this%leg_helper)
    !PERFOFF
    !write(*,"(A,I4,ES20.13)") "after_td:  ",this%nR,sum(real(conjg(dVxBhLM(:,this%nR_Mag))*dVxBhLM(:,this%nR_Mag)))
    !-- Finish calculation of TO variables:
    IF ( this%lTOcalc ) THEN                                   
       CALL getTOfinish(this%nR,dtLast,this%leg_helper%zAS,this%leg_helper%dzAS,this%leg_helper%ddzAS, &
            &           this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM,&
            &           this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM)
    END IF

    !--- Form partial horizontal derivaties of magnetic production and
    !    advection terms:
    IF ( l_dtB ) THEN
       PERFON('dtBLM')
       CALL get_dH_dtBLM(this%nR,this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM,&
            &            this%dtB_arrays%BrVtLM,                      &
            &            this%dtB_arrays%BrVpLM,this%dtB_arrays%BtVpLM,&
            &            this%dtB_arrays%BpVtLM,this%dtB_arrays%BrVZLM,&
            &            this%dtB_arrays%BtVZLM,this%dtB_arrays%BtVpCotLM, &
            &            this%dtB_arrays%BpVtCotLM,this%dtB_arrays%BtVZcotLM,&
            &            this%dtB_arrays%BtVpSn2LM,                &
            &            this%dtB_arrays%BpVtSn2LM,this%dtB_arrays%BtVZsn2LM)
       PERFOFF
    END IF
  END SUBROUTINE do_iteration_ThetaBlocking_OpenMP

END MODULE rIterThetaBlocking_OpenMP_mod
