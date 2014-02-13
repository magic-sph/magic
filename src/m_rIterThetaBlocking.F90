#include "perflib_preproc.cpp"
MODULE rIterThetaBlocking_mod
  USE rIteration_mod, only: rIteration_t

  USE truncation, ONLY: lm_max,lmP_max,nrp,ncp,l_max,lmP_max_dtB,&
       &n_phi_maxStr,n_theta_maxStr,n_r_maxStr
  USE blocking, only: nfs
  USE logic, ONLY: l_mag,l_conv,l_mag_kin,l_heat,l_ht,l_anel,l_mag_LF,&
       & l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic, l_cond_ic, &
       & l_rot_ma, l_cond_ma, l_viscBcCalc, l_dtB, l_store_frame, l_movie_oc
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
  USE grid_space_arrays_mod,ONLY: grid_space_arrays_t
#ifdef WITH_LIKWID
#   include "likwid_f90.h"
#endif
  IMPLICIT NONE

  TYPE :: dtB_arrays_t
     !----- Local dtB output stuff:
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: BtVrLM, BpVrLM, BrVtLM, BrVpLM, BtVpLM, BpVtLM
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: BtVpCotLM, BpVtCotLM, BtVpSn2LM, BpVtSn2LM
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: BrVZLM, BtVZLM, BtVZcotLM, BtVZsn2LM
  END TYPE dtB_arrays_t

  TYPE :: TO_arrays_t
     !----- Local TO output stuff:
     REAL(kind=8),ALLOCATABLE :: dzRstrLM(:),dzAstrLM(:)
     REAL(kind=8),ALLOCATABLE :: dzCorLM(:),dzLFLM(:)
  END TYPE TO_arrays_t

  
  TYPE,PUBLIC,ABSTRACT,EXTENDS(rIteration_t) :: rIterThetaBlocking_t
     ! or with len parameters for the theta block size and number
     !TYPE,PUBLIC,EXTENDS(rIteration_t) :: rIterThetaBlocking_t(sizeThetaB,nThetaBs)
     !integer,len :: sizeThetaB,nThetaBs
     INTEGER :: sizeThetaB, nThetaBs

     !TYPE(nonlinear_lm_t) :: nl_lm
     type(leg_helper_t)   :: leg_helper
     type(dtB_arrays_t)   :: dtB_arrays
     type(TO_arrays_t)    :: TO_arrays

     !CLASS(grid_space_arrays_t),PRIVATE :: gsa

     !----- Saved magnetic field components from last time step:
     !      This is needed for the current TO version. However,
     !      the variables calulated with this don't give any
     !      deep insight. TO should be changes in the future to
     !      eliminate this.
     REAL(kind=8),DIMENSION(:,:,:),ALLOCATABLE :: BsLast, BpLast, BzLast
   CONTAINS
     !PROCEDURE :: initialize => initialize_rIterThetaBlocking
     procedure :: allocate_common_arrays
     procedure :: deallocate_common_arrays
     procedure :: set_ThetaBlocking
     !PROCEDURE,DEFERRED :: do_iteration
     procedure :: transform_to_grid_space
     procedure :: transform_to_lm_space
  END TYPE rIterThetaBlocking_t

CONTAINS
  SUBROUTINE allocate_common_arrays(this)
    class(rIterThetaBlocking_t) :: this

    !----- Nonlinear terms in lm-space: 
    !CALL this%nl_lm%initialize(lmP_max)


    !----- Help arrays for Legendre transform calculated in legPrepG:
    !      Parallelizatio note: these are the R-distributed versions
    !      of the field scalars.
    CALL this%leg_helper%initialize(lm_max,l_max)

    !----- Local dtB output stuff:
    ALLOCATE( this%dtB_arrays%BtVrLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BpVrLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BrVtLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BrVpLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BtVpLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BpVtLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BtVpCotLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BpVtCotLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BtVpSn2LM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BpVtSn2LM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BrVZLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BtVZLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BtVZcotLM(lmP_max_dtB) )
    ALLOCATE( this%dtB_arrays%BtVZsn2LM(lmP_max_dtB) )

    !----- Local TO output stuff:
    ALLOCATE( this%TO_arrays%dzRstrLM(l_max+2),this%TO_arrays%dzAstrLM(l_max+2) )
    ALLOCATE( this%TO_arrays%dzCorLM(l_max+2),this%TO_arrays%dzLFLM(l_max+2) )

    allocate( this%BsLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr) )
    allocate( this%BpLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr) )
    allocate( this%BzLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr) )

  END SUBROUTINE allocate_common_arrays

  SUBROUTINE deallocate_common_arrays(this)
    class(rIterThetaBlocking_t) :: this
    DEALLOCATE( this%dtB_arrays%BtVrLM )
    DEALLOCATE( this%dtB_arrays%BpVrLM )
    DEALLOCATE( this%dtB_arrays%BrVtLM )
    DEALLOCATE( this%dtB_arrays%BrVpLM )
    DEALLOCATE( this%dtB_arrays%BtVpLM )
    DEALLOCATE( this%dtB_arrays%BpVtLM )
    DEALLOCATE( this%dtB_arrays%BtVpCotLM )
    DEALLOCATE( this%dtB_arrays%BpVtCotLM )
    DEALLOCATE( this%dtB_arrays%BtVpSn2LM )
    DEALLOCATE( this%dtB_arrays%BpVtSn2LM )
    DEALLOCATE( this%dtB_arrays%BrVZLM )
    DEALLOCATE( this%dtB_arrays%BtVZLM )
    DEALLOCATE( this%dtB_arrays%BtVZcotLM )
    DEALLOCATE( this%dtB_arrays%BtVZsn2LM )

    !----- Local TO output stuff:
    DEALLOCATE( this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM )
    DEALLOCATE( this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM )

    deallocate( this%BsLast)
    deallocate( this%BpLast)
    deallocate( this%BzLast)
  END SUBROUTINE deallocate_common_arrays

  SUBROUTINE set_ThetaBlocking(this,nThetaBs,sizeThetaB)
    CLASS(rIterThetaBlocking_t) :: this
    INTEGER,INTENT(IN) :: nThetaBs, sizeThetaB
    
    this%nThetaBs = nThetaBs
    this%sizeThetaB = sizeThetaB
  END SUBROUTINE set_ThetaBlocking

  ! Generally usable functions
  SUBROUTINE transform_to_grid_space(this,nThetaStart,nThetaStop,gsa)
    class(rIterThetaBlocking_t) :: this
    INTEGER, INTENT(IN) :: nThetaStart,nThetaStop
    type(grid_space_arrays_t) :: gsa

    ! Local variables
    INTEGER :: nTheta,nPhi,lm
    LOGICAL :: DEBUG_OUTPUT=.FALSE.
    !----- Legendre transform from (r,l,m) to (r,theta,m):
    !      First version with PlmTF needed for first-touch policy  
    IF ( l_mag ) THEN
       !PERFON('legTFG')
       !LIKWID_ON('legTFG')
       CALL legTFG(this%nBc,this%lDeriv,nThetaStart,                &
            &      gsa%vrc,gsa%vtc,gsa%vpc,gsa%dvrdrc,&
            &      gsa%dvtdrc,gsa%dvpdrc,gsa%cvrc, &
            &      gsa%dvrdtc,gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,           &
            &      gsa%brc,gsa%btc,gsa%bpc,gsa%cbrc,&
            &      gsa%cbtc,gsa%cbpc,gsa%sc,gsa%drSc,    &
            &      gsa%dsdtc, gsa%dsdpc, &
            &      this%leg_helper)
       !LIKWID_OFF('legTFG')
       !PERFOFF
       IF (DEBUG_OUTPUT) THEN
          DO nTheta=1,this%sizeThetaB
             WRITE(*,"(2I3,A,6ES20.12)") this%nR,nTheta,": sum v = ",&
                  &SUM(gsa%vrc(:,nTheta))!,SUM(vtc(:,nTheta)),SUM(vpc(:,nTheta))
          END DO
       END IF
    ELSE
       !PERFON('legTFGnm')
       !LIKWID_ON('legTFGnm')
       CALL legTFGnomag(this%nBc,this%lDeriv,nThetaStart,                 &
            &           gsa%vrc,gsa%vtc,gsa%vpc,gsa%dvrdrc,&
            &           gsa%dvtdrc,gsa%dvpdrc,gsa%cvrc,  &
            &           gsa%dvrdtc,gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,&
            &           gsa%sc,gsa%drSc,    &
            &           this%leg_helper%dLhw,this%leg_helper%dLhdw,this%leg_helper%dLhz,&
            &           this%leg_helper%vhG,this%leg_helper%vhC,this%leg_helper%dvhdrG,&
            &           this%leg_helper%dvhdrC,&
            &           this%leg_helper%sR,this%leg_helper%dsR)
       !LIKWID_OFF('legTFGnm')
       !PERFOFF
    END IF

    !------ Fourier transform from (r,theta,m) to (r,theta,phi):
    IF ( l_conv .OR. l_mag_kin ) THEN
       IF ( l_heat ) CALL fft_thetab(gsa%sc,1)
       IF ( l_HT ) CALL fft_thetab(gsa%drSc,1)
       IF ( this%nBc.EQ.0 ) THEN
          CALL fft_thetab(gsa%vrc,1)
          CALL fft_thetab(gsa%vtc,1)
          CALL fft_thetab(gsa%vpc,1)
          IF ( this%lDeriv ) THEN
             CALL fft_thetab(gsa%dvrdrc,1)
             CALL fft_thetab(gsa%dvtdrc,1)
             CALL fft_thetab(gsa%dvpdrc,1)
             CALL fft_thetab(gsa%cvrc,1)
             CALL fft_thetab(gsa%dvrdtc,1)
             CALL fft_thetab(gsa%dvrdpc,1)
             CALL fft_thetab(gsa%dvtdpc,1)
             CALL fft_thetab(gsa%dvpdpc,1)
          END IF
       ELSE IF ( this%nBc.EQ.1 ) THEN ! Stress free
          gsa%vrc = CMPLX(0.D0,0.D0,kind=KIND(gsa%vrc))
          CALL fft_thetab(gsa%vtc,1)
          CALL fft_thetab(gsa%vpc,1)
          IF ( this%lDeriv ) THEN
             gsa%dvrdtc = CMPLX(0.D0,0.D0,kind=KIND(gsa%dvrdtc))
             gsa%dvrdpc = CMPLX(0.D0,0.D0,kind=KIND(gsa%dvrdpc))
             CALL fft_thetab(gsa%dvrdrc,1)
             CALL fft_thetab(gsa%dvtdrc,1)
             CALL fft_thetab(gsa%dvpdrc,1)
             CALL fft_thetab(gsa%cvrc,1)
             CALL fft_thetab(gsa%dvtdpc,1)
             CALL fft_thetab(gsa%dvpdpc,1)
          END IF
       ELSE IF ( this%nBc.EQ.2 ) THEN 
          IF ( this%nR.EQ.n_r_cmb ) THEN
             CALL v_rigid_boundary(this%nR,this%leg_helper%omegaMA,this%lDeriv, &
                  &                gsa%vrc,gsa%vtc,gsa%vpc,&
                  &                gsa%cvrc,gsa%dvrdtc,gsa%dvrdpc,&
                  &                gsa%dvtdpc,gsa%dvpdpc, &
                  &                nThetaStart)
          ELSE IF ( this%nR.EQ.n_r_icb ) THEN
             CALL v_rigid_boundary(this%nR,this%leg_helper%omegaIC,this%lDeriv, &
                  &                gsa%vrc,gsa%vtc,gsa%vpc,      &
                  &                gsa%cvrc,gsa%dvrdtc,gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc, &
                  &                nThetaStart)
          END IF
          IF ( this%lDeriv ) THEN
             CALL fft_thetab(gsa%dvrdrc,1)
             CALL fft_thetab(gsa%dvtdrc,1)
             CALL fft_thetab(gsa%dvpdrc,1)
          END IF
       END IF
    END IF
    IF ( l_mag .OR. l_mag_LF ) THEN
       CALL fft_thetab(gsa%brc,1)
       CALL fft_thetab(gsa%btc,1)
       CALL fft_thetab(gsa%bpc,1)
       IF ( this%lDeriv ) THEN
          CALL fft_thetab(gsa%cbrc,1)
          CALL fft_thetab(gsa%cbtc,1)
          CALL fft_thetab(gsa%cbpc,1)
       END IF
    END IF
  END SUBROUTINE transform_to_grid_space

  SUBROUTINE transform_to_lm_space(this,nThetaStart,nThetaStop,gsa,nl_lm)
    class(rIterThetaBlocking_t) :: this
    INTEGER,INTENT(IN) :: nThetaStart,nThetaStop
    type(grid_space_arrays_t) :: gsa
    type(nonlinear_lm_t) :: nl_lm
    
    ! Local variables
    INTEGER :: nTheta,nPhi

    IF ( (.NOT.this%isRadialBoundaryPoint) .AND. ( l_conv_nl .OR. l_mag_LF ) ) THEN
       !PERFON('inner1')
       IF ( l_conv_nl .AND. l_mag_LF ) THEN
          DO nTheta=1,this%sizeThetaB
             DO nPhi=1,nrp
                gsa%Advr(nPhi,nTheta)=gsa%Advr(nPhi,nTheta) + gsa%LFr(nPhi,nTheta)
                gsa%Advt(nPhi,nTheta)=gsa%Advt(nPhi,nTheta) + gsa%LFt(nPhi,nTheta)
                gsa%Advp(nPhi,nTheta)=gsa%Advp(nPhi,nTheta) + gsa%LFp(nPhi,nTheta)
             END DO
          END DO
       ELSE IF ( l_mag_LF ) THEN
          DO nTheta=1,this%sizeThetaB
             DO nPhi=1,nrp
                gsa%Advr(nPhi,nTheta)=gsa%LFr(nPhi,nTheta)
                gsa%Advt(nPhi,nTheta)=gsa%LFt(nPhi,nTheta)
                gsa%Advp(nPhi,nTheta)=gsa%LFp(nPhi,nTheta)
             END DO
          END DO
       END IF
       CALL fft_thetab(gsa%Advr,-1)
       CALL fft_thetab(gsa%Advt,-1)
       CALL fft_thetab(gsa%Advp,-1)
       CALL legTF3(nThetaStart,nl_lm%AdvrLM,nl_lm%AdvtLM,nl_lm%AdvpLM,    &
            &      gsa%Advr,gsa%Advt,gsa%Advp)
       IF ( this%lRmsCalc .AND. l_mag_LF ) THEN ! LF treated extra:
          CALL fft_thetab(gsa%LFr,-1)
          CALL fft_thetab(gsa%LFt,-1)
          CALL fft_thetab(gsa%LFp,-1)
          CALL legTF3(nThetaStart,nl_lm%LFrLM,nl_lm%LFtLM,nl_lm%LFpLM,    &
               &      gsa%LFr,gsa%LFt,gsa%LFp)
       END IF
       !PERFOFF
    END IF
    IF ( (.NOT.this%isRadialBoundaryPoint) .AND. l_heat ) THEN
       !PERFON('inner2')
       CALL fft_thetab(gsa%VSr,-1)
       CALL fft_thetab(gsa%VSt,-1)
       CALL fft_thetab(gsa%VSp,-1)
       CALL legTF3(nThetaStart,nl_lm%VSrLM,nl_lm%VStLM,nl_lm%VSpLM,       &
            &      gsa%VSr,gsa%VSt,gsa%VSp)
       IF (l_anel) THEN ! anelastic stuff 
          IF (l_mag_nl) THEN
             CALL fft_thetab(gsa%ViscHeat,-1)
             CALL fft_thetab(gsa%OhmLoss,-1)
             CALL legTF2(nThetaStart,nl_lm%OhmLossLM,         &
                  &      nl_lm%ViscHeatLM,gsa%OhmLoss,            &
                  &      gsa%ViscHeat)
          ELSE
             CALL fft_thetab(gsa%ViscHeat,-1)
             CALL legTF1(nThetaStart,nl_lm%ViscHeatLM,gsa%ViscHeat)
          END IF
       END IF
       !PERFOFF
    END IF
    IF ( l_mag_nl ) THEN
       !PERFON('mag_nl')
       IF ( .not.this%isRadialBoundaryPoint ) THEN
          CALL fft_thetab(gsa%VxBr,-1)
          CALL fft_thetab(gsa%VxBt,-1)
          CALL fft_thetab(gsa%VxBp,-1)
          CALL legTF3(nThetaStart,nl_lm%VxBrLM,nl_lm%VxBtLM,nl_lm%VxBpLM,&
               &       gsa%VxBr,gsa%VxBt,gsa%VxBp)
       ELSE
          !write(*,"(I4,A,ES20.13)") this%nR,", VxBt = ",sum(VxBt*VxBt)
          CALL fft_thetab(gsa%VxBt,-1)
          CALL fft_thetab(gsa%VxBp,-1)
          CALL legTF2(nThetaStart,nl_lm%VxBtLM,nl_lm%VxBpLM,&
               &      gsa%VxBt,gsa%VxBp)
       END IF
       !PERFOFF
    END IF

  END SUBROUTINE transform_to_lm_space

END MODULE rIterThetaBlocking_mod
