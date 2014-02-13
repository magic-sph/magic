!$Id$
#include "perflib_preproc.cpp"
MODULE radialLoop
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE torsional_oscillations
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE const
  USE parallel_mod, ONLY: rank, n_procs
  USE radial_data,ONLY: nRstart,nRstop,n_r_cmb
  USE legendre_trafo,only: legTFG
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
#ifdef WITH_LIKWID
#  include "likwid_f90.h"
#endif
  USE rIteration_mod, only: rIteration_t
  USE rIterThetaBlocking_mod, only: rIterThetaBlocking_t
  USE rIterThetaBlocking_seq_mod,only: rIterThetaBlocking_seq_t
  USE rIterThetaBlocking_OpenMP_mod,only: rIterThetaBlocking_OpenMP_t
  IMPLICIT NONE

  PRIVATE
  !---- Nonlinear terms, field components and help arrays
  !     for Legendre transform, field components, and dtB and 
  !     TO output. These are all stored in COMMON BLOCKS
  !     that are thread privat rather than using explicit 
  !     PRIVATE statements in the PARALLEL DO clause.


  ! public elements of the module

  PUBLIC :: initialize_radialLoop,finalize_radialLoop,radialLoopG

  CLASS(rIteration_t),pointer :: this_rIteration

CONTAINS
  SUBROUTINE initialize_radialLoop
    CHARACTER(len=100) :: this_type

#ifdef WITHOMP
    allocate( rIterThetaBlocking_OpenMP_t :: this_rIteration )
#else
    ALLOCATE( rIterThetaBlocking_seq_t :: this_rIteration )
#endif
    this_type = this_rIteration%getType()
    WRITE(*,"(2A)") "Using rIteration type: ",TRIM(this_type)
    CALL this_rIteration%initialize()
    SELECT TYPE (this_rIteration)
    CLASS is (rIterThetaBlocking_t)
       CALL this_rIteration%set_ThetaBlocking(nThetaBs,sizeThetaB)
    CLASS default
       PRINT*,"this_rIteration has no matching type in m_radialLoop.F90"
    END SELECT
  END SUBROUTINE initialize_radialLoop

  SUBROUTINE finalize_radialLoop
    call this_rIteration%finalize()
    DEALLOCATE(this_rIteration)
  END SUBROUTINE finalize_radialLoop

  !***********************************************************************
  SUBROUTINE radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,        &
       &                 lTOCalc,lTONext,lTONext2,lHelCalc,lRmsCalc,   &
       &                 dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM, &
       &                 lorentz_torque_ic,lorentz_torque_ma,          &
       &                 br_vt_lm_cmb,br_vp_lm_cmb,                    &
       &                 br_vt_lm_icb,br_vp_lm_icb,                    &
       &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,&
       &                 gradsLMr,dtrkc,dthkc)
    !***********************************************************************

    !    !------------ This is release 2 level 10  --------------!
    !    !------------ Created on 2/5/02  by JW. -----------

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine performs the actual time-stepping.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+

    !--- Input of variables:
    LOGICAL,intent(IN) :: l_graph,l_cour,l_frame
    LOGICAL,intent(IN) :: lTOcalc,lTONext,lTONext2,lHelCalc
    LOGICAL,intent(IN) :: lRmsCalc
    REAL(kind=8),intent(IN) :: time,dt,dtLast

    !---- Output of explicit time step:
    !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
    !     need a further treatment (radial derivatives required):
    COMPLEX(kind=8),INTENT(OUT),DIMENSION(lm_max,nRstart:nRstop) :: &
         & dwdt,dzdt,dpdt,dsdt,dVSrLM
    COMPLEX(kind=8),INTENT(OUT),DIMENSION(lm_maxMag,nRstartMag:nRstopMag) :: &
         & dbdt,djdt,dVxBhLM
    REAL(kind=8),intent(OUT) :: lorentz_torque_ma,lorentz_torque_ic

    !---- Output for axisymmetric helicity:
    REAL(kind=8),INTENT(OUT),DIMENSION(l_max+1,nRstart:nRstop) :: &
         & HelLMr,Hel2LMr,HelnaLMr,Helna2LMr
    REAL(kind=8),INTENT(OUT),DIMENSION(l_max+1,nRstart:nRstop) :: &
         & uhLMr,duhLMr,gradsLMr

    !---- Output of nonlinear products for nonlinear
    !     magnetic boundary conditions (needed in s_updateB.f):
    COMPLEX(kind=8),intent(OUT) :: br_vt_lm_cmb(lmP_max) ! product br*vt at CMB
    COMPLEX(kind=8),intent(OUT) :: br_vp_lm_cmb(lmP_max) ! product br*vp at CMB
    COMPLEX(kind=8),intent(OUT) :: br_vt_lm_icb(lmP_max) ! product br*vt at ICB
    COMPLEX(kind=8),intent(OUT) :: br_vp_lm_icb(lmP_max) ! product br*vp at ICB

    !---- Output for Courant criteria:
    REAL(kind=8),INTENT(OUT) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)


    !--- Local variables:

    !---- Counter, logicals ...
    INTEGER :: nR,nr_Mag,nBc,lm
    !INTEGER :: nTheta,nThetaB,nThetaLast
    INTEGER :: nThetaStart!,nThetaStop
    LOGICAL :: lDeriv,lOutBc,lMagNlBc
    LOGICAL :: lGraphHeader    ! Write header into graph file
    logical :: isRadialBoundaryPoint

    !---- End of declaration
    !----------------------------------------------------------------------
    PERFON('rloop')
    !LIKWID_ON('rloop')
    lGraphHeader=l_graph
    IF ( lGraphHeader ) THEN
       CALL graphOut_mpi_header(time,nR,ngform,&
            &        nThetaStart,sizeThetaB)
    END IF

    IF ( l_cour ) THEN
       IF (rank.EQ.0) THEN
          dtrkc(n_r_cmb)=1.D10
          dthkc(n_r_cmb)=1.D10
       ELSEIF (rank.EQ.n_procs-1) THEN
          dtrkc(n_r_icb)=1.D10
          dthkc(n_r_icb)=1.D10
       END IF
    END IF

    !------ Set nonlinear terms that are possibly needed at the boundaries.
    !       They may be overwritten by get_td later.
    DO lm=1,lm_max
       IF (rank.EQ.0) THEN
          dVSrLM(lm,n_r_cmb) =zero
          IF ( l_mag ) THEN
             dVxBhLM(lm,n_r_cmb)=zero
          END IF
       ELSEIF (rank.EQ.n_procs-1) then
          dVSrLM(lm,n_r_icb) =zero
          IF ( l_mag ) THEN
             dVxBhLM(lm,n_r_icb)=zero
          END IF
       END IF
    END DO

    !------ Having to calculate non-linear boundary terms?
    lMagNlBc=.FALSE.
    IF ( ( l_mag_nl .OR. l_mag_kin ) .AND.                          &
         &       ( ktopv.EQ.1 .OR. l_cond_ma .OR.                           &
         &          ( ktopv.EQ.2 .AND. l_rot_ma ) ) .OR.                    &
         &       ( kbotv.EQ.1 .OR. l_cond_ic .OR.                           &
         &          ( kbotv.EQ.2 .AND. l_rot_ic ) ) )                       &
         &     lMagNlBc=.TRUE.

    !------ When boundary output, Courant criterion, or non-magnetic 
    !       boundary conditions are required I have to calculate 
    !       the fields at the boundaries. This is done in one thread and 
    !       is triggered by lOutBc=.TRUE.
    lOutBc=.FALSE.
    IF ( lTOCalc .OR. lHelCalc .OR. l_frame .OR.         &
         & l_cour .OR. l_dtB .OR. lMagNlBc .OR. l_graph  &
         &     ) lOutBc=.TRUE.

    !nRstart=n_r_cmb
    !nRstop =n_r_icb-1

    !--- Start the big do loop over the radial threads:

    !nThreadsRmax=1
    DO nR=nRstart,nRstop
       !IF( nTh.GT.nThreadsRmax ) nThreadsRmax=nTh
       IF ( lVerbose ) THEN
          WRITE(*,'(/," ! Starting radial level ",i4)') nR
       END IF

       !nR = nRC
       nBc = 0
       lDeriv = .true.
       isRadialBoundaryPoint=(nR.EQ.n_r_cmb).OR.(nR.EQ.n_r_icb)

       IF ( nR.EQ.n_r_cmb ) THEN 
          IF ( lOutBc ) THEN
             !nR  = n_r_cmb
             nBc = ktopv
             lDeriv= lTOCalc .OR. lHelCalc .OR. l_frame 
          ELSE
             CYCLE   ! Nothing needs to be done by thread one !
          END IF
       ELSEif ( nR.eq.n_r_icb ) then
          IF ( lOutBc ) THEN
             !nR = n_r_icb
             nBc = kbotv
             lDeriv= lTOCalc .OR. lHelCalc .OR. l_frame 
          ELSE
             CYCLE
          END IF
       END IF

       IF ( l_mag .OR. l_mag_LF ) THEN
          nR_Mag=nR
       ELSE
          nR_Mag=1
       END IF

       CALL this_rIteration%set_steering_variables(l_cour,lTOCalc,lTOnext,lTOnext2,&
            & lDeriv,lRmsCalc,lHelCalc,l_frame,lMagNlBc,l_graph)

       CALL this_rIteration%do_iteration(nR,nBc,time,dt,dtLast,&
            & dsdt(:,nR),dwdt(:,nR),dzdt(:,nR),dpdt(:,nR),dbdt(:,nR_Mag),djdt(:,nR_Mag),&
            & dVxBhLM(:,nR_Mag),dVSrLM(:,nR), &
            & br_vt_lm_cmb,br_vp_lm_cmb,&
            & br_vt_lm_icb,br_vp_lm_icb,lorentz_torque_ic,lorentz_torque_ma,&
            & HelLMr(:,nR),Hel2LMr(:,nR),HelnaLMr(:,nR),Helna2LMr(:,nR),&
            & uhLMr(:,nR),duhLMr(:,nR),gradsLMr(:,nR))

       !WRITE(*,"(I3,A,2ES20.12)") nR,": lorentz_torque = ",lorentz_torque_ic,lorentz_torque_ma

       dtrkc(nR)=this_rIteration%dtrkc      
       dthkc(nR)=this_rIteration%dthkc      

    END DO    ! Loop over radial levels 

    !----- Correct sign of mantel Lorentz torque (see above):
    lorentz_torque_ma=-lorentz_torque_ma

    !LIKWID_OFF('rloop')
    PERFOFF
  END SUBROUTINE radialLoopG
END MODULE radialLoop
