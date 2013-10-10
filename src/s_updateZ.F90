!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"
SUBROUTINE updateZ(z,dz,dzdt,dzdtLast,time, &
     &             omega_ma,d_omega_ma_dtLast, &
     &             omega_ic,d_omega_ic_dtLast, &
     &             lorentz_torque_ma,lorentz_torque_maLast, &
     &             lorentz_torque_ic,lorentz_torque_icLast, &
     &             w1,coex,dt,lRmsNext,nTh)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |     Changed by J.Wicht 20.07.2000                                 |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  !  its radial derivatives
  !  adds explicit part to time derivatives of z

  !  Input:  w1 - weight for dbdt-contribution from current time step
  !               (w2=1-w1: weight for contrib. from previous step)
  !          coex - factor depending on weighting alpha of implicit contribution
  !          m1,m2- range of mca-indices in which field is updated
  !                 (harmonic order m=(mca-1)*minc)

  !--------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE torsional_oscillations
  USE init_fields
  USE blocking,ONLY: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map,&
       &lmStartB,lmStopB
  USE horizontal_data
  USE logic
  USE matrices
  USE RMS
  USE const
  use parallel_mod
  USE algebra, ONLY: cgesl,cgeslML
  USE LMLoop_data,ONLY: llm,ulm,llm_real,ulm_real
  USE communications, only:get_global_sum
  IMPLICIT NONE

  !-- Input/output of scalar fields:
  COMPLEX(kind=8),INTENT(INOUT) :: z(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(OUT) :: dz(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN)    :: dzdt(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(INOUT) :: dzdtLast(llm:ulm,n_r_max)
  REAL(kind=8),intent(INOUT) :: d_omega_ma_dtLast,d_omega_ic_dtLast
  REAL(kind=8),intent(OUT) :: omega_ma,omega_ic
  REAL(kind=8),intent(IN) :: lorentz_torque_ma,lorentz_torque_maLast
  REAL(kind=8),intent(IN) :: lorentz_torque_ic,lorentz_torque_icLast
  !-- Output:
  !       Updated z,dz,dzdtLast,d_omega_maLast,d_omega_icLast,ddzAS

  !-- Input of other variables:
  REAL(kind=8),intent(IN) :: time
  REAL(kind=8),intent(IN) :: w1    ! weight for time step !
  REAL(kind=8),intent(IN) :: coex  ! factor depending on alpha
  REAL(kind=8),intent(IN) :: dt
  LOGICAL,intent(IN) :: lRmsNext
  INTEGER,intent(IN) :: nTh ! Thread number

  !-- Input of recycled work arrays:
  COMPLEX(kind=8) :: workA(llm:ulm,n_r_max)
  COMPLEX(kind=8) :: workB(llm:ulm,n_r_max)
  COMPLEX(kind=8) :: workC(llm:ulm,n_r_max)


  !-- local variables:
  REAL(kind=8) :: w2                  ! weight of second time step
  REAL(kind=8) :: O_dt
  REAL(kind=8) :: d_omega_ic_dt,d_omega_ma_dt
  INTEGER :: l1,m1              ! degree and order
  INTEGER :: lm1,lm,lmB         ! position of (l,m) in array
  INTEGER :: lmStart_real       ! range of lm for real array
  INTEGER :: lmStop_real        !
  INTEGER :: lmStart_00         ! excluding l=0,m=0
  INTEGER :: lmStart,lmStop ! max and min number of orders m
  INTEGER :: nLMB2
  INTEGER :: nR                 ! counts radial grid points
  INTEGER :: n_cheb             ! counts cheb modes
  COMPLEX(kind=8) :: rhs(n_r_max)   ! RHS of matrix multiplication
  COMPLEX(kind=8) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! RHS for other modes
  COMPLEX(kind=8) :: z10(n_r_max),z11(n_r_max) ! toroidal flow scalar components
  REAL(kind=8) :: angular_moment(3)   ! total angular momentum
  REAL(kind=8) :: angular_moment_oc(3)! x,y,z component of outer core angular mom.
  REAL(kind=8) :: angular_moment_ic(3)! x,y,z component of inner core angular mom.
  REAL(kind=8) :: angular_moment_ma(3)! x,y,z component of mantle angular mom.
  COMPLEX(kind=8) :: corr_l1m0       ! correction factor for z(l=1,m=0)
  COMPLEX(kind=8) :: corr_l1m1       ! correction factor for z(l=1,m=1)
  REAL(kind=8) :: r_E_2               ! =r**2
  REAL(kind=8) :: nomi                ! nominator for Z10 AM correction
  INTEGER :: l1m0,l1m1          ! position of (l=1,m=0) and (l=1,m=1) in lm.
  INTEGER :: i                  ! counter
  LOGICAL :: l10
  integer :: nLMB

  COMPLEX(kind=8) :: Dif(lm_max)

  INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
  INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
  INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

  logical :: DEBUG_OUTPUT=.false.

  !-- end of declaration
  !-----------------------------------------------------------------------

  !CALL mpi_barrier(MPI_COMM_WORLD,ierr)
  !WRITE(*,"(3(A,2ES20.12))") "begin upZ: dzdt = ",get_global_sum( dzdt ),&
  !     &", z = ",get_global_sum( z ),&
  !     &", dzdtLast = ",get_global_sum( dzdtLast )
  !CALL mpi_barrier(MPI_COMM_WORLD,ierr)

  IF ( .NOT. l_update_v ) RETURN

  nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
  sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
  lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
  lm22l(1:,1:,1:) => lo_sub_map%lm22l
  lm22m(1:,1:,1:) => lo_sub_map%lm22m
  lm2(0:,0:) => lo_map%lm2
  lm2l(1:lm_max) => lo_map%lm2l
  lm2m(1:lm_max) => lo_map%lm2m


  nLMB = 1+rank
  lmStart     =lmStartB(nLMB)
  lmStop      =lmStopB(nLMB)
  lmStart_00  =MAX(2,lmStart)
  lmStart_real=2*lmStart_00-1
  lmStop_real =2*lmStop
  l1m0        =lm2(1,0)

  w2  =1.D0-w1
  O_dt=1.D0/dt

  l10=.FALSE.
  DO nLMB2=1,nLMBs2(nLMB)
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
        IF ( lm1 == l1m0 ) l10= .TRUE. 

        IF ( l_z10mat .AND. lm1 == l1m0 ) THEN
           !WRITE(*,"(A,3I3)") "l_z10mat and lm1=",lm1,l1,m1
           !PERFON('upZ_z10')
           !----- Special treatment of z10 component if ic or mantle
           !      are allowed to rotate about z-axis (l_z10mat=.true.) and
           !      we use no slip boundary condition (ktopv=1,kbotv=1):
           !      Lorentz torque is the explicit part of this time integration
           !      at the boundaries!
           !      Note: no angular momentum correction necessary for this case !
           IF ( .NOT. lZ10mat ) THEN
              CALL get_z10Mat(dt,l1,hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                   z10Mat,z10Pivot)
              lZ10mat=.TRUE.
           END IF
           IF ( l_SRMA ) THEN
              tOmega_ma1=time+tShift_ma1
              tOmega_ma2=time+tShift_ma2
              omega_ma=                                      &
                   omega_ma1*DCOS(omegaOsz_ma1*tOmega_ma1) + &
                   omega_ma2*DCOS(omegaOsz_ma2*tOmega_ma2)
              rhs(1)=omega_ma
           ELSE IF ( ktopv == 2 .AND. l_rot_ma ) THEN  ! time integration
              d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_ma
              rhs(1)=O_dt*c_dt_z10_ma*z(lm1,1) + &
                   w1*d_omega_ma_dt + &
                   w2*d_omega_ma_dtLast
           ELSE
              rhs(1)=0.d0
           END IF
           IF ( l_SRIC ) THEN
              tOmega_ic1=time+tShift_ic1
              tOmega_ic2=time+tShift_ic2
              omega_ic= omega_ic1*COS(omegaOsz_ic1*tOmega_ic1) + &
                   &    omega_ic2*COS(omegaOsz_ic2*tOmega_ic2)
              rhs(n_r_max)=omega_ic
           ELSE if ( kbotv == 2 .AND. l_rot_ic ) then  ! time integration
              d_omega_ic_dt=LFfac*c_lorentz_ic*lorentz_torque_ic
              rhs(n_r_max)=O_dt*c_dt_z10_ic*z(lm1,n_r_max) + &
                   w1*d_omega_ic_dt + &
                   w2*d_omega_ic_dtLast
              !WRITE(*,"(7ES15.8)") d_omega_ic_dt,LFfac,c_lorentz_ic,lorentz_torque_ic
           ELSE
              rhs(n_r_max)=0.d0
           END IF

           !----- This is the normal RHS for the other radial grid points:
           DO nR=2,n_r_max-1
              rhs(nR)=O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*z(lm1,nR)+ &
                   w1*dzdt(lm1,nR)+ &
                   w2*dzdtLast(lm1,nR)
           END DO
           !WRITE(*,"(A,3ES15.8,A,ES20.13)") "div = ",O_dt,dLh(lm1),SUM(ABS(or2)),", rhs=",SUM(REAL(rhs*CONJG(rhs)))
           !WRITE(*,"(I4,2ES20.13)") n_r_max,rhs(n_r_max)
           CALL cgesl(z10Mat,n_r_max,n_r_max,z10Pivot,rhs)
           !PERFOFF

        ELSE IF ( l1 /= 0 ) THEN
           !PERFON('upZ_ln0')
           IF ( .NOT. lZmat(l1) ) THEN
              !PERFON('upZ_mat')
              CALL get_zMat(dt,l1,hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                   zMat(1,1,l1),zPivot(1,l1))
              !WRITE(*,"(A,2I4,ES25.16)") "zMat = ",lm1,l1,SUM(zMat(:,:,l1))
              lZmat(l1)=.TRUE.
              !PERFOFF
           END If
           lmB=lmB+1
           rhs1(1,lmB)      =0.D0
           rhs1(n_r_max,lmB)=0.D0
           DO nR=2,n_r_max-1
              rhs1(nR,lmB)=O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*z(lm1,nR)+ &
                   w1*dzdt(lm1,nR)+ &
                   w2*dzdtLast(lm1,nR)
           END DO
           !PERFOFF
        END IF
        !IF (lmB.GT.0) THEN
        !   IF (lm1.EQ.242) THEN
        !      DO nR=1,n_r_max
        !         WRITE(*,"(4X,A,2I4,2(I4,F20.16))") "rhs1 = ",lm1, nR,EXPONENT(REAL(rhs1(nR,lmB))),FRACTION(REAL(rhs1(nR,lmB))),&
        !              &EXPONENT(aimag(rhs1(nR,lmB))),FRACTION(aimag(rhs1(nR,lmB)))
        !      END DO
        !   END IF
        !END IF
     END DO

     !PERFON('upZ_sol')
     IF ( lmB > 0 ) THEN
        CALL cgeslML(zMat(1,1,l1),n_r_max,n_r_max, &
             &       zPivot(1,l1),rhs1,n_r_max,lmB)
     END IF
     !PERFOFF
     IF ( lRmsNext ) THEN ! Store old z
        DO nR=1,n_r_max
           DO lm=1,sizeLMB2(nLMB2,nLMB)
              lm1=lm22lm(lm,nLMB2,nLMB)
              workB(lm1,nR)=z(lm1,nR)
           END DO
        END DO
     END IF
     !PERFON('upZ_pp')
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
           
        IF ( l_z10mat .AND. lm1 == l1m0 ) THEN
           DO n_cheb=1,n_cheb_max
              z(lm1,n_cheb)=REAL(rhs(n_cheb))
           END DO
        ELSE IF ( l1 /= 0 ) THEN
           lmB=lmB+1
           !IF (lm1.EQ.242) THEN
           !   DO nR=1,n_r_max
           !      WRITE(*,"(4X,A,2I4,2(I4,F20.16))") "after sol, rhs1 = ",lm1,nR,&
           !           &EXPONENT(REAL(rhs1(nR,lmB))),FRACTION(REAL(rhs1(nR,lmB))),&
           !           &EXPONENT(aimag(rhs1(nR,lmB))),FRACTION(aimag(rhs1(nR,lmB)))
           !   END DO
           !END IF
           IF ( m1 > 0 ) THEN
              DO n_cheb=1,n_cheb_max
                 z(lm1,n_cheb)=rhs1(n_cheb,lmB)
              END DO
           ELSE
              DO n_cheb=1,n_cheb_max
                 z(lm1,n_cheb)= &
                      CMPLX(REAL(rhs1(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
              END DO
           END IF
        END IF
     END DO
     !PERFOFF
  END DO       ! end of loop over lm blocks


  !-- set cheb modes > n_cheb_max to zero (dealiazing)
  DO n_cheb=n_cheb_max+1,n_r_max
     DO lm1=lmStart,lmStop
        z(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END DO

  !IF (llm.LE.242.AND.242.LE.ulm) THEN
  !   DO nR=1,n_r_max
        !IF ((nR.EQ.1).OR.(nR.EQ.5)) THEN
        !DO lm=llm,ulm
  !      WRITE(*,"(4X,A,I3,2ES22.14)") "bef_der: ",nR,z(242,nR)
        !END DO
        !   END IF
        !   WRITE(*,"(A,I3,2ES22.14)") "bef_der: ",nR,get_global_SUM( z(:,nR) )
  !   END DO
  !END IF
  !PERFON('upZ_der')
  !-- Get derivatives:
  CALL costf1(z, ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &      dzdtLast, i_costf_init, d_costf_init)
  CALL get_ddr(z, dz, workA, ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &       n_r_max, n_cheb_max, dzdtLast, workC, &
       &       i_costf_init,d_costf_init,drx,ddrx)
  !PERFOFF
  !PERFON('upZ_icma')
  !--- Update of inner core and mantle rotation:
  IF ( l10 ) THEN
     IF ( l_rot_ma .AND. .NOT. l_SRMA ) THEN
        IF ( ktopv == 1 ) THEN  ! free slip, explicit time stepping of omega !
           omega_ma=O_dt*omega_ma + LFfac/c_moi_ma * &
                (w1*lorentz_torque_ma+w2*lorentz_torque_maLast)
           omega_ma=dt*omega_ma
        ELSE IF ( ktopv == 2 ) THEN ! no slip, omega given by z10
           omega_ma=c_z10_omega_ma*REAL(z(l1m0,n_r_cmb))
        END IF
        omega_ma1=omega_ma
     END IF
     IF ( l_rot_ic .AND. .NOT. l_SRIC ) THEN
        IF ( kbotv == 1 ) THEN  ! free slip, explicit time stepping of omega !
           omega_ic=O_dt*omega_ic + LFfac/c_moi_ic * &
                (w1*lorentz_torque_ic+w2*lorentz_torque_icLast)
           omega_ic=dt*omega_ic
        ELSE IF ( kbotv == 2 ) THEN ! no slip, omega given by z10
           omega_ic=c_z10_omega_ic*REAL(z(l1m0,n_r_icb))
        END IF
        omega_ic1=omega_ic
        !write(*,"(A,I4,A,ES20.13)") "after ic update, nLMB = ",nLMB,", omega_ic = ",omega_ic
     END IF
  END IF  ! l=1,m=0 contained in block ?
  !PERFOFF
  !PERFON('upZ_ang')
  !--- We correct so that the angular moment about axis in the equatorial plane
  !    vanish and the angular moment about the (planetary) rotation axis
  !    is kept constant.
  l1m1=lm2(1,1)
  IF ( l_correct_AMz .AND.  l1m0 > 0 .AND. &
       lmStart_00 <= l1m0 .AND. lmStop >= l1m0 ) THEN

     DO nR=1,n_r_max
        z10(nR)=z(l1m0,nR)
     END DO
     CALL get_angular_moment(z10,z11,omega_ic,omega_ma, &
          angular_moment_oc, &
          angular_moment_ic,angular_moment_ma)
     DO i=1,3
        angular_moment(i)=angular_moment_oc(i) + &
             angular_moment_ic(i) + &
             angular_moment_ma(i)
     END DO
     IF ( ( ktopv == 2 .AND. l_rot_ma ) .AND. &
          ( kbotv == 2 .AND. l_rot_ic ) ) THEN
        nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb + &
             c_moi_ic*c_z10_omega_ic*r_icb*r_icb + &
             c_moi_oc*y10_norm
     ELSE IF ( ktopv == 2 .AND. l_rot_ma ) THEN
        nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb+c_moi_oc*y10_norm
     ELSE IF ( kbotv == 2 .AND. l_rot_ic ) THEN
        nomi=c_moi_ic*c_z10_omega_ic*r_icb*r_icb+c_moi_oc*y10_norm
     ELSE
        nomi=c_moi_oc*y10_norm
     END IF
     corr_l1m0=CMPLX(angular_moment(3)-AMstart,0.d0,KIND=KIND(0d0))/nomi

     !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
     !         derivatives:
     DO nR=1,n_r_max
        r_E_2=r(nR)*r(nR)
        z(l1m0,nR)  =z(l1m0,nR)  - rho0(nR)*r_E_2*corr_l1m0
        dz(l1m0,nR) =dz(l1m0,nR) - rho0(nR)*( &
             2.d0*r(nR)+r_E_2*beta(nR))*corr_l1m0
        workA(l1m0,nR)=workA(l1m0,nR)-rho0(nR)*( &
             2.d0+4.d0*beta(nR)*r(nR) + &
             dbeta(nR)*r_E_2 + &
             beta(nR)*beta(nR)*r_E_2 )*corr_l1m0
     END DO
     IF ( ktopv == 2 .AND. l_rot_ma ) &
          omega_ma=c_z10_omega_ma*REAL(z(l1m0,n_r_cmb))
     IF ( kbotv == 2 .AND. l_rot_ic ) &
          omega_ic=c_z10_omega_ic*REAL(z(l1m0,n_r_icb))
     omega_ic1=omega_ic
     omega_ma1=omega_ma
     !write(*,"(A,I4,A,ES20.13)") "after 2nd ic update, nLMB = ",nLMB,", omega_ic = ",omega_ic

  END IF ! l=1,m=0 contained in lm-block ?

  IF ( l_correct_AMe .AND.  l1m1 > 0 .AND. &
       lmStart_00 <= l1m1 .AND. lmStop >= l1m1 ) THEN

     DO nR=1,n_r_max
        z11(nR)=z(l1m1,nR)
     END DO
     CALL get_angular_moment(z10,z11,omega_ic,omega_ma, &
          angular_moment_oc, &
          angular_moment_ic,angular_moment_ma)
     DO i=1,3
        angular_moment(i)=angular_moment_oc(i) + &
             angular_moment_ic(i) + &
             angular_moment_ma(i)
     END DO
     corr_l1m1=CMPLX(angular_moment(1),-angular_moment(2),KIND=KIND(0d0)) / &
          (2.d0*y11_norm*c_moi_oc)

     !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
     !         derivatives:
     DO nR=1,n_r_max
        r_E_2=r(nR)*r(nR)
        z(l1m1,nR)  =z(l1m1,nR)  -  rho0(nR)*r_E_2*corr_l1m1
        dz(l1m1,nR) =dz(l1m1,nR) -  rho0(nR)*( &
             2.d0*r(nR)+r_E_2*beta(nR))*corr_l1m1
        workA(l1m1,nR)=workA(l1m1,nR)-rho0(nR)*( &
             2.d0+4.d0*beta(nR)*r(nR) + &
             dbeta(nR)*r_E_2 + &
             beta(nR)*beta(nR)*r_E_2 )*corr_l1m1
     END DO

  END IF ! l=1,m=1 contained in lm-block ?

  IF (DEBUG_OUTPUT) THEN
     DO nR=1,n_r_max
        IF ((nR.EQ.1).OR.(nR.EQ.5)) THEN
           DO lm=llm,ulm
              WRITE(*,"(4X,A,2I3,4ES22.14)") "upZ_new: ",nR,lm,z(lm,nR),dz(lm,nR)
           END DO
        END IF
        WRITE(*,"(A,I3,4ES22.14)") "upZ_new: ",nR,get_global_SUM( z(:,nR) ),get_global_SUM( dz(:,nR) )
     END DO
  END IF
  !-- Calculate explicit time step part:
  DO nR=n_r_cmb+1,n_r_icb-1
     !WRITE(*,"(A,I4,5ES20.12)") "r-dependent : ",nR,dLvisc(nR),beta(nR),or1(nR),or2(nR),dbeta(nR)
     DO lm1=lmStart_00,lmStop
        Dif(lm1)=hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1)))*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*visc(nR)* &
             & ( workA(lm1,nR) &
             &   +(dLvisc(nR)-beta(nR))*dz(lm1,nR) &
             &   -( dLvisc(nR)*beta(nR) &
             &      + 2.d0*dLvisc(nR)*or1(nR) &
             &      + dLh(st_map%lm2(lm2l(lm1),lm2m(lm1))) * or2(nR)&
             &      + dbeta(nR)&
             &      + 2.d0*beta(nR)*or1(nR) &
             &    ) * z(lm1,nR) &
             & )

!        IF (nR.EQ.2) THEN
!           WRITE(*,"(2I4,8ES20.12)") nR,lm1,workA(lm1,nR),Dif(lm1),z(lm1,nR),dz(lm1,nR)
!        END IF
        dzdtLast(lm1,nR)=dzdt(lm1,nR)-coex*Dif(lm1)
        IF ( lRmsNext ) THEN
           workB(lm1,nR)= &
                O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*(z(lm1,nR)-workB(lm1,nR))
           IF ( l_RMStest ) &
                workB(lm1,nR)= workB(lm1,nR)-Dif(lm1)
        END IF
     END DO
     IF ( lRmsNext ) THEN
        CALL hInt2Tor(Dif,1,lm_max,nR,lmStart_00,lmStop, &
             DifTor2hInt(nR,nTh),DifTorAs2hInt(nR,nTh),lo_map)
        CALL hInt2Tor(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
             dtVTor2hInt(nR,nTh),dtVTorAs2hInt(nR,nTh),lo_map)
     END IF
  END DO
  !--- Note: from ddz=workA only the axisymmetric contributions are needed
  !    beyond this point for the TO calculation.
  !    Parallization note: Very likely, all axisymmetric modes m=0 are
  !    located on the first processor #0.
  IF ( l_TO ) THEN
     DO nR=1,n_r_max
        DO lm1=lmStart_00,lmStop
           l1=lm2l(lm1)
           m1=lm2m(lm1)
           IF ( m1 == 0 ) ddzASL(l1+1,nR)=REAL(workA(lm1,nR))
        END DO
     END DO
  END IF

  !----- Special thing for l=1,m=0 for rigid boundaries and
  !      if IC or mantle are allowed to rotate:
  IF ( l10 .AND. l_z10mat ) THEN ! z10 term !
     lm1=lm2(1,0)
     !----- NOTE opposite sign of visouse torque on ICB and CMB:
     IF ( .NOT. l_SRMA .AND. ktopv == 2 .AND. l_rot_ma ) THEN
        d_omega_ma_dtLast=d_omega_ma_dt -            &
             coex * ( 2.D0*or1(1)*REAL(z(lm1,1))  - &
             REAL(dz(lm1,1)) )
     END IF
     IF ( .NOT. l_SRIC .AND. kbotv == 2 .AND. l_rot_ic ) THEn
        d_omega_ic_dtLast=d_omega_ic_dt +                        &
             coex * ( 2.D0*or1(n_r_max)*REAL(z(lm1,n_r_max))  - &
             REAL(dz(lm1,n_r_max)) )
     END IF
  END IF
  !PERFOFF


  RETURN
end SUBROUTINE updateZ


!--------------------------------------------------------------------
!-- End of subroutine update_z
!--------------------------------------------------------------------
