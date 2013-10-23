!$Id$
#include "perflib_preproc.cpp"
!***********************************************************************
SUBROUTINE updateB(b,db,ddb,aj,dj,ddj,dVxBhLM, &
     &             dbdt,dbdtLast,djdt,djdtLast, &
     &             b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic, &
     &             dbdt_icLast,djdt_icLast, &
     &             b_nl_cmb,aj_nl_cmb,aj_nl_icb,omega_ic, &
     &             w1,coex,dt,time,nLMB,lRmsNext,nTh)
  !***********************************************************************

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Calculated update of magnetic field potential and the time       |
  !  |  stepping arrays dbdtLast, ...                                    |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  !  updates the magnetic field potentials b, aj and
  !  their derivatives,
  !  adds explicit part to time derivatives of b and j

  !  input:  w1 - weight for dbdt-contribution from current time step
  !               (w2=1-w1: weight for contrib. from previous step)
  !          coex - factor depending on weighting alpha of
  !                 implicit contribution
  !          mc_min,mc_max
  !               - range of mca-indices in which field is updated
  !                 (harmonic order m=(mca-1)*minc)
  !          b_nl_cmb = RHS of nonlinear BC for poloidal magnetic field
  !                      potential in the case of stress free CMB
  !          aj_nl_cmb = RHS of nonlinear BC for toroidal magnetic field
  !                      potential in the case of stress free CMB
  !          aj_nl_icb = RHS of nonlinear BC for radial derivative of
  !                      toroidal magnetic field
  !                      potential in the case of stress free ICB

  !-----------------------------------------------------------------------

  USE truncation
  USE radial_functions, ONLY: i_costf_init,d_costf_init,drx,ddrx,or2,r_cmb,&
       & i_costf1_ic_init,d_costf1_ic_init,&
       & i_costf2_ic_init,d_costf2_ic_init,&
       & dr_fac_ic,lambda,dLlambda,o_r_ic
  USE radial_data,ONLY: n_r_cmb,n_r_icb
  USE physical_parameters
  USE num_param
  USE init_fields
  USE blocking,ONLY: nLMBs,st_map,lo_map,st_sub_map,lo_sub_map,lmStartB,lmStopB
  USE horizontal_data
  USE logic
  USE matrices
  USE RMS
  USE const
  USE Bext
  USE algebra, ONLY: cgeslML
  USE LMLoop_data, ONLY: llmMag,ulmMag,llm_realMag,ulm_realMag
  USE parallel_mod, only:rank
  IMPLICIT NONE

  !-- Input/output of scalar potentials and time stepping arrays:
  COMPLEX(kind=8),INTENT(INOUT) :: b(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: db(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: ddb(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(INOUT) :: aj(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: dj(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: ddj(llmMag:ulmMag,n_r_maxMag)

  COMPLEX(kind=8),INTENT(IN) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)

  COMPLEX(kind=8),INTENT(IN) :: dbdt(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(INOUT) :: dbdtLast(llmMag:ulmMag,n_r_maxMag)

  COMPLEX(kind=8),INTENT(INOUT) :: djdt(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(INOUT) :: djdtLast(llmMag:ulmMag,n_r_maxMag)
  
  COMPLEX(kind=8),INTENT(INOUT) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(INOUT) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(OUT) :: ddj_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(INOUT) :: dbdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(INOUT) :: djdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)

  REAL(kind=8),INTENT(IN) :: omega_ic
  !-- Output:
  !       UPDATED b,db,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic,
  !       dbdtLast,djdtLast,dbdt_icLast,djdt_icLast

  !-- Input of non-linear magnetic boundary conditions:
  COMPLEX(kind=8),INTENT(IN) :: b_nl_cmb(*)  ! nonlinear BC for b at CMB
  COMPLEX(kind=8),INTENT(IN) :: aj_nl_cmb(*) ! nonlinear BC for aj at CMB
  COMPLEX(kind=8),INTENT(IN) :: aj_nl_icb(*) ! nonlinear BC for dr aj at ICB

  !-- Input of other variables
  REAL(kind=8),INTENT(IN) :: w1    ! weight for time step !
  REAL(kind=8),INTENT(IN) :: coex  ! factor depending on alpha
  REAL(kind=8),INTENT(IN) :: dt
  REAL(kind=8),INTENT(IN) :: time
  INTEGER,INTENT(IN) :: nLMB
  LOGICAL,INTENT(IN) :: lRmsNext
  INTEGER,INTENT(IN) :: nTh

  !-- Local work arrays:
  COMPLEX(kind=8) :: workA(llmMag:ulmMag,n_r_max)
  COMPLEX(kind=8) :: workB(llmMag:ulmMag,n_r_max)

  !-- Local variables:
  REAL(kind=8) :: w2             ! weight of second time step
  REAL(kind=8) :: O_dt
  INTEGER :: l1,m1               ! degree and order
  INTEGER :: lm1,lm,lmB          ! position of (l,m) in array
  INTEGER :: lmStart,lmStop      ! max and min number of orders m
  INTEGER :: lmStart_real        ! range of lm for real array
  INTEGER :: lmStop_real
  INTEGER :: lmStart_00          ! excluding l=0,m=0
  INTEGER :: nLMB2
  INTEGER :: n_cheb              ! No of cheb polynome (degree+1)
  INTEGER :: nR                 ! No of radial grid point
  INTEGER :: n_r_real            ! total number of used grid points

  COMPLEX(kind=8) :: fac
  COMPLEX(kind=8) :: dbdt_ic,djdt_ic  ! they are calculated here !

  !-- right hand sides:
  COMPLEX(kind=8) :: rhs1(2*n_r_max,lo_sub_map%sizeLMB2max)
  COMPLEX(kind=8) :: rhs2(2*n_r_max,lo_sub_map%sizeLMB2max)

  INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
  INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
  INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

  !-- for feedback
  REAL(kind=8) :: ff,cimp,aimp,b10max

  REAL(kind=8) :: direction
  SAVE direction

  !-- end of declaration
  !-----------------------------------------------------------------------

  IF ( .NOT. l_update_b ) RETURN

  nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
  sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
  lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
  lm22l(1:,1:,1:) => lo_sub_map%lm22l
  lm22m(1:,1:,1:) => lo_sub_map%lm22m
  lm2(0:,0:) => lo_map%lm2
  lm2l(1:lm_max) => lo_map%lm2l
  lm2m(1:lm_max) => lo_map%lm2m

  n_r_real=n_r_max
  if ( l_cond_ic ) n_r_real=n_r_max+n_r_ic_max

  lmStart     =lmStartB(nLMB)
  lmStop      =lmStopB(nLMB)
  lmStart_00  =MAX(2,lmStart)
  lmStart_real=2*lmStart_00-1
  lmStop_real =2*lmStop

  ! output the input 
  !WRITE(*,"(4(A,I4),I4,A,I4)") "nLMB=",nLMB,", from ",lmStart," to ",lmStop,&
  !     &", reals: ",lmStart_real,lmStop_real,", nLMBs2 = ",nLMBs2(nLMB)

  w2  =1.D0-w1
  O_dt=1.D0/dt

  !--- Start with finishing djdt:
  !    dVxBhLM is still in the R-distributed space,
  !    the ouput workA is in the LM-distributed space.
  !IF (2*lmStart-1 - llm_realMag+1.NE.1) THEN
  !   WRITE(*,"(I4,2(A,I6))") rank,": lmStart = ",lmStart,", llm_realMag = ",llm_realMag
  !   STOP
  !END IF
  !IF (lmStop_real .NE. ulm_realMag) THEN
  !   WRITE(*,"(I4,A,2I6)") rank,": ",ulm_realMag,lmStop_real
  !   stop
  !END IF
  !CALL get_drNS( dVxBhLM,workA, &
  !     &         ulm_realMag-llm_realMag+1,(2*lmStart-1)-llm_realMag+1,lmStop_real-llm_realMag+1, &
  !     &         n_r_max,n_cheb_max,workB, &
  !     &         i_costf_init,d_costf_init,drx)
  ! simplified interface
  !PRINT*,rank,": computing for ",ulm_realMag-llm_realMag+1," rows, i_costf_init = ",i_costf_init

  !PERFON('drNS')
  CALL get_drNS( dVxBhLM,workA, &
       &         ulm_realMag-llm_realMag+1,1,ulm_realMag-llm_realMag+1, &
       &         n_r_max,n_cheb_max,workB, &
       &         i_costf_init,d_costf_init,drx)
  !PERFOFF

  DO nR=1,n_r_max
      DO lm=lmStart_00,lmStop
        djdt(lm,nR)=djdt(lm,nR)+or2(nR)*workA(lm,nR)
     END DO
   END DO
 
  DO nLMB2=1,nLMBs2(nLMB)
     !WRITE(*,"(2(A,I4))") "nLMB2 = ",nLMB2,", sizeLMB2 = ",sizeLMB2(nLMB2,nLMB)
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
        IF ( l1 > 0 ) THEN
           lmB=lmB+1

           IF ( .NOT. lBmat(l1) ) THEN
              !PERFON('get_bMat')
              CALL get_bMat(dt,l1,hdif_B(st_map%lm2(l1,m1)), &
                   bMat(1,1,l1),bPivot(1,l1), &
                   jMat(1,1,l1),jPivot(1,l1))
              lBmat(l1)=.TRUE.
              !PERFOFF
           END IF

           !-------- Magnetic boundary conditions, outer core:
           !         Note: the CMB condition is not correct if we assume free slip
           !         and a conducting mantle (conductance_ma>0).
           IF ( l_b_nl_cmb ) THEN ! finitely conducting mantle
              rhs1(1,lmB) =  b_nl_cmb(st_map%lm2(l1,m1))
              rhs2(1,lmB) = aj_nl_cmb(st_map%lm2(l1,m1))
           ELSE
              rhs1(1,lmB) = 0.D0
              rhs2(1,lmB) = 0.D0
           END IF

           rhs1(n_r_max,lmB)=0.D0
           IF ( kbotb == 2 ) rhs1(n_r_max-1,lmB)=0.D0

           rhs2(n_r_max,lmB)=0.D0
           IF ( m1 == 0 ) THEN   ! Magnetoconvection boundary conditions
              IF ( imagcon /= 0 .AND. tmagcon <= time ) THEN
                 IF ( l1 == 2 .AND. imagcon > 0 .AND. imagcon .NE. 12 ) THEN
                    rhs2(1,lmB)      =CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))
                    rhs2(n_r_max,lmB)=CMPLX(bpeakbot,0.D0,KIND=KIND(0d0))
                 ELSE IF( l1 == 1 .AND. imagcon == 12 ) THEN
                    rhs2(1,lmB)      =CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))
                    rhs2(n_r_max,lmB)=CMPLX(bpeakbot,0.D0,KIND=KIND(0d0))
                 ELSE IF( l1 == 1 .AND. imagcon == -1) THEN
                    rhs1(n_r_max,lmB)=CMPLX(bpeakbot,0.D0,KIND=KIND(0d0))
                 ELSE IF( l1 == 1 .AND. imagcon == -2) THEN
                    rhs1(1,lmB)      =CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))
                 ELSE IF( l1 == 3 .AND. imagcon == -10 ) THEN
                    rhs2(1,lmB)      =CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))
                    rhs2(n_r_max,lmB)=CMPLX(bpeakbot,0.D0,KIND=KIND(0d0))
                 END IF
              END IF
              IF ( n_imp > 1 .AND. l1 == 1 ) THEN
                 IF ( n_imp == 2 ) THEN
                    !  Chose external field coefficient so that amp_imp is the amplitude of
                    !  the external dipole field:
                    bpeaktop=3.D0/2.D0*r_cmb/y10_norm*amp_imp
                    rhs1(1,lmB)=CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))
                 ELSE IF ( n_imp == 3 ) THEN
                    !  Chose external field coefficient so that amp_imp is the amplitude of
                    !  the external dipole field:
                    bpeaktop=3.D0/2.D0*r_cmb/y10_norm*amp_imp
                    IF ( REAL(b(2,1)) > 1.D-9 ) &
                         direction=REAL(b(2,1))/DABS(REAL(b(2,1)))
                    rhs1(1,lmB)=CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))*direction
                 ELSE IF ( n_imp == 4 ) THEN
                    !  I have for gotten what this was supposed to do:
                    bpeaktop=3.D0/r_cmb*amp_imp*REAL(b(2,1))**2
                    rhs1(1,lmB)=CMPLX(bpeaktop,0.D0,KIND=KIND(0d0))/b(2,1)

                 ELSE

                    ! Special Heyner feedback functions:
                    ! We don't provide an actual external field strength but rather its
                    ! dependence on the internal dipole via a feedback function ff:
                    !              b10_external=ff(b10_internal)
                    ! where b10_external is the external axial dipole field coeff. and
                    ! b10_internal=b(2,1) is the internal axial dipole field coeff.
                    ! Then rhs1 is always given by:
                    !              rhs1 = (2*l+1)/r_cmb * ff
                    ! because of the special formulation of the CMB matching condition!
                    ! Note that
                    !  B_r_internal(r_cmb) = l*(l+1)/r_cmb**2 * b10_internal*y10_norm*cos(theta)
                    !  B_r_external(r)     = l*(l+1)/r_cmb**2 * b10_external*y10_norm*cos(theta)
                    ! This determins the units of ff=b10_external.
                    ! B itself is given in units sqrt(rho*omega/sigma) so that B**2 is
                    ! the Elsasser number. ff is thus given in units L**2*sqrt(rho*omega/sigma)
                    ! with L=(r_cmb-r_icb). Note that the external dipole field does not depend
                    ! on the radius.
                    ! Here amp_imp provides the relative amplitude so that the maximum of
                    ! the external field is   max(b10_external/b10_internal)=amp_imp

                    ff=0.D0
                    IF ( n_imp == 7 ) THEN
                       ! Using a feedback function of the form
                       !    ff= aimp * b10_internal**expo_imp / (cimp+b10_internal**expo_imp)
                       ! Here expo_imp is an input parameter and aimp and cimp are determined
                       ! such that the maximum of b10_external/b10_internal is amp_imp
                       ! and is located at b10_internal=bmax_imp.
                       ! amp_imp and bmax_imp are two more input parameters.
                       ! Note that bmax_imp on input hat the dimensionless unit of the magnetic
                       ! field B and we convert is to dimensionless units for dipole
                       ! coefficients b10 first by multiplying with r_cmb**2
                       ! Other factors like  l*(l+1)*y10_norm*mean(cos(theta)) may also be
                       ! considered, but like r_cmb**2 they are all of order one here.
                       ! Since amp_imp is dimensionless aimp and thus ff has the dimension of
                       ! b10 as required.
                       b10max=bmax_imp*r_cmb**2
                       cimp=b10max**expo_imp / (expo_imp-1)
                       aimp=amp_imp*cimp*expo_imp / &
                            (cimp*(expo_imp-1))**((expo_imp-1)/expo_imp)

                       ff=  aimp*REAL(b(2,1))**expo_imp/ &
                            (cimp+REAL(b(2,1))**expo_imp)

                    END IF
                    rhs1(1,lmB)=(2*l1+1)/r_cmb*ff

                 END IF
              END IF
           END IF

           DO nR=2,n_r_max-1
              rhs1(nR,lmB)= ( w1*dbdt(lm1,nR) + w2*dbdtLast(lm1,nR) ) &
                   & + O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*b(lm1,nR)
              rhs2(nR,lmB)= ( w1*djdt(lm1,nR) + w2*djdtLast(lm1,nR) ) &
                   & + O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*aj(lm1,nR)
           END DO

           !-------- Magnetic boundary conditions, inner core for radial derivatives
           !         of poloidal and toroidal magnetic potentials:
           IF ( l_cond_ic ) then    ! inner core
              rhs1(n_r_max+1,lmB)=0.d0
              IF ( l_b_nl_icb ) THEN
                 rhs2(n_r_max+1,lmB)=aj_nl_icb(st_map%lm2(l1,m1))
              ELSE
                 rhs2(n_r_max+1,lmB)=0.d0
              END IF

              DO nR=2,n_r_ic_max
                 IF ( omega_ic == 0.D0 .OR. .NOT. l_rot_ic .OR. &
                      .NOT. l_mag_nl ) THEN
                    dbdt_ic=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    djdt_ic=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 ELSE
                    fac=-omega_ic*or2(n_r_max)*dPhi(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))
                    dbdt_ic=fac*b_ic(lm1,nR)
                    djdt_ic=fac*aj_ic(lm1,nR)
                 END IF
                 rhs1(n_r_max+nR,lmB)=            &
                      & ( w1*dbdt_ic + w2*dbdt_icLast(lm1,nR) ) &
                      & +O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * b_ic(lm1,nR)
                 rhs2(n_r_max+nR,lmB)=            &
                      & ( w1*djdt_ic + w2*djdt_icLast(lm1,nR) ) &
                      & +O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * aj_ic(lm1,nR)

                 !--------- Store the IC non-linear terms for the usage below:
                 dbdt_icLast(lm1,nR)=dbdt_ic
                 djdt_icLast(lm1,nR)=djdt_ic
              END DO
           END IF

        END IF ! l>0
        !IF (lmB.GT.0) WRITE(*,"(2I4,4ES22.14)") l1,m1,SUM( rhs1(1:n_r_max+n_r_ic_max,lmB) ),&
        !     & SUM( rhs2(1:n_r_max+n_r_ic_max,lmB) )
     END DO    ! loop over lm in block
     IF ( lmB > 0 ) THEN
        !PERFON('bMat_sol')
        CALL cgeslML(bMat(1,1,l1),n_r_tot,n_r_real, &
             bPivot(1,l1),rhs1,2*n_r_max,lmB)
        CALL cgeslML(jMat(1,1,l1),n_r_tot,n_r_real, &
             jPivot(1,l1),rhs2,2*n_r_max,lmB)
        !PERFOFF
     END IF

     IF ( lRmsNext ) THEN ! Store old b,aj
        DO nR=1,n_r_max
           DO lm=1,sizeLMB2(nLMB2,nLMB)
              lm1=lm22lm(lm,nLMB2,nLMB)
              workA(lm1,nR)=b(lm1,nR)
              workB(lm1,nR)=aj(lm1,nR)
           END DO
        END DO
     END IF

     !----- Update magnetic field in cheb space:
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)

        IF ( l1 > 0 ) THEN
           lmB=lmB+1

           IF ( m1 > 0 ) THEN
              DO n_cheb=1,n_cheb_max  ! outer core
                 b(lm1,n_cheb) =rhs1(n_cheb,lmB)
                 aj(lm1,n_cheb)=rhs2(n_cheb,lmB)
              END DO
              IF ( l_cond_ic ) then   ! inner core
                 DO n_cheb=1,n_cheb_ic_max
                    b_ic(lm1,n_cheb) = &
                         rhs1(n_r_max+n_cheb,lmB)
                    aj_ic(lm1,n_cheb)= &
                         rhs2(n_r_max+n_cheb,lmB)
                 END DO
              END IF
           ELSE
              DO n_cheb=1,n_cheb_max   ! outer core
                 b(lm1,n_cheb) = &
                      CMPLX(REAL(rhs1(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                 aj(lm1,n_cheb)= &
                      CMPLX(REAL(rhs2(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
              END DO
              IF ( l_cond_ic ) THEN    ! inner core
                 DO n_cheb=1,n_cheb_ic_max
                    b_ic(lm1,n_cheb)= &
                         CMPLX(REAL(rhs1(n_r_max+n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                    aj_ic(lm1,n_cheb)= &
                         CMPLX(REAL(rhs2(n_r_max+n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                 END DO
              END IF
           END IF

        END IF
     END DO

  END DO      ! end of do loop over lm1

  !-- Set cheb modes > n_cheb_max to zero (dealiazing)
  !   for inner core modes > 2*n_cheb_ic_max = 0
  DO n_cheb=n_cheb_max+1,n_r_max
     DO lm1=lmStart_00,lmStop
        b(lm1,n_cheb) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        aj(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END DO
  IF ( l_cond_ic ) THEN
     DO n_cheb=n_cheb_ic_max+1,n_r_ic_max
        DO lm1=lmStart_00,lmStop
           b_ic(lm1,n_cheb) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           aj_ic(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
     END DO
  END IF

  !PERFON('rad_der')
  !-- Radial derivatives: dbdtLast and djdtLast used as work arrays
  CALL costf1(b,ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
       &      dbdtLast,i_costf_init,d_costf_init)
  CALL get_ddr(b,db,ddb,ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
       &       n_r_max,n_cheb_max,dbdtLast,djdtLast, &
       &       i_costf_init,d_costf_init,drx,ddrx)
  CALL costf1(aj,ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
       &      dbdtLast,i_costf_init,d_costf_init)
  CALL get_ddr(aj,dj,ddj,ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
       &       n_r_max,n_cheb_max,dbdtLast,djdtLast, &
       &       i_costf_init,d_costf_init,drx,ddrx)
  !CALL costf1(b,lm_max_real,lmStart_real,lmStop_real, &
  !     dbdtLast,i_costf_init,d_costf_init)
  !CALL get_ddr(b,db,ddb,lm_max_real,lmStart_real,lmStop_real, &
  !     n_r_max,n_cheb_max,dbdtLast,djdtLast, &
  !     i_costf_init,d_costf_init,drx,ddrx)
  !CALL costf1(aj,lm_max_real,lmStart_real,lmStop_real, &
  !     dbdtLast,i_costf_init,d_costf_init)
  !CALL get_ddr(aj,dj,ddj,lm_max_real,lmStart_real,lmStop_real, &
  !     n_r_max,n_cheb_max,dbdtLast,djdtLast, &
  !     i_costf_init,d_costf_init,drx,ddrx)

  !-- Same for inner core:
  IF ( l_cond_ic ) THEN
     CALL costf1(b_ic,ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
          dbdtLast,i_costf1_ic_init,d_costf1_ic_init)
     CALL get_ddr_even( b_ic,db_ic,ddb_ic, &
          & ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
          & n_r_ic_max,n_cheb_ic_max,dr_fac_ic,dbdtLast,djdtLast, &
          & i_costf1_ic_init,d_costf1_ic_init, &
          & i_costf2_ic_init,d_costf2_ic_init)
     CALL costf1(aj_ic,ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
          & dbdtLast,i_costf1_ic_init,d_costf1_ic_init)
     CALL get_ddr_even( aj_ic,dj_ic,ddj_ic, &
          & ulm_realMag-llm_realMag+1,lmStart_real-llm_realMag+1,lmStop_real-llm_realMag+1, &
          & n_r_ic_max,n_cheb_ic_max,dr_fac_ic,dbdtLast,djdtLast, &
          & i_costf1_ic_init,d_costf1_ic_init, &
          & i_costf2_ic_init,d_costf2_ic_init)
     !CALL costf1(b_ic,lm_max_real,lmStart_real,lmStop_real, &
     !     dbdtLast,i_costf1_ic_init,d_costf1_ic_init)
     !CALL get_ddr_even( b_ic,db_ic,ddb_ic, &
     !     & lm_max_real,lmStart_real,lmStop_real, &
     !     & n_r_ic_max,n_cheb_ic_max,dr_fac_ic,dbdtLast,djdtLast, &
     !     & i_costf1_ic_init,d_costf1_ic_init, &
     !     & i_costf2_ic_init,d_costf2_ic_init)
     !CALL costf1(aj_ic,lm_max_real,lmStart_real,lmStop_real, &
     !     & dbdtLast,i_costf1_ic_init,d_costf1_ic_init)
     !CALL get_ddr_even( aj_ic,dj_ic,ddj_ic, &
     !     & lm_max_real,lmStart_real,lmStop_real, &
     !     & n_r_ic_max,n_cheb_ic_max,dr_fac_ic,dbdtLast,djdtLast, &
     !     & i_costf1_ic_init,d_costf1_ic_init, &
     !     & i_costf2_ic_init,d_costf2_ic_init)
  END IF
  !PERFOFF
  !-- We are now back in radial space !

  DO nR=n_r_cmb+1,n_r_icb-1
     DO lm1=lmStart_00,lmStop
        l1=lm2l(lm1)
        m1=lm2m(lm1)
        dbdtLast(lm1,nR)= dbdt(lm1,nR) - &
             coex*opm*lambda(nR)*hdif_B(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*or2(nR) * &
             ( ddb(lm1,nR) - dLh(st_map%lm2(l1,m1))*or2(nR)*b(lm1,nR) )
        djdtLast(lm1,nR)= djdt(lm1,nR) - &
             coex*opm*lambda(nR)*hdif_B(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*or2(nR) * &
             ( ddj(lm1,nR) + dLlambda(nR)*dj(lm1,nR) - dLh(st_map%lm2(l1,m1))*or2(nR)*aj(lm1,nR) )
        IF ( lRmsNext ) THEN
           workA(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * (  b(lm1,nR)-workA(lm1,nR) )
           workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * ( aj(lm1,nR)-workB(lm1,nR) )
        END IF
     END DO
     IF ( lRmsNext ) THEN
        !WRITE(*,"(A,2I3,2ES20.12)") "workA = ",nLMB,nR,SUM( workA(lmStart_00:lmStop,nR) )
        !CALL hInt2Pol(workA(1,nR),nR,lmStart_00,lmStop,dtBPolLMr, &
        !     dtBPol2hInt(nR,nTh),dtBPolAs2hInt(nR,nTh),lo_map)
        !WRITE(*,"(A,2I3,ES20.13)") "upB: before dtBPol2hInt = ",nLMB,nR,dtBPol2hInt(nR,1)
        !CALL hInt2Pol(workA(llmMag,nR),ulmMag-llmMag+1,nR,lmStart_00-llmMag+1,lmStop-llmMag+1,&
        !     &dtBPolLMr,dtBPol2hInt(nR,nTh),dtBPolAs2hInt(nR,nTh),lo_map)
        CALL hInt2Pol(workA(llmMag,nR),llmMag,ulmMag,nR,lmStart_00,lmStop,&
             &dtBPolLMr,dtBPol2hInt(nR,nTh),dtBPolAs2hInt(nR,nTh),lo_map)
        !WRITE(*,"(A,2I3,ES20.13)") "upB: after  dtBPol2hInt = ",nLMB,nR,dtBPol2hInt(nR,1)

        !WRITE(*,"(A,2I3,3ES20.13)") "upB: before dtBTor2hInt = ",nLMB,nR,dtBTor2hInt(nR,1),SUM( workB(lmStart:lmStop,nR) )
        !CALL hInt2Tor(workB(llmMag,nR),ulmMag-llmMag+1,nR,lmStart_00-llmMag+1,lmStop-llmMag+1, &
        !     &        dtBTor2hInt(nR,nTh),dtBTorAs2hInt(nR,nTh))
        CALL hInt2Tor(workB(llmMag,nR),llmMag,ulmMag,nR,lmStart_00,lmStop, &
             &        dtBTor2hInt(nR,nTh),dtBTorAs2hInt(nR,nTh),lo_map)
        !WRITE(*,"(A,2I3,ES20.13)") "upB: after  dtBTor2hInt = ",nLMB,nR,dtBTor2hInt(nR,1)
     END IF
  END DO


  !----- equations for inner core are different:
  !      D_lP1(lm1)=l+1, O_sr=sigma/sigma_ic
  !      NOTE: no hyperdiffusion in inner core !
  IF ( l_cond_ic ) THEN
     DO nR=2,n_r_ic_max-1
        DO lm1=lmStart_00,lmStop
           l1=lm2l(lm1)
           m1=lm2m(lm1)
           dbdt_icLast(lm1,nR)=dbdt_icLast(lm1,nR) - &
                coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
                (                        ddb_ic(lm1,nR) + &
                2.d0*D_lP1(st_map%lm2(l1,m1))*O_r_ic(nR)*db_ic(lm1,nR) )
           djdt_icLast(lm1,nR)=djdt_icLast(lm1,nR) - &
                coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
                (                        ddj_ic(lm1,nR) + &
                2.d0*D_lP1(st_map%lm2(l1,m1))*O_r_ic(nR)*dj_ic(lm1,nR) )
        END DO
     END DO
     nR=n_r_ic_max
     DO lm1=lmStart_00,lmStop
        l1=lm2l(lm1)
        m1=lm2m(lm1)
        dbdt_icLast(lm1,nR)=dbdt_icLast(lm1,nR) - &
             coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
             (1.d0+2.d0*D_lP1(st_map%lm2(l1,m1)))*ddb_ic(lm1,nR)
        djdt_icLast(lm1,nR)=djdt_icLast(lm1,nR) - &
             coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
             (1.d0+2.d0*D_lP1(st_map%lm2(l1,m1)))*ddj_ic(lm1,nR)
     END DO
  END IF


  RETURN
end SUBROUTINE updateB

!-------------------------------------------------------------------------
!-- End of subroutine update_b
!-------------------------------------------------------------------------
