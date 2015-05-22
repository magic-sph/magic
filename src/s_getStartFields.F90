!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"
SUBROUTINE getStartFields(time,dt,dtNew,n_time_step)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/15/02  by JW. --------------!

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to initialize the fields and       |
  !  |  other auxiliary parameters.                                      |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  use mpi
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE init_fields
  USE Grenoble
  USE blocking
  USE horizontal_data
  USE logic
  USE fields,ONLY: w,dw,ddw,z,dz,s,ds,p,dp,b,db,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic,&
       &omega_ic,omega_ma,&
       &w_LMloc,p_LMloc,s_LMloc,b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc,&
       &ds_LMloc,dp_LMloc,dw_LMloc,ddw_LMloc,db_LMloc,dj_LMloc,ddb_LMloc,&
       &ddj_LMloc,db_ic_LMloc,dj_ic_LMloc,ddb_ic_LMloc,ddj_ic_LMloc,z_LMloc,dz_LMloc,&
       &s_Rloc,ds_Rloc,z_Rloc,dz_Rloc,w_Rloc,dw_Rloc,ddw_Rloc,p_Rloc,dp_Rloc,&
       &b_Rloc,db_Rloc,ddb_Rloc,aj_Rloc,dj_Rloc,&
       & w_LMloc_container,w_Rloc_container,&
       & s_LMloc_container,s_Rloc_container,&
       & z_LMloc_container,z_Rloc_container,&
       & p_LMloc_container,p_Rloc_container,&
       & b_LMloc_container,b_Rloc_container,&
       & aj_LMloc_container,aj_Rloc_container
  USE fieldsLast
  USE output_data
  USE const
  USE usefull, ONLY: cc2real
  USE LMLoop_data,ONLY: lm_per_rank,lm_on_last_rank,llm_realMag,ulm_realMag,llm_real,ulm_real
  USE parallel_mod,ONLY: rank,n_procs
  USE communications, ONLY: lo2r_redist_start,& !lo2r_redist,&
       & lo2r_s,lo2r_z, lo2r_p,&
       & lo2r_b, lo2r_aj, scatter_from_rank0_to_lo,&
       &get_global_sum, lo2r_w
  IMPLICIT NONE

  !---- Output variables:
  REAL(kind=8) :: time,dt,dtNew
  INTEGER :: n_time_step

  !-- Local variables:
  INTEGER :: nR,l1m0,nLMB,l,m
  INTEGER :: lm
  INTEGER :: lmStart,lmStop,lmStartReal,lmStopReal
  REAL(kind=8) :: coex
  REAL(kind=8) :: d_omega_ma_dt,d_omega_ic_dt
  CHARACTER(len=76) :: message

  REAL(kind=8) :: sEA,sES,sAA

  REAL(kind=8) :: s0(n_r_max),ds0(n_r_max)
  REAL(kind=8) :: w1(n_r_max),w2(n_r_max)

  !COMPLEX(kind=8) :: workA_LMloc(llm:ulm,n_r_max)
  !COMPLEX(kind=8) :: workB_LMloc(llm:ulm,n_r_max)
  COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: workA_LMloc,workB_LMloc
  !COMPLEX(kind=8) :: temp_lo(lm_max,n_r_max)
  !COMPLEX(kind=8) :: temp_lo(lm_max)

  INTEGER :: ierr
  logical :: DEBUG_OUTPUT=.false.
  !-- end of declaration
  !---------------------------------------------------------------
  !PERFON('getFlds')
  !print*,"Starting getStartFields"
  !WRITE(*,"(2(A,L1))") "l_conv=",l_conv,", l_heat=",l_heat
  !---- Computations for the Nusselt number if we are anelastic
  !     Can be done before setting the fields
  IF (l_heat) THEN

     IF ( index(interior_model,'EARTH') /= 0 ) THEN
        topcond=-1.D0/epsS*dtemp0(1)
        botcond=-1.D0/epsS*dtemp0(n_r_max)
     ELSE
        CALL s_cond(s0)
        CALL get_dr(s0,ds0,1,1,1,n_r_max,n_cheb_max, &
               &    w1,w2,i_costf_init,d_costf_init,drx)
        topcond=-1.D0/DSQRT(4.D0*pi)*ds0(1)
        botcond=-1.D0/DSQRT(4.D0*pi)*ds0(n_r_max)
     END IF
  END IF


  !-- Start with setting fields to zero:
  !   Touching the fields with the appropriate processor
  !   for the LM-distribute parallel region (LMLoop) makes
  !   sure that they are located close the individual
  !   processors in memory:

  IF (rank.EQ.0) THEN
     IF ( l_start_file ) THEN
        !PERFON('readFlds')
        CALL readStartFields( w,dwdtLast,z,dzdtLast, &
             &                p,dpdtLast,s,dsdtLast, &
             &                b,dbdtLast,aj,djdtLast, &
             &                b_ic,dbdt_icLast,aj_ic,djdt_icLast, &
             &                omega_ic,omega_ma, &
             &                lorentz_torque_icLast,lorentz_torque_maLast, &
             &                time,dt,dtNew,n_time_step)
        IF ( dt > 0.D0 ) THEN
           WRITE(message,'(''! Using old time step:'',D16.6)') dt
        ELSE
           dt=dtMax
           WRITE(message,'(''! Using dtMax time step:'',D16.6)') dtMax
        END IF
        !PERFOFF
     ELSE
        ! Initialize with zero
        IF ( l_conv ) THEN
           w       =zero
           dwdtLast=zero
           z       =zero
           dzdtLast=zero
           p       =zero
           dpdtLast=zero
        END IF
        IF ( l_heat ) THEN
           s       =zero
           dsdtLast=zero
        END IF
        IF ( l_mag ) THEN
           b       =zero
           dbdtLast=zero
           aj      =zero
           djdtLast=zero
        END IF
        IF ( l_cond_ic ) THEN
           b_ic       =zero
           dbdt_icLast=zero
           aj_ic      =zero
           djdt_icLast=zero
        END IF

        time =0.D0
        dt   =dtMax
        dtNew=dtMax
        n_time_step=0
        WRITE(message,'(''! Using dtMax time step:'',D16.6)') dtMax
     END IF
     CALL logWrite(message)

     !----- Get radial derivatives and initialize:

     DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lmStartReal=2*lmStart-1
        lmStopReal =2*lmStop

        !----- Initialize/add magnetic field:
        IF ( ( imagcon /= 0 .OR. init_b1 /= 0 .OR. lGrenoble ) &
             & .AND. ( l_mag .OR. l_mag_LF ) ) THEN
           !CALL initB(b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc, &
           !     lorentz_torque_icLast, lorentz_torque_maLast, &
           !     lmStart,lmStop)
           CALL initB(b,aj,b_ic,aj_ic, &
                &     lorentz_torque_icLast, lorentz_torque_maLast, &
                &     lmStart,lmStop)
        END IF

        !----- Initialize/add velocity, set IC and ma rotation:
        IF ( l_conv .OR. l_mag_kin .OR. l_SRIC .OR. l_SRMA ) THEN
           !CALL initV(w_LMloc,z_LMloc,omega_ic,omega_ma,lmStart,lmStop)
           CALL initV(w,z,omega_ic,omega_ma,lmStart,lmStop)
        END IF

        !----- Initialize/add entropy:
        IF ( ( init_s1 /= 0 .OR. impS /= 0 ) .AND. l_heat ) THEN
           !CALL initS(s_LMloc,lmStart,lmStop)
           CALL initS(s,lmStart,lmStop)
        END IF

        IF (DEBUG_OUTPUT) THEN
           WRITE(*,"(A,I3,10ES22.15)") "direct after init: w,z,s,b,aj ",nLMB,SUM(w), &
                & SUM(z), SUM(s),SUM(b),SUM(aj)
        END IF

     END DO ! Loop over LM blocks

  END IF

  ! ========== Redistribution of the fields ============
  ! 1. Broadcast the scalars
  CALL MPI_Bcast(omega_ic,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(omega_ma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(lorentz_torque_icLast,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(lorentz_torque_maLast,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(dtNew,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(n_time_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! 2. Scatter the d?dtLast arrays, they are only used in LMLoop
  !write(*,"(4X,A)") "Start Scatter d?dtLast arrays"
  DO nR=1,n_r_max
     !write(*,"(8X,A,I4)") "nR = ",nR
     CALL scatter_from_rank0_to_lo(dwdtLast(1,nR),dwdtLast_LMloc(llm,nR))
     CALL scatter_from_rank0_to_lo(dzdtLast(1,nR),dzdtLast_lo(llm,nR))
     CALL scatter_from_rank0_to_lo(dpdtLast(1,nR),dpdtLast_LMloc(llm,nR))
     CALL scatter_from_rank0_to_lo(dsdtLast(1,nR),dsdtLast_LMloc(llm,nR))

     IF (l_mag) THEN
        CALL scatter_from_rank0_to_lo(dbdtLast(1,nR),dbdtLast_LMloc(llm,nR))
        CALL scatter_from_rank0_to_lo(djdtLast(1,nR),djdtLast_LMloc(llm,nR))
     END IF
  END DO
  IF (l_cond_ic) THEN
     DO nR=1,n_r_ic_max
        CALL scatter_from_rank0_to_lo(dbdt_icLast(1,nR),dbdt_icLast_LMloc(llm,nR))
        CALL scatter_from_rank0_to_lo(djdt_icLast(1,nR),djdt_icLast_LMloc(llm,nR))
     END DO
  END IF

  ! 3. Scatter the fields to the LMloc space
  !write(*,"(4X,A)") "Start Scatter the fields"
  IF (rank.EQ.0) WRITE(*,"(A,2ES20.12)") "init z = ",SUM(z)
  DO nR=1,n_r_max
     CALL scatter_from_rank0_to_lo(w(1,nR),w_LMloc(llm:,nR))
     CALL scatter_from_rank0_to_lo(z(1,nR),z_LMloc(llm:,nR))
     CALL scatter_from_rank0_to_lo(p(1,nR),p_LMloc(llm:,nR))
     CALL scatter_from_rank0_to_lo(s(1,nR),s_LMloc(llm:,nR))
     IF (l_mag) THEN
        CALL scatter_from_rank0_to_lo(b(1,nR),b_LMloc(llmMag:,nR))
        CALL scatter_from_rank0_to_lo(aj(1,nR),aj_LMloc(llmMag:,nR))
     END IF
     !IF (DEBUG_OUTPUT) THEN
     !   IF (rank.EQ.0) THEN
     !      WRITE(*,"(A,I4,6ES22.14)") "full arrays: ",nR,SUM( s(:,nR) ),SUM( b(:,nR) ),SUM( aj(:,nR) )
     !   END IF
     !   WRITE(*,"(A,I4,6ES22.14)") "LMloc arrays: ",nR,SUM( s_LMloc(:,nR) ),SUM( b_LMloc(:,nR) ),SUM( aj_LMloc(:,nR) )
     !END IF

  END DO
  !WRITE(*,"(A,2ES20.12)") "init z_LMloc = ",get_global_sum(z_LMloc)
  IF (l_cond_ic) THEN
     DO nR=1,n_r_ic_max
        CALL scatter_from_rank0_to_lo(b_ic(1,nR),b_ic_LMloc(llm,nR))
        CALL scatter_from_rank0_to_lo(aj_ic(1,nR),aj_ic_LMloc(llm,nR))
     END DO
  END IF

  !IF (DEBUG_OUTPUT) THEN
  !   IF (rank.EQ.0) THEN
  !      WRITE(*,"(A,4ES20.12)") "getStartFields: z,dzdtLast full = ",SUM( z ),SUM( dzdtLast )
  !   END IF!

  !   WRITE(*,"(A,4ES20.12)") "getStartFields: z,dzdtLast = ",SUM( z_LMloc ),SUM( dzdtLast_lo )
  !END IF

  ALLOCATE(workA_LMloc(llm:ulm,n_r_max))
  ALLOCATE(workB_LMloc(llm:ulm,n_r_max))

  !  print*,"Computing derivatives"
  DO nLMB=1+rank*nLMBs_per_rank,MIN((rank+1)*nLMBs_per_rank,nLMBs) ! Blocking of loop over all (l,m)
     lmStart=lmStartB(nLMB)
     lmStop =lmStopB(nLMB)
     lmStartReal=2*lmStart-1
     lmStopReal =2*lmStop

     !IF (DEBUG_OUTPUT) THEN
     !   WRITE(*,"(A,I3,10ES22.15)") "after init: w,z,s,b,aj ",nLMB,SUM(w_LMloc), SUM(z_LMloc), SUM(s_LMloc),SUM(b_LMloc),SUM(aj_LMloc)
     !END IF


     IF ( l_conv .OR. l_mag_kin ) THEN
        CALL get_ddr( w_LMloc,dw_LMloc,ddw_LMloc,ulm_real-llm_real+1, &
             lmStartReal-llm_real+1,lmStopReal-llm_real+1, &
             n_r_max,n_cheb_max,workA_LMloc,workB_LMloc, &
             i_costf_init,d_costf_init,drx,ddrx)
        CALL get_dr( z_LMloc,dz_LMloc,ulm_real-llm_real+1, &
             lmStartReal-llm_real+1,lmStopReal-llm_real+1, &
             n_r_max,n_cheb_max,workA_LMloc,workB_LMloc, &
             i_costf_init,d_costf_init,drx)
     END IF

     IF ( l_mag .OR. l_mag_kin  ) THEN
        CALL get_ddr( b_LMloc,db_LMloc,ddb_LMloc,ulm_realMag-llm_realMag+1, &
             lmStartReal-llm_realMag+1,lmStopReal-llm_realMag+1, &
             n_r_max,n_cheb_max,workA_LMloc,workB_LMloc, &
             i_costf_init,d_costf_init,drx,ddrx)
        CALL get_ddr( aj_LMloc,dj_LMloc,ddj_LMloc,ulm_realMag-llm_realMag+1, &
             lmStartReal-llm_realMag+1,lmStopReal-llm_realMag+1, &
             n_r_max,n_cheb_max,workA_LMloc,workB_LMloc, &
             i_costf_init,d_costf_init,drx,ddrx)
     END IF
     IF ( l_cond_ic ) then
        CALL get_ddr_even(b_ic_LMloc,db_ic_LMLoc,ddb_ic_LMloc,ulm_realMag-llm_realMag+1, &
             lmStartReal-llm_realMag+1,lmStopReal-llm_realMag+1, &
             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA_LMloc,workB_LMloc, &
             i_costf1_ic_init,d_costf1_ic_init, &
             i_costf2_ic_init,d_costf2_ic_init)
        CALL get_ddr_even(aj_ic_LMloc,dj_ic_LMloc,ddj_ic_LMloc,ulm_realMag-llm_realMag+1, &
             lmStartReal-llm_realMag+1,lmStopReal-llm_realMag+1, &
             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA_LMloc,workB_LMloc, &
             i_costf1_ic_init,d_costf1_ic_init, &
             i_costf2_ic_init,d_costf2_ic_init)
     END IF

     IF ( l_LCR ) THEN
       DO nR=n_r_cmb,n_r_icb-1
          IF ( nR<=n_r_LCR ) THEN
             DO lm=lmStart,lmStop
                l=lo_map%lm2l(lm)
                m=lo_map%lm2m(lm)

                b_LMloc(lm,nR)=(r(n_r_LCR)/r(nR))**DBLE(l)* &
                                                b_LMloc(lm,n_r_LCR)
                db_LMloc(lm,nR)=-DBLE(l)*(r(n_r_LCR))**DBLE(l)/ &
                         (r(nR))**DBLE(l+1)*b_LMloc(lm,n_r_LCR)
                ddb_LMloc(lm,nR)=DBLE(l)*DBLE(l+1)*(r(n_r_LCR))**(l)/ &
                         (r(nR))**DBLE(l+2)*b_LMloc(lm,n_r_LCR)
                aj_LMloc(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                dj_LMloc(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                ddj_LMloc(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
             END DO
          END IF
       END DO
     END IF


     IF ( l_heat ) THEN
        !-- Get radial derivatives of entropy:
        !IF (DEBUG_OUTPUT) THEN
         !  DO nR=1,n_r_max
        !      WRITE(*,"(A,I4)") "nR=",nR
        !      DO lm=lmStart,lmStop
        !         WRITE(*,"(4X,A,4I5,2ES22.14)") "s : ", nR,lm,lo_map%lm2l(lm),lo_map%lm2m(lm),&
        !              &s_LMloc(lm,nR)
        !      END DO
        !   END DO
        !END IF
        CALL get_dr( s_LMloc,ds_LMloc,ulm_real-llm_real+1, &
             lmStartReal-llm_real+1,lmStopReal-llm_real+1, &
             n_r_max,n_cheb_max,workA_LMloc,workB_LMloc, &
             i_costf_init,d_costf_init,drx)
     END IF

     IF (DEBUG_OUTPUT) THEN
        !DO nR=1,n_r_max
        !   WRITE(*,"(A,I5,4ES22.14)") "Rdep: s,ds : ", nR,SUM(s_LMloc(lmStart:lmStop,nR)),SUM(ds_LMloc(lmStart:lmStop,nR))
        !END DO
        !DO lm=lmStart,lmStop
        !   WRITE(*,"(A,3I5,4ES22.14)") "s,ds : ", lm,lo_map%lm2l(lm),lo_map%lm2m(lm),&
        !        &SUM(s_LMloc(lm,:)),SUM(ds_LMloc(lm,:))
        !END DO
        WRITE(*,"(A,I3,10ES22.15)") "derivatives: w,z,s,b,aj ",nLMB,SUM(dw_LMloc), SUM(dz_LMloc), &
             & SUM(ds_LMloc),SUM(db_LMloc),SUM(dj_LMloc)
     END IF
     
  END DO

  deallocate(workA_LMloc)
  deallocate(workB_LMloc)
  !--- Get symmetry properties of tops excluding l=m=0:
  sES=0.D0
  sEA=0.D0
  sAA=0.D0
  DO m=0,l_max,minc
     DO l=m,l_max
        IF ( l > 0 ) THEN
           IF ( MOD(l+m,2) == 0 ) THEN
              sES=sES+cc2real(tops(l,m),m)
           ELSE
              sEA=sEA+cc2real(tops(l,m),m)
           END IF
           IF ( m /= 0 ) sAA=sAA+cc2real(tops(l,m),m)
        END IF
     END DO
  END DO
  IF ( sEA+sES == 0 ) THEN
     WRITE(message,'(''! Only l=m=0 comp. in tops:'')')
     CALL logWrite(message)
  ELSE
     sEA=DSQRT(sEA/(sEA+sES))
     sAA=DSQRT(sAA/(sEA+sES))
     WRITE(message,'(''! Rel. RMS equ. asym. tops:'',D16.6)') sEA
     CALL logWrite(message)
     WRITE(message,'(''! Rel. RMS axi. asym. tops:'',D16.6)') sAA
     CALL logWrite(message)
  END IF

  !----- Get changes in mantle and ic rotation rate:
  IF ( .NOT. l_mag_LF ) THEN
     lorentz_torque_icLast=0.D0
     lorentz_torque_maLast=0.D0
  END IF
  IF ( l_z10mat ) THEN
     l1m0=lo_map%lm2(1,0)
     coex=-2.D0*(alpha-1.D0)
     IF ( ( .NOT. l_SRMA .AND. ktopv == 2 .AND. l_rot_ma ).AND.&
          & (l1m0.ge.llm .and.l1m0.le.ulm) ) THEN
        d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_maLast
        d_omega_ma_dtLast=d_omega_ma_dt -           &
             coex * ( 2.d0*or1(1)*REAL(z_LMloc(l1m0,1)) - &
             REAL(dz_LMloc(l1m0,1)) )
     END IF
     IF ( ( .NOT. l_SRIC .AND. kbotv == 2 .AND. l_rot_ic ).AND.&
          & (l1m0.ge.llm .and. l1m0.le.ulm) ) THEN
        d_omega_ic_dt=LFfac*c_lorentz_ic*lorentz_torque_icLast
        d_omega_ic_dtLast= d_omega_ic_dt +                      &
             coex * ( 2.D0*or1(n_r_max)*REAL(z_LMloc(l1m0,n_r_max)) - &
             REAL(dz_LMloc(l1m0,n_r_max)) )
     END IF
  ELSE
     d_omega_ma_dtLast=0.D0
     d_omega_ic_dtLast=0.D0
  END IF


     ! --------------- end of insertion ----------

  !print*,"Start redistribution in getStartfields"
  ! start the redistribution
  IF (l_heat) THEN
     CALL lo2r_redist_start(lo2r_s,s_LMloc_container,s_Rloc_container)
     !CALL lo2r_redist_start(lo2r_s,s_LMloc,s_Rloc)
     !CALL lo2r_redist_start(lo2r_ds,ds_LMloc,ds_Rloc)
  END IF
  IF (l_conv) THEN
     CALL lo2r_redist_start(lo2r_z,z_LMloc_container,z_Rloc_container)
     !CALL lo2r_redist_start(lo2r_z,z_LMloc,z_Rloc)
     !CALL lo2r_redist_start(lo2r_dz,dz_LMloc,dz_Rloc)

     !DO nR=1,n_r_max
     !   WRITE(*,"(A,I2,A,2ES20.12)") "before: dw_LMloc for nR=",nR," is ",SUM( dw_LMloc(llm:ulm,nR) )
     !END DO

     CALL lo2r_redist_start(lo2r_w,w_LMloc_container,w_Rloc_container)
     !CALL lo2r_redist_start(lo2r_w,w_LMloc,w_Rloc)
     !CALL lo2r_redist_start(lo2r_dw,dw_LMloc,dw_Rloc)
     !CALL lo2r_redist_start(lo2r_ddw,ddw_LMloc,ddw_Rloc)
     CALL lo2r_redist_start(lo2r_p,p_LMloc_container,p_Rloc_container)
     !CALL lo2r_redist_start(lo2r_p,p_LMloc,p_Rloc)
     !CALL lo2r_redist_start(lo2r_dp,dp_LMloc,dp_Rloc)
  END IF

  IF (l_mag) THEN
     CALL lo2r_redist_start(lo2r_b,  b_LMloc_container,b_Rloc_container)
     !CALL lo2r_redist_start(lo2r_b,  b_LMloc,b_Rloc)
     !CALL lo2r_redist_start(lo2r_db, db_LMloc,db_Rloc)
     !CALL lo2r_redist_start(lo2r_ddb,ddb_LMloc,ddb_Rloc)
     
     CALL lo2r_redist_start(lo2r_aj, aj_LMloc_container,aj_Rloc_container)
     !CALL lo2r_redist_start(lo2r_aj, aj_LMloc,aj_Rloc)
     !CALL lo2r_redist_start(lo2r_dj, dj_LMloc,dj_Rloc)
  END IF

  !WRITE(*,"(A,10ES22.15)") "end of getStartFields: w,z,s,b,aj ",GET_GLOBAL_SUM(w_LMloc), &
  !     & GET_GLOBAL_SUM(z_LMloc), GET_GLOBAL_SUM(s_LMloc),GET_GLOBAL_SUM(b_LMloc),GET_GLOBAL_SUM(aj_LMloc)

  !print*,"End of getStartFields"
  !PERFOFF
  RETURN
end SUBROUTINE getStartFields


!---------------------------------------------------------------------------
!-- end of subroutine getStartFields
!---------------------------------------------------------------------------
