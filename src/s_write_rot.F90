!$Id$
!***********************************************************************
SUBROUTINE write_rot(time,dt,eKinIC,ekinMA,w,z,dz,b, &
     &               omega_ic,omega_ma, &
     &               lorentz_torque_ic,lorentz_torque_ma)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to write information on the        |
  !  |  outputfile file_rot.                                             |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  use mpi
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking,ONLY: lo_map,st_map,lmStartB,lmStopB
  USE horizontal_data
  USE logic
  USE output_data
  USE const
  USE parallel_mod,only: rank
  USE LMLoop_data,ONLY: llm,ulm,llmMag,ulmMag
  USE integration, ONLY: rInt

  IMPLICIT NONE

  !-- Input of variables:
  REAL(kind=8),INTENT(IN) :: omega_ic,omega_ma
  REAL(kind=8),INTENT(IN) :: lorentz_torque_ma,lorentz_torque_ic
  REAL(kind=8),INTENT(IN) :: time,dt

  !-- Input of scalar fields of toroidal flow and poloidal mag. field:
  COMPLEX(kind=8),INTENT(IN) :: w(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: z(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: dz(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: b(llmMag:ulmMag,n_r_maxMag)

  !-- Output into rot_file
  REAL(kind=8),INTENT(OUT) :: eKinIC,eKinMA

  !-- Local variables:
  REAL(kind=8),PARAMETER :: tolerance=1e-16
  REAL(kind=8) :: eKinOC
  INTEGER :: n_r1,n_r2,n_r3,nR
  INTEGER :: l1m0,l1m1
  REAL(kind=8) :: viscous_torque_ic,viscous_torque_ma
  REAL(kind=8) :: AMz,eKinAMz
  REAL(kind=8) :: angular_moment_oc(3)
  REAL(kind=8) :: angular_moment_ic(3)
  REAL(kind=8) :: angular_moment_ma(3)
  COMPLEX(kind=8) :: z10(n_r_max),z11(n_r_max)
  CHARACTER(len=80) :: filename

  REAL(kind=8) :: powerLor,powerVis

  REAL(kind=8) :: AMzLast,eKinAMzLast
  SAVE AMzLast,eKinAMzLast

  INTEGER, DIMENSION(:,:),POINTER :: lm2
  INTEGER :: i,l,m,ilm,lm_vals(21),n_lm_vals
  COMPLEX(kind=8),DIMENSION(8,3) :: zvals_on_rank0,bvals_on_rank0
  COMPLEX(kind=8),DIMENSION(21) :: vals_on_rank0_1d

  INTEGER :: sr_tag,status(MPI_STATUS_SIZE),ierr
  LOGICAL :: rank_has_l1m0,rank_has_l1m1
  LOGICAL :: DEBUG_OUTPUT=.false.
  !-- end of declaration
  !-----------------------------------------------------------------------
  ! some arbitrary tag for the send and recv
  sr_tag=12345

  lm2(0:,0:) => lo_map%lm2
  l1m0=lm2(1,0)

  IF (DEBUG_OUTPUT) WRITE(*,"(I3,A,3I6)") rank,":lmStartB,lmStopB,l1m0=",lmStartB(rank+1),lmStopB(rank+1),l1m0

  IF ( lmStartB(rank+1).LE.l1m0 .AND. lmStopB(rank+1).GE.l1m0 ) THEN
     !IF (rank.NE.0) THEN
     !   PRINT*,"in s_write_rot, l1m0 not on rank 0"
     !   stop
     !END IF
     !-- Calculating viscous torques:
     IF ( l_rot_ic .AND. kbotv == 2 ) THEN
        CALL get_viscous_torque(viscous_torque_ic, &
             &                  z(l1m0,n_r_max),dz(l1m0,n_r_max),r_icb)
     ELSE
        viscous_torque_ic=0.d0
     END IF
     IF ( l_rot_ma .AND. ktopv == 2 ) THEN
        CALL get_viscous_torque(viscous_torque_ma, &
             &                  z(l1m0,1),dz(l1m0,1),r_cmb)
     ELSE
        viscous_torque_ma=0.d0
     END IF
     rank_has_l1m0=.TRUE.
     IF (rank.NE.0) THEN
        ! send viscous_torque_ic and viscous_torque_ma to rank 0 for 
        ! output
        CALL MPI_Send(viscous_torque_ic,1,MPI_DOUBLE_PRECISION,0,&
             &sr_tag,MPI_COMM_WORLD,ierr)
        CALL MPI_Send(viscous_torque_ma,1,MPI_DOUBLE_PRECISION,0,&
             &sr_tag+1,MPI_COMM_WORLD,ierr)
     END IF
  ELSE
     rank_has_l1m0=.FALSE.
  END IF

  IF (rank.EQ.0) THEN
     IF (.NOT.rank_has_l1m0) THEN
        CALL MPI_Recv(viscous_torque_ic,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
             &sr_tag,MPI_COMM_WORLD,status,ierr)
        CALL MPI_Recv(viscous_torque_ma,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
             &sr_tag+1,MPI_COMM_WORLD,status,ierr)
     END IF
     IF ( l_SRIC ) THEN
        powerLor=lorentz_torque_ic*omega_IC
        powerVis=viscous_torque_ic*omega_IC
        OPEN(n_SRIC_file,file=SRIC_file,status="unknown", &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1p,2x,d20.12,4d17.6)') &
             time*tScale, &
             omega_ic/tScale, &
             (powerLor+powerVis)*vScale*vScale/tScale, &
             powerVis*vScale*vScale/tScale, &
             powerLor*vScale*vScale/tScale
        CLOSE(n_SRIC_file)
     END IF
     IF ( l_SRMA ) THEN
        powerLor=lorentz_torque_ma*omega_ma
        powerVis=viscous_torque_ma*omega_ma
        OPEN(n_SRMA_file,file=SRMA_file,status="unknown", &
             POSITION='APPEND')
        WRITE(n_SRMA_file,'(1p,2x,d20.12,4d17.6)') &
             time*tScale, &
             omega_ma/tScale, &
             (powerLor+powerVis)*vScale*vScale/tScale, &
             powerVis*vScale*vScale/tScale, &
             powerLor*vScale*vScale/tScale
        CLOSE(n_SRMA_file)
     END IF
  END IF

  IF ( l_drift ) THEN
     DO i=1,4
        lm_vals(i)=lm2(i*minc,i*minc)
        lm_vals(4+i)=lm2(i*minc+1,i*minc)
     END DO
     n_r1=INT(1.D0/3.D0*(n_r_max-1))
     n_r2=INT(2.D0/3.D0*(n_r_max-1))
     n_r3=n_r_max-1
     CALL sendvals_to_rank0(z,n_r1,lm_vals(1:8),zvals_on_rank0(:,1))
     CALL sendvals_to_rank0(z,n_r2,lm_vals(1:8),zvals_on_rank0(:,2))
     CALL sendvals_to_rank0(z,n_r3,lm_vals(1:8),zvals_on_rank0(:,3))

     IF (rank.EQ.0) THEN
        filename='driftVD.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,24D12.4)') &
             time, &
             (zvals_on_rank0(ilm,1),ilm=1,4), &
             (zvals_on_rank0(ilm,2),ilm=1,4), &
             (zvals_on_rank0(ilm,3),ilm=1,4)
        CLOSE(n_SRIC_file)
        filename='driftVQ.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,24D12.4)') &
             time, &
             (zvals_on_rank0(ilm,1),ilm=5,8), &
             (zvals_on_rank0(ilm,2),ilm=5,8), &
             (zvals_on_rank0(ilm,3),ilm=5,8)
        CLOSE(n_SRIC_file)
     END IF
     
     IF ( l_mag .OR. l_mag_LF ) THEN
        n_r1=n_r_CMB
        n_r2=n_r_ICB
        CALL sendvals_to_rank0(b,n_r1,lm_vals(1:8),bvals_on_rank0(:,1))
        CALL sendvals_to_rank0(b,n_r2,lm_vals(1:8),bvals_on_rank0(:,2))

        IF (rank.EQ.0) THEN
           filename='driftBD.'//tag
           OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
                POSITION='APPEND')
           WRITE(n_SRIC_file,'(1P,2X,D20.12,16D12.4)') &
                time, &
                (bvals_on_rank0(ilm,1),ilm=5,8), &
                (bvals_on_rank0(ilm,2),ilm=5,8)
           CLOSE(n_SRIC_file)
           filename='driftBQ.'//tag
           OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
                POSITION='APPEND')
           WRITE(n_SRIC_file,'(1P,2X,D20.12,16D12.4)') &
                time, &
                (bvals_on_rank0(ilm,1),ilm=1,4), &
                (bvals_on_rank0(ilm,2),ilm=1,4)
           CLOSE(n_SRIC_file)
        END IF
     END IF ! l_mag
  END IF

  IF ( .NOT. l_SRIC .AND. ( l_rot_ic .OR. l_rot_ma ) ) THEN
     IF (rank.EQ.0) THEN
        IF ( l_save_out ) THEN
           OPEN(n_rot_file,file=rot_file,status='UNKNOWN', &
                POSITION='APPEND')
        END IF
        WRITE(n_rot_file,'(1P,2X,D20.12,6D14.6)') &
             time*tScale, &
             omega_ic/tScale, &
             lScale**2*vScale*lorentz_torque_ic, &
             lScale**2*vScale*viscous_torque_ic, &
             omega_ma/tScale, &
             lScale**2*vScale*lorentz_torque_ma, &
             -lScale**2*vScale*viscous_torque_ma
        IF ( l_save_out ) CLOSE(n_rot_file)
     END IF
  END IF
  
  IF ( l_AM ) THEN
     rank_has_l1m0=.false.
     rank_has_l1m1=.false.
     l1m0=lo_map%lm2(1,0)
     l1m1=lo_map%lm2(1,1)
     IF ( (lmStartB(rank+1).LE.l1m0) .AND. (l1m0.LE.lmStopB(rank+1)) ) THEN
        DO nR=1,n_r_max
           z10(nR)=z(l1m0,nR)
        END DO
        rank_has_l1m0=.true.
        IF (rank.NE.0) THEN
           CALL MPI_Send(z10,n_r_max,MPI_DOUBLE_COMPLEX,0,sr_tag,MPI_COMM_WORLD,ierr)
        END IF
     END IF

     IF ( l1m1 > 0 ) THEN
        IF ( (lmStartB(rank+1).LE.l1m1) .AND. (l1m1.LE.lmStopB(rank+1)) ) THEN
           DO nR=1,n_r_max
              z11(nR)=z(l1m1,nR)
           END DO
           rank_has_l1m1=.TRUE.
           IF (rank.NE.0) THEN
              CALL MPI_Send(z11,n_r_max,MPI_DOUBLE_COMPLEX,0,sr_tag+1,MPI_COMM_WORLD,ierr)
           END IF
        END IF
     ELSE
        DO nR=1,n_r_max
           z11(nR)=CMPLX(0d0,0d0,KIND=KIND(0d0))
        END DO
     END IF
     ! now we have z10 and z11 in the worst case on two different
     ! ranks, which are also different from rank 0
     IF (rank.EQ.0) THEN
        IF (.NOT.rank_has_l1m0) THEN
           CALL MPI_Recv(z10,n_r_max,MPI_DOUBLE_COMPLEX,&
                & MPI_ANY_SOURCE,sr_tag,MPI_COMM_WORLD,status,ierr)
        END IF
        IF ( l1m1 > 0 ) THEN
           IF (.NOT.rank_has_l1m1) THEN
              CALL MPI_Recv(z11,n_r_max,MPI_DOUBLE_COMPLEX,&
                   & MPI_ANY_SOURCE,sr_tag+1,MPI_COMM_WORLD,status,ierr)
           END IF
        END IF

        CALL get_angular_moment(z10,z11,omega_ic,omega_ma, &
             angular_moment_oc, &
             angular_moment_ic,angular_moment_ma)
        IF ( l_save_out ) THEN
           OPEN(n_angular_file,FILE=angular_file,STATUS='UNKNOWN', &
                POSITION='APPEND')
        END IF
        AMz=   angular_moment_oc(3) + &
             & angular_moment_ic(3)+angular_moment_ma(3)
        IF (ABS(AMz).LT.tolerance) AMz=0.0D0
        eKinAMz=0.5d0*(angular_moment_oc(3)**2/c_moi_oc + &
             angular_moment_ic(3)**2/c_moi_ic + &
             angular_moment_ma(3)**2/c_moi_ma )
        if (abs(eKinAMz).lt.tolerance) eKinAMz=0.0D0
        eKinIC=0.5d0*angular_moment_ic(3)**2/c_moi_ic
        eKinOC=0.5d0*angular_moment_oc(3)**2/c_moi_oc
        eKinMA=0.5d0*angular_moment_ma(3)**2/c_moi_ma
        IF ( AMzLast .ne. 0.0D0 ) THEN
           !WRITE(*,"(A,4ES22.15)") "col9 = ",eKinAMz,eKinAMzLast,dt,(eKinAMz-eKinAMzLast)
           WRITE(n_angular_file,'(1p,2x,d20.12,5d14.6,3d20.12)',advance='no') &
                time*tScale, &
                & angular_moment_oc, &
                & angular_moment_ic(3), &
                & angular_moment_ma(3), &
                & AMz,&
                & (AMz-AMzLast)/AMzLast/dt,&
                & eKinAMz
           IF (eKinAMzLast.NE.0.0d0) THEN
              WRITE(n_angular_file,'(1d20.12)',advance='no') &
                & (eKinAMz-eKinAMzLast)/eKinAMzLast/dt
           ELSE
              WRITE(n_angular_file,'(1d20.12)',advance='no') 0.0
           END IF
           WRITE(n_angular_file,'(3d20.12)') &
                & eKinIC,&
                & eKinOC,&
                & eKinMA
        END IF
        IF ( l_save_out ) CLOSE(n_angular_file)
        AMzLast=AMz
        eKinAMzLast=eKinAMz
     END IF
  END IF

  IF ( l_iner ) THEN
     ! l_iner can only be .TRUE. for minc=1
     n_lm_vals=0
     DO l=1,6
        DO m=1,l
           n_lm_vals = n_lm_vals + 1
           lm_vals(n_lm_vals)=lm2(l,m)
        END DO
     END DO
     n_r1=INT(1.D0/2.D0*(n_r_max-1))
     CALL sendvals_to_rank0(w,n_r1,lm_vals(1:n_lm_vals),vals_on_rank0_1d)

     IF (rank.EQ.0) THEN
        filename='inerP.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,21D12.4)') &
             time, ( REAL(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals )
        CLOSE(n_SRIC_file)
     END IF

     n_r1=INT(1.D0/2.D0*(n_r_max-1))
     CALL sendvals_to_rank0(z,n_r1,lm_vals(1:n_lm_vals),vals_on_rank0_1d)

     IF (rank.EQ.0) THEN
        filename='inerT.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,21D12.4)') &
             time, ( REAL(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals ) 
        CLOSE(n_SRIC_file)
     END IF

  END IF

  RETURN
CONTAINS
  SUBROUTINE sendvals_to_rank0(field,n_r,lm_vals,vals_on_rank0)
    COMPLEX(kind=8),DIMENSION(llm:ulm,n_r_max),INTENT(IN) :: field
    INTEGER,INTENT(IN) :: n_r
    INTEGER,DIMENSION(:),INTENT(IN) :: lm_vals
    COMPLEX(kind=8),DIMENSION(:),INTENT(OUT) :: vals_on_rank0

    INTEGER :: ilm,lm,ierr,status(MPI_STATUS_SIZE),tag,n_lm_vals
    
    n_lm_vals=SIZE(lm_vals)
    IF (SIZE(vals_on_rank0).LT.n_lm_vals) THEN
       WRITE(*,"(2(A,I4))") "write_rot: length of vals_on_rank0=",SIZE(vals_on_rank0),&
            &" must be >= size(lm_vals)=",n_lm_vals
       CALL mpi_abort(MPI_COMM_WORLD,43,ierr)
    END IF

    DO ilm=1,n_lm_vals
       lm=lm_vals(ilm)
       IF (lmStartB(1).LE.lm .AND. lm.LE.lmStopB(1)) THEN
          ! the value is already on rank 0
          if (rank.eq.0) vals_on_rank0(ilm)=field(lm,n_r)
       ELSE
          tag=876+ilm
          ! on which process is the lm value?
          IF (lmStartB(rank+1).LE.lm .AND. lm.LE.lmStopB(rank+1)) THEN
             CALL MPI_Send(field(lm,n_r),1,MPI_DOUBLE_COMPLEX,&
                  & 0,tag,MPI_COMM_WORLD,ierr)
          END IF
          IF (rank.EQ.0) THEN
             CALL MPI_Recv(vals_on_rank0(ilm),1,MPI_DOUBLE_COMPLEX,&
                  & MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,status,ierr)
          END IF
       END IF
    END DO
  END SUBROUTINE sendvals_to_rank0
end SUBROUTINE write_rot
!-----------------------------------------------------------------------
