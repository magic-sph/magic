!$Id$
!*************************************************************************
SUBROUTINE dtBrms(time)
  !*************************************************************************

  !------------ This is release 2 level 1  --------------!
  !------------ Created on 1/17/02  by JW. --------------!

  !-------------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE RMS,ONLY: dtBPolLMr,dtBPol2hInt,dtBPolAs2hInt,&
       & dtBTor2hInt,dtBTorAs2hInt
  USE dtB_mod,ONLY: PstrLM,TstrLM,PadvLM,TadvLM,TomeLM,& ! st_map
       & PdifLM,TdifLM,& ! these are in lo order
       &PstrRms,TstrRms,PstrAsRms,TstrAsRms,PadvRms,TadvRms,&
       &PadvAsRms,TadvAsRms,&
       &PdifRms,TdifRms,PdifAsRms,TdifAsRms,&
       &TomeRms,TomeAsRms,dtB_gather_Rloc_on_rank0
  USE output_data
  USE const
  USE parallel_mod
  USE integration, ONLY: rInt_R
  USE communications, ONLY:myallgather
  USE LMLoop_data,ONLY: llm,ulm
  IMPLICIT NONE

  !-- Input of variables:
  REAL(kind=8),intent(IN) :: time

  !-- Local
  INTEGER :: nR,n,l1m0,lm
  CHARACTER(len=80) :: fileName

  REAL(kind=8) :: dtBPolRms,dtBPolAsRms
  REAL(kind=8) :: dtBTorRms,dtBTorAsRms
  REAL(kind=8) :: DstrRms
  REAL(kind=8) :: DadvRms
  REAL(kind=8) :: DdifRms
  REAL(kind=8) :: DdynRms
  REAL(kind=8) :: PdynRms,PdynAsRms
  REAL(kind=8) :: TdynRms,TdynAsRms
  REAL(kind=8) :: dummy1,dummy2,dummy3

  !----- Note: five full additional fields needed here!
  !      This may cause memory problems.
  COMPLEX(kind=8) :: PdynLM(lm_max_dtB,n_r_max_dtB)
  COMPLEX(kind=8) :: drPdynLM(lm_max_dtB,n_r_max_dtB)
  COMPLEX(kind=8) :: TdynLM(lm_max_dtB,n_r_max_dtB)
  COMPLEX(kind=8) :: workA(lm_max_dtB,n_r_max_dtB)
  COMPLEX(kind=8) :: workB(lm_max_dtB,n_r_max_dtB)
  !COMPLEX(kind=8) :: PdifLM_global(lm_max_dtB,n_r_max_dtB)
  !COMPLEX(kind=8) :: TdifLM_global(lm_max_dtB,n_r_max_dtB)

  !-- For new movie output
  INTEGER :: nField,nFields,nFieldSize
  INTEGER :: nTheta,nThetaN,nThetaS,nThetaStart
  INTEGER :: nPos
  REAL(kind=8) :: dumm(12),rS
  REAL(kind=8) :: fOut(n_theta_max*n_r_max)
  REAL(kind=8) :: outBlock(nfs)
  CHARACTER(len=80) :: version
  LOGICAL :: lRmsMov

  REAL(kind=8) :: global_sum(n_r_max)
  !INTEGER,DIMENSION(0:n_procs-1) :: recvcounts,displs
  !INTEGER :: sendcount
  !-- end of declaration
  !----------------------------------------------------------------------
  !CALL dtB_gather_Rloc_on_rank0

  CALL myAllGather(dtBPolLMr,lm_maxMag,n_r_maxMag)
  !CALL myAllGather(dtBTorLMr,lm_maxMag,n_r_maxMag)
  !CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
  !     & dtBPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  !CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
  !     & dtBTorLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

  ! dtBPol2hInt and dtBPolAs2hInt must be added over all ranks.
  CALL mpi_reduce(dtBPol2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  if (rank.eq.0) dtBPol2hInt(:,1)=global_sum
  CALL mpi_reduce(dtBPolAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  if (rank.eq.0) dtBPolAs2hInt(:,1)=global_sum

  !DO nR=1,n_r_max
  !   WRITE(*,"(I3,ES20.12)") nR,dtBTor2hInt(nR,1)
  !END DO
  CALL mpi_reduce(dtBTor2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) THEN
     dtBTor2hInt(:,1)=global_sum
     !   DO nR=1,n_r_max
     !      WRITE(*,"(A,I3,ES20.12)") "global_sum = ",nR,global_sum(nR)
     !   END DO
  END IF


  ! for dtBPolLMr we need a gathering in LM direction as it is set in updateB for
  ! all nR
  !CALL MPI_AllGatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
  !     & dtBPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)


  IF (rank.EQ.0) THEN
     !WRITE(*,"(A,ES20.13)") "dtBPol2hInt = ",SUM( dtBPol2hInt(:,1) )
     l1m0=lm2(1,0)

     !--- Stretching
     CALL get_drNS(PstrLM,workA,lm_max_real,1,lm_max_real, &
          &        n_r_max,n_cheb_max,workB, &
          &        i_costf_init,d_costf_init,drx)
     !--- Finalize rms poloidal and toroidal stretching:
     !WRITE(*,"(A,6ES22.12)") "PstrLM,drPstrLM,TstrLM = ",SUM(PstrLM),SUM(workA),SUM(TstrLM)
     CALL get_PolTorRms(PstrLM,workA,TstrLM,PstrRms,TstrRms,PstrAsRms,TstrAsRms,st_map)
     !WRITE(*,"(A,4ES22.12)") "PstrRms,TstrRMS = ",PstrRms,TstrRms,PstrAsRms,TstrAsRms

     !--- Calculate dipole stretching and copy for total dynamo term:
     DO nR=1,n_r_max
        DO lm=1,lm_max
           PdynLM(lm,nR)  =PstrLM(lm,nR)
           drPdynLM(lm,nR)=workA(lm,nR)
           IF ( lm /= l1m0 ) THEN
              PstrLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              workA(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
        END DO
     END DO
     !--- Get dipole stretching
     CALL get_PolTorRms(PstrLM,workA,TstrLM,DstrRms,dummy1,dummy2,dummy3,st_map)

     !--- Finalize advection
     CALL get_drNS(PadvLM,workA,lm_max_real,1,lm_max_real, &
          &        n_r_max,n_cheb_max,workB, &
          &        i_costf_init,d_costf_init,drx)
     CALL get_PolTorRms(PadvLM,workA,TadvLM, &
          &             PadvRms,TadvRms,PadvAsRms,TadvAsRms,st_map)
     DO nR=1,n_r_max
        DO lm=1,lm_max
           PstrLM(lm,nR)  =PdynLM(lm,nR)
           PdynLM(lm,nR)  =PdynLM(lm,nR)-PadvLM(lm,nR)
           drPdynLM(lm,nR)=drPdynLM(lm,nR)+workA(lm,nR)
           TdynLM(lm,nR)  =TstrLM(lm,nR)-TadvLM(lm,nR)
           IF ( lm /= l1m0 ) THEN
              PadvLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              workA(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
        END DO
     END DO
     !--- Get dipole advection and total dynamo terms:
     CALL get_PolTorRms(PadvLM,workA,TadvLM, DadvRms,dummy1,dummy2,dummy3,st_map)
     CALL get_PolTorRms(PdynLM,drPdynLM,TdynLM, PdynRms,TdynRms,PdynAsRms,TdynAsRms,st_map)
     DO nR=1,n_r_max
        DO lm=1,lm_max
           PadvLM(lm,nR)=PdynLM(lm,nR)
           IF ( lm /= l1m0 ) THEN
              PdynLM(lm,nR)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              drPdynLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
        END DO
     END DO
     !--- Get dipole dynamo terms:
     CALL get_PolTorRms(PdynLM,drPdynLM,TdynLM, DdynRms,dummy1,dummy2,dummy3,st_map)

     !--- Diffusion:
     !WRITE(*,"(A,2ES22.15)") "dtBrms, PdifLM = ",SUM( PdifLM )
     CALL get_drNS(PdifLM,workA,lm_max_real,1,lm_max_real, &
          &        n_r_max,n_cheb_max,workB, &
          &        i_costf_init,d_costf_init,drx)
     CALL get_PolTorRms(PdifLM,workA,TdifLM,PdifRms,TdifRms,PdifAsRms,TdifAsRms,st_map)

     DO nR=1,n_r_max
        DO lm=1,lm_max
           PdynLM(lm,nR)=PdifLM(lm,nR)
           IF ( lm /= st_map%lm2(1,0) ) THEN
              PdifLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              workA(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
        END DO
     END DO
     CALL get_PolTorRms(PdifLM,workA,TdifLM,DdifRms,dummy1,dummy2,dummy3,st_map)


     !--- Omega effect: (different mapping for PdifLM,workA and TomeLM)
     CALL get_PolTorRms(PdifLM,workA,TomeLM,dummy1,TomeRms,dummy2,TomeAsRms,st_map)


     !--- B changes:
     CALL get_drNS(dtBPolLMr,workA,lm_max_real,1,lm_max_real, &
          &        n_r_max,n_cheb_max,workB, &
          &        i_costf_init,d_costf_init,drx)
     DO nR=1,n_r_max
        CALL hInt2dPol(workA(1,nR),2,lm_max, &
             dtBPol2hInt(nR,1),dtBPolAs2hInt(nR,1),lo_map)
        !DO n=2,nThreadsLMmax
        !   dtBPol2hInt(nR,1)  = dtBPol2hInt(nR,1)   + dtBPol2hInt(nR,n)
        !   dtBPolAs2hInt(nR,1)= dtBPolAs2hInt(nR,1) + dtBPolAs2hInt(nR,n)
        !   dtBTor2hInt(nR,1)  = dtBTor2hInt(nR,1)   + dtBTor2hInt(nR,n)
        !   dtBTorAs2hInt(nR,1)= dtBTorAs2hInt(nR,1) + dtBTorAs2hInt(nR,n)
        !END DO
     END DO
     dtBPolRms  =rInt_R(dtBPol2hInt(1,1),n_r_max, &
          &             n_r_max,drx,i_costf_init,d_costf_init)
     dtBPolAsRms=rInt_R(dtBPolAs2hInt(1,1),n_r_max, &
          &             n_r_max,drx,i_costf_init,d_costf_init)
     dtBTorRms  =rInt_R(dtBTor2hInt(1,1),n_r_max, &
          &             n_r_max,drx,i_costf_init,d_costf_init)
     dtBTorAsRms=rInt_R(dtBTorAs2hInt(1,1),n_r_max, &
          &             n_r_max,drx,i_costf_init,d_costf_init)
     dtBPolRms  =SQRT(dtBPolRms  /vol_oc)
     dtBPolAsRms=SQRT(dtBPolAsRms/vol_oc)
     dtBTorRms  =SQRT(dtBTorRms  /vol_oc)
     dtBTorAsRms=SQRT(dtBTorAsRms/vol_oc)

     !-- OUTPUT:
     IF ( l_save_out) THEN
        OPEN(n_dtbrms_file,FILE=dtbrms_file,FORM='FORMATTED', &
             STATUS='UNKNOWN',POSITION='APPEND')
     END IF
     WRITE(n_dtbrms_file,'(1P,D20.12,12D16.8)') &
          time, &
          dtBPolRms,dtBTorRms, & !2,3
          PstrRms,TstrRms, &     !4,5
          PadvRms,TadvRms, &     !6,7
          PdifRms,TdifRms, &     !8,9
          TomeRms/TstrRms,TomeRms, & !10,11
          PdynRms,TdynRms  !12,13
     IF ( l_save_out) THEN
        CLOSE(n_dtbrms_file)
     END IF

     IF ( l_save_out) THEN
        OPEN(n_dtdrms_file,FILE=dtdrms_file,FORM='FORMATTED', &
             STATUS='UNKNOWN',POSITION='APPEND')
     END IF
     WRITE(n_dtdrms_file,'(1P,D20.12,3D16.8)') &
          time,DstrRms,DadvRms,DdifRms
     IF ( l_save_out) THEN
        CLOSE(n_dtdrms_file)
     END IF

     !-- Output of movie files for axisymmetric toroidal field changes:
     !   Tstr,Tome,Tdyn=Tstr+Tadv,
     lRmsMov=.FALSE.
     IF ( lRmsMov ) THEN

        nFieldSize=n_theta_max*n_r_max
        nFields=7
        fileName='dtTas_mov.'//TAG
        OPEN(90,file=fileName,STATUS='UNKNOWN',FORM='UNFORMATTED')

        !------ Write header
        version='JW_Movie_Version_2'
        WRITE(90) version
        dumm(1)=112           ! type of input
        dumm(2)=3             ! marker for constant phi plane
        dumm(3)=0.D0          ! surface constant
        dumm(4)=nFields       ! no of fields
        WRITE(90) (real(dumm(n),4),n=1,4)

        !------ Define marker for output fields stored in movie field
        dumm(1)=101           ! Field marker for AS Br stretching
        dumm(2)=102           ! Field marker for AS Br dynamo term
        dumm(3)=103           ! Field marker for AS Br diffusion
        dumm(4)=104           ! Field marker for AS Bp stretching
        dumm(5)=105           ! Field marker for AS Bp dynamo term
        dumm(6)=106           ! Field marker for AS Bp omega effect
        dumm(7)=107           ! Field marker for AS Bp diffusion
        WRITE(90) (SNGL(dumm(n)),n=1,nFields)

        !------ Now other info about grid and parameters:
        WRITE(90) runid        ! run identifier
        dumm( 1)=n_r_max          ! total number of radial points
        dumm( 2)=n_r_max          ! no of radial point in outer core
        dumm( 3)=n_theta_max      ! no. of theta points
        dumm( 4)=n_phi_max        ! no. of phi points
        dumm( 5)=minc             ! imposed symmetry
        dumm( 6)=ra               ! control parameters
        dumm( 7)=ek               ! (for information only)
        dumm( 8)=pr               !      -"-
        dumm( 9)=prmag            !      -"-
        dumm(10)=radratio         ! ratio of inner / outer core
        dumm(11)=tScale           ! timescale
        WRITE(90) (SNGL(dumm(n)),     n=1,11)
        WRITE(90) (SNGL(r(n)/r_CMB),  n=1,n_r_max)
        WRITE(90) (SNGL(theta_ord(n)),n=1,n_theta_max)
        WRITE(90) (SNGL(phi(n)),      n=1,n_phi_max)

        dumm(1)=1    ! time frame number for movie
        dumm(2)=0.D0 ! time
        dumm(3)=0.D0
        dumm(4)=0.D0
        dumm(5)=0.D0
        dumm(6)=0.D0
        dumm(7)=0.D0
        dumm(8)=0.D0
        WRITE(90) (SNGL(dumm(n)),n=1,8)

        !------ Loop over different output field:
        Do nField=1,nFields

           !------ Loop over r and theta:
           DO nR=1,n_r_max ! Loop over radial points
              rS=r(nR)
              DO n=1,nThetaBs ! Loop over theta blocks
                 nThetaStart=(n-1)*sizeThetaB+1

                 !------ Convert from lm to theta block and store in outBlock:
                 IF ( nField == 1 ) THEN
                    CALL get_RAS(PstrLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 ELSE IF ( nField == 2 ) THEN
                    ! Note that PadvLM stores PdynLM=PstrLM+PadvLM at this point!
                    CALL get_RAS(PadvLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 ELSE IF ( nField == 3 ) THEN
                    ! Note that PdynLM stores PdifLM at this point!
                    CALL get_RAS(PdynLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 ELSE IF ( nField == 4 ) THEN
                    CALL get_PASLM(TstrLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 ELSE IF ( nField == 5 ) THEN
                    CALL get_PASLM(TdynLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 ELSE IF ( nField == 6 ) THEN
                    CALL get_PASLM(TomeLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 ELSE IF ( nField == 7 ) THEN
                    CALL get_PASLM(TdifLM(1,nR),outBlock, &
                         rS,nThetaStart,sizeThetaB)
                 END IF

                 !------ Storage of field in fout for theta block
                 DO nTheta=1,sizeThetaB,2
                    !--------- Convert to correct order in theta grid points and store of fOut:
                    nThetaN=(nThetaStart+nTheta)/2
                    nPos=(nR-1)*n_theta_max+nThetaN
                    fOut(nPos)=outBlock(nTheta)
                    nThetaS=n_theta_max-nThetaN+1
                    nPos=(nR-1)*n_theta_max+nThetaS
                    fOut(nPos)=outBlock(nTheta+1)
                 END DO ! Loop over thetas in block

              END DO ! Loop over theta blocks

           END DO ! Loop over R

           !------ Output of field:
           WRITE(90) (SNGL(fOut(nPos)),nPos=1,nFieldSize)

        END DO ! Loop over different fields

     END IF ! output of mov fields ?
  END IF

  RETURN
end SUBROUTINE dtBrms

!----------------------------------------------------------------------------
