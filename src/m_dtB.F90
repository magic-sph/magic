!$Id$
!---------------------------------------------------------------------------------
!  This module contains magnetic field stretching and advection terms 
!  plus a separate omega-effect.
!  It is used for movie output.
!--------------------------------------------------------------------------------

MODULE dtB_mod
  use truncation
  USE radial_data,ONLY: nRstart,nRstop
  USE parallel_mod,ONLY: nr_per_rank,nr_on_last_rank,MPI_COMM_WORLD,n_procs,&
       &MPI_DOUBLE_COMPLEX,rank
  USE LMLoop_data,ONLY: llmMag,ulmMag
  implicit none
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: PstrLM,TstrLM,PadvLM,TadvLM,TomeLM
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: PstrLM_Rloc,TstrLM_Rloc,&
       &PadvLM_Rloc,TadvLM_Rloc,TomeLM_Rloc

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: PdifLM,TdifLM,&
       &PadvLMIC,PdifLMIC,TadvLMIC,TdifLMIC
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: PdifLM_LMloc,TdifLM_LMloc,&
       &PadvLMIC_LMloc,PdifLMIC_LMloc,TadvLMIC_LMloc,TdifLMIC_LMloc

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: TstrRLM,TadvRLM,TomeRLM
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: TstrRLM_Rloc,TadvRLM_Rloc,TomeRLM_Rloc

  REAL(kind=8) :: PstrRms,PstrAsRms
  REAL(kind=8) :: PadvRms,PadvAsRms
  REAL(kind=8) :: PdifRms,PdifAsRms
  REAL(kind=8) :: TstrRms,TstrAsRms
  REAL(kind=8) :: TadvRms,TadvAsRms
  REAL(kind=8) :: TdifRms,TdifAsRms
  REAL(kind=8) :: TomeRms,TomeAsRms

CONTAINS
  SUBROUTINE initialize_dtB_mod

    IF (rank.EQ.0) THEN
       ALLOCATE( PstrLM(lm_max_dtB,n_r_max_dtB) )
       ALLOCATE( PadvLM(lm_max_dtB,n_r_max_dtB) )
       ALLOCATE( TstrLM(lm_max_dtB,n_r_max_dtB) )
       ALLOCATE( TadvLM(lm_max_dtB,n_r_max_dtB) )
       ALLOCATE( TomeLM(lm_max_dtB,n_r_max_dtB) )
    ELSE
       ALLOCATE( PstrLM(1,1) )
       ALLOCATE( PadvLM(1,1) )
       ALLOCATE( TstrLM(1,1) )
       ALLOCATE( TadvLM(1,1) )
       ALLOCATE( TomeLM(1,1) )
    END IF
    ALLOCATE( PstrLM_Rloc(lm_max_dtB,nRstart:nRstop) )
    ALLOCATE( PadvLM_Rloc(lm_max_dtB,nRstart:nRstop) )
    ALLOCATE( TstrLM_Rloc(lm_max_dtB,nRstart:nRstop) )
    ALLOCATE( TadvLM_Rloc(lm_max_dtB,nRstart:nRstop) )
    ALLOCATE( TomeLM_Rloc(lm_max_dtB,nRstart:nRstop) )

    IF (rank.EQ.0) THEN
       ALLOCATE( PdifLM(lm_max_dtB,n_r_max_dtB) )
       ALLOCATE( TdifLM(lm_max_dtB,n_r_max_dtB) )
       ALLOCATE( PadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
       ALLOCATE( PdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
       ALLOCATE( TadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
       ALLOCATE( TdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
    ELSE
       ALLOCATE( PdifLM(1,1) )
       ALLOCATE( TdifLM(1,1) )
       ALLOCATE( PadvLMIC(1,1) )
       ALLOCATE( PdifLMIC(1,1) )
       ALLOCATE( TadvLMIC(1,1) )
       ALLOCATE( TdifLMIC(1,1) )
    END IF
    ALLOCATE( PdifLM_LMloc(llmMag:ulmMag,n_r_max_dtB) )
    ALLOCATE( TdifLM_LMloc(llmMag:ulmMag,n_r_max_dtB) )
    ALLOCATE( PadvLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
    ALLOCATE( PdifLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
    ALLOCATE( TadvLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
    ALLOCATE( TdifLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )

    IF (rank.EQ.0) THEN
       ALLOCATE(TstrRLM(lm_max_dtB,n_r_max_dtB))
       ALLOCATE(TadvRLM(lm_max_dtB,n_r_max_dtB))
       ALLOCATE(TomeRLM(lm_max_dtB,n_r_max_dtB))
    ELSE
       ALLOCATE(TstrRLM(1,1))
       ALLOCATE(TadvRLM(1,1))
       ALLOCATE(TomeRLM(1,1))
    END IF
    ALLOCATE(TstrRLM_Rloc(lm_max_dtB,nRstart:nRstop))
    ALLOCATE(TadvRLM_Rloc(lm_max_dtB,nRstart:nRstop))
    ALLOCATE(TomeRLM_Rloc(lm_max_dtB,nRstart:nRstop))

  END SUBROUTINE initialize_dtB_mod

SUBROUTINE dtb_gather_Rloc_on_rank0
  INTEGER :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
  INTEGER :: i,ierr

  IF (ldtBmem.EQ.1) THEN
     sendcount  = (nRstop-nRstart+1)*lm_max_dtB
     recvcounts = nr_per_rank*lm_max_dtB
     recvcounts(n_procs-1) = nr_on_last_rank*lm_max_dtB
     DO i=0,n_procs-1
        displs(i) = i*nr_per_rank*lm_max_dtB
     END DO
     CALL MPI_GatherV(TstrRLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & TstrRLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(TadvRLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & TadvRLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(TomeRLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & TomeRLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

     CALL MPI_GatherV(TstrLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & TstrLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(TadvLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & TadvLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(PstrLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & PstrLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(PadvLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & PadvLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

     CALL MPI_GatherV(TomeLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
          & TomeLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
  END IF
END SUBROUTINE dtb_gather_Rloc_on_rank0

END MODULE dtB_mod
