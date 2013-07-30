MODULE LMLoop_data
  use parallel_mod
  use blocking
  implicit none

  private
  INTEGER,PUBLIC :: llm,ulm,llmMag,ulmMag
  INTEGER, PUBLIC :: llm_real,ulm_real,llm_realMag,ulm_realMag
  INTEGER, PUBLIC :: lm_per_rank,lm_on_last_rank

  PUBLIC :: initialize_LMLoop_data
CONTAINS
  SUBROUTINE initialize_LMLoop_data
    
    INTEGER :: nLMB_start,nLMB_end
    logical :: DEBUG_OUTPUT=.false.

    ! set the local lower and upper index for lm
#ifdef WITH_MPI
    ! we have nLMBs LM blocks which are distributed over the ranks
    ! with nLMBs_per_rank blocks per rank (computed in m_blocking.F90)
    
    nLMB_start = 1+rank*nLMBs_per_rank
    nLMB_end   = MIN((rank+1)*nLMBs_per_rank,nLMBs)
    llm = lmStartB(nLMB_start)
    ulm = lmStopB(nLMB_end)
    IF (l_mag) THEN
       llmMag = llm
       ulmMag = ulm
    ELSE
       llmMag = 1
       ulmMag = 1
    END IF
    lm_per_rank=nLMBs_per_rank*sizeLMB
    lm_on_last_rank=lmStopB(MIN(n_procs*nLMBs_per_rank,nLMBs))-lmStartB(1+(n_procs-1)*nLMBs_per_rank)+1
#else
    llm = 1
    ulm = lm_max
    llmMag = 1
    ulmMag = lm_maxMag
    lm_per_rank=lm_max
    lm_on_last_rank=lm_max
#endif

    llm_real = 2*llm-1
    ulm_real = 2*ulm
    llm_realMag = 2*llmMag-1
    ulm_realMag = 2*ulmMag
    IF (DEBUG_OUTPUT) THEN
       WRITE(*,"(4(A,I6))") "ulm = ",ulm,", llm = ",llm,", ulmMag = ",ulmMag,", llmMag = ",llmMag
    END IF
  END SUBROUTINE initialize_LMLoop_data
END MODULE LMLoop_data

