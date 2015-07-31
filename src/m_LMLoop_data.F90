!$Id$
module LMLoop_data

   use parallel_mod, only: rank, nLMBs_per_rank, n_procs
   use blocking, only: nLMBs, sizeLMB, lmStartB, lmStopB
   use logic, only: l_mag
 
   implicit none
 
   private
 
   integer, public :: llm,ulm,llmMag,ulmMag
   integer, public :: llm_real,ulm_real,llm_realMag,ulm_realMag
   integer, public :: lm_per_rank,lm_on_last_rank
 
   public :: initialize_LMLoop_data

contains

   subroutine initialize_LMLoop_data
    
      !-- Local variables:
      integer :: nLMB_start,nLMB_end
      logical :: DEBUG_OUTPUT=.false.

      ! set the local lower and upper index for lm
#ifdef WITH_MPI
      ! we have nLMBs LM blocks which are distributed over the ranks
      ! with nLMBs_per_rank blocks per rank (computed in m_blocking.F90)
      
      nLMB_start = 1+rank*nLMBs_per_rank
      nLMB_end   = min((rank+1)*nLMBs_per_rank,nLMBs)
      llm = lmStartB(nLMB_start)
      ulm = lmStopB(nLMB_end)
      if ( l_mag ) then
         llmMag = llm
         ulmMag = ulm
      else
         llmMag = 1
         ulmMag = 1
      end if
      lm_per_rank=nLMBs_per_rank*sizeLMB
      lm_on_last_rank=lmStopB (min(n_procs*nLMBs_per_rank,nLMBs))- &
                      lmStartB(1+(n_procs-1)*nLMBs_per_rank)+1
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
      if ( DEBUG_OUTPUT ) then
         write(*,"(4(A,I6))") "llm = ",llm,", ulm = ",ulm,", llmMag = ", &
                              llmMag,", ulmMag = ",ulmMag
      end if

   end subroutine initialize_LMLoop_data
!---------------------------------------------------------------------------
end module LMLoop_data

