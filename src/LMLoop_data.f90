module LMLoop_data

   use parallel_mod, only: rank, nLMBs_per_rank, n_procs
   use blocking, only: nLMBs, sizeLMB, lmStartB, lmStopB
   use logic, only: l_mag
#ifndef WITH_MPI
   use truncation, only: lm_max, lm_maxMag 
#endif
 
   implicit none
 
   private
 
   integer, public :: llm,ulm,llmMag,ulmMag
 
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
#else
      llm = 1
      ulm = lm_max
      llmMag = 1
      ulmMag = lm_maxMag
#endif

      if ( DEBUG_OUTPUT ) then
         write(*,"(4(A,I6))") "llm = ",llm,", ulm = ",ulm,", llmMag = ", &
                              llmMag,", ulmMag = ",ulmMag
      end if

   end subroutine initialize_LMLoop_data
!---------------------------------------------------------------------------
end module LMLoop_data

