module LMLoop_data

   use parallel_mod, only: rank, n_procs
   use blocking, only: lm_balance
   use logic, only: l_mag
   use truncation, only: lm_max, lm_maxMag 
 
   implicit none
 
   private
 
   integer, public :: llm,ulm,llmMag,ulmMag
 
   public :: initialize_LMLoop_data

contains

   subroutine initialize_LMLoop_data
    
      llm = lm_balance(rank)%nStart
      ulm = lm_balance(rank)%nStop
      if ( l_mag ) then
         llmMag = llm
         ulmMag = ulm
      else
         llmMag = 1
         ulmMag = 1
      end if

   end subroutine initialize_LMLoop_data
!---------------------------------------------------------------------------
end module LMLoop_data

