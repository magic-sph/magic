!$Id$
module fieldsLast

   use precision_mod, only: cp
   use truncation, only: n_r_max, lm_max, n_r_maxMag, lm_maxMag, &
                         n_r_ic_maxMag
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use parallel_Mod, only: rank


   implicit none

   private

   !--- The following variables labeled Last are provided
   !    by the restart file for the first time step or
   !    calculated here or by the update routines for the
   !    following time step.
   !    These fields remain in the LM-distributed space 
 
   complex(cp), public, allocatable :: dwdtLast(:,:)
   complex(cp), public, allocatable :: dpdtLast(:,:)
   complex(cp), public, allocatable :: dwdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: dpdtLast_LMloc(:,:)
 
   complex(cp), public, allocatable :: dzdtLast(:,:)
   complex(cp), public, allocatable :: dzdtLast_lo(:,:)
 
   complex(cp), public, allocatable :: dsdtLast(:,:)
   complex(cp), public, allocatable :: dsdtLast_LMloc(:,:)
 
   complex(cp), public, allocatable :: dbdtLast(:,:)
   complex(cp), public, allocatable :: djdtLast(:,:)
   complex(cp), public, allocatable :: dbdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: djdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: dbdt_icLast(:,:)
   complex(cp), public, allocatable :: djdt_icLast(:,:)
   complex(cp), public, allocatable :: dbdt_icLast_LMloc(:,:)
   complex(cp), public, allocatable :: djdt_icLast_LMloc(:,:)
 
   real(cp), public :: d_omega_ma_dtLast,d_omega_ic_dtLast
   real(cp), public :: lorentz_torque_maLast,lorentz_torque_icLast

   public :: initialize_fieldsLast

contains

  subroutine initialize_fieldsLast

      if ( rank == 0 ) then
         allocate( dwdtLast(lm_max,n_r_max) )
         allocate( dpdtLast(lm_max,n_r_max) )
         allocate( dzdtLast(lm_max,n_r_max) )
         allocate( dsdtLast(lm_max,n_r_max) )
         allocate( dbdtLast(lm_maxMag,n_r_maxMag) )
         allocate( djdtLast(lm_maxMag,n_r_maxMag) )
         allocate( dbdt_icLast(lm_maxMag,n_r_ic_maxMag) )
         allocate( djdt_icLast(lm_maxMag,n_r_ic_maxMag) )
      else
         allocate( dwdtLast(1,n_r_max) )
         allocate( dpdtLast(1,n_r_max) )
         allocate( dzdtLast(1,n_r_max) )
         allocate( dsdtLast(1,n_r_max) )
         allocate( dbdtLast(1,n_r_max) )
         allocate( djdtLast(1,n_r_max) )
         allocate( dbdt_icLast(1,n_r_max) )
         allocate( djdt_icLast(1,n_r_max) )
      end if
      allocate( dwdtLast_LMloc(llm:ulm,n_r_max) )
      allocate( dpdtLast_LMloc(llm:ulm,n_r_max) )
      allocate( dzdtLast_lo(llm:ulm,n_r_max) )
      allocate( dsdtLast_LMloc(llm:ulm,n_r_max) )

      allocate( dbdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
      allocate( djdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
      allocate( dbdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( djdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    
   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
