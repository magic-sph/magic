!$Id$
module fieldsLast

   implicit none

   !--- The following variables labeled Last are provided
   !    by the restart file for the first time step or
   !    calculated here or by the update routines for the
   !    following time step.
   !    These fields remain in the LM-distributed space 
 
   complex(kind=8), allocatable :: dwdtLast(:,:)
   complex(kind=8), allocatable :: dpdtLast(:,:)
   complex(kind=8), allocatable :: dwdtLast_LMloc(:,:)
   complex(kind=8), allocatable :: dpdtLast_LMloc(:,:)
 
   complex(kind=8), allocatable :: dzdtLast(:,:)
   complex(kind=8), allocatable :: dzdtLast_lo(:,:)
 
   complex(kind=8), allocatable :: dsdtLast(:,:)
   complex(kind=8), allocatable :: dsdtLast_LMloc(:,:)
 
   complex(kind=8), allocatable :: dbdtLast(:,:)
   complex(kind=8), allocatable :: djdtLast(:,:)
   complex(kind=8), allocatable :: dbdtLast_LMloc(:,:)
   complex(kind=8), allocatable :: djdtLast_LMloc(:,:)
   complex(kind=8), allocatable :: dbdt_icLast(:,:)
   complex(kind=8), allocatable :: djdt_icLast(:,:)
   complex(kind=8), allocatable :: dbdt_icLast_LMloc(:,:)
   complex(kind=8), allocatable :: djdt_icLast_LMloc(:,:)
 
   real(kind=8) :: d_omega_ma_dtLast,d_omega_ic_dtLast
   real(kind=8) :: lorentz_torque_maLast,lorentz_torque_icLast

contains

  subroutine initialize_fieldsLast

      use truncation, only: n_r_max, lm_max, n_r_maxMag, lm_maxMag, &
                            n_r_ic_maxMag
      use LMLoop_data, only: llm,ulm,llmMag,ulmMag
      use parallel_Mod, only: rank

      if (rank == 0) then
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

end module fieldsLast
