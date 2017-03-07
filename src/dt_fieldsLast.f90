module fieldsLast
   !
   ! This module contains time-derivaties array of the previous time-step
   ! They are needed in the time-stepping scheme.
   !
   ! The variables labeled with a suffix 'Last' are provided
   ! by the restart file for the first time step or
   ! calculated here or by the update routines for the
   ! following time step.
   ! These fields remain in the LM-distributed space 
 
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, n_r_maxMag, lm_maxMag, &
       &                 n_r_ic_maxMag
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use logic, only: l_chemical_conv
   use parallel_Mod, only: rank

   implicit none

   private

   complex(cp), public, allocatable :: dwdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: dpdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: dzdtLast_lo(:,:)
   complex(cp), public, allocatable :: dsdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: dxidtLast_LMloc(:,:)
 
   complex(cp), public, allocatable :: dbdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: djdtLast_LMloc(:,:)
   complex(cp), public, allocatable :: dbdt_icLast_LMloc(:,:)
   complex(cp), public, allocatable :: djdt_icLast_LMloc(:,:)

   real(cp), public :: d_omega_ma_dtLast,d_omega_ic_dtLast
   real(cp), public :: lorentz_torque_maLast,lorentz_torque_icLast

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast

      allocate( dwdtLast_LMloc(llm:ulm,n_r_max) )
      allocate( dpdtLast_LMloc(llm:ulm,n_r_max) )
      allocate( dzdtLast_lo(llm:ulm,n_r_max) )
      allocate( dsdtLast_LMloc(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated + &
                        4*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( dxidtLast_LMloc(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      else
         allocate( dxidtLast_LMloc(1,1) )
      end if

      allocate( dbdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
      allocate( djdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
      bytes_allocated = bytes_allocated + &
                        2*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      allocate( dbdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( djdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      bytes_allocated = bytes_allocated + &
                        2*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      deallocate( dwdtLast_LMloc, dpdtLast_LMloc, dzdtLast_lo )
      deallocate( dsdtLast_LMloc, dbdtLast_LMloc, djdtLast_LMloc )
      deallocate( dbdt_icLast_LMloc, djdt_icLast_LMloc )
      deallocate( dxidtLast_LMloc )

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
