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
   use blocking, only: llm, ulm, llmMag, ulmMag
   use logic, only: l_chemical_conv, l_heat, l_mag, l_cond_ic
   use parallel_Mod, only: rank
   use time_array

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

   type(type_tarray), public :: dsdt, dwdt, dpdt, dzdt, dxidt
   type(type_tarray), public :: dbdt, djdt, dbdt_ic, djdt_ic
   type(type_tscalar), public :: domega_ma_dt
   type(type_tscalar), public :: domega_ic_dt

   real(cp), public :: d_omega_ma_dtLast,d_omega_ic_dtLast
   real(cp), public :: lorentz_torque_maLast,lorentz_torque_icLast

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(norder_imp, norder_exp,        &
              &                     norder_imp_lin)

      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin


      allocate( dwdtLast_LMloc(llm:ulm,n_r_max) )
      allocate( dpdtLast_LMloc(llm:ulm,n_r_max) )
      allocate( dzdtLast_lo(llm:ulm,n_r_max) )
      allocate( dsdtLast_LMloc(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated + &
      &                 4*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( dxidtLast_LMloc(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      else
         allocate( dxidtLast_LMloc(1,1) )
      end if

      allocate( dbdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
      allocate( djdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
      bytes_allocated = bytes_allocated + &
      &                 2*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      allocate( dbdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( djdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      bytes_allocated = bytes_allocated + &
      &                 2*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      call domega_ma_dt%initialize(norder_imp, norder_exp, norder_imp_lin)
      call domega_ic_dt%initialize(norder_imp, norder_exp, norder_imp_lin)

      if ( l_heat ) call dsdt%initialize(llm, ulm, n_r_max, norder_imp, norder_exp, &
                         &               norder_imp_lin)
      if ( l_chemical_conv ) call dxidt%initialize(llm, ulm, n_r_max, norder_imp, &
                                  &                norder_exp, norder_imp_lin)

      if ( l_mag ) then
         call dbdt%initialize(llmMag, ulmMag, n_r_maxMag, norder_imp, norder_exp, &
              &               norder_imp_lin)
         call djdt%initialize(llmMag, ulmMag, n_r_maxMag, norder_imp, norder_exp, &
              &               norder_imp_lin)
      end if

      if ( l_cond_ic ) then
         call dbdt_ic%initialize(llmMag, ulmMag, n_r_ic_maxMag, norder_imp, &
              &                  norder_exp, norder_imp_lin, l_allocate_exp=.true.)
         call djdt_ic%initialize(llmMag, ulmMag, n_r_ic_maxMag, norder_imp, &
              &                  norder_exp, norder_imp_lin, l_allocate_exp=.true.)
      end if

      call dwdt%initialize(llm, ulm, n_r_max, norder_imp, norder_exp, &
           &               norder_imp_lin)
      call dpdt%initialize(llm, ulm, n_r_max, norder_imp, norder_exp, &
           &               norder_imp_lin)
      call dzdt%initialize(llm, ulm, n_r_max, norder_imp, norder_exp, &
           &               norder_imp_lin)

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      call domega_ma_dt%finalize()
      call domega_ic_dt%finalize()
      call dzdt%finalize()
      call dpdt%finalize()
      call dwdt%finalize()
      if ( l_heat ) call dsdt%finalize()
      if ( l_chemical_conv ) call dxidt%finalize()
      if ( l_mag ) then
         call dbdt%finalize()
         call djdt%finalize()
      end if
      if ( l_cond_ic ) then
         call dbdt_ic%finalize()
         call djdt_ic%finalize()
      end if
      deallocate( dwdtLast_LMloc, dpdtLast_LMloc, dzdtLast_lo )
      deallocate( dsdtLast_LMloc, dbdtLast_LMloc, djdtLast_LMloc )
      deallocate( dbdt_icLast_LMloc, djdt_icLast_LMloc )
      deallocate( dxidtLast_LMloc )

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
