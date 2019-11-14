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
   use truncation, only: n_r_max, lm_max, n_r_maxMag, lm_maxMag, &
       &                 n_r_ic_maxMag
   use blocking, only: llm, ulm, llmMag, ulmMag
   use logic, only: l_chemical_conv, l_heat, l_mag, l_cond_ic, l_double_curl, &
       &            l_RMS
   use time_array

   implicit none

   private

   type(type_tarray), public :: dsdt, dwdt, dpdt, dzdt, dxidt
   type(type_tarray), public :: dbdt, djdt, dbdt_ic, djdt_ic
   type(type_tscalar), public :: domega_ma_dt, domega_ic_dt
   type(type_tscalar), public :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(nold, nexp, nimp)

      integer, intent(in) :: nold
      integer, intent(in) :: nexp
      integer, intent(in) :: nimp

      call domega_ma_dt%initialize(nold, nexp, nimp)
      call domega_ic_dt%initialize(nold, nexp, nimp)

      call lorentz_torque_ic_dt%initialize(nold, nexp, nimp)
      call lorentz_torque_ma_dt%initialize(nold, nexp, nimp)

      if ( l_heat ) call dsdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)
      if ( l_chemical_conv ) call dxidt%initialize(llm, ulm, n_r_max, nold, &
                                  &                nexp, nimp)

      if ( l_mag ) then
         call dbdt%initialize(llmMag, ulmMag, n_r_maxMag, nold, nexp, nimp)
         call djdt%initialize(llmMag, ulmMag, n_r_maxMag, nold, nexp, nimp)
      end if

      if ( l_cond_ic ) then
         call dbdt_ic%initialize(llmMag, ulmMag, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
         call djdt_ic%initialize(llmMag, ulmMag, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
      end if

      call dwdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)
      if ( (.not. l_double_curl) .or. l_RMS ) then
         call dpdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)
      end if
      call dzdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      call lorentz_torque_ma_dt%finalize()
      call lorentz_torque_ic_dt%finalize()
      call domega_ma_dt%finalize()
      call domega_ic_dt%finalize()
      call dzdt%finalize()
      if ( .not. l_double_curl .or. l_RMS ) call dpdt%finalize()
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

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
