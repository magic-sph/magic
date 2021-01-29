module fieldsLast
   !
   ! This module contains all the work arrays of the previous time-steps
   ! needed to time advance the code.
   ! They are needed in the time-stepping scheme.
   !

   use precision_mod
   use truncation, only: n_r_max, lm_max, n_r_maxMag, lm_maxMag, &
       &                 n_r_ic_maxMag, fd_order, fd_order_bound
   use blocking, only: llm, ulm, llmMag, ulmMag
   use logic, only: l_chemical_conv, l_heat, l_mag, l_cond_ic, l_double_curl, &
       &            l_RMS, l_finite_diff, l_parallel_solve, l_mag_par_solve
   use constants, only: zero
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag
   use mem_alloc, only: bytes_allocated
   use time_array

   implicit none

   private

   type(type_tarray), public :: dsdt, dwdt, dpdt, dzdt, dxidt
   type(type_tarray), public :: dbdt, djdt, dbdt_ic, djdt_ic
   type(type_tscalar), public :: domega_ma_dt, domega_ic_dt
   type(type_tscalar), public :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

   !DIR$ ATTRIBUTES ALIGN:64 :: dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dsdt_Rloc,dVSrLM_Rloc,dVXirLM_Rloc
   complex(cp), public, allocatable, target  :: dflowdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_Rloc_container(:,:,:)
   complex(cp), public, pointer :: dwdt_Rloc(:,:),dzdt_Rloc(:,:)
   complex(cp), public, pointer :: dpdt_Rloc(:,:), dsdt_Rloc(:,:), dVSrLM_Rloc(:,:)
   complex(cp), public, pointer :: dxidt_Rloc(:,:), dVXirLM_Rloc(:,:)
   complex(cp), public, pointer :: dVxVhLM_Rloc(:,:)

   !DIR$ ATTRIBUTES ALIGN:64 :: djdt_Rloc,dbdt_Rloc,dVxBhLM_Rloc
   complex(cp), public, pointer :: djdt_Rloc(:,:), dVxBhLM_Rloc(:,:)
   complex(cp), public, pointer :: dbdt_Rloc(:,:)

   ! The same arrays, but now the LM local part
   complex(cp), public, allocatable, target  :: dflowdt_LMloc_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_LMloc_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_LMloc_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_LMloc_container(:,:,:,:)
   complex(cp), public, pointer :: dVSrLM_LMloc(:,:,:), dVXirLM_LMloc(:,:,:)
   complex(cp), public, pointer :: dVxVhLM_LMloc(:,:,:), dVxBhLM_LMloc(:,:,:)

   complex(cp), public, allocatable :: dbdt_CMB_LMloc(:)

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(nold, nexp, nimp)
      !
      ! This routine defines all the arrays needed to time advance MagIC
      !

      integer, intent(in) :: nold ! Number of storage of the old state
      integer, intent(in) :: nexp ! Number of explicit states
      integer, intent(in) :: nimp ! Number of implicit states

      !-- Local variable
      integer :: n_fields

      call domega_ma_dt%initialize(nold, nexp, nimp)
      call domega_ic_dt%initialize(nold, nexp, nimp)

      call lorentz_torque_ic_dt%initialize(nold, nexp, nimp)
      call lorentz_torque_ma_dt%initialize(nold, nexp, nimp)

      if ( l_parallel_solve ) then
         if ( l_heat ) call dsdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, &
                            &               nimp, l_allocate_exp=.true.)
         call dzdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, nimp, &
              &               l_allocate_exp=.true.)
         call dwdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, nimp, &
              &               l_allocate_exp=.true.)
         if ( (.not. l_double_curl) .or. l_RMS ) then
            call dpdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, nimp, &
                 &               l_allocate_exp=.true.)
         else
            allocate( dpdt%expl(1,1,nexp) ) ! For debug
         end if
         if ( l_chemical_conv ) call dxidt%initialize(1, lm_max, nRstart,nRstop, nold, &
                                     &                nexp, nimp, l_allocate_exp=.true.)
         if ( l_mag .and. l_mag_par_solve ) then
            call dbdt%initialize(1, lm_maxMag, nRstartMag, nRstopMag, nold, nexp, nimp, &
                 &               l_allocate_exp=.true.)
            call djdt%initialize(1, lm_maxMag, nRstartMag, nRstopMag, nold, nexp, nimp, &
                 &               l_allocate_exp=.true.)
         else
            call dbdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
            call djdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
         end if
      else
         if ( l_heat ) call dsdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
         call dzdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
         call dwdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
         if ( (.not. l_double_curl) .or. l_RMS ) then
            call dpdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
         end if
         if ( l_mag ) then
            call dbdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
            call djdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
         end if
         if ( l_chemical_conv ) call dxidt%initialize(llm, ulm, 1, n_r_max, nold, &
                                     &                nexp, nimp)
      end if

      if ( l_cond_ic ) then
         call dbdt_ic%initialize(llmMag, ulmMag, 1, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
         call djdt_ic%initialize(llmMag, ulmMag, 1, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
      end if

      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         if ( l_parallel_solve ) then
            if ( l_mag .and. (.not. l_mag_par_solve) ) then
               allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
               dbdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
               djdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            else
               allocate( dbdt_Rloc(1,1), djdt_Rloc(1,1) )
            end if
         else
            n_fields=3
            if ( l_mag ) n_fields=n_fields+2
            allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:n_fields) )
            dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            dsdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
            if ( l_mag .and. (.not. l_mag_par_solve) ) then
               dbdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,4)
               djdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,5)
            end if
            allocate(dpdt_Rloc(lm_max,nRstart:nRstop))
         end if
         allocate(dVxVhLM_Rloc(lm_max,nRstart:nRstop))
         allocate(dVSrLM_Rloc(lm_max,nRstart:nRstop))
         allocate(dVxBhLM_Rloc(lm_maxMag,nRstartMag:nRstopMag))
         bytes_allocated = bytes_allocated+                               &
         &                 6*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX+& 
         &                 3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
      else
         if ( l_double_curl ) then
            allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:4) )
            dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            dpdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
            dVxVhLM_Rloc(1:,nRstart:) => &
            &                         dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,4)
            bytes_allocated = bytes_allocated+ &
            &                 4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         else
            allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:3) )
            dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            dpdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
            allocate( dVxVhLM_Rloc(1:1,1:1) )
            bytes_allocated = bytes_allocated+ &
            &                 3*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         end if

         allocate( dsdt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         dsdt_Rloc(1:,nRstart:)   => dsdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dVSrLM_Rloc(1:,nRstart:) => dsdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated+ &
         &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

         ! the magnetic part
         allocate( dbdt_Rloc_container(lm_maxMag,nRstartMag:nRstopMag,1:3) )
         dbdt_Rloc(1:,nRstartMag:) => &
         &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,1)
         djdt_Rloc(1:,nRstartMag:) => &
         &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,2)
         dVxBhLM_Rloc(1:,nRstartMag:)=> &
         &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,3)
         bytes_allocated = bytes_allocated+ &
         &                 3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         if ( l_parallel_solve ) then
            allocate( dVXirLM_Rloc(lm_max,nRstart:nRstop) )
            bytes_allocated = bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         else
            allocate( dxidt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
            dxidt_Rloc(1:,nRstart:)   => dxidt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dVXirLM_Rloc(1:,nRstart:) => dxidt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            bytes_allocated = bytes_allocated+ &
            &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         end if
      else
         allocate( dxidt_Rloc_container(1,1,1:2) )
         dxidt_Rloc(1:,1:)   => dxidt_Rloc_container(1:1,1:1,1)
         dVXirLM_Rloc(1:,1:) => dxidt_Rloc_container(1:1,1:1,2)
      end if

      !-- Set the initial values to zero
      if ( l_mag ) then
         if ( .not. l_mag_par_solve ) then
            dbdt_Rloc(:,:)=zero
            djdt_Rloc(:,:)=zero
         end if
         dVxBhLM_Rloc(:,:)=zero
      end if
      if ( .not. l_parallel_solve ) then
         dwdt_Rloc(:,:)=zero
         dzdt_Rloc(:,:)=zero
         dsdt_Rloc(:,:)=zero
         dpdt_Rloc(:,:)=zero
      end if
      dVSrLM_Rloc(:,:)=zero
      if ( l_double_curl ) dVxVhLM_Rloc(:,:)=zero
      if ( l_chemical_conv ) then
         if (.not. l_parallel_solve ) dxidt_Rloc(:,:)  =zero
         dVXirLM_Rloc(:,:)=zero
      end if

      ! The same arrays, but now the LM local part
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         if ( l_parallel_solve ) then
            if ( l_mag .and. (.not. l_mag_par_solve) ) then
               allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
               dbdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
               djdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
               bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
               &                 SIZEOF_DEF_COMPLEX
            end if
         else
            n_fields=3
            if ( l_mag ) n_fields=n_fields+2
            !--@> TODO: clean this ugly stuff:
            allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:n_fields,1:nexp))
            dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            dsdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
            bytes_allocated = bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
            if ( l_mag ) then
               dbdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,4,1:nexp)
               djdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,5,1:nexp)
               bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
               &                 SIZEOF_DEF_COMPLEX
            end if
            if ( ((.not. l_double_curl) .or. l_RMS) ) then
               allocate( dpdt%expl(llm:ulm,n_r_max,nexp) )
               bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*nexp* &
               &                 SIZEOF_DEF_COMPLEX
            else
               allocate( dpdt%expl(1,1,1) ) ! To avoid debug
            end if
         end if
      else ! This is either high-order F.D. or Cheb
         if ( l_double_curl ) then
            allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:4,1:nexp))
            dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            dpdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
            dVxVhLM_LMloc(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,4,1:nexp)
            bytes_allocated = bytes_allocated+4*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:3,1:nexp))
            dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            dpdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
            allocate( dVxVhLM_LMloc(1:1,1:1,1:1) )
            bytes_allocated = bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
         end if

         allocate(dsdt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
         if ( .not. l_parallel_solve ) then
            dsdt%expl(llm:,1:,1:) => dsdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
         end if
         dVSrLM_LMloc(llm:,1:,1:) => dsdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
         bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX

         allocate(dbdt_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3,1:nexp))
         dbdt%expl(llmMag:,1:,1:) => dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,1,1:nexp)
         djdt%expl(llmMag:,1:,1:) => dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,2,1:nexp)
         dVxBhLM_LMloc(llmMag:,1:,1:) => &
         &                         dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,3,1:nexp)
         bytes_allocated = bytes_allocated+ &
         &                 3*nexp*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         if ( .not. l_parallel_solve ) then
            allocate(dxidt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
            dxidt%expl(llm:,1:,1:)   => dxidt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dVXirLM_LMloc(llm:,1:,1:) => dxidt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate(dxidt_LMloc_container(1,1,1:2,1))
            !dxidt%expl(1:,1:,1:)   => dxidt_LMloc_container(1:1,1:1,1,1:)
            dVXirLM_LMloc(1:,1:,1:) => dxidt_LMloc_container(1:1,1:1,2,1:)
         end if
      else
         allocate(dxidt_LMloc_container(1,1,1:2,1))
         dxidt%expl(1:,1:,1:)   => dxidt_LMloc_container(1:1,1:1,1,1:)
         dVXirLM_LMloc(1:,1:,1:) => dxidt_LMloc_container(1:1,1:1,2,1:)
      end if

      ! Only when l_dt_cmb_field is requested
      ! There might be a way to allocate only when needed
      allocate ( dbdt_CMB_LMloc(llmMag:ulmMag) )
      bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast
      !
      ! Memory deallocation of d?dt arrays.
      !

      if ( (.not. l_parallel_solve) .and. (.not. l_mag_par_solve) ) then
         deallocate( dflowdt_Rloc_container, dflowdt_LMloc_container )
      end if
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         deallocate( dVxVhLM_Rloc, dVxBhLM_Rloc, dVSrLM_Rloc)
         if (.not. l_parallel_solve ) deallocate( dpdt_Rloc )
      else
         deallocate( dbdt_Rloc_container, dbdt_LMloc_container )
         deallocate( dsdt_Rloc_container, dsdt_LMloc_container )
         if ( .not. l_double_curl ) deallocate( dVxVhLM_Rloc, dVxVhLM_LMloc )
      end if
      deallocate( dbdt_CMB_LMloc )

      if ( l_chemical_conv ) then
         if ( .not. l_parallel_solve ) then
            deallocate( dxidt_Rloc_container, dxidt_LMloc_container )
         else
            deallocate( dVXirLM_Rloc )
         end if
      end if

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
