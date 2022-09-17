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
       &            l_RMS, l_finite_diff, l_parallel_solve, l_mag_par_solve,  &
       &            l_phase_field
   use constants, only: zero
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use time_array

   implicit none

   private

   type(type_tarray), public :: dsdt, dwdt, dpdt, dzdt, dxidt
   type(type_tarray), public :: dbdt, djdt, dbdt_ic, djdt_ic, dphidt
   type(type_tscalar), public :: domega_ma_dt, domega_ic_dt
   type(type_tscalar), public :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

   !DIR$ ATTRIBUTES ALIGN:64 :: dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dsdt_Rloc,dVSrLM_Rloc,dVXirLM_Rloc
   complex(cp), public, allocatable, target  :: dflowdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable  :: dphidt_Rloc(:,:)
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

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: domega_ma_dt, domega_ic_dt)
      !$omp target update to(domega_ma_dt, domega_ic_dt) nowait
      !$omp target enter data map(alloc: lorentz_torque_ic_dt, lorentz_torque_ma_dt)
      !$omp target update to(lorentz_torque_ic_dt, lorentz_torque_ma_dt) nowait
#endif

      if ( l_parallel_solve ) then
         if ( l_heat ) then
            call dsdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, &
                               &               nimp, l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dsdt)
            !$omp target update to(dsdt) nowait
#endif
         end if

         call dzdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, nimp, &
              &               l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dzdt)
         !$omp target update to(dzdt) nowait
#endif
         call dwdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, nimp, &
              &               l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dwdt)
         !$omp target update to(dwdt) nowait
#endif

         if ( (.not. l_double_curl) .or. l_RMS ) then
            call dpdt%initialize(1, lm_max, nRstart, nRstop, nold, nexp, nimp, &
                 &               l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dpdt)
            !$omp target update to(dpdt) nowait
#endif
         else
            allocate( dpdt%expl(1,1,nexp) ) ! For debug
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dpdt%expl)
            !$omp target update to(dpdt%expl) nowait
#endif
         end if

         if ( l_chemical_conv ) then
            call dxidt%initialize(1, lm_max, nRstart,nRstop, nold, &
                                        &                nexp, nimp, l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dxidt)
            !$omp target update to(dxidt) nowait
#endif
         end if

         if ( l_phase_field ) then
            call dphidt%initialize(1, lm_max, nRstart,nRstop, nold, &
                                      &                 nexp, nimp, l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dphidt)
            !$omp target update to(dphidt) nowait
#endif
         end if

         if ( l_mag .and. l_mag_par_solve ) then
            call dbdt%initialize(1, lm_maxMag, nRstartMag, nRstopMag, nold, nexp, nimp, &
                 &               l_allocate_exp=.true.)
            call djdt%initialize(1, lm_maxMag, nRstartMag, nRstopMag, nold, nexp, nimp, &
                 &               l_allocate_exp=.true.)
         else
            call dbdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
            call djdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
         end if
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dbdt, djdt)
         !$omp target update to(dbdt, djdt) nowait
#endif
      else
         if ( l_heat ) then
            call dsdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dsdt)
            !$omp target update to(dsdt) nowait
#endif
         end if

         call dzdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
         call dwdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dzdt, dwdt)
         !$omp target update to(dzdt, dwdt) nowait
#endif

         if ( (.not. l_double_curl) .or. l_RMS ) then
            call dpdt%initialize(llm, ulm, 1, n_r_max, nold, nexp, nimp)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dpdt)
            !$omp target update to(dpdt) nowait
#endif
         end if

         if ( l_mag ) then
            call dbdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
            call djdt%initialize(llmMag, ulmMag, 1, n_r_maxMag, nold, nexp, nimp)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dbdt, djdt)
            !$omp target update to(dbdt, djdt) nowait
#endif
         end if

         if ( l_chemical_conv ) then
            call dxidt%initialize(llm, ulm, 1, n_r_max, nold, &
                                     &                nexp, nimp)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dxidt)
            !$omp target update to(dxidt) nowait
#endif
         end if

         if ( l_phase_field ) then
            call dphidt%initialize(llm, ulm, 1, n_r_max, nold, &
                                   &                 nexp, nimp, l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dphidt)
            !$omp target update to(dphidt) nowait
#endif
         end if
      end if

      if ( l_cond_ic ) then
         call dbdt_ic%initialize(llmMag, ulmMag, 1, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
         call djdt_ic%initialize(llmMag, ulmMag, 1, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dbdt_ic, djdt_ic)
         !$omp target update to(dbdt_ic, djdt_ic) nowait
#endif
      end if

      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         if ( l_parallel_solve ) then
            if ( l_mag .and. (.not. l_mag_par_solve) ) then
               allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
               dflowdt_Rloc_container(:,:,:)=zero
               dbdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
               djdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: dbdt_Rloc, djdt_Rloc)
               !$omp target update to(dbdt_Rloc, djdt_Rloc) nowait
#endif
            else
               allocate( dbdt_Rloc(1,1), djdt_Rloc(1,1) )
               dbdt_Rloc(1,1)=zero; djdt_Rloc(1,1)=zero
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: dbdt_Rloc, djdt_Rloc)
               !$omp target update to(dbdt_Rloc, djdt_Rloc) nowait
#endif
            end if
         else
            n_fields=3
            if ( l_mag ) n_fields=n_fields+2
            allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:n_fields) )
            dflowdt_Rloc_container(:,:,:)=zero
            dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            dsdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dwdt_Rloc, dzdt_Rloc, dsdt_Rloc)
            !$omp target update to(dwdt_Rloc, dzdt_Rloc, dsdt_Rloc) nowait
#endif
            if ( l_mag .and. (.not. l_mag_par_solve) ) then
               dbdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,4)
               djdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,5)
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: dbdt_Rloc, djdt_Rloc)
               !$omp target update to(dbdt_Rloc, djdt_Rloc) nowait
#endif
            end if
            allocate(dpdt_Rloc(lm_max,nRstart:nRstop))
            dpdt_Rloc(:,:)=zero
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dpdt_Rloc)
            !$omp target update to(dpdt_Rloc) nowait
#endif
         end if
         allocate(dVxVhLM_Rloc(lm_max,nRstart:nRstop))
         allocate(dVSrLM_Rloc(lm_max,nRstart:nRstop))
         allocate(dVxBhLM_Rloc(lm_maxMag,nRstartMag:nRstopMag))
         dVxVhLM_Rloc(:,:)=zero
         dVSrLM_Rloc(:,:) =zero
         dVxBhLM_Rloc(:,:)=zero
         bytes_allocated = bytes_allocated+                               &
         &                 6*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX+& 
         &                 3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dVxVhLM_Rloc, dVSrLM_Rloc, dVxBhLM_Rloc)
         !$omp target update to(dVxVhLM_Rloc, dVSrLM_Rloc, dVxBhLM_Rloc) nowait
         gpu_bytes_allocated = gpu_bytes_allocated+                               &
         &                     6*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX+&
         &                     3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
#endif
      else
         if ( l_double_curl ) then
            allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:4) )
            dflowdt_Rloc_container(:,:,:)=zero
            dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            dpdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
            dVxVhLM_Rloc(1:,nRstart:) => &
            &                         dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,4)
            bytes_allocated = bytes_allocated+ &
            &                 4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dwdt_Rloc, dzdt_Rloc, dpdt_Rloc, dVxVhLM_Rloc)
            !$omp target update to(dwdt_Rloc, dzdt_Rloc, dpdt_Rloc, dVxVhLM_Rloc) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+ &
            &                     4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif
         else
            allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:3) )
            dflowdt_Rloc_container(:,:,:)=zero
            dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            dpdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
            !!allocate( dVxVhLM_Rloc(1:1,1:1) )
            allocate( dVxVhLM_Rloc(lm_max,nRstart:nRstop) )
            dVxVhLM_Rloc(:,:) = zero
            bytes_allocated = bytes_allocated+ &
            &                 4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dwdt_Rloc, dzdt_Rloc, dpdt_Rloc, dVxVhLM_Rloc)
            !$omp target update to(dwdt_Rloc, dzdt_Rloc, dpdt_Rloc, dVxVhLM_Rloc) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+ &
            &                     4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif
         end if

         allocate( dsdt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         dsdt_Rloc_container(:,:,:)=zero
         dsdt_Rloc(1:,nRstart:)   => dsdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dVSrLM_Rloc(1:,nRstart:) => dsdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated+ &
         &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dsdt_Rloc, dVSrLM_Rloc)
         !$omp target update to(dsdt_Rloc, dVSrLM_Rloc) nowait
         gpu_bytes_allocated = gpu_bytes_allocated+ &
         &                     2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif

         ! the magnetic part
         allocate( dbdt_Rloc_container(lm_maxMag,nRstartMag:nRstopMag,1:3) )
         dbdt_Rloc_container(:,:,:)=zero
         dbdt_Rloc(1:,nRstartMag:) => &
         &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,1)
         djdt_Rloc(1:,nRstartMag:) => &
         &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,2)
         dVxBhLM_Rloc(1:,nRstartMag:)=> &
         &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,3)
         bytes_allocated = bytes_allocated+ &
         &                 3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dbdt_Rloc, djdt_Rloc, dVxBhLM_Rloc)
         !$omp target update to(dbdt_Rloc, djdt_Rloc, dVxBhLM_Rloc) nowait
         gpu_bytes_allocated = gpu_bytes_allocated+ &
         &                     3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
#endif
      end if

      if ( l_chemical_conv ) then
         if ( l_parallel_solve ) then
            allocate( dVXirLM_Rloc(lm_max,nRstart:nRstop) )
            dVXirLM_Rloc(:,:)=zero
            bytes_allocated = bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dVXirLM_Rloc)
            !$omp target update to(dVXirLM_Rloc) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif
         else
            allocate( dxidt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
            dxidt_Rloc_container(:,:,:)=zero
            dxidt_Rloc(1:,nRstart:)   => dxidt_Rloc_container(1:lm_max,nRstart:nRstop,1)
            dVXirLM_Rloc(1:,nRstart:) => dxidt_Rloc_container(1:lm_max,nRstart:nRstop,2)
            bytes_allocated = bytes_allocated+ &
            &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dxidt_Rloc, dVXirLM_Rloc)
            !$omp target update to(dxidt_Rloc, dVXirLM_Rloc) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+ &
            &                     2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif
         end if
      else
         allocate( dxidt_Rloc_container(1,1,1:2) )
         dxidt_Rloc_container(:,:,:)=zero
         dxidt_Rloc(1:,1:)   => dxidt_Rloc_container(1:1,1:1,1)
         dVXirLM_Rloc(1:,1:) => dxidt_Rloc_container(1:1,1:1,2)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dxidt_Rloc, dVXirLM_Rloc)
         !$omp target update to(dxidt_Rloc, dVXirLM_Rloc) nowait
#endif
      end if

      if ( l_phase_field ) then
         allocate( dphidt_Rloc(lm_max,nRstart:nRstop) )
         dphidt_Rloc(:,:)=zero
         bytes_allocated = bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dphidt_Rloc)
         !$omp target update to(dphidt_Rloc) nowait
         gpu_bytes_allocated = gpu_bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( dphidt_Rloc(1:1,1:1) )
         dphidt_Rloc(:,:)=zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dphidt_Rloc)
         !$omp target update to(dphidt_Rloc) nowait
#endif
      end if

      ! The same arrays, but now the LM local part
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         if ( l_parallel_solve ) then
            if ( l_mag .and. (.not. l_mag_par_solve) ) then
               allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
               dflowdt_LMloc_container(:,:,:,:)=zero
               dbdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
               djdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
               bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
               &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: dflowdt_LMloc_container)
               !$omp target update to(dflowdt_LMloc_container) nowait
               gpu_bytes_allocated = gpu_bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
               &                     SIZEOF_DEF_COMPLEX
#endif
            end if
         else
            n_fields=3
            if ( l_mag ) n_fields=n_fields+2
            !--@> TODO: clean this ugly stuff:
            allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:n_fields,1:nexp))
            dflowdt_LMloc_container(:,:,:,:)=zero
            dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            dsdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
            bytes_allocated = bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dflowdt_LMloc_container)
            !$omp target update to(dflowdt_LMloc_container) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
            &                     SIZEOF_DEF_COMPLEX
#endif
            if ( l_mag ) then
               dbdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,4,1:nexp)
               djdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,5,1:nexp)
               bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
               &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
               gpu_bytes_allocated = gpu_bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
               &                     SIZEOF_DEF_COMPLEX
#endif
            end if
            if ( ((.not. l_double_curl) .or. l_RMS) ) then
               allocate( dpdt%expl(llm:ulm,n_r_max,nexp) )
               dpdt%expl(:,:,:)=zero
               bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*nexp* &
               &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: dpdt%expl)
               !$omp target update to(dpdt%expl) nowait
               gpu_bytes_allocated = gpu_bytes_allocated+(ulm-llm+1)*n_r_max*nexp* &
               &                     SIZEOF_DEF_COMPLEX
#endif
            else
               allocate( dpdt%expl(1,1,nexp) ) ! To avoid debug
               dpdt%expl(:,:,:)=zero
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: dpdt%expl)
               !$omp target update to(dpdt%expl) nowait
#endif
            end if
         end if
      else ! This is either high-order F.D. or Cheb
         if ( l_double_curl ) then
            allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:4,1:nexp))
            dflowdt_LMloc_container(:,:,:,:)=zero
            dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            dpdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
            dVxVhLM_LMloc(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,4,1:nexp)
            bytes_allocated = bytes_allocated+4*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dflowdt_LMloc_container)
            !$omp target update to(dflowdt_LMloc_container) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+4*(ulm-llm+1)*n_r_max*nexp* &
            &                     SIZEOF_DEF_COMPLEX
#endif
         else
            allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:3,1:nexp))
            dflowdt_LMloc_container(:,:,:,:)=zero
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dflowdt_LMloc_container)
            !$omp target update to(dflowdt_LMloc_container) nowait
#endif
            dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            dpdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
            allocate( dVxVhLM_LMloc(1,1,nexp) )
            dVxVhLM_LMloc(:,:,:)=zero
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dVxVhLM_LMloc)
            !$omp target update to(dVxVhLM_LMloc) nowait
#endif
            bytes_allocated = bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            gpu_bytes_allocated = gpu_bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
            &                     SIZEOF_DEF_COMPLEX
#endif
         end if

         allocate(dsdt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
         dsdt_LMloc_container(:,:,:,:)=zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dsdt_LMloc_container)
         !$omp target update to(dsdt_LMloc_container) nowait
#endif
         if ( .not. l_parallel_solve ) then
            dsdt%expl(llm:,1:,1:) => dsdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
         end if
         dVSrLM_LMloc(llm:,1:,1:) => dsdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
         bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
         &                     SIZEOF_DEF_COMPLEX
#endif

         allocate(dbdt_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3,1:nexp))
         dbdt_LMloc_container(:,:,:,:)=zero
         dbdt%expl(llmMag:,1:,1:) => dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,1,1:nexp)
         djdt%expl(llmMag:,1:,1:) => dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,2,1:nexp)
         dVxBhLM_LMloc(llmMag:,1:,1:) => &
         &                         dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,3,1:nexp)
         bytes_allocated = bytes_allocated+ &
         &                 3*nexp*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dbdt_LMloc_container)
         !$omp target update to(dbdt_LMloc_container) nowait
         gpu_bytes_allocated = gpu_bytes_allocated+ &
         &                 3*nexp*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
#endif
      end if

      if ( l_chemical_conv ) then
         if ( .not. l_parallel_solve ) then
            allocate(dxidt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
            dxidt_LMloc_container(:,:,:,:)=zero
            dxidt%expl(llm:,1:,1:)   => dxidt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
            dVXirLM_LMloc(llm:,1:,1:) => dxidt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
            bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU_
            !$omp target enter data map(alloc: dxidt_LMloc_container)
            !$omp target update to(dxidt_LMloc_container) nowait
            gpu_bytes_allocated = gpu_bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
            &                 SIZEOF_DEF_COMPLEX
#endif
         else
            allocate(dxidt_LMloc_container(1,1,1:2,1))
            dxidt_LMloc_container(:,:,:,:)=zero
            !dxidt%expl(1:,1:,1:)   => dxidt_LMloc_container(1:1,1:1,1,1:)
            dVXirLM_LMloc(1:,1:,1:) => dxidt_LMloc_container(1:1,1:1,2,1:)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: dxidt_LMloc_container)
            !$omp target update to(dxidt_LMloc_container) nowait
#endif
         end if
      else
         allocate(dxidt_LMloc_container(1,1,1:2,1:nexp))
         dxidt_LMloc_container(:,:,:,:)=zero
         dxidt%expl(1:,1:,1:)   => dxidt_LMloc_container(1:1,1:1,1,1:nexp)
         dVXirLM_LMloc(1:,1:,1:) => dxidt_LMloc_container(1:1,1:1,2,1:nexp)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dxidt_LMloc_container)
         !$omp target update to(dxidt_LMloc_container) nowait
#endif
      end if

      if ( .not. l_phase_field ) then
         allocate(dphidt%expl(1,1,1:nexp)) ! for debug
         dphidt%expl(:,:,:) = zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dphidt%expl)
         !$omp target update to(dphidt%expl) nowait
#endif
      end if

      ! Only when l_dt_cmb_field is requested
      ! There might be a way to allocate only when needed
      allocate ( dbdt_CMB_LMloc(llmMag:ulmMag) )
      dbdt_CMB_LMloc(:)=zero
      bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: dbdt_CMB_LMloc)
      !$omp target update to(dbdt_CMB_LMloc) nowait
      gpu_bytes_allocated = gpu_bytes_allocated+(ulmMag-llmMag+1)*SIZEOF_DEF_COMPLEX
#endif

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast
      !
      ! Memory deallocation of d?dt arrays.
      !

      if ( (.not. l_parallel_solve) .and. (.not. l_mag_par_solve) ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dflowdt_LMloc_container)
#endif
         deallocate( dflowdt_Rloc_container, dflowdt_LMloc_container )
      end if

      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dVxVhLM_Rloc)
         !$omp target exit data map(delete: dVxBhLM_Rloc)
         !$omp target exit data map(delete: dVSrLM_Rloc)
#endif
         deallocate( dVxVhLM_Rloc, dVxBhLM_Rloc, dVSrLM_Rloc)
         if (.not. l_parallel_solve ) then
#ifdef WITH_OMP_GPU
            !$omp target exit data map(delete: dpdt_Rloc)
#endif
            deallocate( dpdt_Rloc )
         end if
      else
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dbdt_LMloc_container)
#endif
         deallocate( dbdt_Rloc_container, dbdt_LMloc_container )
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dsdt_LMloc_container)
#endif
         deallocate( dsdt_Rloc_container, dsdt_LMloc_container )
         if ( .not. l_double_curl ) then
#ifdef WITH_OMP_GPU
            !$omp target exit data map(delete: dVxVhLM_Rloc)
            !$omp target exit data map(delete: dVxVhLM_LMloc)
#endif
            deallocate( dVxVhLM_Rloc, dVxVhLM_LMloc )
         end if
      end if
#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: dbdt_CMB_LMloc)
#endif
      deallocate( dbdt_CMB_LMloc )

      if ( l_chemical_conv ) then
         if ( .not. l_parallel_solve ) then
#ifdef WITH_OMP_GPU
            !$omp target exit data map(delete: dxidt_LMloc_container)
#endif
            deallocate( dxidt_Rloc_container, dxidt_LMloc_container )
         else
#ifdef WITH_OMP_GPU
            !$omp target exit data map(delete: dVXirLM_Rloc)
#endif
            deallocate( dVXirLM_Rloc )
         end if
      end if

      if ( l_phase_field ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dphidt_Rloc)
#endif
         deallocate( dphidt_Rloc )
      end if

#ifdef WITH_OMP_GPU
         !$omp target exit data map(release: lorentz_torque_ma_dt)
         !$omp target exit data map(release: lorentz_torque_ic_dt)
         !$omp target exit data map(release: domega_ma_dt)
         !$omp target exit data map(release: domega_ic_dt)
!         !$omp target exit data map(release: dzdt) !-- TODO: Error when releasing %expl (for case where expl pointes to *_contennair)
#endif
      call lorentz_torque_ma_dt%finalize()
      call lorentz_torque_ic_dt%finalize()
      call domega_ma_dt%finalize()
      call domega_ic_dt%finalize()
      call dzdt%finalize()
      if ( .not. l_double_curl .or. l_RMS ) then
#ifdef WITH_OMP_GPU
!         !$omp target exit data map(release: dpdt) !-- TODO: Error when releasing %expl
#endif
         call dpdt%finalize()
      end if
#ifdef WITH_OMP_GPU
!         !$omp target exit data map(release: dwdt) !-- TODO: Error when releasing %expl
#endif
      call dwdt%finalize()
      if ( l_heat ) then
#ifdef WITH_OMP_GPU
!         !$omp target exit data map(release: dsdt) !-- TODO: Error when releasing %expl
#endif
         call dsdt%finalize()
      end if
      if ( l_chemical_conv ) then
#ifdef WITH_OMP_GPU
!         !$omp target exit data map(release: dxidt) !-- TODO: Error when releasing %expl
#endif
         call dxidt%finalize()
      end if
      if ( l_phase_field ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(release: dphidt)
#endif
         call dphidt%finalize()
      end if
      if ( l_mag ) then
#ifdef WITH_OMP_GPU
!         !$omp target exit data map(release: dbdt) !-- TODO: Error when releasing %expl
!         !$omp target exit data map(release: djdt) !-- TODO: Error when releasing %expl
#endif
         call dbdt%finalize()
         call djdt%finalize()
      end if
      if ( l_cond_ic ) then
#ifdef WITH_OMP_GPU
!         !$omp target exit data map(release: djdt_ic) !-- TODO: Error when releasing %expl
!         !$omp target exit data map(release: djdt_ic) !-- TODO: Error when releasing %expl
#endif
         call dbdt_ic%finalize()
         call djdt_ic%finalize()
      end if

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
