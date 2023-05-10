module fields
   !
   ! This module contains all the fields used in MagIC in the hybrid (LM,r)
   ! space as well as their radial derivatives. It defines both the
   ! LM-distributed arrays and the R-distributed arrays.
   !
   use precision_mod
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use constants, only: zero
   use physical_parameters, only: ampForce
   use truncation, only: lm_max, n_r_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag, fd_order, fd_order_bound
   use logic, only: l_chemical_conv, l_finite_diff, l_mag, l_parallel_solve, &
       &            l_mag_par_solve, l_phase_field
   use blocking, only: llm, ulm, llmMag, ulmMag
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag
   use parallel_mod, only: rank

   implicit none

   private

   !-- Velocity potentials:
   complex(cp), public, allocatable, target :: flow_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: flow_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target :: press_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: press_Rloc_container(:,:,:)
   complex(cp), public, pointer :: w_LMloc(:,:),dw_LMloc(:,:),ddw_LMloc(:,:)
   complex(cp), public, pointer :: w_Rloc(:,:), dw_Rloc(:,:), ddw_Rloc(:,:)

   complex(cp), public, pointer :: z_LMloc(:,:),dz_LMloc(:,:)
   complex(cp), public, pointer :: z_Rloc(:,:), dz_Rloc(:,:)

   !-- Entropy:
   complex(cp), public, allocatable, target :: s_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: s_Rloc_container(:,:,:)
   complex(cp), public, pointer :: s_LMloc(:,:), ds_LMloc(:,:)
   complex(cp), public, pointer :: s_Rloc(:,:), ds_Rloc(:,:)

   !-- Chemical composition:
   complex(cp), public, allocatable, target :: xi_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: xi_Rloc_container(:,:,:)
   complex(cp), public, pointer :: xi_LMloc(:,:), dxi_LMloc(:,:)
   complex(cp), public, pointer :: xi_Rloc(:,:), dxi_Rloc(:,:)

   !-- Phase field
   complex(cp), public, allocatable :: phi_LMloc(:,:), phi_Rloc(:,:)

   !-- Pressure:
   complex(cp), public, pointer :: p_LMloc(:,:), dp_LMloc(:,:)
   complex(cp), public, pointer :: p_Rloc(:,:), dp_Rloc(:,:)

   !-- Magnetic field potentials:
   complex(cp), public, allocatable :: bICB(:)
   complex(cp), public, allocatable, target :: field_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: field_Rloc_container(:,:,:)
   complex(cp), public, pointer :: b_LMloc(:,:), db_LMloc(:,:), ddb_LMloc(:,:)
   complex(cp), public, pointer :: b_Rloc(:,:), db_Rloc(:,:), ddb_Rloc(:,:)
   complex(cp), public, pointer :: aj_LMloc(:,:), dj_LMloc(:,:), ddj_LMloc(:,:)
   complex(cp), public, pointer :: aj_Rloc(:,:), dj_Rloc(:,:)
   complex(cp), public, allocatable :: ddj_Rloc(:,:)

   !-- Magnetic field potentials in inner core:
   !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
   !         for even chebs
   complex(cp), public, allocatable :: b_ic(:,:)
   complex(cp), public, allocatable :: db_ic(:,:)
   complex(cp), public, allocatable :: ddb_ic(:,:)
   complex(cp), public, allocatable :: aj_ic(:,:)
   complex(cp), public, allocatable :: dj_ic(:,:)
   complex(cp), public, allocatable :: b_ic_LMloc(:,:)
   complex(cp), public, allocatable :: db_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddb_ic_LMloc(:,:)
   complex(cp), public, allocatable :: aj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: dj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: bodyForce(:,:)
   complex(cp), public, allocatable :: work_LMloc(:,:) ! Needed in update routines

   !-- Rotation rates:
   real(cp), public :: omega_ic,omega_ma

   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields
      !
      ! This subroutine allocates the different fields used in MagIC
      !

      integer :: n_fields

      !-- Velocity potentials:
      if ( rank == 0 ) then
         allocate( bICB(lm_maxMag) )
         bytes_allocated = bytes_allocated + lm_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
         &                 5*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + &
         &                 5*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( bICB(1) )
         bytes_allocated = bytes_allocated + n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) )
         allocate( dj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 5*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 5*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
#endif
      end if
      bICB(:)    =zero
      b_ic(:,:)  =zero
      db_ic(:,:) =zero
      ddb_ic(:,:)=zero
      aj_ic(:,:) =zero
      dj_ic(:,:) =zero
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: b_ic, db_ic, ddb_ic, aj_ic, dj_ic)
      !$omp target update to(b_ic, db_ic, ddb_ic, aj_ic, dj_ic) nowait
#endif

      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         if ( l_parallel_solve ) then
            allocate(w_LMloc(llm:ulm,n_r_max), z_LMloc(llm:ulm,n_r_max))
            allocate(s_LMloc(llm:ulm,n_r_max))
            w_LMloc(:,:)=zero
            z_LMloc(:,:)=zero
            s_LMloc(:,:)=zero
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: w_LMloc, z_LMloc, s_LMloc)
            !$omp target update to(w_LMloc, z_LMloc, s_LMloc) nowait
#endif
            if ( l_mag ) then
               if ( l_mag_par_solve ) then
                  allocate(aj_LMloc(llm:ulm,n_r_max), b_LMloc(llm:ulm,n_r_max))
                  aj_LMloc(:,:)=zero
                  b_LMloc(:,:) =zero
#ifdef WITH_OMP_GPU
                  !$omp target enter data map(alloc: aj_LMloc, b_LMloc)
                  !$omp target update to(aj_LMloc, b_LMloc) nowait
#endif
               else
                  allocate( flow_LMloc_container(llm:ulm,n_r_max,1:2) )
                  flow_LMloc_container(:,:,:)=zero
                  b_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,1)
                  aj_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,2)
#ifdef WITH_OMP_GPU
                  !$omp target enter data map(alloc: aj_LMloc, b_LMloc)
                  !$omp target update to(aj_LMloc, b_LMloc) nowait
#endif
               end if
            else
               allocate ( b_LMloc(1,1), aj_LMloc(1,1) )
               b_LMloc(1,1)=zero; aj_LMloc(1,1)=zero
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: b_LMloc, aj_LMloc)
               !$omp target update to(b_LMloc, aj_LMloc) nowait
#endif
            end if
         else
            n_fields = 3
            if ( l_mag ) n_fields = n_fields+2
            allocate( flow_LMloc_container(llm:ulm,n_r_max,1:n_fields) )
            flow_LMloc_container(:,:,:)=zero
            w_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,1)
            z_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,2)
            s_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,3)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: w_LMloc, z_LMloc, s_LMloc)
            !$omp target update to(w_LMloc, z_LMloc, s_LMloc) nowait
#endif
            if ( l_mag ) then
               b_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,4)
               aj_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,5)
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: b_LMloc, aj_LMloc)
               !$omp target update to(b_LMloc, aj_LMloc) nowait
#endif
            end if
         end if

         allocate(dw_LMloc(llm:ulm,n_r_max), ddw_LMloc(llm:ulm,n_r_max))
         dw_LMloc(:,:) =zero
         ddw_LMloc(:,:)=zero
         allocate(dz_LMloc(llm:ulm,n_r_max), ds_LMloc(llm:ulm,n_r_max))
         dz_LMloc(:,:) =zero
         ds_LMloc(:,:) =zero
         allocate(db_LMloc(llmMag:ulmMag,n_r_maxMag))
         db_LMloc(:,:) =zero
         allocate(ddb_LMloc(llmMag:ulmMag,n_r_maxMag))
         ddb_LMloc(:,:)=zero
         allocate(dj_LMloc(llmMag:ulmMag,n_r_maxMag))
         dj_LMloc(:,:) =zero
         allocate(ddj_LMloc(llmMag:ulmMag,n_r_maxMag))
         ddj_LMloc(:,:)=zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dw_LMloc, ddw_LMloc, dz_LMloc, ds_LMloc, &
         !$omp&                             db_LMloc, ddb_LMloc, dj_LMloc, ddj_LMloc)
         !$omp target update to(dw_LMloc, ddw_LMloc, dz_LMloc, ds_LMloc, &
         !$omp&                 db_LMloc, ddb_LMloc, dj_LMloc, ddj_LMloc) nowait
#endif

         if ( l_parallel_solve ) then
            allocate(w_Rloc(lm_max,nRstart:nRstop), z_Rloc(lm_max,nRstart:nRstop))
            allocate(s_Rloc(lm_max,nRstart:nRstop))
            w_Rloc(:,:)=zero
            z_Rloc(:,:)=zero
            s_Rloc(:,:)=zero
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: s_Rloc, w_Rloc, z_Rloc)
            !$omp target update to(s_Rloc, w_Rloc, z_Rloc) nowait
#endif
            if ( l_mag ) then
               if( l_mag_par_solve ) then
                  allocate(b_Rloc(lm_max,nRstart:nRstop), aj_Rloc(lm_max,nRstart:nRstop))
                  b_Rloc(:,:) =zero
                  aj_Rloc(:,:)=zero
#ifdef WITH_OMP_GPU
                  !$omp target enter data map(alloc: b_RLoc, aj_Rloc)
                  !$omp target update to(b_RLoc, aj_Rloc) nowait
#endif
               else
                  allocate( flow_Rloc_container(1:lm_max,nRstart:nRstop,1:2) )
                  flow_Rloc_container(:,:,:)=zero
                  b_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,1)
                  aj_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,2)
#ifdef WITH_OMP_GPU
                  !$omp target enter data map(alloc: b_RLoc, aj_Rloc)
                  !$omp target update to(b_RLoc, aj_Rloc) nowait
#endif
               end if
            else
               allocate ( b_Rloc(1,1), aj_Rloc(1,1) )
               b_Rloc(1,1)=zero; aj_Rloc(1,1)=zero
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc:b_RLoc, aj_Rloc)
               !$omp target update to(b_RLoc, aj_Rloc) nowait
#endif
            end if
         else
            allocate( flow_Rloc_container(1:lm_max,nRstart:nRstop,1:n_fields) )
            flow_Rloc_container(:,:,:)=zero
            w_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,1)
            z_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,2)
            s_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,3)
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: s_Rloc, w_Rloc, z_Rloc)
            !$omp target update to(s_Rloc, w_Rloc, z_Rloc) nowait
#endif
            if ( l_mag ) then
               b_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,4)
               aj_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,5)
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: b_RLoc, aj_Rloc)
               !$omp target update to(b_RLoc, aj_Rloc) nowait
#endif
            end if
         end if
         allocate(dw_Rloc(lm_max,nRstart:nRstop), ddw_Rloc(lm_max,nRstart:nRstop))
         dw_Rloc(:,:) =zero
         ddw_Rloc(:,:)=zero
         allocate(dz_Rloc(lm_max,nRstart:nRstop), ds_Rloc(lm_max,nRstart:nRstop))
         dz_Rloc(:,:) =zero
         ds_Rloc(:,:) =zero
         allocate(db_Rloc(lm_maxMag,nRstartMag:nRstopMag))
         db_Rloc(:,:) =zero
         allocate(ddb_Rloc(lm_maxMag,nRstartMag:nRstopMag))
         ddb_Rloc(:,:)=zero
         allocate(dj_Rloc(lm_maxMag,nRstartMag:nRstopMag))
         dj_Rloc(:,:) =zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: ds_Rloc, dw_Rloc, ddw_Rloc, dz_Rloc, db_Rloc, ddb_Rloc, dj_Rloc)
         !$omp target update to(ds_Rloc, dw_Rloc, ddw_Rloc, dz_Rloc, db_Rloc, ddb_Rloc, dj_Rloc) nowait
#endif
      else
         allocate( flow_LMloc_container(llm:ulm,n_r_max,1:5) )
         flow_LMloc_container(:,:,:)=zero
         w_LMloc(llm:,1:)   => flow_LMloc_container(llm:ulm,1:n_r_max,1)
         dw_LMloc(llm:,1:)  => flow_LMloc_container(llm:ulm,1:n_r_max,2)
         ddw_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,3)
         z_LMloc(llm:,1:)   => flow_LMloc_container(llm:ulm,1:n_r_max,4)
         dz_LMloc(llm:,1:)  => flow_LMloc_container(llm:ulm,1:n_r_max,5)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: w_LMloc, dw_LMloc, ddw_LMloc, z_LMloc, dz_LMloc)
         !$omp target update to(w_LMloc, dw_LMloc, ddw_LMloc, z_LMloc, dz_LMloc) nowait
#endif

         allocate( flow_Rloc_container(lm_max,nRstart:nRstop,1:5) )
         flow_Rloc_container(:,:,:)=zero
         w_Rloc(1:,nRstart:)   => flow_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dw_Rloc(1:,nRstart:)  => flow_Rloc_container(1:lm_max,nRstart:nRstop,2)
         ddw_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,3)
         z_Rloc(1:,nRstart:)   => flow_Rloc_container(1:lm_max,nRstart:nRstop,4)
         dz_Rloc(1:,nRstart:)  => flow_Rloc_container(1:lm_max,nRstart:nRstop,5)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc)
         !$omp target update to(w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc) nowait
#endif

         !-- Entropy:
         allocate( s_LMloc_container(llm:ulm,n_r_max,1:2) )
         s_LMloc_container(:,:,:)=zero
         s_LMloc(llm:,1:)  => s_LMloc_container(llm:ulm,1:n_r_max,1)
         ds_LMloc(llm:,1:) => s_LMloc_container(llm:ulm,1:n_r_max,2)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: s_LMloc, ds_LMloc)
         !$omp target update to(s_LMloc, ds_LMloc) nowait
#endif
         allocate( s_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         s_Rloc_container(:,:,:)=zero
         s_Rloc(1:,nRstart:)  => s_Rloc_container(1:lm_max,nRstart:nRstop,1)
         ds_Rloc(1:,nRstart:) => s_Rloc_container(1:lm_max,nRstart:nRstop,2)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: s_Rloc, ds_Rloc)
         !$omp target update to(s_Rloc, ds_Rloc) nowait
#endif

         !-- Magnetic field potentials:
         allocate( field_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:6) )
         field_LMloc_container(:,:,:)=zero
         b_LMloc(llmMag:,1:)   => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,1)
         db_LMloc(llmMag:,1:)  => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,2)
         ddb_LMloc(llmMag:,1:) => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,3)
         aj_LMloc(llmMag:,1:)  => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,4)
         dj_LMloc(llmMag:,1:)  => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,5)
         ddj_LMloc(llmMag:,1:) => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,6)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: b_LMloc, db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc, ddj_LMloc)
         !$omp target update to(b_LMloc, db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc, ddj_LMloc) nowait
#endif

         allocate( field_Rloc_container(lm_maxMag,nRstart:nRstop,1:5) )
         field_Rloc_container(:,:,:)=zero
         b_Rloc(1:,nRstart:)   => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,1)
         db_Rloc(1:,nRstart:)  => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,2)
         ddb_Rloc(1:,nRstart:) => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,3)
         aj_Rloc(1:,nRstart:)  => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,4)
         dj_Rloc(1:,nRstart:)  => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,5)
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: b_Rloc, db_Rloc, ddb_Rloc, aj_Rloc, dj_Rloc)
         !$omp target update to(b_Rloc, db_Rloc, ddb_Rloc, aj_Rloc, dj_Rloc) nowait
#endif
      end if

      if ( l_mag_par_solve ) then
         allocate(ddj_Rloc(lm_maxMag,nRstartMag:nRstopMag))
         ddj_Rloc(:,:)=zero
         bytes_allocated = bytes_allocated+(nRstopMag-nRstartMag+1)*lm_maxMag* &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: ddj_Rloc)
         !$omp target update to(ddj_Rloc) nowait
         gpu_bytes_allocated = gpu_bytes_allocated+(nRstopMag-nRstartMag+1)*lm_maxMag* &
         &                 SIZEOF_DEF_COMPLEX
#endif
      end if

      allocate( press_LMloc_container(llm:ulm,n_r_max,1:2) )
      press_LMloc_container(:,:,:)=zero
      p_LMloc(llm:,1:)   => press_LMloc_container(llm:ulm,1:n_r_max,1)
      dp_LMloc(llm:,1:)  => press_LMloc_container(llm:ulm,1:n_r_max,2)
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: p_LMloc, dp_LMloc)
      !$omp target update to(p_LMloc, dp_LMloc) nowait
#endif

      allocate( press_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      press_Rloc_container(:,:,:)=zero
      p_Rloc(1:,nRstart:)   => press_Rloc_container(1:lm_max,nRstart:nRstop,1)
      dp_Rloc(1:,nRstart:)  => press_Rloc_container(1:lm_max,nRstart:nRstop,2)
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: p_Rloc, dp_Rloc)
      !$omp target update to(p_Rloc, dp_Rloc) nowait
#endif

      bytes_allocated = bytes_allocated + &
      &                 9*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
      &                 9*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
      &                 6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
      &                 5*lm_maxMag*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
      gpu_bytes_allocated = gpu_bytes_allocated + &
      &                     9*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      gpu_bytes_allocated = gpu_bytes_allocated + &
      &                     9*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      gpu_bytes_allocated = gpu_bytes_allocated + &
      &                     6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
      gpu_bytes_allocated = gpu_bytes_allocated + &
      &                     5*lm_maxMag*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMloc_container(llm:ulm,n_r_max,1:2) )
         xi_LMloc_container(:,:,:)=zero
         xi_LMloc(llm:,1:)  => xi_LMloc_container(llm:ulm,1:n_r_max,1)
         dxi_LMloc(llm:,1:) => xi_LMloc_container(llm:ulm,1:n_r_max,2)
         allocate( xi_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         xi_Rloc_container(:,:,:)=zero
         xi_Rloc(1:,nRstart:)  => xi_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dxi_Rloc(1:,nRstart:) => xi_Rloc_container(1:lm_max,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated + &
         &                 2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
         &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + &
         &                 2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         gpu_bytes_allocated = gpu_bytes_allocated + &
         &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( xi_LMloc_container(1,1,2) ) ! For debugging
         xi_LMloc_container(:,:,:)=zero
         xi_LMloc(1:,1:)  => xi_LMloc_container(1:1,1:1,1)
         dxi_LMloc(1:,1:) => xi_LMloc_container(1:1,1:1,2)
         allocate( xi_Rloc_container(1,1,2) )
         xi_Rloc_container(:,:,:)=zero
         xi_Rloc(1:,1:)   => xi_Rloc_container(1:1,1:1,1)
         dxi_Rloc(1:,1:)  => xi_Rloc_container(1:1,1:1,2)
      end if
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: xi_Rloc, dxi_Rloc, xi_LMloc, dxi_LMloc)
      !$omp target update to(xi_Rloc, dxi_Rloc, xi_LMloc, dxi_LMloc) nowait
#endif

      !-- Phase field
      if ( l_phase_field ) then
         allocate( phi_LMloc(llm:ulm,1:n_r_max) )
         phi_LMloc(:,:)=zero
         bytes_allocated = bytes_allocated + &
         &                 (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         allocate( phi_Rloc(1:lm_max,nRstart:nRstop) )
         phi_Rloc(:,:)=zero
         bytes_allocated = bytes_allocated + &
         &                 (nRstop-nRstart+1)*lm_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + &
         &                 (nRstop-nRstart+1)*lm_max*SIZEOF_DEF_COMPLEX
#endif
      else ! For debugging
         allocate( phi_LMloc(1:1,1:1), phi_Rloc(1:1,1:1) )
         phi_LMloc(:,:)=zero
         phi_Rloc(:,:)=zero
      end if
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: phi_Rloc, phi_LMloc)
      !$omp target update to(phi_Rloc, phi_LMloc) nowait
#endif


      !-- Magnetic field potentials in inner core:
      !   NOTE: n_r-dimension may be smaller once CHEBFT is adopted
      !         for even chebs
      allocate( b_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( db_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddb_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( aj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( dj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      b_ic_LMloc(:,:)  =zero
      db_ic_LMloc(:,:) =zero
      ddb_ic_LMloc(:,:)=zero
      aj_ic_LMloc(:,:) =zero
      dj_ic_LMloc(:,:) =zero
      ddj_ic_LMloc(:,:)=zero
      bytes_allocated = bytes_allocated + &
      &                 6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
      gpu_bytes_allocated = gpu_bytes_allocated + &
      &                 6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
#endif

      allocate( work_LMloc(llm:ulm,1:n_r_max) )
      work_LMloc(:,:)=zero
      bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
      gpu_bytes_allocated = gpu_bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
#endif

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, &
      !$omp&                             aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc, work_LMloc)
      !$omp target update to(b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, &
      !$omp&                 aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc, work_LMloc) nowait
#endif

      if (ampForce /= 0.0_cp) then
         allocate(bodyForce(llm:ulm,n_r_max))
         bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         bodyForce(:,:) = zero
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated=gpu_bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         !$omp target enter data map(alloc: bodyForce)
         !$omp target update to(bodyForce)
#endif
      end if

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields
      !
      ! This subroutine deallocates the field arrays used in MagIC
      !

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: s_Rloc, s_LMloc)
      !$omp target exit data map(delete: p_Rloc, p_LMloc, dp_LMloc, dp_Rloc)
      !$omp target exit data map(delete: xi_Rloc, xi_LMloc)
      !$omp target exit data map(delete: w_Rloc, w_LMloc)
      !$omp target exit data map(delete: z_Rloc, z_LMloc)
      !$omp target exit data map(delete: b_Rloc, b_LMloc)
      !$omp target exit data map(delete: aj_Rloc, aj_LMloc)
#endif

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: b_ic, db_ic, ddb_ic, aj_ic, dj_ic)
#endif
      deallocate( bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic )
      deallocate( press_LMloc_container, press_Rloc_container )
      if ( l_parallel_solve ) then
         deallocate( w_LMloc, z_LMloc, s_LMloc, w_RLoc, z_Rloc, s_Rloc )
         if ( l_mag ) then
            if ( l_mag_par_solve ) then
               deallocate( b_LMloc, aj_LMloc, b_RLoc, aj_Rloc )
            else
               deallocate( flow_Rloc_container, flow_LMloc_container )
            end if
         end if
      else
         deallocate( flow_Rloc_container, flow_LMloc_container )
      end if
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: ds_Rloc, dw_Rloc, ddw_Rloc, dz_Rloc, &
         !$omp&                             db_Rloc, ddb_Rloc, dj_Rloc,          &
         !$omp&                             dw_LMloc, ddw_LMloc, dz_LMloc, ds_LMloc, &
         !$omp&                             db_LMloc, ddb_LMloc, dj_LMloc, ddj_LMloc)
#endif
         deallocate( dw_LMloc, ddw_LMloc, dz_LMloc, ds_LMloc)
         deallocate( db_LMloc, ddb_LMloc, dj_LMloc, ddj_LMloc)
         deallocate( dw_Rloc, ddw_Rloc, dz_Rloc, ds_Rloc)
         deallocate( db_Rloc, ddb_Rloc, dj_Rloc)
      else
#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: ds_Rloc, ds_LMloc)
#endif
         deallocate( s_LMloc_container, s_Rloc_container )
         deallocate( field_LMloc_container, field_Rloc_container )
      end if
#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc, &
      !$omp&                             dj_ic_LMloc, ddj_ic_LMloc, &
      !$omp&                             work_LMloc, &
      !$omp&                             phi_LMloc, phi_Rloc)
#endif
      deallocate( b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc )
      deallocate( dj_ic_LMloc, ddj_ic_LMloc )
      deallocate( xi_LMloc_container, xi_Rloc_container )
      deallocate( work_LMloc )
      deallocate( phi_LMloc, phi_Rloc )
      if ( l_mag_par_solve ) then
#ifdef WITH_OMP_GPU_
      !$omp target exit data map(delete: ddj_Rloc)
#endif
         deallocate(ddj_Rloc)
      end if
      if ( ampForce /= 0.0_cp ) then
#ifdef WITH_OMP_GPU_
      !$omp target exit data map(delete: bodyForce)
#endif
         deallocate(bodyForce)
      end if

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
