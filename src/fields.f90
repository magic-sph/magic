module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_r_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag
   use logic, only: l_chemical_conv, l_mag_hel
   use blocking, only: llm, ulm, llmMag, ulmMag
   use radial_data, only: nRstart, nRstop
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

   !-- Pressure:
   complex(cp), public, pointer :: p_LMloc(:,:), dp_LMloc(:,:)
   complex(cp), public, pointer :: p_Rloc(:,:), dp_Rloc(:,:)

   !-- Magnetic field potentials:
   complex(cp), public, allocatable :: b(:,:)
   complex(cp), public, allocatable, target :: field_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: field_Rloc_container(:,:,:)
   complex(cp), public, pointer :: b_LMloc(:,:), db_LMloc(:,:), ddb_LMloc(:,:)
   complex(cp), public, pointer :: b_Rloc(:,:), db_Rloc(:,:), ddb_Rloc(:,:)
   complex(cp), public, pointer :: aj_LMloc(:,:), dj_LMloc(:,:), ddj_LMloc(:,:)
   complex(cp), public, pointer :: aj_Rloc(:,:), dj_Rloc(:,:)

   !-- Magnetic field potentials in inner core:
   !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
   !         for even chebs
   complex(cp), public, allocatable :: b_ic(:,:)
   complex(cp), public, allocatable :: db_ic(:,:)
   complex(cp), public, allocatable :: ddb_ic(:,:)
   complex(cp), public, allocatable :: aj_ic(:,:)
   complex(cp), public, allocatable :: dj_ic(:,:)
   complex(cp), public, allocatable :: ddj_ic(:,:)
   complex(cp), public, allocatable :: b_ic_LMloc(:,:)
   complex(cp), public, allocatable :: db_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddb_ic_LMloc(:,:)
   complex(cp), public, allocatable :: aj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: dj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddj_ic_LMloc(:,:)

   complex(cp), public, allocatable :: work_LMloc(:,:) ! Needed in update routines

   complex(cp), public, allocatable :: zeros_alike_b_Rloc(:) ! trick magnetic helicity

   !-- Rotation rates:
   real(cp), public :: omega_ic,omega_ma

   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      !-- Velocity potentials:
      if ( rank == 0 ) then
         allocate( b(lm_maxMag,n_r_maxMag) )
         bytes_allocated = bytes_allocated +  &
         &                 lm_maxMag*n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
         &                 6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      else
         allocate( b(1,n_r_maxMag) )
         bytes_allocated = bytes_allocated + n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) )
         allocate( dj_ic(1,n_r_ic_maxMag) )
         allocate( ddj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      end if
      allocate( flow_LMloc_container(llm:ulm,n_r_max,1:5) )
      w_LMloc(llm:,1:)   => flow_LMloc_container(llm:ulm,1:n_r_max,1)
      dw_LMloc(llm:,1:)  => flow_LMloc_container(llm:ulm,1:n_r_max,2)
      ddw_LMloc(llm:,1:) => flow_LMloc_container(llm:ulm,1:n_r_max,3)
      z_LMloc(llm:,1:)   => flow_LMloc_container(llm:ulm,1:n_r_max,4)
      dz_LMloc(llm:,1:)  => flow_LMloc_container(llm:ulm,1:n_r_max,5)

      allocate( press_LMloc_container(llm:ulm,n_r_max,1:2) )
      p_LMloc(llm:,1:)   => press_LMloc_container(llm:ulm,1:n_r_max,1)
      dp_LMloc(llm:,1:)  => press_LMloc_container(llm:ulm,1:n_r_max,2)

      allocate( flow_Rloc_container(lm_max,nRstart:nRstop,1:5) )
      w_Rloc(1:,nRstart:)   => flow_Rloc_container(1:lm_max,nRstart:nRstop,1)
      dw_Rloc(1:,nRstart:)  => flow_Rloc_container(1:lm_max,nRstart:nRstop,2)
      ddw_Rloc(1:,nRstart:) => flow_Rloc_container(1:lm_max,nRstart:nRstop,3)
      z_Rloc(1:,nRstart:)   => flow_Rloc_container(1:lm_max,nRstart:nRstop,4)
      dz_Rloc(1:,nRstart:)  => flow_Rloc_container(1:lm_max,nRstart:nRstop,5)

      allocate( press_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      p_Rloc(1:,nRstart:)   => press_Rloc_container(1:lm_max,nRstart:nRstop,1)
      dp_Rloc(1:,nRstart:)  => press_Rloc_container(1:lm_max,nRstart:nRstop,2)

      !-- Entropy:
      allocate( s_LMloc_container(llm:ulm,n_r_max,1:2) )
      s_LMloc(llm:,1:)  => s_LMloc_container(llm:ulm,1:n_r_max,1)
      ds_LMloc(llm:,1:) => s_LMloc_container(llm:ulm,1:n_r_max,2)
      allocate( s_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      s_Rloc(1:,nRstart:)  => s_Rloc_container(1:lm_max,nRstart:nRstop,1)
      ds_Rloc(1:,nRstart:) => s_Rloc_container(1:lm_max,nRstart:nRstop,2)

      bytes_allocated = bytes_allocated + &
      &                 9*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
      &                 9*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMloc_container(llm:ulm,n_r_max,1:2) )
         xi_LMloc(llm:,1:)  => xi_LMloc_container(llm:ulm,1:n_r_max,1)
         dxi_LMloc(llm:,1:) => xi_LMloc_container(llm:ulm,1:n_r_max,2)
         allocate( xi_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         xi_Rloc(1:,nRstart:)  => xi_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dxi_Rloc(1:,nRstart:) => xi_Rloc_container(1:lm_max,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated + &
         &                 2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
         &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_LMloc_container(1,1,2) ) ! For debugging
         xi_LMloc(1:,1:)  => xi_LMloc_container(1:1,1:1,1)
         dxi_LMloc(1:,1:) => xi_LMloc_container(1:1,1:1,2)
         allocate( xi_Rloc_container(1,1,2) )
         xi_Rloc(1:,1:)   => xi_Rloc_container(1:1,1:1,1)
         dxi_Rloc(1:,1:)  => xi_Rloc_container(1:1,1:1,2)
      end if

      !-- Magnetic field potentials:
      allocate( field_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:6) )
      b_LMloc(llmMag:,1:)   => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,1)
      db_LMloc(llmMag:,1:)  => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,2)
      ddb_LMloc(llmMag:,1:) => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,3)
      aj_LMloc(llmMag:,1:)  => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,4)
      dj_LMloc(llmMag:,1:)  => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,5)
      ddj_LMloc(llmMag:,1:) => field_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,6)

      allocate( field_Rloc_container(lm_maxMag,nRstart:nRstop,1:5) )
      b_Rloc(1:,nRstart:)   => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,1)
      db_Rloc(1:,nRstart:)  => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,2)
      ddb_Rloc(1:,nRstart:) => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,3)
      aj_Rloc(1:,nRstart:)  => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,4)
      dj_Rloc(1:,nRstart:)  => field_Rloc_container(1:lm_maxMag,nRstart:nRstop,5)

      bytes_allocated = bytes_allocated + &
      &                 6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
      &                 5*lm_maxMag*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      if (l_mag_hel) then
         allocate(zeros_alike_b_Rloc(1:lm_maxMag))
         bytes_allocated = bytes_allocated+lm_maxMag*SIZEOF_DEF_COMPLEX

         zeros_alike_b_Rloc(:) = zero
      end if


      !-- Magnetic field potentials in inner core:
      !   NOTE: n_r-dimension may be smaller once CHEBFT is adopted
      !         for even chebs
      allocate( b_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( db_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddb_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( aj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( dj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )

      bytes_allocated = bytes_allocated + &
      &                 6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      allocate( work_LMloc(llm:ulm,1:n_r_max) )
      bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      deallocate( b, b_ic, db_ic, ddb_ic )
      deallocate( aj_ic, dj_ic, ddj_ic, flow_LMloc_container )
      deallocate( press_LMloc_container, press_Rloc_container )
      deallocate( flow_Rloc_container, s_LMloc_container, s_Rloc_container )
      deallocate( field_LMloc_container, field_Rloc_container )
      deallocate( b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc )
      deallocate( dj_ic_LMloc, ddj_ic_LMloc )
      deallocate( xi_LMloc_container, xi_Rloc_container )
      deallocate( work_LMloc )

      if (l_mag_hel) deallocate(zeros_alike_b_Rloc)

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
