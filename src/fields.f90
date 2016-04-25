module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_r_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag
   use logic, only: l_chemical_conv
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use radial_data, only: nRstart, nRstop
   use parallel_mod, only: rank
 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), public, allocatable :: w(:,:)
   complex(cp), public, allocatable, target :: flow_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: flow_Rloc_container(:,:,:)
   complex(cp), public, pointer :: w_LMloc(:,:),dw_LMloc(:,:),ddw_LMloc(:,:)
   complex(cp), public, pointer :: w_Rloc(:,:), dw_Rloc(:,:), ddw_Rloc(:,:)
 
   complex(cp), public, allocatable :: z(:,:)
   complex(cp), public, pointer :: z_LMloc(:,:),dz_LMloc(:,:)
   complex(cp), public, pointer :: z_Rloc(:,:), dz_Rloc(:,:)
 
   !-- Entropy:
   complex(cp), public, allocatable :: s(:,:)
   complex(cp), public, allocatable, target :: s_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: s_Rloc_container(:,:,:)
   complex(cp), public, pointer :: s_LMloc(:,:), ds_LMloc(:,:)
   complex(cp), public, pointer :: s_Rloc(:,:), ds_Rloc(:,:)
 
   !-- Chemical composition:
   complex(cp), public, allocatable :: xi(:,:)
   complex(cp), public, allocatable, target :: xi_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: xi_Rloc_container(:,:,:)
   complex(cp), public, pointer :: xi_LMloc(:,:), dxi_LMloc(:,:)
   complex(cp), public, pointer :: xi_Rloc(:,:), dxi_Rloc(:,:)

   !-- Pressure:
   complex(cp), public, allocatable :: p(:,:)
   complex(cp), public, pointer :: p_LMloc(:,:), dp_LMloc(:,:)
   complex(cp), public, pointer :: p_Rloc(:,:), dp_Rloc(:,:)
 
   !-- Magnetic field potentials:
   complex(cp), public, allocatable :: b(:,:)
   complex(cp), public, allocatable :: db(:,:)
   complex(cp), public, allocatable :: aj(:,:)
   complex(cp), public, allocatable :: dj(:,:)
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
 
   !-- Rotation rates:
   real(cp), public :: omega_ic,omega_ma

   public :: initialize_fields

contains

   subroutine initialize_fields

      !-- Velocity potentials:
      if ( rank == 0 ) then
         allocate( w(lm_max,n_r_max) )
         allocate( z(lm_max,n_r_max) )
         allocate( s(lm_max,n_r_max) )
         allocate( p(lm_max,n_r_max) )
         bytes_allocated = bytes_allocated + 4*lm_max*n_r_max*SIZEOF_DEF_COMPLEX

         if ( l_chemical_conv ) then
            allocate ( xi(lm_max,n_r_max) )
            bytes_allocated = bytes_allocated+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         end if

         allocate( b(lm_maxMag,n_r_maxMag) )
         allocate( db(lm_maxMag,n_r_maxMag) )
         allocate( aj(lm_maxMag,n_r_maxMag) )
         allocate( dj(lm_maxMag,n_r_maxMag) )
         bytes_allocated = bytes_allocated +  &
                           4*lm_maxMag*n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )  
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
                           6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      else
         allocate( w(1,n_r_max) )
         allocate( z(1,n_r_max) )
         allocate( s(1,n_r_max) )
         allocate( p(1,n_r_max) )
         bytes_allocated = bytes_allocated + 4*n_r_max*SIZEOF_DEF_COMPLEX
         if ( l_chemical_conv ) then
            allocate ( xi(1,n_r_max) )
            bytes_allocated = bytes_allocated + n_r_max*SIZEOF_DEF_COMPLEX
         end if
         allocate( b(1,n_r_maxMag) )
         allocate( db(1,n_r_maxMag) )
         allocate( aj(1,n_r_maxMag) )
         allocate( dj(1,n_r_maxMag) )
         bytes_allocated = bytes_allocated + 4*n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )  
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) ) 
         allocate( dj_ic(1,n_r_ic_maxMag) )
         allocate( ddj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      end if
      allocate( flow_LMloc_container(llm:ulm,n_r_max,1:7) )
      w_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,1)
      dw_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,2)
      ddw_LMloc(llm:ulm,1:n_r_max) => flow_LMloc_container(:,:,3)
      z_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,4)
      dz_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,5)
      p_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,6)
      dp_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,7)

      allocate( flow_Rloc_container(lm_max,nRstart:nRstop,1:7) )
      w_Rloc(1:lm_max,nRstart:nRstop)   => flow_Rloc_container(:,:,1)
      dw_Rloc(1:lm_max,nRstart:nRstop)  => flow_Rloc_container(:,:,2)
      ddw_Rloc(1:lm_max,nRstart:nRstop) => flow_Rloc_container(:,:,3)
      z_Rloc(1:lm_max,nRstart:nRstop)   => flow_Rloc_container(:,:,4)
      dz_Rloc(1:lm_max,nRstart:nRstop)  => flow_Rloc_container(:,:,5)
      p_Rloc(1:lm_max,nRstart:nRstop)   => flow_Rloc_container(:,:,6)
      dp_Rloc(1:lm_max,nRstart:nRstop)  => flow_Rloc_container(:,:,7)

      !-- Entropy:
      allocate( s_LMloc_container(llm:ulm,n_r_max,1:2) )
      s_LMloc(llm:ulm,1:n_r_max)  => s_LMloc_container(:,:,1)
      ds_LMloc(llm:ulm,1:n_r_max) => s_LMloc_container(:,:,2)
      allocate( s_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      s_Rloc(1:lm_max,nRstart:nRstop)   => s_Rloc_container(:,:,1)
      ds_Rloc(1:lm_max,nRstart:nRstop)  => s_Rloc_container(:,:,2)

      bytes_allocated = bytes_allocated + &
                        9*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
                        9*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMloc_container(llm:ulm,n_r_max,1:2) )
         xi_LMloc(llm:ulm,1:n_r_max)  => xi_LMloc_container(:,:,1)
         dxi_LMloc(llm:ulm,1:n_r_max) => xi_LMloc_container(:,:,2)
         allocate( xi_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         xi_Rloc(1:lm_max,nRstart:nRstop)   => xi_Rloc_container(:,:,1)
         dxi_Rloc(1:lm_max,nRstart:nRstop)  => xi_Rloc_container(:,:,2)
         bytes_allocated = bytes_allocated + &
                           2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
                           2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_LMloc_container(1,1,2) ) ! For debugging
         xi_LMloc(1:1,1:1)  => xi_LMloc_container(:,:,1)
         dxi_LMloc(1:1,1:1) => xi_LMloc_container(:,:,2)
         allocate( xi_Rloc_container(1,1,2) )
         xi_Rloc(1:1,1:1)   => xi_Rloc_container(:,:,1)
         dxi_Rloc(1:1,1:1)  => xi_Rloc_container(:,:,2)
      end if

      !-- Magnetic field potentials:
      allocate( field_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:6) )
      b_LMloc(llmMag:ulmMag,1:n_r_maxMag)   => field_LMloc_container(:,:,1)
      db_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,2)
      ddb_LMloc(llmMag:ulmMag,1:n_r_maxMag) => field_LMloc_container(:,:,3)
      aj_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,4)
      dj_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,5)
      ddj_LMloc(llmMag:ulmMag,1:n_r_maxMag) => field_LMloc_container(:,:,6)

      allocate( field_Rloc_container(lm_maxMag,nRstart:nRstop,1:5) )
      b_Rloc(1:lm_maxMag,nRstart:nRstop)   => field_Rloc_container(:,:,1)
      db_Rloc(1:lm_maxMag,nRstart:nRstop)  => field_Rloc_container(:,:,2)
      ddb_Rloc(1:lm_maxMag,nRstart:nRstop) => field_Rloc_container(:,:,3)
      aj_Rloc(1:lm_maxMag,nRstart:nRstop)  => field_Rloc_container(:,:,4)
      dj_Rloc(1:lm_maxMag,nRstart:nRstop)  => field_Rloc_container(:,:,5)


      bytes_allocated = bytes_allocated + &
                        6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
                        5*lm_maxMag*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

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
                        6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fields
!----------------------------------------------------------------------------
end module fields
