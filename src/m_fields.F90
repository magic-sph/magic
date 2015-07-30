!$Id$
#include "intrinsic_sizes.h"
module fields
!***************************************************************
!  Common blocks containing the potential fields and their radial
!  derivatives
!***************************************************************
   use truncation
   use LMLoop_data, only: llm,ulm,llmMag,ulmMag
   use radial_data, only: nRstart,nRstop
   use parallel_mod, only: rank
 
   implicit none
 
   !-- Velocity potentials:
   complex(kind=8), allocatable :: w(:,:)
   complex(kind=8), allocatable :: dw(:,:)
   complex(kind=8), allocatable :: ddw(:,:)
   complex(kind=8), allocatable, target :: w_LMloc_container(:,:,:)
   complex(kind=8), allocatable, target :: w_Rloc_container(:,:,:)
   complex(kind=8), pointer :: w_LMloc(:,:),dw_LMloc(:,:),ddw_LMloc(:,:)
   complex(kind=8), pointer :: w_Rloc(:,:), dw_Rloc(:,:), ddw_Rloc(:,:)
 
   complex(kind=8), allocatable :: z(:,:)
   complex(kind=8), allocatable :: dz(:,:)
   complex(kind=8), allocatable, target :: z_LMloc_container(:,:,:)
   complex(kind=8), allocatable, target :: z_Rloc_container(:,:,:)
   complex(kind=8), pointer :: z_LMloc(:,:),dz_LMloc(:,:), z_Rloc(:,:), dz_Rloc(:,:)
 
   !-- Pressure and entropy:
   complex(kind=8), allocatable :: s(:,:)
   complex(kind=8), allocatable :: ds(:,:)
   complex(kind=8), allocatable, target :: s_LMloc_container(:,:,:)
   complex(kind=8), allocatable, target :: s_Rloc_container(:,:,:)
   complex(kind=8), pointer :: s_LMloc(:,:), ds_LMloc(:,:), s_Rloc(:,:), ds_Rloc(:,:)
 
   complex(kind=8), allocatable :: p(:,:)
   complex(kind=8), allocatable :: dp(:,:)
   complex(kind=8), allocatable, target :: p_LMloc_container(:,:,:)
   complex(kind=8), allocatable, target :: p_Rloc_container(:,:,:)
   complex(kind=8), pointer :: p_LMloc(:,:), dp_LMloc(:,:), p_Rloc(:,:), dp_Rloc(:,:)
 
   !-- Magnetic field potentials:
   complex(kind=8), allocatable :: b(:,:)
   complex(kind=8), allocatable :: db(:,:)
   complex(kind=8), allocatable :: ddb(:,:)
   complex(kind=8), allocatable :: aj(:,:)
   complex(kind=8), allocatable :: dj(:,:)
   complex(kind=8), allocatable :: ddj(:,:)
   complex(kind=8), allocatable, target :: b_LMloc_container(:,:,:)
   complex(kind=8), allocatable, target :: b_Rloc_container(:,:,:)
   complex(kind=8), pointer :: b_LMloc(:,:), db_LMloc(:,:), ddb_LMloc(:,:)
   complex(kind=8), pointer :: b_Rloc(:,:), db_Rloc(:,:), ddb_Rloc(:,:)
   complex(kind=8), allocatable, target :: aj_LMloc_container(:,:,:)
   complex(kind=8), allocatable, target :: aj_Rloc_container(:,:,:)
   complex(kind=8), pointer :: aj_LMloc(:,:), dj_LMloc(:,:), ddj_LMloc(:,:)
   complex(kind=8), pointer :: aj_Rloc(:,:), dj_Rloc(:,:)
 
   !-- Magnetic field potentials in inner core:
   !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
   !         for even chebs
   complex(kind=8), allocatable :: b_ic(:,:)  
   complex(kind=8), allocatable :: db_ic(:,:)
   complex(kind=8), allocatable :: ddb_ic(:,:)
   complex(kind=8), allocatable :: aj_ic(:,:) 
   complex(kind=8), allocatable :: dj_ic(:,:)
   complex(kind=8), allocatable :: ddj_ic(:,:)
   complex(kind=8), allocatable :: b_ic_LMloc(:,:)  
   complex(kind=8), allocatable :: db_ic_LMloc(:,:)
   complex(kind=8), allocatable :: ddb_ic_LMloc(:,:)
   complex(kind=8), allocatable :: aj_ic_LMloc(:,:) 
   complex(kind=8), allocatable :: dj_ic_LMloc(:,:)
   complex(kind=8), allocatable :: ddj_ic_LMloc(:,:)
 
   !-- Rotation rates:
   real(kind=8) :: omega_ic,omega_ma

contains

   subroutine initialize_fields

      integer(kind=8) :: bytes_allocated

      bytes_allocated = 0
      !-- Velocity potentials:
      if ( rank == 0 ) then
         allocate( w(lm_max,n_r_max) )
         allocate( z(lm_max,n_r_max) )
         allocate( dw(lm_max,n_r_max) )
         allocate( ddw(lm_max,n_r_max) )
         allocate( dz(lm_max,n_r_max) )
         allocate( s(lm_max,n_r_max) )
         allocate( ds(lm_max,n_r_max) )
         allocate( p(lm_max,n_r_max) )
         allocate( dp(lm_max,n_r_max) )
         bytes_allocated = bytes_allocated + 9*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX
         allocate( b(lm_maxMag,n_r_maxMag) )
         allocate( db(lm_maxMag,n_r_maxMag) )
         allocate( ddb(lm_maxMag,n_r_maxMag) )
         allocate( aj(lm_maxMag,n_r_maxMag) )
         allocate( dj(lm_maxMag,n_r_maxMag) )
         allocate( ddj(lm_maxMag,n_r_maxMag) )
         bytes_allocated = bytes_allocated +  &
                           6*lm_maxMag*n_r_maxMag*SIZEOF_DOUBLE_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )  
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
                           6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DOUBLE_COMPLEX

      else
         allocate( w(1,n_r_max) )
         allocate( z(1,n_r_max) )
         allocate( dw(1,n_r_max) )
         allocate( ddw(1,n_r_max) )
         allocate( dz(1,n_r_max) )
         allocate( s(1,n_r_max) )
         allocate( ds(1,n_r_max) )
         allocate( p(1,n_r_max) )
         allocate( dp(1,n_r_max) )
         bytes_allocated = bytes_allocated + 9*n_r_max*SIZEOF_DOUBLE_COMPLEX
         allocate( b(1,n_r_maxMag) )
         allocate( db(1,n_r_maxMag) )
         allocate( ddb(1,n_r_maxMag) )
         allocate( aj(1,n_r_maxMag) )
         allocate( dj(1,n_r_maxMag) )
         allocate( ddj(1,n_r_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_maxMag*SIZEOF_DOUBLE_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )  
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) ) 
         allocate( dj_ic(1,n_r_ic_maxMag) )
         allocate( ddj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DOUBLE_COMPLEX
      end if
      allocate( w_LMloc_container(llm:ulm,n_r_max,1:3) )
      w_LMloc(llm:,1:)   => w_LMloc_container(llm:ulm,1:n_r_max,1)
      dw_LMloc(llm:,1:)  => w_LMloc_container(llm:ulm,1:n_r_max,2)
      ddw_LMloc(llm:,1:) => w_LMloc_container(llm:ulm,1:n_r_max,3)
      allocate( w_Rloc_container(lm_max,nRstart:nRstop,1:3) )
      w_Rloc(1:,nRstart:)   => w_Rloc_container(1:,nRstart:,1)
      dw_Rloc(1:,nRstart:)  => w_Rloc_container(1:,nRstart:,2)
      ddw_Rloc(1:,nRstart:) => w_Rloc_container(1:,nRstart:,3)

      allocate( z_LMloc_container(llm:ulm,n_r_max,1:2) )
      z_LMloc(llm:,1:)   => z_LMloc_container(llm:ulm,1:n_r_max,1)
      dz_LMloc(llm:,1:)  => z_LMloc_container(llm:ulm,1:n_r_max,2)
      allocate( z_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      z_Rloc(1:,nRstart:)   => z_Rloc_container(1:,nRstart:,1)
      dz_Rloc(1:,nRstart:)  => z_Rloc_container(1:,nRstart:,2)

      !-- Pressure and entropy:
      allocate( s_LMloc_container(llm:ulm,n_r_max,1:2) )
      s_LMloc(llm:,1:)   => s_LMloc_container(llm:ulm,1:n_r_max,1)
      ds_LMloc(llm:,1:)  => s_LMloc_container(llm:ulm,1:n_r_max,2)
      allocate( s_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      s_Rloc(1:,nRstart:)   => s_Rloc_container(1:,nRstart:,1)
      ds_Rloc(1:,nRstart:)  => s_Rloc_container(1:,nRstart:,2)

      allocate( p_LMloc_container(llm:ulm,n_r_max,1:2) )
      p_LMloc(llm:,1:)   => p_LMloc_container(llm:ulm,1:n_r_max,1)
      dp_LMloc(llm:,1:)  => p_LMloc_container(llm:ulm,1:n_r_max,2)
      allocate( p_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      p_Rloc(1:,nRstart:)   => p_Rloc_container(1:,nRstart:,1)
      dp_Rloc(1:,nRstart:)  => p_Rloc_container(1:,nRstart:,2)

      bytes_allocated = bytes_allocated + &
                        9*(ulm-llm+1)*n_r_max*SIZEOF_DOUBLE_COMPLEX
      bytes_allocated = bytes_allocated + &
                        9*lm_max*(nRstop-nRstart+1)*SIZEOF_DOUBLE_COMPLEX

      !-- Magnetic field potentials:
      allocate( b_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3) )
      b_LMloc(llmMag:,1:)   => b_LMloc_container(llmMag:,1:,1)
      db_LMloc(llmMag:,1:)  => b_LMloc_container(llmMag:,1:,2)
      ddb_LMloc(llmMag:,1:) => b_LMloc_container(llmMag:,1:,3)
      allocate( b_Rloc_container(lm_maxMag,nRstart:nRstop,1:3) )
      b_Rloc(1:,nRstart:)   => b_Rloc_container(1:,nRstart:,1)
      db_Rloc(1:,nRstart:)  => b_Rloc_container(1:,nRstart:,2)
      ddb_Rloc(1:,nRstart:) => b_Rloc_container(1:,nRstart:,3)

      allocate( aj_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3) )
      aj_LMloc(llmMag:,1:)  => aj_LMloc_container(llmMag:,1:,1)
      dj_LMloc(llmMag:,1:)  => aj_LMloc_container(llmMag:,1:,2)
      ddj_LMloc(llmMag:,1:) => aj_LMloc_container(llmMag:,1:,3)
      allocate( aj_Rloc_container(lm_maxMag,nRstart:nRstop,1:2) )
      aj_Rloc(1:,nRstart:)  => aj_Rloc_container(1:,nRstart:,1)
      dj_Rloc(1:,nRstart:)  => aj_Rloc_container(1:,nRstart:,2)

      bytes_allocated = bytes_allocated + &
                        6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DOUBLE_COMPLEX
      bytes_allocated = bytes_allocated + &
                        5*lm_maxMag*(nRstop-nRstart+1)*SIZEOF_DOUBLE_COMPLEX

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
                        6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DOUBLE_COMPLEX

      write(*,"(I4,A,I12,A)") rank,": Allocated in m_fields ", &
                              bytes_allocated," bytes."

   end subroutine initialize_fields
!----------------------------------------------------------------------------
end module fields
