!$Id$
!***************************************************************
!  Common blocks containing the potential fields and their radial
!  derivatives
!***************************************************************
#include "intrinsic_sizes.h"
MODULE fields
  use truncation
  USE LMLoop_data,ONLY: llm,ulm,llmMag,ulmMag
  USE radial_data,ONLY: nRstart,nRstop
  USE parallel_mod,only: rank
  !USE cutils
  implicit none

  !-- Velocity potentials:
  COMPLEX(kind=8),ALLOCATABLE :: w(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dw(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddw(:,:)
  COMPLEX(kind=8),ALLOCATABLE,TARGET :: w_LMloc_container(:,:,:),w_Rloc_container(:,:,:)
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: w_LMloc,dw_LMloc,ddw_LMloc
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: w_Rloc, dw_Rloc, ddw_Rloc

  COMPLEX(kind=8),ALLOCATABLE :: z(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dz(:,:)
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:),TARGET :: z_LMloc_container,z_Rloc_container
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: z_LMloc,dz_LMloc,z_Rloc,dz_Rloc

  !-- Pressure and entropy:
  COMPLEX(kind=8),ALLOCATABLE :: s(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ds(:,:)
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:),TARGET :: s_LMloc_container,s_Rloc_container
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: s_LMloc,ds_LMloc,s_Rloc,ds_Rloc

  COMPLEX(kind=8),ALLOCATABLE :: p(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dp(:,:)
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:),TARGET :: p_LMloc_container,p_Rloc_container
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: p_LMloc,dp_LMloc,p_Rloc,dp_Rloc

  !-- Magnetic field potentials:
  COMPLEX(kind=8),ALLOCATABLE :: b(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: db(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj(:,:)
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:),TARGET :: b_LMloc_container,b_Rloc_container
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: b_LMloc,db_LMloc,ddb_LMloc,b_Rloc,db_Rloc,ddb_Rloc
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:),TARGET :: aj_LMloc_container,aj_Rloc_container
  COMPLEX(kind=8),DIMENSION(:,:),POINTER :: aj_LMloc,dj_LMloc,ddj_LMloc,aj_Rloc,dj_Rloc

  !-- Magnetic field potentials in inner core:
  !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
  !         for even chebs
  COMPLEX(kind=8),ALLOCATABLE :: b_ic(:,:)  
  COMPLEX(kind=8),ALLOCATABLE :: db_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_ic(:,:) 
  COMPLEX(kind=8),ALLOCATABLE :: dj_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: b_ic_LMloc(:,:)  
  COMPLEX(kind=8),ALLOCATABLE :: db_ic_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb_ic_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_ic_LMloc(:,:) 
  COMPLEX(kind=8),ALLOCATABLE :: dj_ic_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj_ic_LMloc(:,:)

  !-- Rotation rates:
  REAL(kind=8) :: omega_ic,omega_ma

CONTAINS
  SUBROUTINE initialize_fields
    INTEGER(kind=8) :: bytes_allocated

    bytes_allocated = 0
    !-- Velocity potentials:
    IF (rank.EQ.0) THEN
       ALLOCATE( w(lm_max,n_r_max) )
       ALLOCATE( z(lm_max,n_r_max) )
       ALLOCATE( dw(lm_max,n_r_max) )
       ALLOCATE( ddw(lm_max,n_r_max) )
       ALLOCATE( dz(lm_max,n_r_max) )
       ALLOCATE( s(lm_max,n_r_max) )
       ALLOCATE( ds(lm_max,n_r_max) )
       ALLOCATE( p(lm_max,n_r_max) )
       ALLOCATE( dp(lm_max,n_r_max) )
       bytes_allocated = bytes_allocated + 9*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX
       ALLOCATE( b(lm_maxMag,n_r_maxMag) )
       ALLOCATE( db(lm_maxMag,n_r_maxMag) )
       ALLOCATE( ddb(lm_maxMag,n_r_maxMag) )
       ALLOCATE( aj(lm_maxMag,n_r_maxMag) )
       ALLOCATE( dj(lm_maxMag,n_r_maxMag) )
       ALLOCATE( ddj(lm_maxMag,n_r_maxMag) )
       bytes_allocated = bytes_allocated + 6*lm_maxMag*n_r_maxMag*SIZEOF_DOUBLE_COMPLEX
       ALLOCATE( b_ic(lm_maxMag,n_r_ic_maxMag) )  
       ALLOCATE( db_ic(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
       ALLOCATE( dj_ic(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
       bytes_allocated = bytes_allocated + 6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DOUBLE_COMPLEX

    ELSE
       ALLOCATE( w(1,n_r_max) )
       ALLOCATE( z(1,n_r_max) )
       ALLOCATE( dw(1,n_r_max) )
       ALLOCATE( ddw(1,n_r_max) )
       ALLOCATE( dz(1,n_r_max) )
       ALLOCATE( s(1,n_r_max) )
       ALLOCATE( ds(1,n_r_max) )
       ALLOCATE( p(1,n_r_max) )
       ALLOCATE( dp(1,n_r_max) )
       bytes_allocated = bytes_allocated + 9*n_r_max*SIZEOF_DOUBLE_COMPLEX
       ALLOCATE( b(1,n_r_maxMag) )
       ALLOCATE( db(1,n_r_maxMag) )
       ALLOCATE( ddb(1,n_r_maxMag) )
       ALLOCATE( aj(1,n_r_maxMag) )
       ALLOCATE( dj(1,n_r_maxMag) )
       ALLOCATE( ddj(1,n_r_maxMag) )
       bytes_allocated = bytes_allocated + 6*n_r_maxMag*SIZEOF_DOUBLE_COMPLEX
       ALLOCATE( b_ic(1,n_r_ic_maxMag) )  
       ALLOCATE( db_ic(1,n_r_ic_maxMag) )
       ALLOCATE( ddb_ic(1,n_r_ic_maxMag) )
       ALLOCATE( aj_ic(1,n_r_ic_maxMag) ) 
       ALLOCATE( dj_ic(1,n_r_ic_maxMag) )
       ALLOCATE( ddj_ic(1,n_r_ic_maxMag) )
       bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DOUBLE_COMPLEX
    END IF
    ALLOCATE( w_LMloc_container(llm:ulm,n_r_max,1:3) )
    w_LMloc(llm:,1:)   => w_LMloc_container(llm:ulm,1:n_r_max,1)
    dw_LMloc(llm:,1:)  => w_LMloc_container(llm:ulm,1:n_r_max,2)
    ddw_LMloc(llm:,1:) => w_LMloc_container(llm:ulm,1:n_r_max,3)
    ALLOCATE( w_Rloc_container(lm_max,nRstart:nRstop,1:3) )
    w_Rloc(1:,nRstart:)   => w_Rloc_container(1:,nRstart:,1)
    dw_Rloc(1:,nRstart:)  => w_Rloc_container(1:,nRstart:,2)
    ddw_Rloc(1:,nRstart:) => w_Rloc_container(1:,nRstart:,3)

    ALLOCATE( z_LMloc_container(llm:ulm,n_r_max,1:2) )
    z_LMloc(llm:,1:)   => z_LMloc_container(llm:ulm,1:n_r_max,1)
    dz_LMloc(llm:,1:)  => z_LMloc_container(llm:ulm,1:n_r_max,2)
    ALLOCATE( z_Rloc_container(lm_max,nRstart:nRstop,1:2) )
    z_Rloc(1:,nRstart:)   => z_Rloc_container(1:,nRstart:,1)
    dz_Rloc(1:,nRstart:)  => z_Rloc_container(1:,nRstart:,2)

    !-- Pressure and entropy:
    ALLOCATE( s_LMloc_container(llm:ulm,n_r_max,1:2) )
    s_LMloc(llm:,1:)   => s_LMloc_container(llm:ulm,1:n_r_max,1)
    ds_LMloc(llm:,1:)  => s_LMloc_container(llm:ulm,1:n_r_max,2)
    ALLOCATE( s_Rloc_container(lm_max,nRstart:nRstop,1:2) )
    s_Rloc(1:,nRstart:)   => s_Rloc_container(1:,nRstart:,1)
    ds_Rloc(1:,nRstart:)  => s_Rloc_container(1:,nRstart:,2)

    ALLOCATE( p_LMloc_container(llm:ulm,n_r_max,1:2) )
    p_LMloc(llm:,1:)   => p_LMloc_container(llm:ulm,1:n_r_max,1)
    dp_LMloc(llm:,1:)  => p_LMloc_container(llm:ulm,1:n_r_max,2)
    ALLOCATE( p_Rloc_container(lm_max,nRstart:nRstop,1:2) )
    p_Rloc(1:,nRstart:)   => p_Rloc_container(1:,nRstart:,1)
    dp_Rloc(1:,nRstart:)  => p_Rloc_container(1:,nRstart:,2)

    bytes_allocated = bytes_allocated + 9*(ulm-llm+1)*n_r_max*SIZEOF_DOUBLE_COMPLEX
    bytes_allocated = bytes_allocated + 9*lm_max*(nRstop-nRstart+1)*SIZEOF_DOUBLE_COMPLEX

    !-- Magnetic field potentials:
    ALLOCATE( b_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3) )
    b_LMloc(llmMag:,1:)   => b_LMloc_container(llmMag:,1:,1)
    db_LMloc(llmMag:,1:)  => b_LMloc_container(llmMag:,1:,2)
    ddb_LMloc(llmMag:,1:) => b_LMloc_container(llmMag:,1:,3)
    ALLOCATE( b_Rloc_container(lm_maxMag,nRstart:nRstop,1:3) )
    b_Rloc(1:,nRstart:)   => b_Rloc_container(1:,nRstart:,1)
    db_Rloc(1:,nRstart:)  => b_Rloc_container(1:,nRstart:,2)
    ddb_Rloc(1:,nRstart:) => b_Rloc_container(1:,nRstart:,3)

    ALLOCATE( aj_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3) )
    aj_LMloc(llmMag:,1:)  => aj_LMloc_container(llmMag:,1:,1)
    dj_LMloc(llmMag:,1:)  => aj_LMloc_container(llmMag:,1:,2)
    ddj_LMloc(llmMag:,1:) => aj_LMloc_container(llmMag:,1:,3)
    ALLOCATE( aj_Rloc_container(lm_maxMag,nRstart:nRstop,1:2) )
    aj_Rloc(1:,nRstart:)  => aj_Rloc_container(1:,nRstart:,1)
    dj_Rloc(1:,nRstart:)  => aj_Rloc_container(1:,nRstart:,2)

    bytes_allocated = bytes_allocated + 6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DOUBLE_COMPLEX
    bytes_allocated = bytes_allocated + 5*lm_maxMag*(nRstop-nRstart+1)*SIZEOF_DOUBLE_COMPLEX

    !-- Magnetic field potentials in inner core:
    !   NOTE: n_r-dimension may be smaller once CHEBFT is adopted
    !         for even chebs
    ALLOCATE( b_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )  
    ALLOCATE( db_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( ddb_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( aj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) ) 
    ALLOCATE( dj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( ddj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )

    bytes_allocated = bytes_allocated + 6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DOUBLE_COMPLEX

    WRITE(*,"(I4,A,I12,A)") rank,": Allocated in m_fields ",bytes_allocated," bytes."

  END SUBROUTINE initialize_fields
END MODULE fields
