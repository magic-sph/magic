!$Id$
!***************************************************************
!  Common blocks containing the potential fields and their radial
!  derivatives
!***************************************************************

MODULE fields
  use truncation
  USE LMLoop_data,ONLY: llm,ulm,llmMag,ulmMag
  USE radial_data,ONLY: nRstart,nRstop
  USE parallel_mod,only: rank
  implicit none

  !-- Velocity potentials:
  COMPLEX(kind=8),ALLOCATABLE :: w(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dw(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddw(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: w_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dw_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddw_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: w_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dw_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddw_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: z(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dz(:,:)
  !COMPLEX(kind=8),ALLOCATABLE :: z_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: z_LMloc(:,:)
  !COMPLEX(kind=8),ALLOCATABLE :: dz_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dz_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: z_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dz_Rloc(:,:)

  !-- Pressure and entropy:
  COMPLEX(kind=8),ALLOCATABLE :: s(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ds(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: s_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ds_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: s_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ds_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: p(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dp(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: p_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dp_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: p_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dp_Rloc(:,:)

  !-- Magnetic field potentials:
  COMPLEX(kind=8),ALLOCATABLE :: b(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: db(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: b_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: db_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: b_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: db_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dj_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_Rloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dj_Rloc(:,:)

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
       ALLOCATE( b(lm_maxMag,n_r_maxMag) )
       ALLOCATE( db(lm_maxMag,n_r_maxMag) )
       ALLOCATE( ddb(lm_maxMag,n_r_maxMag) )
       ALLOCATE( aj(lm_maxMag,n_r_maxMag) )
       ALLOCATE( dj(lm_maxMag,n_r_maxMag) )
       ALLOCATE( ddj(lm_maxMag,n_r_maxMag) )
       ALLOCATE( b_ic(lm_maxMag,n_r_ic_maxMag) )  
       ALLOCATE( db_ic(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
       ALLOCATE( dj_ic(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
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
       ALLOCATE( b(1,n_r_maxMag) )
       ALLOCATE( db(1,n_r_maxMag) )
       ALLOCATE( ddb(1,n_r_maxMag) )
       ALLOCATE( aj(1,n_r_maxMag) )
       ALLOCATE( dj(1,n_r_maxMag) )
       ALLOCATE( ddj(1,n_r_maxMag) )
       ALLOCATE( b_ic(1,n_r_ic_maxMag) )  
       ALLOCATE( db_ic(1,n_r_ic_maxMag) )
       ALLOCATE( ddb_ic(1,n_r_ic_maxMag) )
       ALLOCATE( aj_ic(1,n_r_ic_maxMag) ) 
       ALLOCATE( dj_ic(1,n_r_ic_maxMag) )
       ALLOCATE( ddj_ic(1,n_r_ic_maxMag) )
    END IF
    ALLOCATE( w_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( dw_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( ddw_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( w_Rloc(lm_max,nRstart:nRstop) )
    ALLOCATE( dw_Rloc(lm_max,nRstart:nRstop) )
    ALLOCATE( ddw_Rloc(lm_max,nRstart:nRstop) )
    !ALLOCATE( z_LMloc(llm:ulm,n_r_max) )
    !ALLOCATE( dz_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( z_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( dz_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( z_Rloc(lm_max,nRstart:nRstop) )
    ALLOCATE( dz_Rloc(lm_max,nRstart:nRstop) )

    !-- Pressure and entropy:
    ALLOCATE( s_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( s_Rloc(lm_max,nRstart:nRstop) )
    ALLOCATE( ds_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( ds_Rloc(lm_max,nRstart:nRstop) )
    ALLOCATE( p_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( dp_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( p_Rloc(lm_max,nRstart:nRstop) )
    ALLOCATE( dp_Rloc(lm_max,nRstart:nRstop) )

    !-- Magnetic field potentials:
    ALLOCATE( b_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( db_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( ddb_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( aj_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( dj_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( ddj_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( b_Rloc(lm_maxMag,nRstart:nRstop) )
    ALLOCATE( db_Rloc(lm_maxMag,nRstart:nRstop) )
    ALLOCATE( ddb_Rloc(lm_maxMag,nRstart:nRstop) )
    ALLOCATE( aj_Rloc(lm_maxMag,nRstart:nRstop) )
    ALLOCATE( dj_Rloc(lm_maxMag,nRstart:nRstop) )

    !-- Magnetic field potentials in inner core:
    !   NOTE: n_r-dimension may be smaller once CHEBFT is adopted
    !         for even chebs
    ALLOCATE( b_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )  
    ALLOCATE( db_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( ddb_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( aj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) ) 
    ALLOCATE( dj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( ddj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
  END SUBROUTINE initialize_fields
END MODULE fields
