!$Id$
!***************************************************************
!  Common blocks containing the potential fields and their radial
!  derivatives
!***************************************************************

MODULE fields
  use truncation
  implicit none

  !-- Velocity potentials:
  COMPLEX(kind=8),ALLOCATABLE :: w(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dw(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddw(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: z(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dz(:,:)
  !COMMON/v_field/w,dw,ddw,z,dz

  !-- Pressure and entropy:
  COMPLEX(kind=8),ALLOCATABLE :: s(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ds(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: p(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dp(:,:)
  !COMMON/ps_field/s,ds,p,dp

  !-- Magnetic field potentials:
  COMPLEX(kind=8),ALLOCATABLE :: b(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: db(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dj(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj(:,:)
  !COMMON/b_field/b,db,ddb,aj,dj,ddj

  !-- Magnetic field potentials in inner core:
  !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
  !         for even chebs
  COMPLEX(kind=8),ALLOCATABLE :: b_ic(:,:)  
  COMPLEX(kind=8),ALLOCATABLE :: db_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddb_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_ic(:,:) 
  COMPLEX(kind=8),ALLOCATABLE :: dj_ic(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: ddj_ic(:,:)
  !COMMON/b_field_ic/b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic

  !-- Rotation rates:
  REAL(kind=8) :: omega_ic,omega_ma
  !COMMON/rotation/omega_ic,omega_ma

CONTAINS
  SUBROUTINE initialize_fields
    !-- Velocity potentials:
    ALLOCATE( w(lm_max,n_r_max) )
    ALLOCATE( dw(lm_max,n_r_max) )
    ALLOCATE( ddw(lm_max,n_r_max) )
    ALLOCATE( z(lm_max,n_r_max) )
    ALLOCATE( dz(lm_max,n_r_max) )

    !-- Pressure and entropy:
    ALLOCATE( s(lm_max,n_r_max) )
    ALLOCATE( ds(lm_max,n_r_max) )
    ALLOCATE( p(lm_max,n_r_max) )
    ALLOCATE( dp(lm_max,n_r_max) )

    !-- Magnetic field potentials:
    ALLOCATE( b(lm_maxMag,n_r_maxMag) )
    ALLOCATE( db(lm_maxMag,n_r_maxMag) )
    ALLOCATE( ddb(lm_maxMag,n_r_maxMag) )
    ALLOCATE( aj(lm_maxMag,n_r_maxMag) )
    ALLOCATE( dj(lm_maxMag,n_r_maxMag) )
    ALLOCATE( ddj(lm_maxMag,n_r_maxMag) )

    !-- Magnetic field potentials in inner core:
    !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
    !         for even chebs
    ALLOCATE( b_ic(lm_maxMag,n_r_ic_maxMag) )  
    ALLOCATE( db_ic(lm_maxMag,n_r_ic_maxMag) )
    ALLOCATE( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
    ALLOCATE( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
    ALLOCATE( dj_ic(lm_maxMag,n_r_ic_maxMag) )
    ALLOCATE( ddj_ic(lm_maxMag,n_r_ic_maxMag) )

  END SUBROUTINE initialize_fields
END MODULE fields
