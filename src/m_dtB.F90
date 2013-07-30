!$Id$
!---------------------------------------------------------------------------------
!  This module contains magnetic field stretching and advection terms 
!  plus a separate omega-effect.
!  It is used for movie output.
!--------------------------------------------------------------------------------

MODULE dtB_mod
  use truncation
  implicit none
  COMPLEX(kind=8),allocatable :: PstrLM(:,:)
  COMPLEX(kind=8),allocatable :: PadvLM(:,:)
  COMPLEX(kind=8),allocatable :: TstrLM(:,:)
  COMPLEX(kind=8),allocatable :: TadvLM(:,:)
  COMPLEX(kind=8),allocatable :: PdifLM(:,:)
  COMPLEX(kind=8),allocatable :: TdifLM(:,:)
  COMPLEX(kind=8),allocatable :: PadvLMIC(:,:)
  COMPLEX(kind=8),allocatable :: PdifLMIC(:,:)
  COMPLEX(kind=8),allocatable :: TadvLMIC(:,:)
  COMPLEX(kind=8),allocatable :: TdifLMIC(:,:)
  COMPLEX(kind=8),allocatable :: TomeLM(:,:)

  REAL(kind=8) :: PstrRms,PstrAsRms
  REAL(kind=8) :: PadvRms,PadvAsRms
  REAL(kind=8) :: PdifRms,PdifAsRms
  REAL(kind=8) :: TstrRms,TstrAsRms
  REAL(kind=8) :: TadvRms,TadvAsRms
  REAL(kind=8) :: TdifRms,TdifAsRms
  REAL(kind=8) :: TomeRms,TomeAsRms

CONTAINS
  SUBROUTINE initialize_dtB_mod

  ALLOCATE( PstrLM(lm_max_dtB,n_r_max_dtB) )
  ALLOCATE( PadvLM(lm_max_dtB,n_r_max_dtB) )
  ALLOCATE( TstrLM(lm_max_dtB,n_r_max_dtB) )
  ALLOCATE( TadvLM(lm_max_dtB,n_r_max_dtB) )
  ALLOCATE( PdifLM(lm_max_dtB,n_r_max_dtB) )
  ALLOCATE( TdifLM(lm_max_dtB,n_r_max_dtB) )
  ALLOCATE( PadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
  ALLOCATE( PdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
  ALLOCATE( TadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
  ALLOCATE( TdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
  ALLOCATE( TomeLM(lm_max_dtB,n_r_max_dtB) )

END SUBROUTINE initialize_dtB_mod
END MODULE dtB_mod
