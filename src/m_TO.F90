!$Id$
!----------------------------------------------------------------------
!  This module contains information for TO calculation and output
!----------------------------------------------------------------------

MODULE torsional_oscillations
  use truncation

  IMPLICIT NONE

  REAL(kind=8),ALLOCATABLE :: dzStrLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: dzRstrLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: dzAstrLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: dzCorLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: dzLFLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: dzdVpLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: dzddVpLMr(:,:)
  REAL(kind=8),ALLOCATABLE :: V2AS(:,:)
  REAL(kind=8),ALLOCATABLE :: Bs2AS(:,:)
  REAL(kind=8),ALLOCATABLE :: BszAS(:,:)
  REAL(kind=8),ALLOCATABLE :: BspAS(:,:)
  REAL(kind=8),ALLOCATABLE :: BpzAS(:,:)
  REAL(kind=8),ALLOCATABLE :: BspdAS(:,:)
  REAL(kind=8),ALLOCATABLE :: BpsdAS(:,:)
  REAL(kind=8),ALLOCATABLE :: BzpdAS(:,:)
  REAL(kind=8),ALLOCATABLE :: BpzdAS(:,:)
  REAL(kind=8),ALLOCATABLE :: ddzASL(:,:)

contains
  SUBROUTINE initialize_TO
    ALLOCATE( dzStrLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( dzRstrLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( dzAstrLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( dzCorLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( dzLFLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( dzdVpLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( dzddVpLMr(l_max+1,n_r_maxStr) )
    ALLOCATE( V2AS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( Bs2AS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BszAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BspAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BpzAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BspdAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BpsdAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BzpdAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( BpzdAS(n_theta_maxStr,n_r_maxStr) )
    ALLOCATE( ddzASL(l_max+1,n_r_maxStr) ) 

  END SUBROUTINE initialize_TO
END MODULE torsional_oscillations
