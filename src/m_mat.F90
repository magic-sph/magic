!$Id$
!*****************************************************************
!  Common block containing matricies for internal time step
!*****************************************************************

MODULE matrices
  use truncation
  IMPLICIT NONE

  !-- the matrices, already LU-decomposed:
  REAL(kind=8),ALLOCATABLE :: s0Mat(:,:)      ! for l=m=0  
  REAL(kind=8),ALLOCATABLE :: sMat(:,:,:)
  REAL(kind=8),ALLOCATABLE :: zMat(:,:,:) 
  REAL(kind=8),ALLOCATABLE :: z10Mat(:,:)    ! for l=1,m=0 
  REAL(kind=8),ALLOCATABLE :: wpMat(:,:,:) 
  REAL(kind=8),ALLOCATABLE :: bMat(:,:,:)
  REAL(kind=8),ALLOCATABLE :: jMat(:,:,:)

  !-- respecitive pivoting information:
  INTEGER, ALLOCATABLE :: s0Pivot(:)
  INTEGER, ALLOCATABLE :: sPivot(:,:)
  INTEGER, ALLOCATABLE :: z10Pivot(:)
  INTEGER, ALLOCATABLE :: zPivot(:,:)
  INTEGER, ALLOCATABLE :: wpPivot(:,:)
  INTEGER, ALLOCATABLE :: bPivot(:,:)
  INTEGER, ALLOCATABLE :: jPivot(:,:)

  ! -- respective linesums of the matrices
#ifdef WITH_PRECOND_WP
  REAL(kind=8),ALLOCATABLE :: wpMat_fac(:,:,:)
#endif
#ifdef WITH_PRECOND_Z
  REAL(kind=8),ALLOCATABLE :: zMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S
  REAL(kind=8),ALLOCATABLE :: sMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_BJ
  REAL(kind=8),ALLOCATABLE :: bMat_fac(:,:)
  REAL(kind=8),ALLOCATABLE :: jMat_fac(:,:)
#endif
  !--- Logicals that inform whether the respective matrix
  !    has been updated:           
  LOGICAL :: lZ10mat
  LOGICAL,ALLOCATABLE :: lSmat(:)
  LOGICAL,ALLOCATABLE :: lZmat(:)
  LOGICAL,ALLOCATABLE :: lWPmat(:)
  LOGICAL,ALLOCATABLE :: lBmat(:)

CONTAINS
  SUBROUTINE initialize_matrices

    !-- the matrices, already LU-decomposed:
    ALLOCATE( s0Mat(n_r_max,n_r_max) )      ! for l=m=0  
    ALLOCATE( sMat(n_r_max,n_r_max,l_max) )
    ALLOCATE( zMat(n_r_max,n_r_max,l_max) )
    ALLOCATE( z10Mat(n_r_max,n_r_max) )    ! for l=1,m=0 
    ALLOCATE( wpMat(2*n_r_max,2*n_r_max,l_max) )
    ALLOCATE( bMat(n_r_totMag,n_r_totMag,l_maxMag) )
    ALLOCATE( jMat(n_r_totMag,n_r_totMag,l_maxMag) )

    !-- respecitive pivoting information:
    ALLOCATE( s0Pivot(n_r_max) )
    ALLOCATE( sPivot(n_r_max,l_max) )
    ALLOCATE( z10Pivot(n_r_max) )
    ALLOCATE( zPivot(n_r_max,l_max) )
    ALLOCATE( wpPivot(2*n_r_max,l_max) )
    ALLOCATE( bPivot(n_r_totMag,l_maxMag) )
    ALLOCATE( jPivot(n_r_totMag,l_maxMag) )

#ifdef WITH_PRECOND_WP
    ALLOCATE(wpMat_fac(2*n_r_max,2,l_max))
#endif
#ifdef WITH_PRECOND_Z
    ALLOCATE(zMat_fac(n_r_max,l_max))
#endif
#ifdef WITH_PRECOND_S
    ALLOCATE(sMat_fac(n_r_max,l_max))
#endif
#ifdef WITH_PRECOND_BJ
    ALLOCATE(bMat_fac(n_r_tot,l_max))
    ALLOCATE(jMat_fac(n_r_tot,l_max))
#endif

    !--- Logicals that inform whether the respective matrix
    !    has been updated:           
    ALLOCATE( lSmat(0:l_max) )
    ALLOCATE( lZmat(0:l_max) )
    ALLOCATE( lWPmat(0:l_max) )
    ALLOCATE( lBmat(0:l_max) )

  END SUBROUTINE initialize_matrices
END MODULE matrices
