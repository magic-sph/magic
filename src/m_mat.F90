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

    !--- Logicals that inform whether the respective matrix
    !    has been updated:           
    ALLOCATE( lSmat(0:l_max) )
    ALLOCATE( lZmat(0:l_max) )
    ALLOCATE( lWPmat(0:l_max) )
    ALLOCATE( lBmat(0:l_max) )

  END SUBROUTINE initialize_matrices
END MODULE matrices
