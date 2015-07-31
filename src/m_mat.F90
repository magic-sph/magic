!$Id$
module matrices
   !--------------------------------------------------------------
   !  Common block containing matricies for internal time step
   !--------------------------------------------------------------

   implicit none
 
   !-- the matrices, already LU-decomposed:
   real(kind=8),allocatable :: s0Mat(:,:)      ! for l=m=0  
   real(kind=8),allocatable :: sMat(:,:,:)
   real(kind=8),allocatable :: zMat(:,:,:) 
   real(kind=8),allocatable :: z10Mat(:,:)    ! for l=1,m=0 
   real(kind=8),allocatable :: wpMat(:,:,:) 
   real(kind=8),allocatable :: bMat(:,:,:)
   real(kind=8),allocatable :: jMat(:,:,:)
 
   !-- respecitive pivoting information:
   integer, allocatable :: s0Pivot(:)
   integer, allocatable :: sPivot(:,:)
   integer, allocatable :: z10Pivot(:)
   integer, allocatable :: zPivot(:,:)
   integer, allocatable :: wpPivot(:,:)
   integer, allocatable :: bPivot(:,:)
   integer, allocatable :: jPivot(:,:)
 
   ! -- respective linesums of the matrices
   real(kind=8), allocatable :: wpMat_fac(:,:,:)
#ifdef WITH_PRECOND_Z
   real(kind=8), allocatable :: zMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_Z10
   real(kind=8), allocatable :: z10Mat_fac(:)
#endif
#ifdef WITH_PRECOND_S
   real(kind=8), allocatable :: sMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(kind=8), allocatable :: s0Mat_fac(:)
#endif
#ifdef WITH_PRECOND_BJ
   real(kind=8), allocatable :: bMat_fac(:,:)
   real(kind=8), allocatable :: jMat_fac(:,:)
#endif
   !--- Logicals that inform whether the respective matrix
   !    has been updated:           
   logical :: lZ10mat
   logical,allocatable :: lSmat(:)
   logical,allocatable :: lZmat(:)
   logical,allocatable :: lWPmat(:)
   logical,allocatable :: lBmat(:)

contains

   subroutine initialize_matrices

      use truncation, only: n_r_max, l_max, l_maxMag, n_r_totMag, &
                            n_r_tot

      !-- the matrices, already LU-decomposed:
      allocate( s0Mat(n_r_max,n_r_max) )      ! for l=m=0  
      allocate( sMat(n_r_max,n_r_max,l_max) )
      allocate( zMat(n_r_max,n_r_max,l_max) )
      allocate( z10Mat(n_r_max,n_r_max) )    ! for l=1,m=0 
      allocate( wpMat(2*n_r_max,2*n_r_max,l_max) )
      allocate( bMat(n_r_totMag,n_r_totMag,l_maxMag) )
      allocate( jMat(n_r_totMag,n_r_totMag,l_maxMag) )

      !-- respecitive pivoting information:
      allocate( s0Pivot(n_r_max) )
      allocate( sPivot(n_r_max,l_max) )
      allocate( z10Pivot(n_r_max) )
      allocate( zPivot(n_r_max,l_max) )
      allocate( wpPivot(2*n_r_max,l_max) )
      allocate( bPivot(n_r_totMag,l_maxMag) )
      allocate( jPivot(n_r_totMag,l_maxMag) )

      allocate(wpMat_fac(2*n_r_max,2,l_max))
#ifdef WITH_PRECOND_Z10
      allocate(z10Mat_fac(n_r_max))
#endif
#ifdef WITH_PRECOND_Z
      allocate(zMat_fac(n_r_max,l_max))
#endif
#ifdef WITH_PRECOND_S
      allocate(sMat_fac(n_r_max,l_max))
#endif
#ifdef WITH_PRECOND_S0
      allocate(s0Mat_fac(n_r_max))
#endif
#ifdef WITH_PRECOND_BJ
      allocate(bMat_fac(n_r_tot,l_max))
      allocate(jMat_fac(n_r_tot,l_max))
#endif

      !--- Logicals that inform whether the respective matrix
      !    has been updated:           
      allocate( lSmat(0:l_max) )
      allocate( lZmat(0:l_max) )
      allocate( lWPmat(0:l_max) )
      allocate( lBmat(0:l_max) )

   end subroutine initialize_matrices
!------------------------------------------------------------------------------
end module matrices
