module matrices
   !
   !  This module contains matricies for internal time step
   !

   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, l_max, l_maxMag, n_r_totMag, &
                         n_r_tot
   use logic, only: l_single_matrix
   use precision_mod

   implicit none

   private
 
   !-- the matrices, already LU-decomposed:
   real(cp), public, allocatable :: p0Mat(:,:)     ! for l=m=0  
   real(cp), public, allocatable :: s0Mat(:,:)     ! for l=m=0  
   real(cp), public, allocatable :: ps0Mat(:,:)    ! for l=m=0  
   real(cp), public, allocatable :: sMat(:,:,:)
   real(cp), public, allocatable :: zMat(:,:,:) 
   real(cp), public, allocatable :: z10Mat(:,:)    ! for l=1,m=0 
   real(cp), public, allocatable :: wpMat(:,:,:) 
   real(cp), public, allocatable :: wpsMat(:,:,:) 
   real(cp), public, allocatable :: bMat(:,:,:)
   real(cp), public, allocatable :: jMat(:,:,:)
 
   !-- respective pivoting information:
   integer, public, allocatable :: p0Pivot(:)
   integer, public, allocatable :: s0Pivot(:)
   integer, public, allocatable :: ps0Pivot(:)
   integer, public, allocatable :: sPivot(:,:)
   integer, public, allocatable :: z10Pivot(:)
   integer, public, allocatable :: zPivot(:,:)
   integer, public, allocatable :: wpPivot(:,:)
   integer, public, allocatable :: wpsPivot(:,:)
   integer, public, allocatable :: bPivot(:,:)
   integer, public, allocatable :: jPivot(:,:)
 
   ! -- respective linesums of the matrices
   real(cp), public, allocatable :: wpMat_fac(:,:,:)
   real(cp), public, allocatable :: wpsMat_fac(:,:,:)
   real(cp), public, allocatable :: ps0Mat_fac(:)
#ifdef WITH_PRECOND_Z
   real(cp), public, allocatable :: zMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_Z10
   real(cp), public, allocatable :: z10Mat_fac(:)
#endif
#ifdef WITH_PRECOND_S
   real(cp), public, allocatable :: sMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), public, allocatable :: s0Mat_fac(:)
#endif
#ifdef WITH_PRECOND_BJ
   real(cp), public, allocatable :: bMat_fac(:,:)
   real(cp), public, allocatable :: jMat_fac(:,:)
#endif
   !--- Logicals that inform whether the respective matrix
   !    has been updated:           
   logical, public :: lZ10mat
   logical, public, allocatable :: lSmat(:)
   logical, public, allocatable :: lZmat(:)
   logical, public, allocatable :: lWPmat(:)
   logical, public, allocatable :: lWPSmat(:)
   logical, public, allocatable :: lBmat(:)

   public :: initialize_matrices

contains

   subroutine initialize_matrices

      !-- the matrices, already LU-decomposed:
      if ( l_single_matrix ) then
         allocate( ps0Mat(2*n_r_max,2*n_r_max) )      ! for l=m=0  
         allocate( wpsMat(3*n_r_max,3*n_r_max,l_max) )
         allocate( s0Mat(n_r_max,n_r_max) )      ! still needed for scond
      else
         allocate( p0Mat(n_r_max,n_r_max) )      ! for l=m=0  
         allocate( s0Mat(n_r_max,n_r_max) )      ! for l=m=0  
         allocate( sMat(n_r_max,n_r_max,l_max) )
         allocate( wpMat(2*n_r_max,2*n_r_max,l_max) )
      end if
      allocate( zMat(n_r_max,n_r_max,l_max) )
      allocate( z10Mat(n_r_max,n_r_max) )    ! for l=1,m=0 
      allocate( bMat(n_r_totMag,n_r_totMag,l_maxMag) )
      allocate( jMat(n_r_totMag,n_r_totMag,l_maxMag) )
      bytes_allocated = bytes_allocated+(3*n_r_max*n_r_max+ &
                        6*n_r_max*n_r_max*l_max+            &
                        2*n_r_totMag*n_r_totMag*l_maxMag)*SIZEOF_DEF_REAL

      !-- respective pivoting information:
      if ( l_single_matrix ) then
         allocate ( ps0Pivot(2*n_r_max) )
         allocate ( wpsPivot(3*n_r_max,l_max) )
         allocate( s0Pivot(n_r_max) ) ! still needed for scond
      else
         allocate( p0Pivot(n_r_max) )
         allocate( s0Pivot(n_r_max) )
         allocate( sPivot(n_r_max,l_max) )
         allocate( wpPivot(2*n_r_max,l_max) )
      end if
      allocate( z10Pivot(n_r_max) )
      allocate( zPivot(n_r_max,l_max) )
      allocate( bPivot(n_r_totMag,l_maxMag) )
      allocate( jPivot(n_r_totMag,l_maxMag) )
      bytes_allocated = bytes_allocated+(3*n_r_max+4*n_r_max*l_max+ &
                        2*n_r_totMag*l_maxMag)*SIZEOF_INTEGER

      if ( l_single_matrix ) then
         allocate(wpsMat_fac(3*n_r_max,2,l_max))
         allocate(ps0Mat_fac(2*n_r_max))
      else
         allocate(wpMat_fac(2*n_r_max,2,l_max))
      end if
      bytes_allocated = bytes_allocated+4*n_r_max*l_max*SIZEOF_DEF_REAL

#ifdef WITH_PRECOND_Z10
      allocate(z10Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_Z
      allocate(zMat_fac(n_r_max,l_max))
      bytes_allocated = bytes_allocated+n_r_max*l_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S
      allocate(sMat_fac(n_r_max,l_max))
      bytes_allocated = bytes_allocated+n_r_max*l_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S0
      allocate(s0Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_BJ
      allocate(bMat_fac(n_r_totMag,l_maxMag))
      allocate(jMat_fac(n_r_totMag,l_maxMag))
      bytes_allocated = bytes_allocated+2*n_r_totMag*l_maxMag*SIZEOF_DEF_REAL
#endif

      !--- Logicals that inform whether the respective matrix
      !    has been updated:           
      if ( l_single_matrix ) then
         allocate( lWPSmat(0:l_max) )
      else
         allocate( lSmat(0:l_max) )
         allocate( lWPmat(0:l_max) )
      end if
      allocate( lZmat(0:l_max) )
      allocate( lBmat(0:l_max) )
      bytes_allocated = bytes_allocated+4*(l_max+1)*SIZEOF_INTEGER

   end subroutine initialize_matrices
!------------------------------------------------------------------------------
end module matrices
