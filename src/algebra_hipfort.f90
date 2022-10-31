module algebra_hipfort

#ifdef WITH_OMP_GPU

   use precision_mod, only: cp
   use constants, only: one
   use iso_c_binding
   use hipfort_check
   use hipfort_hipblas
   use omp_lib

   implicit none

   private

   real(cp), parameter :: zero_tolerance=1.0e-15_cp

   public :: gpu_prepare_mat, gpu_solve_mat

   interface gpu_solve_mat
      module procedure gpu_solve_mat_real_rhs
      module procedure gpu_solve_mat_complex_rhs
      module procedure gpu_solve_mat_real_rhs_multi
   end interface gpu_solve_mat

contains

   subroutine gpu_solve_mat_complex_rhs(a,len_a,n,pivot,tmpr,tmpi,handle,devInfo)

      !
      !  This routine does the backward substitution into a LU-decomposed real
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side
      !  vector. On return x is stored in bc1.
      !

      !-- Input variables:
      integer,        intent(in) :: n          ! dimension of problem
      integer,        intent(in) :: len_a      ! first dim of a
      integer,        intent(in) :: pivot(n)   ! pivot pointer of legth n
      real(cp),       intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      real(cp),    intent(inout) :: tmpr(n) ! on input RHS of problem
      real(cp),    intent(inout) :: tmpi(n) ! on input RHS of problem
      type(c_ptr), intent(inout) :: handle
      integer,     intent(inout) :: devInfo(:)

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a, pivot, tmpr, tmpi, devInfo)
      call hipblasCheck(hipblasSgetrs(handle, HIPBLAS_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
      &                               c_loc(tmpr), n, c_loc(devInfo)))
      call hipblasCheck(hipblasSgetrs(handle, HIPBLAS_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
      &                               c_loc(tmpi), n, c_loc(devInfo)))
      !$omp end target data
#elif (DEFAULT_PRECISION==dble)
      !-- TODO: Replace these two calls by a single to hipblasZgetrs direcly on input complex rhs
      !$omp target data use_device_addr(a, pivot, tmpr, tmpi, devInfo)
      call hipblasCheck(hipblasDgetrs(handle, HIPBLAS_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
      &                               c_loc(tmpr), n, c_loc(devInfo)))
      call hipblasCheck(hipblasDgetrs(handle, HIPBLAS_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
      &                               c_loc(tmpi), n, c_loc(devInfo)))
      !$omp end target data
#endif

   end subroutine gpu_solve_mat_complex_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_solve_mat_real_rhs_multi(a,len_a,n,pivot,rhs,nRHSs,handle,devInfo)
      !
      !  This routine does the backward substitution into a LU-decomposed real
      !  matrix a (to solve a * x = bc ) simultaneously for nRHSs real
      !  vectors bc. On return the results are stored in the bc.
      !

      !-- Input variables:
      integer,  intent(in)       :: n           ! dimension of problem
      integer,  intent(in)       :: len_a       ! leading dimension of a
      integer,  intent(in)       :: pivot(n)    ! pivot pointer of length n
      real(cp), intent(in)       :: a(len_a,n)  ! real n X n matrix
      integer,  intent(in)       :: nRHSs       ! number of right-hand sides

      real(cp), intent(inout)    :: rhs(:,:) ! on input RHS of problem
      integer,  intent(inout)    :: devInfo(:)
      type(c_ptr), intent(inout) :: handle

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a, pivot, rhs, devInfo)
      call hipblasCheck(hipblasSgetrs(handle, HIPBLAS_OP_N, n, nRHSs, c_loc(a(1:n,1:n)), n, c_loc(pivot(1:n)), &
      &                               c_loc(rhs(1:n,:)), n, c_loc(devInfo)))
      !$omp end target data
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(a, pivot, rhs, devInfo)
      call hipblasCheck(hipblasDgetrs(handle, HIPBLAS_OP_N, n, nRHSs, c_loc(a(1:n,1:n)), n, c_loc(pivot(1:n)), &
      &                               c_loc(rhs(1:n,:)), n, c_loc(devInfo)))
      !$omp end target data
#endif

   end subroutine gpu_solve_mat_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine gpu_solve_mat_real_rhs(a,len_a,n,pivot,rhs,handle,devInfo)
      !
      !  Backward substitution of vector b into lu-decomposed matrix a
      !  to solve  a * x = b for a single real vector b
      !

      !-- Input variables:
      integer,  intent(in)       :: n         ! dim of problem
      integer,  intent(in)       :: len_a     ! first dim of a
      integer,  intent(in)       :: pivot(n)  ! pivot information
      real(cp), intent(in)       :: a(len_a,n)

      !-- Output: solution stored in rhs(n)
      real(cp), intent(inout)    :: rhs(n)
      type(c_ptr), intent(inout) :: handle
      integer, intent(inout)     :: devInfo(:)

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a, pivot, rhs, devInfo)
      call hipblasCheck(hipblasSgetrs(handle, HIPBLAS_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
      &                               c_loc(rhs), n, c_loc(devInfo)))
      !$omp end target data
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(a, pivot, rhs, devInfo)
      call hipblasCheck(hipblasDgetrs(handle, HIPBLAS_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
      &                               c_loc(rhs), n, c_loc(devInfo)))
      !$omp end target data
#endif

   end subroutine gpu_solve_mat_real_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_prepare_mat(a,len_a,n,pivot,info,handle,devInfo)
      !
      ! LU decomposition of the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in)       :: len_a,n
      real(cp), intent(inout)    :: a(len_a,n)
      integer,  intent(inout)    :: devInfo(:)
      type(c_ptr), intent(inout) :: handle

      !-- Output variables:
      integer, intent(out) :: pivot(:)   ! pivoting information
      integer, intent(out) :: info

#ifdef WITH_LIBFLAME
      !$omp critical
#endif

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a, pivot, devInfo)
      call hipblasCheck(hipblasSgetrf(handle, n, c_loc(a(1:n,1:n)), n, c_loc(pivot(1:n)), c_loc(devInfo)))
      !$omp end target data
      !$omp target update from(devInfo)
      info = devInfo(1)
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(a, pivot, devInfo)
      call hipblasCheck(hipblasDgetrf(handle, n, c_loc(a(1:n,1:n)), n, c_loc(pivot(1:n)), c_loc(devInfo)))
      !$omp end target data
      !$omp target update from(devInfo)
      info = devInfo(1)
#endif

#ifdef WITH_LIBFLAME
      !$omp end critical
#endif

   end subroutine gpu_prepare_mat
!-----------------------------------------------------------------------------

#endif

end module algebra_hipfort
