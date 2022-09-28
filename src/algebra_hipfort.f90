module algebra_hipfort

#ifdef WITH_OMP_GPU

   use precision_mod, only: cp
   use constants, only: one
   use iso_c_binding
   use hipfort_check
   use hipfort_hipblas
   use hipfort_hipsolver
   use omp_lib
   use hipfort, only: hipDeviceSynchronize

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

   subroutine gpu_solve_mat_complex_rhs(a,len_a,n,pivot,rhs)

      !
      !  This routine does the backward substitution into a LU-decomposed real
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side
      !  vector. On return x is stored in bc1.
      !

      !-- Input variables:
      integer,  intent(in) :: n          ! dimension of problem
      integer,  intent(in) :: len_a      ! first dim of a
      integer,  intent(in) :: pivot(n)   ! pivot pointer of legth n
      real(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n) ! on input RHS of problem

      !-- Local variables:
      real(cp), allocatable, target :: tmpr(:), tmpi(:)
      type(c_ptr) :: handle = c_null_ptr
      integer, allocatable, target :: devInfo(:)
      real(cp), allocatable, target :: dWork_r(:)
      real(cp), allocatable, target :: dWork_i(:)
      integer(c_int) :: size_work_bytes_r ! size of workspace to pass to getrs
      integer(c_int) :: size_work_bytes_i ! size of workspace to pass to getrs
      integer :: i

      !-- Create handle
      call hipsolverCheck(hipsolverCreate(handle))

      allocate(tmpr(n), tmpi(n))
      tmpi(:) = 0.0_cp; tmpr(:) = 0.0_cp
      !$omp target enter data map(alloc : tmpi, tmpr)
      !$omp target update to(tmpi, tmpr)

      !-- Extract real and imag parts of input rhs matrix
      !$omp target teams distribute simd
      do i=1,n
         tmpr(i) = real(rhs(i))
         tmpi(i) = aimag(rhs(i))
      end do
      !$omp end target teams distribute simd

#if (DEFAULT_PRECISION==sngl)
      !-- Real
      !$omp target data use_device_addr(a, pivot, tmpr)
      call hipsolverCheck(hipsolverSgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, &
                   & c_loc(pivot), c_loc(tmpr), n, size_work_bytes_r))
      !$omp end target data

      !-- Imag
      !$omp target data use_device_addr(a, pivot, tmpi)
      call hipsolverCheck(hipsolverSgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, &
                   & c_loc(pivot), c_loc(tmpi), n, size_work_bytes_i))
      !$omp end target data

      !--
      allocate(dWork_i(size_work_bytes_i), dWork_r(size_work_bytes_r), devInfo(1))
      dWork_i(:) = 0.0_cp
      dWork_r(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork_i, dWork_r, devInfo)
      !$omp target update to(dWork_i, dWork_r, devInfo)

      !--
      !$omp target data use_device_addr(a, pivot, tmpr, tmpi, devInfo, dWork_i, dWork_r)
      call hipsolverCheck(hipsolverSgetrs(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
                   & c_loc(tmpr), n, c_loc(dWork_r), size_work_bytes_r, devInfo(1)))
      call hipsolverCheck(hipsolverSgetrs(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
                   & c_loc(tmpi), n, c_loc(dWork_i), size_work_bytes_i, devInfo(1)))
      !$omp end target data

      !--
      !$omp target exit data map(delete : dWork_i, dWork_r, devInfo)
      deallocate(dWork_i, dWork_r, devInfo)
#elif (DEFAULT_PRECISION==dble)
      !-- Real
      !$omp target data use_device_addr(a, pivot, tmpr)
      call hipsolverCheck(hipsolverDgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, &
                   & c_loc(pivot), c_loc(tmpr), n, size_work_bytes_r))
      !$omp end target data

      !-- Imag
      !$omp target data use_device_addr(a, pivot, tmpi)
      call hipsolverCheck(hipsolverDgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, &
                   & c_loc(pivot), c_loc(tmpi), n, size_work_bytes_i))
      !$omp end target data

      !--
      allocate(dWork_i(size_work_bytes_i), dWork_r(size_work_bytes_r), devInfo(1))
      dWork_i(:) = 0.0_cp
      dWork_r(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork_i, dWork_r, devInfo)
      !$omp target update to(dWork_i, dWork_r, devInfo)

      !-- TODO: Replace these two call by a single to hipsolverZgetrs direcly on input complex rhs
      !$omp target data use_device_addr(a, pivot, tmpr, tmpi, devInfo, dWork_i, dWork_r)
      call hipsolverCheck(hipsolverDgetrs(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
                   & c_loc(tmpr), n, c_loc(dWork_r), size_work_bytes_r, devInfo(1)))
      call hipsolverCheck(hipsolverDgetrs(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
                   & c_loc(tmpi), n, c_loc(dWork_i), size_work_bytes_i, devInfo(1)))
      !$omp end target data

      !--
      !$omp target exit data map(delete : devInfo, dWork_i, dWork_r)
      deallocate(dWork_i, dWork_r, devInfo)
#endif

      call hipCheck(hipDeviceSynchronize())

      !$omp target teams distribute parallel do
      do i=1,n
         rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
      end do
      !$omp end target teams distribute parallel do

      !$omp target exit data map(delete : tmpi, tmpr)
      deallocate(tmpi, tmpr)

      !-- Destroy handle
      call hipsolverCheck(hipsolverDestroy(handle))

   end subroutine gpu_solve_mat_complex_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_solve_mat_real_rhs_multi(a,len_a,n,pivot,rhs,nRHSs)
      !
      !  This routine does the backward substitution into a LU-decomposed real
      !  matrix a (to solve a * x = bc ) simultaneously for nRHSs real
      !  vectors bc. On return the results are stored in the bc.
      !

      !-- Input variables:
      integer,  intent(in) :: n           ! dimension of problem
      integer,  intent(in) :: len_a       ! leading dimension of a
      integer,  intent(in) :: pivot(n)    ! pivot pointer of length n
      real(cp), intent(in) :: a(len_a,n)  ! real n X n matrix
      integer,  intent(in) :: nRHSs       ! number of right-hand sides

      real(cp), intent(inout) :: rhs(:,:) ! on input RHS of problem

      !-- Local variables:
      type(c_ptr) :: handle = c_null_ptr
      integer, allocatable, target :: devInfo(:)
      real(cp), allocatable, target :: dWork(:)
      integer(c_int) :: size_work_bytes ! size of workspace to pass to getrs

      !-- Create handle
      call hipsolverCheck(hipsolverCreate(handle))

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a, pivot, rhs)
      call hipsolverCheck(hipsolverSgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, nRHSs, c_loc(a(1:n,1:n)), n, &
                   & c_loc(pivot(1:n)), c_loc(rhs(1:n,:)), n, size_work_bytes))
      !$omp end target data

      !--
      allocate(dWork(size_work_bytes), devInfo(1))
      dWork(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork, devInfo)
      !$omp target update to(dWork, devInfo)

      !--
      !$omp target data use_device_addr(a, pivot, rhs, dWork, devInfo)
      call hipsolverCheck(hipsolverSgetrs(handle, HIPSOLVER_OP_N, n, nRHSs, c_loc(a(1:n,1:n)), n, c_loc(pivot(1:n)), &
                   & c_loc(rhs(1:n,:)), n, c_loc(dWork), size_work_bytes, devInfo(1)))
      !$omp end target data

      !--
      !$omp target exit data map(delete : dWork, devInfo)
      deallocate(dWork, devInfo)
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(a, pivot, rhs)
      call hipsolverCheck(hipsolverDgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, nRHSs, c_loc(a(1:n,1:n)), n, &
                   & c_loc(pivot(1:n)), c_loc(rhs(1:n,:)), n, size_work_bytes))
      !$omp end target data

      !--
      allocate(dWork(size_work_bytes), devInfo(1))
      dWork(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork, devInfo)
      !$omp target update to(dWork, devInfo)

      !--
      !$omp target data use_device_addr(a, pivot, rhs, dWork, devInfo)
      call hipsolverCheck(hipsolverDgetrs(handle, HIPSOLVER_OP_N, n, nRHSs, c_loc(a(1:n,1:n)), n, c_loc(pivot(1:n)), &
                   & c_loc(rhs(1:n,:)), n, c_loc(dWork), size_work_bytes, devInfo(1)))
      !$omp end target data

      !--
      !$omp target exit data map(delete : dWork, devInfo)
      deallocate(dWork, devInfo)
#endif

      call hipCheck(hipDeviceSynchronize())

      !-- Destroy handle
      call hipsolverCheck(hipsolverDestroy(handle))

   end subroutine gpu_solve_mat_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine gpu_solve_mat_real_rhs(a,len_a,n,pivot,rhs)
      !
      !  Backward substitution of vector b into lu-decomposed matrix a
      !  to solve  a * x = b for a single real vector b
      !

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: len_a     ! first dim of a
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: a(len_a,n)

      !-- Output: solution stored in rhs(n)
      real(cp), intent(inout) :: rhs(n)

      !-- Local variables
      type(c_ptr) :: handle = c_null_ptr
      integer, allocatable, target :: devInfo(:)
      real(cp), allocatable, target :: dWork(:)
      integer(c_int) :: size_work_bytes ! size of workspace to pass to getrs

      !-- Create handle
      call hipsolverCheck(hipsolverCreate(handle))

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a, pivot, rhs)
      call hipsolverCheck(hipsolverSgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, &
                   & c_loc(pivot), c_loc(rhs), n, size_work_bytes))
      !$omp end target data

      !--
      allocate(dWork(size_work_bytes), devInfo(1))
      dWork(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork, devInfo)
      !$omp target update to(dWork, devInfo)

      !--
      !$omp target data use_device_addr(a, pivot, rhs, devInfo, dWork)
      call hipsolverCheck(hipsolverSgetrs(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
                   & c_loc(rhs), n, c_loc(dWork), size_work_bytes, devInfo(1)))
      !$omp end target data

      !--
      !$omp target exit data map(delete : dWork, devInfo)
      deallocate(dWork, devInfo)
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(a, pivot, rhs)
      call hipsolverCheck(hipsolverDgetrs_bufferSize(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, &
                   & c_loc(pivot), c_loc(rhs), n, size_work_bytes))
      !$omp end target data

      !--
      allocate(dWork(size_work_bytes), devInfo(1))
      dWork(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork, devInfo)
      !$omp target update to(dWork, devInfo)

      !--
      !$omp target data use_device_addr(a, pivot, rhs, devInfo, dWork)
      call hipsolverCheck(hipsolverDgetrs(handle, HIPSOLVER_OP_N, n, 1, c_loc(a), len_a, c_loc(pivot), &
                   & c_loc(rhs), n, c_loc(dWork), size_work_bytes, devInfo(1)))
      !$omp end target data

      !--
      !$omp target exit data map(delete : dWork, devInfo)
      deallocate(dWork, devInfo)
#endif

      call hipCheck(hipDeviceSynchronize())

      !-- Destroy handle
      call hipsolverCheck(hipsolverDestroy(handle))

   end subroutine gpu_solve_mat_real_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_prepare_mat(a,len_a,n,pivot,info)
      !
      ! LU decomposition of the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      real(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(:)   ! pivoting information
      integer, intent(out) :: info

      !-- Local variables
      type(c_ptr) :: handle = c_null_ptr
      integer, allocatable, target :: devInfo(:)
      real(cp), allocatable, target :: dWork(:)
      integer(c_int) :: size_work_bytes ! size of workspace to pass to getrf

      !-- Create handle
      call hipsolverCheck(hipsolverCreate(handle))

#ifdef WITH_LIBFLAME
      !$omp critical
#endif

#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(a)
      call hipsolverCheck(hipsolverSgetrf_bufferSize(handle, n, n, c_loc(a(1:n,1:n)), &
                        & n, size_work_bytes))
      !$omp end target data

      allocate(dWork(size_work_bytes), devInfo(1))
      dWork(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork, devInfo)
      !$omp target update to(devInfo, dWork)

      !$omp target data use_device_addr(a, dWork, pivot, devInfo)
      call hipsolverCheck(hipsolverSgetrf(handle, n, n, c_loc(a(1:n,1:n)), n, c_loc(dWork), &
                        & size_work_bytes, c_loc(pivot(1:n)), devInfo(1)))
      !$omp end target data
      !$omp target update from(devInfo)
      info = devInfo(1)

      !$omp target exit data map(delete : dWork, devInfo)
      deallocate(dWork, devInfo)
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(a)
      call hipsolverCheck(hipsolverDgetrf_bufferSize(handle, n, n, c_loc(a(1:n,1:n)), &
                        & n, size_work_bytes))
      !$omp end target data

      allocate(dWork(size_work_bytes), devInfo(1))
      dWork(:) = 0.0_cp
      devInfo(1) = 0
      !$omp target enter data map(alloc : dWork, devInfo)
      !$omp target update to(devInfo, dWork)

      !$omp target data use_device_addr(a, dWork, pivot, devInfo)
      call hipsolverCheck(hipsolverDgetrf(handle, n, n, c_loc(a(1:n,1:n)), n, c_loc(dWork), &
                        & size_work_bytes, c_loc(pivot(1:n)), devInfo(1)))
      !$omp end target data
      !$omp target update from(devInfo)
      info = devInfo(1)

      !$omp target exit data map(delete : dWork, devInfo)
      deallocate(dWork, devInfo)
#endif

#ifdef WITH_LIBFLAME
      !$omp end critical
#endif

      call hipCheck(hipDeviceSynchronize())

      !-- Destroy handle
      call hipsolverCheck(hipsolverDestroy(handle))

   end subroutine gpu_prepare_mat
!-----------------------------------------------------------------------------

#endif

end module algebra_hipfort
