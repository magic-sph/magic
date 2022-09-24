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

   public :: gpu_prepare_mat, gpu_solve_mat, gpu_solve_tridiag, gpu_prepare_tridiag, &
   &         gpu_solve_bordered, gpu_prepare_band, gpu_solve_band, gpu_prepare_bordered

   interface gpu_solve_mat
      module procedure gpu_solve_mat_real_rhs
      module procedure gpu_solve_mat_complex_rhs
      module procedure gpu_solve_mat_real_rhs_multi
   end interface gpu_solve_mat

   interface gpu_solve_tridiag
      module procedure gpu_solve_tridiag_real_rhs
      module procedure gpu_solve_tridiag_complex_rhs
      module procedure gpu_solve_tridiag_real_rhs_multi
   end interface gpu_solve_tridiag

   interface gpu_solve_band
      module procedure gpu_solve_band_real_rhs
      module procedure gpu_solve_band_real_rhs_multi
      module procedure gpu_solve_band_complex_rhs
   end interface gpu_solve_band

   interface gpu_solve_bordered
      module procedure gpu_solve_bordered_real_rhs
      module procedure gpu_solve_bordered_real_rhs_multi
      module procedure gpu_solve_bordered_complex_rhs
   end interface gpu_solve_bordered

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
   subroutine gpu_solve_tridiag_real_rhs(dl,d,du,du2,n,pivot,rhs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower
      real(cp), intent(in) :: du(n-1)   ! Upper
      real(cp), intent(in) :: du2(n-2)  ! For pivot

      !-- Output: solution stored in rhs(n)
      real(cp), intent(inout) :: rhs(n)

#ifndef HIP_SPARSE
      !-- Local variables
      integer :: i, ip
      real(cp) :: temp

      !-- Solve L*x = rhs.
      do i = 1, n-1
         ip = pivot(i)
         temp = rhs(i+1-ip+i)-dl(i)*rhs(ip)
         rhs(i) = rhs(ip)
         rhs(i+1) = temp
      end do

      !-- Solve U*x = rhs.
      rhs(n) = rhs(n)/d(n)
      rhs(n-1) = (rhs(N-1)-du(n-1)*rhs(n)) / d(n-1)
      do  i = n-2,1,-1
         rhs(i) = (rhs(i)-du(i)*rhs(i+1)-du2(i)*rhs(i+2)) / d(i)
      end do
#else
#if (DEFAULT_PRECISION==sngl)
      call sgttrs('N',n,1,dl,d,du,du2,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrs('N',n,1,dl,d,du,du2,pivot,rhs,n,info)
      !call hipsparseDgtsv2StridedBatch_(handle,n,dl,d,du,rhs,batchCount,batchStride,pBuffer)
      !-- Possible to replace by hipsparseXgtsv2StridedBatch (include factorisation step aka sgttrf)
#endif
#endif

   end subroutine gpu_solve_tridiag_real_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_solve_tridiag_complex_rhs(dl,d,du,du2,n,pivot,rhs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower
      real(cp), intent(in) :: du(n-1)   ! Upper
      real(cp), intent(in) :: du2(n-2)  ! For pivot

      !-- Output: solution stored in rhs(n)
      complex(cp), intent(inout) :: rhs(n)

#ifndef HIP_SPARSE
      !-- Local variables
      integer :: i, ip
      complex(cp) :: temp

      !-- Solve L*x = rhs.
      do i = 1, n-1
         ip = pivot(i)
         temp = rhs(i+1-ip+i)-dl(i)*rhs(ip)
         rhs(i) = rhs(ip)
         rhs(i+1) = temp
      end do

      !-- Solve U*x = rhs.
      rhs(n) = rhs(n)/d(n)
      rhs(n-1) = (rhs(N-1)-du(n-1)*rhs(n)) / d(n-1)
      do  i = n-2,1,-1
         rhs(i) = (rhs(i)-du(i)*rhs(i+1)-du2(i)*rhs(i+2)) / d(i)
      end do
#else
      !-- Local variables
      real(cp) :: tmpr(n), tmpi(n)
      integer :: info, i

      do i=1,n
         tmpr(i) = real(rhs(i))
         tmpi(i) = aimag(rhs(i))
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgttrs('N',n,1,dl,d,du,du2,pivot,tmpr,n,info)
      call sgttrs('N',n,1,dl,d,du,du2,pivot,tmpi,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrs('N',n,1,dl,d,du,du2,pivot,tmpr,n,info)
      call dgttrs('N',n,1,dl,d,du,du2,pivot,tmpi,n,info)
#endif

      do i=1,n
         rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
      end do
#endif

   end subroutine gpu_solve_tridiag_complex_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_solve_tridiag_real_rhs_multi(dl,d,du,du2,n,pivot,rhs,nRHSs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower
      real(cp), intent(in) :: du(n-1)   ! Upper
      real(cp), intent(in) :: du2(n-2)  ! For pivot
      integer,  intent(in) :: nRHSs     ! Number of right-hand side

      !-- Output: solution stored in rhs(n)
      real(cp), intent(inout) :: rhs(:,:)

#ifndef HIP_SPARSE
      !-- Local variables
      integer :: i, nRHS
      real(cp) :: temp

      do nRHS = 1, nRHSs
         !-- Solve L*x = rhs.
         do i = 1, n-1
            if ( pivot(i) == i ) then
               rhs(i+1,nRHS) = rhs(i+1,nRHS) - dl(i)*rhs(i,nRHS)
            else
               temp = rhs(i,nRHS)
               rhs(i,nRHS) = rhs(i+1,nRHS)
               rhs(i+1,nRHS) = temp - dl(i)*RHS(i,nRHS)
            end if
         end do

         !-- Solve U*x = rhs.
         rhs(n,nRHS) = rhs(n,nRHS)/d(n)
         rhs(n-1,nRHS) = (rhs(n-1,nRHS)-du(n-1)*rhs(n,nRHS))/d(n-1)
         do i = n-2,1,-1
            rhs(i,nRHS) = (rhs(i,nRHS)-du(i)*rhs(i+1,nRHS)-du2(i)* &
            &              rhs(i+2,nRHS))/d(i)
         end do
      end do
#else
      !-- Local variables:
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgttrs('N',n,nRHSs,dl,d,du,du2,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrs('N',n,nRHSs,dl,d,du,du2,pivot,rhs,n,info)
#endif
#endif

   end subroutine gpu_solve_tridiag_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine gpu_prepare_tridiag(dl,d,du,du2,n,pivot,info)

      !-- Input variables:
      integer,  intent(in) :: n
      real(cp), intent(inout) :: d(n)
      real(cp), intent(inout) :: dl(n-1)
      real(cp), intent(inout) :: du(n-1)

      !-- Output variable:
      real(cp), intent(inout) :: du2(n-2)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

#ifndef HIP_SPARSE
      !--
      real(cp), parameter :: zero_tolerance=1.0e-15_cp

      !-- Local variables
      integer :: i
      real(cp) :: fact, temp

      info = 0
      !-- Initialize pivot(i) = i and du2(I) = 0
      do i = 1, n
         pivot(i) = i
      end do

      du2(:) = 0.0_cp
      do i = 1, n-2
         if ( abs(d(i)) >= abs(dl(i))) then
            !-- No row interchange required, eliminate DL(I)
            if ( d(i) > zero_tolerance ) then
               fact = dl(i)/d(i)
               dl(i) = fact
               d(i+1) = d(i+1) - fact*du(i)
            end if
         else
            !-- Interchange rows I and I+1, eliminate DL(I)
            fact = d(i)/dl(i)
            d(i) = dl(i)
            dl(i) = fact
            temp = du(i)
            du(i) = d(i+1)
            d(i+1) = temp - fact*d(i+1)
            du2(i) = du(i+1)
            du(i+1) = -fact*du(i+1)
            pivot(i) = i + 1
         end if
      end do

      i = n - 1
      if ( abs(d(i)) >= abs(dl(i)) ) then
         if ( d(i) > zero_tolerance ) then
            fact = dl(i)/d(i)
            dl(i) = fact
            d(i+1) = d(i+1) - fact*du(i)
         end if
      else
         fact = d(i) / dl(i)
         d(i) = dl(i)
         dl(i) = fact
         temp = du(i)
         du(i) = d(i+1)
         d(i+1) = temp - fact*d(i+1)
         pivot(i) = i + 1
      end if

      !-- Check for a zero on the diagonal of u.
      outer: do i = 1, n
         if ( d(i) <= zero_tolerance ) then
            info = i
            exit outer
         end if
      end do outer
#else
#ifdef WITH_LIBFLAME
      !$omp critical
#endif
#if (DEFAULT_PRECISION==sngl)
      call sgttrf(n,dl,d,du,du2,pivot(1:n),info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrf(n,dl,d,du,du2,pivot(1:n),info)
#endif
#ifdef WITH_LIBFLAME
      !$omp end critical
#endif
#endif

   end subroutine gpu_prepare_tridiag
!-----------------------------------------------------------------------------
#ifndef HIP_SPARSE
   subroutine gpu_solve_band_real_rhs(abd, n, kl, ku, pivot, rhs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: n
      integer,  intent(in) :: pivot(n)
      real(cp), intent(in) :: abd(2*kl+ku+1, n)

      !-- Output variable
      real(cp), intent(inout) :: rhs(n)

      !-- Local variables
      real(cp) :: t
      integer :: k, kb, l, la, lb, lm, m, nm1

      m = ku + kl + 1
      nm1 = n - 1

      !-- First solve Ly = rhs
      if ( kl /= 0 .and. nm1 >= 1) then
         do k = 1, nm1
            lm = min(kl,n-k)
            l = pivot(k)
            t = rhs(l)
            if (l /= k) then
               rhs(l) = rhs(k)
               rhs(k) = t
            end if
            rhs(k+1:k+lm)=rhs(k+1:k+lm)+t*abd(m+1:m+lm,k)
         end do
      end if

      !-- Solve u*x =y
      do kb = 1, n
         k = n + 1 - kb
         rhs(k) = rhs(k)/abd(m,k)
         lm = min(k,m) - 1
         la = m - lm
         lb = k - lm
         t = -rhs(k)
         rhs(lb:lb+lm-1)=rhs(lb:lb+lm-1)+t*abd(la:la+lm-1,k)
      end do

   end subroutine gpu_solve_band_real_rhs
#else
   subroutine gpu_solve_band_real_rhs(A, lenA, kl, ku, pivotA, rhs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA
      integer,  intent(in) :: pivotA(lenA)
      real(cp), intent(in) :: A(2*kl+ku+1,lenA)

      !-- Output variable
      real(cp), intent(inout) :: rhs(lenA)

      !-- Local variables:
      integer :: n_bands, info

      n_bands = 2*kl+ku+1

#if (DEFAULT_PRECISION==sngl)
      call sgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, rhs(:), lenA, info)
#elif (DEFAULT_PRECISION==dble)
      call dgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, rhs(:), lenA, info)
#endif

   end subroutine gpu_solve_band_real_rhs
#endif
!-----------------------------------------------------------------------------
#ifndef HIP_SPARSE
   subroutine gpu_solve_band_complex_rhs(abd, n, kl, ku, pivot, rhs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: n
      integer,  intent(in) :: pivot(n)
      real(cp), intent(in) :: abd(2*kl+ku+1, n)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(n)

      !-- Local variables
      complex(cp) :: t
      integer :: k, kb, l, la, lb, lm, m, nm1

      m = ku + kl + 1
      nm1 = n - 1

      !-- First solve Ly = rhs
      if ( kl /= 0 .and. nm1 >= 1) then
         do k = 1, nm1
            lm = min(kl,n-k)
            l = pivot(k)
            t = rhs(l)
            if (l /= k) then
               rhs(l) = rhs(k)
               rhs(k) = t
            end if
            rhs(k+1:k+lm)=rhs(k+1:k+lm)+t*abd(m+1:m+lm,k)
         end do
      end if

      !-- Solve u*x =y
      do kb = 1, n
         k = n + 1 - kb
         rhs(k) = rhs(k)/abd(m,k)
         lm = min(k,m) - 1
         la = m - lm
         lb = k - lm
         t = -rhs(k)
         rhs(lb:lb+lm-1)=rhs(lb:lb+lm-1)+t*abd(la:la+lm-1,k)
      end do

   end subroutine gpu_solve_band_complex_rhs
#else
   subroutine gpu_solve_band_complex_rhs(A, lenA, kl, ku, pivotA, rhs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA
      integer,  intent(in) :: pivotA(lenA)
      real(cp), intent(in) :: A(2*kl+ku+1,lenA)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(lenA)

      !-- Local variables:
      real(cp) :: tmpr(lenA), tmpi(lenA)
      integer :: info, i, n_bands

      n_bands = 2*kl+ku+1

      do i=1,lenA
         tmpr(i) = real(rhs(i))
         tmpi(i) = aimag(rhs(i))
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, tmpr(:), lenA, info)
      call sgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, tmpi(:), lenA, info)
#elif (DEFAULT_PRECISION==dble)
      call dgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, tmpr(:), lenA, info)
      call dgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, tmpi(:), lenA, info)
#endif

      do i=1,lenA
         rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
      end do

   end subroutine gpu_solve_band_complex_rhs
#endif
!-----------------------------------------------------------------------------
   subroutine gpu_solve_band_real_rhs_multi(abd, n, kl, ku, pivot, rhs, nRHSs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: n
      integer,  intent(in) :: nRHSs
      integer,  intent(in) :: pivot(n)
      real(cp), intent(in) :: abd(2*kl+ku+1, n)

      !-- Output variable
      real(cp), intent(inout) :: rhs(:,:)

      !-- Local variables
      real(cp) :: t
      integer :: k, kb, l, la, lb, lm, m, nm1, nRHS

      m = ku + kl + 1
      nm1 = n - 1

      !-- First solve Ly = rhs
      if ( kl /= 0 .and. nm1 >= 1) then
         do nRHS=1,nRHSs
            do k = 1, nm1
               lm = min(kl,n-k)
               l = pivot(k)
               t = rhs(l,nRHS)
               if (l /= k) then
                  rhs(l,nRHS) = rhs(k,nRHS)
                  rhs(k,nRHS) = t
               end if
               rhs(k+1:k+lm,nRHS)=rhs(k+1:k+lm,nRHS)+t*abd(m+1:m+lm,k)
            end do
         end do
      end if

      !-- Solve u*x =y
      do nRHS=1,nRHSs
         do kb = 1, n
            k = n + 1 - kb
            rhs(k,nRHS) = rhs(k,nRHS)/abd(m,k)
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -rhs(k,nRHS)
            rhs(lb:lb+lm-1,nRHS)=rhs(lb:lb+lm-1,nRHS)+t*abd(la:la+lm-1,k)
         end do
      end do

   end subroutine gpu_solve_band_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine gpu_prepare_band(abd,n,kl,ku,pivot,info)

      !-- Input variables
      integer, intent(in) :: n
      integer, intent(in) :: kl, ku
      real(cp), intent(inout) :: abd(2*kl+ku+1,n)

      !-- Output variables
      integer, intent(out) :: pivot(n)
      integer, intent(out) :: info

      !-- Local variables
      real(cp) ::  t
      integer :: i, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm, nm1

      m = kl + ku + 1
      info = 0

      j0 = ku + 2
      j1 = min(n,m) - 1
      if ( j1 >= j0 ) then
         do jz = j0, j1
            i0 = m + 1 - jz
            do i = i0, kl
               abd(i,jz) = 0.0_cp
            end do
         end do
      end if
      jz = j1
      ju = 0

      !-- Gaussian elimination
      nm1 = n - 1
      if (nm1 >= 1) then
         do k = 1, nm1
            kp1 = k + 1

            jz = jz + 1
            if ( jz <= n .and. kl >= 1 ) then
               do i = 1, kl
                  abd(i,jz) = 0.0_cp
               end do
            end if

            lm = min(kl,n-k)
            !l = isamax(lm+1,abd(m,k)) + m - 1
            l = maxloc(abs(abd(m:m+lm,k)),dim=1)+m-1

            pivot(k) = l + k - m

            if ( abs(abd(l,k)) > zero_tolerance ) then

               if (l /= m) then
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               end if

               !-- Compute multipliers
               t = -one/abd(m,k)
               abd(m+1:,k)=t*abd(m+1:,k)

               !-- Row elimination
               ju = min(max(ju,ku+pivot(k)),n)
               mm = m
               if ( ju >=  kp1 ) then
                  do j = kp1, ju
                     l = l - 1
                     mm = mm - 1
                     t = abd(l,j)
                     if ( l /= mm) then
                        abd(l,j) = abd(mm,j)
                        abd(mm,j) = t
                     end if
                     abd(mm+1:mm+lm,j)=abd(mm+1:mm+lm,j)+t*abd(m+1:m+lm,k)
                  end do
               end if
            else
               info = k
            end if
         end do
      end if

      pivot(n) = n

      if ( abs(abd(m,n)) <= zero_tolerance ) info = n

   end subroutine gpu_prepare_band
!-----------------------------------------------------------------------------
   subroutine gpu_solve_bordered_real_rhs(A1,A2,A3,A4,lenA1,n_boundaries, kl, &
              &                       ku,pivotA1,pivotA4,rhs,lenRhs)

      !-- Input variables
      integer,  intent(in) :: n_boundaries
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA1
      integer,  intent(in) :: lenRhs
      integer,  intent(in) :: pivotA1(lenA1)
      integer,  intent(in) :: pivotA4(n_boundaries)
      real(cp), intent(in) :: A1(2*kl+ku+1,lenA1)
      real(cp), intent(in) :: A2(lenA1,n_boundaries)
      real(cp), intent(in) :: A3(lenA1)
      real(cp), intent(in) :: A4(n_boundaries,n_boundaries)

      !-- Output variable
      real(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: nStart, k
      type(c_ptr) :: handle = c_null_ptr

      nStart = lenA1+1

      !-- Solve A1*w = rhs1
      call gpu_solve_band_real_rhs(A1, lenA1, kl, ku, pivotA1, rhs(1:lenA1))

      !-- rhs2 <- rhs2-A3*rhs1
      do k=1,lenA1
         rhs(nStart)=rhs(nStart)-A3(k)*rhs(k)
      end do

      !!$omp target update to(rhs, A1, pivotA1)

      !-- Solve A4*y = rhs2
      call gpu_solve_mat_real_rhs(A4, n_boundaries, n_boundaries, pivotA4, rhs(nStart:))

      !-- Create handle
      call hipblasCheck(hipblasCreate(handle))

      !-- Assemble rhs1 <- rhs1-A2*rhs2
#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(A2, rhs)
      call hipblasCheck(hipblasSgemv(handle, HIPBLAS_OP_N, lenA1, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(rhs(nStart:)), 1, one, c_loc(rhs(1:lenA1)), 1))
      !$omp end target data
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(A2, rhs)
      call hipblasCheck(hipblasDgemv(handle, HIPBLAS_OP_N, lenA1, n_boundaries, -one, c_loc(A2), &
                                   & lenA1, c_loc(rhs(nStart:)), 1, one, c_loc(rhs(1:lenA1)), 1))
      !$omp end target data
#endif

      !-- Destroy handle
      call hipblasCheck(hipblasDestroy(handle))

   end subroutine gpu_solve_bordered_real_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_solve_bordered_complex_rhs(A1,A2,A3,A4,lenA1,n_boundaries, kl, &
              &                          ku,pivotA1,pivotA4,rhs,lenRhs)

      !-- Input variables
      integer,  intent(in) :: n_boundaries
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA1
      integer,  intent(in) :: lenRhs
      integer,  intent(in) :: pivotA1(lenA1)
      integer,  intent(in) :: pivotA4(n_boundaries)
      real(cp), intent(in) :: A1(2*kl+ku+1,lenA1)
      real(cp), intent(in) :: A2(lenA1,n_boundaries)
      real(cp), intent(in) :: A3(lenA1)
      real(cp), intent(in) :: A4(n_boundaries,n_boundaries)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: nStart, n_bands_A1, info, k
      real(cp), allocatable, target :: tmpr(:), tmpi(:)
      type(c_ptr) :: handle = c_null_ptr

      allocate(tmpr(lenRhs), tmpi(lenRhs))
      nStart = lenA1+1

      tmpr(:) =  real(rhs)
      tmpi(:) = aimag(rhs)

      !-- Solve A1*w = rhs1
      call gpu_solve_band_real_rhs(A1, lenA1, kl, ku, pivotA1, tmpr(1:lenA1))
      call gpu_solve_band_real_rhs(A1, lenA1, kl, ku, pivotA1, tmpi(1:lenA1))

      !-- rhs2 <- rhs2-A3*rhs1
      do k=1,lenA1
         tmpr(nStart)=tmpr(nStart)-A3(k)*tmpr(k)
         tmpi(nStart)=tmpi(nStart)-A3(k)*tmpi(k)
      end do

      !-- Solve A4*y = rhs2
      call gpu_solve_mat_real_rhs(A4, n_boundaries, n_boundaries, pivotA4, tmpr(nStart:))
      call gpu_solve_mat_real_rhs(A4, n_boundaries, n_boundaries, pivotA4, tmpi(nStart:))

      !-- Create handle
      call hipblasCheck(hipblasCreate(handle))

      !$omp target enter data map(alloc: tmpi, tmpr)
      !$omp target update to(tmpi, tmpr)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(A2, tmpi, tmpr)
      call hipblasCheck(hipblasSgemv(handle, HIPBLAS_OP_N, lenA1, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(tmpr(nStart:)), 1, one, c_loc(tmpr(1:lenA1)), 1))
      call hipblasCheck(hipblasSgemv(handle, HIPBLAS_OP_N, lenA1, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(tmpi(nStart:)), 1, one, c_loc(tmpi(1:lenA1)), 1))
      !$omp end target data
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(A2, tmpi, tmpr)
      call hipblasCheck(hipblasDgemv(handle, HIPBLAS_OP_N, lenA1, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(tmpr(nStart:)), 1, one, c_loc(tmpr(1:lenA1)), 1))
      call hipblasCheck(hipblasDgemv(handle, HIPBLAS_OP_N, lenA1, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(tmpi(nStart:)), 1, one, c_loc(tmpi(1:lenA1)), 1))
      !$omp end target data
#endif
      !$omp target update from(tmpi, tmpr)

      rhs(:)=cmplx(tmpr(:),tmpi(:),kind=cp)

      !$omp target exit data map(delete: tmpi, tmpr)
      deallocate(tmpr, tmpi)

      !-- Destroy handle
      call hipblasCheck(hipblasDestroy(handle))

   end subroutine gpu_solve_bordered_complex_rhs
!-----------------------------------------------------------------------------
   subroutine gpu_solve_bordered_real_rhs_multi(A1,A2,A3,A4,lenA1,n_boundaries,kl, &
              &                             ku,pivotA1,pivotA4,rhs,nRHSs)

      !-- Input variables
      integer,  intent(in) :: n_boundaries
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA1
      integer,  intent(in) :: nRHSs
      integer,  intent(in) :: pivotA1(lenA1)
      integer,  intent(in) :: pivotA4(n_boundaries)
      real(cp), intent(in) :: A1(2*kl+ku+1,lenA1)
      real(cp), intent(in) :: A2(lenA1,n_boundaries)
      real(cp), intent(in) :: A3(lenA1)
      real(cp), intent(in) :: A4(n_boundaries,n_boundaries)

      !-- Output variable
      real(cp), intent(inout) :: rhs(:,:)

      !-- Local variables:
      integer :: nStart, k
      type(c_ptr) :: handle = c_null_ptr

      nStart = lenA1+1

      !-- Solve A1*w = rhs1
      call gpu_solve_band_real_rhs_multi(A1, lenA1, kl, ku, pivotA1, rhs(1:lenA1,:), nRHSs)

      !-- rhs2 <- rhs2-A3*rhs1
      do k=1,lenA1
         rhs(nStart,:)=rhs(nStart,:)-A3(k)*rhs(k,:)
      end do

      !-- Solve A4*y = rhs2
      call gpu_solve_mat_real_rhs_multi(A4, n_boundaries, n_boundaries, pivotA4, rhs(nStart:,:), nRHSs)

      !-- Create handle
      call hipblasCheck(hipblasCreate(handle))

      !-- Assemble rhs1 <- rhs1-A2*rhs2
#if (DEFAULT_PRECISION==sngl)
      !$omp target data use_device_addr(A2, rhs)
      call hipblasCheck(hipblasSgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, lenA1, nRHSs, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(rhs(nStart:,:)), n_boundaries, one, c_loc(rhs(1:lenA1,:)), lenA1))
      !$omp end target data
#elif (DEFAULT_PRECISION==dble)
      !$omp target data use_device_addr(A2, rhs)
      call hipblasCheck(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, lenA1, nRHSs, n_boundaries, -one, c_loc(A2), &
                               & lenA1, c_loc(rhs(nStart:,:)), n_boundaries, one, c_loc(rhs(1:lenA1,:)), lenA1))
      !$omp end target data
#endif

      !-- Destroy handle
      call hipblasCheck(hipblasDestroy(handle))

   end subroutine gpu_solve_bordered_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine gpu_prepare_bordered(A1,A2,A3,A4,lenA1,n_boundaries,kl,ku,pivotA1, &
              &                pivotA4,info)

      !-- Input variables
      integer, intent(in) :: n_boundaries
      integer, intent(in) :: lenA1
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      real(cp), intent(inout) :: A1(2*kl+ku+1,lenA1)
      real(cp), intent(inout) :: A2(lenA1,n_boundaries)
      real(cp), intent(inout) :: A3(lenA1)
      real(cp), intent(inout) :: A4(n_boundaries,n_boundaries)
      integer,  intent(out)   :: pivotA1(lenA1)
      integer,  intent(out)   :: pivotA4(n_boundaries)
      integer,  intent(out)   :: info

      !-- Local variables
      integer :: i, j

      !-- LU factorisation for the banded block
      call gpu_prepare_band(A1, lenA1, kl, ku, pivotA1, info)

      !-- Solve A1*v = A2 (on output v = A2)
      call gpu_solve_band_real_rhs_multi(A1, lenA1, kl, ku, pivotA1, A2, n_boundaries)

      !-- Assemble the Schur complement of A1: A4 <- A4-A3*v
      do i=1,n_boundaries
         do j=1,lenA1
            A4(1,i)=A4(1,i)-A3(j)*A2(j,i)
         end do
      end do

      !-- LU factorisation of the Schur complement
      call gpu_prepare_mat(A4, n_boundaries, n_boundaries, pivotA4, info)

   end subroutine gpu_prepare_bordered
!-----------------------------------------------------------------------------

#endif

end module algebra_hipfort
