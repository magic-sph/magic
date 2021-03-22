#include "perflib_preproc.cpp"
module algebra

   use precision_mod, only: cp
   use constants, only: one

   implicit none

   private

   public :: prepare_mat, solve_mat, solve_tridiag, prepare_tridiag, &
   &         solve_bordered, prepare_band, solve_band, prepare_bordered

   interface solve_mat
      module procedure solve_mat_real_rhs
      module procedure solve_mat_complex_rhs
      module procedure solve_mat_real_rhs_multi
   end interface solve_mat

   interface solve_tridiag
      module procedure solve_tridiag_real_rhs
      module procedure solve_tridiag_complex_rhs
      module procedure solve_tridiag_real_rhs_multi
   end interface solve_tridiag

   interface solve_band
      module procedure solve_band_real_rhs
      module procedure solve_band_real_rhs_multi
      module procedure solve_band_complex_rhs
   end interface solve_band

   interface solve_bordered
      module procedure solve_bordered_real_rhs
      module procedure solve_bordered_real_rhs_multi
      module procedure solve_bordered_complex_rhs
   end interface solve_bordered

contains

   subroutine solve_mat_complex_rhs(a,len_a,n,pivot,rhs)
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
      real(cp) :: tmpr(n), tmpi(n)
      integer :: info, i

      do i=1,n
         tmpr(i) = real(rhs(i))
         tmpi(i) = aimag(rhs(i))
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,1,a,len_a,pivot,tmpr,n,info)
      call sgetrs('N',n,1,a,len_a,pivot,tmpi,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,1,a,len_a,pivot,tmpr,n,info)
      call dgetrs('N',n,1,a,len_a,pivot,tmpi,n,info)
#endif

      do i=1,n
         rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
      end do

   end subroutine solve_mat_complex_rhs
!-----------------------------------------------------------------------------
   subroutine solve_mat_real_rhs_multi(a,len_a,n,pivot,rhs,nRHSs)
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
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),rhs(1:n,:),n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),rhs(1:n,:),n,info)
#endif

   end subroutine solve_mat_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine solve_mat_real_rhs(a,len_a,n,pivot,rhs)
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
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#endif

   end subroutine solve_mat_real_rhs
!-----------------------------------------------------------------------------
   subroutine prepare_mat(a,len_a,n,pivot,info)
      !
      ! LU decomposition of the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      real(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

      PERFON('LUdecomp') 
#ifdef WITH_LIBFLAME
      !$omp critical
#endif
#if (DEFAULT_PRECISION==sngl)
      call sgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#endif
#ifdef WITH_LIBFLAME
      !$omp end critical
#endif
      PERFOFF

   end subroutine prepare_mat
!-----------------------------------------------------------------------------
   subroutine solve_tridiag_real_rhs(dl,d,du,du2,n,pivot,rhs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower
      real(cp), intent(in) :: du(n-1)   ! Upper
      real(cp), intent(in) :: du2(n-2)  ! For pivot

      !-- Output: solution stored in rhs(n)
      real(cp), intent(inout) :: rhs(n)
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgttrs('N',n,1,dl,d,du,du2,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrs('N',n,1,dl,d,du,du2,pivot,rhs,n,info)
#endif

   end subroutine solve_tridiag_real_rhs
!-----------------------------------------------------------------------------
   subroutine solve_tridiag_complex_rhs(dl,d,du,du2,n,pivot,rhs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower
      real(cp), intent(in) :: du(n-1)   ! Upper
      real(cp), intent(in) :: du2(n-2)  ! For pivot

      !-- Output: solution stored in rhs(n)
      complex(cp), intent(inout) :: rhs(n)

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

   end subroutine solve_tridiag_complex_rhs
!-----------------------------------------------------------------------------
   subroutine solve_tridiag_real_rhs_multi(dl,d,du,du2,n,pivot,rhs,nRHSs)

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

      !-- Local variables:
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgttrs('N',n,nRHSs,dl,d,du,du2,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrs('N',n,nRHSs,dl,d,du,du2,pivot,rhs,n,info)
#endif

   end subroutine solve_tridiag_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine prepare_tridiag(dl,d,du,du2,n,pivot,info)

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

   end subroutine prepare_tridiag
!-----------------------------------------------------------------------------
   subroutine solve_band_real_rhs(A, lenA, kl, ku, pivotA, rhs)

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

   end subroutine solve_band_real_rhs
!-----------------------------------------------------------------------------
   subroutine solve_band_complex_rhs(A, lenA, kl, ku, pivotA, rhs)

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

   end subroutine solve_band_complex_rhs
!-----------------------------------------------------------------------------
   subroutine solve_band_real_rhs_multi(A, lenA, kl, ku, pivotA, rhs, nRHSs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA
      integer,  intent(in) :: nRHSs
      integer,  intent(in) :: pivotA(lenA)
      real(cp), intent(in) :: A(2*kl+ku+1,lenA)

      !-- Output variable
      real(cp), intent(inout) :: rhs(:,:)

      !-- Local variables:
      integer :: n_bands, info

      n_bands = 2*kl+ku+1

#if (DEFAULT_PRECISION==sngl)
      call sgbtrs('N', lenA, kl, ku, nRHSs, A, n_bands, pivotA, rhs, lenA, info)
#elif (DEFAULT_PRECISION==dble)
      call dgbtrs('N', lenA, kl, ku, nRHSs, A, n_bands, pivotA, rhs, lenA, info)
#endif

   end subroutine solve_band_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine prepare_band(A,lenA,kl,ku,pivot,info)

      !-- Input variables
      integer, intent(in) :: lenA
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      real(cp), intent(inout) :: A(2*kl+ku+1,lenA)
      integer,  intent(out)   :: pivot(lenA)
      integer,  intent(out)   :: info

      !-- Local variables
      integer :: n_bands

      n_bands = 2*kl+ku+1

#ifdef WITH_LIBFLAME
      !$omp critical
#endif
#if (DEFAULT_PRECISION==sngl)
      call sgbtrf(lenA, lenA, kl, ku, A, n_bands, pivot, info)
#elif (DEFAULT_PRECISION==dble)
      call dgbtrf(lenA, lenA, kl, ku, A, n_bands, pivot, info)
#endif
#ifdef WITH_LIBFLAME
      !$omp end critical
#endif

   end subroutine prepare_band
!-----------------------------------------------------------------------------
   subroutine solve_bordered_real_rhs(A1,A2,A3,A4,lenA1,n_boundaries, kl, &
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
      real(cp), intent(in) :: A3(n_boundaries,lenA1)
      real(cp), intent(in) :: A4(n_boundaries,n_boundaries)

      !-- Output variable
      real(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: nStart, n_bands_A1, info

      nStart = lenA1+1
      n_bands_A1 = 2*kl+ku+1

#if (DEFAULT_PRECISION==sngl)
      !-- Solve A1*w = rhs1
      call sgbtrs('N', lenA1, kl, ku, 1, A1, n_bands_A1, pivotA1, &
           &      rhs(1:lenA1), lenA1, info)

      !-- rhs2 <- rhs2-A3*rhs1
      call sgemv('N', n_boundaries, lenA1, -one, A3, n_boundaries,  &
           &     rhs(1:lenA1), 1, one, rhs(nStart:), 1)

      !-- Solve A4*y = rhs2
      call sgetrs('N', n_boundaries, 1, A4, n_boundaries, pivotA4, &
           &      rhs(nStart:), n_boundaries, info)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
      call sgemv('N', lenA1, n_boundaries, -one, A2, lenA1, &
           &      rhs(nStart:), 1, one, rhs(1:lenA1), 1)
#elif (DEFAULT_PRECISION==dble)
      !-- Solve A1*w = rhs1
      call dgbtrs('N', lenA1, kl, ku, 1, A1, n_bands_A1, pivotA1, &
           &      rhs(1:lenA1), lenA1, info)

      !-- rhs2 <- rhs2-A3*rhs1
      call dgemv('N', n_boundaries, lenA1, -one, A3, n_boundaries,  &
           &     rhs(1:lenA1), 1, one, rhs(nStart:), 1)

      !-- Solve A4*y = rhs2
      call dgetrs('N', n_boundaries, 1, A4, n_boundaries, pivotA4, &
           &      rhs(nStart:), n_boundaries, info)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
      call dgemv('N', lenA1, n_boundaries, -one, A2, lenA1, &
           &      rhs(nStart:), 1, one, rhs(1:lenA1), 1)
#endif

   end subroutine solve_bordered_real_rhs
!-----------------------------------------------------------------------------
   subroutine solve_bordered_complex_rhs(A1,A2,A3,A4,lenA1,n_boundaries, kl, &
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
      real(cp), intent(in) :: A3(n_boundaries,lenA1)
      real(cp), intent(in) :: A4(n_boundaries,n_boundaries)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: nStart, n_bands_A1, info
      real(cp) :: tmpr(lenRhs), tmpi(lenRhs)

      nStart = lenA1+1
      n_bands_A1 = 2*kl+ku+1

      tmpr(:) =  real(rhs)
      tmpi(:) = aimag(rhs)

#if (DEFAULT_PRECISION==sngl)
      !-- Solve A1*w = rhs1
      call sgbtrs('N', lenA1, kl, ku, 1, A1, n_bands_A1, pivotA1, &
           &      tmpr(1:lenA1), lenA1, info)
      call sgbtrs('N', lenA1, kl, ku, 1, A1, n_bands_A1, pivotA1, &
           &      tmpi(1:lenA1), lenA1, info)

      !-- rhs2 <- rhs2-A3*rhs1
      call sgemv('N', n_boundaries, lenA1, -one, A3, n_boundaries,  &
           &     tmpr(1:lenA1), 1, one, tmpr(nStart:), 1)
      call sgemv('N', n_boundaries, lenA1, -one, A3, n_boundaries,  &
           &     tmpi(1:lenA1), 1, one, tmpi(nStart:), 1)

      !-- Solve A4*y = rhs2
      call sgetrs('N', n_boundaries, 1, A4, n_boundaries, pivotA4, &
           &      tmpr(nStart:), n_boundaries, info)
      call sgetrs('N', n_boundaries, 1, A4, n_boundaries, pivotA4, &
           &      tmpi(nStart:), n_boundaries, info)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
      call sgemv('N', lenA1, n_boundaries, -one, A2, lenA1, &
           &      tmpr(nStart:), 1, one, tmpr(1:lenA1), 1)
      call sgemv('N', lenA1, n_boundaries, -one, A2, lenA1, &
           &      tmpi(nStart:), 1, one, tmpi(1:lenA1), 1)
#elif (DEFAULT_PRECISION==dble)
      !-- Solve A1*w = rhs1
      call dgbtrs('N', lenA1, kl, ku, 1, A1, n_bands_A1, pivotA1, &
           &      tmpr(1:lenA1), lenA1, info)
      call dgbtrs('N', lenA1, kl, ku, 1, A1, n_bands_A1, pivotA1, &
           &      tmpi(1:lenA1), lenA1, info)

      !-- rhs2 <- rhs2-A3*rhs1
      call dgemv('N', n_boundaries, lenA1, -one, A3, n_boundaries,  &
           &     tmpr(1:lenA1), 1, one, tmpr(nStart:), 1)
      call dgemv('N', n_boundaries, lenA1, -one, A3, n_boundaries,  &
           &     tmpi(1:lenA1), 1, one, tmpi(nStart:), 1)

      !-- Solve A4*y = rhs2
      call dgetrs('N', n_boundaries, 1, A4, n_boundaries, pivotA4, &
           &      tmpr(nStart:), n_boundaries, info)
      call dgetrs('N', n_boundaries, 1, A4, n_boundaries, pivotA4, &
           &      tmpi(nStart:), n_boundaries, info)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
      call dgemv('N', lenA1, n_boundaries, -one, A2, lenA1, &
           &      tmpr(nStart:), 1, one, tmpr(1:lenA1), 1)
      call dgemv('N', lenA1, n_boundaries, -one, A2, lenA1, &
           &      tmpi(nStart:), 1, one, tmpi(1:lenA1), 1)
#endif

      rhs(:)=cmplx(tmpr(:),tmpi(:),kind=cp)

   end subroutine solve_bordered_complex_rhs
!-----------------------------------------------------------------------------
   subroutine solve_bordered_real_rhs_multi(A1,A2,A3,A4,lenA1,n_boundaries,kl, &
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
      real(cp), intent(in) :: A3(n_boundaries,lenA1)
      real(cp), intent(in) :: A4(n_boundaries,n_boundaries)

      !-- Output variable
      real(cp), intent(inout) :: rhs(:,:)

      !-- Local variables:
      integer :: nStart, n_bands_A1, info

      nStart = lenA1+1
      n_bands_A1 = 2*kl+ku+1

#if (DEFAULT_PRECISION==sngl)
      !-- Solve A1*w = rhs1
      call sgbtrs('N', lenA1, kl, ku, nRHSs, A1, n_bands_A1, pivotA1, &
           &      rhs(1:lenA1,:), lenA1, info)

      !-- rhs2 <- rhs2-A3*rhs1
      call sgemm('N', 'N', n_boundaries, nRHSs, lenA1, -one, A3,  &
           &     n_boundaries, rhs(1:lenA1,:), lenA1, one, rhs(nStart:,:), &
           &     n_boundaries)

      !-- Solve A4*y = rhs2
      call sgetrs('N', n_boundaries, nRHSs, A4, n_boundaries, pivotA4, &
           &      rhs(nStart:,:), n_boundaries, info)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
      call sgemm('N', 'N', lenA1, nRHSs, n_boundaries, -one, A2, lenA1, &
           &      rhs(nStart:,:), n_boundaries, one, rhs(1:lenA1,:), lenA1)
#elif (DEFAULT_PRECISION==dble)
      !-- Solve A1*w = rhs1
      call dgbtrs('N', lenA1, kl, ku, nRHSs, A1, n_bands_A1, pivotA1, &
           &      rhs(1:lenA1,:), lenA1, info)

      !-- rhs2 <- rhs2-A3*rhs1
      call dgemm('N', 'N', n_boundaries, nRHSs, lenA1, -one, A3,  &
           &     n_boundaries, rhs(1:lenA1,:), lenA1, one, rhs(nStart:,:), &
           &     n_boundaries)

      !-- Solve A4*y = rhs2
      call dgetrs('N', n_boundaries, nRHSs, A4, n_boundaries, pivotA4, &
           &      rhs(nStart:,:), n_boundaries, info)

      !-- Assemble rhs1 <- rhs1-A2*rhs2
      call dgemm('N', 'N', lenA1, nRHSs, n_boundaries, -one, A2, lenA1, &
           &      rhs(nStart:,:), n_boundaries, one, rhs(1:lenA1,:), lenA1)
#endif

   end subroutine solve_bordered_real_rhs_multi
!-----------------------------------------------------------------------------
   subroutine prepare_bordered(A1,A2,A3,A4,lenA1,n_boundaries,kl,ku,pivotA1, &
              &                pivotA4,info)

      !-- Input variables
      integer, intent(in) :: n_boundaries
      integer, intent(in) :: lenA1
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      real(cp), intent(inout) :: A1(2*kl+ku+1,lenA1)
      real(cp), intent(inout) :: A2(lenA1,n_boundaries)
      real(cp), intent(inout) :: A3(n_boundaries,lenA1)
      real(cp), intent(inout) :: A4(n_boundaries,n_boundaries)
      integer,  intent(out)   :: pivotA1(lenA1)
      integer,  intent(out)   :: pivotA4(n_boundaries)
      integer,  intent(out)   :: info

      !-- Local variables
      integer :: n_bands_A1

      n_bands_A1 = 2*kl+ku+1

#ifdef WITH_LIBFLAME
      !$omp critical
#endif
#if (DEFAULT_PRECISION==sngl)
      !-- LU factorisation for the banded block
      call sgbtrf(lenA1, lenA1, kl, ku, A1, n_bands_A1, pivotA1, info)

      !-- Solve A1*v = A2 (on output v = A2)
      call sgbtrs('N', lenA1, kl, ku, n_boundaries, A1, n_bands_A1, pivotA1, &
           &      A2, lenA1, info)

      !-- Assemble the Schur complement of A1: A4 <- A4-A3*v
      call sgemm('N', 'N', n_boundaries, n_boundaries, lenA1, -one, A3,  &
           &     n_boundaries, A2, lenA1, one, A4,  n_boundaries)

      !-- LU factorisation of the Schur complement
      call sgetrf(n_boundaries, n_boundaries, A4, n_boundaries, pivotA4, info)
#elif (DEFAULT_PRECISION==dble)
      !-- LU factorisation for the banded block
      call dgbtrf(lenA1, lenA1, kl, ku, A1, n_bands_A1, pivotA1, info)

      !-- Solve A1*v = A2 (on output v = A2)
      call dgbtrs('N', lenA1, kl, ku, n_boundaries, A1, n_bands_A1, pivotA1, &
           &      A2, lenA1, info)

      !-- Assemble the Schur complement of A1: A4 <- A4-A3*v
      call dgemm('N', 'N', n_boundaries, n_boundaries, lenA1, -one, A3,  &
           &     n_boundaries, A2, lenA1, one, A4,  n_boundaries)

      !-- LU factorisation of the Schur complement
      call dgetrf(n_boundaries, n_boundaries, A4, n_boundaries, pivotA4, info)
#endif
#ifdef WITH_LIBFLAME
      !$omp end critical
#endif

   end subroutine prepare_bordered
!-----------------------------------------------------------------------------
end module algebra
