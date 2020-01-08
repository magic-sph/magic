module algebra

   use precision_mod, only: cp

   implicit none

   private

   public :: prepare_mat, solve_mat, solve_tridiag, prepare_tridiag, &
   &         prepare_band, solve_band

   interface solve_mat
      module procedure solve_mat_real_rhs
      module procedure solve_mat_complex_rhs
      module procedure solve_mat_complex_rhs_multi
   end interface solve_mat

   interface solve_tridiag
      module procedure solve_tridiag_real_rhs
      module procedure solve_tridiag_complex_rhs
      module procedure solve_tridiag_complex_rhs_multi
   end interface solve_tridiag

   interface solve_band
      module procedure solve_band_real_rhs
      module procedure solve_band_complex_rhs_multi
      module procedure solve_band_complex_rhs
   end interface solve_band

contains

   subroutine solve_mat_complex_rhs(a,len_a,n,pivot,rhs)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
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
   subroutine solve_mat_complex_rhs_multi(a,len_a,n,pivot,rhs,nRHSs)
      !
      !  This routine does the backward substitution into a lu-decomposed real
      !  matrix a (to solve a * x = bc ) simultaneously for nRHSs complex 
      !  vectors bc. On return the results are stored in the bc.                  
      !

      !-- Input variables:
      integer,  intent(in) :: n           ! dimension of problem
      integer,  intent(in) :: len_a       ! leading dimension of a
      integer,  intent(in) :: pivot(n)       ! pivot pointer of length n
      real(cp), intent(in) :: a(len_a,n)  ! real n X n matrix
      integer,  intent(in) :: nRHSs       ! number of right-hand sides

      complex(cp), intent(inout) :: rhs(:,:) ! on input RHS of problem

      !-- Local variables:
      real(cp) :: tmpr(n,nRHSs), tmpi(n,nRHSs)
      integer :: info, i, j

      do j=1,nRHSs
         do i=1,n
            tmpr(i,j) = real(rhs(i,j))
            tmpi(i,j) = aimag(rhs(i,j))
         end do
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,:),n,info)
      call sgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpi(1:n,:),n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,:),n,info)
      call dgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpi(1:n,:),n,info)
#endif

      do j=1,nRHSs
         do i=1,n
            rhs(i,j)=cmplx(tmpr(i,j),tmpi(i,j),kind=cp)
         end do
      end do

   end subroutine solve_mat_complex_rhs_multi
!-----------------------------------------------------------------------------
   subroutine solve_mat_real_rhs(a,len_a,n,pivot,rhs)
      !
      !     like the linpack routine
      !     backward substitution of vector b into lu-decomposed matrix a
      !     to solve  a * x = b for a single real vector b
      !
      !     sub sgefa must be called once first to initialize a and pivot
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
      !     like the linpack routine
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      real(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

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

   end subroutine prepare_mat
!-----------------------------------------------------------------------------
   subroutine solve_tridiag_real_rhs(dl,d,du,du2,n,pivot,rhs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower 
      real(cp), intent(in) :: du(n-1)   ! Lower 
      real(cp), intent(in) :: du2(n-2)  ! Upper

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
      real(cp), intent(in) :: du(n-1)   ! Lower 
      real(cp), intent(in) :: du2(n-2)  ! Upper

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
   subroutine solve_tridiag_complex_rhs_multi(dl,d,du,du2,n,pivot,rhs,nRHSs)

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: d(n)      ! Diagonal
      real(cp), intent(in) :: dl(n-1)   ! Lower 
      real(cp), intent(in) :: du(n-1)   ! Lower 
      real(cp), intent(in) :: du2(n-2)  ! Upper
      integer,  intent(in) :: nRHSs     ! Number of right-hand side

      !-- Output: solution stored in rhs(n)
      complex(cp), intent(inout) :: rhs(:,:)

      !-- Local variables:
      real(cp) :: tmpr(n,nRHSs), tmpi(n,nRHSs)
      integer :: info, i, j

      do j=1,nRHSs
         do i=1,n
            tmpr(i,j) = real(rhs(i,j))
            tmpi(i,j) = aimag(rhs(i,j))
         end do
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgttrs('N',n,nRHSs,dl,d,du,du2,pivot,tmpr,n,info)
      call sgttrs('N',n,nRHSs,dl,d,du,du2,pivot,tmpi,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgttrs('N',n,nRHSs,dl,d,du,du2,pivot,tmpr,n,info)
      call dgttrs('N',n,nRHSs,dl,d,du,du2,pivot,tmpi,n,info)
#endif

      do j=1,nRHSs
         do i=1,n
            rhs(i,j)=cmplx(tmpr(i,j),tmpi(i,j),kind=cp)
         end do
      end do

   end subroutine solve_tridiag_complex_rhs_multi
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
   subroutine solve_band_complex_rhs_multi(A, lenA, kl, ku, pivotA, rhs, nRHSs)
   
      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA
      integer,  intent(in) :: nRHSs
      integer,  intent(in) :: pivotA(lenA)
      real(cp), intent(in) :: A(2*kl+ku+1,lenA)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(:,:)

      !-- Local variables:
      integer :: n_bands, info, i, j
      real(cp) :: tmpr(lenA,nRHSs), tmpi(lenA,nRHSs)

      n_bands = 2*kl+ku+1

      do j=1,nRHSs
         do i=1,lenA
            tmpr(i,j) = real(rhs(i,j))
            tmpi(i,j) = aimag(rhs(i,j))
         end do
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgbtrs('N', lenA, kl, ku, nRHSs, A, n_bands, pivotA, tmpr, lenA, info)
      call sgbtrs('N', lenA, kl, ku, nRHSs, A, n_bands, pivotA, tmpi, lenA, info)
#elif (DEFAULT_PRECISION==dble)
      call dgbtrs('N', lenA, kl, ku, nRHSs, A, n_bands, pivotA, tmpr, lenA, info)
      call dgbtrs('N', lenA, kl, ku, nRHSs, A, n_bands, pivotA, tmpi, lenA, info)
#endif

      do j=1,nRHSs
         do i=1,lenA
            rhs(i,j)=cmplx(tmpr(i,j),tmpi(i,j),kind=cp)
         end do
      end do

   end subroutine solve_band_complex_rhs_multi
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
end module algebra
