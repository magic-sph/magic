module algebra_loops
   !
   ! This module is defined to handle LU factorisations and linear solves for a
   ! series of nmat matrices. Those can be dense, banded or tridiagonal. Linear
   ! solves are written for complex right hand sides.
   !

   use constants, only: one
   use precision_mod

   implicit none

   private

   real(cp), parameter :: zero_tolerance=10.0_cp*epsilon(0.0_cp)

   public :: prepare_dense_all, prepare_tridiag_all, prepare_band_all, &
   &         solve_dense_all, solve_tridiag_all, solve_band_all,       &
   &         solve_band_real_all

contains

   subroutine prepare_dense_all(dat, pivot, n, nmat, info)
      !
      ! This subroutines preforms the LU decomposition of a series of nmat dense
      ! matrices of size (n,n).
      !

      !-- Input variables
      integer, intent(in) :: n    ! Matrix size
      integer, intent(in) :: nmat ! Number of matrices

      !-- Output variables
      real(cp), intent(inout) :: dat(n,n,nmat) ! Matrix
      integer,  intent(out) :: pivot(n,nmat)   ! Pivot
      integer,  intent(out) :: info ! An integer to return the success of the pivoting

      !-- Local variables:
      integer :: nm1,k,kp1,l,i,j,idx
      real(cp) :: help

      info=0
      nm1 =n-1

      !-- This external loop should be put on GPU
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(k,kp1,l,i,j,help) &
      !$omp reduction(max:info)
#endif
      do idx=1,nmat
         do k=1,nm1

            kp1=k+1
            l  =k

            do i=kp1,n
               if ( abs(dat(i,k,idx)) > abs(dat(l,k,idx)) ) l=i
            end do
            pivot(k,idx)=l

            if ( abs(dat(l,k,idx)) > zero_tolerance ) then
               if ( l /= k ) then
                  do i=1,n
                     help        =dat(k,i,idx)
                     dat(k,i,idx)=dat(l,i,idx)
                     dat(l,i,idx)=help
                  end do
               end if

               help=one/dat(k,k,idx)
               do i=kp1,n
                  dat(i,k,idx)=help*dat(i,k,idx)
               end do

               do j=kp1,n
                  do i=kp1,n
                     dat(i,j,idx)=dat(i,j,idx)-dat(k,j,idx)*dat(i,k,idx)
                  end do
               end do
            else
               info=k
            end if

         end do

         pivot(n,idx)=n
         if ( abs(dat(n,n,idx)) <= zero_tolerance ) info=n

         do i=1,n
            dat(i,i,idx)=one/dat(i,i,idx)
         end do

      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

   end subroutine prepare_dense_all
!------------------------------------------------------------------------------
   subroutine prepare_tridiag_all(dat, du2, pivot, n, nmat, info)
      !
      ! This subroutines preforms the LU decomposition of a series of nmat tridiag
      ! matrices of size (n).
      !

      !-- Input variables
      integer, intent(in) :: n    ! Matrix size
      integer, intent(in) :: nmat ! Number of matrices

      !-- Output variables
      real(cp), intent(inout) :: dat(3,n,nmat) ! Matrix
      real(cp), intent(out) :: du2(n-2,nmat) ! Help array
      integer,  intent(out) :: pivot(n,nmat)   ! Pivot
      integer,  intent(out) :: info ! An integer to return the success of the pivoting

      !-- Local variables
      integer :: i, idx
      real(cp) :: fact, temp

      print*, 'hello'
      info=0

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(i,fact,temp) &
      !$omp reduction(max:info)
#endif
      do idx=1,nmat
         !-- Initialize pivot(i) = i and du2(I) = 0
         do i = 1,n
            pivot(i,idx) = i
         end do

         du2(:,idx)=0.0_cp
         do i = 1,n-2
            if ( abs(dat(2,i,idx)) >= abs(dat(3,i,idx))) then
               !-- No row interchange required, eliminate DL(I)
               if ( dat(2,i,idx) > zero_tolerance ) then
                  fact = dat(3,i,idx)/dat(2,i,idx)
                  dat(3,i,idx) = fact
                  dat(2,i+1,idx) = dat(2,i+1,idx) - fact*dat(1,i+1,idx)
               end if
            else
               !-- Interchange rows I and I+1, eliminate DL(I)
               fact = dat(2,i,idx)/dat(3,i,idx)
               dat(2,i,idx) = dat(3,i,idx)
               dat(3,i,idx) = fact
               temp = dat(1,i+1,idx)
               dat(1,i+1,idx) = dat(2,i+1,idx)
               dat(2,i+1,idx) = temp - fact*dat(2,i+1,idx)
               du2(i,idx) = dat(1,i+2,idx)
               dat(1,i+2,idx) = -fact*dat(1,i+2,idx)
               pivot(i,idx) = i + 1
            end if
         end do

         i = n - 1
         if ( abs(dat(2,i,idx)) >= abs(dat(3,i,idx)) ) then
            if ( dat(2,i,idx) > zero_tolerance ) then
               fact = dat(3,i,idx)/dat(2,i,idx)
               dat(3,i,idx) = fact
               dat(2,i+1,idx) = dat(2,i+1,idx) - fact*dat(1,i,idx)
            end if
         else
            fact = dat(2,i,idx) / dat(3,i,idx)
            dat(2,i,idx) = dat(3,i,idx)
            dat(3,i,idx) = fact
            temp = dat(1,i,idx)
            dat(1,i,idx) = dat(2,i+1,idx)
            dat(2,i+1,idx) = temp - fact*dat(2,i+1,idx)
            pivot(i,idx) = i + 1
         end if

         !-- Check for a zero on the diagonal of u.
         outer: do i = 1,n
            if ( dat(2,i,idx) <= zero_tolerance ) then
               info = i
               exit outer
            end if
         end do outer

      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      print*, 'value of info is', info

   end subroutine prepare_tridiag_all
!------------------------------------------------------------------------------
   subroutine prepare_band_all(dat, pivot, n, kl, ku, nmat, info)
      !
      ! This subroutines preforms the LU decomposition of a series of nmat band
      ! matrices of size (kl+ku+1,n).
      !

      !-- Input variables
      integer, intent(in) :: n    ! Matrix size
      integer, intent(in) :: kl   ! Number of lower bands
      integer, intent(in) :: ku   ! Number of upper bands
      integer, intent(in) :: nmat ! Number of matrices

      !-- Output variables
      real(cp), intent(inout) :: dat(2*kl+ku+1,n,nmat) ! Matrix
      integer,  intent(out) :: pivot(n,nmat)   ! Pivot
      integer,  intent(out) :: info ! An integer to return the success of the pivoting

      !-- Local variables
      integer :: i, idx, j1, jz, i0, lm, l, k, j, nm1, j0, ju, kp1, mm, m
      real(cp) :: temp

      m = kl + ku + 1
      j0 = ku + 2
      j1 = min(n,m) - 1
      info = 0

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do &
      !$omp private(i,temp,j1,jz,i0,lm,l,k,j,nm1,j0,ju,kp1,mm,m) &
      !$omp reduction(max:info)
#endif
      do idx=1,nmat
         if ( j1 >= j0 ) then
            do jz = j0, j1
               i0 = m + 1 - jz
               do i = i0, kl
                  dat(i,jz,idx)=0.0_cp
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
                     dat(i,jz,idx)=0.0_cp
                  end do
               end if

               lm = min(kl,n-k)
               l = maxloc(abs(dat(m:m+lm,k,idx)),dim=1)+m-1

               pivot(k,idx) = l + k - m

               if ( abs(dat(l,k,idx)) > zero_tolerance ) then

                  if (l /= m) then
                     temp = dat(l,k,idx)
                     dat(l,k,idx) = dat(m,k,idx)
                     dat(m,k,idx) = temp
                  end if

                  !-- Compute multipliers
                  temp = -1.0_cp/dat(m,k,idx)
                  dat(m+1:,k,idx)=temp*dat(m+1:,k,idx)

                  !-- Row elimination
                  ju = min(max(ju,ku+pivot(k,idx)),n)
                  mm = m
                  if ( ju >=  kp1 ) then
                     do j = kp1, ju
                        l = l - 1
                        mm = mm - 1
                        temp = dat(l,j,idx)
                        if ( l /= mm) then
                           dat(l,j,idx) = dat(mm,j,idx)
                           dat(mm,j,idx) = temp
                        end if
                        dat(mm+1:mm+lm,j,idx)=dat(mm+1:mm+lm,j,idx) + &
                        &                          temp*dat(m+1:m+lm,k,idx)
                     end do
                  end if
               else
                  info = k
               end if
            end do
         end if

         pivot(n,idx) = n

         if ( abs(dat(m,n,idx)) <= zero_tolerance ) info = n

      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

   end subroutine prepare_band_all
!------------------------------------------------------------------------------
   subroutine solve_dense_all(dat, pivot, rhs, lm2l, l2nLMB2, n, nmat, llm, ulm)
      !
      ! This subroutine performs a series of linear solves for llm:ulm right hand
      ! sides and for a corresponding series of nmat dense matrices of size (n,n).
      !

      !-- Input variables
      integer,  intent(in) :: n ! Matrix size
      integer,  intent(in) :: nmat ! Number of matrices
      integer,  intent(in) :: llm ! Lower lm index
      integer,  intent(in) :: ulm ! Upper lm index
      real(cp), intent(in) :: dat(n,n,nmat) ! Matrix
      integer,  intent(in) :: pivot(n,nmat) ! Pivot
      integer,  intent(in) :: lm2l(:) ! lm to degree converter
      integer,  intent(in) :: l2nLMB2(0:) ! degree to matrix index converter

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n,llm:ulm) ! r.h.s. and solution of the linear problem

      !-- Local variables:
      integer :: i, m, k, k1, nLMB2, l, nm1, nodd, nRHS
      complex(cp) :: help

      nm1=n-1
      nodd=mod(n,2)

      !-- Single loop over lm's
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(i,m,k,k1,nLMB2,l,help)
#endif
      do nRHS=llm,ulm

         l=lm2l(nRHS)
         nLMB2=l2nLMB2(l)

         !-- Permute vectors rhs
         do k=1,nm1
            m=pivot(k,nLMB2)
            help       =rhs(m,nRHS)
            rhs(m,nRHS)=rhs(k,nRHS)
            rhs(k,nRHS)=help
         end do

         !-- Solve  l * y = b
         do k=1,n-2,2
            k1=k+1
            rhs(k1,nRHS)=rhs(k1,nRHS)-rhs(k,nRHS)*dat(k1,k,nLMB2)
            !DIR$ CONCURRENT
            do i=k+2,n
               rhs(i,nRHS)=rhs(i,nRHS)-(rhs(k,nRHS)*dat(i,k,nLMB2) + &
               &                      rhs(k1,nRHS)*dat(i,k1,nLMB2))
            end do
         end do
         if ( nodd == 0 ) then
            rhs(n,nRHS)=rhs(n,nRHS)-rhs(nm1,nRHS)*dat(n,nm1,nLMB2)
         end if

         !-- Solve  u * x = y
         do k=n,3,-2
            k1=k-1
            rhs(k,nRHS) =rhs(k,nRHS)*dat(k,k,nLMB2)
            rhs(k1,nRHS)=(rhs(k1,nRHS)-rhs(k,nRHS)*dat(k1,k,nLMB2)) * &
            &            dat(k1,k1,nLMB2)
            !DIR$ CONCURRENT
            do i=1,k-2
               rhs(i,nRHS)=rhs(i,nRHS)-rhs(k,nRHS)*dat(i,k,nLMB2) - &
               &           rhs(k1,nRHS)*dat(i,k1,nLMB2)
            end do
         end do
         if ( nodd == 0 ) then
            rhs(2,nRHS)=rhs(2,nRHS)*dat(2,2,nLMB2)
            rhs(1,nRHS)=(rhs(1,nRHS)-rhs(2,nRHS)*dat(1,2,nLMB2))*dat(1,1,nLMB2)
         else
            rhs(1,nRHS)=rhs(1,nRHS)*dat(1,1,nLMB2)
         end if

      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

   end subroutine solve_dense_all
!------------------------------------------------------------------------------
   subroutine solve_tridiag_all(dat, pivot, du2, rhs, lm2l, l2nLMB2, n, nmat, &
              &                 llm, ulm)
      !
      ! This subroutine performs a series of linear solves for llm:ulm right hand
      ! sides and for a corresponding series of nmat tridiagonal matrices of size n.
      !

      !-- Input variables
      integer,  intent(in) :: n ! Matrix size
      integer,  intent(in) :: nmat ! Number of matrices
      integer,  intent(in) :: llm ! Lower lm index
      integer,  intent(in) :: ulm ! Upper lm index
      real(cp), intent(in) :: dat(3,n,nmat) ! Matrix
      real(cp), intent(in) :: du2(n-2,nmat) ! Help vector
      integer,  intent(in) :: pivot(n,nmat) ! Pivot
      integer,  intent(in) :: lm2l(:) ! lm to degree converter
      integer,  intent(in) :: l2nLMB2(0:) ! degree to matrix index converter

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n,llm:ulm) ! r.h.s. and solution of the linear problem

      !-- Local variables
      complex(cp) :: temp
      integer :: i,nLMB2,l,nRHS

      !-- Single loop over lm's
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(i,nLMB2,l,temp)
#endif
      do nRHS=llm,ulm

         l=lm2l(nRHS)
         nLMB2=l2nLMB2(l)

         !-- Solve L*x = rhs.
         do i=1,n-1
            if ( pivot(i,nLMB2) == i ) then
               rhs(i+1,nRHS)=rhs(i+1,nRHS)-dat(3,i,nLMB2)*rhs(i,nRHS)
            else
               temp=rhs(i,nRHS)
               rhs(i,nRHS)  =rhs(i+1,nRHS)
               rhs(i+1,nRHS)=temp-dat(3,i,nLMB2)*rhs(i,nRHS)
            end if
         end do

         !-- Solve U*x = rhs.
         rhs(n,nRHS)  =rhs(n,nRHS)/dat(2,n,nLMB2)
         rhs(n-1,nRHS)=(rhs(n-1,nRHS)-dat(1,n,nLMB2)*rhs(n,nRHS))/dat(2,n-1,nLMB2)
         do i = n-2,1,-1
            rhs(i,nRHS)=(rhs(i,nRHS)-dat(1,i+1,nLMB2)*rhs(i+1,nRHS) - &
            &            du2(i,nLMB2)*rhs(i+2,nRHS))/dat(2,i,nLMB2)
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif


   end subroutine solve_tridiag_all
!------------------------------------------------------------------------------
   subroutine solve_band_all(dat, pivot, rhs, lm2l, l2nLMB2, n, kl, ku, nmat, &
              &              llm, ulm)
      !
      ! This subroutine performs a series of linear solves for llm:ulm right hand
      ! sides and for a corresponding series of nmat band matrices of size (kl+ku+1,n).
      !

      !-- Input variables
      integer,  intent(in) :: n ! Matrix size
      integer,  intent(in) :: kl ! Number of lower bands
      integer,  intent(in) :: ku ! Number of upper bands
      integer,  intent(in) :: nmat ! Number of matrices
      integer,  intent(in) :: llm ! Lower lm index
      integer,  intent(in) :: ulm ! Upper lm index
      real(cp), intent(in) :: dat(2*kl+ku+1,n,nmat) ! Matrix
      integer,  intent(in) :: pivot(n,nmat) ! Pivot
      integer,  intent(in) :: lm2l(:) ! lm to degree converter
      integer,  intent(in) :: l2nLMB2(0:) ! degree to matrix index converter

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n,llm:ulm) ! r.h.s. and solution of the linear problem

      !-- Local variables
      integer :: k,l,nLMB2,ll,lm,kb,la,lb,m,nm1,nRHS
      complex(cp) :: temp

      m = ku + kl + 1
      nm1 = n - 1

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(k,l,nLMB2,ll,temp,lm,kb,la,lb)
#endif
      do nRHS=llm,ulm

         ll=lm2l(nRHS)
         nLMB2=l2nLMB2(ll)

         !-- First solve Ly = rhs
         do k = 1, nm1
            lm = min(kl,n-k)
            l = pivot(k,nLMB2)
            temp = rhs(l,nRHS)
            if (l /= k) then
               rhs(l,nRHS) = rhs(k,nRHS)
               rhs(k,nRHS) = temp
            end if
            rhs(k+1:k+lm,nRHS)=rhs(k+1:k+lm,nRHS)+temp*dat(m+1:m+lm,k,nLMB2)
         end do

         !-- Solve u*x =y
         do kb = 1, n
            k = n + 1 - kb
            rhs(k,nRHS) = rhs(k,nRHS)/dat(m,k,nLMB2)
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            temp = -rhs(k,nRHS)
            rhs(lb:lb+lm-1,nRHS)=rhs(lb:lb+lm-1,nRHS)+temp*dat(la:la+lm-1,k,nLMB2)
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

   end subroutine solve_band_all
!------------------------------------------------------------------------------
   subroutine solve_band_real_all(dat, pivot, rhs, n, kl, ku, nmat, nrhs)
      !
      ! This subroutine performs a series of linear solves for llm:ulm right hand
      ! sides and for a corresponding series of nmat band matrices of size (kl+ku+1,n).
      !

      !-- Input variables
      integer,  intent(in) :: n ! Matrix size
      integer,  intent(in) :: kl ! Number of lower bands
      integer,  intent(in) :: ku ! Number of upper bands
      integer,  intent(in) :: nmat ! Number of matrices
      integer,  intent(in) :: nrhs ! Number of right hand sides
      real(cp), intent(in) :: dat(2*kl+ku+1,n,nmat) ! Matrix
      integer,  intent(in) :: pivot(n,nmat) ! Pivot

      !-- Output variables
      real(cp), intent(inout) :: rhs(n,nrhs,nmat) ! r.h.s. and solution of the linear problem

      !-- Local variables
      integer :: k,l,idx,lm,kb,la,lb,m,nm1,iRHS
      real(cp) :: temp

      m = ku + kl + 1
      nm1 = n - 1

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2) &
      !$omp private(k,l,temp,lm,kb,la,lb)
#endif
      do idx=1,nmat
         do iRHS=1,nrhs

            !-- First solve Ly = rhs
            do k = 1, nm1
               lm = min(kl,n-k)
               l = pivot(k,idx)
               temp = rhs(l,iRHS,idx)
               if (l /= k) then
                  rhs(l,iRHS,idx) = rhs(k,iRHS,idx)
                  rhs(k,iRHS,idx) = temp
               end if
               rhs(k+1:k+lm,iRHS,idx)=rhs(k+1:k+lm,iRHS,idx)+temp*dat(m+1:m+lm,k,idx)
            end do

            !-- Solve u*x =y
            do kb = 1, n
               k = n + 1 - kb
               rhs(k,iRHS,idx) = rhs(k,iRHS,idx)/dat(m,k,idx)
               lm = min(k,m) - 1
               la = m - lm
               lb = k - lm
               temp = -rhs(k,iRHS,idx)
               rhs(lb:lb+lm-1,iRHS,idx)=rhs(lb:lb+lm-1,iRHS,idx) + &
               &                        temp*dat(la:la+lm-1,k,idx)
            end do
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

   end subroutine solve_band_real_all
!------------------------------------------------------------------------------
end module algebra_loops
