#include "perflib_preproc.cpp"
module algebra

   use omp_lib
   use precision_mod
   use constants, only: one

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   implicit none

   private

   real(cp), parameter :: zero_tolerance=1.0e-15_cp

   public :: sgefa, sgesl, cgesl, cgeslML

contains

   subroutine cgesl(a,ia,n,ip,bc1)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side  
      !  vector. On return x is stored in bc1.                            
      !                                                                     

      !-- Input variables:
      integer,  intent(in) :: n         ! dimension of problem
      integer,  intent(in) :: ia        ! first dim of a
      integer,  intent(in) :: ip(*)     ! pivot pointer of legth n
      real(cp), intent(in) :: a(ia,*)    ! real n X n matrix

      !-- Output: solution stored in bc1(*)
      complex(cp), intent(inout) :: bc1(*) ! on input RHS of problem

      !-- Local variables:
      integer :: nm1,nodd,i,m
      integer :: k,k1
      complex(cp) :: c1

      nm1=n-1
      nodd=mod(n,2)

      !--   permute vectors b1
      do k=1,nm1
         m=ip(k)
         c1=bc1(m)
         bc1(m)=bc1(k)
         bc1(k)=c1
      end do

      !--   solve  l * y = b
      do k=1,n-2,2
         k1=k+1
         bc1(k1)=bc1(k1)-bc1(k)*a(k1,k)
         do i=k+2,n
            bc1(i)=bc1(i)-(bc1(k)*a(i,k)+bc1(k1)*a(i,k1))
         end do
      end do

      if ( nodd == 0 ) then
         bc1(n)=bc1(n)-bc1(nm1)*a(n,nm1)
      end if

      !--   solve  u * x = y
      do k=n,3,-2
         k1=k-1
         bc1(k)=bc1(k)*a(k,k)
         bc1(k1)=(bc1(k1)-bc1(k)*a(k1,k))*a(k1,k1)
         do i=1,k-2
            bc1(i)=bc1(i)-bc1(k)*a(i,k)-bc1(k1)*a(i,k1)
         end do
      end do

      if ( nodd == 0 ) then
         bc1(2)=bc1(2)*a(2,2)
         bc1(1)=(bc1(1)-bc1(2)*a(1,2))*a(1,1)
      else
         bc1(1)=bc1(1)*a(1,1)
      end if

   end subroutine cgesl
!-----------------------------------------------------------------------------
   subroutine cgeslML(a,ia,n,ip,bc,ldBc,nRHSs)
      !
      !  This routine does the backward substitution into a lu-decomposed real
      !  matrix a (to solve a * x = bc ) simultaneously for nRHSs complex 
      !  vectors bc. On return the results are stored in the bc.                  
      !

      !-- Input variables:
      integer,  intent(in) :: n           ! dimension of problem
      integer,  intent(in) :: ia          ! leading dimension of a
      integer,  intent(in) :: ip(n)       ! pivot pointer of length n
      integer,  intent(in) :: ldBc        ! leading dimension of bc
      real(cp), intent(in) :: a(ia,n)     ! real n X n matrix
      integer,  intent(in) :: nRHSs       ! number of right-hand sides

      complex(cp), intent(inout) :: bc(ldBc,nRHSs) ! on input RHS of problem

      !-- Local variables:
      integer :: nm1,nodd,i,m
      integer :: k,k1,nRHS,nRHS2,noddRHS
      complex(cp) :: help

      nm1    = n-1
      nodd   = mod(n,2)
      noddRHS= mod(nRHSs,2)

      !     permute vectors bc
      LIKWID_ON('perm')
      do nRHS=1,nRHSs
         do k=1,nm1
            m=ip(k)
            help       =bc(m,nRHS)
            bc(m,nRHS) =bc(k,nRHS)
            bc(k,nRHS) =help
         end do
      end do
      LIKWID_OFF('perm')

      !     solve  l * y = b

      LIKWID_ON('cgeslML_1')
      !write(*,"(A,I4,A,I2,A)") "OpenMP loop over ",(nRHSs-1)/2,&
      !     &" iterations on ",omp_get_num_threads()," threads"
      do nRHS=1,nRHSs-1,2
         nRHS2=nRHS+1

         !PERFON('sol_1')
         do k=1,n-2,2
            k1=k+1
            bc(k1,nRHS) =bc(k1,nRHS)-bc(k,nRHS)*a(k1,k)
            bc(k1,nRHS2)=bc(k1,nRHS2)-bc(k,nRHS2)*a(k1,k)
            do i=k+2,n
               bc(i,nRHS) =bc(i,nRHS) - (bc(k,nRHS)*a(i,k)+bc(k1,nRHS)*a(i,k1))
               bc(i,nRHS2)=bc(i,nRHS2) - (bc(k,nRHS2)*a(i,k)+bc(k1,nRHS2)*a(i,k1))
            end do
         end do
         if ( nodd == 0 ) then
            bc(n,nRHS) =bc(n,nRHS) -bc(nm1,nRHS)*a(n,nm1)
            bc(n,nRHS2)=bc(n,nRHS2)-bc(nm1,nRHS2)*a(n,nm1)
         end if
         !PERFOFF
         !     solve  u * x = y
         !PERFON('sol_2')
         do k=n,3,-2
            k1=k-1
            bc(k,nRHS)  =bc(k,nRHS)*a(k,k)
            bc(k1,nRHS) =(bc(k1,nRHS)-bc(k,nRHS)*a(k1,k))*a(k1,k1)
            bc(k,nRHS2) =bc(k,nRHS2)*a(k,k)
            bc(k1,nRHS2)=(bc(k1,nRHS2)-bc(k,nRHS2)*a(k1,k))*a(k1,k1)
            do i=1,k-2
               bc(i,nRHS)=bc(i,nRHS) - bc(k,nRHS)*a(i,k)-bc(k1,nRHS)*a(i,k1)
               bc(i,nRHS2)=bc(i,nRHS2) - bc(k,nRHS2)*a(i,k)-bc(k1,nRHS2)*a(i,k1)
            end do
         end do
         if ( nodd == 0 ) then
            bc(2,nRHS)=bc(2,nRHS)*a(2,2)
            bc(1,nRHS)=(bc(1,nRHS)-bc(2,nRHS)*a(1,2))*a(1,1)
            bc(2,nRHS2)=bc(2,nRHS2)*a(2,2)
            bc(1,nRHS2)=(bc(1,nRHS2)-bc(2,nRHS2)*a(1,2))*a(1,1)
         else
            bc(1,nRHS)=bc(1,nRHS)*a(1,1)
            bc(1,nRHS2)=bc(1,nRHS2)*a(1,1)
         end if
         !PERFOFF

      end do

      if ( noddRHS == 1 ) then
         nRHS=nRHSs

         do k=1,n-2,2
            k1=k+1
            bc(k1,nRHS)=bc(k1,nRHS)-bc(k,nRHS)*a(k1,k)
            do i=k+2,n
               bc(i,nRHS)=bc(i,nRHS) - (bc(k,nRHS)*a(i,k)+bc(k1,nRHS)*a(i,k1))
            end do
         end do
         if ( nodd == 0 ) bc(n,nRHS)=bc(n,nRHS)-bc(nm1,nRHS)*a(n,nm1)
         do k=n,3,-2
            k1=k-1
            bc(k,nRHS) =bc(k,nRHS)*a(k,k)
            bc(k1,nRHS)=(bc(k1,nRHS)-bc(k,nRHS)*a(k1,k))*a(k1,k1)
            do i=1,k-2
               bc(i,nRHS)=bc(i,nRHS) - bc(k,nRHS)*a(i,k)-bc(k1,nRHS)*a(i,k1)
            end do
         end do
         if ( nodd == 0 ) then
            bc(2,nRHS)=bc(2,nRHS)*a(2,2)
            bc(1,nRHS)=(bc(1,nRHS)-bc(2,nRHS)*a(1,2))*a(1,1)
         else
            bc(1,nRHS)=bc(1,nRHS)*a(1,1)
         end if

      end if
      LIKWID_OFF('cgeslML_1')

   end subroutine cgeslML
!-----------------------------------------------------------------------------
   subroutine sgesl(a,ia,n,ip,b)
      !
      !     like the linpack routine
      !     backward substitution of vector b into lu-decomposed matrix a
      !     to solve  a * x = b for a single real vector b
      !
      !     sub sgefa must be called once first to initialize a and ip
      !
      !     a: (input)  nxn real matrix
      !     n: (input)  size of a and b
      !     ip: (input) pivot pointer array of length n
      !     b: (in/output) rhs-vector on input, solution on output
      !

      !-- Input variables:
      integer,  intent(in) :: n      ! dim of problem
      integer,  intent(in) :: ia     ! first dim of a
      integer,  intent(in) :: ip(*)  ! pivot information
      real(cp), intent(in) :: a(ia,*)

      !-- Output: solution stored in b(n)
      real(cp), intent(inout) :: b(*)

      !-- Local variables:
      integer :: nm1,i
      integer :: k,k1,m,nodd
      real(cp) :: help


      nm1 =n-1
      nodd=mod(n,2)

      !--   permute vector b
      do k=1,nm1
         m   =ip(k)
         help=b(m)
         b(m)=b(k)
         b(k)=help
      end do

      !--   solve  l * y = b
      do k=1,n-2,2
         k1=k+1
         b(k1)=b(k1)-b(k)*a(k1,k)
         do i=k+2,n
            b(i)=b(i)-(b(k)*a(i,k)+b(k1)*a(i,k1))
         end do
      end do
      if ( nodd == 0 ) b(n)=b(n)-b(nm1)*a(n,nm1)

      !--   solve  u * x = y
      do k=n,3,-2
         k1=k-1
         b(k) =b(k)*a(k,k)
         b(k1)=(b(k1)-b(k)*a(k1,k))*a(k1,k1)
         do i=1,k-2
            b(i)=b(i)-(b(k)*a(i,k)+b(k1)*a(i,k1))
         end do
      end do
      if ( nodd == 0 ) then
         b(2)=b(2)*a(2,2)
         b(1)=(b(1)-b(2)*a(1,2))*a(1,1)
      else
         b(1)=b(1)*a(1,1)
      end if

   end subroutine sgesl
!-----------------------------------------------------------------------------
   subroutine sgefa(a,ia,n,ip,info)
      !
      !     like the linpack routine
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !
      !     a: (in/output) real nxn matrix on input, lu-decomposed matrix on output
      !     ia: (input) first dimension of a (must be >= n)
      !     n: (input) 2nd dimension and rank of a
      !     ip: (output) pivot pointer array
      !     info: (output) error message when  /=  0
      !

      !-- Input variables:
      integer,  intent(in) :: ia,n
      real(cp), intent(inout) :: a(ia,*)

      !-- Output variables:
      integer, intent(out) :: ip(*)   ! pivoting information
      integer, intent(out) :: info

      !-- Local variables:
      integer :: nm1,k,kp1,l,i,j
      real(cp) :: help

      if ( n <= 1 ) stop '45'

      info=0
      nm1 =n-1

      do k=1,nm1
         kp1=k+1
         l  =k

         do i=kp1,n
            if ( abs(a(i,k)) > abs(a(l,k)) ) l=i
         end do

         ip(k)=l

         if ( abs(a(l,k))  >  zero_tolerance ) then

            if ( l /= k ) then
               do i=1,n
                  help  =a(k,i)
                  a(k,i)=a(l,i)
                  a(l,i)=help
               end do
            end if

            help=one/a(k,k)
            do i=kp1,n
               a(i,k)=help*a(i,k)
            end do

            do j=kp1,n
               do i=kp1,n
                  a(i,j)=a(i,j)-a(k,j)*a(i,k)
               end do
            end do

         else
            info=k
         end if

      end do

      ip(n)=n
      if ( abs(a(n,n))  <=  zero_tolerance ) info=n
      if ( info > 0 ) return

      do i=1,n
         a(i,i)=one/a(i,i)
      end do

   end subroutine sgefa
!-----------------------------------------------------------------------------
end module algebra
