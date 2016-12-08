module algebra

   use precision_mod, only: cp

   implicit none

   private

   public :: sgefa, sgesl, cgesl, cgeslML

contains

   subroutine cgesl(a,len_a,n,pivot,rhs)
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
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call cgetrs('N',n,1,cmplx(a,0.0_cp,kind=cp),len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call zgetrs('N',n,1,cmplx(a,0.0_cp,kind=cp),len_a,pivot,rhs,n,info)
#endif

   end subroutine cgesl
!-----------------------------------------------------------------------------
   subroutine cgeslML(a,len_a,n,pivot,rhs,nRHSs)
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
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call cgetrs('N',n,nRHSs,cmplx(a(1:n,1:n),0.0_cp,kind=cp),n,pivot(1:n), &
                  rhs(1:n,:),n,info)
#elif (DEFAULT_PRECISION==dble)
      call zgetrs('N',n,nRHSs,cmplx(a(1:n,1:n),0.0_cp,kind=cp),n,pivot(1:n), &
                  rhs(1:n,:),n,info)
#endif

   end subroutine cgeslML
!-----------------------------------------------------------------------------
   subroutine sgesl(a,len_a,n,pivot,rhs)
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

   end subroutine sgesl
!-----------------------------------------------------------------------------
   subroutine sgefa(a,len_a,n,pivot,info)
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

#if (DEFAULT_PRECISION==sngl)
      call sgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#endif

   end subroutine sgefa
!-----------------------------------------------------------------------------
end module algebra
