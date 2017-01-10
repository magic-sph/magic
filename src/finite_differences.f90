module finite_differences
   !
   ! This module is used to calculate the radial grid when finite differences
   ! are used
   !

   use precision_mod
   use parallel_mod, only: rank
   use constants, only: one, two
   use truncation, only: n_r_max
   use useful, only: logWrite
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_stencil
      integer :: order
      real(cp), allocatable :: dr(:,:)
      real(cp), allocatable :: ddr(:,:)
      real(cp), allocatable :: dddr(:,:)
      real(cp), allocatable :: dr_top(:,:)
      real(cp), allocatable :: dr_bot(:,:)
      real(cp), allocatable :: ddr_top(:,:)
      real(cp), allocatable :: ddr_bot(:,:)
      real(cp), allocatable :: dddr_top(:,:)
      real(cp), allocatable :: dddr_bot(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: get_FD_coeffs
      procedure :: get_FD_matder
   end type type_stencil


   public :: get_FD_grid

contains

   subroutine initialize(this,n_r_max,order)
      !
      ! This subroutine allocates the arrays used when finite difference are used
      !

      class(type_stencil) :: this
      integer, intent(in) :: n_r_max ! Number of radial grid points
      integer, intent(in) :: order   ! FD order

      this%order = order

      allocate( this%dr(n_r_max,0:order) )
      allocate( this%ddr(n_r_max,0:order) )
      allocate( this%dddr(n_r_max,0:order+2) )
      allocate( this%dr_top(order/2,0:order) )
      allocate( this%dr_bot(order/2,0:order) )
      allocate( this%ddr_top(order/2,0:order+1) )
      allocate( this%ddr_bot(order/2,0:order+1) )
      allocate( this%dddr_top(order/2+1,0:order+2) )
      allocate( this%dddr_bot(order/2+1,0:order+2) )

      bytes_allocated=bytes_allocated+(n_r_max*(3*order+5)+         &
      &               order/2*(4*order+6)+(order/2+1)*(2*order+6))* &
      &               SIZEOF_DEF_REAL

   end subroutine initialize
!---------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! This subroutine deallocates the arrays used in FD
      !

      class(type_stencil) :: this

      deallocate( this%dr, this%ddr, this%dddr )
      deallocate( this%dr_top, this%dr_bot )
      deallocate( this%ddr_top, this%ddr_bot )
      deallocate( this%dddr_top, this%dddr_bot )

   end subroutine finalize
!---------------------------------------------------------------------------
   subroutine get_FD_grid(ratio1, ratio2, ricb, rcmb, r)
      !
      ! This subroutine constructs the radial grid
      !

      !-- Input quantities:
      real(cp), intent(inout) :: ratio1  ! Nboudary/Nbulk
      real(cp), intent(in) :: ratio2     ! drMin/drMax
      real(cp), intent(in) :: ricb       ! inner boundary
      real(cp), intent(in) :: rcmb       ! outer boundary

      !-- Output quantities:
      real(cp), intent(out) :: r(:) ! radius

      !-- Local quantities:
      real(cp) :: dr_before, dr_after
      character(len=80) :: message
      integer :: n_r
      integer :: n_bulk_points, n_boundary_points

      r(1)=rcmb ! start with the outer boundary

      if ( ratio2 == one ) then ! Regular grid

         dr_before = one/(real(n_r_max,cp)-one)
         dr_after = dr_before
         do n_r=2,n_r_max
            r(n_r)=r(n_r-1)-dr_before
         end do

      else  ! irregular grid

         n_boundary_points = int( real(n_r_max-1,cp)/(two*(one+ratio1)) )
         ratio1 = real(n_r_max-1, cp)/real(two*n_boundary_points)-one

         n_bulk_points = n_r_max-1-2*n_boundary_points
         dr_after  = exp( log(ratio2)/real(n_boundary_points,cp) )
         dr_before = one

         do n_r=1,n_boundary_points
            dr_before=dr_before*dr_after
         end do
         dr_before=one/(real(n_bulk_points,cp)+two*dr_after*((one-dr_before) &
         &         /(one-dr_after)))

         write(message,'(''!      drMax='',ES16.6)') dr_before
         call logWrite(message)

         do n_r=1,n_boundary_points
            dr_before = dr_before*dr_after
         end do

         write(message,'(''!      drMax='',ES16.6)') dr_before
         call logWrite(message)

         do n_r=2,n_boundary_points+1
            r(n_r) = r(n_r-1)-dr_before
            dr_before = dr_before/dr_after
         end do

         do n_r=1,n_bulk_points
            r(n_r+n_boundary_points+1)=r(n_r+n_boundary_points)-dr_before
         end do

         do n_r=1,n_boundary_points
            dr_before = dr_before*dr_after
            r(n_r+n_boundary_points+n_bulk_points+1)=         &
            &        r(n_r+n_boundary_points+n_bulk_points)-dr_before
         end do

      end if

      if ( abs(r(n_r_max)-ricb) > 10.0_cp*epsilon(one) ) then
         if ( rank == 0 ) then
            write(*,*) 'Wrong internal radius!'
            stop
         end if
      end if

   end subroutine get_FD_grid
!---------------------------------------------------------------------------
   subroutine get_FD_coeffs(this, r)

      class(type_stencil) :: this

      !-- Input quantities:
      real(cp), intent(in) :: r(:) ! Radius

      !-- Local quantities:
      real(cp), allocatable :: dr_spacing(:)
      real(cp), allocatable :: taylor_exp(:,:)
      real(cp), allocatable :: taylor_exp_inv(:,:)
      real(cp) :: weight
      integer :: n_r, od, od_in

      allocate( dr_spacing(this%order+1) )
      allocate( taylor_exp(0:this%order,0:this%order) )
      allocate( taylor_exp_inv(0:this%order,0:this%order) )

      !
      !-- Step 1: First and 2nd derivatives in the bulk
      !
      do n_r=1+this%order/2,n_r_max-this%order/2
         do od=0,this%order
            dr_spacing(od+1)=r(n_r-this%order/2+od)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+1)

         do od=0,this%order
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+1)

         do od_in=0,this%order
            !-- Preconditioning
            do od=0,this%order
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%dr(n_r,od_in) =taylor_exp_inv(1,od_in)
            this%ddr(n_r,od_in)=taylor_exp_inv(2,od_in)
         end do
      end do

      !
      !-- Step 2: First derivative for the outer points
      !
      do n_r=1,this%order/2
         do od=0,this%order
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+1)

         do od=0,this%order
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+1)

         do od_in=0,this%order
            !-- Preconditioning
            do od=0,this%order
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%dr_top(n_r,od_in) =taylor_exp_inv(1,od_in)
         end do
      end do

      !
      !-- Step 3: First derivative for the inner points
      !
      do n_r=1,this%order/2
         do od=0,this%order
            dr_spacing(od+1)=r(n_r_max-od)-r(n_r_max-n_r+1)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+1)

         do od=0,this%order
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+1)

         do od_in=0,this%order
            !-- Preconditioning
            do od=0,this%order
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%dr_bot(n_r,od_in) =taylor_exp_inv(1,od_in)
         end do
      end do

      deallocate( dr_spacing, taylor_exp, taylor_exp_inv )

      !
      !-- Step 4: 2nd derivative for the outer points
      !
      allocate( dr_spacing(this%order+2) )
      allocate( taylor_exp(0:this%order+1,0:this%order+1) )
      allocate( taylor_exp_inv(0:this%order+1,0:this%order+1) )

      do n_r=1,this%order/2
         do od=0,this%order+1
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+2)

         do od=0,this%order+1
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+2)

         do od_in=0,this%order+1
            !-- Preconditioning
            do od=0,this%order+1
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%ddr_top(n_r,od_in) =taylor_exp_inv(2,od_in)
         end do
      end do

      !
      !-- Step 5: 2nd derivative for the inner points
      !
      do n_r=1,this%order/2
         do od=0,this%order+1
            dr_spacing(od+1)=r(n_r_max-od)-r(n_r_max-n_r+1)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+2)

         do od=0,this%order+1
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+2)

         do od_in=0,this%order+1
            !-- Preconditioning
            do od=0,this%order+1
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%ddr_bot(n_r,od_in) =taylor_exp_inv(2,od_in)
         end do
      end do

      deallocate( dr_spacing, taylor_exp, taylor_exp_inv )

      !
      !-- Step 6: 3rd derivative in the bulk
      !
      allocate( dr_spacing(this%order+3) )
      allocate( taylor_exp(0:this%order+2,0:this%order+2) )
      allocate( taylor_exp_inv(0:this%order+2,0:this%order+2) )

      do n_r=2+this%order/2,n_r_max-this%order/2-1
         do od=0,this%order+2
            dr_spacing(od+1)=r(n_r-this%order/2-1+od)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+3)

         do od=0,this%order+2
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+3)

         do od_in=0,this%order+2
            !-- Preconditioning
            do od=0,this%order+2
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%dddr(n_r,od_in) =taylor_exp_inv(3,od_in)
         end do
      end do

      !
      !-- Step 7: 3rd derivative for the outer points
      !
      do n_r=1,this%order/2+1
         do od=0,this%order+2
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+3)

         do od=0,this%order+2
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+3)

         do od_in=0,this%order+2
            !-- Preconditioning
            do od=0,this%order+2
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%dddr_top(n_r,od_in) =taylor_exp_inv(3,od_in)
         end do
      end do

      !
      !-- Step 8: 3rd derivative for the inner points
      !
      do n_r=1,this%order/2+1
         do od=0,this%order+2
            dr_spacing(od+1)=r(n_r_max-od)-r(n_r_max-n_r+1)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(this%order+3)

         do od=0,this%order+2
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, this%order+3)

         do od_in=0,this%order+2
            !-- Preconditioning
            do od=0,this%order+2
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            this%dddr_bot(n_r,od_in) =taylor_exp_inv(3,od_in)
         end do
      end do

      deallocate( dr_spacing, taylor_exp, taylor_exp_inv )

   end subroutine get_FD_coeffs
!---------------------------------------------------------------------------
   subroutine get_FD_matder(this, n_r_max, eye, drMat, d2rMat, d3rMat)

      class(type_stencil) :: this

      !-- Input variable
      integer,  intent(in) :: n_r_max

      !-- Output variable
      real(cp), intent(out) :: eye(n_r_max, n_r_max)    ! Identity matrix
      real(cp), intent(out) :: drMat(n_r_max, n_r_max)  ! First derivative
      real(cp), intent(out) :: d2rMat(n_r_max, n_r_max) ! Second derivative
      real(cp), intent(out) :: d3rMat(n_r_max, n_r_max) ! Third derivative

      !-- Local variables
      integer :: i, j, n_r

      !-- Set everything to zero
      do j=1,n_r_max
         do i=1,n_r_max
            drMat(i,j) =0.0_cp
            d2rMat(i,j)=0.0_cp
            d3rMat(i,j)=0.0_cp
            eye(i,j)   =0.0_cp
         end do
         eye(j,j)=1.0_cp
      end do

      !-- Bulk points for 1st and 2nd derivatives
      do n_r=this%order/2+1,n_r_max-this%order/2
         drMat(n_r,n_r-this%order/2:n_r+this%order/2) = this%dr(n_r,:)
         d2rMat(n_r,n_r-this%order/2:n_r+this%order/2)=this%ddr(n_r,:)
      end do

      !-- Bulk points for 3rd derivative
      do n_r=2+this%order/2,n_r_max-this%order/2-1
         d3rMat(n_r,n_r-this%order/2-1:n_r+this%order/2+1)=this%dddr(n_r,:)
      end do

      !-- Boundary points for 1st derivative
      do n_r=1,this%order/2
         drMat(n_r,1:this%order+1)                         =this%dr_top(n_r,:)
         drMat(n_r_max-n_r+1,n_r_max:n_r_max-this%order:-1)=this%dr_bot(n_r,:)
      end do

      !-- Boundary points for 2nd derivative
      do n_r=1,this%order/2
         d2rMat(n_r,1:this%order+2)                           =this%ddr_top(n_r,:)
         d2rMat(n_r_max-n_r+1,n_r_max:n_r_max-this%order-1:-1)=this%ddr_bot(n_r,:)
      end do

      !-- Boundary points for 3rd derivative
      do n_r=1,this%order/2+1
         d3rMat(n_r,1:this%order+3)                           =this%dddr_top(n_r,:)
         d3rMat(n_r_max-n_r+1,n_r_max:n_r_max-this%order-2:-1)=this%dddr_bot(n_r,:)
      end do

   end subroutine get_FD_matder
!---------------------------------------------------------------------------
   integer function factorial(n)
      !
      ! Compute the factorial of n
      !

      !-- Input variable
      integer, intent(inout) :: n

      integer :: i

      factorial = 1
      do i=1,n
         factorial = factorial*i
      end do

   end function factorial
!---------------------------------------------------------------------------
   subroutine inverse(a,c,n)
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      !

      !-- Input:
      integer,  intent(in) :: n
      real(cp), intent(inout) :: a(n,n)

      !-- Output
      real(cp), intent(out) :: c(n,n)

      !-- Local variables
      real(cp) :: L(n,n), U(n,n), b(n), d(n), x(n)
      real(cp) :: coeff
      integer :: i, j, k

      ! step 0: initialization for matrices L and U and b
      L(:,:)=0.0_cp
      U(:,:)=0.0_cp
      b(:)  =0.0_cp

      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = one
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
         do i=1,j
            U(i,j) = a(i,j)
         end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
         b(k)=one
         d(1) = b(1)
         ! Step 3a: Solve Ld=b using the forward substitution
         do i=2,n
            d(i)=b(i)
            do j=1,i-1
               d(i) = d(i) - L(i,j)*d(j)
            end do
         end do
         ! Step 3b: Solve Ux=d using the back substitution
         x(n)=d(n)/U(n,n)
         do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
               x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
         end do
         ! Step 3c: fill the solutions x(n) into column k of C
         do i=1,n
            c(i,k) = x(i)
         end do
         b(k)=0.0_cp
      end do

   end subroutine inverse
!---------------------------------------------------------------------------
end module finite_differences
