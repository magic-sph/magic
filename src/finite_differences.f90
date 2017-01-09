module finite_differences
   !
   ! This module is used to calculate the radial grid when finite differences
   ! are used
   !

   use precision_mod
   use parallel_mod, only: rank
   use constants, only: one, two, zero
   use truncation, only: n_r_max
   use useful, only: logWrite
   use mem_alloc, only: bytes_allocated


   implicit none

   private

   real(cp), allocatable :: dr_stencil(:,:)
   real(cp), allocatable :: ddr_stencil(:,:)
   real(cp), allocatable :: dddr_stencil(:,:)
   real(cp), allocatable :: dr_stencil_top(:,:)
   real(cp), allocatable :: dr_stencil_bot(:,:)
   real(cp), allocatable :: ddr_stencil_top(:,:)
   real(cp), allocatable :: ddr_stencil_bot(:,:)
   real(cp), allocatable :: dddr_stencil_top(:,:)
   real(cp), allocatable :: dddr_stencil_bot(:,:)

   public :: get_FD_grid, initialize_FD_arrays, finalize_FD_arrays

contains

   subroutine initialize_FD_arrays(n_r_max,order)
      !
      ! This subroutine allocates the arrays used when finite difference are used
      !

      integer, intent(in) :: n_r_max ! Number of radial grid points
      integer, intent(in) :: order   ! FD order

      allocate( dr_stencil(n_r_max,0:order) )
      allocate( ddr_stencil(n_r_max,0:order) )
      allocate( dddr_stencil(n_r_max,0:order+2) )
      allocate( dr_stencil_top(order/2,0:order) )
      allocate( dr_stencil_bot(order/2,0:order) )
      allocate( ddr_stencil_top(order/2,0:order+1) )
      allocate( ddr_stencil_bot(order/2,0:order+1) )
      allocate( dddr_stencil_top(order/2+1,0:order+2) )
      allocate( dddr_stencil_bot(order/2+1,0:order+2) )

      bytes_allocated=bytes_allocated+(n_r_max*(3*order+5)+         &
      &               order/2*(4*order+6)+(order/2+1)*(2*order+6))* &
      &               SIZEOF_DEF_REAL

   end subroutine initialize_FD_arrays
!---------------------------------------------------------------------------
   subroutine finalize_FD_arrays
      !
      ! This subroutine deallocates the arrays used in FD
      !

      deallocate( dr_stencil, ddr_stencil, dddr_stencil )
      deallocate( dr_stencil_top, dr_stencil_bot )
      deallocate( ddr_stencil_top, ddr_stencil_bot )
      deallocate( dddr_stencil_top, dddr_stencil_bot )

   end subroutine finalize_FD_arrays
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

      ! To be removed:
      real(cp) :: f(n_r_max), df(n_r_max)
      complex(cp) :: f_c(1,n_r_max), df_c(1,n_r_max), ddf_c(1,n_r_max), dddf_c(1,n_r_max)
      character(len=72) :: fileName
      integer :: fileHandle

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

      f = sin(r)
      do n_r=1,n_r_max
         f_c(1,n_r)=sin(r(n_r))
      end do
      call get_FD_coeffs(2, r)
      ! call get_dr_real_1d_fd(f, df, n_r_max, 2)
      call get_dr_complex_fd(f_c, df_c, 1, 1, 1, n_r_max, 2)
      call get_ddr_fd(f_c, df_c, ddf_c, 1, 1, 1, n_r_max, 2)
      call get_dddr_fd(f_c, df_c, ddf_c, dddf_c, 1, 1, 1, n_r_max, 2)

      open(newunit=fileHandle, file='test.txt', status='unknown')

      do n_r=1,n_r_max
         write(fileHandle, '(5ES20.12)') r(n_r), real(f_c(1,n_r)), &
         &                               real(df_c(1,n_r)), real(ddf_c(1,n_r)), &
         &                               real(dddf_c(1,n_r))
      end do

      close(fileHandle)
      stop

   end subroutine get_FD_grid
!---------------------------------------------------------------------------
   subroutine get_FD_coeffs(order, r)

      !-- Input quantities:
      integer, intent(in) :: order ! Order of the finite difference scheme
      real(cp), intent(in) :: r(:) ! Radius

      !-- Local quantities:
      real(cp), allocatable :: dr_spacing(:)
      real(cp), allocatable :: taylor_exp(:,:)
      real(cp), allocatable :: taylor_exp_inv(:,:)
      real(cp) :: weight
      integer :: n_r, od, od_in

      allocate( dr_spacing(order+1) )
      allocate( taylor_exp(0:order,0:order) )
      allocate( taylor_exp_inv(0:order,0:order) )

      !
      !-- Step 1: First and 2nd derivatives in the bulk
      !
      do n_r=1+order/2,n_r_max-order/2
         do od=0,order
            dr_spacing(od+1)=r(n_r-order/2+od)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+1)

         do od=0,order
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+1)

         do od_in=0,order
            !-- Preconditioning
            do od=0,order
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            dr_stencil(n_r,od_in) =taylor_exp_inv(1,od_in)
            ddr_stencil(n_r,od_in)=taylor_exp_inv(2,od_in)
         end do
      end do

      !
      !-- Step 2: First derivative for the outer points
      !
      do n_r=1,order/2
         do od=0,order
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+1)

         do od=0,order
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+1)

         do od_in=0,order
            !-- Preconditioning
            do od=0,order
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            dr_stencil_top(n_r,od_in) =taylor_exp_inv(1,od_in)
         end do
      end do

      !
      !-- Step 3: First derivative for the inner points
      !
      do n_r=1,order/2
         do od=0,order
            dr_spacing(od+1)=r(n_r_max-od)-r(n_r_max-n_r+1)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+1)

         do od=0,order
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+1)

         do od_in=0,order
            !-- Preconditioning
            do od=0,order
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            dr_stencil_bot(n_r,od_in) =taylor_exp_inv(1,od_in)
         end do
      end do

      deallocate( dr_spacing, taylor_exp, taylor_exp_inv )

      !
      !-- Step 4: 2nd derivative for the outer points
      !
      allocate( dr_spacing(order+2) )
      allocate( taylor_exp(0:order+1,0:order+1) )
      allocate( taylor_exp_inv(0:order+1,0:order+1) )

      do n_r=1,order/2
         do od=0,order+1
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+2)

         do od=0,order+1
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+2)

         do od_in=0,order+1
            !-- Preconditioning
            do od=0,order+1
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            ddr_stencil_top(n_r,od_in) =taylor_exp_inv(2,od_in)
         end do
      end do

      !
      !-- Step 5: 2nd derivative for the inner points
      !
      do n_r=1,order/2
         do od=0,order+1
            dr_spacing(od+1)=r(n_r_max-od)-r(n_r_max-n_r+1)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+2)

         do od=0,order+1
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+2)

         do od_in=0,order+1
            !-- Preconditioning
            do od=0,order+1
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            ddr_stencil_bot(n_r,od_in) =taylor_exp_inv(2,od_in)
         end do
      end do

      deallocate( dr_spacing, taylor_exp, taylor_exp_inv )

      !
      !-- Step 6: 3rd derivative in the bulk
      !
      allocate( dr_spacing(order+3) )
      allocate( taylor_exp(0:order+2,0:order+2) )
      allocate( taylor_exp_inv(0:order+2,0:order+2) )

      do n_r=2+order/2,n_r_max-order/2-1
         do od=0,order+2
            dr_spacing(od+1)=r(n_r-order/2-1+od)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+3)

         do od=0,order+2
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+3)

         do od_in=0,order+2
            !-- Preconditioning
            do od=0,order+2
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            dddr_stencil(n_r,od_in) =taylor_exp_inv(3,od_in)
         end do
      end do

      !
      !-- Step 7: 3rd derivative for the outer points
      !
      do n_r=1,order/2+1
         do od=0,order+2
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+3)

         do od=0,order+2
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+3)

         do od_in=0,order+2
            !-- Preconditioning
            do od=0,order+2
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            dddr_stencil_top(n_r,od_in) =taylor_exp_inv(3,od_in)
         end do
      end do

      !
      !-- Step 8: 3rd derivative for the inner points
      !
      do n_r=1,order/2+1
         do od=0,order+2
            dr_spacing(od+1)=r(n_r_max-od)-r(n_r_max-n_r+1)
         end do

         !-- This is a weight for matrix preconditioning
         weight = sum(abs(dr_spacing))/(order+3)

         do od=0,order+2
            taylor_exp(:, od) = (dr_spacing(:)/weight)**od
         end do

         !-- Invert the matrix to get the FD coeffs
         call inverse(taylor_exp, taylor_exp_inv, order+3)

         do od_in=0,order+2
            !-- Preconditioning
            do od=0,order+2
               taylor_exp_inv(od, od_in) = taylor_exp_inv(od,od_in)* &
               &                           factorial(od)/weight**od
            end do
            dddr_stencil_bot(n_r,od_in) =taylor_exp_inv(3,od_in)
         end do
      end do

      deallocate( dr_spacing, taylor_exp, taylor_exp_inv )

   end subroutine get_FD_coeffs
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
   subroutine get_dr_real_1d_fd(f,df,n_r_max, order)

      !-- Input variables
      integer,  intent(in) :: n_r_max     ! Number of radial grid points
      integer,  intent(in) :: order       ! Order of the FD scheme
      real(cp), intent(in) :: f(n_r_max)  ! Input array

      !-- Output variable
      real(cp), intent(out) :: df(n_r_max)! Output array

      !-- Local variables
      integer :: n_r, od

      df(:) = 0.0_cp

      print*, dr_stencil_bot

      do od=0,order
         !-- Bulk points
         do n_r=1+order/2,n_r_max-order/2
            df(n_r) = df(n_r)+dr_stencil(n_r, od) * f(n_r-order/2+od)
         end do

         !-- Boundary points
         do n_r=1,order/2
            df(n_r) = df(n_r)+dr_stencil_top(n_r,od) * f(od+1)
         end do
         do n_r=1,order/2
            df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+dr_stencil_bot(n_r,od)*f(n_r_max-od)
         end do
      end do

   end subroutine get_dr_real_1d_fd
!---------------------------------------------------------------------------
   subroutine get_dr_complex_fd(f,df,n_f_max,n_f_start,n_f_stop,n_r_max,order)

      !-- Input variables
      integer,     intent(in) :: n_r_max            ! Number of radial grid points
      integer,     intent(in) :: n_f_max            ! first dimension of f
      integer,     intent(in) :: n_f_start          ! first function to be treated
      integer,     intent(in) :: n_f_stop           ! last function to be treated
      integer,     intent(in) :: order              ! Order of the FD scheme
      complex(cp), intent(in) :: f(n_f_max,n_r_max) ! Input array

      !-- Output variable
      complex(cp), intent(out) :: df(n_f_max,n_r_max) ! First derivative

      !-- Local variables
      integer :: n_r, n_f, od

      !-- Initiatise to zero:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r) =zero
         end do
      end do

      !-- Bulk points for 1st derivative
      do od=0,order
         do n_r=1+order/2,n_r_max-order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)=df(n_f,n_r)+dr_stencil(n_r,od)*f(n_f,n_r-order/2+od)
            end do
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,order
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) = df(n_f,n_r)+dr_stencil_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
               &                       dr_stencil_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

   end subroutine get_dr_complex_fd
!---------------------------------------------------------------------------
   subroutine get_ddr_fd(f,df,ddf,n_f_max,n_f_start,n_f_stop,n_r_max, order)

      !-- Input variables
      integer,     intent(in) :: n_r_max            ! Number of radial grid points
      integer,     intent(in) :: n_f_max            ! first dimension of f
      integer,     intent(in) :: n_f_start          ! first function to be treated
      integer,     intent(in) :: n_f_stop           ! last function to be treated
      integer,     intent(in) :: order              ! Order of the FD scheme
      complex(cp), intent(in) :: f(n_f_max,n_r_max) ! Input array

      !-- Output variable
      complex(cp), intent(out) :: df(n_f_max,n_r_max) ! First derivative
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max)! Second derivative

      !-- Local variables
      integer :: n_r, n_f, od

      !-- Initiatise to zero:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r) =zero
            ddf(n_f,n_r)=zero
         end do
      end do

      !-- Bulk points for 1st and 2nd derivatives
      do od=0,order
         do n_r=1+order/2,n_r_max-order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)  = df(n_f,n_r) + dr_stencil(n_r,od) * f(n_f,n_r-order/2+od)
               ddf(n_f,n_r) = ddf(n_f,n_r)+ddr_stencil(n_r,od) * f(n_f,n_r-order/2+od)
            end do
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,order
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) = df(n_f,n_r)+dr_stencil_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
               &                       dr_stencil_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

      !-- Boundary points for 2nd derivative
      do od=0,order+1
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r) = ddf(n_f,n_r)+ddr_stencil_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
               &                       ddr_stencil_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

   end subroutine get_ddr_fd
!---------------------------------------------------------------------------
   subroutine get_dddr_fd(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop,n_r_max,order)

      !-- Input variables
      integer,     intent(in) :: n_r_max            ! Number of radial grid points
      integer,     intent(in) :: n_f_max            ! first dimension of f
      integer,     intent(in) :: n_f_start          ! first function to be treated
      integer,     intent(in) :: n_f_stop           ! last function to be treated
      integer,     intent(in) :: order              ! Order of the FD scheme
      complex(cp), intent(in) :: f(n_f_max,n_r_max) ! Input array

      !-- Output variable
      complex(cp), intent(out) :: df(n_f_max,n_r_max)  ! First derivative
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max) ! Second derivative
      complex(cp), intent(out) :: dddf(n_f_max,n_r_max)! Third derivative

      !-- Local variables
      integer :: n_r, n_f, od

      !-- Initiatise to zero:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r)  =zero
            ddf(n_f,n_r) =zero
            dddf(n_f,n_r)=zero
         end do
      end do

      !-- Bulk points for 1st and 2nd derivatives
      do od=0,order
         do n_r=1+order/2,n_r_max-order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)  = df(n_f,n_r) + dr_stencil(n_r,od) * f(n_f,n_r-order/2+od)
               ddf(n_f,n_r) = ddf(n_f,n_r)+ddr_stencil(n_r,od) * f(n_f,n_r-order/2+od)
            end do
         end do
      end do

      print*, dddr_stencil(16, :)
      !-- Bulk points for 3rd derivative
      do od=0,order+2
         do n_r=2+order/2,n_r_max-order/2-1
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r)=dddf(n_f,n_r)+dddr_stencil(n_r,od)*f(n_f,n_r-order/2-1+od)
            end do
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,order
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) = df(n_f,n_r)+dr_stencil_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
               &                       dr_stencil_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

      !-- Boundary points for 2nd derivative
      do od=0,order+1
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r) = ddf(n_f,n_r)+ddr_stencil_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
               &                       ddr_stencil_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

      !-- Boundary points for 3rd derivative
      do od=0,order+2
         do n_r=1,order/2+1
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r) = dddf(n_f,n_r)+dddr_stencil_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,order/2+1
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r_max-n_r+1) = dddf(n_f,n_r_max-n_r+1)+               &
               &                       dddr_stencil_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

   end subroutine get_dddr_fd
!---------------------------------------------------------------------------
end module finite_differences
