module radial_der
   !
   ! Radial derivatives functions
   !

   use constants, only: zero, one, three
   use precision_mod
   use mem_alloc
   use cosine_transform_odd
   use radial_scheme, only: type_rscheme
   use logic, only: l_finite_diff
   use parallel_mod
   use useful, only: abortRun

   implicit none

   private

   interface get_dr
      module procedure get_dr_real_1d
      module procedure get_dr_complex
   end interface get_dr

   interface get_dcheb
      module procedure get_dcheb_real_1d
      module procedure get_dcheb_complex
   end interface get_dcheb

   public :: get_ddr, get_dddr, get_dcheb, get_dr, initialize_der_arrays, &
   &         finalize_der_arrays, get_dr_Rloc, get_ddr_Rloc

   complex(cp), allocatable :: work(:,:)
   real(cp), allocatable :: work_1d_real(:)

contains

!------------------------------------------------------------------------------
   subroutine initialize_der_arrays(n_r_max,n_mlo_loc)
      !
      ! Allocate work arrays to compute derivatives
      !

      integer, intent(in) :: n_r_max
      integer, intent(in) :: n_mlo_loc

      if ( .not. l_finite_diff ) then
         allocate( work_1d_real(n_r_max) )
         allocate( work(n_mlo_loc,n_r_max) )
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL+&
         &                 n_r_max*n_mlo_loc*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_der_arrays
!------------------------------------------------------------------------------
   subroutine finalize_der_arrays
      !
      ! Deallocate work arrays
      !

      if ( .not. l_finite_diff ) deallocate( work_1d_real, work )

   end subroutine finalize_der_arrays
!------------------------------------------------------------------------------
   subroutine get_dcheb_complex(f,df,n_f_max,n_f_start,n_f_stop, &
              &                 n_r_max,n_cheb_max,d_fac)
      !
      !  Returns Chebyshev coeffitients of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_f_max,n_r_max).                             
      !

      !-- Input variables:
      integer,     intent(in) :: n_f_start  ! No of function to start with
      integer,     intent(in) :: n_f_stop   ! No of function to stop with
      integer,     intent(in) :: n_f_max    ! Max no of functions
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(n_f_max,n_r_max)
      real(cp),    intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      complex(cp), intent(out) ::  df(n_f_max,n_r_max)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(cp) :: fac_cheb


      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=zero
         end do
      end do

      !-- First Coefficient
      n_cheb  =n_cheb_max-1
      if ( n_r_max == n_cheb_max ) then
         fac_cheb=d_fac*real(n_cheb,kind=cp)
      else
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      end if
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2) + fac_cheb*f(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_dcheb_complex
!------------------------------------------------------------------------------
   subroutine get_dcheb_real_1d(f,df, n_r_max,n_cheb_max,d_fac)

      !-- Input variables:
      integer,  intent(in) :: n_r_max    ! Dimension of f,df,ddf
      integer,  intent(in) :: n_cheb_max ! Number of cheb modes
      real(cp), intent(in) :: f(n_r_max)
      real(cp), intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      real(cp), intent(out) ::  df(n_r_max)

      !-- Local variables:
      integer :: n_cheb
      real(cp) :: fac_cheb

      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         df(n_cheb)=0.0_cp
      end do
      n_cheb  =n_cheb_max-1

      !-- First coefficient
      if ( n_r_max == n_cheb_max ) then
         fac_cheb=d_fac*real(n_cheb,kind=cp)
      else
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      end if
      df(n_cheb)=fac_cheb*f(n_cheb+1)

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         df(n_cheb)=df(n_cheb+2)+fac_cheb*f(n_cheb+1)
      end do

   end subroutine get_dcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddcheb(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
              &          n_r_max,n_cheb_max,d_fac)
      !
      !  Returns Chebyshev coefficients of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_c_tot,n_r_max).                             
      !
    
      !-- Input variables:
      integer,     intent(in) :: n_f_start  ! No of column to start with
      integer,     intent(in) :: n_f_stop   ! No of column to stop with
      integer,     intent(in) :: n_f_max    ! First dimension of f,df,ddf
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(n_f_max,n_r_max)
      real(cp),    intent(in) :: d_fac      ! factor for interval mapping
    
      !-- Output variables:
      complex(cp), intent(out) ::  df(n_f_max,n_r_max)
      complex(cp), intent(out) ::  ddf(n_f_max,n_r_max)
    
      !-- local variables:
      integer :: n_f,n_cheb
      real(cp) :: fac_cheb
    
      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=zero
            ddf(n_f,n_cheb)=zero
         end do
      end do

      !-- First coefficients:
      n_cheb=n_cheb_max-1
      if ( n_cheb_max == n_r_max ) then
         fac_cheb=d_fac*real(n_cheb,kind=cp)
      else
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      end if
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb) =fac_cheb*f(n_f,n_cheb+1)
         ddf(n_f,n_cheb)=zero
      end do
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb) = df(n_f,n_cheb+2) + fac_cheb* f(n_f,n_cheb+1)
            ddf(n_f,n_cheb)=ddf(n_f,n_cheb+2) + fac_cheb*df(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_ddcheb
!------------------------------------------------------------------------------
   subroutine get_dddcheb(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop, &
              &           n_r_max,n_cheb_max,d_fac)
      !
      !  Returns chebychev coeffitiens of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_c_tot,n_r_max).                             
      !
         
      !-- Input variables:
      integer,     intent(in) :: n_f_start  ! No of column to start with
      integer,     intent(in) :: n_f_stop   ! No of column to stop with
      integer,     intent(in) :: n_f_max    ! First dimension of f,df,ddf
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(n_f_max,n_r_max)
      real(cp),    intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,n_r_max)
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max)
      complex(cp), intent(out) :: dddf(n_f_max,n_r_max)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(cp) :: fac_cheb

      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)  =zero
            ddf(n_f,n_cheb) =zero
            dddf(n_f,n_cheb)=zero
         end do
      end do

      !-- First coefficients
      n_cheb=n_cheb_max-1
      if ( n_cheb_max == n_r_max ) then
         fac_cheb=d_fac*real(n_cheb,kind=cp)
      else
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      end if
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)  =fac_cheb*f(n_f,n_cheb+1)
         ddf(n_f,n_cheb) =zero
         dddf(n_f,n_cheb)=zero
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)  =  df(n_f,n_cheb+2) + fac_cheb*  f(n_f,n_cheb+1)
            ddf(n_f,n_cheb) = ddf(n_f,n_cheb+2) + fac_cheb* df(n_f,n_cheb+1)
            dddf(n_f,n_cheb)=dddf(n_f,n_cheb+2) + fac_cheb*ddf(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_dddcheb
!---------------------------------------------------------------------------
   subroutine get_dr_real_1d(f,df,n_r_max,r_scheme)
 
      !-- Input variables:
      integer,             intent(in) :: n_r_max    ! number of radial grid points
      real(cp),            intent(in) :: f(n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
    
      !-- Output variables:
      real(cp), intent(out) :: df(n_r_max)   ! first derivative of f
    
      !-- Local:
      integer :: n_r, od
    
    
      if ( r_scheme%version == 'cheb' ) then

         !-- Copy input functions:
         do n_r=1,n_r_max
            work_1d_real(n_r)=f(n_r)
         end do
    
         !-- Transform f to cheb space:
         call r_scheme%costf1(work_1d_real)
    
         !-- Get derivatives:
         call get_dcheb(work_1d_real,df,n_r_max,r_scheme%n_max,one)
    
         !-- Transform back:
         call r_scheme%costf1(df)
    
         !-- New map:
         do n_r=1,n_r_max
            df(n_r)=r_scheme%drx(n_r)*df(n_r)
         end do

      else

         df(:) = 0.0_cp

         do od=0,r_scheme%order
            !-- Bulk points
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               df(n_r) = df(n_r)+r_scheme%dr(n_r, od) * f(n_r-r_scheme%order/2+od)
            end do
         end do

         do od=0,r_scheme%order_boundary
            !-- Boundary points
            do n_r=1,r_scheme%order/2
               df(n_r) = df(n_r)+r_scheme%dr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+r_scheme%dr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

      end if

   end subroutine get_dr_real_1d
!------------------------------------------------------------------------------
   subroutine get_dr_complex(f,df,n_f_max,n_f_start,n_f_stop, &
              &              n_r_max,r_scheme,nocopy,l_dct_in)
      !
      !  Returns first radial derivative df of the input function f.      
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points . 
      !
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max          ! number of radial grid points
      integer,             intent(in) :: n_f_max          ! first dim of f
      complex(cp),         intent(inout) :: f(n_f_max,n_r_max)
      integer,             intent(in) :: n_f_start        ! first function to be treated
      integer,             intent(in) :: n_f_stop         ! last function to be treated
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: nocopy
      logical, optional,   intent(in) :: l_dct_in
    
      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,n_r_max)   ! first derivative of f
    
      !-- Local:
      integer :: n_r,n_f,od
      logical :: copy_array, l_dct_in_loc
    
      if ( r_scheme%version == 'cheb' ) then

         if ( present(l_dct_in) ) then
            l_dct_in_loc=l_dct_in
         else
            l_dct_in_loc=.true.
         end if

         if ( present(nocopy) ) then
            copy_array=.false.
         else
            copy_array=.true.
         end if
    
         if ( copy_array )  then
            do n_r=1,n_r_max
               do n_f=n_f_start,n_f_stop
                  work(n_f,n_r)=f(n_f,n_r)
               end do
            end do
       
            !-- Transform f to cheb space:
            if ( l_dct_in_loc ) then
               call r_scheme%costf1(work(1:n_f_max,1:n_r_max),n_f_max,n_f_start,n_f_stop)
            end if
          
            !-- Get derivatives:
            call get_dcheb(work(1:n_f_max,1:n_r_max),df,n_f_max,n_f_start,n_f_stop,n_r_max, &
                 &         r_scheme%n_max,one)
          
            !-- Transform back:
            call r_scheme%costf1(df,n_f_max,n_f_start,n_f_stop)

         else

            !-- Transform f to cheb space:
            if ( l_dct_in_loc ) then
               call r_scheme%costf1(f,n_f_max,n_f_start,n_f_stop)
            end if
          
            !-- Get derivatives:
            call get_dcheb(f,df,n_f_max,n_f_start,n_f_stop,n_r_max, &
                 &         r_scheme%n_max,one)
          
            !-- Transform back:
            if ( l_dct_in_loc ) then
               call r_scheme%costf1(f,n_f_max,n_f_start,n_f_stop)
            end if
            call r_scheme%costf1(df,n_f_max,n_f_start,n_f_stop)

         end if
       
         !-- New map:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)=r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) =zero
            end do
         end do

         !-- Bulk points for 1st derivative
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r)=df(n_f,n_r)+r_scheme%dr(n_r,od)*f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_dr_complex
!------------------------------------------------------------------------------
   subroutine get_ddr(f,df,ddf,n_f_max,n_f_start,n_f_stop,n_r_max,r_scheme, &
              &       l_dct_in)
      !
      !  Returns first radial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivatives of    
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max       ! number of radial grid points
      integer,             intent(in) :: n_f_max       ! first dim of f
      complex(cp),         intent(in) :: f(n_f_max,n_r_max)
      integer,             intent(in) :: n_f_start     ! first function to be treated
      integer,             intent(in) :: n_f_stop      ! last function to be treated
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: l_dct_in
    
      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,n_r_max)   ! first derivative of f
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max)  ! second derivative of f
    
      !-- Local variables:
      integer :: n_r,n_f,od
      logical :: l_dct_in_loc

      if ( r_scheme%version == 'cheb' ) then

         if ( present(l_dct_in) ) then
            l_dct_in_loc=l_dct_in
         else
            l_dct_in_loc=.true.
         end if
    
         !-- Copy input functions:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               work(n_f,n_r)=f(n_f,n_r)
            end do
         end do
    
         !-- Transform f to cheb space:
         if ( l_dct_in_loc ) then
            call r_scheme%costf1(work(1:n_f_max,1:n_r_max),n_f_max,n_f_start,n_f_stop)
         end if
    
         !-- Get derivatives:
         call get_ddcheb(work(1:n_f_max,1:n_r_max),df,ddf,n_f_max,n_f_start,n_f_stop, &
              &          n_r_max,r_scheme%n_max,one)
    
         !-- Transform back:
         call r_scheme%costf1(df,n_f_max,n_f_start,n_f_stop)
         call r_scheme%costf1(ddf,n_f_max,n_f_start,n_f_stop)
    
         !-- New map:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r)=r_scheme%ddrx(n_r)*df(n_f,n_r)+r_scheme%drx(n_r)* &
               &             r_scheme%drx(n_r)*ddf(n_f,n_r)
               df(n_f,n_r) = r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) =zero
               ddf(n_f,n_r)=zero
            end do
         end do

         !-- Bulk points for 1st and 2nd derivatives
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r)  = df(n_f,n_r) + r_scheme%dr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

         !-- Boundary points for 2nd derivative
         do od=0,r_scheme%order_boundary+1
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%ddr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_ddr
!------------------------------------------------------------------------------
   subroutine get_dddr(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop, &
              &        n_r_max,r_scheme,l_dct_in)
      !
      !  Returns first radial derivative df, the second radial deriv. ddf,
      !  and the third radial derivative dddf of the input function f.    
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !

      !-- Input variables:
      integer,             intent(in) :: n_r_max         ! number of radial grid points
      integer,             intent(in) :: n_f_max         ! first dim of f
      complex(cp),         intent(in) :: f(n_f_max,n_r_max)
      integer,             intent(in) :: n_f_start       ! first function to be treated
      integer,             intent(in) :: n_f_stop        ! last function to be treated
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: l_dct_in
    
      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,n_r_max)    ! first derivative of f
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max)   ! second derivative of f
      complex(cp), intent(out) :: dddf(n_f_max,n_r_max)  ! third derivative of f

      !-- Local variables
      integer :: n_r,n_f,od
      logical :: l_dct_in_loc

      if ( r_scheme%version == 'cheb' ) then

         if ( present(l_dct_in) ) then
            l_dct_in_loc=l_dct_in
         else
            l_dct_in_loc=.true.
         end if

         !-- Copy input functions:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               work(n_f,n_r)=f(n_f,n_r)
            end do
         end do

         !-- Transform f to cheb space:
         if ( l_dct_in_loc ) then
            call r_scheme%costf1(work(1:n_f_max,1:n_r_max),n_f_max,n_f_start,n_f_stop)
         end if

         !-- Get derivatives:
         call get_dddcheb(work(1:n_f_max,1:n_r_max),df,ddf,dddf,n_f_max,n_f_start,n_f_stop,  &
              &           n_r_max,r_scheme%n_max,one)

         !-- Transform back:
         call r_scheme%costf1(df,n_f_max,n_f_start,n_f_stop)
         call r_scheme%costf1(ddf,n_f_max,n_f_start,n_f_stop)
         call r_scheme%costf1(dddf,n_f_max,n_f_start,n_f_stop)

         !-- New map:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r)=        r_scheme%dddrx(n_r)*df(n_f,n_r) +   &
               &             three*r_scheme%ddrx(n_r)*r_scheme%drx(n_r)* &
               &                                        ddf(n_f,n_r) +   &
               &             r_scheme%drx(n_r)*r_scheme%drx(n_r)*        &
               &             r_scheme%drx(n_r)*        dddf(n_f,n_r)
               ddf(n_f,n_r) =r_scheme%ddrx(n_r)*df(n_f,n_r) + r_scheme%drx(n_r)* &
               &             r_scheme%drx(n_r)*ddf(n_f,n_r)
               df(n_f,n_r)  = r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)  =zero
               ddf(n_f,n_r) =zero
               dddf(n_f,n_r)=zero
            end do
         end do

         !-- Bulk points for 1st and 2nd derivatives
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r)  = df(n_f,n_r) + r_scheme%dr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Bulk points for 3rd derivative
         do od=0,r_scheme%order+2
            do n_r=2+r_scheme%order/2,n_r_max-r_scheme%order/2-1
               do n_f=n_f_start,n_f_stop
                  dddf(n_f,n_r)=dddf(n_f,n_r)+r_scheme%dddr(n_r,od)*f(n_f,n_r-r_scheme%order/2-1+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

         !-- Boundary points for 2nd derivative
         do od=0,r_scheme%order_boundary+1
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=n_f_start,n_f_stop
                  ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%ddr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

         !-- Boundary points for 3rd derivative
         do od=0,r_scheme%order_boundary+2
            do n_r=1,r_scheme%order/2+1
               do n_f=n_f_start,n_f_stop
                  dddf(n_f,n_r) = dddf(n_f,n_r)+r_scheme%dddr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2+1
               do n_f=n_f_start,n_f_stop
                  dddf(n_f,n_r_max-n_r+1) = dddf(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dddr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_dddr
!------------------------------------------------------------------------------
   subroutine get_dr_Rloc(f_Rloc, df_Rloc, n_lm_loc, nRstart, nRstop, n_r_max, r_scheme)
      !
      ! Purpose of this subroutine is to take the first radial derivative of an input
      ! complex array distributed over radius. This can only be used with
      ! finite differences.
      !

      !-- Input variables
      integer,             intent(in) :: n_lm_loc, nRstart, nRstop, n_r_max
      class(type_rscheme), intent(in) :: r_scheme
      complex(cp),         intent(in) :: f_Rloc(n_lm_loc,nRstart:nRstop)

      !-- Output variable
      complex(cp), intent(out) ::  df_Rloc(n_lm_loc,nRstart:nRstop)

      !-- Local variables:
      complex(cp) :: work_ghost(n_lm_loc,nRstart-r_scheme%order/2:nRstop+r_scheme%order/2)
      complex(cp) :: ftop(n_lm_loc,r_scheme%order_boundary+1)
      complex(cp) :: fbot(n_lm_loc,n_r_max-r_scheme%order_boundary:n_r_max)
      integer :: n_r, od, start_lm, stop_lm

      if ( (r_scheme%order>2 .or. r_scheme%order_boundary>2) .and. &
       &   (nRstop-nRstart+1)<r_scheme%order ) then
         call abortRun('Distributed r-der not implemented in this case yet!')
      end if

      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=1; stop_lm=n_lm_loc
      call get_openmp_blocks(start_lm,stop_lm)

      !-- Copy input array
      work_ghost(start_lm:stop_lm,nRstart:nRstop)=f_Rloc(start_lm:stop_lm,:)
      do n_r=1,r_scheme%order_boundary+1
         if (n_r >= nRstart .and. n_r <= nRstop) then
            ftop(start_lm:stop_lm,n_r)=f_Rloc(start_lm:stop_lm,n_r)
         end if
      end do

      do n_r=n_r_max-r_scheme%order_boundary,n_r_max
         if (n_r >= nRstart .and. n_r <= nRstop) then
            fbot(start_lm:stop_lm,n_r)=f_Rloc(start_lm:stop_lm,n_r)
         end if
      end do

      !-- Exchange the ghost zones
      !$omp barrier
      !$omp master
      call exch_ghosts(work_ghost, n_lm_loc, nRstart, nRstop, r_scheme%order/2)
      !$omp end master
      !$omp barrier

      !-- Initialize to zero:
      do n_r=nRstart,nRstop
         df_Rloc(start_lm:stop_lm,n_r) =zero
      end do

      !-- Bulk points for 1st derivative
      do od=0,r_scheme%order
         do n_r=nRstart,nRstop
            df_Rloc(start_lm:stop_lm,n_r)=df_Rloc(start_lm:stop_lm,n_r) + &
            &                             r_scheme%dr(n_r,od)*            &
            &              work_ghost(start_lm:stop_lm,n_r-r_scheme%order/2+od)
         end do
      end do

      !-- Exchange boundary values
      !$omp barrier
      !$omp master
      call get_bound_vals(fbot, ftop, n_lm_loc, nRstart, nRstop, n_r_max, &
           &              r_scheme%order_boundary+1)
      !$omp end master
      !$omp barrier

      !-- Boundary points for 1st derivative
      if ( coord_r == 0 ) then
         df_Rloc(start_lm:stop_lm,1)=zero
         do od=0,r_scheme%order_boundary
            df_Rloc(start_lm:stop_lm,1)=df_Rloc(start_lm:stop_lm,1) + &
            &            r_scheme%dr_top(1,od)*ftop(start_lm:stop_lm,od+1)
         end do
      end if

      if ( coord_r == n_ranks_r -1 ) then
         df_Rloc(start_lm:stop_lm,n_r_max)=zero
         do od=0,r_scheme%order_boundary
            df_Rloc(start_lm:stop_lm,n_r_max)=df_Rloc(start_lm:stop_lm,n_r_max)+  &
            &           r_scheme%dr_bot(1,od)*fbot(start_lm:stop_lm,n_r_max-od)
         end do
      end if

      !$omp end parallel

   end subroutine get_dr_Rloc
!------------------------------------------------------------------------------
   subroutine get_ddr_Rloc(f_Rloc, df_Rloc, ddf_Rloc, n_lm_loc, nRstart, nRstop, &
              &            n_r_max, r_scheme)
      !
      ! Purpose of this subroutine is to take the first and second
      ! radial derivatives of an input complex array distributed over radius.
      ! This can only be used with finite differences.
      !

      !-- Input variables
      integer,             intent(in) :: n_lm_loc, nRstart, nRstop, n_r_max
      class(type_rscheme), intent(in) :: r_scheme
      complex(cp),         intent(in) :: f_Rloc(n_lm_loc,nRstart:nRstop)

      !-- Output variable
      complex(cp), intent(out) ::  df_Rloc(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) ::  ddf_Rloc(n_lm_loc,nRstart:nRstop)

      !-- Local variables:
      complex(cp) :: work_ghost(n_lm_loc,nRstart-r_scheme%order/2:nRstop+r_scheme%order/2)
      complex(cp) :: ftop(n_lm_loc,r_scheme%order_boundary+2)
      complex(cp) :: fbot(n_lm_loc,n_r_max-r_scheme%order_boundary-1:n_r_max)
      integer :: n_r, od, start_lm, stop_lm

      if ( (r_scheme%order>2 .or. r_scheme%order_boundary>2) .and. &
       &   (nRstop-nRstart+1)<r_scheme%order ) then
         call abortRun('Distributed r-der not implemented in this case yet!')
      end if

      !$omp parallel default(shared) private(start_lm,stop_lm,n_r,od)
      start_lm=1; stop_lm=n_lm_loc
      call get_openmp_blocks(start_lm,stop_lm)

      !-- Copy input array
      work_ghost(start_lm:stop_lm,nRstart:nRstop)=f_Rloc(start_lm:stop_lm,:)
      do n_r=1,r_scheme%order_boundary+2
         if (n_r >= nRstart .and. n_r <= nRstop) then
            ftop(start_lm:stop_lm,n_r)=f_Rloc(start_lm:stop_lm,n_r)
         end if
      end do

      do n_r=n_r_max-r_scheme%order_boundary-1,n_r_max
         if (n_r >= nRstart .and. n_r <= nRstop) then
            fbot(start_lm:stop_lm,n_r)=f_Rloc(start_lm:stop_lm,n_r)
         end if
      end do

      !-- Exchange the ghost zones
      !$omp barrier
      !$omp master
      call exch_ghosts(work_ghost, n_lm_loc, nRstart, nRstop, r_scheme%order/2)
      !$omp end master
      !$omp barrier

      !-- Initialize to zero:
      do n_r=nRstart,nRstop
         df_Rloc(start_lm:stop_lm,n_r) =zero
         ddf_Rloc(start_lm:stop_lm,n_r)=zero
      end do

      !-- Bulk points for 1st and 2nd derivatives
      do od=0,r_scheme%order
         do n_r=nRstart,nRstop
            df_Rloc(start_lm:stop_lm,n_r)=df_Rloc(start_lm:stop_lm,n_r)+   &
            &                             r_scheme%dr(n_r,od)*             &
            &              work_ghost(start_lm:stop_lm,n_r-r_scheme%order/2+od)
            ddf_Rloc(start_lm:stop_lm,n_r)=ddf_Rloc(start_lm:stop_lm,n_r)+ &
            &                            r_scheme%ddr(n_r,od)*             &
            &              work_ghost(start_lm:stop_lm,n_r-r_scheme%order/2+od)
         end do
      end do

      !-- Exchange boundary values
      !$omp barrier
      !$omp master
      call get_bound_vals(fbot, ftop, n_lm_loc, nRstart, nRstop, n_r_max, &
           &              r_scheme%order_boundary+2)
      !$omp end master
      !$omp barrier

      !-- Boundary points
      if ( coord_r == 0 ) then
         df_Rloc(start_lm:stop_lm,1) =zero
         ddf_Rloc(start_lm:stop_lm,1)=zero
         do od=0,r_scheme%order_boundary
            df_Rloc(start_lm:stop_lm,1)=df_Rloc(start_lm:stop_lm,1) + &
            &                 r_scheme%dr_top(1,od)*ftop(start_lm:stop_lm,od+1)
         end do
         do od=0,r_scheme%order_boundary+1
            ddf_Rloc(start_lm:stop_lm,1) = ddf_Rloc(start_lm:stop_lm,1) + &
            &                 r_scheme%ddr_top(1,od)*ftop(start_lm:stop_lm,od+1)
         end do
      end if

      if ( coord_r == n_ranks_r-1 ) then
         df_Rloc(start_lm:stop_lm,n_r_max) =zero
         ddf_Rloc(start_lm:stop_lm,n_r_max)=zero
         do od=0,r_scheme%order_boundary
            df_Rloc(start_lm:stop_lm,n_r_max)=df_Rloc(start_lm:stop_lm,n_r_max)+  &
            &             r_scheme%dr_bot(1,od)*fbot(start_lm:stop_lm,n_r_max-od)
         end do
         do od=0,r_scheme%order_boundary+1
            ddf_Rloc(start_lm:stop_lm,n_r_max)=ddf_Rloc(start_lm:stop_lm,n_r_max)+  &
            &             r_scheme%ddr_bot(1,od)*fbot(start_lm:stop_lm,n_r_max-od)
         end do
      end if

      !$omp end parallel

   end subroutine get_ddr_Rloc
!------------------------------------------------------------------------------
   subroutine exch_ghosts(f, n_lm_loc, nRstart, nRstop, nghosts)

      integer, intent(in) :: n_lm_loc, nRstart, nRstop, nghosts

      complex(cp), intent(inout) :: f(n_lm_loc, nRstart-nghosts:nRstop+nghosts)

#ifdef WITH_MPI
      integer :: n_counts, rightProc, leftProc, st(MPI_STATUS_SIZE)
      !integer :: tagsend, tagrecv, req(4)

      !req(:)=MPI_REQUEST_NULL
      !tagsend = coord_r
      !tagrecv = MPI_ANY_TAG
      !n_counts = n_lm_loc*nghosts
      !if ( coord_r < n_ranks_r-1 ) then
      !   call MPI_ISend(f(:,nRstop-nghosts+1:nRstop), n_counts, MPI_DEF_COMPLEX, &
      !        &         coord_r+1, tagsend, comm_r, req(1), ierr)
      !   call MPI_IRecv(f(:,nRstop+1:nRstop+nghosts), n_counts, MPI_DEF_COMPLEX, &
      !        &         coord_r+1, tagrecv, comm_r, req(2), ierr)
      !end if
      !if ( coord_r > 0 ) then
      !   call MPI_ISend(f(:,nRstart:nRstart+nghosts-1), n_counts, MPI_DEF_COMPLEX, &
      !        &         coord_r-1, tagsend, comm_r, req(3), ierr)
      !   call MPI_IRecv(f(:,nRstart-nghosts:nRstart-1), n_counts, MPI_DEF_COMPLEX, &
      !        &         coord_r-1, tagrecv, comm_r, req(4), ierr)
      !end if

      !call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)

      leftProc=coord_r-1
      if ( leftProc < 0 ) leftProc = MPI_PROC_NULL
      rightProc=coord_r+1
      if ( rightProc >= n_ranks_r ) rightProc = MPI_PROC_NULL
      n_counts = n_lm_loc*nghosts

      !-- Halo exchange in forward direction
      call MPI_SendRecv(f(:,nRstop-nghosts+1:nRstop), n_counts, MPI_DEF_COMPLEX,    &
           &            rightProc, 12345, f(:,nRstart-nghosts:nRstart-1), n_counts, &
           &            MPI_DEF_COMPLEX, leftProc, 12345, comm_r, st, ierr)

      !-- Halo exchange in backward direction
      call MPI_SendRecv(f(:,nRstart:nRstart+nghosts-1), n_counts, MPI_DEF_COMPLEX, &
           &            leftProc, 12345, f(:,nRstop+1:nRstop+nghosts), n_counts,   &
           &            MPI_DEF_COMPLEX, rightProc, 12345, comm_r, st, ierr)
#endif

   end subroutine exch_ghosts
!------------------------------------------------------------------------------
   subroutine get_bound_vals(fbot, ftop, n_lm_loc, nRstart, nRstop, n_r_max, nbounds)

      !-- Input variables
      integer, intent(in) :: n_lm_loc, nRstart, nRstop
      integer, intent(in) :: nbounds, n_r_max

      !-- Output boundary values
      complex(cp), intent(inout) :: ftop(n_lm_loc,nbounds)
      complex(cp), intent(inout) :: fbot(n_lm_loc,n_r_max-nbounds+1:n_r_max)

#ifdef WITH_MPI
      !-- Local variables
      integer :: nreq, nR, tag, req(2*nbounds)

      if ( nRstart > nbounds .and. nRstop < n_r_max-nbounds ) return

      req(:)=MPI_REQUEST_NULL
      nreq=0
      do nR=1,nbounds
         tag = 754432+nR
         if ( coord_r /= 0 .and. nRstart<=nR .and. nRstop>=nR) then
            call MPI_Send(ftop(:,nR), n_lm_loc, MPI_DEF_COMPLEX, 0, tag, &
                 &        comm_r, ierr)
         else if ( coord_r == 0 .and. nRstop < nR ) then
            nreq=nreq+1
            call MPI_IRecv(ftop(:,nR), n_lm_loc, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, tag, &
                 &         comm_r, req(nreq), ierr)
         end if
      end do
      call MPI_Waitall(nreq, req(1:nreq), MPI_STATUSES_IGNORE, ierr)

      req(:)=MPI_REQUEST_NULL
      nreq=0
      do nR=n_r_max-nbounds+1,n_r_max
         tag = 92113+nR
         if ( coord_r /= n_ranks_r-1 .and. nRstart<=nR .and. nRstop>=nR) then
            call MPI_Send(fbot(:,nR), n_lm_loc, MPI_DEF_COMPLEX, n_ranks_r-1, tag, &
                 &        comm_r, ierr)
         else if ( coord_r == n_ranks_r-1 .and. nR < nRstart ) then
            nreq=nreq+1
            call MPI_IRecv(fbot(:,nR), n_lm_loc, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, tag, &
                 &         comm_r, req(nreq), ierr)
         end if
      end do
      call MPI_Waitall(nreq, req(1:nreq), MPI_STATUSES_IGNORE, ierr)
#endif

   end subroutine get_bound_vals
!------------------------------------------------------------------------------
end module radial_der
