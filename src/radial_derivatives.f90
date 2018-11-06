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
   &         finalize_der_arrays

   complex(cp), allocatable :: work(:,:)
   real(cp), allocatable :: work_1d_real(:)

contains

!------------------------------------------------------------------------------
   subroutine initialize_der_arrays(n_r_max,llm,ulm)
      !
      ! Allocate work arrays to compute derivatives
      !

      integer, intent(in) :: n_r_max
      integer, intent(in) :: llm
      integer, intent(in) :: ulm

      if ( .not. l_finite_diff ) then
         allocate( work_1d_real(n_r_max) )
         allocate( work(1:ulm-llm+1,n_r_max) )
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL+&
         &                 n_r_max*(ulm-llm+1)*SIZEOF_DEF_COMPLEX
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
               call r_scheme%costf1(work,n_f_max,n_f_start,n_f_stop)
            end if
          
            !-- Get derivatives:
            call get_dcheb(work,df,n_f_max,n_f_start,n_f_stop,n_r_max, &
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
            call r_scheme%costf1(work,n_f_max,n_f_start,n_f_stop)
         end if
    
         !-- Get derivatives:
         call get_ddcheb(work,df,ddf,n_f_max,n_f_start,n_f_stop, &
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
            call r_scheme%costf1(work,n_f_max,n_f_start,n_f_stop)
         end if

         !-- Get derivatives:
         call get_dddcheb(work,df,ddf,dddf,n_f_max,n_f_start,n_f_stop,  &
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
end module radial_der
