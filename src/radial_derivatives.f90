module radial_der
   !
   ! Radial derivatives functions
   !

   use constants, only: zero, one, three
   use precision_mod
   use cosine_transform_odd
   use finite_differences, only: type_stencil

   implicit none

   private

   interface get_dr
      module procedure get_dr_real_1d
      module procedure get_dr_complex
   end interface get_dr

   interface get_dr_fd
      module procedure get_dr_real_1d_fd
      module procedure get_dr_complex_fd
   end interface get_dr_fd

   interface get_dcheb
      module procedure get_dcheb_real_1d
      module procedure get_dcheb_complex
   end interface get_dcheb

   public :: get_dr, get_drNS, get_ddr, get_dddr, get_dcheb, &
   &         get_dr_fd, get_ddr_fd, get_dddr_fd

contains

!------------------------------------------------------------------------------
   subroutine get_dr_complex(f,df,n_f_max,n_f_start,n_f_stop, &
        &            n_r_max,n_cheb_max,work1,work2,chebt,drx)
      !
      !  Returns first radial derivative df of the input function f.      
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points . 
      !
    
      !-- Input variables:
      integer,           intent(in) :: n_r_max          ! number of radial grid points
      integer,           intent(in) :: n_f_max          ! first dim of f
      complex(cp),       intent(in) :: f(n_f_max,n_r_max)
      integer,           intent(in) :: n_f_start        ! first function to be treated
      integer,           intent(in) :: n_f_stop         ! last function to be treated
      integer,           intent(in) :: n_cheb_max       ! max number of cheb modes
      real(cp),          intent(in) :: drx(n_r_max)     ! first derivatives of x(r)
      type(costf_odd_t), intent(in) :: chebt
    
      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,n_r_max)   ! first derivative of f
      complex(cp), intent(out) :: work1(n_f_max,n_r_max)! work array needed for costf
      complex(cp), intent(out) :: work2(n_f_max,n_r_max)! work array for f transfer
    
      !-- Local:
      integer :: n_r,n_f
    
    
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do
    
      !-- Transform f to cheb space:
      call chebt%costf1(work2,n_f_max,n_f_start,n_f_stop,work1)
    
      !-- Get derivatives:
      call get_dcheb(work2,df,n_f_max,n_f_start,n_f_stop, n_r_max,n_cheb_max,one)
    
      !-- Transform back:
      call chebt%costf1(df,n_f_max,n_f_start,n_f_stop,work1)
    
      !-- New map:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r)=drx(n_r)*df(n_f,n_r)
         end do
      end do

   end subroutine get_dr_complex
!------------------------------------------------------------------------------
   subroutine get_dr_real_1d(f,df,n_r_max,n_cheb_max,work1,work2,  &
              &              chebt,drx)
 
      !-- Input variables:
      integer,           intent(in) :: n_r_max          ! number of radial grid points
      real(cp),          intent(in) :: f(n_r_max)
      integer,           intent(in) :: n_cheb_max       ! max number of cheb modes
      real(cp),          intent(in) :: drx(n_r_max)     ! first derivatives of x(r)
      type(costf_odd_t), intent(in) :: chebt
    
      !-- Output variables:
      real(cp), intent(out) :: df(n_r_max)   ! first derivative of f
      real(cp), intent(out) :: work1(n_r_max)! work array needed for costf
      real(cp), intent(out) :: work2(n_r_max)! work array for f transfer
    
      !-- Local:
      integer :: n_r
    
    
      !-- Copy input functions:
      do n_r=1,n_r_max
         work2(n_r)=f(n_r)
      end do
    
      !-- Transform f to cheb space:
      call chebt%costf1(work2,work1)
    
      !-- Get derivatives:
      call get_dcheb(work2,df,n_r_max,n_cheb_max,one)
    
      !-- Transform back:
      call chebt%costf1(df,work1)
    
      !-- New map:
      do n_r=1,n_r_max
         df(n_r)=drx(n_r)*df(n_r)
      end do

   end subroutine get_dr_real_1d
!------------------------------------------------------------------------------
   subroutine get_drNS(f,df,n_f_max,n_f_start,n_f_stop, &
        &              n_r_max,n_cheb_max,work1,        &
        &              chebt,drx)
      !
      !  Returns first radial derivative df of the input function f.      
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivatives of    
      !  the functions f(n_f_start,*) to f(n_f_stop,*) by transforming    
      !  to a Chebychev representation using n_r_max radial grid points . 
      !  Note: when using this function the input field f is slightly     
      !  changed by the back and forth transform. Use s_get_dr.f to       
      !  avoid this.                                                      
      !
    
      !-- Input variables:
      integer,           intent(in) :: n_r_max          ! number of radial grid points
      integer,           intent(in) :: n_f_max          ! first dim of f
      complex(cp),       intent(inout) :: f(n_f_max,n_r_max)
      integer,           intent(in) :: n_f_start        ! first function to be treated
      integer,           intent(in) :: n_f_stop         ! last function to be treated
      integer,           intent(in) :: n_cheb_max       ! max number of cheb modes
      real(cp),          intent(in) :: drx(*)           ! first derivatives of x(r)
      type(costf_odd_t), intent(in) :: chebt
    
      !-- Output variables:
      complex(cp), intent(out) :: work1(n_f_max,n_r_max) ! work array needed for costf
      complex(cp), intent(out) :: df(n_f_max,n_r_max)    ! first derivative of f
    
      !-- Local variables:
      integer :: n_r,n_f
    
      !-- Transform f to cheb space:
      call chebt%costf1(f,n_f_max,n_f_start,n_f_stop,work1)
    
      !-- Get derivatives:
      call get_dcheb(f,df,n_f_max,n_f_start,n_f_stop,n_r_max,n_cheb_max,one)
    
      !-- Transform back:
      call chebt%costf1(f,n_f_max,n_f_start,n_f_stop,work1)
      call chebt%costf1(df,n_f_max,n_f_start,n_f_stop,work1)
    
      !-- New map:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r)=drx(n_r)*df(n_f,n_r)
         end do
      end do

   end subroutine get_drNS
!------------------------------------------------------------------------------
   subroutine get_ddr(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
        &             n_r_max,n_cheb_max,work1,work2,      &
        &             chebt,drx,ddrx)
      !
      !  Returns first radial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivatives of    
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !
    
      !-- Input variables:
      integer,           intent(in) :: n_r_max         ! number of radial grid points
      integer,           intent(in) :: n_f_max         ! first dim of f
      complex(cp),       intent(in) :: f(n_f_max,n_r_max)
      integer,           intent(in) :: n_f_start       ! first function to be treated
      integer,           intent(in) :: n_f_stop        ! last function to be treated
      integer,           intent(in) :: n_cheb_max      ! number of cheb modes
      real(cp),          intent(in) :: drx(n_r_max)    ! first derivatives of x(r)
      real(cp),          intent(in) :: ddrx(n_r_max)   ! second derivatives of x(r)
      type(costf_odd_t), intent(in) :: chebt
    
      !-- Output variables:
      complex(cp), intent(out) :: work1(n_f_max,n_r_max)! work array needed for costf
      complex(cp), intent(out) :: work2(n_f_max,n_r_max)! work array for f transfer
      complex(cp), intent(out) :: df(n_f_max,n_r_max)   ! first derivative of f
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max)  ! second derivative of f
    
      !-- Local variables:
      integer :: n_r,n_f
    
      !-- Copy input functions:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do
    
      !-- Transform f to cheb space:
      call chebt%costf1(work2,n_f_max,n_f_start,n_f_stop,work1)
    
      !-- Get derivatives:
      call get_ddcheb(work2,df,ddf,n_f_max,n_f_start,n_f_stop, &
                      n_r_max,n_cheb_max,one)
    
      !-- Transform back:
      call chebt%costf1(df,n_f_max,n_f_start,n_f_stop,work1)
      call chebt%costf1(ddf,n_f_max,n_f_start,n_f_stop,work1)
    
      !-- New map:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            ddf(n_f,n_r)=   ddrx(n_r)*df(n_f,n_r) + drx(n_r)**2*ddf(n_f,n_r)
            df(n_f,n_r) =    drx(n_r)*df(n_f,n_r)
         end do
      end do

   end subroutine get_ddr
!------------------------------------------------------------------------------
   subroutine get_dddr(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop, &
                       n_r_max,n_cheb_max,work1,work2,           &
                       chebt,drx,ddrx,dddrx)
      !
      !  Returns first radial derivative df, the second radial deriv. ddf,
      !  and the third radial derivative dddf of the input function f.    
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !

      !-- Input variables:
      integer,           intent(in) :: n_r_max         ! number of radial grid points
      integer,           intent(in) :: n_f_max         ! first dim of f
      complex(cp),       intent(in) :: f(n_f_max,n_r_max)
      integer,           intent(in) :: n_f_start       ! first function to be treated
      integer,           intent(in) :: n_f_stop        ! last function to be treated
      integer,           intent(in) :: n_cheb_max      ! number of cheb_modes
      real(cp),          intent(in) :: drx(n_r_max)    ! first derivatives of x(r)
      real(cp),          intent(in) :: ddrx(n_r_max)   ! second derivatives of x(r)
      real(cp),          intent(in) :: dddrx(n_r_max)  ! third derivatives of x(r)
      type(costf_odd_t), intent(in) :: chebt
    

      !-- Work arrays:
      complex(cp), intent(out) :: work1(n_f_max,n_r_max) ! work array needed for costf
      complex(cp), intent(out):: work2(n_f_max,n_r_max)  ! work array needed for costf

      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,n_r_max)    ! first derivative of f
      complex(cp), intent(out) :: ddf(n_f_max,n_r_max)   ! second derivative of f
      complex(cp), intent(out) :: dddf(n_f_max,n_r_max)  ! third derivative of f

      !-- Local variables
      integer :: n_r,n_f

      !-- Copy input functions:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do

      !-- Transform f to cheb space:
      call chebt%costf1(work2,n_f_max,n_f_start,n_f_stop,work1)

      !-- Get derivatives:
      call get_dddcheb(work2,df,ddf,dddf,n_f_max,n_f_start,n_f_stop,  &
                       n_r_max,n_cheb_max,one)

      !-- Transform back:
      call chebt%costf1(df,n_f_max,n_f_start,n_f_stop,work1)
      call chebt%costf1(ddf,n_f_max,n_f_start,n_f_stop,work1)
      call chebt%costf1(dddf,n_f_max,n_f_start,n_f_stop,work1)

      !-- New map:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            dddf(n_f,n_r)=               dddrx(n_r)*df(n_f,n_r) + &
                          three*ddrx(n_r)*drx(n_r)*ddf(n_f,n_r) + &
                                      drx(n_r)**3*dddf(n_f,n_r)
            ddf(n_f,n_r) =ddrx(n_r)*df(n_f,n_r) + drx(n_r)**2*ddf(n_f,n_r)
            df(n_f,n_r)  = drx(n_r)*df(n_f,n_r)
         end do
      end do

   end subroutine get_dddr
!------------------------------------------------------------------------------
   subroutine get_dcheb_complex(f,df,n_f_max,n_f_start,n_f_stop, &
                                n_r_max,n_cheb_max,d_fac)
      !
      !  Returns chebychev coeffitiens of first derivative df and second  
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
      n_cheb  =n_cheb_max-1
      fac_cheb=d_fac*real(2*n_cheb,kind=cp)
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
      fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      df(n_cheb)=fac_cheb*f(n_cheb+1)

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         df(n_cheb)=df(n_cheb+2)+fac_cheb*f(n_cheb+1)
      end do

   end subroutine get_dcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddcheb(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
        &                n_r_max,n_cheb_max,d_fac)
      !
      !  Returns chebychev coefficents of first derivative df and second  
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
      n_cheb=n_cheb_max-1
      fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
         ddf(n_f,n_cheb)=zero
      end do
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2) + fac_cheb*f(n_f,n_cheb+1)
            ddf(n_f,n_cheb)=ddf(n_f,n_cheb+2) + fac_cheb*df(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_ddcheb
!------------------------------------------------------------------------------
   subroutine get_dddcheb(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop, &
                          n_r_max,n_cheb_max,d_fac)
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
            df(n_f,n_cheb)=zero
            ddf(n_f,n_cheb)=zero
            dddf(n_f,n_cheb)=zero
         end do
      end do
      n_cheb=n_cheb_max-1
      fac_cheb=d_fac*real(2*n_cheb,kind=cp)
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
         ddf(n_f,n_cheb)=zero
         dddf(n_f,n_cheb)=zero
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*real(2*n_cheb,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2)    + fac_cheb*f(n_f,n_cheb+1)
            ddf(n_f,n_cheb)=ddf(n_f,n_cheb+2)  + fac_cheb*df(n_f,n_cheb+1)
            dddf(n_f,n_cheb)=dddf(n_f,n_cheb+2)+ fac_cheb*ddf(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_dddcheb
!------------------------------------------------------------------------------
   subroutine get_dr_real_1d_fd(f,df,n_r_max,stencil)

      !-- Input variables
      integer,            intent(in) :: n_r_max   ! Number of radial grid points
      type(type_stencil), intent(in) :: stencil   ! Structure that contains stencils
      real(cp),           intent(in) :: f(n_r_max)! Input array

      !-- Output variable
      real(cp), intent(out) :: df(n_r_max)! Output array

      !-- Local variables
      integer :: n_r, od

      df(:) = 0.0_cp

      do od=0,stencil%order
         !-- Bulk points
         do n_r=1+stencil%order/2,n_r_max-stencil%order/2
            df(n_r) = df(n_r)+stencil%dr(n_r, od) * f(n_r-stencil%order/2+od)
         end do

         !-- Boundary points
         do n_r=1,stencil%order/2
            df(n_r) = df(n_r)+stencil%dr_top(n_r,od) * f(od+1)
         end do
         do n_r=1,stencil%order/2
            df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+stencil%dr_bot(n_r,od)*f(n_r_max-od)
         end do
      end do

   end subroutine get_dr_real_1d_fd
!---------------------------------------------------------------------------
   subroutine get_dr_complex_fd(f,df,n_f_max,n_f_start,n_f_stop,n_r_max, &
              &                 stencil)

      !-- Input variables
      integer,            intent(in) :: n_r_max           ! Number of radial grid points
      integer,            intent(in) :: n_f_max           ! first dimension of f
      integer,            intent(in) :: n_f_start         ! first function to be treated
      integer,            intent(in) :: n_f_stop          ! last function to be treated
      type(type_stencil), intent(in) :: stencil           ! Structure that contains stencils
      complex(cp),        intent(in) :: f(n_f_max,n_r_max)! Input array

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
      do od=0,stencil%order
         do n_r=1+stencil%order/2,n_r_max-stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)=df(n_f,n_r)+stencil%dr(n_r,od)*f(n_f,n_r-stencil%order/2+od)
            end do
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,stencil%order
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) = df(n_f,n_r)+stencil%dr_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
               &                       stencil%dr_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

   end subroutine get_dr_complex_fd
!---------------------------------------------------------------------------
   subroutine get_ddr_fd(f,df,ddf,n_f_max,n_f_start,n_f_stop,n_r_max,stencil)

      !-- Input variables
      integer,            intent(in) :: n_r_max           ! Number of radial grid points
      integer,            intent(in) :: n_f_max           ! first dimension of f
      integer,            intent(in) :: n_f_start         ! first function to be treated
      integer,            intent(in) :: n_f_stop          ! last function to be treated
      type(type_stencil), intent(in) :: stencil           ! Structure that contains stencils
      complex(cp),        intent(in) :: f(n_f_max,n_r_max)! Input array

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
      do od=0,stencil%order
         do n_r=1+stencil%order/2,n_r_max-stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)  = df(n_f,n_r) + stencil%dr(n_r,od) * f(n_f,n_r-stencil%order/2+od)
               ddf(n_f,n_r) = ddf(n_f,n_r)+stencil%ddr(n_r,od) * f(n_f,n_r-stencil%order/2+od)
            end do
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,stencil%order
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) = df(n_f,n_r)+stencil%dr_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
               &                       stencil%dr_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

      !-- Boundary points for 2nd derivative
      do od=0,stencil%order+1
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r) = ddf(n_f,n_r)+stencil%ddr_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
               &                       stencil%ddr_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

   end subroutine get_ddr_fd
!---------------------------------------------------------------------------
   subroutine get_dddr_fd(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop,n_r_max, &
              &           stencil)

      !-- Input variables
      integer,            intent(in) :: n_r_max           ! Number of radial grid points
      integer,            intent(in) :: n_f_max           ! first dimension of f
      integer,            intent(in) :: n_f_start         ! first function to be treated
      integer,            intent(in) :: n_f_stop          ! last function to be treated
      type(type_stencil), intent(in) :: stencil           ! Structure that contains stencils
      complex(cp),        intent(in) :: f(n_f_max,n_r_max)! Input array

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
      do od=0,stencil%order
         do n_r=1+stencil%order/2,n_r_max-stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r)  = df(n_f,n_r) + stencil%dr(n_r,od) * f(n_f,n_r-stencil%order/2+od)
               ddf(n_f,n_r) = ddf(n_f,n_r)+stencil%ddr(n_r,od) * f(n_f,n_r-stencil%order/2+od)
            end do
         end do
      end do

      !-- Bulk points for 3rd derivative
      do od=0,stencil%order+2
         do n_r=2+stencil%order/2,n_r_max-stencil%order/2-1
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r)=dddf(n_f,n_r)+stencil%dddr(n_r,od)*f(n_f,n_r-stencil%order/2-1+od)
            end do
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,stencil%order
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r) = df(n_f,n_r)+stencil%dr_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
               &                       stencil%dr_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

      !-- Boundary points for 2nd derivative
      do od=0,stencil%order+1
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r) = ddf(n_f,n_r)+stencil%ddr_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,stencil%order/2
            do n_f=n_f_start,n_f_stop
               ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
               &                       stencil%ddr_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

      !-- Boundary points for 3rd derivative
      do od=0,stencil%order+2
         do n_r=1,stencil%order/2+1
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r) = dddf(n_f,n_r)+stencil%dddr_top(n_r,od) * f(n_f,od+1)
            end do
         end do
         do n_r=1,stencil%order/2+1
            do n_f=n_f_start,n_f_stop
               dddf(n_f,n_r_max-n_r+1) = dddf(n_f,n_r_max-n_r+1)+               &
               &                       stencil%dddr_bot(n_r,od)*f(n_f,n_r_max-od)
            end do
         end do
      end do

   end subroutine get_dddr_fd
!---------------------------------------------------------------------------
end module radial_der
