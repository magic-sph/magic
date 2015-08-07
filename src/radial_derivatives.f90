!$Id$
module radial_der

   use cosine_transform, only: costf1

   implicit none

   private

   interface get_dr
      module procedure get_dr_real_1d
      module procedure get_dr_complex
   end interface

   interface get_dcheb
      module procedure get_dcheb_real_1d
      module procedure get_dcheb_complex
   end interface

   public :: get_dr, get_drNS, get_ddr, get_dddr, get_dcheb

contains

!------------------------------------------------------------------------------
   subroutine get_dr_complex(f,df,n_f_max,n_f_start,n_f_stop, &
        &            n_r_max,n_cheb_max,work1,work2,          &
        &            i_costf_init,d_costf_init,drx)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns first radial derivative df of the input function f.      |
      !  |  Array f(n_f_max,*) may contain several functions numbered by     |
      !  |  the first index. The subroutine calculates the derivaties of     |
      !  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
      !  |  to a Chebychev representation using n_r_max radial grid points . |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
    
      !-- Input variables:
      integer,         intent(in) :: n_f_max          ! first dim of f
      complex(kind=8), intent(in) :: f(n_f_max,*)
      integer,         intent(in) :: n_f_start        ! first function to be treated
      integer,         intent(in) :: n_f_stop         ! last function to be treated
      integer,         intent(in) :: n_r_max          ! number of radial grid points
      integer,         intent(in) :: n_cheb_max       ! max number of cheb modes
      real(kind=8),    intent(in) :: drx(*)           ! first derivatives of x(r)
      integer,         intent(in) :: i_costf_init(*)  ! info for costf
      real(kind=8),    intent(in) :: d_costf_init(*)  ! info for costf
    
      !-- Output variables:
      complex(kind=8), intent(out) :: df(n_f_max,*)       ! first derivative of f
      complex(kind=8), intent(out) :: work1(n_f_max,*)    ! work array needed for costf
      complex(kind=8), intent(out) :: work2(n_f_max,n_r_max) ! work array for f transfer
    
      !-- Local:
      integer :: n_r,n_f
    
    
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do
    
      !-- Transform f to cheb space:
      call costf1(work2,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
    
      !-- Get derivatives:
      call get_dcheb(work2,df,n_f_max,n_f_start,n_f_stop, n_r_max,n_cheb_max,1.D0)
    
      !-- Transform back:
      call costf1(df,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
    
      !-- New map:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r)=drx(n_r)*df(n_f,n_r)
         end do
      end do

   end subroutine get_dr_complex
!------------------------------------------------------------------------------
   subroutine get_dr_real_1d(f,df,n_r_max,n_cheb_max,work1,work2,  &
              &              i_costf_init,d_costf_init,drx)
 
      !-- Input variables:
      real(kind=8), intent(in) :: f(*)
      integer,      intent(in) :: n_r_max          ! number of radial grid points
      integer,      intent(in) :: n_cheb_max       ! max number of cheb modes
      real(kind=8), intent(in) :: drx(*)           ! first derivatives of x(r)
      integer,      intent(in) :: i_costf_init(*)  ! info for costf
      real(kind=8), intent(in) :: d_costf_init(*)  ! info for costf
    
      !-- Output variables:
      real(kind=8), intent(out) :: df(*)          ! first derivative of f
      real(kind=8), intent(out) :: work1(*)       ! work array needed for costf
      real(kind=8), intent(out) :: work2(n_r_max) ! work array for f transfer
    
      !-- Local:
      integer :: n_r,n_f
    
    
      !-- Copy input functions:
      do n_r=1,n_r_max
         work2(n_r)=f(n_r)
      end do
    
      !-- Transform f to cheb space:
      call costf1(work2,work1,i_costf_init,d_costf_init)
    
      !-- Get derivatives:
      call get_dcheb(work2,df,n_r_max,n_cheb_max,1.D0)
    
      !-- Transform back:
      call costf1(df,work1,i_costf_init,d_costf_init)
    
      !-- New map:
      do n_r=1,n_r_max
         df(n_r)=drx(n_r)*df(n_r)
      end do

   end subroutine get_dr_real_1d
!------------------------------------------------------------------------------
   subroutine get_drNS(f,df,n_f_max,n_f_start,n_f_stop, &
        &              n_r_max,n_cheb_max,work1,        &
        &              i_costf_init,d_costf_init,drx)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns first radial derivative df of the input function f.      |
      !  |  Array f(n_f_max,*) may contain several functions numbered by     |
      !  |  the first index. The subroutine calculates the derivatives of    |
      !  |  the functions f(n_f_start,*) to f(n_f_stop,*) by transforming    |
      !  |  to a Chebychev representation using n_r_max radial grid points . |
      !  |  Note: when using this function the input field f is slightly     |
      !  |  changed by the back and forth transform. Use s_get_dr.f to       |
      !  |  avoid this.                                                      |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
    
      !-- Input variables:
      integer,         intent(in) :: n_f_max          ! first dim of f
      complex(kind=8), intent(inout) :: f(n_f_max,*)
      integer,         intent(in) :: n_f_start        ! first function to be treated
      integer,         intent(in) :: n_f_stop         ! last function to be treated
      integer,         intent(in) :: n_r_max          ! number of radial grid points
      integer,         intent(in) :: n_cheb_max       ! max number of cheb modes
      real(kind=8),    intent(in) :: drx(*)           ! first derivatives of x(r)
    
      integer,      intent(in) :: i_costf_init(*)  ! info for costf
      real(kind=8), intent(in) :: d_costf_init(*)  ! info for costf
    
      !-- Output variables:
      complex(kind=8), intent(out) :: work1(n_f_max,*) ! work array needed for costf
      complex(kind=8), intent(out) :: df(n_f_max,*)    ! first derivative of f
    
      !-- Local variables:
      integer :: n_r,n_f
    
      !-- Transform f to cheb space:
      call costf1(f,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
    
      !-- Get derivatives:
      call get_dcheb(f,df,n_f_max,n_f_start,n_f_stop,n_r_max,n_cheb_max,1.D0)
    
      !-- Transform back:
      call costf1(f,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
      call costf1(df,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
    
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
        &             i_costf_init,d_costf_init,drx,ddrx)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns first radial derivative df and second radial             |
      !  |  derivative ddf of the input function f.                          |
      !  |  Array f(n_f_max,*) may contain several functions numbered by     |
      !  |  the first index. The subroutine calculates the derivatives of    |
      !  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
      !  |  to a Chebychev representation using n_r_max radial grid points.  |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
    
      !-- Input variables:
      integer,         intent(in) :: n_f_max         ! first dim of f
      complex(kind=8), intent(in) :: f(n_f_max,*)
      integer,         intent(in) :: n_f_start       ! first function to be treated
      integer,         intent(in) :: n_f_stop        ! last function to be treated
      integer,         intent(in) :: n_r_max         ! number of radial grid points
      integer,         intent(in) :: n_cheb_max      ! number of cheb modes
      real(kind=8),    intent(in) :: drx(*)          ! first derivatives of x(r)
      real(kind=8),    intent(in) :: ddrx(*)         ! second derivatives of x(r)
      integer,         intent(in) :: i_costf_init(*) ! info for costf
      real(kind=8),    intent(in) :: d_costf_init(*) ! info for costf
    
      !-- Output variables:
      complex(kind=8), intent(out) :: work1(n_f_max,*)  ! work array needed for costf
      complex(kind=8), intent(out) :: work2(n_f_max,*)  ! work array for f transfer
      complex(kind=8), intent(out) :: df(n_f_max,*)  ! first derivative of f
      complex(kind=8), intent(out) :: ddf(n_f_max,*) ! second derivative of f
    
      !-- Local variables:
      integer :: n_r,n_f
    
      !-- Copy input functions:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do
    
      !-- Transform f to cheb space:
      call costf1(work2,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
    
      !-- Get derivatives:
      call get_ddcheb(work2,df,ddf,n_f_max,n_f_start,n_f_stop, &
                      n_r_max,n_cheb_max,1.D0)
    
      !-- Transform back:
      call costf1(df,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
      call costf1(ddf,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
    
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
                       i_costf_init,d_costf_init,drx,ddrx,dddrx)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns first radial derivative df, the second radial deriv. ddf,|
      !  |  and the third radial derivative dddf of the input function f.    |
      !  |  Array f(n_f_max,*) may contain several functions numbered by     |
      !  |  the first index. The subroutine calculates the derivaties of     |
      !  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
      !  |  to a Chebychev representation using n_r_max radial grid points.  |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,         intent(in) :: n_f_max         ! first dim of f
      complex(kind=8), intent(in) :: f(n_f_max,*)
      integer,         intent(in) :: n_f_start       ! first function to be treated
      integer,         intent(in) :: n_f_stop        ! last function to be treated
      integer,         intent(in) :: n_r_max         ! number of radial grid points
      integer,         intent(in) :: n_cheb_max      ! number of cheb_modes
      real(kind=8),    intent(in) :: drx(*)          ! first derivatives of x(r)
      real(kind=8),    intent(in) :: ddrx(*)         ! second derivatives of x(r)
      real(kind=8),    intent(in) :: dddrx(*)        ! third derivatives of x(r)
      integer,         intent(in) :: i_costf_init(*) ! info for costf
      real(kind=8),    intent(in) :: d_costf_init(*) ! info for costf

      !-- Work arrays:
      complex(kind=8), intent(out) :: work1(n_f_max,*) ! work array needed for costf
      complex(kind=8) , intent(out):: work2(n_f_max,*) ! work array needed for costf

      !-- Output variables:
      complex(kind=8), intent(out) :: df(n_f_max,*)    ! first derivative of f
      complex(kind=8), intent(out) :: ddf(n_f_max,*)   ! second derivative of f
      complex(kind=8), intent(out) :: dddf(n_f_max,*)  ! third derivative of f

      !-- Local variables
      integer :: n_r,n_f

      !-- Copy input functions:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do

      !-- Transform f to cheb space:
      call costf1(work2,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)

      !-- Get derivatives:
      call get_dddcheb(work2,df,ddf,dddf,n_f_max,n_f_start,n_f_stop,  &
                       n_r_max,n_cheb_max,1.D0)

      !-- Transform back:
      call costf1(df,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
      call costf1(ddf,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)
      call costf1(dddf,n_f_max,n_f_start,n_f_stop,work1,i_costf_init,d_costf_init)

      !-- New map:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            dddf(n_f,n_r)=           dddrx(n_r)*df(n_f,n_r) + &
                          3*ddrx(n_r)*drx(n_r)*ddf(n_f,n_r) + &
                                  drx(n_r)**3*dddf(n_f,n_r)
            ddf(n_f,n_r) =ddrx(n_r)*df(n_f,n_r) + drx(n_r)**2*ddf(n_f,n_r)
            df(n_f,n_r)  = drx(n_r)*df(n_f,n_r)
         end do
      end do

   end subroutine get_dddr
!------------------------------------------------------------------------------
   subroutine get_dcheb_complex(f,df,n_f_max,n_f_start,n_f_stop, &
                                n_r_max,n_cheb_max,d_fac)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns chebychev coeffitiens of first derivative df and second  |
      !  |  derivative ddf for a function whose cheb-coeff. are given as     |
      !  |  columns in array f(n_f_max,n_r_max).                             |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,         intent(in) :: n_f_start  ! No of function to start with
      integer,         intent(in) :: n_f_stop   ! No of function to stop with
      integer,         intent(in) :: n_f_max    ! Max no of functions
      integer,         intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,         intent(in) :: n_cheb_max ! Number of cheb modes
      complex(kind=8), intent(in) :: f(n_f_max,*)
      real(kind=8),    intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      complex(kind=8), intent(out) ::  df(n_f_max,*)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(kind=8) :: fac_cheb


      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=0.d0
         end do
      end do
      n_cheb  =n_cheb_max-1
      fac_cheb=d_fac*dble(2*n_cheb)
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*dble(2*n_cheb)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2) + fac_cheb*f(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_dcheb_complex
!------------------------------------------------------------------------------
   subroutine get_dcheb_real_1d(f,df, n_r_max,n_cheb_max,d_fac)

      !-- Input variables:
      integer,      intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,      intent(in) :: n_cheb_max ! Number of cheb modes
      real(kind=8), intent(in) :: f(*)
      real(kind=8), intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      real(kind=8), intent(out) ::  df(*)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(kind=8) :: fac_cheb

      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         df(n_cheb)=0.d0
      end do
      n_cheb  =n_cheb_max-1
      fac_cheb=d_fac*dble(2*n_cheb)
      df(n_cheb)=fac_cheb*f(n_cheb+1)

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*dble(2*n_cheb)
         df(n_cheb)=df(n_cheb+2)+fac_cheb*f(n_cheb+1)
      end do

   end subroutine get_dcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddcheb(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
        &                n_r_max,n_cheb_max,d_fac)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns chebychev coefficents of first derivative df and second  |
      !  |  derivative ddf for a function whose cheb-coeff. are given as     |
      !  |  columns in array f(n_c_tot,n_r_max).                             |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
    
      !-- Input variables:
      integer,         intent(in) :: n_f_start  ! No of column to start with
      integer,         intent(in) :: n_f_stop   ! No of column to stop with
      integer,         intent(in) :: n_f_max    ! First dimension of f,df,ddf
      integer,         intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,         intent(in) :: n_cheb_max ! Number of cheb modes
      complex(kind=8), intent(in) :: f(n_f_max,*)
      real(kind=8),    intent(in) :: d_fac      ! factor for interval mapping
    
      !-- Output variables:
      complex(kind=8), intent(out) ::  df(n_f_max,*)
      complex(kind=8), intent(out) ::  ddf(n_f_max,*)
    
      !-- local variables:
      integer :: n_f,n_cheb
      real(kind=8) :: fac_cheb
    
      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=0.d0
            ddf(n_f,n_cheb)=0.d0
         end do
      end do
      n_cheb=n_cheb_max-1
      fac_cheb=d_fac*dble(2*n_cheb)
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
         ddf(n_f,n_cheb)=0.d0
      end do
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*dble(2*n_cheb)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2) + fac_cheb*f(n_f,n_cheb+1)
            ddf(n_f,n_cheb)=ddf(n_f,n_cheb+2) + fac_cheb*df(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_ddcheb
!------------------------------------------------------------------------------
   subroutine get_dddcheb(f,df,ddf,dddf,n_f_max,n_f_start,n_f_stop, &
                          n_r_max,n_cheb_max,d_fac)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Returns chebychev coeffitiens of first derivative df and second  |
      !  |  derivative ddf for a function whose cheb-coeff. are given as     |
      !  |  columns in array f(n_c_tot,n_r_max).                             |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
         
      !-- Input variables:
      integer,      intent(in) :: n_f_start  ! No of column to start with
      integer,      intent(in) :: n_f_stop   ! No of column to stop with
      integer,      intent(in) :: n_f_max    ! First dimension of f,df,ddf
      integer,      intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,      intent(in) :: n_cheb_max ! Number of cheb modes
      complex(kind=8), intent(in) :: f(n_f_max,*)
      real(kind=8), intent(in) :: d_fac      ! factor for intervall mapping

      !-- Output variables:
      complex(kind=8), intent(out) :: df(n_f_max,*)
      complex(kind=8), intent(out) :: ddf(n_f_max,*)
      complex(kind=8), intent(out) :: dddf(n_f_max,*)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(kind=8) :: fac_cheb

      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=0.d0
            ddf(n_f,n_cheb)=0.d0
            dddf(n_f,n_cheb)=0.d0
         end do
      end do
      n_cheb=n_cheb_max-1
      fac_cheb=d_fac*dble(2*n_cheb)
      do n_f=n_f_start,n_f_stop
         df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
         ddf(n_f,n_cheb)=0.d0
         dddf(n_f,n_cheb)=0.d0
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=d_fac*dble(2*n_cheb)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2) + fac_cheb*f(n_f,n_cheb+1)
            ddf(n_f,n_cheb)=ddf(n_f,n_cheb+2) + fac_cheb*df(n_f,n_cheb+1)
            dddf(n_f,n_cheb)=dddf(n_f,n_cheb+2) + fac_cheb*ddf(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_dddcheb
!------------------------------------------------------------------------------
end module radial_der
