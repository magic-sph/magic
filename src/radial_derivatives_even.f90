module radial_der_even

   use const, only: zero
   use precision_mod
   use cosine_transform, only: costf1, costf2

   implicit none

   private

   public :: get_drNS_even, get_ddrNS_even, get_ddr_even

contains

   subroutine get_ddr_even(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
                           n_r_max,n_cheb_max,dr_fac,work1,work2, &
                                     i_costf1_init,d_costf1_init, &
                                      i_costf2_init,d_costf2_init)
      !
      !  Returns first rarial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !  The cheb transforms have to be initialized by calling            
      !   init_costf1 and init_costf2.                                    
      !

      !-- Input variables:
      integer,     intent(in) :: n_f_max          ! first dim of f
      complex(cp), intent(in) :: f(n_f_max,*)
      integer,     intent(in) :: n_f_start        ! first function to be treated
      integer,     intent(in) :: n_f_stop         ! last function to be treated
      integer,     intent(in) :: n_r_max          ! number of radial grid points
      integer,     intent(in) :: n_cheb_max       ! number of cheb modes
      real(cp),    intent(in) :: dr_fac           ! mapping factor
      integer,     intent(in) :: i_costf1_init(*)
      integer,     intent(in) :: i_costf2_init(*)
      real(cp),    intent(in) :: d_costf1_init(*)
      real(cp),    intent(in) :: d_costf2_init(*)

      !-- Output variables:
      complex(cp), intent(out) :: work1(n_f_max,*)  ! work array needed for costf
      complex(cp), intent(out) :: work2(n_f_max,*)  ! work array needed for costf
      complex(cp), intent(out) :: df(n_f_max,*)  ! first derivative of f
      complex(cp), intent(out) :: ddf(n_f_max,*) ! second derivative of f

      !-- Local variables:
      integer :: n_r,n_f


      !-- Copy input functions:
      do n_r=1,n_r_max
         do n_f=n_f_start,n_f_stop
            work2(n_f,n_r)=f(n_f,n_r)
         end do
      end do

      !-- Transform f to cheb space:
      call costf1(work2,n_f_max,n_f_start,n_f_stop,work1, &
                  i_costf1_init,d_costf1_init)

      !-- Get derivatives:
      call get_ddcheb_even(work2,df,ddf,n_f_max,n_f_start,n_f_stop, &
                           n_r_max,n_cheb_max,dr_fac)

      !-- Transform back, note the different transform used for df,
      !   cause df is odd:
      call costf2(df,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf2_init,d_costf2_init,1)
      call costf1(ddf,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf1_init,d_costf1_init)

   end subroutine get_ddr_even
!------------------------------------------------------------------------------
   subroutine get_drNS_even(f,df,n_f_max,n_f_start,n_f_stop, &
                            n_r_max,n_cheb_max,dr_fac,work1, &
                                i_costf1_init,d_costf1_init, &
                                i_costf2_init,d_costf2_init)
      !
      !  Returns first rarial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !  The cheb transforms have to be initialized by calling            
      !   init_costf1 and init_costf2.                                    
      !

      !-- Input variables:
      integer,  intent(in) :: n_f_max           ! first dim of f
      integer,  intent(in) :: n_f_start         ! first function to be treated
      integer,  intent(in) :: n_f_stop          ! last function to be treated
      integer,  intent(in) :: n_r_max           ! number of radial grid points
      integer,  intent(in) :: n_cheb_max        ! number of cheb modes
      real(cp), intent(in) :: dr_fac            ! mapping factor
      integer,  intent(in) :: i_costf1_init(*)
      integer,  intent(in) :: i_costf2_init(*)
      real(cp), intent(in) :: d_costf1_init(*)
      real(cp), intent(in) :: d_costf2_init(*)

      !-- Output variables:
      complex(cp), intent(inout) :: f(n_f_max,*)
      complex(cp), intent(out) :: work1(n_f_max,*)   ! work array needed for costf
      complex(cp), intent(out) :: df(n_f_max,*)  ! first derivative of f


      !-- Transform f to cheb space:
      call costf1(f,n_f_max,n_f_start,n_f_stop,work1, &
                  i_costf1_init,d_costf1_init)

      !-- Get derivatives:
      call get_dcheb_even(f,df,n_f_max,n_f_start,n_f_stop, &
                          n_r_max,n_cheb_max,dr_fac)

      !-- Transform back, note the different transform used for df,
      !   cause df is odd:
      call costf1(f,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf1_init,d_costf1_init)
      call costf2(df,n_f_max,n_f_start,n_f_stop,work1,&
                  i_costf2_init,d_costf2_init,1)


   end subroutine get_drNS_even
!------------------------------------------------------------------------------
   subroutine get_ddrNS_even(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
                                 n_r_max,n_cheb_max,dr_fac,work1, &
                                     i_costf1_init,d_costf1_init, &
                                     i_costf2_init,d_costf2_init)
      !
      !  Returns first rarial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(n_f_max,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(n_f_start,*) to f(n_f_stop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !  The cheb transforms have to be initialized by calling            
      !   init_costf1 and init_costf2.                                    
      !

      !-- Input variables:
      integer,  intent(in) :: n_f_max          ! first dim of f
      integer,  intent(in) :: n_f_start        ! first function to be treated
      integer,  intent(in) :: n_f_stop         ! last function to be treated
      integer,  intent(in) :: n_r_max          ! number of radial grid points
      integer,  intent(in) :: n_cheb_max       ! number of cheb modes
      real(cp), intent(in) :: dr_fac            ! mapping factor
      integer,  intent(in) :: i_costf1_init(*)
      integer,  intent(in) :: i_costf2_init(*)
      real(cp), intent(in) :: d_costf1_init(*)
      real(cp), intent(in) :: d_costf2_init(*)

      complex(cp), intent(inout) :: f(n_f_max,*)

      !-- Output variables:
      complex(cp), intent(out) :: work1(n_f_max,*)  ! work array needed for costf
      complex(cp), intent(out) :: df(n_f_max,*)  ! first derivative of f
      complex(cp), intent(out) :: ddf(n_f_max,*) ! second derivative of f


      !-- Transform f to cheb space:
      call costf1(f,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf1_init,d_costf1_init)

      !-- Get derivatives:
      call get_ddcheb_even(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
                           n_r_max,n_cheb_max,dr_fac)

      !-- Transform back, note the different transform used for df,
      !   cause df is odd:
      call costf1(f,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf1_init,d_costf1_init)
      call costf2(df,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf2_init,d_costf2_init,1)
      call costf1(ddf,n_f_max,n_f_start,n_f_stop, &
                  work1,i_costf1_init,d_costf1_init)

   end subroutine get_ddrNS_even
!------------------------------------------------------------------------------
   subroutine get_dcheb_even(f,df,n_f_max,n_f_start,n_f_stop, &
                             n_r_max,n_cheb_max,d_fac)
      !
      !  Returns chebychev coeffitiens of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_f_max,n_r_max).                             
      !

      !-- Input variables:
      integer,     intent(in) :: n_f_start  ! No of function to start with
      integer,     intent(in) :: n_f_stop   ! No of function to stop with
      integer,     intent(in) :: n_f_max    ! First dimension of f,df
      integer,     intent(in) :: n_r_max    ! second dimension of f,df
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(n_f_max,*)
      real(cp),    intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      complex(cp), intent(out) ::  df(n_f_max,*)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(cp) :: fac_cheb_odd

      !-- initialize inner core derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=zero
         end do
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-1,1,-1
         fac_cheb_odd=d_fac*real(4*n_cheb,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+1) + fac_cheb_odd*f(n_f,n_cheb+1)
         end do
      end do

   end subroutine get_dcheb_even
!------------------------------------------------------------------------------
   subroutine get_ddcheb_even(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
                              n_r_max,n_cheb_max,d_fac)
      !
      !  Returns chebychev coeffitiens of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_f_max,n_r_max).                             
      !

      !-- Input variables:
      integer,     intent(in) :: n_f_start  ! No of function to start with
      integer,     intent(in) :: n_f_stop   ! No of function to stop with
      integer,     intent(in) :: n_f_max    ! First dimension of f,df,ddf
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(n_f_max,*)
      real(cp),    intent(in) :: d_fac      ! factor for interval mapping

      !-- Output variables:
      complex(cp), intent(out) ::  df(n_f_max,*)
      complex(cp), intent(out) ::  ddf(n_f_max,*)

      !-- Local variables:
      integer :: n_f,n_cheb
      real(cp) :: fac_cheb_even
      real(cp) :: fac_cheb_odd

      !----- initialize inner core derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=zero
            ddf(n_f,n_cheb)=zero
         end do
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-1,1,-1
         fac_cheb_odd=d_fac*real(4*n_cheb,kind=cp)
         fac_cheb_even=d_fac*real(4*n_cheb-2,kind=cp)
         do n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+1) + fac_cheb_odd*f(n_f,n_cheb+1)
            ddf(n_f,n_cheb)=ddf(n_f,n_cheb+1) + fac_cheb_even*df(n_f,n_cheb)
         end do
      end do

   end subroutine get_ddcheb_even
!------------------------------------------------------------------------------
end module radial_der_even
