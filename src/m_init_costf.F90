!$Id$
module init_costf
 
   use const, only: pi
   use useful, only: factorise

   implicit none

   private

   public :: init_costf1, init_costf2

contains

   subroutine init_costf1(n,i_costf_init,ni,d_costf_init,nd)
      !  +-------------+----------------+------------------------------------+
      !  |  Purpose of this subroutine is to calculate and store several     |
      !  |  values that will be needed for a fast cosine transform of the    |
      !  |  frist kind. The actual transform is performed by the             |
      !  |  subroutine costf1.                                               |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer, intent(in) :: n          ! No of grid points !
      integer, intent(in) :: ni         ! dimension of i_costf_init
      integer, intent(in) :: nd         ! dimension of d_costf_init

      !-- Output variables:
      integer,      intent(out) :: i_costf_init(ni) ! array for integers
      real(kind=8), intent(out) :: d_costf_init(nd) ! array for integers

      !-- Local variables:
      integer :: j,k
      real(kind=8) :: theta
      real(kind=8) :: wr,wi,wpr,wpi,wtemp
      integer :: n_facs,fac(20),n_factors,factor(40)

      real(kind=8) :: rad


      !-- Checking number of datapoints:
      if ( n <= 3 ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Sorry, I need more than 3 grid points!'
         stop
      end if

      if ( mod(n-1,4) /= 0 ) then
         write(*,*) '! Note from subroutine init_costf1:'
         write(*,*) '! Number of data points -1 has to be'
         write(*,*) '! a mutiple of 4!'
         stop
      end if

      if ( nd < 2*n+5 ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array d_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',2*n+5
         stop
      end if

      if ( ni < n+1 ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array i_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n+1
         stop
      end if

      !-- first information stored in i_costf_init is the dimension:
      i_costf_init(1)=n

      !-- Resorting: ??????????????????
      !   second thing stored in i_costf_init
      i_costf_init(2)=1
      i_costf_init(3)=2
      i_costf_init(n+1)=n
      do k=3,n-2,2
         i_costf_init(k+1)=n+1-k
         i_costf_init(k+2)=i_costf_init(k+1)+1
      end do
       

      !-- Factorisation of n for FFT:
      !-- Factors to be checked:
      n_facs=4    ! No of factors to be checked
      fac(1)=4    ! List of factors
      fac(2)=2
      fac(3)=3
      fac(4)=5

      call factorise((n-1)/2,n_facs,fac,n_factors,factor)

      !-- third info stored in i_costf_init:
      if ( ni <= n+2+n_factors ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array i_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n+1+n_factors
         stop
      end if
      i_costf_init(n+2)=n_factors
      do j=1,n_factors
         i_costf_init(n+2+j)=factor(j)
      end do

      !-- Recurrence to get trigonometric auxiliary functions for cos TF
      theta=pi/dble(n-1)
      wr=1.d0
      wi=0.d0
      wpr=-2.d0*dsin(0.5d0*theta)**2  ! = cos(theta)-1
      wpi=dsin(theta)
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !-- Normal way of using these:
         !           wr_costf(j-1)=2.d0*wr  ! = 2*cos((j-1)*theta), factor 2 is special !
         !           wi_costf(j-1)=2.d0*wi  ! = 2*sin((j-1)*theta)
         !-- Storage in one array to make subroutine calls simpler:
         d_costf_init(j-1)=2.d0*wr
      end do

      !-- Recurrence to get trigonometric auxiliary functions for real TF
      theta=pi/dble((n-1)/2)
      wr=1.d0
      wi=0.d0
      wpr=-2.d0*dsin(0.5d0*theta)**2  ! = cos(theta)-1
      wpi=dsin(theta)
      do j=2,n/4
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         d_costf_init(n/2+2*j-1)=wr  ! = cos((j-1)*theta)
         d_costf_init(n/2+2*j)=wi    ! = sin((j-1)*theta)
      end do

      !-- And this is the way they are needed in fft_fac:
      theta=2.d0*pi/dble((n-1)/2)
      wr=1.d0
      wi=0.d0
      wpr=-2.d0*dsin(0.5d0*theta)**2  ! = cos(theta)-1
      wpi=dsin(theta)
      d_costf_init(n+1)=wr           ! = cos(0)
      d_costf_init(n+2)=wi           ! = sin(0)
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         d_costf_init(n+2*j-1)=wr    ! = cos((j-1)*theta)
         d_costf_init(n+2*j  )=wi    ! = sin((j-1)*theta)
      end do

      rad=pi/180.D0
      d_costf_init(2*n+1)=dsin(36.D0*rad)
      d_costf_init(2*n+2)=dcos(36.D0*rad)
      d_costf_init(2*n+3)=dsin(72.D0*rad)
      d_costf_init(2*n+4)=dcos(72.D0*rad)
      d_costf_init(2*n+5)=dsin(60.D0*rad)

   end subroutine init_costf1
!------------------------------------------------------------------------------
   subroutine init_costf2(n,i_costf_init,ni,d_costf_init,nd)
      !  +-------------+----------------+------------------------------------+
      !  |  Purpose of this subroutine is to calculate several things        |
      !  |  needed for the cheb transform.                                   |
      !  |  Prepares costf2 for even number of grid points.                  |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer, intent(in) :: n           ! No of grid points !
      integer, intent(in) :: ni          ! dimension of i_costf_init
      integer, intent(in) :: nd          ! dimension of i_costf_init

      !-- Output variables:
      integer,      intent(out) :: i_costf_init(ni) ! array for integers
      real(kind=8), intent(out) :: d_costf_init(nd) ! array for integers

      !-- Local variables:
      integer :: j,k
      real(kind=8) :: theta
      real(kind=8) :: wr,wi,wpr,wpi,wtemp
      integer :: n_facs,fac(20),n_factors,factor(40)

      real(kind=8) :: rad


      !-- Checking number of datapoints:
      if ( n <= 3 ) then
         write(*,*) '! Note from subroutine init_costf2:'
         write(*,*) '! Sorry, I need more than 3 grid points!'
         stop
      end if

      if ( mod(n,4) /= 0 ) then
         write(*,*) '! Note from subroutine init_costf2:'
         write(*,*) '! Number of data points -1 has to be'
         write(*,*) '! a mutiple of 4!'
         stop
      end if

      if ( nd < 2*n+n/2+5 ) then
         write(*,*) '! Message from subroutine init_costf2:'
         write(*,*) '! Increase dimension of array d_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',2*n+n/2+5
         stop
      end if

      if ( ni < n+1 ) then
         write(*,*) '! Message from subroutine init_costf2:'
         write(*,*) '! Increase dimension of array i_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n
         stop
      end if


      !-- first information stored in i_costf_init is the dimension:
      i_costf_init(1)=n

      !-- Resorting: ??????????????????
      !   second thing stored in i_costf_init
      i_costf_init(2)=1
      i_costf_init(3)=2
      do k=3,n-1,2
         i_costf_init(k+1)=n+2-k
         i_costf_init(k+2)=i_costf_init(k+1)+1
      end do

      !-- Factorisation of n for FFT:
      !-- Factors to be checked:
      n_facs=4    ! No of factors to be checked
      fac(1)=4    ! List of factors
      fac(2)=2
      fac(3)=3
      fac(4)=5

      call factorise(n/2,n_facs,fac,n_factors,factor)

      !-- Third info stored in i_costf_init:
      if ( ni < n+2+n_factors ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array i_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n+2+n_factors
         stop
      end if
      i_costf_init(n+2)=n_factors
      do j=1,n_factors
         i_costf_init(n+2+j)=factor(j)
      end do

      !-- Recurrencies to get trigonometric auxiliary functions for second cos TF
      !   only needed if IC included
      theta=0.5d0*pi/dble(n)
      wr=dcos(theta)
      wi=dsin(theta)
      wpr=-2.d0*wi**2              ! = cos(2*theta)-1
      wpi=dsin(2.d0*theta)
      !wi_costf_2(1)=2.d0*wi        ! = 2*sin(theta)
      !w_save(n+1)=2.d0*wi
      d_costf_init(2*n+1)=2.d0*wi
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !wi_costf_2(j)=2.d0*wi     ! = 2*sin((2*j+1)*theta)
         !w_save(n+j)=2.d0*wi
         d_costf_init(2*n+j)=2.d0*wi
      end do

      wr=wpr+1.d0
      wi=wpi
      !wr_costf(1)=wr             ! = cos(2*theta)
      !wi_costf(1)=wi             ! = sin(2*theta)
      !w_save(1)=wr
      d_costf_init(1)=wr
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !wr_costf(j)=wr          ! = cos(2*j*theta)
         !wi_costf(j)=wi          ! = sin(2*j*theta)
         !w_save(j)=wr
         d_costf_init(j)=wr
      end do

      !-- Recurrence to get trigonometric auxiliary functions for real TF
      theta=pi/dble(n/2)
      wr=1.d0
      wi=0.d0
      wpr=-2.d0*dsin(0.5d0*theta)**2  ! Thats cos(theta)-1
      wpi=dsin(theta)
      do j=2,n/4
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !wr_realtf(j)=wr          ! Thats cos((j-1)*theta)
         !wi_realtf(j)=wi          ! Thats sin((j-1)*theta)
         !w_save(n/2+2*j-1)=wr
         !w_save(n/2+2*j)=wi
         d_costf_init(n/2+2*j-1)=wr  ! = cos((j-1)*theta)
         d_costf_init(n/2+2*j)=wi    ! = sin((j-1)*theta)
      end do


      !-- And this is the way they are needed in fft_fac:
      theta=2.d0*pi/dble(n/2)
      wr=1.d0
      wi=0.d0
      wpr=-2.d0*dsin(0.5d0*theta)**2  ! Thats cos(theta)-1
      wpi=dsin(theta)
      d_costf_init(n+1)=wr           ! = cos(0)
      d_costf_init(n+2)=wi           ! = sin(0)
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         d_costf_init(n+2*j-1)=wr    ! = cos((j-1)*theta)
         d_costf_init(n+2*j  )=wi    ! = sin((j-1)*theta)
      end do

      !-- JW addition to NumRes routine:
      rad=pi/180.D0
      d_costf_init(2*n+n/2+1)=dsin(36.D0*rad)
      d_costf_init(2*n+n/2+2)=dcos(36.D0*rad)
      d_costf_init(2*n+n/2+3)=dsin(72.D0*rad)
      d_costf_init(2*n+n/2+4)=dcos(72.D0*rad)
      d_costf_init(2*n+n/2+5)=dsin(60.D0*rad)

   end subroutine init_costf2
!------------------------------------------------------------------------------
end module init_costf
