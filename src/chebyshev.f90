module chebyshev

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, lm_max_real
   use fft_fac_mod, only: fft_fac_complex, fft_fac_real
   use constants, only: half, one, two, three, four, pi, sin36, &
       &                cos36, sin60, sin72, cos72
   use LMLoop_data, only: llm,ulm
   use radial_scheme, only: type_rscheme
   use logic, only: l_newmap, l_PV, l_TO
   use useful, only: factorise

   implicit none

   private

   type, public, extends(type_rscheme) :: type_cheb_odd
      real(cp) :: alpha1 !Input parameter for non-linear map to define degree of spacing (0.0:2.0)
      real(cp) :: alpha2 !Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
      integer, allocatable  :: i_costf_init(:)
      real(cp), allocatable :: d_costf_init(:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: get_der_mat
      procedure :: get_grid => initialize_mapping
      procedure :: costf1_complex_1d
      procedure :: costf1_complex
      procedure :: costf1_real
      procedure :: costf1_real_1d
      !generic :: costf1 => costf1_complex,costf1_real,costf1_real_1d,costf1_complex_1d
   end type type_cheb_odd

   real(cp), allocatable :: r_cheb(:)

   !-- Work arrays that are needed to compute the cosine transforms
   complex(cp), allocatable :: work(:,:)
   complex(cp), allocatable :: work1d(:)
   real(cp), allocatable :: work1d_real(:)

contains

   subroutine initialize(this, n_r_max, order)
      !
      !  Purpose of this subroutine is to calculate and store several     
      !  values that will be needed for a fast cosine transform of the    
      !  first kind. The actual transform is performed by the             
      !  subroutine costf1.                                               
      !

      class(type_cheb_odd) :: this
      
      integer, intent(in) :: n_r_max
      integer, intent(in) :: order ! This is going to be n_cheb_max

      !-- Local variables:
      integer :: j,k,ni,nd
      real(cp) :: theta
      real(cp) :: wr,wi,wpr,wpi,wtemp
      integer :: n_facs,fac(20),n_factors,factor(40)

      this%rnorm = sqrt(two/real(n_r_max-1,kind=cp))
      this%n_max = order  ! n_cheb_max
      this%boundary_fac = half
      this%version = 'cheb'

      allocate( r_cheb(n_r_max) )
      allocate( this%rMat(n_r_max,n_r_max) )
      allocate( this%drMat(n_r_max,n_r_max) )
      allocate( this%d2rMat(n_r_max,n_r_max) )
      allocate( this%d3rMat(n_r_max,n_r_max) )
      bytes_allocated=bytes_allocated+5*n_r_max*n_r_max*SIZEOF_DEF_REAL

      ! if ( l_PV .or. l_TO ) then
         ! allocate( work(1:lm_max,n_r_max) )
      ! else
         ! allocate( work(1:lm_max,n_r_max) )
         !allocate( work(llm:ulm,n_r_max) )
      allocate( work(1:ulm-llm+1,n_r_max) )
      ! end if
      allocate( work1d(n_r_max), work1d_real(n_r_max) )

      ni = 2*n_r_max+2
      nd = 2*n_r_max+5

      allocate( this%d_costf_init(nd) )
      allocate( this%i_costf_init(ni) )
      bytes_allocated = bytes_allocated+nd*SIZEOF_DEF_REAL+ni*SIZEOF_INTEGER

      !-- Checking number of datapoints:
      if ( n_r_max <= 3 ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Sorry, I need more than 3 grid points!'
         stop
      end if

      if ( mod(n_r_max-1,4) /= 0 ) then
         write(*,*) '! Note from subroutine init_costf1:'
         write(*,*) '! Number of data points -1 has to be'
         write(*,*) '! a mutiple of 4!'
         stop
      end if

      if ( nd < 2*n_r_max+5 ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array d_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',2*n_r_max+5
         stop
      end if

      if ( ni < n_r_max+1 ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array i_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n_r_max+1
         stop
      end if

      !-- first information stored in i_costf_init is the dimension:
      this%i_costf_init(1)=n_r_max

      !-- Resorting: ??????????????????
      !   second thing stored in i_costf_init
      this%i_costf_init(2)=1
      this%i_costf_init(3)=2
      this%i_costf_init(n_r_max+1)=n_r_max
      do k=3,n_r_max-2,2
         this%i_costf_init(k+1)=n_r_max+1-k
         this%i_costf_init(k+2)=this%i_costf_init(k+1)+1
      end do
       

      !-- Factorisation of n_r_max for FFT:
      !-- Factors to be checked:
      n_facs=4    ! No of factors to be checked
      fac(1)=4    ! List of factors
      fac(2)=2
      fac(3)=3
      fac(4)=5

      call factorise((n_r_max-1)/2,n_facs,fac,n_factors,factor)

      !-- third info stored in i_costf_init:
      if ( ni <= n_r_max+2+n_factors ) then
         write(*,*) '! Message from subroutine init_costf1:'
         write(*,*) '! Increase dimension of array i_costf_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n_r_max+1+n_factors
         stop
      end if
      this%i_costf_init(n_r_max+2)=n_factors
      do j=1,n_factors
         this%i_costf_init(n_r_max+2+j)=factor(j)
      end do

      !-- Recurrence to get trigonometric auxiliary functions for cos TF
      theta=pi/real(n_r_max-1,cp)
      wr=one
      wi=0.0_cp
      wpr=-two*sin(half*theta)**2  ! = cos(theta)-1
      wpi=sin(theta)
      do j=2,n_r_max/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !-- Normal way of using these:
         !           wr_costf(j-1)=two*wr  ! = 2*cos((j-1)*theta), factor 2 is special !
         !           wi_costf(j-1)=two*wi  ! = 2*sin((j-1)*theta)
         !-- Storage in one array to make subroutine calls simpler:
         this%d_costf_init(j-1)=two*wr
      end do

      !-- Recurrence to get trigonometric auxiliary functions for real TF
      theta=pi/real((n_r_max-1)/2,cp)
      wr=one
      wi=0.0_cp
      wpr=-two*sin(half*theta)**2  ! = cos(theta)-1
      wpi=sin(theta)
      do j=2,n_r_max/4
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         this%d_costf_init(n_r_max/2+2*j-1)=wr  ! = cos((j-1)*theta)
         this%d_costf_init(n_r_max/2+2*j)=wi    ! = sin((j-1)*theta)
      end do

      !-- And this is the way they are needed in fft_fac:
      theta=two*pi/real((n_r_max-1)/2,cp)
      wr=one
      wi=0.0_cp
      wpr=-two*sin(half*theta)**2  ! = cos(theta)-1
      wpi=sin(theta)
      this%d_costf_init(n_r_max+1)=wr           ! = cos(0)
      this%d_costf_init(n_r_max+2)=wi           ! = sin(0)
      do j=2,n_r_max/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         this%d_costf_init(n_r_max+2*j-1)=wr    ! = cos((j-1)*theta)
         this%d_costf_init(n_r_max+2*j  )=wi    ! = sin((j-1)*theta)
      end do

      this%d_costf_init(2*n_r_max+1)=sin36
      this%d_costf_init(2*n_r_max+2)=cos36
      this%d_costf_init(2*n_r_max+3)=sin72
      this%d_costf_init(2*n_r_max+4)=cos72
      this%d_costf_init(2*n_r_max+5)=sin60

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine initialize_mapping(this, n_r_max, ricb, rcmb, ratio1, ratio2, r)

      class(type_cheb_odd) :: this

      !-- Input variables:
      integer,  intent(in) :: n_r_max
      real(cp), intent(in) :: ricb
      real(cp), intent(in) :: rcmb
      real(cp), intent(inout) :: ratio1
      real(cp), intent(in) :: ratio2

      !-- Output variable:
      real(cp), intent(out) :: r(n_r_max)

      !-- Local variables:
      integer :: n_r
      real(cp) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

      allocate ( this%drx(n_r_max), this%ddrx(n_r_max), this%dddrx(n_r_max) )
      bytes_allocated=bytes_allocated+3*n_r_max*SIZEOF_DEF_REAL

      
      if ( l_newmap ) then
         this%alpha1=ratio1
         this%alpha2=ratio2
         paraK=atan(this%alpha1*(one+this%alpha2))/atan(this%alpha1*(one-this%alpha2))
         paraX0=(paraK-one)/(paraK+one)
         lambd=atan(this%alpha1*(one-this%alpha2))/(one-paraX0)
      else
         this%alpha1=0.0_cp
         this%alpha2=0.0_cp
      end if

      call cheb_grid(ricb,rcmb,n_r_max-1,r,r_cheb,this%alpha1,this%alpha2, &
           &         paraX0,lambd)

      if ( l_newmap ) then

         do n_r=1,n_r_max
            this%drx(n_r) =                          (two*this%alpha1) /   &
            &    ((one+this%alpha1**2*(two*r(n_r)-ricb-rcmb-this%alpha2)**2)* &
            &    lambd)
            this%ddrx(n_r) = -(8.0_cp*this%alpha1**3*(two*r(n_r)-ricb-rcmb-this%alpha2)) / &
            &    ((one+this%alpha1**2*(-two*r(n_r)+ricb+rcmb+this%alpha2)**2)**2*     &
            &    lambd)
            this%dddrx(n_r) =        (16.0_cp*this%alpha1**3*(-one+three*this%alpha1**2* &
            &                     (-two*r(n_r)+ricb+rcmb+this%alpha2)**2)) / &
            &    ((one+this%alpha1**2*(-two*r(n_r)+ricb+rcmb+this%alpha2)**2)**3* &
            &    lambd)
         end do

      else

         do n_r=1,n_r_max
            this%drx(n_r)  =two/(rcmb-ricb)
            this%ddrx(n_r) =0.0_cp
            this%dddrx(n_r)=0.0_cp
         end do

      end if

   end subroutine initialize_mapping
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_cheb_odd) :: this

      deallocate( this%d_costf_init )
      deallocate( this%i_costf_init )
      deallocate( this%rMat, this%drMat, this%d2rMat, this%d3rMat )
      deallocate( this%drx, this%ddrx, this%dddrx )
      deallocate( r_cheb )
      deallocate( work, work1d, work1d_real )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine cheb_grid(a,b,n,x,y,a1,a2,x0,lbd)
      !
      !   Given the interval [a,b] the routine returns the
      !   n+1 points that should be used to support a
      !   Chebychev expansion. These are the n+1 extrema y(i) of
      !   the Chebychev polynomial of degree n in the
      !   interval [-1,1].
      !   The respective points mapped into the interval of
      !   question [a,b] are the x(i).
      !
      !   .. note:: x(i) and y(i) are stored in the reversed order:
      !             x(1)=b, x(n+1)=a, y(1)=1, y(n+1)=-1
      !
       
      !-- Input variables
      real(cp), intent(in) :: a,b   ! interval boundaries
      integer,  intent(in) :: n ! degree of Cheb polynomial to be represented by the grid points
      real(cp), intent(in) :: a1,a2,x0,lbd

      !-- Output variables
      real(cp), intent(out) :: x(*) ! grid points in interval [a,b]
      real(cp), intent(out) :: y(*) ! grid points in interval [-1,1]
       
      !-- Local variables:
      real(cp) :: bpa,bma
      integer :: k
       
      bma=half*(b-a)
      bpa=half*(a+b)
       
      do k=1,n+1
         y(k)=cos( pi*real(k-1,cp)/real(n,cp) )
         if ( l_newmap ) then
            x(k)=(a2+tan(lbd*(y(k)-x0))/a1)/2 + bpa
         else
            x(k)=bma * y(k) + bpa
         end if
      end do
        
   end subroutine cheb_grid
!------------------------------------------------------------------------------
   subroutine get_der_mat(this, n_r_max)
      !
      !  Construct Chebychev polynomials and their first, second,
      !  and third derivative up to degree n_r at n_r points x
      !  in the interval [a,b]. Since the Chebs are only defined
      !  in [-1,1] we have to use a map, mapping the points x
      !  points y in the interval [-1,1]. This map is executed
      !  by the subroutine cheb_grid and has to be done
      !  before calling this program.
      !

      class(type_cheb_odd) :: this

      !-- Input variables:
      integer, intent(in) :: n_r_max

      !-- Local variables:
      integer :: n,k   ! counter
      ! in [-1,1] to x-derivatives in [a,b]

      !-- definition of map_fac:
      !   d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
      !                   = map_fac * d Cheb(y) / d y

      !-- construction of chebs and derivatives with recursion:
      do k=1,n_r_max  ! do loop over the n_r grid points !

         !----- set first two chebs:
         this%rMat(1,k)=one
         this%rMat(2,k)=r_cheb(k)
         this%drMat(1,k)=0.0_cp
         this%drMat(2,k)=this%drx(k)
         this%d2rMat(1,k)=0.0_cp
         this%d2rMat(2,k)=this%ddrx(k)
         this%d3rMat(1,k)=0.0_cp
         this%d3rMat(2,k)=this%dddrx(k)

         !----- now construct the rest with a recursion:
         do n=3,n_r_max ! do loop over the (n-1) order of the chebs

            this%rMat(n,k)=    two*r_cheb(k)*this%rMat(n-1,k)-this%rMat(n-2,k)
            this%drMat(n,k)=   two*this%drx(k)*this%rMat(n-1,k) + &
            &                        two*r_cheb(k)*this%drMat(n-1,k) - &
            &                                 this%drMat(n-2,k)
            this%d2rMat(n,k)=  two*this%ddrx(k)*this%rMat(n-1,k) + &
            &                four*this%drx(k)*this%drMat(n-1,k) + &
            &                       two*r_cheb(k)*this%d2rMat(n-1,k) - &
            &                                this%d2rMat(n-2,k)
            this%d3rMat(n,k)=  two*this%dddrx(k)*this%rMat(n-1,k) + &
            &              6.0_cp*this%ddrx(k)*this%drMat(n-1,k) + &
            &             6.0_cp*this%drx(k)*this%d2rMat(n-1,k) + &
            &                       two*r_cheb(k)*this%d3rMat(n-1,k) - &
            &                                this%d3rMat(n-2,k)

         end do

      end do

      this%rMat  =transpose(this%rMat)
      this%drMat =transpose(this%drMat)
      this%d2rMat=transpose(this%d2rMat)
      this%d3rMat=transpose(this%d3rMat)

   end subroutine get_der_mat
!------------------------------------------------------------------------------
   subroutine costf1_complex(this,f,n_f_max,n_f_start,n_f_stop)
      !
      !  Purpose of this subroutine is to perform a multiple
      !  cosine transforms for n+1 datapoints
      !  on the columns numbered n_f_start to n_f_stop in the array
      !  f(n_f_max,n+1)
      !  Depending whether the input f contains data or coeff arrays
      !  coeffs or data are returned in f.
      !
      class(type_cheb_odd) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
    
      !-- Output variables:
      complex(cp), intent(inout) :: f(n_f_max,*) ! data/coeff input
    
      !-- Local variables:
      integer :: n
    
      logical :: l_f2_data
      integer :: n_f
      integer :: j,j1,j2,j3,j4
      integer :: i,i1,i2
      integer :: n_P1  ! n+1
      integer :: n_P2  ! n+2
      integer :: n_P3  ! n+3
      integer :: n_O2  ! n/2
      integer :: n_O2_P1 ! n/2+1
      integer :: n_O2_P2 ! n/2+2
    
    
      complex(cp) :: f_h1,f_h2,f_h3,f_h4 ! help variables
      complex(cp) :: w_h1,w_h2
      real(cp) :: fac_norm,facn
      real(cp) :: wr_j,wi_j,wr_i,wi_i
    
      integer :: n_factors,n_fac,fac,fac_tot
    
      complex(cp) :: tot_sum(lm_max)
    
      if ( n_f_start < 1 ) then
         write(*,*) '! Message from costf1:'
         write(*,*) '! n_f_start should be >=1'
         write(*,*) '! but is:',n_f_start
         stop
      end if
      if ( n_f_stop > n_f_max ) then
         write(*,*) '! Message from costf1:'
         write(*,*) '! n_f_stop > n_f_max !'
         stop
      end if
    
      n=this%i_costf_init(1)-1
    
      n_P1=n+1
      n_P2=n+2
      n_P3=n+3
      n_O2=n/2
      n_O2_P1=n/2+1
      n_O2_P2=n/2+2
    
      !-- Normalisation factor:
      !   The actual normalization factor for the cos transform is the
      !   usual sqrt(2/n). We have in addition a factor 1/2 from the
      !   pre-processing (sum over two f's) and another factor 1/2 from
      !   post-processing (again sum over two f2's).
      fac_norm=one/sqrt(8.0_cp*real(n,cp))
    
    
      !-- Build auxiliary function for cos transform
      !   and shuffle data according to k2k in the process:
      j1=this%i_costf_init(2)
      j2=this%i_costf_init(n_O2_P2)
      do n_f=n_f_start,n_f_stop
         tot_sum(n_f)=f(n_f,1)-f(n_f,n_P1)
         work(n_f,j1)=f(n_f,1)+f(n_f,n_P1)
         work(n_f,j2)=two*f(n_f,n_O2_P1)
      end do
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=this%i_costf_init(j+2)    ! first step
         j2=this%i_costf_init(n_P2-j)
         wr_j=this%d_costf_init(j)
         wi_j=this%d_costf_init(n_O2-j)
    
         i=j+1
         i1=this%i_costf_init(i+2)   ! second step
         i2=this%i_costf_init(n_P2-i)
         wr_i=this%d_costf_init(i)
         wi_i=this%d_costf_init(n_O2-i)

         do n_f=n_f_start,n_f_stop
            f_h1=f(n_f,j+1)+f(n_f,n_P1-j)
            f_h2=f(n_f,j+1)-f(n_f,n_P1-j)
            work(n_f,j1)=f_h1-wi_j*f_h2
            work(n_f,j2)=f_h1+wi_j*f_h2
            tot_sum(n_f)=tot_sum(n_f)+wr_j*f_h2
    
            f_h1=f(n_f,i+1)+f(n_f,n_P1-i)
            f_h2=f(n_f,i+1)-f(n_f,n_P1-i)
            work(n_f,i1)=f_h1-wi_i*f_h2
            work(n_f,i2)=f_h1+wi_i*f_h2
            tot_sum(n_f)=tot_sum(n_f)+wr_i*f_h2
         end do
      end do
    
      j=n/2-1   ! last step
      j1=this%i_costf_init(j+2)
      j2=this%i_costf_init(n_P2-j)
      wr_j=this%d_costf_init(j)
      wi_j=this%d_costf_init(n_O2-j)
    
      do n_f=n_f_start,n_f_stop
         f_h1=f(n_f,j+1)+f(n_f,n_P1-j)
         f_h2=f(n_f,j+1)-f(n_f,n_P1-j)
         work(n_f,j1)=f_h1-wi_j*f_h2
         work(n_f,j2)=f_h1+wi_j*f_h2
         tot_sum(n_f)=tot_sum(n_f)+wr_j*f_h2
      end do
    
      !-- Perform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2
    
      n_factors=this%i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=this%i_costf_init(n+3+n_fac)
         if ( l_f2_data ) then
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_complex(work(1,1),work(1,2),f(1,1),f(1,2),    &
                 &       this%d_costf_init(n+2),n_f_max,n_f_start,      &
                 &       n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_complex(f(1,1),f(1,2),work(1,1),work(1,2),   &
                 &       this%d_costf_init(n+2),n_f_max,n_f_start,     &
                 &       n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.true.
         end if
         fac_tot=fac_tot*fac
      end do
    
      if ( l_f2_data ) then
         !----- Copy data on f2:
         do j1=1,n,4       ! Step size 4: Loop unrolling
            j2=j1+1
            j3=j2+1
            j4=j3+1
            do n_f=n_f_start,n_f_stop
               f(n_f,j1)=work(n_f,j1)
               f(n_f,j2)=work(n_f,j2)
               f(n_f,j3)=work(n_f,j3)
               f(n_f,j4)=work(n_f,j4)
            end do
         end do
      end if
    
      !-- Postprocessing:
      !----- Unscramble two real transforms from the complex transform:
      facn=two*fac_norm
      do n_f=n_f_start,n_f_stop
         f_h1=f(n_f,1)
         f(n_f,1)      =facn*(f(n_f,1)+f(n_f,2))
         f(n_f,2)      =facn*(f_h1-f(n_f,2))
         f(n_f,n_O2_P1)=facn*f(n_f,n_O2_P1)
         f(n_f,n_O2_P2)=facn*f(n_f,n_O2_P2)
      end do
    
      do j=2,n/4
         j2=2*j
         j1=j2-1
         j3=n_P3-j2
         j4=j3+1
    
         wr_j=this%d_costf_init(n_O2+2*j-1)
         wi_j=this%d_costf_init(n_O2+2*j)
    
         do n_f=n_f_start,n_f_stop
            f_h1=fac_norm*(f(n_f,j1)+f(n_f,j3))
            f_h2=fac_norm*(f(n_f,j2)-f(n_f,j4))
            f_h3=fac_norm*(f(n_f,j2)+f(n_f,j4))
            f_h4=fac_norm*(f(n_f,j3)-f(n_f,j1))
    
            w_h1=-wr_j*f_h4+wi_j*f_h3
            w_h2= wr_j*f_h3+wi_j*f_h4
    
            f(n_f,j1)= f_h1+w_h2
            f(n_f,j2)=-f_h2+w_h1
            f(n_f,j3)= f_h1-w_h2
            f(n_f,j4)= f_h2+w_h1
         end do
      end do
    
    
      !---- Extract auxiliary function for costf using recurrence:
      do n_f=n_f_start,n_f_stop
         f(n_f,n+1)=f(n_f,2)
         tot_sum(n_f)=facn*tot_sum(n_f)
         f(n_f,2)=tot_sum(n_f)
      end do
    
      do j=4,n,2
         do n_f=n_f_start,n_f_stop
            tot_sum(n_f)=tot_sum(n_f)+f(n_f,j)
            f(n_f,j)=tot_sum(n_f)
         end do
      end do

   end subroutine costf1_complex
!------------------------------------------------------------------------------
   subroutine costf1_complex_1d(this,f)
    
      class(type_cheb_odd) :: this

      !-- Output variables:
      complex(cp), intent(inout) :: f(*)   ! data/coeff input
    
      !-- Local variables:
      integer :: n
    
      logical :: l_f2_data
      integer :: j,j1,j2,j3,j4
      integer :: i,i1,i2
      integer :: n_P1  ! n+1
      integer :: n_P2  ! n+2
      integer :: n_P3  ! n+3
      integer :: n_O2  ! n/2
      integer :: n_O2_P1 ! n/2+1
      integer :: n_O2_P2 ! n/2+2
    
      complex(cp) :: f_h1,f_h2,f_h3,f_h4 ! help variables
      complex(cp) :: w_h1,w_h2
      complex(cp) :: tot_sum
      real(cp) :: fac_norm,facn
      real(cp) :: wr_j,wi_j,wr_i,wi_i
      integer :: n_factors,n_fac,fac,fac_tot
    
      n=this%i_costf_init(1)-1
    
      n_P1=n+1
      n_P2=n+2
      n_P3=n+3
      n_O2=n/2
      n_O2_P1=n/2+1
      n_O2_P2=n/2+2
    
      !-- Normalisation factor:
      !   The actual normalization factor for the cos transform is the
      !   usual sqrt(2/n). We have in addition a factor 1/2 from the
      !   pre-processing (sum over two f's) and another factor 1/2 from
      !   post-processing (again sum over two f2's).
      fac_norm=one/sqrt(8.0_cp*real(n,cp))
    
      !-- Build auxiliary function for cos transform
      !   and shuffle data according to k2k in the process:
    
      j1=this%i_costf_init(2)
      j2=this%i_costf_init(n_O2_P2)
      tot_sum=f(1)-f(n_P1)
      work1d(j1)=f(1)+f(n_P1)
      work1d(j2)=two*f(n_O2_P1)
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=this%i_costf_init(j+2)    ! first step
         j2=this%i_costf_init(n_P2-j)
         wr_j=this%d_costf_init(j)
         wi_j=this%d_costf_init(n_O2-j)
    
         i=j+1
         i1=this%i_costf_init(i+2)   ! second step
         i2=this%i_costf_init(n_P2-i)
         wr_i=this%d_costf_init(i)
         wi_i=this%d_costf_init(n_O2-i)
    
         f_h1=f(j+1)+f(n_P1-j)
         f_h2=f(j+1)-f(n_P1-j)
         work1d(j1)=f_h1-wi_j*f_h2
         work1d(j2)=f_h1+wi_j*f_h2
         tot_sum=tot_sum+wr_j*f_h2

         f_h1=f(i+1)+f(n_P1-i)
         f_h2=f(i+1)-f(n_P1-i)
         work1d(i1)=f_h1-wi_i*f_h2
         work1d(i2)=f_h1+wi_i*f_h2
         tot_sum=tot_sum+wr_i*f_h2
      end do
    
      j=n/2-1   ! last step
      j1=this%i_costf_init(j+2)
      j2=this%i_costf_init(n_P2-j)
      wr_j=this%d_costf_init(j)
      wi_j=this%d_costf_init(n_O2-j)
    
      f_h1=f(j+1)+f(n_P1-j)
      f_h2=f(j+1)-f(n_P1-j)
      work1d(j1)=f_h1-wi_j*f_h2
      work1d(j2)=f_h1+wi_j*f_h2
      tot_sum=tot_sum+wr_j*f_h2
    
      !-- Perform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2
    
      n_factors=this%i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=this%i_costf_init(n+3+n_fac)
    
         if ( l_f2_data ) then
    
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_complex(work1d(1),work1d(2),f(1),f(2), &
                 &       this%d_costf_init(n+2),1,1,1,  &
                 &       n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_complex(f(1),f(2),work1d(1),work1d(2), &
                 &       this%d_costf_init(n+2),1,1,1,  &
                 &       n_O2,fac,fac_tot)
            l_f2_data=.true.
         end if
         fac_tot=fac_tot*fac
      end do
    
      if ( l_f2_data ) then
         !----- Copy data on work1d:
         do j1=1,n,4       ! Step size 4: Loop unrolling
            j2=j1+1
            j3=j2+1
            j4=j3+1
            f(j1)=work1d(j1)
            f(j2)=work1d(j2)
            f(j3)=work1d(j3)
            f(j4)=work1d(j4)
         end do
      end if
    
      !-- Postprocessing:
      !----- Unscramble two real transforms from the complex transform:
      facn=two*fac_norm
      f_h1=f(1)
      f(1)      =facn*(f(1)+f(2))
      f(2)      =facn*(f_h1-f(2))
      f(n_O2_P1)=facn*f(n_O2_P1)
      f(n_O2_P2)=facn*f(n_O2_P2)
    
      do j=2,n/4
         j2=2*j
         j1=j2-1
         j3=n_P3-j2
         j4=j3+1
    
         wr_j=this%d_costf_init(n_O2+2*j-1)
         wi_j=this%d_costf_init(n_O2+2*j)
    
         f_h1=fac_norm*(f(j1)+f(j3))
         f_h2=fac_norm*(f(j2)-f(j4))
         f_h3=fac_norm*(f(j2)+f(j4))
         f_h4=fac_norm*(f(j3)-f(j1))

         w_h1=-wr_j*f_h4+wi_j*f_h3
         w_h2= wr_j*f_h3+wi_j*f_h4

         f(j1)= f_h1+w_h2
         f(j2)=-f_h2+w_h1
         f(j3)= f_h1-w_h2
         f(j4)= f_h2+w_h1
      end do
    
      !---- Extract auxiliary function for costf using recurrence:
      f(n+1)=f(2)
      tot_sum=facn*tot_sum
      f(2)=tot_sum
    
      do j=4,n,2
         tot_sum=tot_sum+f(j)
         f(j)=tot_sum
      end do

   end subroutine costf1_complex_1d
!------------------------------------------------------------------------------
   subroutine costf1_real(this,f,n_f_max,n_f_start,n_f_stop,f2)

      class(type_cheb_odd) :: this
    
      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
    
      !-- Output variables:
      real(cp), intent(inout) :: f(n_f_max,*)   ! data/coeff input
      real(cp), intent(out) :: f2(n_f_max,*)    ! work array of the same size as f
    
      !-- Local variables:
      integer :: n
    
      logical :: l_f2_data
      integer :: n_f
      integer :: j,j1,j2,j3,j4
      integer :: i,i1,i2
      integer :: n_P1  ! n+1
      integer :: n_P2  ! n+2
      integer :: n_P3  ! n+3
      integer :: n_O2  ! n/2
      integer :: n_O2_P1 ! n/2+1
      integer :: n_O2_P2 ! n/2+2
    
      real(cp) :: f_h1,f_h2,f_h3,f_h4 ! help variables
      real(cp) :: fac_norm,facn
      real(cp) :: w_h1,w_h2
      real(cp) :: wr_j,wi_j,wr_i,wi_i
    
      integer :: n_factors,n_fac,fac,fac_tot
    
      real(cp) :: tot_sum(lm_max_real)
    
      if ( n_f_start < 1 ) then
         write(*,*) '! Message from costf1:'
         write(*,*) '! n_f_start should be >=1'
         write(*,*) '! but is:',n_f_start
         stop
      end if
      if ( n_f_stop > n_f_max ) then
         write(*,*) '! Message from costf1:'
         write(*,*) '! n_f_stop > n_f_max !'
         stop
      end if
    
      n=this%i_costf_init(1)-1
    
      n_P1=n+1
      n_P2=n+2
      n_P3=n+3
      n_O2=n/2
      n_O2_P1=n/2+1
      n_O2_P2=n/2+2
    
      !-- Normalisation factor:
      !   The actual normalization factor for the cos transform is the
      !   usual sqrt(2/n). We have in addition a factor 1/2 from the
      !   pre-processing (sum over two f's) and another factor 1/2 from
      !   post-processing (again sum over two f2's).
      fac_norm=one/sqrt(8.0_cp*real(n,cp))
    
      !-- Build auxiliary function for cos transform
      !   and shuffle data according to k2k in the process:
      j1=this%i_costf_init(2)
      j2=this%i_costf_init(n_O2_P2)
      do n_f=n_f_start,n_f_stop
         tot_sum(n_f)=f(n_f,1)-f(n_f,n_P1)
         f2(n_f,j1)=f(n_f,1)+f(n_f,n_P1)
         f2(n_f,j2)=two*f(n_f,n_O2_P1)
      end do
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=this%i_costf_init(j+2)    ! first step
         j2=this%i_costf_init(n_P2-j)
         wr_j=this%d_costf_init(j)
         wi_j=this%d_costf_init(n_O2-j)
    
         i=j+1
         i1=this%i_costf_init(i+2)   ! second step
         i2=this%i_costf_init(n_P2-i)
         wr_i=this%d_costf_init(i)
         wi_i=this%d_costf_init(n_O2-i)
    
         do n_f=n_f_start,n_f_stop
            f_h1=f(n_f,j+1)+f(n_f,n_P1-j)
            f_h2=f(n_f,j+1)-f(n_f,n_P1-j)
            f2(n_f,j1)=f_h1-wi_j*f_h2
            f2(n_f,j2)=f_h1+wi_j*f_h2
            tot_sum(n_f)=tot_sum(n_f)+wr_j*f_h2
    
            f_h1=f(n_f,i+1)+f(n_f,n_P1-i)
            f_h2=f(n_f,i+1)-f(n_f,n_P1-i)
            f2(n_f,i1)=f_h1-wi_i*f_h2
            f2(n_f,i2)=f_h1+wi_i*f_h2
            tot_sum(n_f)=tot_sum(n_f)+wr_i*f_h2
         end do
      end do
    
      j=n/2-1   ! last step
      j1=this%i_costf_init(j+2)
      j2=this%i_costf_init(n_P2-j)
      wr_j=this%d_costf_init(j)
      wi_j=this%d_costf_init(n_O2-j)
    
      do n_f=n_f_start,n_f_stop
         f_h1=f(n_f,j+1)+f(n_f,n_P1-j)
         f_h2=f(n_f,j+1)-f(n_f,n_P1-j)
         f2(n_f,j1)=f_h1-wi_j*f_h2
         f2(n_f,j2)=f_h1+wi_j*f_h2
         tot_sum(n_f)=tot_sum(n_f)+wr_j*f_h2
      end do
    
      !-- Perform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2
    
      n_factors=this%i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=this%i_costf_init(n+3+n_fac)
         if ( l_f2_data ) then
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_real(f2(1,1),f2(1,2),f(1,1),f(1,2),  &
                 &       this%d_costf_init(n+2),n_f_max,n_f_start, &
                 &       n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_real(f(1,1),f(1,2),f2(1,1),f2(1,2),  &
                 &       this%d_costf_init(n+2),n_f_max,n_f_start, &
                 &       n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.true.
         end if
         fac_tot=fac_tot*fac
      end do
    
      if ( l_f2_data ) then
         !----- Copy data on f2:
         do j1=1,n,4       ! Step size 4: Loop unrolling
            j2=j1+1
            j3=j2+1
            j4=j3+1
            do n_f=n_f_start,n_f_stop
               f(n_f,j1)=f2(n_f,j1)
               f(n_f,j2)=f2(n_f,j2)
               f(n_f,j3)=f2(n_f,j3)
               f(n_f,j4)=f2(n_f,j4)
            end do
         end do
      end if
    
      !-- Postprocessing:
      !----- Unscramble two real transforms from the complex transform:
      facn=two*fac_norm
      do n_f=n_f_start,n_f_stop
         f_h1=f(n_f,1)
         f(n_f,1)      =facn*(f(n_f,1)+f(n_f,2))
         f(n_f,2)      =facn*(f_h1-f(n_f,2))
         f(n_f,n_O2_P1)=facn*f(n_f,n_O2_P1)
         f(n_f,n_O2_P2)=facn*f(n_f,n_O2_P2)
      end do
    
      do j=2,n/4
         j2=2*j
         j1=j2-1
         j3=n_P3-j2
         j4=j3+1
    
         wr_j=this%d_costf_init(n_O2+2*j-1)
         wi_j=this%d_costf_init(n_O2+2*j)
         do n_f=n_f_start,n_f_stop
            f_h1=fac_norm*(f(n_f,j1)+f(n_f,j3))
            f_h2=fac_norm*(f(n_f,j2)-f(n_f,j4))
            f_h3=fac_norm*(f(n_f,j2)+f(n_f,j4))
            f_h4=fac_norm*(f(n_f,j3)-f(n_f,j1))
    
            w_h1=-wr_j*f_h4+wi_j*f_h3
            w_h2= wr_j*f_h3+wi_j*f_h4
    
            f(n_f,j1)= f_h1+w_h2
            f(n_f,j2)=-f_h2+w_h1
            f(n_f,j3)= f_h1-w_h2
            f(n_f,j4)= f_h2+w_h1
         end do
      end do
    
      !---- Extract auxiliary function for costf using recurrence:
      do n_f=n_f_start,n_f_stop
         f(n_f,n+1)=f(n_f,2)
         tot_sum(n_f)=facn*tot_sum(n_f)
         f(n_f,2)=tot_sum(n_f)
      end do
    
      do j=4,n,2
         do n_f=n_f_start,n_f_stop
            tot_sum(n_f)=tot_sum(n_f)+f(n_f,j)
            f(n_f,j)=tot_sum(n_f)
         end do
      end do

   end subroutine costf1_real
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this,f)
    
      class(type_cheb_odd) :: this

      !-- Output variables:
      real(cp), intent(inout) :: f(*)   ! data/coeff input
    
      !-- Local variables:
      integer :: n
    
      logical :: l_f2_data
      integer :: j,j1,j2,j3,j4
      integer :: i,i1,i2
      integer :: n_P1  ! n+1
      integer :: n_P2  ! n+2
      integer :: n_P3  ! n+3
      integer :: n_O2  ! n/2
      integer :: n_O2_P1 ! n/2+1
      integer :: n_O2_P2 ! n/2+2
    
      real(cp) :: f_h1,f_h2,f_h3,f_h4 ! help variables
      real(cp) :: fac_norm,facn
      real(cp) :: w_h1,w_h2
      real(cp) :: wr_j,wi_j,wr_i,wi_i
      real(cp) :: tot_sum
      integer :: n_factors,n_fac,fac,fac_tot

      n=this%i_costf_init(1)-1
    
      n_P1=n+1
      n_P2=n+2
      n_P3=n+3
      n_O2=n/2
      n_O2_P1=n/2+1
      n_O2_P2=n/2+2
    
      !-- Normalisation factor:
      !   The actual normalization factor for the cos transform is the
      !   usual sqrt(2/n). We have in addition a factor 1/2 from the
      !   pre-processing (sum over two f's) and another factor 1/2 from
      !   post-processing (again sum over two f2's).
      fac_norm=one/sqrt(8.0_cp*real(n,cp))
    
      !-- Build auxiliary function for cos transform
      !   and shuffle data according to k2k in the process:
    
      j1=this%i_costf_init(2)
      j2=this%i_costf_init(n_O2_P2)
      tot_sum=f(1)-f(n_P1)
      work1d_real(j1)=f(1)+f(n_P1)
      work1d_real(j2)=two*f(n_O2_P1)
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=this%i_costf_init(j+2)    ! first step
         j2=this%i_costf_init(n_P2-j)
         wr_j=this%d_costf_init(j)
         wi_j=this%d_costf_init(n_O2-j)
    
         i=j+1
         i1=this%i_costf_init(i+2)   ! second step
         i2=this%i_costf_init(n_P2-i)
         wr_i=this%d_costf_init(i)
         wi_i=this%d_costf_init(n_O2-i)
    
         f_h1=f(j+1)+f(n_P1-j)
         f_h2=f(j+1)-f(n_P1-j)
         work1d_real(j1)=f_h1-wi_j*f_h2
         work1d_real(j2)=f_h1+wi_j*f_h2
         tot_sum=tot_sum+wr_j*f_h2

         f_h1=f(i+1)+f(n_P1-i)
         f_h2=f(i+1)-f(n_P1-i)
         work1d_real(i1)=f_h1-wi_i*f_h2
         work1d_real(i2)=f_h1+wi_i*f_h2
         tot_sum=tot_sum+wr_i*f_h2
      end do
    
      j=n/2-1   ! last step
      j1=this%i_costf_init(j+2)
      j2=this%i_costf_init(n_P2-j)
      wr_j=this%d_costf_init(j)
      wi_j=this%d_costf_init(n_O2-j)
    
      f_h1=f(j+1)+f(n_P1-j)
      f_h2=f(j+1)-f(n_P1-j)
      work1d_real(j1)=f_h1-wi_j*f_h2
      work1d_real(j2)=f_h1+wi_j*f_h2
      tot_sum=tot_sum+wr_j*f_h2
    
      !-- Perform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2
    
      n_factors=this%i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=this%i_costf_init(n+3+n_fac)
    
         if ( l_f2_data ) then
    
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_real(work1d_real(1),work1d_real(2),f(1),f(2),  &
                 &       this%d_costf_init(n+2),1,1,1,     &
                 &       n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_real(f(1),f(2),work1d_real(1),work1d_real(2), &
                 &       this%d_costf_init(n+2),1,1,1,    &
                 &       n_O2,fac,fac_tot)
            l_f2_data=.true.
         end if
         fac_tot=fac_tot*fac
      end do
    
      if ( l_f2_data ) then
         !----- Copy data on f2:
         do j1=1,n,4       ! Step size 4: Loop unrolling
            j2=j1+1
            j3=j2+1
            j4=j3+1
            f(j1)=work1d_real(j1)
            f(j2)=work1d_real(j2)
            f(j3)=work1d_real(j3)
            f(j4)=work1d_real(j4)
         end do
      end if
    
      !-- Postprocessing:
      !----- Unscramble two real transforms from the complex transform:
      facn=two*fac_norm
      f_h1=f(1)
      f(1)      =facn*(f(1)+f(2))
      f(2)      =facn*(f_h1-f(2))
      f(n_O2_P1)=facn*f(n_O2_P1)
      f(n_O2_P2)=facn*f(n_O2_P2)
    
      do j=2,n/4
         j2=2*j
         j1=j2-1
         j3=n_P3-j2
         j4=j3+1
    
         wr_j=this%d_costf_init(n_O2+2*j-1)
         wi_j=this%d_costf_init(n_O2+2*j)
    
         f_h1=fac_norm*(f(j1)+f(j3))
         f_h2=fac_norm*(f(j2)-f(j4))
         f_h3=fac_norm*(f(j2)+f(j4))
         f_h4=fac_norm*(f(j3)-f(j1))

         w_h1=-wr_j*f_h4+wi_j*f_h3
         w_h2= wr_j*f_h3+wi_j*f_h4

         f(j1)= f_h1+w_h2
         f(j2)=-f_h2+w_h1
         f(j3)= f_h1-w_h2
         f(j4)= f_h2+w_h1
      end do
    
      !---- Extract auxiliary function for costf using recurrence:
      f(n+1)=f(2)
      tot_sum=facn*tot_sum
      f(2)=tot_sum
    
      do j=4,n,2
         tot_sum=tot_sum+f(j)
         f(j)=tot_sum
      end do

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
end module chebyshev
