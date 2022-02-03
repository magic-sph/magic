module cosine_transform_even

   use iso_fortran_env, only: output_unit
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max
   use fft_fac_mod, only: fft_fac_complex
   use constants, only: half, one, two, pi, sin36, cos36, sin60, sin72, cos72
   use useful, only: factorise, abortRun

   implicit none

   private

   type, public :: costf_even_t
      integer, allocatable  :: i_costf_init(:)
      real(cp), allocatable :: d_costf_init(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: costf2
   end type costf_even_t

contains

   subroutine initialize(this, n, ni, nd)
      !
      !  Purpose of this subroutine is to calculate several things
      !  needed for the Chebyshev transform.
      !  Prepares ``costf2`` for even number of grid points.
      !

      class(costf_even_t) :: this

      !-- Input variables:
      integer, intent(in) :: n           ! No of grid points !
      integer, intent(in) :: ni          ! dimension of i_costf_init
      integer, intent(in) :: nd          ! dimension of i_costf_init

      !-- Local variables:
      integer :: j,k
      real(cp) :: theta
      real(cp) :: wr,wi,wpr,wpi,wtemp
      integer :: n_facs,fac(20),n_factors,factor(40)

      allocate( this%d_costf_init(nd), this%i_costf_init(ni) )
      bytes_allocated = bytes_allocated+nd*SIZEOF_DEF_REAL+ni*SIZEOF_INTEGER

      !-- Checking number of datapoints:
      if ( n <= 3 ) call abortRun('At least 3 points: stop run in costf_even')
      if ( mod(n,4) /= 0 ) call abortRun('n-1 has to be a multiple of 4: stop run in costf_even')
      if ( nd < 2*n+n/2+5 ) call abortRun('d_costf_init size should be larger')
      if ( ni < n+1 ) call abortRun('i_costf_init size should be larger')

      !-- first information stored in i_costf_init is the dimension:
      this%i_costf_init(1)=n

      !-- Resorting: ??????????????????
      !   second thing stored in i_costf_init
      this%i_costf_init(2)=1
      this%i_costf_init(3)=2
      do k=3,n-1,2
         this%i_costf_init(k+1)=n+2-k
         this%i_costf_init(k+2)=this%i_costf_init(k+1)+1
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
         write(output_unit,*) '! Message from subroutine init_costf1:'
         write(output_unit,*) '! Increase dimension of array i_costf_init'
         write(output_unit,*) '! in calling routine.'
         write(output_unit,*) '! Should be at least:',n+2+n_factors
         call abortRun('Stop run in costf_even')
      end if
      this%i_costf_init(n+2)=n_factors
      do j=1,n_factors
         this%i_costf_init(n+2+j)=factor(j)
      end do

      !-- Recurrencies to get trigonometric auxiliary functions for second cos TF
      !   only needed if IC included
      theta=half*pi/real(n,cp)
      wr=cos(theta)
      wi=sin(theta)
      wpr=-two*wi**2              ! = cos(2*theta)-1
      wpi=sin(two*theta)
      !wi_costf_2(1)=two*wi        ! = 2*sin(theta)
      !w_save(n+1)=two*wi
      this%d_costf_init(2*n+1)=two*wi
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !wi_costf_2(j)=two*wi     ! = 2*sin((2*j+1)*theta)
         !w_save(n+j)=two*wi
         this%d_costf_init(2*n+j)=two*wi
      end do

      wr=wpr+one
      wi=wpi
      !wr_costf(1)=wr             ! = cos(2*theta)
      !wi_costf(1)=wi             ! = sin(2*theta)
      !w_save(1)=wr
      this%d_costf_init(1)=wr
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !wr_costf(j)=wr          ! = cos(2*j*theta)
         !wi_costf(j)=wi          ! = sin(2*j*theta)
         !w_save(j)=wr
         this%d_costf_init(j)=wr
      end do

      !-- Recurrence to get trigonometric auxiliary functions for real TF
      theta=pi/real(n/2,cp)
      wr=one
      wi=0.0_cp
      wpr=-two*sin(half*theta)**2  ! Thats cos(theta)-1
      wpi=sin(theta)
      do j=2,n/4
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         !wr_realtf(j)=wr          ! Thats cos((j-1)*theta)
         !wi_realtf(j)=wi          ! Thats sin((j-1)*theta)
         !w_save(n/2+2*j-1)=wr
         !w_save(n/2+2*j)=wi
         this%d_costf_init(n/2+2*j-1)=wr  ! = cos((j-1)*theta)
         this%d_costf_init(n/2+2*j)=wi    ! = sin((j-1)*theta)
      end do


      !-- And this is the way they are needed in fft_fac:
      theta=two*pi/real(n/2,cp)
      wr=one
      wi=0.0_cp
      wpr=-two*sin(half*theta)**2  ! Thats cos(theta)-1
      wpi=sin(theta)
      this%d_costf_init(n+1)=wr           ! = cos(0)
      this%d_costf_init(n+2)=wi           ! = sin(0)
      do j=2,n/2
         wtemp=wr
         wr=wr*wpr-wi*wpi+wr
         wi=wi*wpr+wtemp*wpi+wi
         this%d_costf_init(n+2*j-1)=wr    ! = cos((j-1)*theta)
         this%d_costf_init(n+2*j  )=wi    ! = sin((j-1)*theta)
      end do

      !-- JW addition to NumRes routine:
      this%d_costf_init(2*n+n/2+1)=sin36
      this%d_costf_init(2*n+n/2+2)=cos36
      this%d_costf_init(2*n+n/2+3)=sin72
      this%d_costf_init(2*n+n/2+4)=cos72
      this%d_costf_init(2*n+n/2+5)=sin60

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(costf_even_t) :: this

      deallocate( this%d_costf_init, this%i_costf_init )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine costf2(this,f,n_f_max,n_f_start,n_f_stop,f2)
      !
      !  Purpose of this subroutine is to perform a multiple
      !  cosine transforms for n+1 datapoints
      !  on the columns numbered n_f_start to n_f_stop in the array
      !  ``y(n_f_max,n+1)``
      !  Depending whether the input y contains data or coeff arrays
      !  coeffs or data are returned in y.
      !

      class(costf_even_t) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in y,y2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed

      !-- Output variables
      complex(cp), intent(inout) :: f(n_f_max,*)   ! data/coeff input
      complex(cp), intent(out) :: f2(n_f_max,*)    ! work array of the same size as y

      !-- Local variables:
      logical :: l_f2_data
      integer :: n, n_f, n_fac, fac, fac_tot, n_factors
      integer :: i, j, j1, j2, j3, j4, k1, k2, k3, k4
      integer :: n_P1, n_P2, n_P3, n_O2, n_O2_P1, n_O2_P2, n_O2_P3
      real(cp) :: fac_norm,facn,wr_j,wi_j,wr_i,wi_i
      complex(cp) :: f_h1, f_h2, f_h3, f_h4 ! help variables
      complex(cp) :: w_h1, w_h2
      complex(cp) :: sum(lm_max)

      n=this%i_costf_init(1)

      n_P1=n+1
      n_P2=n+2
      n_P3=n+3
      n_O2=n/2
      n_O2_P1=n/2+1
      n_O2_P2=n/2+2
      n_O2_P3=n/2+3

      !-- Normalisation factor:
      !   The actual normalisition factor for the cos transform is the
      !   usual sqrt(2/n). We have in addition a factor 1/2 from the
      !   pre-processing (sum over two f's) and another factor 1/2 from
      !   post-processing (again sum over two f2's).
      fac_norm=one/sqrt(real(8*n,cp))

      !-- Build auxiliary function for cos transform
      !   and shuffle data according to k2k in the process:
      do j=1,n/2-1,2    ! step 2 unrolling

         j1=this%i_costf_init(j+1)    ! first step
         j2=this%i_costf_init(n_P2-j)
         wi_j=this%d_costf_init(2*n+j)

         i=j+1
         j3=this%i_costf_init(i+1)   ! second step
         j4=this%i_costf_init(n_P2-i)
         wi_i=this%d_costf_init(2*n+i)

         do n_f=n_f_start,n_f_stop
            f_h1=f(n_f,j)+f(n_f,n_P1-j)
            f_h2=wi_j*(f(n_f,j)-f(n_f,n_P1-j))
            f2(n_f,j1)=f_h1+f_h2
            f2(n_f,j2)=f_h1-f_h2

            f_h1=f(n_f,i)+f(n_f,n_P1-i)
            f_h2=wi_i*(f(n_f,i)-f(n_f,n_P1-i))
            f2(n_f,j3)=f_h1+f_h2
            f2(n_f,j4)=f_h1-f_h2
         end do
      end do

      !----- Preform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2

      n_factors=this%i_costf_init(n+2)

      do n_fac=1,n_factors   ! loop over factors
         fac=this%i_costf_init(n+2+n_fac)
         if ( l_f2_data ) then
            !-- Vpassm returns complex transform of f2's on f's:
            call fft_fac_complex(f2(1,1),f2(1,2),f(1,1),f(1,2),      &
                 &               this%d_costf_init(n+1),n_f_max,     &
                 &               n_f_start,n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !-- Vpassm returns complex transform of f's on f2's:
            call fft_fac_complex(f(1,1),f(1,2),f2(1,1),f2(1,2),      &
                 &               this%d_costf_init(n+1),n_f_max,     &
                 &               n_f_start,n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.true.
         end if
         fac_tot=fac_tot*fac
      end do

      if ( l_f2_data ) then
         !----- Copy data on f:
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

      !----- Postprocessing:
      !----- Unscramble two real transforms from the complex transform:
      facn=two*fac_norm
      do n_f=n_f_start,n_f_stop
         f_h1=f(n_f,1)
         f(n_f,1)      =facn*(f_h1+f(n_f,2))
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

      !----- Extract auxiliary function for cos TF:
      do j1=3,n,2
         j2=j1+1
         i=(j1-1)/2

         wr_i=this%d_costf_init(i)
         wi_i=this%d_costf_init(n_O2-i)

         do n_f=n_f_start,n_f_stop
            f_h1=  wr_i*f(n_f,j1) - wi_i*f(n_f,j2)
            f_h2=  wr_i*f(n_f,j2) + wi_i*f(n_f,j1)
            f(n_f,j1)=f_h1
            f(n_f,j2)=f_h2
         end do
      end do

      !----- Initialize recurrence:
      do n_f=n_f_start,n_f_stop
         sum(n_f)=half*f(n_f,2)
      end do

      !----- Carry out recurrence for odd terms, even terms unchanged:
      do j=n,2,-2
         do n_f=n_f_start,n_f_stop
            f_h1=sum(n_f)
            sum(n_f)=sum(n_f)+f(n_f,j)
            f(n_f,j)=f_h1
         end do
      end do

   end subroutine costf2
!------------------------------------------------------------------------------
end module cosine_transform_even
