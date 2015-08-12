module cosine_transform

   use precision_mod, only: cp
   use truncation, only: lm_max
   use fft_fac_mod, only: fft_fac_complex, fft_fac_real
   use const, only: half, one, two

   implicit none

   private

   interface costf1
      module procedure costf1_complex
      module procedure costf1_complex_1d
      module procedure costf1_real
      module procedure costf1_real_1d
   end interface costf1

   public :: costf1, costf2

contains

   subroutine costf1_complex(f,n_f_max,n_f_start,n_f_stop,f2, &
        &                    i_costf_init,d_costf_init)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to perform a multiple              |
      !  |  cosine transforms for n+1 datapoints                             |
      !  |  on the columns numbered n_f_start to n_f_stop in the array       |
      !  |   f(n_f_max,n+1)                                                  |
      !  |  Depending whether the input f contains data or coeff arrays      |
      !  |  coeffs or data are returned in f.                                |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
      integer,  intent(in) :: i_costf_init(*)    ! prestored integers
      real(cp), intent(in) :: d_costf_init(*)    ! prestored dble numbers
    
      !-- Output variables:
      complex(cp), intent(inout) :: f(n_f_max,*)   ! data/coeff input
      complex(cp), intent(out) :: f2(n_f_max,*)  ! work array of the same size as f
    
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
    
      n=i_costf_init(1)-1
    
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
      j1=i_costf_init(2)
      j2=i_costf_init(n_O2_P2)
      do n_f=n_f_start,n_f_stop
         tot_sum(n_f)=f(n_f,1)-f(n_f,n_P1)
         f2(n_f,j1)=f(n_f,1)+f(n_f,n_P1)
         f2(n_f,j2)=two*f(n_f,n_O2_P1)
      end do
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=i_costf_init(j+2)    ! first step
         j2=i_costf_init(n_P2-j)
         wr_j=d_costf_init(j)
         wi_j=d_costf_init(n_O2-j)
    
         i=j+1
         i1=i_costf_init(i+2)   ! second step
         i2=i_costf_init(n_P2-i)
         wr_i=d_costf_init(i)
         wi_i=d_costf_init(n_O2-i)

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
      j1=i_costf_init(j+2)
      j2=i_costf_init(n_P2-j)
      wr_j=d_costf_init(j)
      wi_j=d_costf_init(n_O2-j)
    
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
    
      n_factors=i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=i_costf_init(n+3+n_fac)
         if ( l_f2_data ) then
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_complex(f2(1,1),f2(1,2),f(1,1),f(1,2),    &
                 &       d_costf_init(n+2),n_f_max,n_f_start,      &
                 &       n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_complex(f(1,1),f(1,2),f2(1,1),f2(1,2),   &
                 &       d_costf_init(n+2),n_f_max,n_f_start,     &
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
    
         wr_j=d_costf_init(n_O2+2*j-1)
         wi_j=d_costf_init(n_O2+2*j)
    
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
   subroutine costf1_complex_1d(f,f2,i_costf_init,d_costf_init)
    
      !-- Input variables:
      integer,  intent(in) :: i_costf_init(*)    ! prestored integers
      real(cp), intent(in) :: d_costf_init(*)    ! prestored dble numbers
    
      !-- Output variables:
      complex(cp), intent(inout) :: f(*)   ! data/coeff input
      complex(cp), intent(out) :: f2(*)    ! work array of the same size as f
    
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
    
      n=i_costf_init(1)-1
    
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
    
      j1=i_costf_init(2)
      j2=i_costf_init(n_O2_P2)
      tot_sum=f(1)-f(n_P1)
      f2(j1)=f(1)+f(n_P1)
      f2(j2)=two*f(n_O2_P1)
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=i_costf_init(j+2)    ! first step
         j2=i_costf_init(n_P2-j)
         wr_j=d_costf_init(j)
         wi_j=d_costf_init(n_O2-j)
    
         i=j+1
         i1=i_costf_init(i+2)   ! second step
         i2=i_costf_init(n_P2-i)
         wr_i=d_costf_init(i)
         wi_i=d_costf_init(n_O2-i)
    
         f_h1=f(j+1)+f(n_P1-j)
         f_h2=f(j+1)-f(n_P1-j)
         f2(j1)=f_h1-wi_j*f_h2
         f2(j2)=f_h1+wi_j*f_h2
         tot_sum=tot_sum+wr_j*f_h2

         f_h1=f(i+1)+f(n_P1-i)
         f_h2=f(i+1)-f(n_P1-i)
         f2(i1)=f_h1-wi_i*f_h2
         f2(i2)=f_h1+wi_i*f_h2
         tot_sum=tot_sum+wr_i*f_h2
      end do
    
      j=n/2-1   ! last step
      j1=i_costf_init(j+2)
      j2=i_costf_init(n_P2-j)
      wr_j=d_costf_init(j)
      wi_j=d_costf_init(n_O2-j)
    
      f_h1=f(j+1)+f(n_P1-j)
      f_h2=f(j+1)-f(n_P1-j)
      f2(j1)=f_h1-wi_j*f_h2
      f2(j2)=f_h1+wi_j*f_h2
      tot_sum=tot_sum+wr_j*f_h2
    
      !-- Perform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2
    
      n_factors=i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=i_costf_init(n+3+n_fac)
    
         if ( l_f2_data ) then
    
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_complex(f2(1),f2(2),f(1),f(2), &
                 &       d_costf_init(n+2),1,1,1,       &
                 &       n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_complex(f(1),f(2),f2(1),f2(2), &
                 &       d_costf_init(n+2),1,1,1,       &
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
            f(j1)=f2(j1)
            f(j2)=f2(j2)
            f(j3)=f2(j3)
            f(j4)=f2(j4)
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
    
         wr_j=d_costf_init(n_O2+2*j-1)
         wi_j=d_costf_init(n_O2+2*j)
    
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
   subroutine costf1_real(f,n_f_max,n_f_start,n_f_stop,f2, &
        &                 i_costf_init,d_costf_init)

      use truncation,  only: lm_max_real
    
      implicit none
    
      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
      integer,  intent(in) :: i_costf_init(*)    ! prestored integers
      real(cp), intent(in) :: d_costf_init(*)    ! prestored dble numbers
    
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
    
      n=i_costf_init(1)-1
    
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
      j1=i_costf_init(2)
      j2=i_costf_init(n_O2_P2)
      do n_f=n_f_start,n_f_stop
         tot_sum(n_f)=f(n_f,1)-f(n_f,n_P1)
         f2(n_f,j1)=f(n_f,1)+f(n_f,n_P1)
         f2(n_f,j2)=two*f(n_f,n_O2_P1)
      end do
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=i_costf_init(j+2)    ! first step
         j2=i_costf_init(n_P2-j)
         wr_j=d_costf_init(j)
         wi_j=d_costf_init(n_O2-j)
    
         i=j+1
         i1=i_costf_init(i+2)   ! second step
         i2=i_costf_init(n_P2-i)
         wr_i=d_costf_init(i)
         wi_i=d_costf_init(n_O2-i)
    
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
      j1=i_costf_init(j+2)
      j2=i_costf_init(n_P2-j)
      wr_j=d_costf_init(j)
      wi_j=d_costf_init(n_O2-j)
    
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
    
      n_factors=i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=i_costf_init(n+3+n_fac)
         if ( l_f2_data ) then
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_real(f2(1,1),f2(1,2),f(1,1),f(1,2),  &
                 &       d_costf_init(n+2),n_f_max,n_f_start, &
                 &       n_f_stop,n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_real(f(1,1),f(1,2),f2(1,1),f2(1,2),  &
                 &       d_costf_init(n+2),n_f_max,n_f_start, &
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
    
         wr_j=d_costf_init(n_O2+2*j-1)
         wi_j=d_costf_init(n_O2+2*j)
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
   subroutine costf1_real_1d(f,f2,i_costf_init,d_costf_init)
    
      !-- Input variables:
      integer,  intent(in) :: i_costf_init(*)    ! prestored integers
      real(cp), intent(in) :: d_costf_init(*)    ! prestored dble numbers
    
      !-- Output variables:
      real(cp), intent(inout) :: f(*)   ! data/coeff input
      real(cp), intent(out) :: f2(*)    ! work array of the same size as f
    
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
    
      n=i_costf_init(1)-1
    
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
    
      j1=i_costf_init(2)
      j2=i_costf_init(n_O2_P2)
      tot_sum=f(1)-f(n_P1)
      f2(j1)=f(1)+f(n_P1)
      f2(j2)=two*f(n_O2_P1)
    
      do j=1,n/2-3,2    ! step 2 unrolling
         j1=i_costf_init(j+2)    ! first step
         j2=i_costf_init(n_P2-j)
         wr_j=d_costf_init(j)
         wi_j=d_costf_init(n_O2-j)
    
         i=j+1
         i1=i_costf_init(i+2)   ! second step
         i2=i_costf_init(n_P2-i)
         wr_i=d_costf_init(i)
         wi_i=d_costf_init(n_O2-i)
    
         f_h1=f(j+1)+f(n_P1-j)
         f_h2=f(j+1)-f(n_P1-j)
         f2(j1)=f_h1-wi_j*f_h2
         f2(j2)=f_h1+wi_j*f_h2
         tot_sum=tot_sum+wr_j*f_h2

         f_h1=f(i+1)+f(n_P1-i)
         f_h2=f(i+1)-f(n_P1-i)
         f2(i1)=f_h1-wi_i*f_h2
         f2(i2)=f_h1+wi_i*f_h2
         tot_sum=tot_sum+wr_i*f_h2
      end do
    
      j=n/2-1   ! last step
      j1=i_costf_init(j+2)
      j2=i_costf_init(n_P2-j)
      wr_j=d_costf_init(j)
      wi_j=d_costf_init(n_O2-j)
    
      f_h1=f(j+1)+f(n_P1-j)
      f_h2=f(j+1)-f(n_P1-j)
      f2(j1)=f_h1-wi_j*f_h2
      f2(j2)=f_h1+wi_j*f_h2
      tot_sum=tot_sum+wr_j*f_h2
    
      !-- Perform transform for n_fac factors:
      fac_tot=1              ! total factor
      l_f2_data=.true.   ! data are on f2
    
      n_factors=i_costf_init(n+3)
    
      do n_fac=1,n_factors   ! loop over factors
         fac=i_costf_init(n+3+n_fac)
    
         if ( l_f2_data ) then
    
            !----- fft_fac returns complex transform of f2's on f's:
            call fft_fac_real(f2(1),f2(2),f(1),f(2),  &
                 &       d_costf_init(n+2),1,1,1,     &
                 &       n_O2,fac,fac_tot)
            l_f2_data=.false.
         else
            !----- fft_fac returns complex transform of f's on f2's:
            call fft_fac_real(f(1),f(2),f2(1),f2(2), &
                 &       d_costf_init(n+2),1,1,1,    &
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
            f(j1)=f2(j1)
            f(j2)=f2(j2)
            f(j3)=f2(j3)
            f(j4)=f2(j4)
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
    
         wr_j=d_costf_init(n_O2+2*j-1)
         wi_j=d_costf_init(n_O2+2*j)
    
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
   subroutine costf2(f,n_f_max,n_f_start,n_f_stop,f2, &
                     i_costf_init,d_costf_init,isign)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to perform a multiple              |
      !  |  cosine transforms for n+1 datapoints                             |
      !  |  on the columns numbered n_f_start to n_f_stop in the array       |
      !  |   y(n_f_max,n+1)                                                  |
      !  |  Depending whether the input y contains data or coeff arrays      |
      !  |  coeffs or data are returned in y.                                |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in y,y2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
      integer,  intent(in) :: i_costf_init(*)    ! prestored integers
      real(cp), intent(in) :: d_costf_init(*)    ! prestored dble numbers
      integer,  intent(in) :: isign     !  = +1 (-1) for forward (backward) transform

      !-- Output variables
      complex(cp), intent(inout) :: f(n_f_max,*)   ! data/coeff input
      complex(cp), intent(out) :: f2(n_f_max,*)    ! work array of the same size as y

      !-- Local variables:
      logical :: l_f2_data
      integer :: n, n_f
      integer :: i,j,j1,j2,j3,j4
      integer :: k1,k2,k3,k4
      integer :: n_P1  ! n+1
      integer :: n_P2  ! n+2
      integer :: n_P3  ! n+3
      integer :: n_O2  ! n/2
      integer :: n_O2_P1 ! n/2+1
      integer :: n_O2_P2 ! n/2+2
      integer :: n_O2_P3 ! n/2+3

      complex(cp) :: f_h1,f_h2,f_h3,f_h4 ! help variables
      complex(cp) :: w_h1,w_h2
      real(cp) :: fac_norm,facn
      real(cp) :: wr_j,wi_j,wr_i,wi_i

      integer :: n_fac,fac,fac_tot,n_factors

      complex(cp) :: sum(lm_max)

      n=i_costf_init(1)

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
      if ( isign == 1 ) then

         do j=1,n/2-1,2    ! step 2 unrolling

            j1=i_costf_init(j+1)    ! first step
            j2=i_costf_init(n_P2-j)
            wi_j=d_costf_init(2*n+j)

            i=j+1
            j3=i_costf_init(i+1)   ! second step
            j4=i_costf_init(n_P2-i)
            wi_i=d_costf_init(2*n+i)
             
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

         n_factors=i_costf_init(n+2)

         do n_fac=1,n_factors   ! loop over factors
            fac=i_costf_init(n+2+n_fac)
            if ( l_f2_data ) then
               !-- Vpassm returns complex transform of f2's on f's:
               call fft_fac_complex(f2(1,1),f2(1,2),f(1,1),f(1,2),  &
                            d_costf_init(n+1),n_f_max,              &
                            n_f_start,n_f_stop,n_O2,fac,fac_tot)
               l_f2_data=.false.
            else
               !-- Vpassm returns complex transform of f's on f2's:
               call fft_fac_complex(f(1,1),f(1,2),f2(1,1),f2(1,2),  &
                            d_costf_init(n+1),n_f_max,              &
                            n_f_start,n_f_stop,n_O2,fac,fac_tot)
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

            wr_j=d_costf_init(n_O2+2*j-1)
            wi_j=d_costf_init(n_O2+2*j)

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

            wr_i=d_costf_init(i)
            wi_i=d_costf_init(n_O2-i)

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

      else if ( isign == -1 ) then  ! Inverse transform:

         !-- Calculation of auxiliary function:
         !----- Save f(n):
         do n_f=n_f_start,n_f_stop
            sum(n_f)=f(n_f,n)          ! save f(n)
         end do

         do j1=n,4,-2
            j2=j1-2
            do n_f=n_f_start,n_f_stop
               f(n_f,j1)=f(n_f,j2)-f(n_f,j1)  ! Difference of odd terms
            end do
         end do

         do n_f=n_f_start,n_f_stop
            f(n_f,2)=two*sum(n_f)      ! Write saved f(n) to f(2)
         end do

         do j1=3,n,2
            j2=j1+1
            i=(j1-1)/2

            wr_i=d_costf_init(i)
            wi_i=d_costf_init(n_O2-i)

            do n_f=n_f_start,n_f_stop
               f_h1=f(n_f,j1)*wr_i + f(n_f,j2)*wi_i
               f_h2=f(n_f,j2)*wr_i - f(n_f,j1)*wi_i
               f(n_f,j1)=f_h1
               f(n_f,j2)=f_h2
            end do
         end do


         !-- Preprocessing for realtf, copying on f2:
         k1=i_costf_init(2)
         k2=i_costf_init(3)
         k3=i_costf_init(n_O2_P2)
         k4=i_costf_init(n_O2_P3)
         facn=two*fac_norm
         do n_f=n_f_start,n_f_stop
            f_h1=f(n_f,1)
            f2(n_f,k1)=fac_norm*(f(n_f,1)+f(n_f,2))
            f2(n_f,k2)=fac_norm*(f_h1-f(n_f,2))
            f2(n_f,k3)=facn*f(n_f,n_O2_P1)
            f2(n_f,k4)=facn*f(n_f,n_O2_P2)
         end do

         do j=2,n/4
            j2=2*j
            j1=j2-1
            j3=n_P3-j2
            j4=j3+1

            wr_j=d_costf_init(n_O2+2*j-1)
            wi_j=d_costf_init(n_O2+2*j)

            k1=i_costf_init(j1+1)
            k2=i_costf_init(j2+1)
            k3=i_costf_init(j3+1)
            k4=i_costf_init(j4+1)

            do n_f=n_f_start,n_f_stop
               f_h1=fac_norm*(f(n_f,j1)+f(n_f,j3))
               f_h2=fac_norm*(f(n_f,j2)-f(n_f,j4))
               f_h3=-fac_norm*(f(n_f,j2)+f(n_f,j4))
               f_h4=fac_norm*(f(n_f,j1)-f(n_f,j3))

               w_h1=-wr_j*f_h4+wi_j*f_h3
               w_h2= wr_j*f_h3+wi_j*f_h4

               f2(n_f,k1)= f_h1+w_h2
               f2(n_f,k2)= f_h2-w_h1
               f2(n_f,k3)= f_h1-w_h2
               f2(n_f,k4)=-f_h2-w_h1
            end do
         end do

         !-- Perform transform for n_fac factors:
         fac_tot=1              ! total factor
         l_f2_data=.true.       ! data are on f2

         n_factors=i_costf_init(n+2)

         do n_fac=1,n_factors   ! loop over factors
            fac=i_costf_init(n+2+n_fac)

            if ( l_f2_data ) then
               !-- Vpassm returns complex transform of f2's on f's:
               call fft_fac_complex(f2(1,1),f2(1,2),f(1,1),f(1,2),  &
                            d_costf_init(n+1),n_f_max,n_f_start,    &
                            n_f_stop,n_O2,fac,fac_tot)
               l_f2_data=.false.
            else
               !-- Vpassm returns complex transform of f's on f2's:
               call fft_fac_complex(f(1,1),f(1,2),f2(1,1),f2(1,2),  &
                            d_costf_init(n+1),n_f_max,n_f_start,    &
                            n_f_stop,n_O2,fac,fac_tot)
               l_f2_data=.true.
            end if
            fac_tot=fac_tot*fac
         end do

         if ( l_f2_data ) then
            !-- Copy data on f2:
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

         !------- Extract auxiliary function for costf using recurrence:
         do j1=1,n_O2
            j2=n_P1-j1
            wi_j=d_costf_init(2*n+j1)
            do n_f=n_f_start,n_f_stop
               f_h1=f(n_f,j1)+f(n_f,j2)
               f_h2=(f(n_f,j1)-f(n_f,j2))/wi_j
               f(n_f,j1)=f_h1+f_h2
               f(n_f,j2)=f_h1-f_h2
            end do
         end do

      end if
        
   end subroutine costf2
!------------------------------------------------------------------------------
end module cosine_transform
