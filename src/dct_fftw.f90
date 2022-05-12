#define dct_many 0
#define dct_loop 1
#define dft_loop 2
#define DCT_VERSION dft_loop
module cosine_transform_odd
   !
   ! This module contains the FFTW wrappers for the discrete Cosine Transforms
   ! (DCT-I). Unfortunately, the MKL has no support for the many_r2r variants
   ! such that one has to either manually loop over lm's and treat the real and the
   ! imaginary parts separately or use DFTs instead. Both strategies still seems
   ! to outperform the built-in transforms but this is not always the case. The
   ! latter approach using DFTs seems to offer the best performance.
   !

   use iso_c_binding
   use blocking, only: llm, ulm
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: half, pi, one, two
#ifdef WITHOMP
   use omp_lib
#endif

   implicit none

   include 'fftw3.f03'

   private

   !-- For type-I DCT, FFTW_EXHAUSTIVE yields a speed-up
   integer(c_int), parameter :: fftw_plan_flag=FFTW_EXHAUSTIVE
   integer(c_int), parameter :: fft_plan_flag=FFTW_PATIENT

   type, public :: costf_odd_t
      integer :: n_r_max                ! Number of radial grid points
      real(cp) :: cheb_fac              ! Normalisation factor
      type(c_ptr) :: plan               ! FFTW many plan
      type(c_ptr) :: plan_1d            ! FFTW single plan for DCT
      type(c_ptr) :: plan_fft_1d        ! FFTW single plan for FFT
      complex(cp), pointer :: work(:,:) ! Complex work array
      real(cp), pointer :: work_r(:,:)  ! Real work array
   contains
      procedure :: initialize
      procedure :: finalize
      procedure, private :: costf1_complex
      procedure, private :: costf1_real_1d
      procedure, private :: costf1_complex_1d
      generic :: costf1 => costf1_real_1d, costf1_complex, costf1_complex_1d
   end type costf_odd_t

contains

   subroutine initialize(this, n_r_max, n_in, n_in2)
      !
      ! Definition of FFTW plans for type I DCTs. 
      !

      class(costf_odd_t) :: this
      
      !-- Input variables
      integer, intent(in) :: n_in    ! Not used here, only for compatibility
      integer, intent(in) :: n_in2   ! Not used here, only for compatibility
      integer, intent(in) :: n_r_max ! Number of radial grid points

      !--Local variables
#ifdef WITHOMP
      integer :: ier
#endif
      integer :: inembed(1), istride, idist, plan_size(1)
      integer :: onembed(1), ostride, odist, isize, howmany
      integer(C_INT) :: plan_type(1)
#if (DCT_VERSION==dft_loop)
      complex(cp) :: array_cplx_1d(2*n_r_max-2), array_cplx_out_1d(2*n_r_max-2)
#elif (DCT_VERSION==dct_many)
      real(cp) :: array_in(2*(ulm-llm+1),n_r_max), array_out(2*(ulm-llm+1),n_r_max)
#endif
      real(cp) :: array_in_1d(n_r_max), array_out_1d(n_r_max)

#ifdef WITHOMP
      ier =  fftw_init_threads()
      call fftw_plan_with_nthreads(1) ! No OMP for those plans
#endif

      this%n_r_max = n_r_max
      plan_type(1) = FFTW_REDFT00

#if (DCT_VERSION==dct_many)
      plan_size(1) = n_r_max
      howmany = 2*(ulm-llm+1)
      inembed(1) = 0
      onembed(1) = 0
      istride = 2*(ulm-llm+1)
      ostride = 2*(ulm-llm+1)
      isize   = 2*(ulm-llm+1)
      idist = 1
      odist = 1

      this%plan = fftw_plan_many_r2r(1, plan_size, isize, array_in,         &
                  &                  inembed, istride, idist, array_out,    &
                  &                  onembed, ostride, odist,               &
                  &                  plan_type, fftw_plan_flag)

      allocate( this%work(1:ulm-llm+1,n_r_max) )
      call c_f_pointer(c_loc(this%work), this%work_r, [isize, n_r_max])

      bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max* &
      &                 SIZEOF_DEF_COMPLEX
#elif (DCT_VERSION==dft_loop)
      plan_size(1) = 2*n_r_max-2
      this%plan_fft_1d = fftw_plan_dft(1, plan_size, array_cplx_1d, array_cplx_out_1d, &
                         &             1, fft_plan_flag)
#endif
      
      plan_size(1) = n_r_max
      this%plan_1d = fftw_plan_r2r(1, plan_size, array_in_1d, array_out_1d, &
                     &             plan_type, fftw_plan_flag)


      this%cheb_fac = sqrt(half/(n_r_max-1))

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Desctruction of FFTW plans for DCT-I and deallocation of help arrays
      !

      class(costf_odd_t) :: this

#if (DCT_VERSION==dct_many)
      deallocate( this%work )
      call fftw_destroy_plan(this%plan)
#elif (DCT_VERSION==dft_loop)
      call fftw_destroy_plan(this%plan_fft_1d)
#endif
#ifdef WITHOMP
      !call fftw_cleanup_threads()
#endif
      call fftw_destroy_plan(this%plan_1d)

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine costf1_complex(this, array_in, n_f_max, n_f_start, n_f_stop, work_2d)
      !
      ! Multiple DCT's for 2-D complex input array.
      !

      class(costf_odd_t), intent(in) :: this

      !-- Input variables
      integer, intent(in) :: n_f_start ! Starting index (OMP)
      integer, intent(in) :: n_f_stop  ! Stopping index (OMP)
      integer, intent(in) :: n_f_max   ! Number of vectors

      !-- Output variables:
      complex(cp), intent(inout) :: array_in(n_f_max,*) ! Array to be transformed
      complex(cp), intent(inout) :: work_2d(n_f_max,*)  ! Help array (not needed)

#if (DCT_VERSION==dft_loop)
      !-- Local variables:
      integer :: n_f
      complex(cp) :: work_1d(2*this%n_r_max-2), work_1d_out(2*this%n_r_max-2)

      do n_f=n_f_start,n_f_stop
         work_1d(1:this%n_r_max) = array_in(n_f,1:this%n_r_max)
         work_1d(this%n_r_max+1:) = array_in(n_f,this%n_r_max-1:2:-1)
         call fftw_execute_dft(this%plan_fft_1d, work_1d, work_1d_out)
         array_in(n_f,1:this%n_r_max)=this%cheb_fac* work_1d_out(1:this%n_r_max)
      end do

#elif (DCT_VERSION==dct_loop)
      !-- Local variables:
      integer :: n_f
      real(cp) :: r_input(this%n_r_max), i_input(this%n_r_max), work_1d(this%n_r_max)

      !- Uncomment in case OpenMP is moved inwards
      !!$omp parallel do default(shared) private(n_f,r_input,work_1d,i_input)
      do n_f=n_f_start,n_f_stop
         work_1d(:) = real(array_in(n_f,1:this%n_r_max))
         call fftw_execute_r2r(this%plan_1d, work_1d, r_input)
         work_1d(:) = aimag(array_in(n_f,1:this%n_r_max))
         call fftw_execute_r2r(this%plan_1d, work_1d, i_input)
         array_in(n_f,1:this%n_r_max)=this%cheb_fac*cmplx(r_input, i_input, kind=cp)
      end do
      !!$omp end parallel do

#elif (DCT_VERSION==dct_many)
      ! This should be the fastest but unfortunately MKL has no support for it:
      !https://software.intel.com/content/www/us/en/develop/documentation/
      !mkl-developer-reference-c/top/
      !appendix-d-fftw-interface-to-intel-math-kernel-library/
      !fftw3-interface-to-intel-math-kernel-library/using-fftw3-wrappers.html
      !-- Local variables
      real(cp), pointer :: r_input(:,:)
      integer :: n_r

      call c_f_pointer(c_loc(array_in), r_input, [2*n_f_max, this%n_r_max])
      call fftw_execute_r2r(this%plan, r_input, this%work_r)

      !$omp parallel do
      do n_r=1,this%n_r_max
         array_in(:,n_r)=this%cheb_fac*this%work(:,n_r)
      end do
      !$omp end parallel do
#endif

   end subroutine costf1_complex
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this, array_in, work_1d)
      !
      ! DCT for one single real vector of dimension ``n_r_max``
      !

      class(costf_odd_t), intent(in) :: this

      !-- Output variables:
      real(cp), intent(inout) :: array_in(:) ! data to be transformed
      real(cp), intent(out) :: work_1d(:)    ! work array (not used)

      call fftw_execute_r2r(this%plan_1d, array_in, work_1d)
      array_in(:)=this%cheb_fac*work_1d(:)

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
   subroutine costf1_complex_1d(this, array_in, work_1d)
      !
      ! DCT for one single complex vector of dimension ``n_r_max``
      !

      class(costf_odd_t), intent(in) :: this

      !-- Output variables:
      complex(cp), intent(inout) :: array_in(:) ! data to be transformed
      complex(cp), intent(out) :: work_1d(:)    ! work array (not needed)

      !-- Local variables:
      real(cp) :: tmpr(this%n_r_max), tmpi(this%n_r_max)
      real(cp) :: outr(this%n_r_max), outi(this%n_r_max)

      tmpr(:)= real(array_in(:))
      tmpi(:)=aimag(array_in(:))

      call fftw_execute_r2r(this%plan_1d, tmpr, outr)
      call fftw_execute_r2r(this%plan_1d, tmpi, outi)

      array_in(:)=this%cheb_fac*cmplx(outr(:), outi(:), cp)

   end subroutine costf1_complex_1d
!------------------------------------------------------------------------------
end module cosine_transform_odd
