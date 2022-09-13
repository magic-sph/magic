module cosine_transform_gpu
   !
   ! This module contains the Hipfort_hipfftw wrappers for the discrete Cosine Transforms
   ! (DCT-I).
   !

#ifdef WITH_OMP_GPU

   use iso_c_binding
   use precision_mod
   use constants, only: half, pi, one, two, zero
   use blocking, only: llm, ulm
   use omp_lib
   use hipfort_hipfft
   use hipfort_check

   implicit none

   private

   type, public :: gpu_costf_odd_t
      integer :: n_r_max                ! Number of radial grid points
      real(cp) :: cheb_fac              ! Normalisation factor
      type(c_ptr) :: plan_many          ! FFTW many plan
      type(c_ptr) :: plan_1d            ! FFTW single plan for DCT
      integer :: plan_1d_size           ! Size of plan_1d
      integer :: plan_many_size(1)      ! Size of plan_many
   contains
      procedure :: initialize
      procedure :: finalize
      procedure, private :: costf1_complex
      procedure, private :: costf1_real_1d
      procedure, private :: costf1_complex_1d
      generic :: costf1 => costf1_real_1d, costf1_complex, costf1_complex_1d
   end type gpu_costf_odd_t

contains

   subroutine initialize(this, n_r_max, n_in, n_in2)
      !
      ! Definition of FFTW plans for type I DCTs.
      !

      class(gpu_costf_odd_t) :: this

      !-- Input variables
      integer, intent(in) :: n_in    ! Not used here, only for compatibility
      integer, intent(in) :: n_in2   ! Not used here, only for compatibility
      integer, intent(in) :: n_r_max ! Number of radial grid points
      integer ::start_lm, stop_lm
      integer :: inembed(1), istride, idist
      integer :: onembed(1), ostride, odist, howmany

      this%n_r_max = n_r_max
      this%plan_1d_size = 2*(n_r_max-1)
      this%cheb_fac = sqrt(half/(n_r_max-1))

      !-- Create a single 1d plan for HIPFFT_Z2Z
      call hipfftCheck(hipfftPlan1d(this%plan_1d, this%plan_1d_size, HIPFFT_Z2Z, 1))

      !-- Create a many plan for many Z2Z : plan_many
      start_lm = llm
      stop_lm  = ulm
      this%plan_many_size(1) = 2*(n_r_max-1)
      idist = 1
      odist = 1
      inembed = this%plan_many_size(1)
      onembed = this%plan_many_size(1)
      howmany = stop_lm-start_lm+1
      istride = stop_lm-start_lm+1
      ostride = stop_lm-start_lm+1
      call hipfftCheck(hipfftPlanMany(this%plan_many, 1, c_loc(this%plan_many_size), c_loc(inembed), istride, &
                                    & idist, c_loc(onembed), ostride, odist, HIPFFT_Z2Z, howmany))

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Desctruction of FFTW plans for DCT-I
      !

      class(gpu_costf_odd_t) :: this

      call hipfftcheck( hipfftDestroy(this%plan_1d))
      call hipfftcheck( hipfftDestroy(this%plan_many))

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this, array_in, work_1d)
      !
      ! DCT for one single real vector of dimension ``n_r_max``
      !

      class(gpu_costf_odd_t), intent(in) :: this

      !-- Output variables:
      real(cp), intent(inout) :: array_in(:) ! data to be transformed
      real(cp), intent(out) :: work_1d(:)    ! work array (not used)

      !-- Local variables:
      complex(cp), allocatable, target :: temp_in(:), temp_out(:)
      real(cp), allocatable :: temp_resu(:)

      !-- Allocate
      allocate(temp_in(this%plan_1d_size), temp_out(this%plan_1d_size))
      allocate(temp_resu(this%plan_1d_size))

      !--
      !-- https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-I
      !-- http://www.fftw.org/fftw3_doc/Real-even_002fodd-DFTs-_0028cosine_002fsine-transforms_0029.html
      temp_in(:) = zero
      temp_in(1:this%n_r_max) = array_in(1:this%n_r_max)
      temp_in(this%n_r_max+1:this%plan_1d_size) = array_in(this%n_r_max-1:2:-1)
      temp_out(:) = zero

      !--
      !$omp target enter data map(alloc : temp_in, temp_out)
      !$omp target update to(temp_in, temp_out)
      !$omp target data use_device_addr(temp_in, temp_out)
      call hipfftCheck(hipfftExecZ2Z(this%plan_1d, c_loc(temp_in), c_loc(temp_out), HIPFFT_FORWARD))
      !$omp end target data
      !$omp target update from(temp_out)
      !$omp target exit data map(delete : temp_in, temp_out)

      !--
      temp_resu(:) = real(temp_out(:))

      !--
      work_1d(1:this%n_r_max) = temp_resu(1:this%n_r_max)
      array_in(:) = this%cheb_fac * work_1d(:)

      !--
      deallocate(temp_in, temp_out)
      deallocate(temp_resu)

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
   subroutine costf1_complex_1d(this, array_in, work_1d)
      !
      ! DCT for one single complex vector of dimension ``n_r_max``
      !

      class(gpu_costf_odd_t), intent(in) :: this

      !-- Output variables:
      complex(cp), intent(inout) :: array_in(:) ! data to be transformed
      complex(cp), intent(out) :: work_1d(:)    ! work array (not needed)

      !-- local variables
      complex(cp), allocatable, target :: real_(:), aimag_(:), out_real(:), out_aim(:)
      real(cp), allocatable :: resu_rr(:), resu_ir(:)

      !-- Allocate
      allocate(real_(this%plan_1d_size), aimag_(this%plan_1d_size), out_real(this%plan_1d_size), out_aim(this%plan_1d_size))
      allocate(resu_rr(this%n_r_max), resu_ir(this%n_r_max))

      !--
      real_(:)  = zero
      real_(1:this%n_r_max)                    = cmplx(real(array_in(1:this%n_r_max)), 0.0_cp, kind=cp)
      real_(this%n_r_max+1:this%plan_1d_size)  = cmplx(real(array_in(this%n_r_max-1:2:-1)), 0.0_cp, kind=cp)
      aimag_(:) = zero
      aimag_(1:this%n_r_max)                   = cmplx(aimag(array_in(1:this%n_r_max)), 0.0_cp, kind=cp)
      aimag_(this%n_r_max+1:this%plan_1d_size) = cmplx(aimag(array_in(this%n_r_max-1:2:-1)), 0.0_cp, kind=cp)

      !--
      out_real(:) = zero
      out_aim(:)  = zero
      !$omp target enter data map(alloc : real_, aimag_, out_real, out_aim)
      !$omp target update to(real_, aimag_, out_real, out_aim)
      !$omp target data use_device_addr(real_, aimag_, out_real, out_aim)
      call hipfftCheck(hipfftExecZ2Z(this%plan_1d, c_loc(real_), c_loc(out_real), HIPFFT_FORWARD))
      call hipfftCheck(hipfftExecZ2Z(this%plan_1d, c_loc(aimag_), c_loc(out_aim), HIPFFT_FORWARD))
      !$omp end target data
      !$omp target update from(out_real, out_aim)
      !$omp target exit data map(delete : real_, aimag_, out_real, out_aim)

      resu_rr(:) = real(out_real(1:this%n_r_max))
      resu_ir(:) = real(out_aim(1:this%n_r_max))

      array_in(:) = this%cheb_fac * cmplx(resu_rr(:), resu_ir(:), cp)

      !--
      deallocate(real_, aimag_, out_real, out_aim, resu_rr, resu_ir)

   end subroutine costf1_complex_1d
!------------------------------------------------------------------------------
   subroutine costf1_complex(this, array_in, n_f_max, n_f_start, n_f_stop, work_2d)
      !
      ! Multiple DCT's for 2-D complex input array.
      !

      class(gpu_costf_odd_t), intent(in) :: this

      !-- Input variables
      integer, intent(in) :: n_f_start ! Starting index
      integer, intent(in) :: n_f_stop  ! Stopping index
      integer, intent(in) :: n_f_max   ! Number of vectors

      !-- Output variables:
      complex(cp), intent(inout) :: array_in(:,:) ! Array to be transformed
      complex(cp), intent(inout) :: work_2d(n_f_max,*)  ! Help array (not needed)
      !-- ftn-7032 ftn: ERROR COSTF1_COMPLEX, File = ../../magic/src/dct_hipfort_hipfftw.f90, Line = 226, 231, 243
      !-- Unsupported OpenMP construct Assumed size arrays -- array_in

      !-- Local variables:
      complex(cp), allocatable, target :: tmp_in(:,:), tmp_out(:,:)
      integer :: n_r, n_f

      !-- Allocate temp arrays
      allocate(tmp_in(n_f_start:n_f_stop,2*this%n_r_max-2), tmp_out(n_f_start:n_f_stop,2*this%n_r_max-2))
      !$omp target enter data map(alloc : tmp_in, tmp_out)

      !-- initialize to zero
      !$omp target teams distribute parallel do collapse(2)
      do n_f=n_f_start,n_f_stop
         do n_r=1,2*this%n_r_max-2
            tmp_in(n_f, n_r) = zero
            tmp_out(n_f, n_r) = zero
         end do
      end do
      !$omp end target teams distribute parallel do

      !-- Prepare array for dft many
      !$omp target teams distribute parallel do
      do n_r=1,this%n_r_max
         tmp_in(:,n_r)=array_in(n_f_start:n_f_stop,n_r)
      end do
      !$omp end target teams distribute parallel do
      !$omp target teams distribute parallel do
      do n_r=this%n_r_max+1,2*this%n_r_max-2
         tmp_in(:,n_r)=array_in(n_f_start:n_f_stop,2*this%n_r_max-n_r)
      end do
      !$omp end target teams distribute parallel do

      !-- Perform DFT many
      !$omp target data use_device_addr(tmp_in, tmp_out)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_in), c_loc(tmp_out), HIPFFT_FORWARD))
      !$omp end target data

      !-- Copy output onto array_in
      !$omp target teams distribute parallel do
      do n_r=1,this%n_r_max
         array_in(n_f_start:n_f_stop,n_r)=this%cheb_fac*tmp_out(:,n_r)
      end do
      !$omp end target teams distribute parallel do

      !-- Free allocated arrays
      !$omp target exit data map(delete : tmp_in, tmp_out)
      deallocate(tmp_in, tmp_out)

   end subroutine costf1_complex
!------------------------------------------------------------------------------

#endif

end module cosine_transform_gpu
