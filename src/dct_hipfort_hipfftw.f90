module cosine_transform_gpu
   !
   ! This module contains the Hipfort_hipfftw wrappers for the discrete Cosine Transforms
   ! (DCT-I).
   !

#ifdef WITH_OMP_GPU

   use iso_c_binding
   use precision_mod
   use constants, only: half, pi, one, two, zero, ci, third, three
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
      procedure :: get_dr_fft
      procedure :: get_ddr_fft
      procedure :: get_dddr_fft
      generic :: costf1 => costf1_real_1d, costf1_complex, costf1_complex_1d
   end type gpu_costf_odd_t

   complex(cp), allocatable, target :: tmp_in(:,:), tmp_out(:,:)
   complex(cp), allocatable, target :: work_2(:,:), work_3(:,:)
   integer, allocatable :: der(:), der2(:)

contains

   subroutine initialize(this, n_r_max, n_cheb_max, n_in)
      !
      ! Definition of FFTW plans for type I DCTs.
      !

      class(gpu_costf_odd_t) :: this

      !-- Input variables
      integer, intent(in) :: n_cheb_max  ! Max number of Chebyshev modes
      integer, intent(in) :: n_in    ! Not used here, only for compatibility
      integer, intent(in) :: n_r_max ! Number of radial grid points
      integer ::start_lm, stop_lm
      integer :: inembed(1), istride, idist, j, k
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
      call hipfftCheck(hipfftPlanMany(this%plan_many, 1, c_loc(this%plan_many_size), &
           &                          c_loc(inembed), istride, idist, c_loc(onembed),&
           &                          ostride, odist, HIPFFT_Z2Z, howmany))

      !-- Allocate tmp_* arrays
      allocate(tmp_in(1:ulm-llm+1,2*this%n_r_max-2), tmp_out(1:ulm-llm+1,2*this%n_r_max-2))
      allocate(work_2(1:ulm-llm+1,2*this%n_r_max-2), work_3(1:ulm-llm+1,2*this%n_r_max-2))
      tmp_in(:,:) = zero; tmp_out(:,:) = zero
      work_2(:,:) = zero; work_3(:,:) = zero

      allocate ( der(2*n_r_max-2), der2(2*n_r_max-2) )
      der(:)=0
      !bytes_allocated=bytes_allocated+(2*n_r_max-2)*SIZEOF_INTEGER
      do k=1,n_cheb_max!-1
         der(k)=k-1
      end do
      j=1
      do k=2*n_r_max-2,n_r_max+1,-1
         if ( j < n_cheb_max ) der(k)=-j
         j=j+1
      end do
      der2(:)=der(:)*der(:)
      der(n_r_max)=0

      !$omp target enter data map(alloc : tmp_in,tmp_out,work_2,work_3,der,der2)
      !$omp target update to(tmp_in,tmp_out,der,der2,work_2,work_3)

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Desctruction of FFTW plans for DCT-I
      !

      class(gpu_costf_odd_t) :: this

      call hipfftcheck( hipfftDestroy(this%plan_1d))
      call hipfftcheck( hipfftDestroy(this%plan_many))

      !-- Free tmp_* arrays
      !$omp target exit data map(delete : tmp_in, tmp_out, der, der2, work_2, work_3)
      deallocate(tmp_in, tmp_out, der, der2, work_2, work_3)

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
      call hipfftCheck(hipfftExecZ2Z(this%plan_1d, c_loc(temp_in), c_loc(temp_out), &
           &           HIPFFT_FORWARD))
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
      call hipfftCheck(hipfftExecZ2Z(this%plan_1d, c_loc(real_), c_loc(out_real), &
           &           HIPFFT_FORWARD))
      call hipfftCheck(hipfftExecZ2Z(this%plan_1d, c_loc(aimag_), c_loc(out_aim), &
           &           HIPFFT_FORWARD))
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

      !-- Local variables:
      integer :: n_r, n_f, tmp_n_r_max
      real(cp) :: tmp_fac_cheb

      n_r = 0; n_f = 0; tmp_n_r_max = this%n_r_max

      !-- Prepare array for dft many
      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,2*tmp_n_r_max-2
         do n_f=n_f_start,n_f_stop
            if(n_r <= tmp_n_r_max) then
               tmp_in(n_f,n_r)=array_in(n_f,n_r)
            else
               tmp_in(n_f,n_r)=array_in(n_f,2*tmp_n_r_max-n_r)
            end if
         end do
      end do
      !$omp end target teams distribute parallel do

      !-- Perform DFT many
      !$omp target data use_device_addr(tmp_in, tmp_out)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_in), c_loc(tmp_out), &
           &           HIPFFT_FORWARD))
      !$omp end target data

      !-- Copy output onto array_in
      tmp_fac_cheb = this%cheb_fac
      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,tmp_n_r_max
         do n_f=n_f_start,n_f_stop
            array_in(n_f,n_r)=tmp_fac_cheb*tmp_out(n_f,n_r)
         end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine costf1_complex
!------------------------------------------------------------------------------
   subroutine get_dr_fft(this, f, df, xcheb, n_f_max, n_f_start, &
              &          n_f_stop, n_cheb_max, l_dct_in)

      class(gpu_costf_odd_t), intent(in) :: this

      !-- Input variables
      integer,     intent(in) :: n_f_start ! Starting index (OMP)
      integer,     intent(in) :: n_f_stop  ! Stopping index (OMP)
      integer,     intent(in) :: n_f_max   ! Number of vectors
      integer,     intent(in) :: n_cheb_max  ! Max cheb
      real(cp),    intent(in) :: xcheb(:) ! Gauss-Lobatto grid
      complex(cp), intent(in) :: f(n_f_max,this%n_r_max) ! Array to be transformed
      logical,     intent(in) :: l_dct_in ! Do we need a DCT for the input array?

      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,this%n_r_max)  ! Radial derivative

      !-- Local variables:
      integer :: n_r, n_f, tmp_n_r_max, k
      real(cp) :: tmp_fac_cheb
      complex(cp) :: tot

      n_r = 0; n_f = 0; tmp_n_r_max = this%n_r_max
      tmp_fac_cheb = this%cheb_fac

      if ( l_dct_in ) then
         !-- Prepare array for dft many
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,2*tmp_n_r_max-2
            do n_f=n_f_start,n_f_stop
               if (n_r <= tmp_n_r_max) then
                  tmp_in(n_f,n_r)=f(n_f,n_r)
               else
                  tmp_in(n_f,n_r)=f(n_f,2*tmp_n_r_max-n_r)
               end if
            end do
         end do
         !$omp end target teams distribute parallel do
         !-- Perform DFT many
         !$omp target data use_device_addr(tmp_in, tmp_out)
         call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_in), c_loc(tmp_out), &
              &           HIPFFT_FORWARD))
         !$omp end target data
      else
         !-- Prepare array for dft many
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,2*tmp_n_r_max-2
            do n_f=n_f_start,n_f_stop
               if (n_r <= tmp_n_r_max) then
                  tmp_out(n_f,n_r)=f(n_f,n_r)/tmp_fac_cheb
               else
                  tmp_out(n_f,n_r)=f(n_f,2*tmp_n_r_max-n_r)/tmp_fac_cheb
               end if
            end do
         end do
         !$omp end target teams distribute parallel do
      end if

      !-- Boundary point
      !$omp target teams distribute parallel do private(k,tot)
      do n_f=n_f_start,n_f_stop
         tot=zero
         tmp_out(n_f,tmp_n_r_max)=half*tmp_out(n_f,tmp_n_r_max)
         do k=1,tmp_n_r_max-1
            tot=tot+k**2 * tmp_out(n_f,k+1)
         end do
         df(n_f,1)=tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+(-1)**(k+1)*k**2*tmp_out(n_f,k+1)
         end do
         df(n_f,tmp_n_r_max)=tot/(tmp_n_r_max-1)
         tmp_out(n_f,tmp_n_r_max)=two*tmp_out(n_f,tmp_n_r_max)
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,2*tmp_n_r_max-2
         do n_f=n_f_start,n_f_stop
            tmp_out(n_f,n_r)=ci*der(n_r)*tmp_out(n_f,n_r)
         end do
      end do
      !$omp end target teams distribute parallel do

      !-- Perform DFT many
      !$omp target data use_device_addr(tmp_in, tmp_out)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_out), c_loc(tmp_in), &
           &           HIPFFT_BACKWARD))
      !$omp end target data

      !$omp target teams distribute parallel do collapse(2)
      do n_r=2,tmp_n_r_max-1
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r)=-tmp_in(n_f,n_r)/sqrt(one-xcheb(n_r)**2)/(2*tmp_n_r_max-2)
         end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine get_dr_fft
!------------------------------------------------------------------------------
   subroutine get_ddr_fft(this, f, df, ddf, xcheb, n_f_max, n_f_start, &
              &           n_f_stop, n_cheb_max, l_dct_in)

      class(gpu_costf_odd_t), intent(in) :: this

      !-- Input variables
      integer,     intent(in) :: n_f_start ! Starting index (OMP)
      integer,     intent(in) :: n_f_stop  ! Stopping index (OMP)
      integer,     intent(in) :: n_f_max   ! Number of vectors
      integer,     intent(in) :: n_cheb_max  ! Max cheb
      real(cp),    intent(in) :: xcheb(:) ! Gauss-Lobatto grid
      complex(cp), intent(in) :: f(n_f_max,this%n_r_max) ! Array to be transformed
      logical,     intent(in) :: l_dct_in ! Do we need a DCT for the input array?

      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,this%n_r_max)  ! 1st Radial derivative
      complex(cp), intent(out) :: ddf(n_f_max,this%n_r_max) ! 2nd Radial derivative

      !-- Local variables:
      integer :: n_r, n_f, tmp_n_r_max, k
      real(cp) :: tmp_fac_cheb
      complex(cp) :: tot

      n_r = 0; n_f = 0; tmp_n_r_max = this%n_r_max
      tmp_fac_cheb = this%cheb_fac

      if ( l_dct_in ) then
         !-- Prepare array for dft many
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,2*tmp_n_r_max-2
            do n_f=n_f_start,n_f_stop
               if (n_r <= tmp_n_r_max) then
                  tmp_in(n_f,n_r)=f(n_f,n_r)
               else
                  tmp_in(n_f,n_r)=f(n_f,2*tmp_n_r_max-n_r)
               end if
            end do
         end do
         !$omp end target teams distribute parallel do
         !-- Perform DFT many
         !$omp target data use_device_addr(tmp_in, tmp_out)
         call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_in), c_loc(tmp_out), &
              &           HIPFFT_FORWARD))
         !$omp end target data
      else
         !-- Prepare array for dft many
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,2*tmp_n_r_max-2
            do n_f=n_f_start,n_f_stop
               if (n_r <= tmp_n_r_max) then
                  tmp_out(n_f,n_r)=f(n_f,n_r)/tmp_fac_cheb
               else
                  tmp_out(n_f,n_r)=f(n_f,2*tmp_n_r_max-n_r)/tmp_fac_cheb
               end if
            end do
         end do
         !$omp end target teams distribute parallel do
      end if

      !-- Boundary points
      !$omp target teams distribute parallel do private(k,tot)
      do n_f=n_f_start,n_f_stop
         tot=zero
         tmp_out(n_f,tmp_n_r_max)=half*tmp_out(n_f,tmp_n_r_max)
         do k=1,tmp_n_r_max-1
            tot=tot+k**2 * tmp_out(n_f,k+1)
         end do
         df(n_f,1)=tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+k**2*(k**2-1) * tmp_out(n_f,k+1)
         end do
         ddf(n_f,1)=third * tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+(-1)**(k+1)*k**2*tmp_out(n_f,k+1)
         end do
         df(n_f,tmp_n_r_max)=tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+(-1)**k*k**2*(k**2-1)*tmp_out(n_f,k+1)
         end do
         ddf(n_f,tmp_n_r_max)=third*tot/(tmp_n_r_max-1)
         tmp_out(n_f,tmp_n_r_max)=two*tmp_out(n_f,tmp_n_r_max)
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,2*tmp_n_r_max-2
         do n_f=n_f_start,n_f_stop
            work_2(n_f,n_r)=ci*der(n_r)*tmp_out(n_f,n_r)
            work_3(n_f,n_r)=-der2(n_r)*tmp_out(n_f,n_r)
         end do
      end do
      !$omp end target teams distribute parallel do

      !-- Perform DFT many
      !$omp target data use_device_addr(tmp_in, work_2)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(work_2), c_loc(tmp_in), &
           &           HIPFFT_BACKWARD))
      !$omp end target data

      !$omp target data use_device_addr(work_2, work_3)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(work_3), c_loc(work_2), &
           &           HIPFFT_BACKWARD))
      !$omp end target data

      !$omp target teams distribute parallel do collapse(2)
      do n_r=2,tmp_n_r_max-1
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r) =-tmp_in(n_f,n_r)/sqrt(one-xcheb(n_r)**2)/(2*tmp_n_r_max-2)
            ddf(n_f,n_r)=-tmp_in(n_f,n_r)*xcheb(n_r)/(one-xcheb(n_r)**2)**1.5_cp / &
            &            (2*tmp_n_r_max-2)+work_2(n_f,n_r)/(one-xcheb(n_r)**2)/    &
            &            (2*tmp_n_r_max-2)
         end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine get_ddr_fft
!------------------------------------------------------------------------------
   subroutine get_dddr_fft(this, f, df, ddf, dddf, xcheb, n_f_max, n_f_start, &
              &           n_f_stop, n_cheb_max, l_dct_in)

      class(gpu_costf_odd_t), intent(in) :: this

      !-- Input variables
      integer,     intent(in) :: n_f_start ! Starting index (OMP)
      integer,     intent(in) :: n_f_stop  ! Stopping index (OMP)
      integer,     intent(in) :: n_f_max   ! Number of vectors
      integer,     intent(in) :: n_cheb_max  ! Max cheb
      real(cp),    intent(in) :: xcheb(:) ! Gauss-Lobatto grid
      complex(cp), intent(in) :: f(n_f_max,this%n_r_max) ! Array to be transformed
      logical,     intent(in) :: l_dct_in ! Do we need a DCT for the input array?

      !-- Output variables:
      complex(cp), intent(out) :: df(n_f_max,this%n_r_max)  ! 1st Radial derivative
      complex(cp), intent(out) :: ddf(n_f_max,this%n_r_max) ! 2nd Radial derivative
      complex(cp), intent(out) :: dddf(n_f_max,this%n_r_max) ! 3rd Radial derivative

      !-- Local variables:
      integer :: n_r, n_f, tmp_n_r_max, k
      real(cp) :: tmp_fac_cheb
      complex(cp) :: tot

      n_r = 0; n_f = 0; tmp_n_r_max = this%n_r_max
      tmp_fac_cheb = this%cheb_fac

      if ( l_dct_in ) then
         !-- Prepare array for dft many
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,2*tmp_n_r_max-2
            do n_f=n_f_start,n_f_stop
               if (n_r <= tmp_n_r_max) then
                  tmp_in(n_f,n_r)=f(n_f,n_r)
               else
                  tmp_in(n_f,n_r)=f(n_f,2*tmp_n_r_max-n_r)
               end if
            end do
         end do
         !$omp end target teams distribute parallel do
         !-- Perform DFT many
         !$omp target data use_device_addr(tmp_in, tmp_out)
         call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_in), c_loc(tmp_out), &
              &           HIPFFT_FORWARD))
         !$omp end target data
      else
         !-- Prepare array for dft many
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,2*tmp_n_r_max-2
            do n_f=n_f_start,n_f_stop
               if (n_r <= tmp_n_r_max) then
                  tmp_out(n_f,n_r)=f(n_f,n_r)/tmp_fac_cheb
               else
                  tmp_out(n_f,n_r)=f(n_f,2*tmp_n_r_max-n_r)/tmp_fac_cheb
               end if
            end do
         end do
         !$omp end target teams distribute parallel do
      end if

      !-- Boundary points
      !$omp target teams distribute parallel do private(k,tot)
      do n_f=n_f_start,n_f_stop
         tot=zero
         tmp_out(n_f,tmp_n_r_max)=half*tmp_out(n_f,tmp_n_r_max)
         do k=1,tmp_n_r_max-1
            tot=tot+k**2 * tmp_out(n_f,k+1)
         end do
         df(n_f,1)=tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+k**2*(k**2-1) * tmp_out(n_f,k+1)
         end do
         ddf(n_f,1)=third * tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+k**2*(k**2-1)*(k**2-4) * tmp_out(n_f,k+1)
         end do
         dddf(n_f,1)=1.0_cp/15.0_cp*tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+(-1)**(k+1)*k**2*tmp_out(n_f,k+1)
         end do
         df(n_f,tmp_n_r_max)=tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+(-1)**k*k**2*(k**2-1)*tmp_out(n_f,k+1)
         end do
         ddf(n_f,tmp_n_r_max)=third*tot/(tmp_n_r_max-1)
         tot=zero
         do k=1,tmp_n_r_max-1
            tot=tot+(-1)**(k+1)*k**2*(k**2-1)*(k**2-4)*tmp_out(n_f,k+1)
         end do
         dddf(n_f,tmp_n_r_max)=1.0_cp/15.0_cp*tot/(tmp_n_r_max-1)
         tmp_out(n_f,tmp_n_r_max)=two*tmp_out(n_f,tmp_n_r_max)
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,2*tmp_n_r_max-2
         do n_f=n_f_start,n_f_stop
            work_3(n_f,n_r)=ci*der(n_r)*tmp_out(n_f,n_r)
         end do
      end do
      !$omp end target teams distribute parallel do

      !-- Perform DFT many
      !$omp target data use_device_addr(tmp_in, work_3)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(work_3), c_loc(tmp_in), &
           &           HIPFFT_BACKWARD))
      !$omp end target data

      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,2*tmp_n_r_max-2
         do n_f=n_f_start,n_f_stop
            work_3(n_f,n_r)=-der2(n_r)*tmp_out(n_f,n_r)
         end do
      end do
      !$omp end target teams distribute parallel do

      !$omp target data use_device_addr(work_2, work_3)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(work_3), c_loc(work_2), &
           &           HIPFFT_BACKWARD))
      !$omp end target data

      !$omp target teams distribute parallel do collapse(2)
      do n_r=1,2*tmp_n_r_max-2
         do n_f=n_f_start,n_f_stop
            tmp_out(n_f,n_r)=-ci*der(n_r)**3*tmp_out(n_f,n_r)
         end do
      end do
      !$omp end target teams distribute parallel do

      !$omp target data use_device_addr(tmp_out, work_3)
      call hipfftCheck(hipfftExecZ2Z(this%plan_many, c_loc(tmp_out), c_loc(work_3), &
           &           HIPFFT_BACKWARD))
      !$omp end target data

      !$omp target teams distribute parallel do collapse(2)
      do n_r=2,tmp_n_r_max-1
         do n_f=n_f_start,n_f_stop
            df(n_f,n_r)  =-tmp_in(n_f,n_r)/sqrt(one-xcheb(n_r)**2)/(2*tmp_n_r_max-2)
            ddf(n_f,n_r) =-tmp_in(n_f,n_r)*xcheb(n_r)/(one-xcheb(n_r)**2)**1.5_cp / &
            &             (2*tmp_n_r_max-2)+work_2(n_f,n_r)/(one-xcheb(n_r)**2)/    &
            &             (2*tmp_n_r_max-2)
            dddf(n_f,n_r)=-tmp_in(n_f,n_r)*(one+two*xcheb(n_r)**2)/                 &
            &             (one-xcheb(n_r)**2)**2.5_cp/(2*tmp_n_r_max-2)+            &
            &             work_2(n_f,n_r)*three*xcheb(n_r)/(one-xcheb(n_r)**2)**2/  &
            &             (2*tmp_n_r_max-2)-work_3(n_f,n_r)/                        &
            &             (one-xcheb(n_r)**2)**1.5_cp/(2*tmp_n_r_max-2)
         end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine get_dddr_fft
!------------------------------------------------------------------------------
#endif
end module cosine_transform_gpu
