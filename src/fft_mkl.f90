#include "perflib_preproc.cpp"
#include "mkl_dfti.f90"

module fft

   use precision_mod
   use truncation, only: nrp, ncp, n_phi_max
   use blocking, only: nfs
   use mkl_dfti
 
   implicit none
   
   private
 
   !----------- MKL specific variables -------------
   integer :: status
   type(DFTI_DESCRIPTOR), pointer :: c2r_handle, r2c_handle
   !----------- END MKL specific variables
 
   interface fft_thetab
      module procedure fft_thetab_real
      module procedure fft_thetab_cmplx
   end interface
 
   public :: fft_thetab, init_fft, fft_to_real

contains

   subroutine init_fft(number_of_points)
      ! 
      !  MKL equivalent of init_fft
      !

      integer, intent(in) :: number_of_points ! number of points
      
      ! Fourier transformation complex->REAL with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( c2r_handle, DFTI_DOUBLE, DFTI_REAL, &
                                     1, number_of_points )
      status = DftiSetValue( c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
      status = DftiSetValue( c2r_handle, DFTI_INPUT_DISTANCE, ncp )
      status = DftiSetValue( c2r_handle, DFTI_OUTPUT_DISTANCE, nrp )
      status = DftiSetValue( c2r_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
                             DFTI_COMPLEX_COMPLEX )
      status = DftiSetValue( c2r_handle, DFTI_PLACEMENT, DFTI_INPLACE )
      status = DftiCommitDescriptor( c2r_handle )
  
      ! Fourier transformation REAL->complex with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( r2c_handle, DFTI_DOUBLE, DFTI_REAL, &
                                     1, number_of_points )
      status = DftiSetValue( r2c_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
      status = DftiSetValue( r2c_handle, DFTI_INPUT_DISTANCE, nrp )
      status = DftiSetValue( r2c_handle, DFTI_OUTPUT_DISTANCE, ncp )
      status = DftiSetValue( r2c_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
                             DFTI_COMPLEX_COMPLEX )
      status = DftiSetValue( r2c_handle, DFTI_PLACEMENT, DFTI_INPLACE )
      status = DftiSetValue( r2c_handle, DFTI_FORWARD_SCALE, &
                             1.0/real(number_of_points,cp) )
      status = DftiCommitDescriptor( r2c_handle )

   end subroutine init_fft
!------------------------------------------------------------------------------
   subroutine fft_thetab_cmplx(f,dir)

      complex(cp), intent(inout) :: f(nrp/2,nfs)
      integer,     intent(in) :: dir            ! back or forth transform

      PERFON('fft_thc')
      if (dir == 1) then
         ! run FFT
         status = DftiComputeBackward( c2r_handle, f(:,1) )
      else if (dir == -1) then
         ! run FFT
         status = DftiComputeForward( r2c_handle, f(:,1) )
         !PRINT*,"Calling fft_thetab with complex array and dir /= 1. &
         !       Don't know what to do!"
         !call TRACEBACKQQ
         !stop
      end if
      PERFOFF

   end subroutine fft_thetab_cmplx
!------------------------------------------------------------------------------
   subroutine fft_thetab_real(f,dir)

      real(cp), intent(inout) :: f(nrp,nfs)
      integer,  intent(in) :: dir            ! back or forth transform

      PERFON('fft_thr')
      if (dir == -1) then
         ! run FFT
         status = DftiComputeForward( r2c_handle, f(:,1) )
      else if (dir == 1) then
         ! we want a backward transform, and the real array f 
         ! is to be interpreted as complex array

         ! run FFT
         status = DftiComputeBackward( c2r_handle, f(:,1) )

         !PRINT*,"Calling fft_thetab with real array and dir /= -1. & 
         !        Don't know what to do!"
         !call TRACEBACKQQ
         !stop
      end if
      PERFOFF

   end subroutine fft_thetab_real
!------------------------------------------------------------------------------
   subroutine fft_to_real(f,ld_f,nrep)

      integer,  intent(in) :: ld_f,nrep
      real(cp), intent(inout) :: f(ld_f,nrep)

      type(DFTI_DESCRIPTOR), pointer :: local_c2r_handle

      PERFON('fft2r')
      ! Fourier transformation complex->REAL with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( local_c2r_handle, DFTI_DOUBLE, &
                                     DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( local_c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nrep )
      status = DftiSetValue( local_c2r_handle, DFTI_INPUT_DISTANCE, ld_f )
      status = DftiSetValue( local_c2r_handle, DFTI_OUTPUT_DISTANCE, ld_f )
      status = DftiSetValue( local_c2r_handle, DFTI_PLACEMENT, DFTI_INPLACE )
      status = DftiCommitDescriptor( local_c2r_handle )

      ! run FFT
      status = DftiComputeBackward( local_c2r_handle, f(:,1) )

      status = DftiFreeDescriptor( local_c2r_handle )
      PERFOFF

   end subroutine fft_to_real
!------------------------------------------------------------------------------
end module fft
