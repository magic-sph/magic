!$Id: m_fftJW.F90 394 2013-02-13 09:40:47Z dannert $
!*******************************************************************************
#include "perflib_preproc.cpp"
#include "mkl_dfti.f90"

MODULE fft_mkl
  use truncation
  use blocking
  use mkl_dfti
  IMPLICIT NONE
  
  private

  !----------- MKL specific variables -------------
  INTEGER :: status
  TYPE(DFTI_DESCRIPTOR),POINTER :: c2r_handle, r2c_handle
  !$OMP THREADPRIVATE( c2r_handle, r2c_handle )
  !----------- END MKL specific variables


  INTERFACE fft_thetab
     module procedure fft_thetab_real
     module procedure fft_thetab_cmplx
  END INTERFACE

  PUBLIC :: fft_thetab, init_fft,fft_to_real

CONTAINS

  SUBROUTINE init_fft(number_of_points)
    integer :: number_of_points
    
    ! Fourier transformation COMPLEX->REAL with MKL DFTI interface
    ! init FFT
    status = DftiCreateDescriptor( c2r_handle, DFTI_DOUBLE, DFTI_REAL, 1, number_of_points )
    status = DftiSetValue( c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
    status = DftiSetValue( c2r_handle, DFTI_INPUT_DISTANCE, ncp )
    status = DftiSetValue( c2r_handle, DFTI_OUTPUT_DISTANCE, nrp )
    status = DftiSetValue( c2r_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX )
    status = DftiSetValue( c2r_handle, DFTI_PLACEMENT, DFTI_INPLACE )
    status = DftiCommitDescriptor( c2r_handle )

    ! Fourier transformation REAL->COMPLEX with MKL DFTI interface
    ! init FFT
    status = DftiCreateDescriptor( r2c_handle, DFTI_DOUBLE, DFTI_REAL, 1, number_of_points )
    status = DftiSetValue( r2c_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
    status = DftiSetValue( r2c_handle, DFTI_INPUT_DISTANCE, nrp )
    status = DftiSetValue( r2c_handle, DFTI_OUTPUT_DISTANCE, ncp )
    status = DftiSetValue( r2c_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX )
    status = DftiSetValue( r2c_handle, DFTI_PLACEMENT, DFTI_INPLACE )
    status = DftiSetValue( r2c_handle, DFTI_FORWARD_SCALE, 1.0/DBLE(number_of_points) )
    status = DftiCommitDescriptor( r2c_handle )
  END SUBROUTINE init_fft

  SUBROUTINE fft_thetab_cmplx(f,dir)
    COMPLEX(kind=8),intent(INOUT) :: f(nrp/2,nfs)

    INTEGER,intent(IN) :: dir            ! back or forth transform

    PERFON('fft_thc')
    if (dir.eq.1) then
       ! run FFT
       status = DftiComputeBackward( c2r_handle, f(:,1) )
    ELSE IF (dir.EQ.-1) THEN
       ! run FFT
       status = DftiComputeForward( r2c_handle, f(:,1) )
       !PRINT*,"Calling fft_thetab with complex array and dir.ne.1. Don't know what to do!"
       !CALL TRACEBACKQQ
       !stop
    END IF
    PERFOFF
  END SUBROUTINE fft_thetab_cmplx

  SUBROUTINE fft_thetab_real(f,dir)
    REAL(kind=8),intent(INOUT) :: f(nrp,nfs)
    INTEGER,intent(IN) :: dir            ! back or forth transform

    PERFON('fft_thr')
    IF (dir.EQ.-1) then
       ! run FFT
       status = DftiComputeForward( r2c_handle, f(:,1) )
    ELSE if (dir.eq.1) THEN
       ! we want a backward transform, and the real array f is to be interpreted as complex array

       ! run FFT
       status = DftiComputeBackward( c2r_handle, f(:,1) )

       !PRINT*,"Calling fft_thetab with real array and dir.ne.-1. Don't know what to do!"
       !CALL TRACEBACKQQ
       !stop
    END IF
    PERFOFF

  END SUBROUTINE fft_thetab_real

  SUBROUTINE fft_to_real(f,ld_f,nrep)
    INTEGER,INTENT(IN) :: ld_f,nrep
    COMPLEX(kind=8), INTENT(INOUT) :: f(ld_f/2,nrep)

    TYPE(DFTI_DESCRIPTOR),POINTER :: local_c2r_handle

    PERFON('fft2r')
    ! Fourier transformation COMPLEX->REAL with MKL DFTI interface
    ! init FFT
    status = DftiCreateDescriptor( local_c2r_handle, DFTI_DOUBLE, DFTI_REAL, 1, n_phi_max )
    status = DftiSetValue( local_c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nrep )
    status = DftiSetValue( local_c2r_handle, DFTI_INPUT_DISTANCE, ld_f/2 )
    status = DftiSetValue( local_c2r_handle, DFTI_OUTPUT_DISTANCE, ld_f )
    status = DftiSetValue( local_c2r_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX )
    status = DftiSetValue( local_c2r_handle, DFTI_PLACEMENT, DFTI_INPLACE )
    status = DftiCommitDescriptor( local_c2r_handle )

    ! run FFT
    status = DftiComputeBackward( local_c2r_handle, f(:,1) )

    status = DftiFreeDescriptor( local_c2r_handle )
    PERFOFF
  END SUBROUTINE fft_to_real

END MODULE fft_mkl
