module hipfort_check

contains

#ifdef USE_CUDA_NAMES
   subroutine hipCheck(cudaError_t)
      use hipfort_cuda_errors
      use hipfort_enums
      implicit none
      integer(kind(cudaSuccess)) :: cudaError_t
      integer(kind(hipSuccess)) :: hipError_t

      hipError_t = hipCUDAErrorTohipError(cudaError_t)

      if (hipError_t /= hipSuccess) then
         write (*, *) "HIP ERROR: Error code = ", hipError_t, ", CUDA error code = ", cudaError_t
         call exit(hipError_t)
      end if
   end subroutine hipCheck
#else
   subroutine hipCheck(hipError_t)
      use hipfort_enums
      implicit none

      integer(kind(hipSuccess)) :: hipError_t

      if (hipError_t /= hipSuccess) then
         write (*, *) "HIP ERROR: Error code = ", hipError_t
         call exit(hipError_t)
      end if
   end subroutine hipCheck
#endif

  ! HIP math libs
  ! TODO: Currently, only AMDGPU is supported
  
  subroutine hipblasCheck(hipblasError_t)
    use hipfort_hipblas_enums

    implicit none

    integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasError_t

    if(hipblasError_t /= HIPBLAS_STATUS_SUCCESS)then
       write(*,*) "HIPBLAS ERROR: Error code = ", hipblasError_t
       call exit(hipblasError_t)
    end if
  end subroutine hipblasCheck

  subroutine hipfftCheck(hipfft_status)
    use hipfort_hipfft_enums

    implicit none

    integer(kind(hipfft_success)) :: hipfft_status

    if(hipfft_status /= hipfft_success)then
       write(*,*) "HIPFFT ERROR: Error code = ", hipfft_status
       call exit(hipfft_status)
    end if
  end subroutine hipfftCheck

  subroutine hipsolverCheck(hipsolverError_t)
    use hipfort_hipsolver_enums

    implicit none

    integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverError_t

    if(hipsolverError_t /= HIPSOLVER_STATUS_SUCCESS)then
       write(*,*) "HIPSOLVER ERROR: Error code = ", hipsolverError_t
       call exit(hipsolverError_t)
    end if
  end subroutine hipsolverCheck
  
end module hipfort_check
