#include "perflib_preproc.cpp"
#include "mkl_dfti.f90"

module fft

   use precision_mod
   use constants, only: one, zero
   use truncation, only: nrp, ncp, n_phi_max, n_theta_loc, n_m_max, n_r_loc
   use blocking, only: nfs
   !use parallel_mod, only: nThreads
   use mkl_dfti
 
   implicit none
   
   private
 
   !----------- MKL specific variables -------------
   integer :: status
   type(DFTI_DESCRIPTOR), pointer :: c2r_handle, r2c_handle
   type(DFTI_DESCRIPTOR), pointer :: phi2m_dhandle, m2phi_dhandle
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_1, m2p_handle_1
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_2, m2p_handle_2
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_3, m2p_handle_3
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_4, m2p_handle_4
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_5, m2p_handle_5
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_6, m2p_handle_6
   type(DFTI_DESCRIPTOR), pointer :: p2m_handle_7, m2p_handle_7
   !----------- END MKL specific variables
 
   public :: fft_thetab, init_fft, fft_to_real, finalize_fft
   public :: finalize_fft_phi, fft_phi_loc
   public :: initialize_fft_phi
   public :: finalize_fft_phi_many, fft_phi_many
   public :: initialize_fft_phi_many, ifft_phi, fft_phi

contains

   subroutine init_fft(number_of_points)
      !
      ! MKL equivalent of init_fft
      !

      !-- Input variable
      integer, intent(in) :: number_of_points ! number of points

      ! Fourier transformation complex->REAL with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( c2r_handle, DFTI_DOUBLE, DFTI_REAL, &
                                     1, number_of_points )
      status = DftiSetValue( c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
      status = DftiSetValue( c2r_handle, DFTI_INPUT_DISTANCE, nrp )
      status = DftiSetValue( c2r_handle, DFTI_OUTPUT_DISTANCE, nrp )
      !status = DftiSetValue( c2r_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
      !                       DFTI_COMPLEX_COMPLEX )
#ifdef WITH_SHTNS
      status = DftiSetValue( c2r_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
#else
      status = DftiSetValue( c2r_handle, DFTI_PLACEMENT, DFTI_INPLACE )
#endif
      ! This is supposed to fix the issue for ifort < 14 but it increases the
      ! walltime
      !status = DftiSetValue( c2r_handle, DFTI_NUMBER_OF_USER_THREADS, nThreads)
      status = DftiCommitDescriptor( c2r_handle )
  
      ! Fourier transformation REAL->complex with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( r2c_handle, DFTI_DOUBLE, DFTI_REAL, &
                                     1, number_of_points )
      status = DftiSetValue( r2c_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
      status = DftiSetValue( r2c_handle, DFTI_INPUT_DISTANCE, nrp )
      status = DftiSetValue( r2c_handle, DFTI_OUTPUT_DISTANCE, nrp )
      !status = DftiSetValue( r2c_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
      !                       DFTI_COMPLEX_COMPLEX )
#ifdef WITH_SHTNS
      status = DftiSetValue( r2c_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
#else
      status = DftiSetValue( r2c_handle, DFTI_PLACEMENT, DFTI_INPLACE )
#endif
      status = DftiSetValue( r2c_handle, DFTI_FORWARD_SCALE, &
                             one/real(number_of_points,cp) )
      !status = DftiSetValue( r2c_handle, DFTI_NUMBER_OF_USER_THREADS, nThreads)
      status = DftiCommitDescriptor( r2c_handle )

   end subroutine init_fft
!------------------------------------------------------------------------------
   subroutine finalize_fft

      status = DftiFreeDescriptor( c2r_handle )
      status = DftiFreeDescriptor( r2c_handle )

   end subroutine finalize_fft
!------------------------------------------------------------------------------
   subroutine fft_thetab(f,dir)

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
      end if
      PERFOFF

   end subroutine fft_thetab
!------------------------------------------------------------------------------
   subroutine fft_to_real(f,ld_f,nrep)

      integer,  intent(in) :: ld_f,nrep
      real(cp), intent(inout) :: f(ld_f,nrep)

      type(DFTI_DESCRIPTOR), pointer :: local_c2r_handle
#ifdef WITH_SHTNS
      real(cp), allocatable :: work_array(:, :)
      allocate(work_array(ld_f+2, nrep))
#endif

      PERFON('fft2r')
      ! Fourier transformation complex->REAL with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( local_c2r_handle, DFTI_DOUBLE, &
                                     DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( local_c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nrep )
      status = DftiSetValue( local_c2r_handle, DFTI_OUTPUT_DISTANCE, ld_f )
#ifdef WITH_SHTNS
      status = DftiSetValue( local_c2r_handle, DFTI_INPUT_DISTANCE, ld_f+2 )
      status = DftiSetValue( local_c2r_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
#else
      status = DftiSetValue( local_c2r_handle, DFTI_INPUT_DISTANCE, ld_f )
      status = DftiSetValue( local_c2r_handle, DFTI_PLACEMENT, DFTI_INPLACE )
#endif
      status = DftiCommitDescriptor( local_c2r_handle )

      ! run FFT
#ifdef WITH_SHTNS
      work_array(1:n_phi_max, :) = f(1:n_phi_max, :)
      work_array(n_phi_max+1:n_phi_max+2, :) = 0.0_cp
      status = DftiComputeBackward( local_c2r_handle, work_array(:, 1), f(:,1) )
      deallocate(work_array)
#else
      status = DftiComputeBackward( local_c2r_handle, f(:,1) )
#endif

      status = DftiFreeDescriptor( local_c2r_handle )
      PERFOFF

   end subroutine fft_to_real
!------------------------------------------------------------------------------


   subroutine initialize_fft_phi
      !
      !   The phi2m_dhandle is the handle for the real transform.
      !   It will store only Nφ/2+1 modes, because the others are the conjugate 
      !   of these. This is for the distributed case as well (and thus, will 
      !   repeat n_theta_loc times)
      !   
      !   Author: Rafael Lago, MPCDF, July 2017
      !   
      integer :: status
      
      status = DftiCreateDescriptor( phi2m_dhandle, DFTI_DOUBLE, DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( phi2m_dhandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      status = DftiSetValue( phi2m_dhandle, DFTI_NUMBER_OF_TRANSFORMS, n_theta_loc)
      status = DftiSetValue( phi2m_dhandle, DFTI_INPUT_DISTANCE, n_phi_max )
      status = DftiSetValue( phi2m_dhandle, DFTI_OUTPUT_DISTANCE, n_phi_max/2+1 )
      status = DftiSetValue( phi2m_dhandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiSetValue( phi2m_dhandle, DFTI_FORWARD_SCALE, 1.0_cp/real(n_phi_max,cp) )
      status = DftiCommitDescriptor( phi2m_dhandle )
      
      status = DftiCreateDescriptor( m2phi_dhandle, DFTI_DOUBLE, DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( m2phi_dhandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      status = DftiSetValue( m2phi_dhandle, DFTI_NUMBER_OF_TRANSFORMS, n_theta_loc)
      status = DftiSetValue( m2phi_dhandle, DFTI_INPUT_DISTANCE, n_phi_max/2+1 )
      status = DftiSetValue( m2phi_dhandle, DFTI_OUTPUT_DISTANCE, n_phi_max )
      status = DftiSetValue( m2phi_dhandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiSetValue( m2phi_dhandle, DFTI_FORWARD_SCALE, 1.0_cp/real(n_phi_max,cp) )
      status = DftiCommitDescriptor( m2phi_dhandle )
    
   end subroutine initialize_fft_phi
   !----------------------------------------------------------------------------
   subroutine finalize_fft_phi
      !
      !   Author: Rafael Lago, MPCDF, July 2017
      !   
   
      integer :: status
      
      status = DftiFreeDescriptor(phi2m_dhandle)
      status = DftiFreeDescriptor(m2phi_dhandle)
   end subroutine
  
   !----------------------------------------------------------------------------
   subroutine fft_phi_loc(f, g, dir)
      !   FFT meant to be used for θ-distributed vectors. It does not perform 
      !   the transposition. That can be achieved with a call to
      !   transpose_m_theta or transpose_theta_m.
      !   
      !   Forward:  converts f(φ,θ_loc) to g(m,θ_loc) 
      !   Backward: converts g(m,θ_loc) to f(φ,θ_loc)
      !   
      !   OBS: g has (n_phi_max/2+1)*n_theta_loc points, not m_max*n_theta_loc points.
      !       Only the first m_max rows are required by MagIC.
      !  
      !   Author: Rafael Lago, MPCDF, November 2017
      ! 
      real(cp),    intent(inout)  :: f(n_phi_max*n_theta_loc)
      complex(cp), intent(inout)  :: g((n_phi_max/2+1)*n_theta_loc)
      integer,     intent(in)     :: dir
      integer :: status
      
      PERFON('fftphi')      
      if (dir == 1) then
         status = DftiComputeForward( phi2m_dhandle, f(:), g(:) )
      else if (dir == -1) then
         status = DftiComputeBackward( m2phi_dhandle, g(:), f(:) )
      else
         print *, "Unknown direction in fft_phi_loc: ", dir
      end if
      PERFOFF
    
   end subroutine fft_phi_loc
!-----------------------------------------------------------------------------------
   subroutine initialize_fft_phi_many

      call initialize_plan(p2m_handle_1, m2p_handle_1, 1)
      call initialize_plan(p2m_handle_2, m2p_handle_2, 2)
      call initialize_plan(p2m_handle_3, m2p_handle_3, 3)
      call initialize_plan(p2m_handle_4, m2p_handle_4, 4)
      call initialize_plan(p2m_handle_5, m2p_handle_5, 5)
      call initialize_plan(p2m_handle_6, m2p_handle_6, 6)
      call initialize_plan(p2m_handle_7, m2p_handle_7, 7)

   end subroutine initialize_fft_phi_many
!-----------------------------------------------------------------------------------
   subroutine finalize_fft_phi_many

      call finalize_plan(p2m_handle_1, m2p_handle_1)
      call finalize_plan(p2m_handle_2, m2p_handle_2)
      call finalize_plan(p2m_handle_3, m2p_handle_3)
      call finalize_plan(p2m_handle_4, m2p_handle_4)
      call finalize_plan(p2m_handle_5, m2p_handle_5)
      call finalize_plan(p2m_handle_6, m2p_handle_6)
      call finalize_plan(p2m_handle_7, m2p_handle_7)

   end subroutine finalize_fft_phi_many
!-----------------------------------------------------------------------------------
   subroutine initialize_plan(p2m_handle, m2p_handle, n_f)

      integer, intent(in) :: n_f
      type(DFTI_DESCRIPTOR), pointer, intent(inout) ::  p2m_handle
      type(DFTI_DESCRIPTOR), pointer, intent(inout) ::  m2p_handle

      integer :: st

      st = DftiCreateDescriptor( p2m_handle, DFTI_DOUBLE, DFTI_REAL, &
           &                     1, n_phi_max )
      st = DftiSetValue( p2m_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
           &             DFTI_COMPLEX_COMPLEX)
      st = DftiSetValue( p2m_handle, DFTI_NUMBER_OF_TRANSFORMS, &
           &             n_theta_loc*n_r_loc*n_f )
      st = DftiSetValue( p2m_handle, DFTI_INPUT_DISTANCE, n_phi_max )
      st = DftiSetValue( p2m_handle, DFTI_OUTPUT_DISTANCE, &
           &             n_phi_max/2+1 )
      st = DftiSetValue( p2m_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      st = DftiSetValue( p2m_handle, DFTI_FORWARD_SCALE, &
           &             1.0_cp/real(n_phi_max,cp) )
      st = DftiCommitDescriptor( p2m_handle )
         
      st = DftiCreateDescriptor( m2p_handle, DFTI_DOUBLE, DFTI_REAL, 1, &
           &                     n_phi_max )
      st = DftiSetValue( m2p_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
           &             DFTI_COMPLEX_COMPLEX)
      st = DftiSetValue( m2p_handle, DFTI_NUMBER_OF_TRANSFORMS, &
           &             n_theta_loc*n_r_loc*n_f )
      st = DftiSetValue( m2p_handle, DFTI_INPUT_DISTANCE, &
           &             n_phi_max/2+1 )
      st = DftiSetValue( m2p_handle, DFTI_OUTPUT_DISTANCE, n_phi_max )
      st = DftiSetValue( m2p_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      st = DftiSetValue( m2p_handle, DFTI_FORWARD_SCALE, &
           &             1.0_cp/real(n_phi_max,cp) )
      st = DftiCommitDescriptor( m2p_handle )
    
   end subroutine initialize_plan
!-----------------------------------------------------------------------------------
   subroutine finalize_plan(p2m_handle, m2p_handle)
   
      type(DFTI_DESCRIPTOR), pointer, intent(inout) :: p2m_handle
      type(DFTI_DESCRIPTOR), pointer, intent(inout) :: m2p_handle

      integer :: st
      
      st = DftiFreeDescriptor(p2m_handle)
      st = DftiFreeDescriptor(m2p_handle)

   end subroutine finalize_plan
!-----------------------------------------------------------------------------------
   subroutine fft_phi_many(f, g, dir)

      real(cp),    intent(inout)  :: f(n_phi_max,n_theta_loc,n_r_loc)
      complex(cp), intent(inout)  :: g(n_m_max,n_theta_loc,n_r_loc)
      integer,     intent(in)     :: dir

      !-- Local variables
      complex(cp) :: tmp((n_phi_max/2+1),n_theta_loc,n_r_loc)
      integer :: st
      
      if (dir == 1) then
         st = DftiComputeForward( p2m_handle_1, f(:,1,1), tmp(:,1,1) )
         g(:,:,:)=tmp(1:n_m_max,:,:)
      else if (dir == -1) then
         tmp(1:n_m_max,:,:) =g(:,:,:)
         tmp(n_m_max+1:,:,:)=zero
         st = DftiComputeBackward( m2p_handle_1, tmp(:,1,1), f(:,1,1) )
      else
         print *, "Unknown direction in fft_phi_loc: ", dir
      end if
    
  end subroutine fft_phi_many
!-----------------------------------------------------------------------------------
   subroutine ifft_phi(g,f,n_fields)

      !-- Input variables
      integer,     intent(in) :: n_fields
      complex(cp), intent(in)  :: g(n_phi_max/2+1,n_theta_loc,n_r_loc,n_fields)

      !-- Output variables
      real(cp),    intent(out)  :: f(n_phi_max,n_theta_loc,n_r_loc,n_fields)

      !-- Local variables
      integer :: st

      if ( n_fields == 1 ) then
         st = DftiComputeBackward( m2p_handle_1, g(:,1,1,1), f(:,1,1,1) )
      else if ( n_fields == 2 ) then
         st = DftiComputeBackward( m2p_handle_2, g(:,1,1,1), f(:,1,1,1) )
      else if ( n_fields == 3 ) then
         st = DftiComputeBackward( m2p_handle_3, g(:,1,1,1), f(:,1,1,1) )
      else if ( n_fields == 4 ) then
         st = DftiComputeBackward( m2p_handle_4, g(:,1,1,1), f(:,1,1,1) )
      else if ( n_fields == 5 ) then
         st = DftiComputeBackward( m2p_handle_5, g(:,1,1,1), f(:,1,1,1) )
      else if ( n_fields == 6 ) then
         st = DftiComputeBackward( m2p_handle_6, g(:,1,1,1), f(:,1,1,1) )
      else if ( n_fields == 7 ) then
         st = DftiComputeBackward( m2p_handle_7, g(:,1,1,1), f(:,1,1,1) )
      end if

   end subroutine ifft_phi
!-----------------------------------------------------------------------------------
   subroutine fft_phi(f,g,n_fields)

      !-- Input variables
      integer,     intent(in) :: n_fields
      real(cp),    intent(in)  :: f(n_phi_max,n_theta_loc,n_r_loc,n_fields)

      !-- Output variables
      complex(cp), intent(out)  :: g(n_phi_max/2+1,n_theta_loc,n_r_loc,n_fields)

      !-- Local variables
      integer :: st

      if ( n_fields == 1 ) then
         st = DftiComputeForward( p2m_handle_1, f(:,1,1,1), g(:,1,1,1) )
      else if ( n_fields == 2 ) then
         st = DftiComputeForward( p2m_handle_2, f(:,1,1,1), g(:,1,1,1) )
      else if ( n_fields == 3 ) then
         st = DftiComputeForward( p2m_handle_3, f(:,1,1,1), g(:,1,1,1) )
      else if ( n_fields == 4 ) then
         st = DftiComputeForward( p2m_handle_4, f(:,1,1,1), g(:,1,1,1) )
      else if ( n_fields == 5 ) then
         st = DftiComputeForward( p2m_handle_5, f(:,1,1,1), g(:,1,1,1) )
      else if ( n_fields == 6 ) then
         st = DftiComputeForward( p2m_handle_6, f(:,1,1,1), g(:,1,1,1) )
      else if ( n_fields == 7 ) then
         st = DftiComputeForward( p2m_handle_7, f(:,1,1,1), g(:,1,1,1) )
      end if

   end subroutine fft_phi
!-----------------------------------------------------------------------------------
end module fft
