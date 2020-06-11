module fft
   !
   ! This is the FFTW version of the fft module
   !

   use iso_c_binding
   use precision_mod
   use constants, only: zero
   use truncation, only: n_phi_max, n_theta_loc, n_r_loc

   implicit none

   include 'fftw3.f03'

   !type(c_ptr) :: r2c_handle, c2r_handle
   type(c_ptr) :: plan_forward, plan_backward
   type(c_ptr) :: plan_forward_1, plan_backward_1
   type(c_ptr) :: plan_forward_2, plan_backward_2
   type(c_ptr) :: plan_forward_3, plan_backward_3
   type(c_ptr) :: plan_forward_4, plan_backward_4
   type(c_ptr) :: plan_forward_5, plan_backward_5
   type(c_ptr) :: plan_forward_6, plan_backward_6
   type(c_ptr) :: plan_forward_7, plan_backward_7

   integer(c_int), parameter :: fftw_plan_flag=FFTW_PATIENT

   public :: initialize_fft_phi, finalize_fft_phi, fft_phi_loc, fft_to_real,    &
   &         initialize_fft_phi_many, finalize_fft_phi_many, fft_phi, ifft_phi, &
   &         init_fft, finalize_fft, fft_thetab

contains 

   subroutine init_fft(n_points)
      !@> TODO: not certain this is needed since when n_ranks_theta=1
      ! the fft_phi_loc should do the job. Keep it for now

      integer, intent(in) :: n_points

      !call initialize_plan(r2c_handle, c2r_handle, n_points, n_theta_loc)

   end subroutine init_fft
!-----------------------------------------------------------------------------------
   subroutine finalize_fft
      !@> TODO: not certain this is needed since when n_ranks_theta=1
      ! the fft_phi_loc should do the job. Keep it for now

      !call finalize_plan(r2c_handle, c2r_handle)

   end subroutine finalize_fft
!-----------------------------------------------------------------------------------
   subroutine fft_thetab(f,dir)

      use truncation, only: nrp, n_theta_max
      use useful, only: abortRun

      real(cp), intent(inout) :: f(nrp,n_theta_max)
      integer,  intent(in) :: dir       ! back or forth transform

      call abortRun('Not ported: likely not needed!')

   end subroutine fft_thetab
!-----------------------------------------------------------------------------------
   subroutine fft_to_real(f,ld_f,nrep)
      !@ TODO> Right now this works only for SHTns config like (i.e. out of
      !place with n_phi_max points and not n_phi_max+2)

      !-- Input variables
      integer,  intent(in) :: ld_f,nrep

      !-- Inout variable
      real(cp), intent(inout) :: f(ld_f,nrep)

      !-- Local variables
      type(c_ptr) :: plan_bwd
      complex(cp) :: fcplx(ld_f/2+1,nrep)
      integer :: inembed(1), onembed(1), istride, ostride, odist, idist, n(1), np, nr

      inembed(1) = 0
      onembed(1) = 0
      istride = 1
      ostride = 1
      odist = ld_f ! or 2*(ld_f/2+1) if in-place
      idist = ld_f/2+1
      n(1)=ld_f

      plan_bwd = fftw_plan_many_dft_c2r(1, n, nrep, fcplx, inembed, istride, idist, &
                 &                      f, onembed, ostride, odist, FFTW_ESTIMATE)

      do nr=1,nrep
         do np=1,ld_f/2
            fcplx(np,nr)=cmplx(f(2*np-1,nr),f(2*np,nr),cp)
         end do
         fcplx(ld_f/2+1,nr)=zero
      end do

      call fftw_execute_dft_c2r(plan_bwd, fcplx, f)
      call fftw_destroy_plan(plan_bwd)

   end subroutine fft_to_real
!-----------------------------------------------------------------------------------
   subroutine initialize_fft_phi

      call initialize_plan(plan_forward, plan_backward, n_phi_max, n_theta_loc)

   end subroutine initialize_fft_phi
!----------------------------------------------------------------------------------
   subroutine finalize_fft_phi

      call fftw_destroy_plan(plan_forward)
      call fftw_destroy_plan(plan_backward)

   end subroutine finalize_fft_phi
!-----------------------------------------------------------------------------------
   subroutine initialize_plan(plan_fwd, plan_bwd, n_points, n_fields)

      !-- Input variables
      integer, intent(in) :: n_points
      integer, intent(in) :: n_fields

      !-- Ouput variables
      type(c_ptr), intent(out) :: plan_bwd
      type(c_ptr), intent(out) :: plan_fwd

      !-- Local variables
      complex(cp) :: dat_in(n_points/2+1, n_theta_loc, n_fields)
      real(cp) :: dat_out(n_points, n_theta_loc, n_fields)
      integer :: inembed(1), onembed(1), istride, ostride, odist, idist, n(1)

      inembed(1) = 0
      onembed(1) = 0
      istride = 1
      ostride = 1
      odist = n_points ! or 2*(n_points/2+1) if in-place
      idist = n_points/2+1
      n(1)=n_points

      plan_bwd = fftw_plan_many_dft_c2r(1, n, n_fields, dat_in,             &
                 &                      inembed, istride, idist, dat_out,   &
                 &                      onembed, ostride, odist, fftw_plan_flag)

      idist = n_points
      odist = n_points/2+1
      plan_fwd = fftw_plan_many_dft_r2c(1, n, n_fields, dat_out,              &
                 &                      inembed, istride, idist, dat_in,      &
                 &                      onembed, ostride, odist, fftw_plan_flag)

   end subroutine initialize_plan
!-----------------------------------------------------------------------------------
   subroutine finalize_plan(plan_bwd, plan_fwd)

      type(c_ptr), intent(inout) :: plan_bwd
      type(c_ptr), intent(inout) :: plan_fwd

      call fftw_destroy_plan(plan_fwd)
      call fftw_destroy_plan(plan_bwd)

   end subroutine finalize_plan
!-----------------------------------------------------------------------------------
   subroutine fft_phi_loc(f, g, dir)

      real(cp),    intent(inout)  :: f(n_phi_max,n_theta_loc)
      complex(cp), intent(inout)  :: g((n_phi_max/2+1),n_theta_loc)
      integer,     intent(in)     :: dir
      
      if (dir == 1) then
         call fftw_execute_dft_r2c(plan_forward, f, g)
         g(:,:)=g(:,:)/n_phi_max
      else if (dir == -1) then
         call fftw_execute_dft_c2r(plan_backward, g, f)
      end if
    
   end subroutine fft_phi_loc
!-----------------------------------------------------------------------------------
   subroutine initialize_fft_phi_many

      call initialize_plan(plan_forward_1, plan_backward_1, n_phi_max, n_theta_loc*n_r_loc)
      call initialize_plan(plan_forward_2, plan_backward_2, n_phi_max, 2*n_theta_loc*n_r_loc)
      call initialize_plan(plan_forward_3, plan_backward_3, n_phi_max, 3*n_theta_loc*n_r_loc)
      call initialize_plan(plan_forward_4, plan_backward_4, n_phi_max, 4*n_theta_loc*n_r_loc)
      call initialize_plan(plan_forward_5, plan_backward_5, n_phi_max, 5*n_theta_loc*n_r_loc)
      call initialize_plan(plan_forward_6, plan_backward_6, n_phi_max, 6*n_theta_loc*n_r_loc)
      call initialize_plan(plan_forward_7, plan_backward_7, n_phi_max, 7*n_theta_loc*n_r_loc)

   end subroutine initialize_fft_phi_many
!-----------------------------------------------------------------------------------
   subroutine finalize_fft_phi_many

      call finalize_plan(plan_forward_1, plan_backward_1)
      call finalize_plan(plan_forward_2, plan_backward_2)
      call finalize_plan(plan_forward_3, plan_backward_3)
      call finalize_plan(plan_forward_4, plan_backward_4)
      call finalize_plan(plan_forward_5, plan_backward_5)
      call finalize_plan(plan_forward_6, plan_backward_6)
      call finalize_plan(plan_forward_7, plan_backward_7)

   end subroutine finalize_fft_phi_many
!-----------------------------------------------------------------------------------
   subroutine ifft_phi(g,f,n_fields)

      !-- Input variables
      integer,     intent(in) :: n_fields
      complex(cp), intent(inout)  :: g(n_phi_max/2+1,n_theta_loc,n_r_loc,n_fields)

      !-- Output variables
      real(cp),    intent(out)  :: f(n_phi_max,n_theta_loc,n_r_loc,n_fields)

      if ( n_fields == 1 ) then
         call fftw_execute_dft_c2r(plan_backward_1, g, f)
      else if ( n_fields == 2 ) then
         call fftw_execute_dft_c2r(plan_backward_2, g, f)
      else if ( n_fields == 3 ) then
         call fftw_execute_dft_c2r(plan_backward_3, g, f)
      else if ( n_fields == 4 ) then
         call fftw_execute_dft_c2r(plan_backward_4, g, f)
      else if ( n_fields == 5 ) then
         call fftw_execute_dft_c2r(plan_backward_5, g, f)
      else if ( n_fields == 6 ) then
         call fftw_execute_dft_c2r(plan_backward_6, g, f)
      else if ( n_fields == 7 ) then
         call fftw_execute_dft_c2r(plan_backward_7, g, f)
      end if

   end subroutine ifft_phi
!-----------------------------------------------------------------------------------
   subroutine fft_phi(f,g,n_fields)

      !-- Input variables
      integer,     intent(in) :: n_fields
      real(cp),    intent(inout)  :: f(n_phi_max,n_theta_loc,n_r_loc,n_fields)

      !-- Output variables
      complex(cp), intent(out)  :: g(n_phi_max/2+1,n_theta_loc,n_r_loc,n_fields)

      if ( n_fields == 1 ) then
         call fftw_execute_dft_r2c(plan_forward_1, f, g)
      else if ( n_fields == 2 ) then
         call fftw_execute_dft_r2c(plan_forward_2, f, g)
      else if ( n_fields == 3 ) then
         call fftw_execute_dft_r2c(plan_forward_3, f, g)
      else if ( n_fields == 4 ) then
         call fftw_execute_dft_r2c(plan_forward_4, f, g)
      else if ( n_fields == 5 ) then
         call fftw_execute_dft_r2c(plan_forward_5, f, g)
      else if ( n_fields == 6 ) then
         call fftw_execute_dft_r2c(plan_forward_6, f, g)
      else if ( n_fields == 7 ) then
         call fftw_execute_dft_r2c(plan_forward_7, f, g)
      end if
      g=g/n_phi_max

   end subroutine fft_phi
!-----------------------------------------------------------------------------------
end module fft
