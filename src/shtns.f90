!#define EXPLICITE_SYNCHRO
module sht
   !
   ! This module contains is a wrapper of the SHTns routines used in MagIC
   !

   use iso_c_binding
   use iso_fortran_env, only: output_unit
   use precision_mod, only: cp, lip
   use blocking, only: st_map
   use constants, only: one, zero
   use truncation, only: m_max, l_max, n_theta_max, n_phi_max, &
       &                 minc, lm_max, nlat_padded
   use horizontal_data, only: dLh
   use radial_data, only: nRstart, nRstop
   use parallel_mod
#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
   use hipfort_check, only: hipCheck
   use hipfort, only: hipDeviceSynchronize
#endif
#endif

   implicit none

   include "shtns.f03"
#ifdef WITH_OMP_GPU
   include "shtns_cuda.f03"
#endif

   private

   public :: initialize_sht, scal_to_spat, scal_to_grad_spat, torpol_to_spat, &
   &         scal_to_SH, spat_to_sphertor, torpol_to_spat_IC, axi_to_spat,    &
   &         torpol_to_curl_spat_IC, spat_to_qst, sphtor_to_spat,             &
   &         toraxi_to_spat, finalize_sht, torpol_to_spat_single

#ifdef WITH_OMP_GPU
   type(c_ptr), public :: sht_l, sht_l_single, sht_l_gpu, sht_l_single_gpu
#else
   type(c_ptr), public :: sht_l, sht_l_single
#endif

   logical :: l_batched_sh

contains

   subroutine initialize_sht(l_scrambled_theta, l_batched)

      !-- Input variable
      logical, intent(in) :: l_batched

      !-- Output variable
      logical, intent(out) :: l_scrambled_theta

      !-- Local variables
      integer :: norm, layout
#ifdef WITH_OMP_GPU
      integer :: layout_gpu
#endif
      real(cp) :: eps_polar
      type(shtns_info), pointer :: sht_info
      complex(cp) :: tmp(l_max+1)
      real(cp), allocatable :: tmpr(:)
      integer :: howmany
      integer(lip) :: dist

      l_batched_sh = l_batched

      if ( l_batched_sh ) then
         howmany = nRstop-nRstart+1
      else
         howmany = 1
      end if
      dist = int(lm_max, kind=8)

      if ( rank == 0 ) then
         write(output_unit,*) ''
         call shtns_verbose(1)
      end if

      nthreads =  shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
#ifdef SHT_PADDING
#ifdef WITH_OMP_GPU
      if ( .not. l_batched_sh ) then
         layout     = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
         layout_gpu = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING + SHT_ALLOW_GPU
      else
         layout     = SHT_GAUSS + SHT_THETA_CONTIGUOUS
         layout_gpu = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_GPU
      end if
#else
      if ( .not. l_batched_sh ) then
         !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
         layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
      else
         !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS
         layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS
      end if
#endif
#else
#ifdef WITH_OMP_GPU
      layout     = SHT_GAUSS + SHT_THETA_CONTIGUOUS
      layout_gpu = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_GPU
#else
      !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS
      layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS
#endif
#endif

      eps_polar = 1.e-10_cp

      sht_l = shtns_create(l_max, m_max/minc, minc, norm)
      if ( l_batched_sh ) call shtns_set_batch(sht_l, howmany, dist)
      call shtns_robert_form(sht_l, 1) ! Use Robert's form
      call shtns_set_grid(sht_l, layout, eps_polar, n_theta_max, n_phi_max)
#ifdef WITH_OMP_GPU
      sht_l_gpu = shtns_create(l_max, m_max/minc, minc, norm)
      if ( l_batched_sh ) call shtns_set_batch(sht_l_gpu, howmany, dist)
      call shtns_robert_form(sht_l_gpu, 1) ! Use Robert's form
      call shtns_set_grid(sht_l_gpu, layout_gpu, eps_polar, n_theta_max, n_phi_max)
#endif

      call c_f_pointer(cptr=sht_l, fptr=sht_info)
#ifdef SHT_PADDING
      if ( .not. l_batched_sh ) then
         nlat_padded = sht_info%nlat_padded
         if ( nlat_padded /= n_theta_max .and. rank == 0 ) then
            write(output_unit,*) '! SHTns uses theta padding with nlat_padded=', &
            &                    nlat_padded
         end if
      else
         nlat_padded = n_theta_max
      end if
#else
      nlat_padded = n_theta_max
#endif

      if ( rank == 0 ) then
         call shtns_verbose(0)
         write(output_unit,*) ''
      end if

      if ( l_batched_sh ) then
         sht_l_single = shtns_create(l_max, m_max/minc, minc, norm)
         call shtns_set_grid(sht_l_single, layout, eps_polar, n_theta_max, n_phi_max)
         call shtns_robert_form(sht_l_single, 1) ! Use Robert's form
#ifdef WITH_OMP_GPU
         sht_l_single_gpu = shtns_create(l_max, m_max/minc, minc, norm)
         call shtns_set_grid(sht_l_single_gpu, layout_gpu, eps_polar, n_theta_max, n_phi_max)
         call shtns_robert_form(sht_l_single_gpu, 1) ! Use Robert's form
#endif
      else
         sht_l_single = sht_l
#ifdef WITH_OMP_GPU
         sht_l_single_gpu = sht_l_gpu
#endif
      end if

      !-- Set a l=1, m=0 mode to determine whether scrambling is there or not
      allocate( tmpr(nlat_padded) )
      tmp(:)=zero
      tmp(2)=(1.0_cp, 0.0_cp)
      call axi_to_spat(tmp, tmpr)
      if ( abs(tmpr(2)+tmpr(1)) <= 10.0_cp * epsilon(1.0_cp) ) then
         l_scrambled_theta=.true.
      else
         l_scrambled_theta=.false.
      end if
      deallocate( tmpr )

   end subroutine initialize_sht
!------------------------------------------------------------------------------
   subroutine finalize_sht

      call shtns_unset_grid(sht_l)
      call shtns_destroy(sht_l)
      if ( l_batched_sh ) then
         call shtns_unset_grid(sht_l_single)
         call shtns_destroy(sht_l_single)
      end if
#ifdef WITH_OMP_GPU
      call shtns_unset_grid(sht_l_gpu)
      call shtns_destroy(sht_l_gpu)
      if ( l_batched_sh ) then
         call shtns_unset_grid(sht_l_single_gpu)
         call shtns_destroy(sht_l_single_gpu)
      end if
#endif

   end subroutine finalize_sht
!------------------------------------------------------------------------------
   subroutine scal_to_spat(sh, Slm, fieldc, lcut, use_gpu)
      ! transform a spherical harmonic field into grid space

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      complex(cp), intent(inout) :: Slm(*)
      integer,     intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variable
      real(cp), intent(out) :: fieldc(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(Slm, fieldc)
         call cu_SH_to_spat(sh, Slm, fieldc, lcut)
         !$omp end target data
      else
         call SH_to_spat_l(sh, Slm, fieldc, lcut)
      end if
#else
      call SH_to_spat_l(sh, Slm, fieldc, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif


   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc, lcut, use_gpu)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variables
      complex(cp), intent(inout) :: Slm(*)
      integer,     intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variables
      real(cp), intent(out) :: gradtc(*)
      real(cp), intent(out) :: gradpc(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(Slm, gradtc, gradpc)
         call cu_SHsph_to_spat(sht_l_gpu, Slm, gradtc, gradpc, lcut)
         !$omp end target data
      else
         call SHsph_to_spat_l(sht_l, Slm, gradtc, gradpc, lcut)
      endif
#else
      call SHsph_to_spat_l(sht_l, Slm, gradtc, gradpc, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif

   end subroutine scal_to_grad_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut, use_gpu)

      !-- Input variables
      complex(cp), intent(inout) :: Wlm(*), dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variables
      real(cp), intent(out) :: vrc(*)
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(Wlm, dWlm, Zlm, vrc, vtc, vpc)
         call cu_SHqst_to_spat(sht_l_gpu, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)
         !$omp end target data
      else
         call SHqst_to_spat_l(sht_l, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)
      end if
#else
      call SHqst_to_spat_l(sht_l, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif

   end subroutine torpol_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_single(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(*)
      complex(cp), intent(inout) :: dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(*)
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(l)
      do lm = 1, lm_max
         l = st_map%lm2l(lm)
         if ( l <= lcut ) then
            Qlm(lm) = dLh(lm) * Wlm(lm)
         else
            Qlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call SHqst_to_spat_l(sht_l_single, Qlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

   end subroutine torpol_to_spat_single
!------------------------------------------------------------------------------
   subroutine sphtor_to_spat(sh, dWlm, Zlm, vtc, vpc, lcut, use_gpu)

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      complex(cp), intent(inout) :: dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variables
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(dWlm, Zlm, vtc, vpc)
         call cu_SHsphtor_to_spat(sh, dWlm, Zlm, vtc, vpc, lcut)
         !$omp end target data
      else
         call SHsphtor_to_spat_l(sh, dWlm, Zlm, vtc, vpc, lcut)
      end if
#else
      call SHsphtor_to_spat_l(sh, dWlm, Zlm, vtc, vpc, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif

   end subroutine sphtor_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_IC(r, r_ICB, dBlm, ddBlm, Jlm, dJlm, &
              &                      cbr, cbt, cbp)
      !
      ! This is a QST transform that contains the transform for the
      ! inner core to compute the three components of the curl of B.
      !

      !-- Input variables
      real(cp),    intent(in) :: r, r_ICB
      complex(cp), intent(in) :: dBlm(:), ddBlm(:)
      complex(cp), intent(in) :: Jlm(:), dJlm(:)

      !-- Output variables
      real(cp), intent(out) :: cbr(:,:)
      real(cp), intent(out) :: cbt(:,:)
      real(cp), intent(out) :: cbp(:,:)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Slm(lm_max), Tlm(lm_max)
      real(cp) :: rDep(0:l_max), rDep2(0:l_max)
      real(cp) :: rRatio
      integer :: lm, l, m

      rRatio = r/r_ICB
      rDep(0) = rRatio
      rDep2(0)= one/r_ICB
      do l=1,l_max
         rDep(l) =rDep(l-1) *rRatio
         rDep2(l)=rDep2(l-1)*rRatio
      end do

      lm = 0
      do m=0,m_max,minc
         do l=m,l_max
            lm = lm+1
            Qlm(lm) = rDep(l) * dLh(lm) * Jlm(lm)
            Slm(lm) = rDep2(l) * ((l+1)*Jlm(lm)+r*dJlm(lm))
            Tlm(lm) = -rDep2(l) * ( 2*(l+1)*dBlm(lm)+r*ddBlm(lm) )
         end do
      end do

      call SHqst_to_spat(sht_l_single, Qlm, Slm, Tlm, cbr, cbt, cbp)

   end subroutine torpol_to_curl_spat_IC
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_IC(r, r_ICB, Wlm, dWlm, Zlm, Br, Bt, Bp)
      !
      ! This is a QST transform that contains the transform for the
      ! inner core.
      !

      !-- Input variables
      real(cp),    intent(in) :: r, r_ICB
      complex(cp), intent(in) :: Wlm(:), dWlm(:), Zlm(:)

      !-- Output variables
      real(cp), intent(out) :: Br(:,:)
      real(cp), intent(out) :: Bt(:,:)
      real(cp), intent(out) :: Bp(:,:)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Slm(lm_max), Tlm(lm_max)
      real(cp) :: rDep(0:l_max), rDep2(0:l_max)
      real(cp) :: rRatio
      integer :: lm, l, m

      rRatio = r/r_ICB
      rDep(0) = rRatio
      rDep2(0)= one/r_ICB
      do l=1,l_max
         rDep(l) =rDep(l-1) *rRatio
         rDep2(l)=rDep2(l-1)*rRatio
      end do

      lm = 0
      do m=0,m_max,minc
         do l=m,l_max
            lm = lm+1
            Qlm(lm) = rDep(l) * dLh(lm) * Wlm(lm)
            Slm(lm) = rDep2(l) * ((l+1)*Wlm(lm)+r*dWlm(lm))
            Tlm(lm) = rDep(l) * Zlm(lm)
         end do
      end do

      call SHqst_to_spat(sht_l_single, Qlm, Slm, Tlm, Br, Bt, Bp)

   end subroutine torpol_to_spat_IC
!------------------------------------------------------------------------------
   subroutine scal_to_SH(sh, f, fLM, lcut, use_gpu)

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      real(cp),    intent(inout) :: f(*)
      integer,     intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variable
      complex(cp), intent(out) :: fLM(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(f, fLM)
         call cu_spat_to_SH(sh, f, fLM, lcut)
         !$omp end target data
      else
         call spat_to_SH_l(sh, f, fLM, lcut)
      end if
#else
      call spat_to_SH_l(sh, f, fLM, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif

   end subroutine scal_to_SH
!------------------------------------------------------------------------------
   subroutine spat_to_qst(f, g, h, qLM, sLM, tLM, lcut, use_gpu)

      !-- Input variables
      real(cp), intent(inout) :: f(*)
      real(cp), intent(inout) :: g(*)
      real(cp), intent(inout) :: h(*)
      integer,  intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variables
      complex(cp), intent(out) :: qLM(*)
      complex(cp), intent(out) :: sLM(*)
      complex(cp), intent(out) :: tLM(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(f, g, h, qLM, sLM, tLM)
         call cu_spat_to_SHqst(sht_l_gpu, f, g, h, qLM, sLM, tLM, lcut)
         !$omp end target data
      else
         call spat_to_SHqst_l(sht_l, f, g, h, qLM, sLM, tLM, lcut)
      end if
#else
      call spat_to_SHqst_l(sht_l, f, g, h, qLM, sLM, tLM, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif

   end subroutine spat_to_qst
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor(sh, f, g, fLM, gLM, lcut, use_gpu)

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      real(cp),    intent(inout) :: f(*)
      real(cp),    intent(inout) :: g(*)
      integer,     intent(in) :: lcut
      logical, optional, intent(in) :: use_gpu

      !-- Output variables
      complex(cp), intent(out) :: fLM(*)
      complex(cp), intent(out) :: gLM(*)

#ifdef WITH_OMP_GPU
      logical :: loc_use_gpu

      if( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      else
         loc_use_gpu = .false.
      end if
#endif

#ifdef WITH_OMP_GPU
      if(loc_use_gpu) then
         !$omp target data use_device_addr(f, g, fLM, gLM)
         call cu_spat_to_SHsphtor(sh, f, g, fLM, gLM, lcut)
         !$omp end target data
      else
         call spat_to_SHsphtor_l(sh, f, g, fLM, gLM, lcut)
      end if
#else
      call spat_to_SHsphtor_l(sh, f, g, fLM, gLM, lcut)
#endif

#ifdef WITH_OMP_GPU
#ifdef EXPLICITE_SYNCHRO
      call hipCheck(hipDeviceSynchronize())
#endif
#endif

   end subroutine spat_to_sphertor
!------------------------------------------------------------------------------
   subroutine axi_to_spat(fl_ax, f)

      !-- Input field
      complex(cp), intent(inout) :: fl_ax(l_max+1) !-- Axi-sym toroidal

      !-- Output field on grid
      real(cp), intent(out) :: f(:)

      !-- Local arrays
      complex(cp) :: tmp(nlat_padded)

      call SH_to_spat_ml(sht_l_single, 0, fl_ax, tmp, l_max)
      f(:)=real(tmp(:))

   end subroutine axi_to_spat
!------------------------------------------------------------------------------
   subroutine toraxi_to_spat(fl_ax, ft, fp)

      !-- Input field
      complex(cp), intent(inout) :: fl_ax(l_max+1) !-- Axi-sym toroidal

      !-- Output fields on grid
      real(cp), intent(out) :: ft(:)
      real(cp), intent(out) :: fp(:)

      !-- Local arrays
      complex(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)

      call SHtor_to_spat_ml(sht_l_single, 0, fl_ax, tmpt, tmpp, l_max)
      ft(:)=real(tmpt(:))
      fp(:)=real(tmpp(:))

   end subroutine toraxi_to_spat
!------------------------------------------------------------------------------
end module sht
