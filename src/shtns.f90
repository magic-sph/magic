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
       &                 minc, lm_max, lmP_max, nlat_padded
   use horizontal_data, only: dLh
   use radial_data, only: nRstart, nRstop
   use parallel_mod

   implicit none

   include "shtns.f03"

   private

   public :: initialize_sht, scal_to_spat, scal_to_grad_spat, torpol_to_spat, &
   &         scal_to_SH, spat_to_sphertor, torpol_to_spat_IC, axi_to_spat,    &
   &         torpol_to_curl_spat_IC, spat_to_SH_axi, spat_to_qst,             &
   &         sphtor_to_spat, toraxi_to_spat, finalize_sht, torpol_to_spat_single

   type(c_ptr), public :: sht_l, sht_lP, sht_l_single, sht_lP_single

   logical :: l_batched_sh

contains

   subroutine initialize_sht(l_scrambled_theta, l_batched)

      !-- Input variable
      logical, intent(in) :: l_batched

      !-- Output variable
      logical, intent(out) :: l_scrambled_theta

      !-- Local variables
      integer :: norm, layout
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
      if ( .not. l_batched_sh ) then
         !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
         layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
      else
         !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS
         layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS
      end if
#else
      !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS
      layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS
#endif
      eps_polar = 1.e-10_cp

      sht_l = shtns_create(l_max, m_max/minc, minc, norm)
      if ( l_batched_sh ) call shtns_set_batch(sht_l, howmany, dist)
      call shtns_set_grid(sht_l, layout, eps_polar, n_theta_max, n_phi_max)

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
      call shtns_robert_form(sht_l, 1) ! Use Robert's form

      if ( rank == 0 ) then
         call shtns_verbose(0)
         write(output_unit,*) ''
      end if

      if ( l_batched_sh ) then
         sht_l_single = shtns_create(l_max, m_max/minc, minc, norm)
         call shtns_set_grid(sht_l_single, layout, eps_polar, n_theta_max, n_phi_max)
         call shtns_robert_form(sht_l_single, 1) ! Use Robert's form
      else
         sht_l_single = sht_l
      end if

      sht_lP = shtns_create(l_max+1, m_max/minc, minc, norm)
      dist = int(lmP_max, kind=8)
      if ( l_batched_sh ) call shtns_set_batch(sht_lP, howmany, dist)
      call shtns_set_grid(sht_lP, layout, eps_polar, n_theta_max, n_phi_max)
      call shtns_robert_form(sht_lP, 1) ! Use Robert's form


      if ( l_batched_sh ) then
         sht_lP_single = shtns_create(l_max+1, m_max/minc, minc, norm)
         call shtns_set_grid(sht_lP_single, layout, eps_polar, n_theta_max, n_phi_max)
         call shtns_robert_form(sht_lP_single, 1) ! Use Robert's form
      else
         sht_lP_single = sht_lP
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
      call shtns_unset_grid(sht_lP)
      call shtns_destroy(sht_lP)
      if ( l_batched_sh ) then
         call shtns_unset_grid(sht_l_single)
         call shtns_destroy(sht_l_single)
         call shtns_unset_grid(sht_lP_single)
         call shtns_destroy(sht_lP_single)
      end if

   end subroutine finalize_sht
!------------------------------------------------------------------------------
   subroutine scal_to_spat(sh, Slm, fieldc, lcut)
      ! transform a spherical harmonic field into grid space

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      complex(cp), intent(inout) :: Slm(*)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: fieldc(*)

      call SH_to_spat_l(sh, Slm, fieldc, lcut)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc, lcut)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variables
      complex(cp), intent(inout) :: Slm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(*)
      real(cp), intent(out) :: gradpc(*)

      call SHsph_to_spat_l(sht_l, Slm, gradtc, gradpc, lcut)

   end subroutine scal_to_grad_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(inout) :: Wlm(*), dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(*)
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

      call SHqst_to_spat_l(sht_l, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

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
   subroutine sphtor_to_spat(sh, dWlm, Zlm, vtc, vpc, lcut)

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      complex(cp), intent(inout) :: dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

      call SHsphtor_to_spat_l(sh, dWlm, Zlm, vtc, vpc, lcut)

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
   subroutine scal_to_SH(sh, f, fLM, lcut)

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      real(cp),    intent(inout) :: f(*)
      integer,     intent(in) :: lcut

      !-- Output variable
      complex(cp), intent(out) :: fLM(*)

      call spat_to_SH_l(sh, f, fLM, lcut+1)

   end subroutine scal_to_SH
!------------------------------------------------------------------------------
   subroutine spat_to_qst(f, g, h, qLM, sLM, tLM, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f(*)
      real(cp), intent(inout) :: g(*)
      real(cp), intent(inout) :: h(*)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: qLM(*)
      complex(cp), intent(out) :: sLM(*)
      complex(cp), intent(out) :: tLM(*)

      call spat_to_SHqst_l(sht_lP, f, g, h, qLM, sLM, tLM, lcut+1)

   end subroutine spat_to_qst
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor(sh, f, g, fLM, gLM, lcut)

      !-- Input variables
      type(c_ptr), intent(in) :: sh
      real(cp),    intent(inout) :: f(*)
      real(cp),    intent(inout) :: g(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLM(*)
      complex(cp), intent(out) :: gLM(*)

      call spat_to_SHsphtor_l(sh, f, g, fLM, gLM, lcut+1)

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
   subroutine spat_to_SH_axi(f, fLM)

      real(cp), intent(in) :: f(:)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp(nlat_padded)
      complex(cp) :: tmpLM(size(fLM))

      tmp(:)=cmplx(f(:),0.0_cp,kind=cp)
      if ( size(fLM) == l_max+2 ) then
         call spat_to_SH_ml(sht_lP_single, 0, tmp, tmpLM, l_max+1)
      else if ( size(fLM) == l_max+1 ) then
         call spat_to_SH_ml(sht_l_single, 0, tmp, tmpLM, l_max)
      end if
      fLM(:)=real(tmpLM(:))

   end subroutine spat_to_SH_axi
!------------------------------------------------------------------------------
end module sht
