module sht

   use iso_fortran_env, only: output_unit
   use iso_c_binding, only: c_ptr
   use precision_mod, only: cp
   use blocking, only: st_map
   use constants, only: one, zero
   use truncation, only: m_max, l_max, n_theta_max, minc, lm_max
   use horizontal_data, only: dLh
   use parallel_mod
   use shtransforms

   implicit none

   private

   public :: initialize_sht, scal_to_spat, scal_to_grad_spat,             &
   &         torpol_to_spat, scal_to_SH, spat_to_sphertor,                &
   &         torpol_to_spat_IC, torpol_to_curl_spat_IC, axi_to_spat,      &
   &         spat_to_qst, sphtor_to_spat, toraxi_to_spat, finalize_sht,   &
   &         torpol_to_spat_single

   type(c_ptr), public :: sht_l, sht_lP, sht_l_single, sht_lP_single

contains

   subroutine initialize_sht(l_scrambled_theta,l_batched_shts)

      !-- Input variable
      logical, intent(in) :: l_batched_shts ! Dummy, only for API compatibility

      !-- Output variable
      logical, intent(out) :: l_scrambled_theta

      l_scrambled_theta=.true.
      call initialize_transforms()

   end subroutine initialize_sht
!------------------------------------------------------------------------------
   subroutine finalize_sht()

      call finalize_transforms()

   end subroutine finalize_sht
!------------------------------------------------------------------------------
   subroutine scal_to_spat(sh, Slm, fieldc, lcut)
      ! transform a spherical harmonic field into grid space

      !-- Input variables
      type(c_ptr), intent(in) :: sh ! Dummy: only for API compatibility with SHTns
      complex(cp), intent(in) :: Slm(*)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: fieldc(*)

      call native_sph_to_spat(Slm, fieldc, lcut)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc, lcut)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variables
      complex(cp), intent(in) :: Slm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(*)
      real(cp), intent(out) :: gradpc(*)

      call native_sph_to_grad_spat(Slm, gradtc, gradpc, lcut)

   end subroutine scal_to_grad_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(*), dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(*)
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

      call native_qst_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

   end subroutine torpol_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_single(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(*), dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(*)
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = st_map%lm2l(lm)
         if ( l <= lcut ) then
            Qlm(lm) = dLh(lm) * Wlm(lm)
         else
            Qlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call native_qst_to_spat(Qlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

   end subroutine torpol_to_spat_single
!------------------------------------------------------------------------------
   subroutine sphtor_to_spat(sh, dWlm, Zlm, vtc, vpc, lcut)

      !-- Input variables
      type(c_ptr), intent(in) :: sh ! Dummy: only for API compatibility with SHTns
      complex(cp), intent(in) :: dWlm(*), Zlm(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vtc(*)
      real(cp), intent(out) :: vpc(*)

      call native_sphtor_to_spat(dWlm, Zlm, vtc, vpc, lcut)

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

      call native_qst_to_spat(Qlm, Slm, Tlm, cbr, cbt, cbp, l_max)

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

      call native_qst_to_spat(Qlm, Slm, Tlm, Br, Bt, Bp, l_max)

   end subroutine torpol_to_spat_IC
!------------------------------------------------------------------------------
   subroutine scal_to_SH(sh, f, fLM, lcut)

      type(c_ptr), intent(in) :: sh ! Dummy, only here for API compatibility
      real(cp),    intent(inout) :: f(*)
      integer,     intent(in) :: lcut
      complex(cp), intent(out) :: fLM(*)

      call native_spat_to_sph(f, fLM, lcut)

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

      call native_spat_to_sph(f, qLM, lcut)
      call native_spat_to_sph_tor(g, h, sLM, tLM, lcut)

   end subroutine spat_to_qst
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor(sh, f, g, fLM, gLM, lcut)

      !-- Input variables
      type(c_ptr), intent(in) :: sh ! Dummy, only here for API compatibility
      real(cp),    intent(inout) :: f(*)
      real(cp),    intent(inout) :: g(*)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLM(*)
      complex(cp), intent(out) :: gLM(*)

      call native_spat_to_sph_tor(f, g, fLM, gLM, lcut)

   end subroutine spat_to_sphertor
!------------------------------------------------------------------------------
   subroutine axi_to_spat(fl_ax, f)

      !-- Input field
      complex(cp), intent(in) :: fl_ax(l_max+1) !-- Axi-sym field

      !-- Output field on grid
      real(cp), intent(out) :: f(:)

      call native_axi_to_spat(fl_ax, f)

   end subroutine axi_to_spat
!------------------------------------------------------------------------------
   subroutine toraxi_to_spat(fl_ax, ft, fp, lcut)

      !-- Input field
      integer, intent(in) :: lcut ! cut-off spherical harmonic degree
      complex(cp), intent(in) :: fl_ax(l_max+1) !-- Axi-sym toroidal

      !-- Output fields on grid
      real(cp), intent(out) :: ft(:)
      real(cp), intent(out) :: fp(:)

      call native_toraxi_to_spat(fl_ax, ft, fp, lcut)

   end subroutine toraxi_to_spat
!------------------------------------------------------------------------------
end module sht
