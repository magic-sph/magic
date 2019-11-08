module shtns

   use iso_fortran_env, only: output_unit
   use precision_mod, only: cp
   use blocking, only: st_map
   use constants, only: ci, one, zero
   use truncation, only: m_max, l_max, n_theta_max, n_phi_max, &
       &                 minc, lm_max, lmP_max
   use horizontal_data, only: dLh, D_m, O_sin_theta_E2, O_sin_theta
   use parallel_mod

   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, scal_to_grad_spat, pol_to_grad_spat, &
   &         torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
   &         torpol_to_dphspat, spat_to_SH, spat_to_sphertor,               &
   &         torpol_to_spat_IC, torpol_to_curl_spat_IC, spat_to_SH_axi,     &
   &         axi_to_spat, spat_to_qst

contains

   subroutine init_shtns()

      integer :: norm

      if ( rank == 0 ) then
         call shtns_verbose(1)
      end if

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max/minc, minc, norm)
      call shtns_precompute(SHT_QUICK_INIT, SHT_PHI_CONTIGUOUS, &
           &                1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(0)

      if ( rank == 0 ) then
         call shtns_verbose(0)
         write(output_unit,*) ''
      end if

      call shtns_set_size(l_max+1, m_max/minc, minc, norm)
      call shtns_precompute(SHT_QUICK_INIT, SHT_PHI_CONTIGUOUS, &
           &                1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(1)

      call shtns_load_cfg(0)

   end subroutine
!------------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc, lcut)
      ! transform a spherical harmonic field into grid space

      !-- Input variables
      complex(cp), intent(in) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: fieldc(n_phi_max, n_theta_max)

      call shtns_SH_to_spat_l(Slm, fieldc, lcut)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc, lcut)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variables
      complex(cp), intent(in) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      call shtns_sph_to_spat_l(Slm, gradtc, gradpc, lcut)

   end subroutine scal_to_grad_spat
!------------------------------------------------------------------------------
   subroutine pol_to_grad_spat(Slm, gradtc, gradpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = st_map%lm2l(lm)
         if ( l <= lcut ) then
            Qlm(lm) = dLh(lm) * Slm(lm)
         else
            Qlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call shtns_sph_to_spat_l(Qlm, gradtc, gradpc, lcut)

   end subroutine pol_to_grad_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(lm_max), dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: vtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: vpc(n_phi_max, n_theta_max)

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

      call shtns_qst_to_spat_l(Qlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

   end subroutine torpol_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_IC(r, r_ICB, dBlm, ddBlm, Jlm, dJlm, &
              &                      cbr, cbt, cbp)
      !
      ! This is a QST transform that contains the transform for the
      ! inner core to compute the three components of the curl of B.
      !

      !-- Input variables
      real(cp),    intent(in) :: r, r_ICB
      complex(cp), intent(in) :: dBlm(lm_max), ddBlm(lm_max)
      complex(cp), intent(in) :: Jlm(lm_max), dJlm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: cbr(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cbt(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cbp(n_phi_max, n_theta_max)

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

      call shtns_qst_to_spat(Qlm, Slm, Tlm, cbr, cbt, cbp)

   end subroutine torpol_to_curl_spat_IC
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_IC(r, r_ICB, Wlm, dWlm, Zlm, Br, Bt, Bp)
      !
      ! This is a QST transform that contains the transform for the
      ! inner core.
      !

      !-- Input variables
      real(cp),    intent(in) :: r, r_ICB
      complex(cp), intent(in) :: Wlm(lm_max), dWlm(lm_max), Zlm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: Br(n_phi_max, n_theta_max)
      real(cp), intent(out) :: Bt(n_phi_max, n_theta_max)
      real(cp), intent(out) :: Bp(n_phi_max, n_theta_max)

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

      call shtns_qst_to_spat(Qlm, Slm, Tlm, Br, Bt, Bp)

   end subroutine torpol_to_spat_IC
!------------------------------------------------------------------------------
   subroutine torpol_to_dphspat(dWlm, Zlm, dvtdp, dvpdp, lcut)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !

      !-- Input variables
      complex(cp), intent(in) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: dvtdp(n_phi_max, n_theta_max)
      real(cp), intent(out) :: dvpdp(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: Slm(lm_max), Tlm(lm_max)
      integer :: lm, it, ip, l
      real(cp) :: m

      !$omp parallel do default(shared) private(lm, m, l)
      do lm = 1, lm_max
         l = st_map%lm2l(lm)
         m = D_m(lm)
         if ( l <= lcut ) then
            Slm(lm) = ci*m*dWlm(lm)
            Tlm(lm) = ci*m*Zlm(lm)
         else
            Slm(lm) = 0.0_cp
            Tlm(lm) = 0.0_cp
         end if
      end do
      !$omp end parallel do

      call shtns_sphtor_to_spat_l(Slm, Tlm, dvtdp, dvpdp, lcut)

      !$omp parallel do default(shared) private(it,ip)
      do it=1, n_theta_max
         do ip=1, n_phi_max
            dvtdp(ip, it) = dvtdp(ip, it) * O_sin_theta_E2(it)
            dvpdp(ip, it) = dvpdp(ip, it) * O_sin_theta_E2(it)
         end do
      end do
      !$omp end parallel do

   end subroutine torpol_to_dphspat
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat(Qlm, cvrc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Qlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: dQlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = st_map%lm2l(lm)
         if ( l <= lcut ) then
            dQlm(lm) = dLh(lm) * Qlm(lm)
         else
            dQlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call shtns_SH_to_spat_l(dQlm, cvrc, lcut)

   end subroutine pol_to_curlr_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat(or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, cvpc, &
              &                   lcut)

      !-- Input variables
      complex(cp), intent(in) :: Blm(lm_max), ddBlm(lm_max)
      complex(cp), intent(in) :: Jlm(lm_max), dJlm(lm_max)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cvtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cvpc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = st_map%lm2l(lm)
         if ( l <= lcut ) then
            Qlm(lm) = dLh(lm) * Jlm(lm)
            Tlm(lm) = or2 * dLh(lm) * Blm(lm) - ddBlm(lm)
         else
            Qlm(lm) = zero
            Tlm(lm) = zero
         end if
      end do
      !
      !$omp end parallel do

      call shtns_qst_to_spat_l(Qlm, dJlm, Tlm, cvrc, cvtc, cvpc, lcut)

   end subroutine torpol_to_curl_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM, lcut)

      real(cp), intent(in) :: f(n_phi_max, n_theta_max)
      integer,  intent(in) :: lcut
      complex(cp), intent(out) :: fLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sh_l(f, fLM, lcut+1)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH
!------------------------------------------------------------------------------
   subroutine spat_to_qst(f, g, h, qLM, sLM, tLM, lcut)

      !-- Input variables
      real(cp), intent(in) :: f(n_phi_max,n_theta_max)
      real(cp), intent(in) :: g(n_phi_max,n_theta_max)
      real(cp), intent(in) :: h(n_phi_max,n_theta_max)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: qLM(lmP_max)
      complex(cp), intent(out) :: sLM(lmP_max)
      complex(cp), intent(out) :: tLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_qst_l(f, g, h, qLM, sLM, tLM, lcut+1)
      call shtns_load_cfg(0)

   end subroutine spat_to_qst
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor(f, g, fLM, gLM, lcut)

      !-- Input variables
      real(cp), intent(in) :: f(n_phi_max,n_theta_max)
      real(cp), intent(in) :: g(n_phi_max,n_theta_max)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLM(lmP_max)
      complex(cp), intent(out) :: gLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sphtor_l(f, g, fLM, gLM, lcut+1)
      call shtns_load_cfg(0)

   end subroutine spat_to_sphertor
!------------------------------------------------------------------------------
   subroutine axi_to_spat(fl_ax, f)

      real(cp), intent(in) :: fl_ax(l_max+1)
      real(cp), intent(out) :: f(n_theta_max)

      !-- Local arrays
      complex(cp) :: tmp(n_theta_max)
      complex(cp) :: tmp_ax(l_max+1)

      tmp_ax(:)=cmplx(fl_ax(:),0.0_cp,kind=cp)
      call shtns_sh_to_spat_ml(0, tmp_ax, tmp, l_max)
      f(:)=real(tmp(:))

   end subroutine axi_to_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH_axi(f, fLM)

      real(cp), intent(in) :: f(n_theta_max)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp(n_theta_max)
      complex(cp) :: tmpLM(size(fLM))

      if ( size(fLM) == l_max+2 ) call shtns_load_cfg(1)
      tmp(:)=cmplx(f(:),0.0_cp,kind=cp)
      call shtns_spat_to_sh_ml(0, tmp, tmpLM, size(fLM)-1)
      if ( size(fLM) == l_max+2 ) call shtns_load_cfg(0)
      fLM(:)=real(tmpLM(:))

   end subroutine spat_to_SH_axi
!------------------------------------------------------------------------------
end module shtns
