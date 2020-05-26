module shtns

   use iso_fortran_env, only: output_unit
   use precision_mod, only: cp, MPI_DEF_REAL
   use blocking, only: st_map
   use constants, only: ci, one, zero
   use truncation, only: m_max, l_max, n_theta_max, n_phi_max, &
       &                 minc, lm_max, lmP_max, n_lm_loc, n_theta_loc, &
       &                 n_m_loc, n_lmP_loc, n_m_max, dist_m,          &
       &                 nThetaStart, nThetaStop, coord_m, dist_theta
   use horizontal_data, only: dLh_loc, dLh, O_sin_theta_E2, O_sin_theta
   use parallel_mod
   use fft, only: fft_phi_loc
   use LMmapping, only: map_dist_st, map_glbl_st
   use mpi_thetap_mod

   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, scal_to_grad_spat, pol_to_grad_spat, &
   &         torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
   &         torpol_to_dphspat, spat_to_SH, spat_to_sphertor,               &
   &         torpol_to_spat_IC, torpol_to_curl_spat_IC, spat_to_SH_axi,     &
   &         axi_to_spat, spat_to_qst,                                      &
   &         spat_to_SH_dist, spat_to_qst_dist, spat_to_sphertor_dist,      &
   &         spat_to_SH_axi_dist, scal_to_spat_dist, scal_to_grad_spat_dist,&
   &         torpol_to_spat_dist, torpol_to_curl_spat_dist,                 &
   &         pol_to_grad_spat_dist, torpol_to_dphspat_dist,                 &
   &         pol_to_curlr_spat_dist

   public :: scal_to_hyb, scal_to_grad_hyb, torpol_to_hyb, torpol_to_curl_hyb, &
   &         pol_to_grad_hyb, torpol_to_dphhyb, pol_to_curlr_hyb, hyb_to_SH,   &
   &         hyb_to_qst, hyb_to_sphertor

contains

   subroutine init_shtns()

      integer :: norm

      if ( l_master_rank ) then
         call shtns_verbose(1)
      end if

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max/minc, minc, norm)
      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, &
           &                1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(0)

      if ( l_master_rank ) then
         call shtns_verbose(0)
         write(output_unit,*) ''
      end if

      call shtns_set_size(l_max+1, m_max/minc, minc, norm)
      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, &
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
         m = st_map%lm2m(lm)
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


!------------------------------------------------------------------------------
!   
!   Θ-Distributed Forward Transforms
!   
!------------------------------------------------------------------------------

   subroutine scal_to_spat_dist(fLM_loc, fieldc_loc, lcut)
      !   
      !   Transform the spherical harmonic coefficients Qlm into its spatial 
      !   representation Vr
      !  
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      
      !-- Input variables
      complex(cp),  intent(inout) :: fLM_loc(n_lm_loc)
      
      !-- Output variables
      real(cp),     intent(out)   :: fieldc_loc(n_phi_max, n_theta_loc)
      integer,      intent(in)    :: lcut

      
      !-- Local variables
      !complex(cp) :: tmp1(n_lm_loc)
      complex(cp) :: fL_loc(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        if (m>lcut) then
           fL_loc(:,i)=zero
           cycle
        end if
        l_lm = map_dist_st%lm2(m, m)
        u_lm = map_dist_st%lm2(l_max, m)
        !tmp1=fLM_loc

        call shtns_sh_to_spat_ml(m/minc, fLM_loc(l_lm:u_lm), fL_loc(:,i),lcut)

        !block
        !
        !   use legendre_spec_to_grid, only: sh_to_spat_ml
        !   complex(cp) :: tmp(n_theta_max)
        !
        !  call sh_to_spat_ml(m, tmp1(l_lm:u_lm), tmp, lcut)
        !
        !   if( maxval(abs(tmp-fL_loc(:,i))) >= 2.0_cp*epsilon(1.0_cp) ) then
        !      print*, maxval(abs(tmp-fL_loc(:,i)))
        !   end if
        !
        !end block

      end do
      !$omp end parallel do
      
      call transform_m2phi(fL_loc, fieldc_loc)
      
   end subroutine scal_to_spat_dist
   
   !------------------------------------------------------------------------------
   subroutine scal_to_grad_spat_dist(Slm, gradtc, gradpc, lcut)
      !
      !   Transform a scalar spherical harmonic field into it's gradient
      !   on the grid
      !
      !   Author: Rafael Lago, MPCDF, April 2020
      !

      !-- Input variables
      complex(cp), intent(in) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_loc)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      integer :: i, l_lm, u_lm, m

      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         call shtns_sph_to_spat_ml(m/minc, Slm(l_lm:u_lm), fL(:,i), gL(:,i), lcut)
      end do
      !$omp end parallel do
      
      call transform_m2phi(fL, gradtc)
      call transform_m2phi(gL, gradpc)

   end subroutine scal_to_grad_spat_dist
   
   !------------------------------------------------------------------------------
   subroutine torpol_to_spat_dist(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(n_lm_loc), dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(n_phi_max, n_theta_loc)
      real(cp), intent(out) :: vtc(n_phi_max, n_theta_loc)
      real(cp), intent(out) :: vpc(n_phi_max, n_theta_loc)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      complex(cp) :: hL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            hL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Qlm(m:l_max) = dLh_loc(l_lm:u_lm) * Wlm(l_lm:u_lm)
         !if (lcut<l_max) Qlm(lcut+1:u_lm) = zero
         call shtns_qst_to_spat_ml(m/minc, Qlm(m:l_max), dWlm(l_lm:u_lm), &
              &                    Zlm(l_lm:u_lm), fL(:,i), gL(:,i), hL(:,i), lcut)
      end do
      !$omp end parallel do
      
      call transform_m2phi(fL, vrc)
      call transform_m2phi(gL, vtc)
      call transform_m2phi(hL, vpc)

   end subroutine torpol_to_spat_dist
   
   !------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_dist(or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, &
              &                        cvpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc)
      complex(cp), intent(in) :: Jlm(n_lm_loc), dJlm(n_lm_loc)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_loc)
      real(cp), intent(out) :: cvtc(n_phi_max, n_theta_loc)
      real(cp), intent(out) :: cvpc(n_phi_max, n_theta_loc)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max), Tlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      complex(cp) :: hL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            hL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Qlm(m:l_max) = dLh_loc(l_lm:u_lm) * Jlm(l_lm:u_lm)
         Tlm(m:l_max) = or2 * dLh_loc(l_lm:u_lm) * Blm(l_lm:u_lm) - ddBlm(l_lm:u_lm)
         call shtns_qst_to_spat_ml(m/minc, Qlm(m:l_max), dJlm(l_lm:u_lm), &
              &                    Tlm(m:l_max), fL(:,i), gL(:,i), hL(:,i), lcut)
      end do
      !$omp end parallel do
      
      call transform_m2phi(fL, cvrc)
      call transform_m2phi(gL, cvtc)
      call transform_m2phi(hL, cvpc)

   end subroutine torpol_to_curl_spat_dist
   
   !------------------------------------------------------------------------------
   subroutine pol_to_grad_spat_dist(Slm, gradtc, gradpc, lcut)

       !-- Input variables
      complex(cp), intent(in) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_loc)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_loc)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Qlm(m:l_max) = dLh_loc(l_lm:u_lm) * Slm(l_lm:u_lm)
         call shtns_sph_to_spat_ml(m/minc, Qlm(m:l_max), fL(:,i), gL(:,i), lcut)
      end do
      !$omp end parallel do
      
      call transform_m2phi(fL, gradtc)
      call transform_m2phi(gL, gradpc)

   end subroutine pol_to_grad_spat_dist
   
   
   !------------------------------------------------------------------------------
   subroutine torpol_to_dphspat_dist(dWlm, Zlm, dvtdp, dvpdp, lcut)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !

      !-- Input variables
      complex(cp), intent(in) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: dvtdp(n_phi_max, nThetaStart:nThetaStop) ! Careful with dimensions here!
      real(cp), intent(out) :: dvpdp(n_phi_max, nThetaStart:nThetaStop) ! Careful with dimensions here!

      !-- Local variables
      complex(cp) :: Slm(0:l_max), Tlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m, it, ip
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Slm(m:l_max) = ci*m*dWlm(l_lm:u_lm)
         Tlm(m:l_max) = ci*m*Zlm(l_lm:u_lm)
         call shtns_sphtor_to_spat_ml(m/minc, Slm(m:l_max), Tlm(m:l_max), &
              &                       fL(:,i), gL(:,i), lcut)
      end do
      !$omp end parallel do
      
      call transform_m2phi(fL, dvtdp)
      call transform_m2phi(gL, dvpdp)
      
      !$omp parallel do default(shared) private(it,ip)
      do it=nThetaStart, nThetaStop
         do ip=1, n_phi_max
            dvtdp(ip, it) = dvtdp(ip, it) * O_sin_theta_E2(it)
            dvpdp(ip, it) = dvpdp(ip, it) * O_sin_theta_E2(it)
         end do
      end do
      !$omp end parallel do

   end subroutine torpol_to_dphspat_dist
   
   !------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_dist(Qlm, cvrc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Qlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: dQlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         dQlm(m:l_max) = dLh_loc(l_lm:u_lm) * Qlm(l_lm:u_lm)
         call shtns_SH_to_spat_ml(m/minc, dQlm(m:l_max), fL(:,i), lcut)
      end do
      !$omp end parallel do
      
      call transform_m2phi(fL, cvrc)

   end subroutine pol_to_curlr_spat_dist
   
   
   
!------------------------------------------------------------------------------
!   
!   Θ-Distributed Backward Transforms
!   
!------------------------------------------------------------------------------
  
   !----------------------------------------------------------------------------
   subroutine spat_to_SH_dist(f_loc, fLMP_loc, lcut)

      !-- Input variables
      real(cp),    intent(inout) :: f_loc(n_phi_max,n_theta_loc)
      integer,     intent(in)    :: lcut

      !-- Output variable
      complex(cp), intent(out)   :: fLMP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      integer :: m, l_lm, u_lm, i
      
      call shtns_load_cfg(1) ! l_max + 1
      
      call transform_phi2m(f_loc, fL_loc)
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            fLMP_loc(l_lm:u_lm)=zero
            cycle
         end if 
         call shtns_spat_to_sh_ml(m/minc,fL_loc(:,i),fLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do
      
      call shtns_load_cfg(0) ! l_max

   end subroutine spat_to_SH_dist
   
   !------------------------------------------------------------------------------
   subroutine spat_to_qst_dist(f_loc, g_loc, h_loc, qLMP_loc, sLMP_loc, tLMP_loc, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f_loc(n_phi_max,n_theta_loc)
      real(cp), intent(inout) :: g_loc(n_phi_max,n_theta_loc)
      real(cp), intent(inout) :: h_loc(n_phi_max,n_theta_loc)
      integer,  intent(in)    :: lcut

      !-- Output variables
      complex(cp), intent(out) :: qLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: sLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: tLMP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  gL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  hL_loc(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m

      call shtns_load_cfg(1)
      
      !>@TODO Vectorial FFT and transpose (f,g,h at once)
      call transform_phi2m(f_loc, fL_loc)
      call transform_phi2m(g_loc, gL_loc)
      call transform_phi2m(h_loc, hL_loc)
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            qLMP_loc(l_lm:u_lm)=zero
            sLMP_loc(l_lm:u_lm)=zero
            tLMP_loc(l_lm:u_lm)=zero
            cycle
         end if
         call shtns_spat_to_qst_ml(m/minc,fL_loc(:,i),gL_loc(:,i),hL_loc(:,i),&
              &                    qLMP_loc(l_lm:u_lm),sLMP_loc(l_lm:u_lm),   &
              &                    tLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do
      
      call shtns_load_cfg(0)

   end subroutine spat_to_qst_dist
   
   !------------------------------------------------------------------------------
   subroutine spat_to_sphertor_dist(f_loc, g_loc, fLMP_loc, gLMP_loc, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f_loc(n_phi_max,n_theta_loc)
      real(cp), intent(inout) :: g_loc(n_phi_max,n_theta_loc)
      integer,  intent(in)    :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: gLMP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  gL_loc(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m

      call shtns_load_cfg(1)
      
      !>@TODO Vectorial FFT and transpose (f,g,h at once)
      call transform_phi2m(f_loc, fL_loc)
      call transform_phi2m(g_loc, gL_loc)
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            fLMP_loc(l_lm:u_lm)=zero
            gLMP_loc(l_lm:u_lm)=zero
            cycle
         end if
         call shtns_spat_to_sphtor_ml(m/minc,fL_loc(:,i),gL_loc(:,i),&
              &                       fLMP_loc(l_lm:u_lm),           &
              &                       gLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do
      
      call shtns_load_cfg(0)

   end subroutine spat_to_sphertor_dist
   
   !------------------------------------------------------------------------------
   subroutine spat_to_SH_axi_dist(f, fLM)

      real(cp), intent(in)  :: f(nThetaStart:nThetaStop)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp_c(n_theta_max)
      real(cp)    :: tmp_r(n_theta_max)
      complex(cp) :: tmpLM(size(fLM))
      integer     :: ierr

      if ( size(fLM) == l_max+2 ) call shtns_load_cfg(1)
      tmp_r(nThetaStart:nThetaStop)=f(nThetaStart:nThetaStop)

#ifdef WITH_MPI
      !@TODO I think that this needs to be done only for coord_theta=0; if so, this is a gatherv 
      !      instead of an allgatherv
      !@TODO it may be determined beforehand if the number of theta points is the same in every 
      !      rank. If yes, then this would be a gather instead of a gatherall
      if ( n_ranks_theta>1 ) then
         call MPI_ALLGATHERV(MPI_IN_PLACE, 0, 0, tmp_r, dist_theta(:,0), &
              &              dist_theta(:,1)-1, MPI_DEF_REAL, comm_theta, ierr)
      end if
#endif
      tmp_c = cmplx(tmp_r, 0.0, kind=cp)
      call shtns_spat_to_sh_ml(0, tmp_c, tmpLM, size(fLM)-1)
      
      if ( size(fLM) == l_max+2 ) call shtns_load_cfg(0)
      fLM(:)=real(tmpLM(:))

   end subroutine spat_to_SH_axi_dist
!------------------------------------------------------------------------------
   

!--------------------------------
!
!-- Legendre only
!
!--------------------------------
   subroutine scal_to_hyb(fLM_loc, field_hyb, lcut)
      !   
      !   Transform the spherical harmonic coefficients Qlm into its hybrid
      !   representation Vr(n_theta,n_m)
      !  
      
      !-- Input variables
      integer,     intent(in) :: lcut
      complex(cp), intent(in) :: fLM_loc(n_lm_loc)
      
      !-- Output variables
      complex(cp), intent(out) :: field_hyb(n_theta_max, n_m_loc)

      !-- Local variables
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        if (m>lcut) then
           field_hyb(:,i)=zero
           cycle
        end if
        l_lm = map_dist_st%lm2(m, m)
        u_lm = map_dist_st%lm2(l_max, m)

        call shtns_sh_to_spat_ml(m/minc, fLM_loc(l_lm:u_lm), field_hyb(:,i),lcut)

      end do
      !$omp end parallel do
      
   end subroutine scal_to_hyb
!------------------------------------------------------------------------------
   subroutine scal_to_grad_hyb(Slm, gradtL, gradpL, lcut)
      !
      !   Transform a scalar spherical harmonic field into it's gradient
      !   on hybrid (n_theta,n_m) space
      !

      !-- Input variables
      complex(cp), intent(in) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: gradtL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: gradpL(n_theta_max,n_m_loc)
      
      !-- Local variables
      
      integer :: i, l_lm, u_lm, m

      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            gradtL(:,i)=zero
            gradpL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         call shtns_sph_to_spat_ml(m/minc, Slm(l_lm:u_lm), gradtL(:,i), gradpL(:,i), &
              &                    lcut)
      end do
      !$omp end parallel do
      
   end subroutine scal_to_grad_hyb
!------------------------------------------------------------------------------
   subroutine torpol_to_hyb(Wlm, dWlm, Zlm, fL, gL, hL, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(n_lm_loc), dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: gL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: hL(n_theta_max,n_m_loc)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            hL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Qlm(m:l_max) = dLh_loc(l_lm:u_lm) * Wlm(l_lm:u_lm)
         !if (lcut<l_max) Qlm(lcut+1:u_lm) = zero
         call shtns_qst_to_spat_ml(m/minc, Qlm(m:l_max), dWlm(l_lm:u_lm), &
              &                    Zlm(l_lm:u_lm), fL(:,i), gL(:,i), hL(:,i), lcut)
      end do
      !$omp end parallel do
      
   end subroutine torpol_to_hyb
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_hyb(or2, Blm, ddBlm, Jlm, dJlm, fL, gL, hL, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc)
      complex(cp), intent(in) :: Jlm(n_lm_loc), dJlm(n_lm_loc)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: gL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: hL(n_theta_max,n_m_loc)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max), Tlm(0:l_max)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            hL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Qlm(m:l_max) = dLh_loc(l_lm:u_lm) * Jlm(l_lm:u_lm)
         Tlm(m:l_max) = or2 * dLh_loc(l_lm:u_lm) * Blm(l_lm:u_lm) - ddBlm(l_lm:u_lm)
         call shtns_qst_to_spat_ml(m/minc, Qlm(m:l_max), dJlm(l_lm:u_lm), &
              &                    Tlm(m:l_max), fL(:,i), gL(:,i), hL(:,i), lcut)
      end do
      !$omp end parallel do
      
   end subroutine torpol_to_curl_hyb
!------------------------------------------------------------------------------
   subroutine pol_to_grad_hyb(Slm, fL, gL, lcut)

       !-- Input variables
      complex(cp), intent(in) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: gL(n_theta_max,n_m_loc)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Qlm(m:l_max) = dLh_loc(l_lm:u_lm) * Slm(l_lm:u_lm)
         call shtns_sph_to_spat_ml(m/minc, Qlm(m:l_max), fL(:,i), gL(:,i), lcut)
      end do
      !$omp end parallel do
      
   end subroutine pol_to_grad_hyb
!------------------------------------------------------------------------------
   subroutine torpol_to_dphhyb(dWlm, Zlm, fL, gL, lcut)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !

      !-- Input variables
      complex(cp), intent(in) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      complex(cp), intent(out) :: gL(n_theta_max,n_m_loc)

      !-- Local variables
      complex(cp) :: Slm(0:l_max), Tlm(0:l_max)
      integer :: i, l_lm, u_lm, m, it
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            gL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         Slm(m:l_max) = ci*m*dWlm(l_lm:u_lm)
         Tlm(m:l_max) = ci*m*Zlm(l_lm:u_lm)
         call shtns_sphtor_to_spat_ml(m/minc, Slm(m:l_max), Tlm(m:l_max), &
              &                       fL(:,i), gL(:,i), lcut)
         do it=1,n_theta_max
            fL(it,i)=fL(it,i)*O_sin_theta_E2(it)
            gL(it,i)=gL(it,i)*O_sin_theta_E2(it)
         end do
      end do
      !$omp end parallel do
      
   end subroutine torpol_to_dphhyb
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_hyb(Qlm, fL, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Qlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variable
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)

      !-- Local variables
      complex(cp) :: dQlm(0:l_max)
      integer :: i, l_lm, u_lm, m
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         if (m>lcut) then
            fL(:,i)=zero
            cycle
         end if
         l_lm = map_dist_st%lm2(m, m)
         u_lm = map_dist_st%lm2(l_max, m)
         
         dQlm(m:l_max) = dLh_loc(l_lm:u_lm) * Qlm(l_lm:u_lm)
         call shtns_SH_to_spat_ml(m/minc, dQlm(m:l_max), fL(:,i), lcut)
      end do
      !$omp end parallel do
      
   end subroutine pol_to_curlr_hyb
!------------------------------------------------------------------------------
   subroutine hyb_to_SH(fL_loc, fLMP_loc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: fL_loc(n_theta_max,n_m_loc)
      integer,     intent(in) :: lcut

      !-- Output variable
      complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
      
      !-- Local variables
      integer :: m, l_lm, u_lm, i
      
      call shtns_load_cfg(1) ! l_max + 1
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            fLMP_loc(l_lm:u_lm)=zero
            cycle
         end if 
         call shtns_spat_to_sh_ml(m/minc,fL_loc(:,i),fLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do
      
      call shtns_load_cfg(0) ! l_max

   end subroutine hyb_to_SH
!------------------------------------------------------------------------------
   subroutine hyb_to_qst(fL_loc, gL_loc, hL_loc, qLMP_loc, sLMP_loc, tLMP_loc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: fL_loc(n_theta_max,n_m_loc)
      complex(cp), intent(in) :: gL_loc(n_theta_max,n_m_loc)
      complex(cp), intent(in) :: hL_loc(n_theta_max,n_m_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: qLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: sLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: tLMP_loc(n_lmP_loc)
      
      !-- Local variables
      integer :: i, l_lm, u_lm, m

      call shtns_load_cfg(1)
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            qLMP_loc(l_lm:u_lm)=zero
            sLMP_loc(l_lm:u_lm)=zero
            tLMP_loc(l_lm:u_lm)=zero
            cycle
         end if
         call shtns_spat_to_qst_ml(m/minc,fL_loc(:,i),gL_loc(:,i),hL_loc(:,i),&
              &                    qLMP_loc(l_lm:u_lm),sLMP_loc(l_lm:u_lm),   &
              &                    tLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do
      
      call shtns_load_cfg(0)

   end subroutine hyb_to_qst
!------------------------------------------------------------------------------
   subroutine hyb_to_sphertor(fL_loc, gL_loc, fLMP_loc, gLMP_loc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: fL_loc(n_theta_max,n_m_loc)
      complex(cp), intent(in) :: gL_loc(n_theta_max,n_m_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: gLMP_loc(n_lmP_loc)
      
      !-- Local variables
      integer :: i, l_lm, u_lm, m

      call shtns_load_cfg(1)
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            fLMP_loc(l_lm:u_lm)=zero
            gLMP_loc(l_lm:u_lm)=zero
            cycle
         end if
         call shtns_spat_to_sphtor_ml(m/minc,fL_loc(:,i),gL_loc(:,i),&
              &                       fLMP_loc(l_lm:u_lm),           &
              &                       gLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do
      
      call shtns_load_cfg(0)

   end subroutine hyb_to_sphertor
!------------------------------------------------------------------------------
end module shtns
