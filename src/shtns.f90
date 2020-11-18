module sht
   !
   ! This module contains is a wrapper of the SHTns routines used in MagIC
   !

   use iso_c_binding
   use iso_fortran_env, only: output_unit
   use precision_mod, only: cp, MPI_DEF_REAL
   use constants, only: ci, one, zero
   use truncation, only: m_max, l_max, n_theta_max, n_phi_max, &
       &                 minc, lm_max, lmP_max, n_lm_loc, n_theta_loc, &
       &                 n_m_loc, n_lmP_loc, n_m_max, dist_m,          &
       &                 nThetaStart, nThetaStop, coord_m, dist_theta, &
       &                 nlat_padded
   use horizontal_data, only: dLh_loc, dLh, O_sin_theta_E2, O_sin_theta
   use parallel_mod
   use LMmapping, only: map_dist_st, map_glbl_st
   use mpi_thetap_mod
   
   use communications, only: theta_transp
   use mpi_transpose_theta !@>TODO replaces mpi_thetap_mod

   implicit none

   include "shtns.f03"

   private

   !@>TODO: the following functions do not have a "_dist" version; the existing
   ! implementation is _loc. Is the _dist version needed?
   public :: initialize_sht, torpol_to_spat_IC, torpol_to_curl_spat_IC, toraxi_to_spat, &
   &         sphtor_to_spat, finalize_sht
   public :: torpol_to_spat_loc, scal_to_spat_loc

   !@>TODO: the following functions do not have a "_loc" version; the existing
   ! implementation is _dist. Is the _loc version needed?
   public :: scal_to_hyb, scal_to_grad_hyb, torpol_to_hyb, torpol_to_curl_hyb, &
   &         pol_to_grad_hyb, torpol_to_dphhyb, pol_to_curlr_hyb, hyb_to_SH,   &
   &         hyb_to_qst, hyb_to_sphertor

   type(c_ptr) :: sht_l, sht_lP
   
   interface
      subroutine scal_to_spat_if(Slm, fieldc, lcut)
         import
         complex(cp), intent(inout) :: Slm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp),   intent(out) :: fieldc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine scal_to_grad_spat_if(Slm, gradtc, gradpc, lcut)
         import
         complex(cp), intent(inout) :: Slm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp),   intent(out) :: gradtc(n_theta_loc,n_phi_max)
         real(cp),   intent(out) :: gradpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine spat_to_SH_if(f_loc, fLMP_loc, lcut)
         import
         real(cp),    intent(inout) :: f_loc(n_theta_loc,n_phi_max)
         integer,     intent(in)    :: lcut
         complex(cp), intent(out)   :: fLMP_loc(n_lmP_loc)
      end subroutine
      
      subroutine spat_to_qst_if(f_loc, g_loc, h_loc, qLMP_loc, sLMP_loc, tLMP_loc, lcut)
         import
         real(cp), intent(inout) :: f_loc(n_theta_loc,n_phi_max)
         real(cp), intent(inout) :: g_loc(n_theta_loc,n_phi_max)
         real(cp), intent(inout) :: h_loc(n_theta_loc,n_phi_max)
         integer,  intent(in)    :: lcut
         complex(cp), intent(out) :: qLMP_loc(n_lmP_loc)
         complex(cp), intent(out) :: sLMP_loc(n_lmP_loc)
         complex(cp), intent(out) :: tLMP_loc(n_lmP_loc)
      end subroutine
      
      subroutine spat_to_sphertor_if(f_loc, g_loc, fLMP_loc, gLMP_loc, lcut)
         import
         real(cp), intent(inout) :: f_loc(n_theta_loc,n_phi_max)
         real(cp), intent(inout) :: g_loc(n_theta_loc,n_phi_max)
         integer,  intent(in)    :: lcut
         complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
         complex(cp), intent(out) :: gLMP_loc(n_lmP_loc)
      end subroutine
      
      subroutine torpol_to_spat_if(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)
         import
         complex(cp), intent(in) :: Wlm(n_lm_loc)
         complex(cp), intent(inout) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp), intent(out) :: vrc(n_theta_loc,n_phi_max)
         real(cp), intent(out) :: vtc(n_theta_loc,n_phi_max)
         real(cp), intent(out) :: vpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine torpol_to_curl_spat_if(or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, &
              &                        cvpc, lcut)
         import
         complex(cp), intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc)
         complex(cp), intent(in) :: Jlm(n_lm_loc)
         complex(cp), intent(inout) :: dJlm(n_lm_loc)
         real(cp),    intent(in) :: or2
         integer,     intent(in) :: lcut
         real(cp), intent(out) :: cvrc(n_theta_loc,n_phi_max)
         real(cp), intent(out) :: cvtc(n_theta_loc,n_phi_max)
         real(cp), intent(out) :: cvpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine pol_to_curlr_spat_if(Qlm, cvrc, lcut)
         import
         complex(cp), intent(in) :: Qlm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp),    intent(out) :: cvrc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine pol_to_grad_spat_if(Slm, gradtc, gradpc, lcut)
         import
         complex(cp), intent(in) :: Slm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp), intent(out) :: gradtc(n_theta_loc,n_phi_max)
         real(cp), intent(out) :: gradpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine torpol_to_dphspat_if(dWlm, Zlm, dvtdp, dvpdp, lcut)
         import
         complex(cp), intent(in) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp), intent(out) :: dvtdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
         real(cp), intent(out) :: dvpdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
      end subroutine
      
      subroutine spat_to_SH_axi_if(f, fLM)
         import
         real(cp), intent(in)  :: f(nThetaStart:nThetaStop)
         real(cp), intent(out) :: fLM(:)
      end subroutine 

   end interface
   
   public :: scal_to_spat, scal_to_grad_spat, scal_to_SH, spat_to_qst, &
   & spat_to_sphertor, torpol_to_spat, torpol_to_curl_spat, pol_to_curlr_spat, &
   & pol_to_grad_spat, torpol_to_dphspat, spat_to_SH_axi
   
   procedure (scal_to_spat_if), pointer :: scal_to_spat => null ()
   procedure (scal_to_grad_spat_if), pointer :: scal_to_grad_spat => null ()
   procedure (spat_to_SH_if), pointer :: scal_to_SH => null ()
   procedure (spat_to_qst_if), pointer :: spat_to_qst => null ()
   procedure (spat_to_sphertor_if), pointer :: spat_to_sphertor => null ()
   procedure (torpol_to_spat_if), pointer :: torpol_to_spat => null ()
   procedure (torpol_to_curl_spat_if), pointer :: torpol_to_curl_spat => null ()
   procedure (pol_to_curlr_spat_if), pointer :: pol_to_curlr_spat => null ()
   procedure (pol_to_grad_spat_if), pointer :: pol_to_grad_spat => null ()
   procedure (torpol_to_dphspat_if), pointer :: torpol_to_dphspat => null ()
   procedure (spat_to_SH_axi_if), pointer :: spat_to_SH_axi => null ()
   
contains

   subroutine initialize_sht()

      integer :: norm, layout, nthreads
      real(cp) :: eps_polar
      type(shtns_info), pointer :: sht_info

      if ( l_master_rank ) then
         write(output_unit,*) ''
         call shtns_verbose(1)
      end if

      nthreads =  shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

#ifdef SHT_PADDING
      !layout = SHT_QUICK_INIT + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
      if ( n_ranks_theta == 1 ) then
         layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_ALLOW_PADDING
      else
         layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS
      end if
#else
      layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS
#endif
      eps_polar = 1.e-10_cp

      sht_l = shtns_create(l_max, m_max/minc, minc, norm)
      call shtns_set_grid(sht_l, layout, eps_polar, n_theta_max, n_phi_max)

      call c_f_pointer(cptr=sht_l, fptr=sht_info)
#ifdef SHT_PADDING
      if ( n_ranks_theta == 1 ) then
         nlat_padded = sht_info%nlat_padded
         nThetaStart = 1
         nThetaStop  = nlat_padded
         n_theta_loc = nlat_padded
      else
         nlat_padded = n_theta_max
      end if
      if ( nlat_padded /= n_theta_max .and. rank == 0 ) then
         write(output_unit,*) '! SHTns uses theta padding with nlat_padded=', nlat_padded
      end if
#else
      nlat_padded = n_theta_max
#endif

      if ( l_master_rank ) then
         call shtns_verbose(0)
         write(output_unit,*) ''
      end if

      sht_lP = shtns_create(l_max+1, m_max/minc, minc, norm)
      call shtns_set_grid(sht_lP, layout, eps_polar, n_theta_max, n_phi_max)
      
      ! Assigns all pointers to the correct functions
      if ( n_ranks_theta>1 ) then
         scal_to_spat => scal_to_spat_dist
         scal_to_grad_spat => scal_to_grad_spat_dist
         scal_to_SH => spat_to_SH_dist
         spat_to_qst => spat_to_qst_dist
         spat_to_sphertor => spat_to_sphertor_dist
         torpol_to_spat => torpol_to_spat_dist
         torpol_to_curl_spat => torpol_to_curl_spat_dist
         pol_to_curlr_spat => pol_to_curlr_spat_dist
         pol_to_grad_spat => pol_to_grad_spat_dist
         torpol_to_dphspat => torpol_to_dphspat_dist
         spat_to_SH_axi => spat_to_SH_axi_dist
      else
         scal_to_spat => scal_to_spat_loc
         scal_to_grad_spat => scal_to_grad_spat_loc
         scal_to_SH => spat_to_SH_loc
         spat_to_qst => spat_to_qst_loc
         spat_to_sphertor => spat_to_sphertor_loc
         torpol_to_spat => torpol_to_spat_loc
         torpol_to_curl_spat => torpol_to_curl_spat_loc
         pol_to_curlr_spat => pol_to_curlr_spat_loc
         pol_to_grad_spat => pol_to_grad_spat_loc
         torpol_to_dphspat => torpol_to_dphspat_loc
         spat_to_SH_axi => spat_to_SH_axi_loc
      end if

   end subroutine initialize_sht
!------------------------------------------------------------------------------
   subroutine finalize_sht

      call shtns_unset_grid(sht_l)
      call shtns_destroy(sht_l)
      call shtns_unset_grid(sht_lP)
      call shtns_destroy(sht_lP)

   end subroutine finalize_sht
!------------------------------------------------------------------------------
   subroutine scal_to_spat_loc(Slm, fieldc, lcut)
      ! transform a spherical harmonic field into grid space

      !-- Input variables
      complex(cp), intent(inout) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: fieldc(nlat_padded,n_phi_max)

      call SH_to_spat_l(sht_l, Slm, fieldc, lcut)

   end subroutine scal_to_spat_loc
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat_loc(Slm, gradtc, gradpc, lcut)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variables
      complex(cp), intent(inout) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: gradpc(nlat_padded,n_phi_max)

      call SHsph_to_spat_l(sht_l, Slm, gradtc, gradpc, lcut)

   end subroutine scal_to_grad_spat_loc
!------------------------------------------------------------------------------
   subroutine pol_to_grad_spat_loc(Slm, gradtc, gradpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: gradpc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = map_glbl_st%lm2l(lm)
         if ( l <= lcut ) then
            Qlm(lm) = dLh(lm) * Slm(lm)
         else
            Qlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call SHsph_to_spat_l(sht_l, Qlm, gradtc, gradpc, lcut)

   end subroutine pol_to_grad_spat_loc
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_loc(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(lm_max)
      complex(cp), intent(inout) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: vtc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: vpc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = map_glbl_st%lm2l(lm)
         if ( l <= lcut ) then
            Qlm(lm) = dLh(lm) * Wlm(lm)
         else
            Qlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call SHqst_to_spat_l(sht_l, Qlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

   end subroutine torpol_to_spat_loc
!------------------------------------------------------------------------------
   subroutine sphtor_to_spat(dWlm, Zlm, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(inout) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vtc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: vpc(nlat_padded,n_phi_max)

      call SHsphtor_to_spat_l(sht_l, dWlm, Zlm, vtc, vpc, lcut)

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
      real(cp), intent(out) :: cbr(nlat_padded,n_phi_max)
      real(cp), intent(out) :: cbt(nlat_padded,n_phi_max)
      real(cp), intent(out) :: cbp(nlat_padded,n_phi_max)

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

      call SHqst_to_spat(sht_l, Qlm, Slm, Tlm, cbr, cbt, cbp)

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
      real(cp), intent(out) :: Br(nlat_padded,n_phi_max)
      real(cp), intent(out) :: Bt(nlat_padded,n_phi_max)
      real(cp), intent(out) :: Bp(nlat_padded,n_phi_max)

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

      call SHqst_to_spat(sht_l, Qlm, Slm, Tlm, Br, Bt, Bp)

   end subroutine torpol_to_spat_IC
!------------------------------------------------------------------------------
   subroutine torpol_to_dphspat_loc(dWlm, Zlm, dvtdp, dvpdp, lcut)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !

      !-- Input variables
      complex(cp), intent(in) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: dvtdp(nlat_padded,n_phi_max)
      real(cp), intent(out) :: dvpdp(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Slm(lm_max), Tlm(lm_max)
      integer :: lm, ip, l
      real(cp) :: m

      !$omp parallel do default(shared) private(lm, m, l)
      do lm = 1, lm_max
         l = map_glbl_st%lm2l(lm)
         m = map_glbl_st%lm2m(lm)
         if ( l <= lcut ) then
            Slm(lm) = ci*m*dWlm(lm)
            Tlm(lm) = ci*m*Zlm(lm)
         else
            Slm(lm) = 0.0_cp
            Tlm(lm) = 0.0_cp
         end if
      end do
      !$omp end parallel do

      call SHsphtor_to_spat_l(sht_l, Slm, Tlm, dvtdp, dvpdp, lcut)

      !$omp parallel do default(shared) private(ip)
      do ip=1, n_phi_max
         dvtdp(:, ip) = dvtdp(:, ip) * O_sin_theta_E2(:)
         dvpdp(:, ip) = dvpdp(:, ip) * O_sin_theta_E2(:)
      end do
      !$omp end parallel do

   end subroutine torpol_to_dphspat_loc
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_loc(Qlm, cvrc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Qlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: cvrc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: dQlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = map_glbl_st%lm2l(lm)
         if ( l <= lcut ) then
            dQlm(lm) = dLh(lm) * Qlm(lm)
         else
            dQlm(lm) = zero
         end if
      end do
      !$omp end parallel do

      call SH_to_spat_l(sht_l, dQlm, cvrc, lcut)

   end subroutine pol_to_curlr_spat_loc
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_loc(or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, cvpc, &
              &                   lcut)

      !-- Input variables
      complex(cp), intent(in) :: Blm(lm_max), ddBlm(lm_max), Jlm(lm_max)
      complex(cp), intent(inout) :: dJlm(lm_max)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: cvrc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: cvtc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: cvpc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm, l

      !$omp parallel do default(shared) private(lm, l)
      do lm = 1, lm_max
         l = map_glbl_st%lm2l(lm)
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

      call SHqst_to_spat_l(sht_l, Qlm, dJlm, Tlm, cvrc, cvtc, cvpc, lcut)

   end subroutine torpol_to_curl_spat_loc
!------------------------------------------------------------------------------
   subroutine spat_to_SH_loc(f, fLM, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f(nlat_padded,n_phi_max)
      integer,  intent(in) :: lcut

      !-- Output variable
      complex(cp), intent(out) :: fLM(lmP_max)

      call spat_to_SH_l(sht_lP, f, fLM, lcut+1)

   end subroutine spat_to_SH_loc
!------------------------------------------------------------------------------
   subroutine spat_to_qst_loc(f, g, h, qLM, sLM, tLM, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f(nlat_padded,n_phi_max)
      real(cp), intent(inout) :: g(nlat_padded,n_phi_max)
      real(cp), intent(inout) :: h(nlat_padded,n_phi_max)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: qLM(lmP_max)
      complex(cp), intent(out) :: sLM(lmP_max)
      complex(cp), intent(out) :: tLM(lmP_max)

      call spat_to_SHqst_l(sht_lP, f, g, h, qLM, sLM, tLM, lcut+1)

   end subroutine spat_to_qst_loc
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor_loc(f, g, fLM, gLM, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f(nlat_padded,n_phi_max)
      real(cp), intent(inout) :: g(nlat_padded,n_phi_max)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLM(lmP_max)
      complex(cp), intent(out) :: gLM(lmP_max)

      call spat_to_SHsphtor_l(sht_lP, f, g, fLM, gLM, lcut+1)

   end subroutine spat_to_sphertor_loc
!------------------------------------------------------------------------------
   subroutine toraxi_to_spat(fl_ax, ft, fp)

      !-- Input field
      complex(cp), intent(inout) :: fl_ax(l_max+1) !-- Axi-sym toroidal

      !-- Output fields on grid
      real(cp), intent(out) :: ft(:)
      real(cp), intent(out) :: fp(:)

      !-- Local arrays
      complex(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)

      call SHtor_to_spat_ml(sht_l, 0, fl_ax, tmpt, tmpp, l_max)
      ft(:)=real(tmpt(:))
      fp(:)=real(tmpp(:))

   end subroutine toraxi_to_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH_axi_loc(f, fLM)

      real(cp), intent(in) :: f(:)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp(nlat_padded)
      complex(cp) :: tmpLM(size(fLM))

      tmp(:)=cmplx(f(:),0.0_cp,kind=cp)
      if ( size(fLM) == l_max+2 ) then
         call spat_to_SH_ml(sht_lP, 0, tmp, tmpLM, l_max+1)
      else if ( size(fLM) == l_max+1 ) then
         call spat_to_SH_ml(sht_l, 0, tmp, tmpLM, l_max)
      end if
      fLM(:)=real(tmpLM(:))

   end subroutine spat_to_SH_axi_loc
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
      real(cp),     intent(out)   :: fieldc_loc(n_theta_loc,n_phi_max)
      integer,      intent(in)    :: lcut

      !-- Local variables
      complex(cp) :: fL_loc(n_theta_max,n_m_loc)

      call scal_to_hyb(fLM_loc, fL_loc, lcut)
!       call transform_m2phi(fL_loc, fieldc_loc)
      call transform_m2phi_new(theta_transp, fL_loc, fieldc_loc) 
      
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
      complex(cp), intent(inout) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_theta_loc,n_phi_max)
      real(cp), intent(out) :: gradpc(n_theta_loc,n_phi_max)
      
      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      call scal_to_grad_hyb(Slm, fL, gL, lcut)
      call transform_m2phi(fL, gradtc)
      call transform_m2phi(gL, gradpc)

   end subroutine scal_to_grad_spat_dist
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_dist(Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(n_lm_loc)
      complex(cp), intent(inout) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vrc(n_theta_loc,n_phi_max)
      real(cp), intent(out) :: vtc(n_theta_loc,n_phi_max)
      real(cp), intent(out) :: vpc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      complex(cp) :: hL(n_theta_max,n_m_loc)
      
      call torpol_to_hyb(Wlm, dWlm, Zlm, fL, gL, hL, lcut)
      call transform_m2phi(fL, vrc)
      call transform_m2phi(gL, vtc)
      call transform_m2phi(hL, vpc)

   end subroutine torpol_to_spat_dist
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_dist(or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, &
              &                        cvpc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc), Jlm(n_lm_loc)
      complex(cp), intent(inout) :: dJlm(n_lm_loc)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: cvrc(n_theta_loc,n_phi_max)
      real(cp), intent(out) :: cvtc(n_theta_loc,n_phi_max)
      real(cp), intent(out) :: cvpc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max), Tlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      complex(cp) :: hL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      call torpol_to_curl_hyb(or2, Blm, ddBlm, Jlm, dJlm, fL, gL, hL, lcut)
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
      real(cp), intent(out) :: gradtc(n_theta_loc,n_phi_max)
      real(cp), intent(out) :: gradpc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      call pol_to_grad_hyb(Slm, fL, gL, lcut)
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
      real(cp), intent(out) :: dvtdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
      real(cp), intent(out) :: dvpdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      call torpol_to_dphhyb(dWlm, Zlm, fL, gL, lcut)      
      call transform_m2phi(fL, dvtdp)
      call transform_m2phi(gL, dvpdp)
      
   end subroutine torpol_to_dphspat_dist
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_dist(Qlm, cvrc, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Qlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), intent(out) :: cvrc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      
      call pol_to_curlr_hyb(Qlm, fL, lcut)      
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
      real(cp),    intent(inout) :: f_loc(n_theta_loc,n_phi_max)
      integer,     intent(in)    :: lcut

      !-- Output variable
      complex(cp), intent(out)   :: fLMP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      integer :: m, l_lm, u_lm, i
      
!       call transform_phi2m_new(theta_transp, f_loc, fL_loc)
      call transform_phi2m(f_loc, fL_loc)
      call hyb_to_SH(fL_loc, fLMP_loc, lcut)
      
   end subroutine spat_to_SH_dist
!------------------------------------------------------------------------------
   subroutine spat_to_qst_dist(f_loc, g_loc, h_loc, qLMP_loc, sLMP_loc, tLMP_loc, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f_loc(n_theta_loc,n_phi_max)
      real(cp), intent(inout) :: g_loc(n_theta_loc,n_phi_max)
      real(cp), intent(inout) :: h_loc(n_theta_loc,n_phi_max)
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

      !>@TODO Vectorial FFT and transpose (f,g,h at once)
      call transform_phi2m(f_loc, fL_loc)
      call transform_phi2m(g_loc, gL_loc)
      call transform_phi2m(h_loc, hL_loc)
      call hyb_to_qst(fL_loc, gL_loc, hL_loc, qLMP_loc, sLMP_loc, tLMP_loc, lcut)

   end subroutine spat_to_qst_dist
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor_dist(f_loc, g_loc, fLMP_loc, gLMP_loc, lcut)

      !-- Input variables
      real(cp), intent(inout) :: f_loc(n_theta_loc,n_phi_max)
      real(cp), intent(inout) :: g_loc(n_theta_loc,n_phi_max)
      integer,  intent(in)    :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: gLMP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  gL_loc(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m

      !>@TODO Vectorial FFT and transpose (f,g,h at once)
      call transform_phi2m(f_loc, fL_loc)
      call transform_phi2m(g_loc, gL_loc)
      call hyb_to_sphertor(fL_loc, gL_loc, fLMP_loc, gLMP_loc, lcut)

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
      if ( size(fLM) == l_max+2 ) then
         call spat_to_SH_ml(sht_lP, 0, tmp_c, tmpLM, size(fLM)-1)
      else if ( size(fLM) == l_max+1 ) then
         call spat_to_SH_ml(sht_l, 0, tmp_c, tmpLM, size(fLM)-1)
      end if
       
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
      complex(cp), intent(inout) :: fLM_loc(n_lm_loc)
      
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

        call SH_to_spat_ml(sht_l, m/minc, fLM_loc(l_lm:u_lm), field_hyb(:,i),lcut)

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
      complex(cp), intent(inout) :: Slm(n_lm_loc)
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
         call SHsph_to_spat_ml(sht_l, m/minc, Slm(l_lm:u_lm), gradtL(:,i), gradpL(:,i), &
              &                lcut)
      end do
      !$omp end parallel do
      
   end subroutine scal_to_grad_hyb
!------------------------------------------------------------------------------
   subroutine torpol_to_hyb(Wlm, dWlm, Zlm, fL, gL, hL, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Wlm(n_lm_loc)
      complex(cp), intent(inout) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
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
         call SHqst_to_spat_ml(sht_l, m/minc, Qlm(m:l_max), dWlm(l_lm:u_lm), &
              &                Zlm(l_lm:u_lm), fL(:,i), gL(:,i), hL(:,i), lcut)
      end do
      !$omp end parallel do
      
   end subroutine torpol_to_hyb
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_hyb(or2, Blm, ddBlm, Jlm, dJlm, fL, gL, hL, lcut)

      !-- Input variables
      complex(cp), intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc), Jlm(n_lm_loc)
      complex(cp), intent(inout) :: dJlm(n_lm_loc)
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
         call SHqst_to_spat_ml(sht_l, m/minc, Qlm(m:l_max), dJlm(l_lm:u_lm), &
              &                Tlm(m:l_max), fL(:,i), gL(:,i), hL(:,i), lcut)
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
         call SHsph_to_spat_ml(sht_l, m/minc, Qlm(m:l_max), fL(:,i), gL(:,i), lcut)
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
         call SHsphtor_to_spat_ml(sht_l, m/minc, Slm(m:l_max), Tlm(m:l_max), &
              &                   fL(:,i), gL(:,i), lcut)
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
         call SH_to_spat_ml(sht_l, m/minc, dQlm(m:l_max), fL(:,i), lcut)
      end do
      !$omp end parallel do
      
   end subroutine pol_to_curlr_hyb
!------------------------------------------------------------------------------
   subroutine hyb_to_SH(fL_loc, fLMP_loc, lcut)

      !-- Input variables
      complex(cp), intent(inout) :: fL_loc(n_theta_max,n_m_loc)
      integer,     intent(in) :: lcut

      !-- Output variable
      complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
      
      !-- Local variables
      integer :: m, l_lm, u_lm, i
      
      !$omp parallel do default(shared) private(i,m,l_lm,u_lm)
      do i = 1, n_m_loc
         m = dist_m(coord_m, i)
         l_lm = map_dist_st%lmP2(m, m)
         u_lm = map_dist_st%lmP2(l_max+1, m)
         if ( m>min(m_max,lcut+1) ) then
            fLMP_loc(l_lm:u_lm)=zero
            cycle
         end if 
         call spat_to_SH_ml(sht_lP, m/minc,fL_loc(:,i),fLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do

   end subroutine hyb_to_SH
!------------------------------------------------------------------------------
   subroutine hyb_to_qst(fL_loc, gL_loc, hL_loc, qLMP_loc, sLMP_loc, tLMP_loc, lcut)

      !-- Input variables
      complex(cp), intent(inout) :: fL_loc(n_theta_max,n_m_loc)
      complex(cp), intent(inout) :: gL_loc(n_theta_max,n_m_loc)
      complex(cp), intent(inout) :: hL_loc(n_theta_max,n_m_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: qLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: sLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: tLMP_loc(n_lmP_loc)
      
      !-- Local variables
      integer :: i, l_lm, u_lm, m

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
         call spat_to_SHqst_ml(sht_lP, m/minc,fL_loc(:,i),gL_loc(:,i),hL_loc(:,i),&
              &                qLMP_loc(l_lm:u_lm),sLMP_loc(l_lm:u_lm),           &
              &                tLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do

   end subroutine hyb_to_qst
!------------------------------------------------------------------------------
   subroutine hyb_to_sphertor(fL_loc, gL_loc, fLMP_loc, gLMP_loc, lcut)

      !-- Input variables
      complex(cp), intent(inout) :: fL_loc(n_theta_max,n_m_loc)
      complex(cp), intent(inout) :: gL_loc(n_theta_max,n_m_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      complex(cp), intent(out) :: fLMP_loc(n_lmP_loc)
      complex(cp), intent(out) :: gLMP_loc(n_lmP_loc)
      
      !-- Local variables
      integer :: i, l_lm, u_lm, m

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
         call spat_to_SHsphtor_ml(sht_lP,m/minc,fL_loc(:,i),gL_loc(:,i),&
              &                   fLMP_loc(l_lm:u_lm),                  &
              &                   gLMP_loc(l_lm:u_lm),lcut+1)
      end do
      !$omp end parallel do

   end subroutine hyb_to_sphertor
!------------------------------------------------------------------------------
end module sht
