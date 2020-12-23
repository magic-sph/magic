#include "perflib_preproc.cpp"
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
   use num_param, only: sht_buffer_size
   use fft, only: fft_many, ifft_many

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
   
   
!    public :: scal_to_spat, scal_to_grad_spat, scal_to_SH, spat_to_qst, &
!    & spat_to_sphertor, torpol_to_spat, torpol_to_curl_spat, pol_to_curlr_spat, &
!    & pol_to_grad_spat, torpol_to_dphspat, spat_to_SH_axi

   type(c_ptr) :: sht_l, sht_lP
   integer, allocatable :: test_ptr(:)
   
   ! --------------------------------------------------------------------------------------
   type, abstract, public :: type_shtns_if
      !
      !   Author: Rafael Lago (MPCDF) December 2020
      ! 
   contains
      procedure :: commit_backward => commit_dummy
      procedure :: commit_forward  => commit_dummy
      procedure(scal_to_spat_if), deferred        :: scal_to_spat      
      procedure(scal_to_grad_spat_if), deferred   :: scal_to_grad_spat 
      procedure(torpol_to_spat_if), deferred      :: torpol_to_spat    
      procedure(torpol_to_curl_spat_if), deferred :: torpol_to_curl_spat
      procedure(pol_to_grad_spat_if), deferred    :: pol_to_grad_spat  
      procedure(pol_to_curlr_spat_if), deferred   :: pol_to_curlr_spat 
      procedure(torpol_to_dphspat_if), deferred   :: torpol_to_dphspat 
      procedure(scal_to_SH_if), deferred          :: scal_to_SH        
      procedure(spat_to_qst_if), deferred         :: spat_to_qst       
      procedure(spat_to_sphertor_if), deferred    :: spat_to_sphertor
      procedure(spat_to_SH_axi_if), deferred      :: spat_to_SH_axi
   end type type_shtns_if
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_shtns_if) :: type_shtns_buff
      !
      !   Author: Rafael Lago (MPCDF) December 2020
      ! 
      type(real_pointer_wrapper),  allocatable :: spat_ptr(:)
      type(cmplx_pointer_wrapper), allocatable :: sh_ptr(:)
      integer, allocatable :: operations(:), lcut(:)
      real(cp), allocatable :: or2(:)
      integer :: max_buff, n_sh, n_spat
   contains
      procedure :: initialize  => initialize_buff
      procedure :: commit_backward  => commit_backward_buff
      procedure :: commit_forward   => commit_forward_buff
      procedure :: scal_to_spat        => scal_to_spat_buff
      procedure :: scal_to_grad_spat   => scal_to_grad_spat_buff
      procedure :: torpol_to_spat      => torpol_to_spat_buff
      procedure :: torpol_to_curl_spat => torpol_to_curl_spat_buff
      procedure :: pol_to_grad_spat    => pol_to_grad_spat_buff
      procedure :: pol_to_curlr_spat   => pol_to_curlr_spat_buff
      procedure :: torpol_to_dphspat   => torpol_to_dphspat_buff
      procedure :: scal_to_SH          => scal_to_SH_buff
      procedure :: spat_to_qst         => spat_to_qst_buff
      procedure :: spat_to_sphertor    => spat_to_sphertor_buff
      procedure :: spat_to_SH_axi      => spat_to_SH_axi_buff
      final :: finalize_buff
   end type type_shtns_buff
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_shtns_if) :: type_shtns_loc
      !
      !   Author: Rafael Lago (MPCDF) December 2020
      ! 
   contains
      procedure :: scal_to_spat        => scal_to_spat_loc
      procedure :: scal_to_grad_spat   => scal_to_grad_spat_loc
      procedure :: torpol_to_spat      => torpol_to_spat_loc
      procedure :: torpol_to_curl_spat => torpol_to_curl_spat_loc
      procedure :: pol_to_grad_spat    => pol_to_grad_spat_loc
      procedure :: pol_to_curlr_spat   => pol_to_curlr_spat_loc
      procedure :: torpol_to_dphspat   => torpol_to_dphspat_loc
      procedure :: scal_to_SH          => scal_to_SH_loc
      procedure :: spat_to_qst         => spat_to_qst_loc
      procedure :: spat_to_sphertor    => spat_to_sphertor_loc
      procedure :: spat_to_SH_axi      => spat_to_SH_axi_loc
   end type type_shtns_loc
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_shtns_if) :: type_shtns_dist
      !
      !   Author: Rafael Lago (MPCDF) December 2020
      ! 
   contains
      procedure :: scal_to_spat        => scal_to_spat_dist
      procedure :: scal_to_grad_spat   => scal_to_grad_spat_dist
      procedure :: torpol_to_spat      => torpol_to_spat_dist
      procedure :: torpol_to_curl_spat => torpol_to_curl_spat_dist
      procedure :: pol_to_grad_spat    => pol_to_grad_spat_dist
      procedure :: pol_to_curlr_spat   => pol_to_curlr_spat_dist
      procedure :: torpol_to_dphspat   => torpol_to_dphspat_dist
      procedure :: scal_to_SH          => scal_to_SH_dist
      procedure :: spat_to_qst         => spat_to_qst_dist
      procedure :: spat_to_sphertor    => spat_to_sphertor_dist
      procedure :: spat_to_SH_axi      => spat_to_SH_axi_dist
   end type type_shtns_dist
   
   !   Workarounds the lack of array of pointers in fortran
   ! --------------------------------------------------------------------------------------
   type cmplx_pointer_wrapper
      complex(cp), contiguous, pointer :: p(:)
   end type cmplx_pointer_wrapper
   type real_pointer_wrapper
      real(cp), contiguous, pointer :: p(:,:)
   end type real_pointer_wrapper

   integer, parameter :: OP_NONE             = 0
   integer, parameter :: OP_SPAT2SH          = 1
   integer, parameter :: OP_SPAT2QST         = 2
   integer, parameter :: OP_SPAT2SPHERTOR    = 3
   integer, parameter :: OP_SPAT2SH_AXI      = 4
   integer, parameter :: OP_SCAL2SPAT        = -1
   integer, parameter :: OP_SCAL2GRAD_SPAT   = -2
   integer, parameter :: OP_TORPOL2SPAT      = -3
   integer, parameter :: OP_TORPOL2CURL_SPAT = -4
   integer, parameter :: OP_POL2CURLR_SPAT   = -5
   integer, parameter :: OP_POL2GRAD_SPAT    = -6
   integer, parameter :: OP_TORPOL2DPHSPAT   = -7
   
   interface
      subroutine scal_to_spat_if(this, Slm, fieldc, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(inout) :: Slm(n_lm_loc)
         real(cp),    target, intent(out) :: fieldc(n_theta_loc,n_phi_max)
         integer,             intent(in) :: lcut
      end subroutine
      
      subroutine scal_to_grad_spat_if(this, Slm, gradtc, gradpc, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(inout) :: Slm(n_lm_loc)
         real(cp),   target, intent(out) :: gradtc(n_theta_loc,n_phi_max)
         real(cp),   target, intent(out) :: gradpc(n_theta_loc,n_phi_max)
         integer,            intent(in) :: lcut
      end subroutine
      
      subroutine scal_to_SH_if(this, f, fLMP, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         real(cp),    target, intent(inout) :: f(n_theta_loc,n_phi_max)
         complex(cp), target, intent(out)   :: fLMP(n_lmP_loc)
         integer,     intent(in)    :: lcut
      end subroutine
      
      subroutine spat_to_qst_if(this, f, g, h, qLMP, sLMP, tLMP, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         real(cp), target, intent(inout) :: f(n_theta_loc,n_phi_max)
         real(cp), target, intent(inout) :: g(n_theta_loc,n_phi_max)
         real(cp), target, intent(inout) :: h(n_theta_loc,n_phi_max)
         integer,  intent(in)    :: lcut
         complex(cp), target, intent(out) :: qLMP(n_lmP_loc)
         complex(cp), target, intent(out) :: sLMP(n_lmP_loc)
         complex(cp), target, intent(out) :: tLMP(n_lmP_loc)
      end subroutine
      
      subroutine spat_to_sphertor_if(this, f, g, fLMP, gLMP, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         real(cp), target, intent(inout) :: f(n_theta_loc,n_phi_max)
         real(cp), target, intent(inout) :: g(n_theta_loc,n_phi_max)
         integer,  intent(in)    :: lcut
         complex(cp), target, intent(out) :: fLMP(n_lmP_loc)
         complex(cp), target, intent(out) :: gLMP(n_lmP_loc)
      end subroutine
      
      subroutine torpol_to_spat_if(this, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(in) :: Wlm(n_lm_loc)
         complex(cp), target, intent(inout) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp), target, intent(out) :: vrc(n_theta_loc,n_phi_max)
         real(cp), target, intent(out) :: vtc(n_theta_loc,n_phi_max)
         real(cp), target, intent(out) :: vpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine torpol_to_curl_spat_if(this, or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, &
              &                        cvpc, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc)
         complex(cp), target, intent(in) :: Jlm(n_lm_loc)
         complex(cp), target, intent(inout) :: dJlm(n_lm_loc)
         real(cp),    intent(in) :: or2
         integer,     intent(in) :: lcut
         real(cp), target, intent(out) :: cvrc(n_theta_loc,n_phi_max)
         real(cp), target, intent(out) :: cvtc(n_theta_loc,n_phi_max)
         real(cp), target, intent(out) :: cvpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine pol_to_curlr_spat_if(this, Qlm, cvrc, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(in) :: Qlm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp),    target, intent(out) :: cvrc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine pol_to_grad_spat_if(this, Slm, gradtc, gradpc, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(in) :: Slm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp), target, intent(out) :: gradtc(n_theta_loc,n_phi_max)
         real(cp), target, intent(out) :: gradpc(n_theta_loc,n_phi_max)
      end subroutine
      
      subroutine torpol_to_dphspat_if(this, dWlm, Zlm, dvtdp, dvpdp, lcut)
         import
         class(type_shtns_if), intent(inout) :: this
         complex(cp), target, intent(in) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
         integer,     intent(in) :: lcut
         real(cp), target, intent(out) :: dvtdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
         real(cp), target, intent(out) :: dvpdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
      end subroutine
      
      subroutine spat_to_SH_axi_if(this, f, fLM)
         import
         class(type_shtns_if), intent(inout) :: this
         real(cp), intent(in)  :: f(nThetaStart:nThetaStop)
         real(cp), intent(out) :: fLM(:)
      end subroutine 
      
   end interface
   
!    procedure (scal_to_spat_if), pointer :: scal_to_spat => null ()
!    procedure (scal_to_grad_spat_if), pointer :: scal_to_grad_spat => null ()
!    procedure (spat_to_qst_if), pointer :: spat_to_qst => null ()
!    procedure (spat_to_sphertor_if), pointer :: spat_to_sphertor => null ()
!    procedure (torpol_to_spat_if), pointer :: torpol_to_spat => null ()
!    procedure (torpol_to_curl_spat_if), pointer :: torpol_to_curl_spat => null ()
!    procedure (pol_to_curlr_spat_if), pointer :: pol_to_curlr_spat => null ()
!    procedure (pol_to_grad_spat_if), pointer :: pol_to_grad_spat => null ()
!    procedure (torpol_to_dphspat_if), pointer :: torpol_to_dphspat => null ()
!    procedure (scal_to_SH_if), pointer :: scal_to_SH => null ()
!    procedure (spat_to_SH_axi_if), pointer :: spat_to_SH_axi => null ()
   class(type_shtns_if), allocatable, public :: SHtransf
   type(type_shtns_loc), public :: SHloc  ! used for outputs, for instance
   
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
      
      if (n_ranks_theta==1) then
         allocate(type_shtns_loc :: SHtransf)
         if (rank==0) print *, '! SHT buffer size = None (local version)' 
      else if (sht_buffer_size==1) then
         allocate(type_shtns_dist :: SHtransf)
         if (rank==0) print *, '! SHT buffer size = 1 (distributed version)' 
      else
         allocate(type_shtns_buff :: SHtransf)
         if (rank==0) print *, '! SHT buffer size =',sht_buffer_size, &
         &    ' (buffered distributed version)' 
         select type(SHtransf)
            class is (type_shtns_buff)
               call SHtransf%initialize(sht_buffer_size)
         end select
      end if

   end subroutine initialize_sht
!------------------------------------------------------------------------------
   subroutine finalize_sht

      call shtns_unset_grid(sht_l)
      call shtns_destroy(sht_l)
      call shtns_unset_grid(sht_lP)
      call shtns_destroy(sht_lP)
      
      if (allocated(SHtransf)) deallocate(SHtransf)

   end subroutine finalize_sht
!------------------------------------------------------------------------------
   subroutine scal_to_spat_loc(this, Slm, fieldc, lcut)
      ! transform a spherical harmonic field into grid space

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(inout) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), target, intent(out) :: fieldc(nlat_padded,n_phi_max)

      PERFON("sht_bwd")
      call SH_to_spat_l(sht_l, Slm, fieldc, lcut)
      PERFOFF

   end subroutine scal_to_spat_loc
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat_loc(this, Slm, gradtc, gradpc, lcut)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(inout) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: gradtc(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: gradpc(nlat_padded,n_phi_max)

      PERFON("sht_bwd")
      call SHsph_to_spat_l(sht_l, Slm, gradtc, gradpc, lcut)
      PERFOFF

   end subroutine scal_to_grad_spat_loc
!------------------------------------------------------------------------------
   subroutine pol_to_grad_spat_loc(this, Slm, gradtc, gradpc, lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(in) :: Slm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: gradtc(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: gradpc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l

      PERFON("sht_bwd")
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
      PERFOFF

   end subroutine pol_to_grad_spat_loc
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_loc(this, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(in) :: Wlm(lm_max)
      complex(cp), target, intent(inout) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: vrc(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: vtc(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: vpc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm, l
      
      PERFON("sht_bwd")
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
      PERFOFF

   end subroutine torpol_to_spat_loc
!------------------------------------------------------------------------------
   subroutine sphtor_to_spat(dWlm, Zlm, vtc, vpc, lcut)

      !-- Input variables
      complex(cp), intent(inout) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), intent(out) :: vtc(nlat_padded,n_phi_max)
      real(cp), intent(out) :: vpc(nlat_padded,n_phi_max)

      PERFON("sht_bwd")
      call SHsphtor_to_spat_l(sht_l, dWlm, Zlm, vtc, vpc, lcut)
      PERFOFF

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

      PERFON("sht_bwd")
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
      PERFOFF

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

      PERFON("sht_bwd")
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
      PERFOFF

   end subroutine torpol_to_spat_IC
!------------------------------------------------------------------------------
   subroutine torpol_to_dphspat_loc(this, dWlm, Zlm, dvtdp, dvpdp, lcut)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(in) :: dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: dvtdp(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: dvpdp(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Slm(lm_max), Tlm(lm_max)
      integer :: lm, ip, l
      real(cp) :: m

      PERFON("sht_bwd")
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
      PERFOFF

   end subroutine torpol_to_dphspat_loc
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_loc(this, Qlm, cvrc, lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(in) :: Qlm(lm_max)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), target, intent(out) :: cvrc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: dQlm(lm_max)
      integer :: lm, l

      PERFON("sht_bwd")
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
      PERFOFF

   end subroutine pol_to_curlr_spat_loc
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_loc(this, or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, cvpc, &
              &                   lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      complex(cp), target, intent(in) :: Blm(lm_max), ddBlm(lm_max), Jlm(lm_max)
      complex(cp), target, intent(inout) :: dJlm(lm_max)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: cvrc(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: cvtc(nlat_padded,n_phi_max)
      real(cp), target, intent(out) :: cvpc(nlat_padded,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm, l

      PERFON("sht_bwd")
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
      PERFOFF

   end subroutine torpol_to_curl_spat_loc
!------------------------------------------------------------------------------
   subroutine scal_to_SH_loc(this, f, fLMP, lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      real(cp), target, intent(inout) :: f(nlat_padded,n_phi_max)
      integer,  intent(in) :: lcut

      !-- Output variable
      complex(cp), target, intent(out) :: fLMP(lmP_max)

      PERFON("sht_fwd")
      call spat_to_SH_l(sht_lP, f, fLMP, lcut+1)
      PERFOFF

   end subroutine scal_to_SH_loc
!------------------------------------------------------------------------------
   subroutine spat_to_qst_loc(this, f, g, h, qLMP, sLMP, tLMP, lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      real(cp), target, intent(inout) :: f(nlat_padded,n_phi_max)
      real(cp), target, intent(inout) :: g(nlat_padded,n_phi_max)
      real(cp), target, intent(inout) :: h(nlat_padded,n_phi_max)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), target, intent(out) :: qLMP(lmP_max)
      complex(cp), target, intent(out) :: sLMP(lmP_max)
      complex(cp), target, intent(out) :: tLMP(lmP_max)

      PERFON("sht_fwd")
      call spat_to_SHqst_l(sht_lP, f, g, h, qLMP, sLMP, tLMP, lcut+1)
      PERFOFF

   end subroutine spat_to_qst_loc
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor_loc(this, f, g, fLMP, gLMP, lcut)

      !-- Input variables
      class(type_shtns_loc), intent(inout) :: this
      real(cp), target, intent(inout) :: f(nlat_padded,n_phi_max)
      real(cp), target, intent(inout) :: g(nlat_padded,n_phi_max)
      integer,  intent(in) :: lcut

      !-- Output variables
      complex(cp), target, intent(out) :: fLMP(lmP_max)
      complex(cp), target, intent(out) :: gLMP(lmP_max)

      PERFON("sht_fwd")
      call spat_to_SHsphtor_l(sht_lP, f, g, fLMP, gLMP, lcut+1)
      PERFOFF

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

       PERFON("sht_fwd")
      call SHtor_to_spat_ml(sht_l, 0, fl_ax, tmpt, tmpp, l_max)
      ft(:)=real(tmpt(:))
      fp(:)=real(tmpp(:))
      PERFOFF

   end subroutine toraxi_to_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH_axi_loc(this, f, fLM)

      class(type_shtns_loc), intent(inout) :: this
      real(cp), intent(in) :: f(nlat_padded)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp(nlat_padded)
      complex(cp) :: tmpLM(size(fLM))

      PERFON("sht_fwd")
      tmp(:)=cmplx(f(:),0.0_cp,kind=cp)
      if ( size(fLM) == l_max+2 ) then
         call spat_to_SH_ml(sht_lP, 0, tmp, tmpLM, l_max+1)
      else if ( size(fLM) == l_max+1 ) then
         call spat_to_SH_ml(sht_l, 0, tmp, tmpLM, l_max)
      end if
      fLM(:)=real(tmpLM(:))
      PERFOFF

   end subroutine spat_to_SH_axi_loc
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!   
!   Θ-Distributed Forward Transforms
!   
!------------------------------------------------------------------------------

   subroutine scal_to_spat_dist(this, Slm, fieldc, lcut)
      !   
      !   Transform the spherical harmonic coefficients Qlm into its spatial 
      !   representation Vr
      !  
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      
      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp),  target, intent(inout) :: Slm(n_lm_loc)
      
      !-- Output variables
      real(cp),     target, intent(out)   :: fieldc(n_theta_loc,n_phi_max)
      integer,      intent(in)    :: lcut

      !-- Local variables
      complex(cp) :: fL_loc(n_theta_max,n_m_loc)
      
      PERFON("sht_bwd")
      call scal_to_hyb(Slm, fL_loc, lcut)
      call transform_m2phi(fL_loc, fieldc)
      PERFOFF
      
   end subroutine scal_to_spat_dist
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat_dist(this, Slm, gradtc, gradpc, lcut)
      !
      !   Transform a scalar spherical harmonic field into it's gradient
      !   on the grid
      !
      !   Author: Rafael Lago, MPCDF, April 2020
      !

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp), target, intent(inout) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: gradtc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: gradpc(n_theta_loc,n_phi_max)
      
      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      PERFON("sht_bwd")
      call scal_to_grad_hyb(Slm, fL, gL, lcut)
      call transform_m2phi(fL, gradtc)
      call transform_m2phi(gL, gradpc)
      PERFOFF

   end subroutine scal_to_grad_spat_dist
!------------------------------------------------------------------------------
   subroutine torpol_to_spat_dist(this, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp), target, intent(in) :: Wlm(n_lm_loc)
      complex(cp), target, intent(inout) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: vrc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: vtc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: vpc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      complex(cp) :: hL(n_theta_max,n_m_loc)
      
      PERFON("sht_bwd")
      call torpol_to_hyb(Wlm, dWlm, Zlm, fL, gL, hL, lcut)
      call transform_m2phi(fL, vrc)
      call transform_m2phi(gL, vtc)
      call transform_m2phi(hL, vpc)
      PERFOFF

   end subroutine torpol_to_spat_dist
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_dist(this, or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, &
              &                        cvpc, lcut)

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp), target, intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc), Jlm(n_lm_loc)
      complex(cp), target, intent(inout) :: dJlm(n_lm_loc)
      real(cp),    intent(in) :: or2
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: cvrc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: cvtc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: cvpc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: Qlm(0:l_max), Tlm(0:l_max)
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      complex(cp) :: hL(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      PERFON("sht_bwd")
      call torpol_to_curl_hyb(or2, Blm, ddBlm, Jlm, dJlm, fL, gL, hL, lcut)
      call transform_m2phi(fL, cvrc)
      call transform_m2phi(gL, cvtc)
      call transform_m2phi(hL, cvpc)
      PERFOFF
      
   end subroutine torpol_to_curl_spat_dist
!------------------------------------------------------------------------------
   subroutine pol_to_grad_spat_dist(this, Slm, gradtc, gradpc, lcut)

       !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp), target, intent(in) :: Slm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: gradtc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: gradpc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      PERFON("sht_bwd")
      call pol_to_grad_hyb(Slm, fL, gL, lcut)
      call transform_m2phi(fL, gradtc)
      call transform_m2phi(gL, gradpc)
      PERFOFF

   end subroutine pol_to_grad_spat_dist
!------------------------------------------------------------------------------
   subroutine torpol_to_dphspat_dist(this, dWlm, Zlm, dvtdp, dvpdp, lcut)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp), target, intent(in) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variables
      real(cp), target, intent(out) :: dvtdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
      real(cp), target, intent(out) :: dvpdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      complex(cp) :: gL(n_theta_max,n_m_loc)
      
      PERFON("sht_bwd")
      call torpol_to_dphhyb(dWlm, Zlm, fL, gL, lcut)      
      call transform_m2phi(fL, dvtdp)
      call transform_m2phi(gL, dvpdp)
      PERFOFF
      
   end subroutine torpol_to_dphspat_dist
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_dist(this, Qlm, cvrc, lcut)

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      complex(cp), target, intent(in) :: Qlm(n_lm_loc)
      integer,     intent(in) :: lcut

      !-- Output variable
      real(cp), target, intent(out) :: cvrc(n_theta_loc,n_phi_max)

      !-- Local variables
      complex(cp) :: fL(n_theta_max,n_m_loc)
      
      PERFON("sht_bwd")
      call pol_to_curlr_hyb(Qlm, fL, lcut)
      call transform_m2phi(fL, cvrc)
      PERFOFF
      
   end subroutine pol_to_curlr_spat_dist

   
!------------------------------------------------------------------------------
!   
!   Θ-Distributed Backward Transforms
!   
!------------------------------------------------------------------------------
  
!----------------------------------------------------------------------------
   subroutine scal_to_SH_dist(this, f, fLMP, lcut)

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      real(cp),    target, intent(inout) :: f(n_theta_loc,n_phi_max)
      integer,     intent(in)    :: lcut

      !-- Output variable
      complex(cp), target, intent(out)   :: fLMP(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      integer :: m, l_lm, u_lm, i
      
      PERFON("sht_fwd")
      call transform_phi2m(f, fL_loc)
      call hyb_to_SH(fL_loc, fLMP, lcut)
      PERFOFF
      
   end subroutine scal_to_SH_dist
!------------------------------------------------------------------------------
   subroutine spat_to_qst_dist(this, f, g, h, qLMP, sLMP, tLMP, lcut)

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      real(cp), target, intent(inout) :: f(n_theta_loc,n_phi_max)
      real(cp), target, intent(inout) :: g(n_theta_loc,n_phi_max)
      real(cp), target, intent(inout) :: h(n_theta_loc,n_phi_max)
      integer,  intent(in)    :: lcut

      !-- Output variables
      complex(cp), target, intent(out) :: qLMP(n_lmP_loc)
      complex(cp), target, intent(out) :: sLMP(n_lmP_loc)
      complex(cp), target, intent(out) :: tLMP(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  gL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  hL_loc(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m

      PERFON("sht_fwd")
      !>@TODO Vectorial FFT and transpose (f,g,h at once)
      call transform_phi2m(f, fL_loc)
      call transform_phi2m(g, gL_loc)
      call transform_phi2m(h, hL_loc)
      call hyb_to_qst(fL_loc, gL_loc, hL_loc, qLMP, sLMP, tLMP, lcut)
      PERFOFF

   end subroutine spat_to_qst_dist
!------------------------------------------------------------------------------
   subroutine spat_to_sphertor_dist(this, f, g, fLMP, gLMP, lcut)

      !-- Input variables
      class(type_shtns_dist), intent(inout) :: this
      real(cp), target, intent(inout) :: f(n_theta_loc,n_phi_max)
      real(cp), target, intent(inout) :: g(n_theta_loc,n_phi_max)
      integer,  intent(in)    :: lcut

      !-- Output variables
      complex(cp), target, intent(out) :: fLMP(n_lmP_loc)
      complex(cp), target, intent(out) :: gLMP(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  fL_loc(n_theta_max,n_m_loc)
      complex(cp) ::  gL_loc(n_theta_max,n_m_loc)
      integer :: i, l_lm, u_lm, m
      
      PERFON("sht_fwd")
      call transform_phi2m(f, fL_loc)
      call transform_phi2m(g, gL_loc)
      call hyb_to_sphertor(fL_loc, gL_loc, fLMP, gLMP, lcut)
      PERFOFF

   end subroutine spat_to_sphertor_dist
!------------------------------------------------------------------------------
   subroutine spat_to_SH_axi_dist(this, f, fLM)
      
      class(type_shtns_dist), intent(inout) :: this
      real(cp), intent(in)  :: f(nThetaStart:nThetaStop)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp_c(n_theta_max)
      real(cp)    :: tmp_r(n_theta_max)
      complex(cp) :: tmpLM(size(fLM))
      integer     :: ierr

      PERFON("sht_fwd")
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
      PERFOFF

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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_bwd")
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
      PERFOFF
      
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
      
      PERFON("leg_fwd")
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
      PERFOFF

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

      PERFON("leg_fwd")
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
      PERFOFF

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

      PERFON("leg_fwd")
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
      PERFOFF

   end subroutine hyb_to_sphertor
!------------------------------------------------------------------------------


!--------------------------------
!
!-- type_sht procedures
!
!--------------------------------

!------------------------------------------------------------------------------
   subroutine initialize_buff(this, max_buff)
      class(type_shtns_buff) :: this
      integer, intent(in) :: max_buff
      
      this%max_buff = max_buff
      this%n_sh = 0
      this%n_spat = 0
      
      ! TODO: figure out the largest buffer needed...
      allocate(this%spat_ptr(2*max_buff))
      allocate(this%sh_ptr(2*max_buff))
      allocate(this%operations(2*max_buff))
      allocate(this%lcut(2*max_buff))
      allocate(this%or2(2*max_buff))
      
      this%operations = OP_NONE
      
   end subroutine initialize_buff

!------------------------------------------------------------------------------   
   subroutine finalize_buff(this)
      type(type_shtns_buff) :: this
      integer :: i
      
      do i=1, size(this%spat_ptr)
         if (associated(this%spat_ptr(i)%p)) nullify(this%spat_ptr(i)%p)
         if (associated(this%sh_ptr(i)%p)) nullify(this%sh_ptr(i)%p)
      end do
      deallocate(this%spat_ptr)
      deallocate(this%sh_ptr)
      deallocate(this%operations)
      deallocate(this%lcut)
      deallocate(this%or2)
      
      this%max_buff = 0
      this%n_sh = 0
      this%n_spat = 0
      
   end subroutine finalize_buff

!------------------------------------------------------------------------------   
#define __CHECK_BUFFERSIZE(__X, __NEEDED)\
   if (this%operations(1)>0) then;\
      print *, "Error in "//__X//"; commit transform before changing direction!";\
      test_ptr(5) = 5;\
      stop;\
   end if;\
   if ((this%n_spat+__NEEDED)>this%max_buff) call commit_backward_buff(this);\

!------------------------------------------------------------------------------      
   subroutine scal_to_spat_buff(this, Slm, fieldc, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp), target, intent(inout) :: Slm(n_lm_loc)
      real(cp),    target, intent(out)   :: fieldc(n_theta_loc,n_phi_max)
      integer,             intent(in)    :: lcut
      
      __CHECK_BUFFERSIZE("scal_to_spat_buff", 1)
      
      this%sh_ptr(this%n_sh+1)%p     => Slm
      this%spat_ptr(this%n_spat+1)%p => fieldc
      this%operations(this%n_spat+1) =  OP_SCAL2SPAT
      this%lcut(this%n_spat+1)       =  lcut
      
      this%n_sh = this%n_sh + 1
      this%n_spat  = this%n_spat + 1
      
   end subroutine
   
!------------------------------------------------------------------------------   
   subroutine scal_to_grad_spat_buff(this, Slm, gradtc, gradpc, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp), target, intent(inout) :: Slm(n_lm_loc)
      real(cp),    target, intent(out)   :: gradtc(n_theta_loc,n_phi_max)
      real(cp),    target, intent(out)   :: gradpc(n_theta_loc,n_phi_max)
      integer,             intent(in)    :: lcut
      
      __CHECK_BUFFERSIZE("scal_to_grad_spat_buff", 2)
      
      this%sh_ptr(this%n_sh+1)%p     => Slm
      this%spat_ptr(this%n_spat+1)%p => gradtc
      this%spat_ptr(this%n_spat+2)%p => gradpc
      this%operations(this%n_spat+1) =  OP_SCAL2GRAD_SPAT
      this%lcut(this%n_spat+1)       =  lcut
      
      this%n_sh = this%n_sh + 1
      this%n_spat  = this%n_spat + 2
      
   end subroutine
   
!------------------------------------------------------------------------------   
   subroutine torpol_to_spat_buff(this, Wlm, dWlm, Zlm, vrc, vtc, vpc, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp),    target, intent(in)    :: Wlm(n_lm_loc)
      complex(cp),    target, intent(inout) :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      real(cp), target, intent(out) :: vrc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: vtc(n_theta_loc,n_phi_max)
      real(cp), target, intent(out) :: vpc(n_theta_loc,n_phi_max)
      integer,          intent(in)  :: lcut
      
      __CHECK_BUFFERSIZE("torpol_to_spat_buff", 3)
      
      this%sh_ptr(this%n_sh+1)%p     => Wlm
      this%sh_ptr(this%n_sh+2)%p     => dWlm
      this%sh_ptr(this%n_sh+3)%p     => Zlm
      this%spat_ptr(this%n_spat+1)%p => vrc
      this%spat_ptr(this%n_spat+2)%p => vtc
      this%spat_ptr(this%n_spat+3)%p => vpc
      this%operations(this%n_spat+1) =  OP_TORPOL2SPAT
      this%lcut(this%n_spat+1)       =  lcut
      
      this%n_sh = this%n_sh + 3
      this%n_spat  = this%n_spat + 3
      
   end subroutine
   
!------------------------------------------------------------------------------   
   subroutine torpol_to_curl_spat_buff(this, or2, Blm, ddBlm, Jlm, dJlm, cvrc, cvtc, &
            &                        cvpc, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp),    target, intent(in) :: Blm(n_lm_loc), ddBlm(n_lm_loc)
      complex(cp),    target, intent(in) :: Jlm(n_lm_loc)
      complex(cp),    target, intent(inout) :: dJlm(n_lm_loc)
      real(cp),       target, intent(out) :: cvrc(n_theta_loc,n_phi_max)
      real(cp),       target, intent(out) :: cvtc(n_theta_loc,n_phi_max)
      real(cp),       target, intent(out) :: cvpc(n_theta_loc,n_phi_max)
      real(cp),                intent(in) :: or2
      integer,                 intent(in) :: lcut
      
      __CHECK_BUFFERSIZE("torpol_to_curl_spat_buff", 3)
      
      this%sh_ptr(this%n_sh+1)%p     => Blm
      this%sh_ptr(this%n_sh+2)%p     => ddBlm
      this%sh_ptr(this%n_sh+3)%p     => Jlm
      this%sh_ptr(this%n_sh+4)%p     => dJlm
      this%spat_ptr(this%n_spat+1)%p => cvrc
      this%spat_ptr(this%n_spat+2)%p => cvtc
      this%spat_ptr(this%n_spat+3)%p => cvpc
      this%operations(this%n_spat+1) =  OP_TORPOL2CURL_SPAT
      this%lcut(this%n_spat+1)       =  lcut
      this%or2(this%n_spat+1)        =  or2
      
      this%n_sh = this%n_sh + 4
      this%n_spat  = this%n_spat + 3
      
   end subroutine

!------------------------------------------------------------------------------      
   subroutine pol_to_grad_spat_buff(this, Slm, gradtc, gradpc, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp), target, intent(in)  :: Slm(n_lm_loc)
      real(cp),    target, intent(out) :: gradtc(n_theta_loc,n_phi_max)
      real(cp),    target, intent(out) :: gradpc(n_theta_loc,n_phi_max)
      integer,             intent(in)  :: lcut
      
      __CHECK_BUFFERSIZE("pol_to_grad_spat_buff", 2)
      
      this%sh_ptr(this%n_sh+1)%p     => Slm
      this%spat_ptr(this%n_spat+1)%p => gradtc
      this%spat_ptr(this%n_spat+2)%p => gradpc
      this%operations(this%n_spat+1) =  OP_POL2GRAD_SPAT
      this%lcut(this%n_spat+1)       =  lcut
      
      this%n_sh   = this%n_sh + 1
      this%n_spat = this%n_spat + 2
   end subroutine
   
!------------------------------------------------------------------------------   
   subroutine torpol_to_dphspat_buff(this, dWlm, Zlm, dvtdp, dvpdp, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp), target, intent(in)  :: dWlm(n_lm_loc), Zlm(n_lm_loc)
      real(cp),    target, intent(out) :: dvtdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
      real(cp),    target, intent(out) :: dvpdp(nThetaStart:nThetaStop,n_phi_max) ! Careful with dimensions here!
      integer,     intent(in) :: lcut
      
      __CHECK_BUFFERSIZE("torpol_to_dphspat_buff", 2)
      
      this%sh_ptr(this%n_sh+1)%p     => dWlm
      this%sh_ptr(this%n_sh+2)%p     => Zlm
      this%spat_ptr(this%n_spat+1)%p => dvtdp
      this%spat_ptr(this%n_spat+2)%p => dvpdp
      this%operations(this%n_spat+1) =  OP_TORPOL2DPHSPAT
      this%lcut(this%n_spat+1)       =  lcut
      
      this%n_sh   = this%n_sh + 2
      this%n_spat = this%n_spat + 2
      
   end subroutine
      
!------------------------------------------------------------------------------   
   subroutine pol_to_curlr_spat_buff(this, Qlm, cvrc, lcut)
      class(type_shtns_buff), intent(inout) :: this
      complex(cp), target, intent(in)  :: Qlm(n_lm_loc)
      real(cp),    target, intent(out) :: cvrc(n_theta_loc,n_phi_max)
      integer,             intent(in)  :: lcut
      
      __CHECK_BUFFERSIZE("pol_to_curlr_spat_buff", 1)
      
      this%sh_ptr(this%n_sh+1)%p     => Qlm
      this%spat_ptr(this%n_spat+1)%p => cvrc
      this%operations(this%n_spat+1) =  OP_POL2CURLR_SPAT
      this%lcut(this%n_spat+1)       =  lcut
      
      this%n_sh   = this%n_sh + 1
      this%n_spat = this%n_spat + 1
      
   end subroutine
#undef __CHECK_BUFFERSIZE

!------------------------------------------------------------------------------   
#define __CHECK_BUFFERSIZE(__X, __NEEDED)\
   if (this%operations(1)<0) then;\
      print *, "Error in "//__X//"; commit transform before changing direction!";\
      test_ptr(5) = 5;\
      stop;\
   end if;\
   if ((this%n_sh+__NEEDED)>this%max_buff) call commit_forward_buff(this);\

   
!------------------------------------------------------------------------------   
   subroutine scal_to_SH_buff(this, f, fLMP, lcut)
      class(type_shtns_buff), intent(inout) :: this
      real(cp),    target, intent(inout) :: f(n_theta_loc,n_phi_max)
      complex(cp), target, intent(out)   :: fLMP(n_lmP_loc)
      integer,             intent(in)    :: lcut
      
      __CHECK_BUFFERSIZE("scal_to_SH_buff", 1)
      this%spat_ptr(this%n_spat+1)%p => f
      this%sh_ptr(this%n_sh+1)%p     => fLMP
      this%operations(this%n_sh+1)   =  OP_SPAT2SH
      this%lcut(this%n_sh+1)         =  lcut
      
      this%n_spat  = this%n_spat + 1
      this%n_sh = this%n_sh + 1
      
   end subroutine

!------------------------------------------------------------------------------   
   subroutine spat_to_qst_buff(this, f, g, h, qLMP, sLMP, tLMP, lcut)
      class(type_shtns_buff), intent(inout) :: this
      real(cp),    target, intent(inout) :: f(n_theta_loc,n_phi_max)
      real(cp),    target, intent(inout) :: g(n_theta_loc,n_phi_max)
      real(cp),    target, intent(inout) :: h(n_theta_loc,n_phi_max)
      complex(cp), target, intent(out)   :: qLMP(n_lmP_loc)
      complex(cp), target, intent(out)   :: sLMP(n_lmP_loc)
      complex(cp), target, intent(out)   :: tLMP(n_lmP_loc)
      integer,  intent(in)    :: lcut
      
      __CHECK_BUFFERSIZE("scal_to_SH_buff", 3)
      this%spat_ptr(this%n_spat+1)%p => f
      this%spat_ptr(this%n_spat+2)%p => g
      this%spat_ptr(this%n_spat+3)%p => h
      this%sh_ptr(this%n_sh+1)%p     => qLMP
      this%sh_ptr(this%n_sh+2)%p     => sLMP
      this%sh_ptr(this%n_sh+3)%p     => tLMP
      this%operations(this%n_sh+1)   =  OP_SPAT2QST
      this%lcut(this%n_sh+1)         =  lcut
      
      this%n_spat  = this%n_spat + 3
      this%n_sh = this%n_sh + 3
   end subroutine
   
!------------------------------------------------------------------------------   
   subroutine spat_to_sphertor_buff(this, f, g, fLMP, gLMP, lcut)
      class(type_shtns_buff), intent(inout) :: this
      real(cp),    target, intent(inout) :: f(n_theta_loc,n_phi_max)
      real(cp),    target, intent(inout) :: g(n_theta_loc,n_phi_max)
      complex(cp), target, intent(out)   :: fLMP(n_lmP_loc)
      complex(cp), target, intent(out)   :: gLMP(n_lmP_loc)
      integer,  intent(in)    :: lcut
      
      __CHECK_BUFFERSIZE("scal_to_SH_buff", 2)
      this%spat_ptr(this%n_spat+1)%p => f
      this%spat_ptr(this%n_spat+2)%p => g
      this%sh_ptr(this%n_sh+1)%p     => fLMP
      this%sh_ptr(this%n_sh+2)%p     => gLMP
      this%operations(this%n_sh+1)   =  OP_SPAT2QST
      this%lcut(this%n_sh+1)         =  lcut
      
      this%n_spat  = this%n_spat + 2
      this%n_sh = this%n_sh + 2
   end subroutine
   
!------------------------------------------------------------------------------
   subroutine spat_to_SH_axi_buff(this, f, fLM)
      
      class(type_shtns_buff), intent(inout) :: this
      real(cp), intent(in)  :: f(nThetaStart:nThetaStop)
      real(cp), intent(out) :: fLM(:)

      !-- Local arrays
      complex(cp) :: tmp_c(n_theta_max)
      real(cp)    :: tmp_r(n_theta_max)
      complex(cp) :: tmpLM(size(fLM))
      integer     :: ierr

      PERFON("sht_fwd")
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
      PERFOFF

   end subroutine spat_to_SH_axi_buff
#undef __CHECK_BUFFERSIZE

!------------------------------------------------------------------------------   
   subroutine commit_backward_buff(this)
      class(type_shtns_buff) :: this
      complex(cp), allocatable :: hyb_buffer(:,:,:), fft_buffer(:,:,:)
      integer :: ic, ir, fftlen, ierr
      
      PERFON("sht_bwd")
      
      fftlen = max(n_m_max, n_phi_max/2+1)
      if (this%n_spat==0) return
      
      !-- Only allocates what is needed
      allocate(hyb_buffer(n_theta_max, n_m_loc, this%n_spat))
      allocate(fft_buffer(n_theta_loc, fftlen, this%n_spat))
      
      !-- Legendre transform
      ir = 1
      ic = 1
      do while (.TRUE.)
         select case(this%operations(ir))
            ! --------------------------------------------
            case (OP_SCAL2SPAT)
               call scal_to_hyb(this%sh_ptr(ic)%p, hyb_buffer(:,:,ir), this%lcut(ir))
               ic = ic + 1
               ir = ir + 1
            
            ! --------------------------------------------
            case (OP_SCAL2GRAD_SPAT)
               call scal_to_grad_hyb(this%sh_ptr(ic)%p, hyb_buffer(:,:,ir), hyb_buffer(:,:,ir+1), this%lcut(ir))
               ic = ic + 1
               ir = ir + 2
            
            ! --------------------------------------------
            case (OP_TORPOL2SPAT)
               call torpol_to_hyb(this%sh_ptr(ic)%p, this%sh_ptr(ic+1)%p, this%sh_ptr(ic+2)%p, &
                  &               hyb_buffer(:,:,ir), hyb_buffer(:,:,ir+1), hyb_buffer(:,:,ir+2), this%lcut(ir))
               ic = ic + 3
               ir = ir + 3
               
            ! --------------------------------------------
            case (OP_TORPOL2CURL_SPAT)
               call torpol_to_curl_hyb(this%or2(ir), this%sh_ptr(ic)%p, this%sh_ptr(ic+1)%p, &
                  & this%sh_ptr(ic+2)%p, this%sh_ptr(ic+3)%p, &
                  & hyb_buffer(:,:,ir),  hyb_buffer(:,:,ir+1),  hyb_buffer(:,:,ir+2), this%lcut(ir))
               ic = ic + 4
               ir = ir + 3
               
            ! --------------------------------------------
            case (OP_POL2GRAD_SPAT)
               call pol_to_grad_hyb(this%sh_ptr(ic)%p, &
                  & hyb_buffer(:,:,ir), hyb_buffer(:,:,ir+1), this%lcut(ir))
               ic = ic + 1
               ir = ir + 2
            
            ! --------------------------------------------   
            case (OP_POL2CURLR_SPAT)
               call pol_to_curlr_hyb(this%sh_ptr(ic)%p, hyb_buffer(:,:,ir), this%lcut(ir))
               ic = ic + 1
               ir = ir + 1
               
            ! --------------------------------------------   
            case (OP_TORPOL2DPHSPAT)
               call torpol_to_dphhyb(this%sh_ptr(ic)%p, this%sh_ptr(ic+1)%p, &
                  & hyb_buffer(:,:,ir), hyb_buffer(:,:,ir+1), this%lcut(ir))
               ic = ic + 2
               ir = ir + 2
               
            case default
               print *, " * Unknown operation in commit_buff: ", this%operations(ir), ". Aborting..."
               stop
         end select
         
         if (ir>this%n_spat) exit
      end do
      
      !-- Transpose from m_loc to th_loc
      !   TODO this inplace!
      call transpose_m2th(hyb_buffer, fft_buffer, this%n_spat)
      
      !-- FFT from hyb to phi
      do ir=1,this%n_spat
         fft_buffer(:,n_m_max+1:,ir) = zero
         call ifft_many(fft_buffer(:,:,ir), this%spat_ptr(ir)%p)
         nullify( this%spat_ptr(ir)%p )
      end do
      
      !-- Cleanup of complex (input) pointers
      do ic=1,this%n_sh
         nullify( this%sh_ptr(ic)%p )
      end do
      
      this%n_sh = 0
      this%n_spat = 0
      this%operations = OP_NONE
      deallocate(hyb_buffer)
      deallocate(fft_buffer)
      
      PERFOFF
      
   end subroutine commit_backward_buff

   
   !------------------------------------------------------------------------------   
   subroutine commit_forward_buff(this)
      !   
      !   Commits the forward transforms in the queue.
      !   So far, the forward transform always has n_spat = n_sh, so we use only n_sh here.
      !   
      !
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      class(type_shtns_buff) :: this
      complex(cp), allocatable :: hyb_buffer(:,:,:), fft_buffer(:,:,:)
      integer :: i, fftlen
      
      PERFON("sht_fwd")
      
      fftlen = max(n_m_max, n_phi_max/2+1)
      if (this%n_sh==0) return
      
      !-- Only allocates what is needed
      allocate(hyb_buffer(n_theta_max, n_m_loc, this%n_sh))
      allocate(fft_buffer(n_theta_loc, fftlen, this%n_sh))
      
      !-- FFT from phi to hyb
      do i=1,this%n_sh
         fft_buffer(:,n_m_max+1:,i) = zero
         call fft_many(this%spat_ptr(i)%p, fft_buffer(:,:,i))
      end do
      
      !-- Transpose from th_loc to m_loc
      call transpose_th2m(fft_buffer, hyb_buffer, this%n_sh)
      
      !-- Legendre transform
      i = 1
      do while (.TRUE.)
         select case(this%operations(i))
            ! --------------------------------------------
            case (OP_SPAT2SH)
               call hyb_to_SH(hyb_buffer(:,:,i), this%sh_ptr(i)%p, this%lcut(i))
               i = i + 1

            ! --------------------------------------------
            case (OP_SPAT2QST)
               call hyb_to_qst(hyb_buffer(:,:,i), hyb_buffer(:,:,i+1), hyb_buffer(:,:,i+2), &
                & this%sh_ptr(i)%p, this%sh_ptr(i+1)%p, this%sh_ptr(i+2)%p, this%lcut(i))
               i = i + 3
            
            ! --------------------------------------------
            case (OP_SPAT2SPHERTOR)
               call hyb_to_sphertor(hyb_buffer(:,:,i), hyb_buffer(:,:,i+1), &
                & this%sh_ptr(i)%p, this%sh_ptr(i+1)%p, this%lcut(i))
               i = i + 2
!                
!             ! --------------------------------------------
!             ! TODO: is it worth to do it here? I don't think so...
!             case (OP_SPAT2SH_AXI) 
               
            case default
               print *, " * Unknown operation in commit_buff: ", this%operations(i), ". Aborting..."
               stop
         end select
         
         if (i>this%n_spat) exit
      end do
      
      !-- Cleanup of (input & output) pointers
      do i=1,this%n_sh
         nullify( this%sh_ptr(i)%p )
         nullify( this%spat_ptr(i)%p )
      end do
      
      this%n_sh = 0
      this%n_spat = 0
      this%operations = OP_NONE
      deallocate(hyb_buffer)
      deallocate(fft_buffer)
      PERFOFF
      
   end subroutine commit_forward_buff
   
   !------------------------------------------------------------------------------   
   subroutine commit_dummy(this)
      !   
      !   Does nothing
      !
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      class(type_shtns_if) :: this
      return
   end subroutine commit_dummy

end module sht
