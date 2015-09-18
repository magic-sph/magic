#include "perflib_preproc.cpp"
module rIterThetaBlocking_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
#ifdef WITH_SHTNS
  use truncation, only: n_phi_max, n_theta_max
  use shtns
#endif

   use rIteration_mod, only: rIteration_t
   use precision_mod
   use truncation, only: lm_max,lmP_max,nrp,l_max,lmP_max_dtB, &
        & n_phi_maxStr,n_theta_maxStr,n_r_maxStr,lm_maxMag
   use blocking, only: nfs
   use logic, only: l_mag,l_conv,l_mag_kin,l_heat,l_ht,l_anel,l_mag_LF,        &
        & l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic, l_cond_ic,    &
        & l_rot_ma, l_cond_ma, l_dtB, l_store_frame, l_movie_oc, l_TO
   use radial_data,only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: or2, orho1
   use fft
   use legendre_spec_to_grid, only: legTFG, legTFGnomag
   use leg_helper_mod, only: leg_helper_t
   use nonlinear_lm_mod, only:nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use physical_parameters, only: kbots,ktops,n_r_LCR
   use nonlinear_bcs, only: v_rigid_boundary
   use legendre_grid_to_spec

   implicit none

   private

   type, public :: dtB_arrays_t
      !----- Local dtB output stuff:
      complex(cp), allocatable :: BtVrLM(:),BpVrLM(:),BrVtLM(:),BrVpLM(:), &
                                  &               BtVpLM(:), BpVtLM(:)
      complex(cp), allocatable :: BtVpCotLM(:),BpVtCotLM(:),BtVpSn2LM(:), &
                                  &               BpVtSn2LM(:)
      complex(cp), allocatable :: BrVZLM(:),BtVZLM(:),BtVZcotLM(:),       &
                                  &               BtVZsn2LM(:)
   end type dtB_arrays_t

   type, public :: TO_arrays_t
      !----- Local TO output stuff:
      real(cp), allocatable :: dzRstrLM(:),dzAstrLM(:)
      real(cp), allocatable :: dzCorLM(:),dzLFLM(:)
   end type TO_arrays_t

   type, public, abstract, extends(rIteration_t) :: rIterThetaBlocking_t
      ! or with len parameters for the theta block size and number
      !type,public,extends(rIteration_t) :: rIterThetaBlocking_t(sizeThetaB,nThetaBs)
      !integer,len :: sizeThetaB,nThetaBs
      integer :: sizeThetaB, nThetaBs

      !type(nonlinear_lm_t) :: nl_lm
      type(leg_helper_t) :: leg_helper
      type(dtB_arrays_t) :: dtB_arrays
      type(TO_arrays_t)  :: TO_arrays

      !class(grid_space_arrays_t),private :: gsa

      !----- Saved magnetic field components from last time step:
      !      This is needed for the current TO version. However,
      !      the variables calulated with this don't give any
      !      deep insight. TO should be changes in the future to
      !      eliminate this.
      real(cp), allocatable :: BsLast(:,:,:), BpLast(:,:,:), BzLast(:,:,:)

   contains

      !procedure :: initialize => initialize_rIterThetaBlocking
      procedure :: allocate_common_arrays
      procedure :: deallocate_common_arrays
      procedure :: set_ThetaBlocking
      !procedure,deferred :: do_iteration
      procedure :: transform_to_grid_space
      procedure :: transform_to_lm_space

   end type rIterThetaBlocking_t

contains

   subroutine allocate_common_arrays(this)

      class(rIterThetaBlocking_t) :: this

      !----- Nonlinear terms in lm-space:
      !call this%nl_lm%initialize(lmP_max)


      !----- Help arrays for Legendre transform calculated in legPrepG:
      !      Parallelizatio note: these are the R-distributed versions
      !      of the field scalars.
      call this%leg_helper%initialize(lm_max,lm_maxMag,l_max)

      !----- Local dtB output stuff:
      allocate( this%dtB_arrays%BtVrLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BpVrLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BrVtLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BrVpLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BtVpLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BpVtLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BtVpCotLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BpVtCotLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BtVpSn2LM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BpVtSn2LM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BrVZLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BtVZLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BtVZcotLM(lmP_max_dtB) )
      allocate( this%dtB_arrays%BtVZsn2LM(lmP_max_dtB) )

      !----- Local TO output stuff:
      if ( l_TO ) then
        allocate( this%TO_arrays%dzRstrLM(l_max+2),this%TO_arrays%dzAstrLM(l_max+2) )
        allocate( this%TO_arrays%dzCorLM(l_max+2),this%TO_arrays%dzLFLM(l_max+2) )
      end if

      allocate( this%BsLast(n_phi_maxStr,n_theta_maxStr,nRstart:nRstop) )
      allocate( this%BpLast(n_phi_maxStr,n_theta_maxStr,nRstart:nRstop) )
      allocate( this%BzLast(n_phi_maxStr,n_theta_maxStr,nRstart:nRstop) )

   end subroutine allocate_common_arrays
!-------------------------------------------------------------------------------
   subroutine deallocate_common_arrays(this)

      class(rIterThetaBlocking_t) :: this

      deallocate( this%dtB_arrays%BtVrLM )
      deallocate( this%dtB_arrays%BpVrLM )
      deallocate( this%dtB_arrays%BrVtLM )
      deallocate( this%dtB_arrays%BrVpLM )
      deallocate( this%dtB_arrays%BtVpLM )
      deallocate( this%dtB_arrays%BpVtLM )
      deallocate( this%dtB_arrays%BtVpCotLM )
      deallocate( this%dtB_arrays%BpVtCotLM )
      deallocate( this%dtB_arrays%BtVpSn2LM )
      deallocate( this%dtB_arrays%BpVtSn2LM )
      deallocate( this%dtB_arrays%BrVZLM )
      deallocate( this%dtB_arrays%BtVZLM )
      deallocate( this%dtB_arrays%BtVZcotLM )
      deallocate( this%dtB_arrays%BtVZsn2LM )

      !----- Local TO output stuff:
      if ( l_TO ) then
        deallocate( this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM )
        deallocate( this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM )
      end if

      deallocate( this%BsLast)
      deallocate( this%BpLast)
      deallocate( this%BzLast)

   end subroutine deallocate_common_arrays
!-------------------------------------------------------------------------------
   subroutine set_ThetaBlocking(this,nThetaBs,sizeThetaB)

      class(rIterThetaBlocking_t) :: this
      integer,intent(in) :: nThetaBs, sizeThetaB

      this%nThetaBs = nThetaBs

      this%sizeThetaB = sizeThetaB

   end subroutine set_ThetaBlocking
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space(this,nThetaStart,nThetaStop,gsa)

#ifdef WITH_SHTNS
      use blocking, only: sizeThetaB
      use truncation, only: n_m_max
      use horizontal_data, only: D_mc2m, osn2
#endif

      class(rIterThetaBlocking_t), target :: this
      integer, intent(in) :: nThetaStart,nThetaStop
      type(grid_space_arrays_t) :: gsa

      ! Local variables
      integer :: nTheta
      logical :: DEBUG_OUTPUT=.false.
#ifdef WITH_SHTNS
      integer :: nThetaNHS, nThetaN, nThetaS, mc
      real(cp) :: dmT, swap
#endif

      !----- Legendre transform from (r,l,m) to (r,theta,m):
      !      First version with PlmTF needed for first-touch policy
#ifndef WITH_SHTNS
      if ( l_mag ) then
         !PERFON('legTFG')
         !LIKWID_ON('legTFG')
         call legTFG(this%nBc,this%lDeriv,this%lViscBcCalc,           &
              &      this%lFluxProfCalc,nThetaStart,                  &
              &      gsa%vrc,gsa%vtc,gsa%vpc,gsa%dvrdrc,              &
              &      gsa%dvtdrc,gsa%dvpdrc,gsa%cvrc,                  &
              &      gsa%dvrdtc,gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,     &
              &      gsa%brc,gsa%btc,gsa%bpc,gsa%cbrc,                &
              &      gsa%cbtc,gsa%cbpc,gsa%sc,gsa%drSc,               &
              &      gsa%dsdtc, gsa%dsdpc, gsa%pc,                    &
              &      this%leg_helper)
         !LIKWID_OFF('legTFG')
         !PERFOFF
         if (DEBUG_OUTPUT) then
            do nTheta=1,this%sizeThetaB
               write(*,"(2I3,A,6ES20.12)") this%nR,nTheta,": sum v = ",&
                    &sum(gsa%vrc(:,nTheta))!,sum(vtc(:,nTheta)),sum(vpc(:,nTheta))
            end do
         end if
      else
         !PERFON('legTFGnm')
         !LIKWID_ON('legTFGnm')
         call legTFGnomag(this%nBc,this%lDeriv,this%lViscBcCalc,            & 
              &           this%lFluxProfCalc,nThetaStart,                   &
              &           gsa%vrc,gsa%vtc,gsa%vpc,gsa%dvrdrc,               &
              &           gsa%dvtdrc,gsa%dvpdrc,gsa%cvrc,                   &
              &           gsa%dvrdtc,gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,      &
              &           gsa%sc,gsa%drSc,                                  &
              &           gsa%dsdtc, gsa%dsdpc,gsa%pc,                      &
              &           this%leg_helper)
         !LIKWID_OFF('legTFGnm')
         !PERFOFF
      end if
#endif

      !------ Fourier transform from (r,theta,m) to (r,theta,phi):
      if ( l_conv .or. l_mag_kin ) then
         if ( l_heat ) then
#ifdef WITH_SHTNS
            gsa%sc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_s((nThetaStart-1)*n_phi_max+1:)
#else
            call fft_thetab(gsa%sc,1)
#endif
            if ( this%lViscBcCalc ) then
#ifdef WITH_SHTNS
               gsa%dsdtc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dsdt((nThetaStart-1)*n_phi_max+1:)
               gsa%dsdpc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dsdp((nThetaStart-1)*n_phi_max+1:)
#else
               call fft_thetab(gsa%dsdtc,1)
               call fft_thetab(gsa%dsdpc,1)
#endif
               if (this%nR == n_r_cmb .and. ktops==1) then
                  gsa%dsdtc=0.0_cp
                  gsa%dsdpc=0.0_cp
               end if
               if (this%nR == n_r_icb .and. kbots==1) then
                  gsa%dsdtc=0.0_cp
                  gsa%dsdpc=0.0_cp
               end if
            end if
            if ( this%lFluxProfCalc ) then
#ifdef WITH_SHTNS
               gsa%pc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_p((nThetaStart-1)*n_phi_max+1:)
#else
               call fft_thetab(gsa%pc,1)
#endif
            end if
         end if
         if ( l_HT .or. this%lViscBcCalc ) then
#ifdef WITH_SHTNS
            gsa%drSc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_drs((nThetaStart-1)*n_phi_max+1:)
#else
            call fft_thetab(gsa%drSc,1)
#endif
         endif
         if ( this%nBc == 0 ) then
#ifdef WITH_SHTNS
            gsa%vrc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_vr((nThetaStart-1)*n_phi_max+1:)
            gsa%vtc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_vt((nThetaStart-1)*n_phi_max+1:)
            gsa%vpc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_vp((nThetaStart-1)*n_phi_max+1:)
#else
            call fft_thetab(gsa%vrc,1)
            call fft_thetab(gsa%vtc,1)
            call fft_thetab(gsa%vpc,1)
#endif
            if ( this%lDeriv ) then
#ifdef WITH_SHTNS
               gsa%dvrdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvrdr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvtdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvtdr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvpdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvpdr((nThetaStart-1)*n_phi_max+1:)
#else
               call fft_thetab(gsa%dvrdrc,1)
               call fft_thetab(gsa%dvtdrc,1)
               call fft_thetab(gsa%dvpdrc,1)
#endif
#ifdef WITH_SHTNS
               gsa%cvrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_cvr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvrdtc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvrdt((nThetaStart-1)*n_phi_max+1:)
               gsa%dvrdpc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvrdp((nThetaStart-1)*n_phi_max+1:)
#else
               call fft_thetab(gsa%cvrc,1)
               call fft_thetab(gsa%dvrdtc,1)
               call fft_thetab(gsa%dvrdpc,1)
#endif

#ifdef WITH_SHTNS
               gsa%dvtdpc = gsa%vtc
               call fft_thetab(gsa%dvtdpc, -1)
               gsa%dvpdpc(:, :) = gsa%vpc(:, :)
               call fft_thetab(gsa%dvpdpc, -1)

               nThetaNHS = (nThetaStart-1)/2
               do nThetaN=1, sizeThetaB, 2   ! Loop over thetas for north HS
                  nThetaS   = nThetaN+1       ! same theta but for southern HS
                  nThetaNHS = nThetaNHS+1     ! theta-index of northern hemisph. point
                  do mc=1, n_m_max
                     dmT = D_mc2m(mc) * osn2(nThetaNHS)
                     swap = -dmT*gsa%dvtdpc(2*mc, nThetaN)
                     gsa%dvtdpc(2*mc  , nThetaN) =  dmT*gsa%dvtdpc(2*mc-1, nThetaN)
                     gsa%dvtdpc(2*mc-1, nThetaN) = swap
                     swap = -dmT*gsa%dvtdpc(2*mc, nThetaS)
                     gsa%dvtdpc(2*mc  , nThetaS) =  dmT*gsa%dvtdpc(2*mc-1, nThetaS)
                     gsa%dvtdpc(2*mc-1, nThetaS) = swap

                     swap = -dmT*gsa%dvpdpc(2*mc, nThetaN)
                     gsa%dvpdpc(2*mc  , nThetaN) =  dmT*gsa%dvpdpc(2*mc-1, nThetaN)
                     gsa%dvpdpc(2*mc-1, nThetaN) = swap
                     swap = -dmT*gsa%dvpdpc(2*mc, nThetaS)
                     gsa%dvpdpc(2*mc  , nThetaS) =  dmT*gsa%dvpdpc(2*mc-1, nThetaS)
                     gsa%dvpdpc(2*mc-1, nThetaS) = swap
                  end do
               end do
               !-- Zero out terms with index mc > n_m_max:
               if ( n_m_max < nrp/2 ) then
                  do nThetaN=1, sizeThetaB
                     do mc=2*n_m_max+1, nrp
                        gsa%dvtdpc(mc, nThetaN) = 0.0_cp
                        gsa%dvpdpc(mc, nThetaN) = 0.0_cp
                     end do
                  end do  ! loop over nThetaN (theta)
               end if
#endif
               call fft_thetab(gsa%dvtdpc,1)
               call fft_thetab(gsa%dvpdpc,1)
            end if
         else if ( this%nBc == 1 ) then ! Stress free
            gsa%vrc = 0.0_cp
#ifdef WITH_SHTNS
            gsa%vtc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_vt((nThetaStart-1)*n_phi_max+1:)
            gsa%vpc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_vp((nThetaStart-1)*n_phi_max+1:)
#else
            call fft_thetab(gsa%vtc,1)
            call fft_thetab(gsa%vpc,1)
#endif
            if ( this%lDeriv ) then
               gsa%dvrdtc = 0.0_cp
               gsa%dvrdpc = 0.0_cp
#ifdef WITH_SHTNS
               gsa%dvrdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvrdr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvtdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvtdr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvpdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvpdr((nThetaStart-1)*n_phi_max+1:)
               gsa%cvrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_cvr((nThetaStart-1)*n_phi_max+1:)
#else
               call fft_thetab(gsa%dvrdrc,1)
               call fft_thetab(gsa%dvtdrc,1)
               call fft_thetab(gsa%dvpdrc,1)
               call fft_thetab(gsa%cvrc,1)
#endif
               call fft_thetab(gsa%dvtdpc,1)
               call fft_thetab(gsa%dvpdpc,1)
            end if
         else if ( this%nBc == 2 ) then 
            if ( this%nR == n_r_cmb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaMA,this%lDeriv, &
                    &                gsa%vrc,gsa%vtc,gsa%vpc,gsa%cvrc,gsa%dvrdtc, &
                    &                gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,            &
                    &                nThetaStart)
            else if ( this%nR == n_r_icb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaIC,this%lDeriv, &
                    &                gsa%vrc,gsa%vtc,gsa%vpc,gsa%cvrc,gsa%dvrdtc, &
                    &                gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,            &
                    &                nThetaStart)
            end if
            if ( this%lDeriv ) then
#ifdef WITH_SHTNS
               gsa%dvrdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvrdr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvtdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvtdr((nThetaStart-1)*n_phi_max+1:)
               gsa%dvpdrc(1:n_phi_max, 1:nfs) => &
                   this%leg_helper%shtns_dvpdr((nThetaStart-1)*n_phi_max+1:)
#else
               call fft_thetab(gsa%dvrdrc,1)
               call fft_thetab(gsa%dvtdrc,1)
               call fft_thetab(gsa%dvpdrc,1)
#endif
            end if
         end if
      end if
      if ( l_mag .or. l_mag_LF ) then
#ifdef WITH_SHTNS
         gsa%brc(1:n_phi_max, 1:nfs) => &
             this%leg_helper%shtns_br((nThetaStart-1)*n_phi_max+1:)
         gsa%btc(1:n_phi_max, 1:nfs) => &
             this%leg_helper%shtns_bt((nThetaStart-1)*n_phi_max+1:)
         gsa%bpc(1:n_phi_max, 1:nfs) => &
             this%leg_helper%shtns_bp((nThetaStart-1)*n_phi_max+1:)
#else
         call fft_thetab(gsa%brc,1)
         call fft_thetab(gsa%btc,1)
         call fft_thetab(gsa%bpc,1)
#endif
         if ( this%lDeriv ) then
#ifdef WITH_SHTNS
            gsa%cbrc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_cbr((nThetaStart-1)*n_phi_max+1:)
            gsa%cbtc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_cbt((nThetaStart-1)*n_phi_max+1:)
            gsa%cbpc(1:n_phi_max, 1:nfs) => &
                this%leg_helper%shtns_cbp((nThetaStart-1)*n_phi_max+1:)
#else
            call fft_thetab(gsa%cbrc,1)
            call fft_thetab(gsa%cbtc,1)
            call fft_thetab(gsa%cbpc,1)
#endif
         end if
      end if

   end subroutine transform_to_grid_space
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space(this,nThetaStart,nThetaStop,gsa,nl_lm)

      class(rIterThetaBlocking_t) :: this
      integer,intent(in) :: nThetaStart,nThetaStop
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm
      
      ! Local variables
      integer :: nTheta,nPhi
  
      if ( (.not.this%isRadialBoundaryPoint) .and. ( l_conv_nl .or. l_mag_LF ) ) then
         !PERFON('inner1')
         if ( l_conv_nl .and. l_mag_LF ) then
            if ( this%nR>n_r_LCR ) then
               do nTheta=1,this%sizeThetaB
                  do nPhi=1,nrp
                     gsa%Advr(nPhi,nTheta)=gsa%Advr(nPhi,nTheta) + gsa%LFr(nPhi,nTheta)
                     gsa%Advt(nPhi,nTheta)=gsa%Advt(nPhi,nTheta) + gsa%LFt(nPhi,nTheta)
                     gsa%Advp(nPhi,nTheta)=gsa%Advp(nPhi,nTheta) + gsa%LFp(nPhi,nTheta)
                  end do
               end do
            end if
         else if ( l_mag_LF ) then
            if ( this%nR>n_r_LCR ) then
               do nTheta=1,this%sizeThetaB
                  do nPhi=1,nrp
                     gsa%Advr(nPhi,nTheta)=gsa%LFr(nPhi,nTheta)
                     gsa%Advt(nPhi,nTheta)=gsa%LFt(nPhi,nTheta)
                     gsa%Advp(nPhi,nTheta)=gsa%LFp(nPhi,nTheta)
                  end do
               end do
            else
               do nTheta=1,this%sizeThetaB
                  do nPhi=1,nrp
                     gsa%Advr(nPhi,nTheta)=0.0_cp
                     gsa%Advt(nPhi,nTheta)=0.0_cp
                     gsa%Advp(nPhi,nTheta)=0.0_cp
                  end do
               end do
            end if
         end if
         call fft_thetab(gsa%Advr,-1)
         call fft_thetab(gsa%Advt,-1)
         call fft_thetab(gsa%Advp,-1)
         call legTF3(nThetaStart,nl_lm%AdvrLM,nl_lm%AdvtLM,nl_lm%AdvpLM,    &
              &      gsa%Advr,gsa%Advt,gsa%Advp)
         if ( this%lRmsCalc .and. l_mag_LF .and. this%nR>n_r_LCR ) then 
            ! LF treated extra:
            call fft_thetab(gsa%LFr,-1)
            call fft_thetab(gsa%LFt,-1)
            call fft_thetab(gsa%LFp,-1)
            call legTF3(nThetaStart,nl_lm%LFrLM,nl_lm%LFtLM,nl_lm%LFpLM,    &
                 &      gsa%LFr,gsa%LFt,gsa%LFp)
         end if
         !PERFOFF
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_heat ) then
         !PERFON('inner2')
         call fft_thetab(gsa%VSr,-1)
         call fft_thetab(gsa%VSt,-1)
         call fft_thetab(gsa%VSp,-1)
         call legTF3(nThetaStart,nl_lm%VSrLM,nl_lm%VStLM,nl_lm%VSpLM,       &
              &      gsa%VSr,gsa%VSt,gsa%VSp)
         if (l_anel) then ! anelastic stuff 
            if ( l_mag_nl .and. this%nR>n_r_LCR ) then
               call fft_thetab(gsa%ViscHeat,-1)
               call fft_thetab(gsa%OhmLoss,-1)
               call legTF2(nThetaStart,nl_lm%OhmLossLM,nl_lm%ViscHeatLM,    &
                    &      gsa%OhmLoss,gsa%ViscHeat)
            else
               call fft_thetab(gsa%ViscHeat,-1)
               call legTF1(nThetaStart,nl_lm%ViscHeatLM,gsa%ViscHeat)
            end if
         end if
         !PERFOFF
      end if
      if ( l_mag_nl ) then
         !PERFON('mag_nl')
         if ( .not.this%isRadialBoundaryPoint .and. this%nR>n_r_LCR ) then
            call fft_thetab(gsa%VxBr,-1)
            call fft_thetab(gsa%VxBt,-1)
            call fft_thetab(gsa%VxBp,-1)
            call legTF3(nThetaStart,nl_lm%VxBrLM,nl_lm%VxBtLM,nl_lm%VxBpLM, &
                 &       gsa%VxBr,gsa%VxBt,gsa%VxBp)
         else
            !write(*,"(I4,A,ES20.13)") this%nR,", VxBt = ",sum(VxBt*VxBt)
            call fft_thetab(gsa%VxBt,-1)
            call fft_thetab(gsa%VxBp,-1)
            call legTF2(nThetaStart,nl_lm%VxBtLM,nl_lm%VxBpLM,              &
                 &      gsa%VxBt,gsa%VxBp)
         end if
         !PERFOFF
      end if

   end subroutine transform_to_lm_space
!-------------------------------------------------------------------------------
end module rIterThetaBlocking_mod
