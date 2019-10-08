#include "perflib_preproc.cpp"
module rIterThetaBlocking_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use rIteration_mod, only: rIteration_t
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max,lmP_max,nrp,l_max,lmP_max_dtB,   &
       &                 n_phi_maxStr,n_theta_maxStr,n_r_maxStr, &
       &                 lm_maxMag,l_axi
   use blocking, only: nfs
   use logic, only: l_mag,l_conv,l_mag_kin,l_heat,l_HT,l_anel,l_mag_LF,    &
       &            l_conv_nl, l_mag_nl, l_b_nl_cmb, l_b_nl_icb, l_rot_ic, &
       &            l_cond_ic, l_rot_ma, l_cond_ma, l_dtB, l_store_frame,  &
       &            l_movie_oc, l_chemical_conv, l_precession,             &
       &            l_centrifuge, l_TO, l_adv_curl
   use radial_data,only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: or2, orho1
   use fft
   use legendre_spec_to_grid, only: leg_scal_to_grad_spat, leg_scal_to_spat,    &
       &                            leg_polsphtor_to_spat, leg_pol_to_grad_spat,&
       &                            leg_dphi_vec
   use leg_helper_mod, only: leg_helper_t
   use fields, only: s_Rloc, ds_Rloc, xi_Rloc, p_Rloc
   use nonlinear_lm_mod, only:nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use physical_parameters, only: kbots,ktops,n_r_LCR
   use nonlinear_bcs, only: v_rigid_boundary
   use legendre_grid_to_spec

   implicit none

   private

   type, public, abstract, extends(rIteration_t) :: rIterThetaBlocking_t
      ! or with len parameters for the theta block size and number
      !type,public,extends(rIteration_t) :: rIterThetaBlocking_t(sizeThetaB,nThetaBs)
      !integer,len :: sizeThetaB,nThetaBs
      integer :: sizeThetaB, nThetaBs

      !type(nonlinear_lm_t) :: nl_lm
      type(leg_helper_t) :: leg_helper

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

      !----- Help arrays for Legendre transform calculated in legPrepG:
      !      Parallelizatio note: these are the R-distributed versions
      !      of the field scalars.
      call this%leg_helper%initialize(lm_max,lm_maxMag,l_max)

      if ( l_TO ) then
         allocate( this%BsLast(n_phi_maxStr,n_theta_maxStr,nRstart:nRstop) )
         allocate( this%BpLast(n_phi_maxStr,n_theta_maxStr,nRstart:nRstop) )
         allocate( this%BzLast(n_phi_maxStr,n_theta_maxStr,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+ &
         &                3*n_phi_maxStr*n_theta_maxStr*(nRstop-nRstart+1)*&
         &                SIZEOF_DEF_REAL
      end if

   end subroutine allocate_common_arrays
!-------------------------------------------------------------------------------
   subroutine deallocate_common_arrays(this)

      class(rIterThetaBlocking_t) :: this

      call this%leg_helper%finalize()
      if ( l_TO ) deallocate( this%BsLast,this%BpLast,this%BzLast )

   end subroutine deallocate_common_arrays
!-------------------------------------------------------------------------------
   subroutine set_ThetaBlocking(this,nThetaBs,sizeThetaB)

      class(rIterThetaBlocking_t) :: this
      integer,intent(in) :: nThetaBs, sizeThetaB

      this%nThetaBs = nThetaBs

      this%sizeThetaB = sizeThetaB

   end subroutine set_ThetaBlocking
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space(this,nThetaStart,nThetaStop,gsa,time)

      class(rIterThetaBlocking_t), target :: this
      integer, intent(in) :: nThetaStart,nThetaStop
      type(grid_space_arrays_t) :: gsa
      real(cp), intent(in) :: time

      ! Local variables
      logical :: l_calc

      l_calc = (this%nBc == 0) .or. this%lDeriv 

      !----- Legendre transform from (r,l,m) to (r,theta,m):
      call leg_polsphtor_to_spat(l_calc,nThetaStart,this%leg_helper%dLhw,  &
           &                     this%leg_helper%vhG,this%leg_helper%vhC,  &
           &                     gsa%vrc,gsa%vtc,gsa%vpc)
      if ( l_adv_curl ) then
         if ( this%lDeriv ) then
            call leg_polsphtor_to_spat(l_calc,nThetaStart,this%leg_helper%dLhz,   &
                 &                     this%leg_helper%cvhG,this%leg_helper%cvhC, &
                 &                     gsa%cvrc,gsa%cvtc,gsa%cvpc)
         end if

         !-- For some outputs one still need the other terms
         if ( this%lViscBcCalc .or. this%lPowerCalc .or.  &
         &    this%lFluxProfCalc .or. this%lTOCalc .or.   &
         &    ( this%l_frame .and. l_movie_oc .and. l_store_frame) ) then
            call leg_pol_to_grad_spat(l_calc, nThetaStart, this%leg_helper%dLhw, &
                 &                    gsa%dvrdtc)
            call leg_dphi_vec(l_calc, nThetaStart, gsa%vrc, gsa%vtc, gsa%vpc, &
                 &            gsa%dvrdpc, gsa%dvtdpc, gsa%dvpdpc)
            call leg_polsphtor_to_spat(l_calc,nThetaStart,this%leg_helper%dLhdw, &
                 &                     this%leg_helper%dvhdrG,                   &
                 &                     this%leg_helper%dvhdrC,gsa%dvrdrc,        &
                 &                     gsa%dvtdrc,gsa%dvpdrc)
         end if

      else

         if ( this%lDeriv ) then
            call leg_scal_to_spat(nThetaStart, this%leg_helper%dLhz, gsa%cvrc)
            call leg_pol_to_grad_spat(l_calc, nThetaStart, this%leg_helper%dLhw, &
                 &                    gsa%dvrdtc)
            call leg_dphi_vec(l_calc, nThetaStart, gsa%vrc, gsa%vtc, gsa%vpc, &
                 &            gsa%dvrdpc, gsa%dvtdpc, gsa%dvpdpc)
            call leg_polsphtor_to_spat(l_calc,nThetaStart,this%leg_helper%dLhdw, &
                 &                     this%leg_helper%dvhdrG,                   &
                 &                     this%leg_helper%dvhdrC,gsa%dvrdrc,        &
                 &                     gsa%dvtdrc,gsa%dvpdrc)
         end if

      end if

      if ( l_mag ) then
         call leg_polsphtor_to_spat(l_calc,nThetaStart,this%leg_helper%dLhb,   &
              &                     this%leg_helper%bhG,this%leg_helper%bhC,   &
              &                     gsa%brc,gsa%btc,gsa%bpc)
         if ( this%lDeriv ) then
            call leg_polsphtor_to_spat(l_calc,nThetaStart,this%leg_helper%dLhj,   &
                 &                     this%leg_helper%cbhG,this%leg_helper%cbhC, &
                 &                     gsa%cbrc,gsa%cbtc,gsa%cbpc)
         end if
      end if

      if ( l_heat ) then
         call leg_scal_to_spat(nThetaStart, s_Rloc(:,this%nR), gsa%sc)
      end if

      if ( l_chemical_conv ) then
         call leg_scal_to_spat(nThetaStart, xi_Rloc(:,this%nR), gsa%xic)
      end if

      if ( this%lPressCalc ) then
         call leg_scal_to_spat(nThetaStart, p_Rloc(:,this%nR), gsa%pc)
      end if

      if ( l_HT .or. this%lViscBcCalc ) then
         call leg_scal_to_spat(nThetaStart, ds_Rloc(:,this%nR), gsa%drSc)
      end if

      if ( this%lViscBcCalc ) then
         call leg_scal_to_grad_spat(nThetaStart, s_Rloc(:,this%nR), &
              &                     gsa%dsdtc, gsa%dsdpc)
      end if

      if ( this%lRmsCalc ) then
         call leg_scal_to_grad_spat(nThetaStart, p_Rloc(:,this%nR), &
              &                     gsa%dpdtc, gsa%dpdpc)
      end if

      !------ Fourier transform from (r,theta,m) to (r,theta,phi):
      if ( l_conv .or. l_mag_kin ) then
         if ( l_heat ) then
            if ( .not. l_axi ) call fft_thetab(gsa%sc,1)
            if ( this%lViscBcCalc ) then
               if ( .not. l_axi ) then
                  call fft_thetab(gsa%dsdtc,1)
                  call fft_thetab(gsa%dsdpc,1)
               end if
               if (this%nR == n_r_cmb .and. ktops==1) then
                  gsa%dsdtc=0.0_cp
                  gsa%dsdpc=0.0_cp
               end if
               if (this%nR == n_r_icb .and. kbots==1) then
                  gsa%dsdtc=0.0_cp
                  gsa%dsdpc=0.0_cp
               end if
            end if
         end if
         if ( l_chemical_conv .and. (.not. l_axi) ) then
            call fft_thetab(gsa%xic,1)
         end if
         if ( this%lPressCalc  .and. (.not. l_axi) ) then
            call fft_thetab(gsa%pc,1)
         end if
         if ( l_HT .or. this%lViscBcCalc ) then
            if ( .not. l_axi ) call fft_thetab(gsa%drSc,1)
         endif
         if ( thIs%lRmsCalc ) then
            if ( .not. l_axi ) then
               call fft_thetab(gsa%dpdtc,1)
               call fft_thetab(gsa%dpdpc,1)
            end if
         end if
         if ( this%nBc == 0 ) then
            if ( .not. l_axi ) then
               call fft_thetab(gsa%vrc,1)
               call fft_thetab(gsa%vtc,1)
               call fft_thetab(gsa%vpc,1)
            end if
            if ( this%lDeriv .and. ( .not. l_axi ) ) then

               if ( l_adv_curl ) then
                  call fft_thetab(gsa%cvrc,1)
                  call fft_thetab(gsa%cvtc,1)
                  call fft_thetab(gsa%cvpc,1)

                  if ( this%lViscBcCalc .or. this%lPowerCalc .or.  &
                  &    this%lFluxProfCalc .or. this%lTOCalc .or.   &
                  &    ( this%l_frame .and. l_movie_oc .and. l_store_frame) ) then
                     call fft_thetab(gsa%dvrdrc,1)
                     call fft_thetab(gsa%dvtdrc,1)
                     call fft_thetab(gsa%dvpdrc,1)
                     call fft_thetab(gsa%dvrdtc,1)
                     call fft_thetab(gsa%dvrdpc,1)
                     call fft_thetab(gsa%dvtdpc,1)
                     call fft_thetab(gsa%dvpdpc,1)
                  end if
               else
                  call fft_thetab(gsa%dvrdrc,1)
                  call fft_thetab(gsa%dvtdrc,1)
                  call fft_thetab(gsa%dvpdrc,1)
                  call fft_thetab(gsa%cvrc,1)
                  call fft_thetab(gsa%dvrdtc,1)
                  call fft_thetab(gsa%dvrdpc,1)
                  call fft_thetab(gsa%dvtdpc,1)
                  call fft_thetab(gsa%dvpdpc,1)
               end if
            end if
         else if ( this%nBc == 1 ) then ! Stress free
            gsa%vrc = 0.0_cp
            if ( .not. l_axi ) then
               call fft_thetab(gsa%vtc,1)
               call fft_thetab(gsa%vpc,1)
            end if
            if ( this%lDeriv ) then
               gsa%dvrdtc = 0.0_cp
               gsa%dvrdpc = 0.0_cp
               if ( .not. l_axi ) then
                  if ( l_adv_curl ) then
                     call fft_thetab(gsa%cvrc,1)
                     call fft_thetab(gsa%cvtc,1)
                     call fft_thetab(gsa%cvpc,1)
                     if ( this%lViscBcCalc .or. this%lPowerCalc .or.  &
                     &    this%lFluxProfCalc .or. this%lTOCalc .or.   &
                     &    ( this%l_frame .and. l_movie_oc .and. l_store_frame) ) then
                        call fft_thetab(gsa%dvrdrc,1)
                        call fft_thetab(gsa%dvtdrc,1)
                        call fft_thetab(gsa%dvpdrc,1)
                        call fft_thetab(gsa%dvtdpc,1)
                        call fft_thetab(gsa%dvpdpc,1)
                     end if
                  else
                     call fft_thetab(gsa%dvrdrc,1)
                     call fft_thetab(gsa%dvtdrc,1)
                     call fft_thetab(gsa%dvpdrc,1)
                     call fft_thetab(gsa%cvrc,1)
                     call fft_thetab(gsa%dvtdpc,1)
                     call fft_thetab(gsa%dvpdpc,1)
                  end if
               end if
            end if
         else if ( this%nBc == 2 ) then
            if ( this%nR == n_r_cmb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaMA,this%lDeriv, &
                    &                gsa%vrc,gsa%vtc,gsa%vpc,gsa%cvrc,gsa%dvrdtc, &
                    &                gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,            &
                    &                nThetaStart,time)
            else if ( this%nR == n_r_icb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaIC,this%lDeriv, &
                    &                gsa%vrc,gsa%vtc,gsa%vpc,gsa%cvrc,gsa%dvrdtc, &
                    &                gsa%dvrdpc,gsa%dvtdpc,gsa%dvpdpc,            &
                    &                nThetaStart,time)
            end if
            if ( this%lDeriv .and. ( .not. l_axi ) ) then
               call fft_thetab(gsa%dvrdrc,1)
               call fft_thetab(gsa%dvtdrc,1)
               call fft_thetab(gsa%dvpdrc,1)
            end if
         end if
      end if
      if ( l_mag .or. l_mag_LF ) then
         if ( .not. l_axi ) then
            call fft_thetab(gsa%brc,1)
            call fft_thetab(gsa%btc,1)
            call fft_thetab(gsa%bpc,1)
         end if
         if ( this%lDeriv .and. ( .not. l_axi ) ) then
            call fft_thetab(gsa%cbrc,1)
            call fft_thetab(gsa%cbtc,1)
            call fft_thetab(gsa%cbpc,1)
         end if
      end if

   end subroutine transform_to_grid_space
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space(this,nThetaStart,nThetaStop,gsa,nl_lm)

      class(rIterThetaBlocking_t) :: this
      integer,intent(in) :: nThetaStart, nThetaStop
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

      ! Local variables
      integer :: nTheta,nPhi

      if ( (.not.this%isRadialBoundaryPoint .or. this%lRmsCalc) .and. &
            ( l_conv_nl .or. l_mag_LF ) ) then
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

         if ( l_precession ) then
            do nTheta=1,this%sizeThetaB
               do nPhi=1,nrp
                  gsa%Advr(nPhi,nTheta)=gsa%Advr(nPhi,nTheta) + gsa%PCr(nPhi,nTheta)
                  gsa%Advt(nPhi,nTheta)=gsa%Advt(nPhi,nTheta) + gsa%PCt(nPhi,nTheta)
                  gsa%Advp(nPhi,nTheta)=gsa%Advp(nPhi,nTheta) + gsa%PCp(nPhi,nTheta)
               end do
            end do
         end if

         if ( l_centrifuge ) then
            do nTheta=1,this%sizeThetaB
               do nPhi=1,nrp
                  gsa%Advr(nPhi,nTheta)=gsa%Advr(nPhi,nTheta) + gsa%CAr(nPhi,nTheta)
                  gsa%Advt(nPhi,nTheta)=gsa%Advt(nPhi,nTheta) + gsa%CAt(nPhi,nTheta)
               end do
            end do
         end if ! centrifuge

         if ( .not. l_axi ) then
            call fft_thetab(gsa%Advr,-1)
            call fft_thetab(gsa%Advt,-1)
            call fft_thetab(gsa%Advp,-1)
         end if
         call legTF3(nThetaStart,nl_lm%AdvrLM,nl_lm%AdvtLM,nl_lm%AdvpLM,    &
              &      gsa%Advr,gsa%Advt,gsa%Advp)
         if ( this%lRmsCalc .and. l_mag_LF .and. this%nR>n_r_LCR ) then
            ! LF treated extra:
            if ( .not. l_axi ) then
               call fft_thetab(gsa%LFr,-1)
               call fft_thetab(gsa%LFt,-1)
               call fft_thetab(gsa%LFp,-1)
            end if
            call legTF3(nThetaStart,nl_lm%LFrLM,nl_lm%LFtLM,nl_lm%LFpLM,    &
                 &      gsa%LFr,gsa%LFt,gsa%LFp)
         end if
         !PERFOFF
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_heat ) then
         !PERFON('inner2')
         if ( .not. l_axi ) then
            call fft_thetab(gsa%VSr,-1)
            call fft_thetab(gsa%VSt,-1)
            call fft_thetab(gsa%VSp,-1)
         end if
         call legTF3(nThetaStart,nl_lm%VSrLM,nl_lm%VStLM,nl_lm%VSpLM,       &
              &      gsa%VSr,gsa%VSt,gsa%VSp)
         if (l_anel) then ! anelastic stuff
            if ( l_mag_nl .and. this%nR>n_r_LCR ) then
               if ( .not. l_axi ) then
                  call fft_thetab(gsa%ViscHeat,-1)
                  call fft_thetab(gsa%OhmLoss,-1)
               end if
               call legTF2(nThetaStart,nl_lm%OhmLossLM,nl_lm%ViscHeatLM,    &
                    &      gsa%OhmLoss,gsa%ViscHeat)
            else
               if ( .not. l_axi ) call fft_thetab(gsa%ViscHeat,-1)
               call legTF1(nThetaStart,nl_lm%ViscHeatLM,gsa%ViscHeat)
            end if
         end if
         !PERFOFF
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_chemical_conv ) then
         if ( .not. l_axi ) then
            call fft_thetab(gsa%VXir,-1)
            call fft_thetab(gsa%VXit,-1)
            call fft_thetab(gsa%VXip,-1)
         end if
         call legTF3(nThetaStart,nl_lm%VXirLM,nl_lm%VXitLM,nl_lm%VXipLM,    &
              &      gsa%VXir,gsa%VXit,gsa%VXip)
      end if
      if ( l_mag_nl ) then
         !PERFON('mag_nl')
         if ( .not.this%isRadialBoundaryPoint .and. this%nR>n_r_LCR ) then
            if ( .not. l_axi ) then
               call fft_thetab(gsa%VxBr,-1)
               call fft_thetab(gsa%VxBt,-1)
               call fft_thetab(gsa%VxBp,-1)
            end if
            call legTF3(nThetaStart,nl_lm%VxBrLM,nl_lm%VxBtLM,nl_lm%VxBpLM, &
                 &       gsa%VxBr,gsa%VxBt,gsa%VxBp)
         else
            if ( .not. l_axi ) then
               call fft_thetab(gsa%VxBt,-1)
               call fft_thetab(gsa%VxBp,-1)
            end if
            call legTF2(nThetaStart,nl_lm%VxBtLM,nl_lm%VxBpLM,   &
                 &      gsa%VxBt,gsa%VxBp)
         end if
         !PERFOFF
      end if

      if ( this%lRmsCalc ) then
         if ( .not. l_axi ) then
            call fft_thetab(gsa%dpdtc,-1)
            call fft_thetab(gsa%dpdpc,-1)
         end if
         call legTF_spher_tor(nThetaStart,nl_lm%PFp2LM,nl_lm%PFt2LM,gsa%dpdpc, &
              &               gsa%dpdtc)
         if ( .not. l_axi ) then
            call fft_thetab(gsa%CFt2,-1)
            call fft_thetab(gsa%CFp2,-1)
         end if
         call legTF_spher_tor(nThetaStart,nl_lm%CFp2LM,nl_lm%CFt2LM, &
              &               gsa%CFp2,gsa%CFt2)
         if ( l_conv_nl ) then
            if ( .not. l_axi ) then
               call fft_thetab(gsa%Advt2,-1)
               call fft_thetab(gsa%Advp2,-1)
            end if
            call legTF_spher_tor(nThetaStart,nl_lm%Advp2LM,nl_lm%Advt2LM, &
                 &               gsa%Advp2,gsa%Advt2)
         end if
         if ( l_mag_nl .and. this%nR>n_r_LCR ) then
            if ( .not. l_axi ) then
               call fft_thetab(gsa%LFt2,-1)
               call fft_thetab(gsa%LFp2,-1)
            end if
            call legTF_spher_tor(nThetaStart,nl_lm%LFp2LM,nl_lm%LFt2LM, &
                 &               gsa%LFp2,gsa%LFt2)
         end if
         if ( .not. l_axi ) then
            call fft_thetab(gsa%dtVr,-1)
            call fft_thetab(gsa%dtVt,-1)
            call fft_thetab(gsa%dtVp,-1)
         end if
         call legTF1(nThetaStart,nl_lm%dtVrLM,gsa%dtVr)
         call legTF_spher_tor(nThetaStart,nl_lm%dtVpLM,nl_lm%dtVtLM, &
              &               gsa%dtVp,gsa%dtVt)
      end if

   end subroutine transform_to_lm_space
!-------------------------------------------------------------------------------
end module rIterThetaBlocking_mod
