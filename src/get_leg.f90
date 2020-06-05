module hybrid_space_mod

   use precision_mod
   use nonlinear_3D_lm_mod, only: nonlinear_3D_lm_t
   use constants, only: zero
   use truncation, only: n_theta_max, n_m_loc, nRstart, nRstop, n_r_loc,     &
       &                 n_r_cmb, n_r_icb, n_m_max, n_theta_loc, nThetaStart,&
       &                 nThetaStop, n_lm_loc, n_phi_max
   use logic, only: l_mag, l_adv_curl, l_chemical_conv, l_viscBcCalc, l_RMS,   &
       &            l_conv, l_mag_kin, l_heat, l_HT, l_movie_oc, l_store_frame,&
       &            l_mag_LF, l_anel, l_mag_nl, l_conv_nl
   use mem_alloc, only: bytes_allocated
   use physical_parameters, only: ktops, kbots, n_r_LCR, ktopv, kbotv
   use radial_functions, only: l_R, or2
   use shtns, only: scal_to_hyb, scal_to_grad_hyb, torpol_to_hyb,          &
       &            torpol_to_curl_hyb, pol_to_grad_hyb, torpol_to_dphhyb, &
       &            pol_to_curlr_hyb, hyb_to_SH, hyb_to_qst, hyb_to_sphertor
   use mpi_thetap_mod, only: transpose_theta_m_many, transpose_m_theta_many

   implicit none

   private

   type, public :: hybrid_3D_arrays_t
      complex(cp), pointer :: vr_Mloc(:,:,:), vt_Mloc(:,:,:), vp_Mloc(:,:,:)
      complex(cp), pointer :: dvrdr_Mloc(:,:,:), dvtdr_Mloc(:,:,:), dvpdr_Mloc(:,:,:)
      complex(cp), pointer :: cvr_Mloc(:,:,:), s_Mloc(:,:,:), dsdr_Mloc(:,:,:)
      complex(cp), pointer :: dvrdt_Mloc(:,:,:), dvrdp_Mloc(:,:,:)
      complex(cp), pointer :: dvtdp_Mloc(:,:,:), dvpdp_Mloc(:,:,:)
      complex(cp), pointer :: br_Mloc(:,:,:), bt_Mloc(:,:,:), bp_Mloc(:,:,:)
      complex(cp), pointer :: cbr_Mloc(:,:,:), cbt_Mloc(:,:,:), cbp_Mloc(:,:,:)
      complex(cp), pointer :: p_Mloc(:,:,:), xi_Mloc(:,:,:), cvt_Mloc(:,:,:), cvp_Mloc(:,:,:)
      complex(cp), pointer :: dsdt_Mloc(:,:,:), dsdp_Mloc(:,:,:)
      complex(cp), pointer :: dpdt_Mloc(:,:,:), dpdp_Mloc(:,:,:)

      complex(cp), pointer :: vel_pThloc(:,:,:,:), gradvel_pThloc(:,:,:,:)
      complex(cp), pointer :: s_pThloc(:,:,:), grads_pThloc(:,:,:,:)
      complex(cp), pointer :: p_pThloc(:,:,:), gradp_pThloc(:,:,:,:)
      complex(cp), pointer :: xi_pThloc(:,:,:), mag_pThloc(:,:,:,:)
      complex(cp), pointer :: vr_Thloc(:,:,:), vt_Thloc(:,:,:), vp_Thloc(:,:,:)
      complex(cp), pointer :: dvrdr_Thloc(:,:,:), dvtdr_Thloc(:,:,:), dvpdr_Thloc(:,:,:)
      complex(cp), pointer :: cvr_Thloc(:,:,:), s_Thloc(:,:,:), dsdr_Thloc(:,:,:)
      complex(cp), pointer :: dvrdt_Thloc(:,:,:), dvrdp_Thloc(:,:,:)
      complex(cp), pointer :: dvtdp_Thloc(:,:,:), dvpdp_Thloc(:,:,:)
      complex(cp), pointer :: br_Thloc(:,:,:), bt_Thloc(:,:,:), bp_Thloc(:,:,:)
      complex(cp), pointer :: cbr_Thloc(:,:,:), cbt_Thloc(:,:,:), cbp_Thloc(:,:,:)
      complex(cp), pointer :: p_Thloc(:,:,:), xi_Thloc(:,:,:), cvt_Thloc(:,:,:), cvp_Thloc(:,:,:)
      complex(cp), pointer :: dsdt_Thloc(:,:,:), dsdp_Thloc(:,:,:)
      complex(cp), pointer :: dpdt_Thloc(:,:,:), dpdp_Thloc(:,:,:)

      complex(cp), pointer :: Advr_Mloc(:,:,:), Advt_Mloc(:,:,:), Advp_Mloc(:,:,:)
      complex(cp), pointer :: LFr_Mloc(:,:,:), LFt_Mloc(:,:,:), LFp_Mloc(:,:,:)
      complex(cp), pointer :: VSr_Mloc(:,:,:), VSt_Mloc(:,:,:), VSp_Mloc(:,:,:)
      complex(cp), pointer :: VXir_Mloc(:,:,:), VXit_Mloc(:,:,:), VXip_Mloc(:,:,:)
      complex(cp), pointer :: VxBr_Mloc(:,:,:), VxBt_Mloc(:,:,:), VxBp_Mloc(:,:,:)
      complex(cp), pointer :: ViscHeat_Mloc(:,:,:), OhmLoss_Mloc(:,:,:)
      complex(cp), pointer :: dpkindr_Mloc(:,:,:)
      complex(cp), pointer :: dtVr_Mloc(:,:,:), dtVt_Mloc(:,:,:), dtVp_Mloc(:,:,:)
      complex(cp), pointer :: PFt2_Mloc(:,:,:), PFp2_Mloc(:,:,:)
      complex(cp), pointer :: CFt2_Mloc(:,:,:), CFp2_Mloc(:,:,:)
      complex(cp), pointer :: Advt2_Mloc(:,:,:), Advp2_Mloc(:,:,:)
      complex(cp), pointer :: LFt2_Mloc(:,:,:), LFp2_Mloc(:,:,:)

      complex(cp), pointer :: NSadv_pThloc(:,:,:,:), heatadv_pThloc(:,:,:,:)
      complex(cp), pointer :: compadv_pThloc(:,:,:,:), emf_pThloc(:,:,:,:)
      complex(cp), pointer :: RMS_pThloc(:,:,:,:), LF_pThloc(:,:,:,:)
      complex(cp), pointer :: LF2_pThloc(:,:,:,:), PF2_pThloc(:,:,:,:)
      complex(cp), pointer :: dtV_pThloc(:,:,:,:), anel_pThloc(:,:,:,:)
      complex(cp), pointer :: Advr_Thloc(:,:,:), Advt_Thloc(:,:,:), Advp_Thloc(:,:,:)
      complex(cp), pointer :: LFr_Thloc(:,:,:), LFt_Thloc(:,:,:), LFp_Thloc(:,:,:)
      complex(cp), pointer :: VSr_Thloc(:,:,:), VSt_Thloc(:,:,:), VSp_Thloc(:,:,:)
      complex(cp), pointer :: VXir_Thloc(:,:,:), VXit_Thloc(:,:,:), VXip_Thloc(:,:,:)
      complex(cp), pointer :: VxBr_Thloc(:,:,:), VxBt_Thloc(:,:,:), VxBp_Thloc(:,:,:)
      complex(cp), pointer :: ViscHeat_Thloc(:,:,:), OhmLoss_Thloc(:,:,:)
      complex(cp), pointer :: dpkindr_Thloc(:,:,:)
      complex(cp), pointer :: dtVr_Thloc(:,:,:), dtVt_Thloc(:,:,:), dtVp_Thloc(:,:,:)
      complex(cp), pointer :: PFt2_Thloc(:,:,:), PFp2_Thloc(:,:,:)
      complex(cp), pointer :: CFt2_Thloc(:,:,:), CFp2_Thloc(:,:,:)
      complex(cp), pointer :: Advt2_Thloc(:,:,:), Advp2_Thloc(:,:,:)
      complex(cp), pointer :: LFt2_Thloc(:,:,:), LFp2_Thloc(:,:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: leg_spec_to_hyb
      procedure :: leg_hyb_to_spec
      procedure :: transp_Mloc_to_Thloc
      procedure :: transp_Thloc_to_Mloc
   end type hybrid_3D_arrays_t

contains

   subroutine initialize(this)
      !
      ! Memory allocation for arrays in hybrid space, i.e. (n_theta, n_m, n_r)
      !

      class(hybrid_3D_arrays_t) :: this


      allocate( this%vr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%vt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%vp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvrdr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvtdr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvpdr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%cvr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvrdt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvrdp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvtdp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dvpdp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%p_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%s_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( this%dsdr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+14*n_theta_max*n_m_loc*n_r_loc*&
      &                 SIZEOF_DEF_COMPLEX


      if ( l_adv_curl ) then
         allocate( this%vel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop, 6) )
      else
         allocate( this%vel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop, 4) )
      end if
      this%vr_Thloc(1:,nThetaStart:,nRstart:) => &
      &          this%vel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%vt_Thloc(1:,nThetaStart:,nRstart:) => &
      &          this%vel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
      this%vp_Thloc(1:,nThetaStart:,nRstart:) => &
      &          this%vel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
      this%cvr_Thloc(1:,nThetaStart:,nRstart:) => &
      &          this%vel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,4)
      if ( l_adv_curl ) then
         this%cvt_Thloc(1:,nThetaStart:,nRstart:) => &
         &          this%vel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,5)
         this%cvp_Thloc(1:,nThetaStart:,nRstart:) => &
         &          this%vel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,6)
      end if
      this%vel_pThloc(:,:,:,:)=zero


      allocate( this%p_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop) )
      this%p_Thloc(1:,nThetaStart:,nRstart:) => &
      &          this%p_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop)
      allocate( this%s_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop) )
      this%s_Thloc(1:,nThetaStart:,nRstart:) => &
      &          this%s_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop)
      this%p_pThloc(:,:,:)=zero
      this%s_pThloc(:,:,:)=zero

      allocate( this%gradvel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop, 7) )
      this%dvrdr_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%dvtdr_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
      this%dvpdr_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
      this%dvrdp_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,4)
      this%dvtdp_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,5)
      this%dvpdp_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,6)
      this%dvrdt_Thloc(1:,nThetaStart:,nRstart:) => &
      &        this%gradvel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,7)
      this%gradvel_pThloc(:,:,:,:)=zero
      bytes_allocated = bytes_allocated+14*n_theta_loc*n_m_max*n_r_loc*&
      &                 SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( this%xi_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+1*n_theta_max*n_m_loc*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         allocate( this%xi_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop) )
         this%xi_Thloc(1:,nThetaStart:,nRstart:) => &
         &            this%xi_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop)
         this%xi_pThloc(:,:,:)=zero
         bytes_allocated = bytes_allocated+1*n_theta_loc*n_m_max*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
      end if

      if ( l_mag ) then
         allocate( this%br_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%bt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%bp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%cbr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%cbt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%cbp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+6*n_theta_max*n_m_loc*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         allocate( this%mag_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,6) )
         bytes_allocated = bytes_allocated+6*n_theta_loc*n_m_max*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         this%br_Thloc(1:,nThetaStart:,nRstart:) => &
         &         this%mag_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%bt_Thloc(1:,nThetaStart:,nRstart:) => &
         &         this%mag_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%bp_Thloc(1:,nThetaStart:,nRstart:) => &
         &         this%mag_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
         this%cbr_Thloc(1:,nThetaStart:,nRstart:) => &
         &         this%mag_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,4)
         this%cbt_Thloc(1:,nThetaStart:,nRstart:) => &
         &         this%mag_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,5)
         this%cbp_Thloc(1:,nThetaStart:,nRstart:) => &
         &         this%mag_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,6)
         this%mag_pThloc(:,:,:,:)=zero
      end if

      if ( l_adv_curl ) then
         allocate( this%cvt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%cvp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+2*n_theta_max*n_m_loc*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
      end if

      if ( l_viscBcCalc ) then
         allocate( this%dsdt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%dsdp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+2*n_theta_max*n_m_loc*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         allocate( this%grads_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+2*n_theta_loc*n_m_max*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%dsdr_Thloc(1:,nThetaStart:,nRstart:) => &
         &       this%grads_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%dsdt_Thloc(1:,nThetaStart:,nRstart:) => &
         &       this%grads_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%dsdp_Thloc(1:,nThetaStart:,nRstart:) => &
         &       this%grads_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
      else
         allocate( this%grads_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,1) )
         this%dsdr_Thloc(1:,nThetaStart:,nRstart:) => &
         &       this%grads_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
      end if
      this%grads_pThloc(:,:,:,:)=zero

      if ( l_RMS ) then
         allocate( this%dpdt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%dpdp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+2*n_theta_max*n_m_loc*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         allocate( this%gradp_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
         bytes_allocated = bytes_allocated+2*n_theta_loc*n_m_max*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%dpdt_Thloc(1:,nThetaStart:,nRstart:) => &
         &          this%gradp_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%dpdp_Thloc(1:,nThetaStart:,nRstart:) => &
         &          this%gradp_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%gradp_pThloc(:,:,:,:)=zero
      end if

      allocate( this%NSadv_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
      this%NSadv_pThloc(:,:,:,:)=zero
      this%Advr_Thloc(1:,nThetaStart:,nRstart:) => &
      &      this%NSadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%Advt_Thloc(1:,nThetaStart:,nRstart:) => &
      &      this%NSadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
      this%Advp_Thloc(1:,nThetaStart:,nRstart:) => &
      &      this%NSadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)

      if ( l_heat ) then
         allocate( this%heatadv_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%heatadv_pThloc(:,:,:,:)=zero
         this%VSr_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%heatadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%VSt_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%heatadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%VSp_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%heatadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
      end if 

      if ( l_chemical_conv ) then
         allocate( this%compadv_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%compadv_pThloc(:,:,:,:)=zero
         this%VXir_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%compadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%VXit_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%compadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%VXip_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%compadv_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
      end if

      if ( l_anel ) then
         if ( l_mag ) then
            allocate( this%anel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
         else
            allocate( this%anel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,1) )
         end if
         this%ViscHeat_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%anel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         if ( l_mag ) then
            this%OhmLoss_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%anel_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         end if
      end if

      if ( l_mag_nl ) then
         allocate( this%emf_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%emf_pThloc(:,:,:,:)=zero
         this%VxBr_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%emf_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%VxBt_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%emf_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%VxBp_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%emf_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
      end if

      allocate( this%Advr_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
      allocate( this%Advt_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
      allocate( this%Advp_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )

      if ( l_heat ) then
         allocate( this%VSr_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
         allocate( this%VSt_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
         allocate( this%VSp_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXir_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
         allocate( this%VXit_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
         allocate( this%VXip_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
      end if

      if ( l_anel ) then
         if ( l_mag ) then
            allocate( this%OhmLoss_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop))
         end if
         allocate( this%ViscHeat_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
      end if

      if ( l_mag_nl ) then
         allocate( this%VxBr_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
         allocate( this%VxBt_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
         allocate( this%VxBp_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop) )
      end if

      if ( l_RMS ) then
         allocate( this%PF2_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
         this%PF2_pThloc(:,:,:,:)=zero
         this%PFt2_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%PF2_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%PFp2_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%PF2_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         if ( l_adv_curl ) then
            allocate( this%RMS_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,5) )
         else
            allocate( this%RMS_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,4) )
         end if
         this%RMS_pThloc(:,:,:,:)=zero
         this%CFt2_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%RMS_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%CFp2_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%RMS_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%Advt2_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%RMS_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
         this%Advp2_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%RMS_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,4)
         if ( l_adv_curl ) then
            this%dpkindr_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%RMS_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,5)
         end if
         allocate( this%dtV_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%dtV_pThloc(:,:,:,:)=zero
         this%dtVr_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%dtV_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%dtVt_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%dtV_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%dtVp_Thloc(1:,nThetaStart:,nRstart:) => &
         &      this%dtV_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)

         if ( l_mag ) then
            allocate( this%LF_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
            allocate( this%LF2_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
            this%LF_pThloc(:,:,:,:)=zero
            this%LF2_pThloc(:,:,:,:)=zero
            this%LFr_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%LF_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
            this%LFt_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%LF_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
            this%LFp_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%LF_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,3)
            this%LFt2_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%LF2_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,1)
            this%LFp2_Thloc(1:,nThetaStart:,nRstart:) => &
            &      this%LF2_pThloc(1:n_m_max,nThetaStart:nThetaStop,nRstart:nRstop,2)
         end if

         allocate( this%PFt2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%PFp2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%CFt2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%CFp2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%Advt2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%Advp2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         if ( l_adv_curl ) then
            allocate( this%dpkindr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         end if

         allocate( this%dtVr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%dtVt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         allocate( this%dtVp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         if ( l_mag ) then
            allocate( this%LFr_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
            allocate( this%LFt_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
            allocate( this%LFp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
            allocate( this%LFt2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
            allocate( this%LFp2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         end if
         bytes_allocated = bytes_allocated+5*n_m_loc*n_theta_max*n_r_loc* &
         &                 SIZEOF_DEF_COMPLEX

         bytes_allocated = bytes_allocated+5*n_m_max*n_theta_loc*n_r_loc* &
         &                 SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize
!-----------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation for arrays in hybrid space
      !

      class(hybrid_3D_arrays_t) :: this

      if ( l_RMS ) then
         deallocate( this%dtVr_Mloc, this%dtVt_Mloc, this%dtVp_Mloc ) 
         deallocate( this%LFt2_Mloc, this%LFp2_Mloc)
         deallocate( this%dpdt_Mloc, this%dpdp_Mloc )
         deallocate( this%RMS_pThloc, this%dtV_pThloc )
         if ( l_mag ) deallocate( this%LF_pThloc )
      end if
      if ( l_adv_curl ) then
         deallocate( this%cvt_Mloc, this%cvp_Mloc )
      end if
      if ( l_mag ) then
         deallocate( this%br_Mloc, this%bt_Mloc, this%bp_Mloc, &
         &           this%cbr_Mloc, this%cbt_Mloc, this%cbp_Mloc )
         deallocate( this%mag_pThloc )
      end if
      if ( l_chemical_conv ) then
         deallocate( this%xi_Mloc, this%xi_pThloc )
      end if

      deallocate( this%vr_Mloc, this%vt_Mloc, this%vp_Mloc )
      deallocate( this%dvrdr_Mloc, this%dvtdr_Mloc, this%dvpdr_Mloc )
      deallocate( this%cvr_Mloc, this%dvrdt_Mloc, this%dvrdp_Mloc, this%dvtdp_Mloc )
      deallocate( this%dvpdp_Mloc, this%p_Mloc, this%s_Mloc, this%dsdr_Mloc )

      deallocate( this%vel_pThloc, this%gradvel_pThloc )
      deallocate( this%s_pThloc,this%p_pThloc, this%grads_pThloc)

   end subroutine finalize
!-----------------------------------------------------------------------------------
   subroutine leg_spec_to_hyb(this, w, dw, ddw, z, dz, b, db, ddb, aj, dj, s, ds, &
              &               p, xi, lVisc, lRmsCalc, lPressCalc, lTOCalc,        &
              &               lPowerCalc, lFluxProfCalc, lPerpParCalc, lHelCalc,  &
              &               l_frame)

      class(hybrid_3D_arrays_t) :: this

      !-- Input variables
      complex(cp), intent(in) :: w(n_lm_loc,nRstart:nRstop)   ! Poloidal potential
      complex(cp), intent(in) :: dw(n_lm_loc,nRstart:nRstop)  ! dw/dr
      complex(cp), intent(in) :: ddw(n_lm_loc,nRstart:nRstop) ! d^2w/dr^2
      complex(cp), intent(in) :: z(n_lm_loc,nRstart:nRstop)   ! Toroidal potential
      complex(cp), intent(in) :: dz(n_lm_loc,nRstart:nRstop)  ! dz/dr
      complex(cp), intent(in) :: b(n_lm_loc,nRstart:nRstop)   ! Magnetic poloidal potential
      complex(cp), intent(in) :: db(n_lm_loc,nRstart:nRstop)  ! db/dr
      complex(cp), intent(in) :: ddb(n_lm_loc,nRstart:nRstop) ! d^2b/dr^2
      complex(cp), intent(in) :: aj(n_lm_loc,nRstart:nRstop)  ! Magnetic toroidal potential
      complex(cp), intent(in) :: dj(n_lm_loc,nRstart:nRstop)  ! dj/dr
      complex(cp), intent(in) :: s(n_lm_loc,nRstart:nRstop)   ! Entropy
      complex(cp), intent(in) :: ds(n_lm_loc,nRstart:nRstop)  ! ds/dr
      complex(cp), intent(in) :: p(n_lm_loc,nRstart:nRstop)   ! Pressure
      complex(cp), intent(in) :: xi(n_lm_loc,nRstart:nRstop)  ! Chemical composition
      logical,     intent(in) :: lVisc, lRmsCalc, lPressCalc, lPowerCalc
      logical,     intent(in) :: lTOCalc, lFluxProfCalc, l_frame, lHelCalc
      logical,     intent(in) :: lPerpParCalc

      !-- Local variables
      logical :: lDeriv
      integer :: nR, nBc

      do nR=nRstart,nRstop

         nBc = 0
         lDeriv = .true.
         if ( nR == n_r_cmb ) then
            nBc = ktopv
            lDeriv= lTOCalc .or. lHelCalc .or. l_frame .or. lPerpParCalc   &
            &       .or. lVisc .or. lFluxProfCalc .or. lRmsCalc .or. lPowerCalc
         else if ( nR == n_r_icb ) then
            nBc = kbotv
            lDeriv= lTOCalc .or. lHelCalc .or. l_frame  .or. lPerpParCalc  &
            &       .or. lVisc .or. lFluxProfCalc .or. lRmsCalc .or. lPowerCalc
         end if

         if ( l_conv .or. l_mag_kin ) then
            if ( l_heat ) then
               call scal_to_hyb(s(:,nR), this%s_Mloc(:,:,nR), l_R(nR))
               if ( lVisc ) then
                  call scal_to_grad_hyb(s(:,nR), this%dsdt_Mloc(:,:,nR), &
                       &                this%dsdp_Mloc(:,:,nR),&
                       &                l_R(nR))
                  if ( nR == n_r_cmb .and. ktops==1) then
                     this%dsdt_Mloc(:,:,nR)=zero
                     this%dsdp_Mloc(:,:,nR)=zero
                  end if
                  if ( nR == n_r_icb .and. kbots==1) then
                     this%dsdt_Mloc(:,:,nR)=zero
                     this%dsdp_Mloc(:,:,nR)=zero
                  end if
               end if
            end if

            if ( lRmsCalc ) call scal_to_grad_hyb(p(:,nR), this%dpdt_Mloc(:,:,nR), &
                                 &                this%dpdp_Mloc(:,:,nR), l_R(nR))

            !-- Pressure
            if ( lPressCalc ) call scal_to_hyb(p(:,nR), this%p_Mloc(:,:,nR), l_R(nR))

            !-- Composition
            if ( l_chemical_conv ) call scal_to_hyb(xi(:,nR), this%xi_Mloc(:,:,nR),&
                                        &           l_R(nR))

            if ( l_HT .or. lVisc ) then
               call scal_to_hyb(ds(:,nR), this%dsdr_Mloc(:,:,nR), l_R(nR))
            endif

            if ( nBc == 0 ) then ! Bulk points
               !-- pol, sph, tor > ur,ut,up
               call torpol_to_hyb(w(:,nR), dw(:,nR), z(:,nR), this%vr_Mloc(:,:,nR),&
                    &             this%vt_Mloc(:,:,nR), this%vp_Mloc(:,:,nR), l_R(nR))

               !-- Advection is treated as u \times \curl u
               if ( l_adv_curl ) then
                  !-- z,dz,w,dd< -> wr,wt,wp
                  call torpol_to_curl_hyb(or2(nR), w(:,nR), ddw(:,nR), z(:,nR), &
                       &                  dz(:,nR), this%cvr_Mloc(:,:,nR),      &
                       &                  this%cvt_Mloc(:,:,nR),                &
                       &                  this%cvp_Mloc(:,:,nR), l_R(nR))

                  !-- For some outputs one still need the other terms
                  if ( lVisc .or. lPowerCalc .or. lRmsCalc             &
                  &    .or. lFluxProfCalc .or. lTOCalc .or.            &
                  &    ( l_frame .and. l_movie_oc .and. l_store_frame) ) then

                     call torpol_to_hyb(dw(:,nR), ddw(:,nR), dz(:,nR),    &
                          &             this%dvrdr_Mloc(:,:,nR),          &
                          &             this%dvtdr_Mloc(:,:,nR),          &
                          &             this%dvpdr_Mloc(:,:,nR), l_R(nR))
                     call pol_to_grad_hyb(w(:,nR), this%dvrdt_Mloc(:,:,nR), &
                          &               this%dvrdp_Mloc(:,:,nR), l_R(nR))
                     call torpol_to_dphhyb(dw(:,nR), z(:,nR),          &
                          &                this%dvtdp_Mloc(:,:,nR),    &
                          &                this%dvpdp_Mloc(:,:,nR), l_R(nR))
                  end if

               else ! Advection is treated as u\grad u

                  call torpol_to_hyb(dw(:,nR), ddw(:,nR), dz(:,nR),     &
                       &             this%dvrdr_Mloc(:,:,nR),           &
                       &             this%dvtdr_Mloc(:,:,nR),           &
                       &             this%dvpdr_Mloc(:,:,nR), l_R(nR))
                  call pol_to_curlr_hyb(z(:,nR), this%cvr_Mloc(:,:,nR), l_R(nR))
                  call pol_to_grad_hyb(w(:,nR), this%dvrdt_Mloc(:,:,nR),  &
                       &               this%dvrdp_Mloc(:,:,nR), l_R(nR))
                  call torpol_to_dphhyb(dw(:,nR), z(:,nR), this%dvtdp_Mloc(:,:,nR), &
                       &                this%dvpdp_Mloc(:,:,nR), l_R(nR))
               end if

            else if ( nBc == 1 ) then ! Stress free
                ! TODO don't compute vr_Mloc(:,:,nR) as it is set to 0 afterward
               call torpol_to_hyb(w(:,nR), dw(:,nR), z(:,nR), this%vr_Mloc(:,:,nR),&
                    &             this%vt_Mloc(:,:,nR), this%vp_Mloc(:,:,nR), l_R(nR))

               this%vr_Mloc(:,:,nR) = zero
               if ( lDeriv ) then
                  this%dvrdt_Mloc(:,:,nR) = zero
                  this%dvrdp_Mloc(:,:,nR) = zero
                  call torpol_to_hyb(dw(:,nR), ddw(:,nR), dz(:,nR),   &
                       &             this%dvrdr_Mloc(:,:,nR),         &
                       &             this%dvtdr_Mloc(:,:,nR),         &
                       &             this%dvpdr_Mloc(:,:,nR), l_R(nR))
                  call pol_to_curlr_hyb(z(:,nR), this%cvr_Mloc(:,:,nR), l_R(nR))
                  call torpol_to_dphhyb(dw(:,nR), z(:,nR), this%dvtdp_Mloc(:,:,nR), &
                       &                this%dvpdp_Mloc(:,:,nR), l_R(nR))
               end if
            else if ( nBc == 2 ) then
               if ( lDeriv ) then
                  call torpol_to_hyb(dw(:,nR), ddw(:,nR), dz(:,nR),  &
                       &             this%dvrdr_Mloc(:,:,nR),        &
                       &             this%dvtdr_Mloc(:,:,nR),        &
                       &             this%dvpdr_Mloc(:,:,nR), l_R(nR))
               end if
            end if

         end if

         if ( l_mag .or. l_mag_LF ) then
            call torpol_to_hyb(b(:,nR), db(:,nR), aj(:,nR), this%br_Mloc(:,:,nR), &
                 &             this%bt_Mloc(:,:,nR), this%bp_Mloc(:,:,nR), l_R(nR))
            if ( lDeriv ) then
               call torpol_to_curl_hyb(or2(nR), b(:,nR), ddb(:,nR), aj(:,nR),       &
                    &                  dj(:,nR), this%cbr_Mloc(:,:,nR),             &
                    &                  this%cbt_Mloc(:,:,nR), this%cbp_Mloc(:,:,nR),&
                    &                  l_R(nR))
            end if
         end if

      end do

   end subroutine leg_spec_to_hyb
!-----------------------------------------------------------------------------------
   subroutine transp_Mloc_to_Thloc(this, lVisc, lRmsCalc, lPressCalc, lTOCalc, &
              &                    lPowerCalc, lFluxProfCalc, lPerpParCalc,    &
              &                    lHelCalc, l_frame)
      !
      ! This subroutine performs the MPI transposition from (n_theta,
      ! n_m_loc,n_r_loc) to (n_theta_loc, n_m_max, n_r_loc)
      !

      class(hybrid_3D_arrays_t) :: this
      logical, intent(in) :: lVisc, lRmsCalc, lPressCalc, lTOCalc, lPowerCalc
      logical, intent(in) :: lFluxProfCalc, l_frame, lPerpParCalc, lHelCalc
       
      if ( l_conv .or. l_mag_kin ) then
         if ( l_heat ) then
            call transpose_theta_m_many(this%s_Mloc, this%s_Thloc, 1)
            if ( lVisc ) then
               call transpose_theta_m_many(this%dsdt_Mloc, this%dsdt_Thloc, 1)
               call transpose_theta_m_many(this%dsdp_Mloc, this%dsdp_Thloc, 1)
            end if
            if ( l_HT .or. lVisc ) then
               call transpose_theta_m_many(this%dsdr_Mloc, this%dsdr_Thloc, 1)
            endif
         end if

         !-- Pressure
         if ( lPressCalc ) call transpose_theta_m_many(this%p_Mloc, &
                                &                 this%p_Thloc, 1)

         !-- Composition
         if ( l_chemical_conv ) call transpose_theta_m_many(this%xi_Mloc, &
                                     &                 this%xi_Thloc, 1)

         call transpose_theta_m_many(this%vr_Mloc, this%vr_Thloc, 1)
         call transpose_theta_m_many(this%vt_Mloc, this%vt_Thloc, 1)
         call transpose_theta_m_many(this%vp_Mloc, this%vp_Thloc, 1)
         call transpose_theta_m_many(this%cvr_Mloc, this%cvr_Thloc, 1)

         if ( l_adv_curl ) then
            call transpose_theta_m_many(this%cvt_Mloc, this%cvt_Thloc, 1)
            call transpose_theta_m_many(this%cvp_Mloc, this%cvp_Thloc, 1)

            !@> TODO: I thinks one also needs lHelCalc, lPerpParCalc, ...
            if ( lVisc .or. lPowerCalc .or. lRmsCalc .or. lFluxProfCalc &
            &    .or. lTOCalc .or. lHelCalc .or. lPerpParCalc .or.      &
            &    ( l_frame .and. l_movie_oc .and. l_store_frame) ) then
               call transpose_theta_m_many(this%dvrdr_Mloc, this%dvrdr_Thloc, 1)
               call transpose_theta_m_many(this%dvtdr_Mloc, this%dvtdr_Thloc, 1)
               call transpose_theta_m_many(this%dvpdr_Mloc, this%dvpdr_Thloc, 1)
               call transpose_theta_m_many(this%dvrdp_Mloc, this%dvrdp_Thloc, 1)
               call transpose_theta_m_many(this%dvtdp_Mloc, this%dvtdp_Thloc, 1)
               call transpose_theta_m_many(this%dvpdp_Mloc, this%dvpdp_Thloc, 1)
               call transpose_theta_m_many(this%dvrdt_Mloc, this%dvrdt_Thloc, 1)

            end if
         else
            call transpose_theta_m_many(this%dvrdr_Mloc, this%dvrdr_Thloc, 1)
            call transpose_theta_m_many(this%dvtdr_Mloc, this%dvtdr_Thloc, 1)
            call transpose_theta_m_many(this%dvpdr_Mloc, this%dvpdr_Thloc, 1)
            call transpose_theta_m_many(this%dvrdp_Mloc, this%dvrdp_Thloc, 1)
            call transpose_theta_m_many(this%dvtdp_Mloc, this%dvtdp_Thloc, 1)
            call transpose_theta_m_many(this%dvpdp_Mloc, this%dvpdp_Thloc, 1)
            call transpose_theta_m_many(this%dvrdt_Mloc, this%dvrdt_Thloc, 1)
         end if
      end if

      if ( l_mag .or. l_mag_LF ) then
         call transpose_theta_m_many(this%br_Mloc, this%br_Thloc, 1)
         call transpose_theta_m_many(this%bt_Mloc, this%bt_Thloc, 1)
         call transpose_theta_m_many(this%bp_Mloc, this%bp_Thloc, 1)
         call transpose_theta_m_many(this%cbr_Mloc, this%cbr_Thloc, 1)
         call transpose_theta_m_many(this%cbt_Mloc, this%cbt_Thloc, 1)
         call transpose_theta_m_many(this%cbp_Mloc, this%cbp_Thloc, 1)
      end if

   end subroutine transp_Mloc_to_Thloc
!-----------------------------------------------------------------------------------
   subroutine leg_hyb_to_spec(this, nl_lm, lRmsCalc)

      class(hybrid_3D_arrays_t) :: this

      !-- Input variable
      logical, intent(in) :: lRmsCalc

      !-- Output variables
      type(nonlinear_3D_lm_t) :: nl_lm

      !-- Local variables
      logical :: l_bound
      integer :: nR

      call shtns_load_cfg(1)
      do nR=nRstart,nRstop
         l_bound = (nR == n_r_cmb ) .or. ( nR == n_r_icb )

         if ( (.not.l_bound .or. lRmsCalc) .and. ( l_conv_nl .or. l_mag_LF ) ) then
            call hyb_to_SH(this%Advr_Mloc(:,:,nR), nl_lm%AdvrLM(:,nR), l_R(nR))
            call hyb_to_SH(this%Advt_Mloc(:,:,nR), nl_lm%AdvtLM(:,nR), l_R(nR))
            call hyb_to_SH(this%Advp_Mloc(:,:,nR), nl_lm%AdvpLM(:,nR), l_R(nR))

            if ( lRmsCalc .and. l_mag_LF .and. nR>n_r_LCR ) then
               ! LF treated extra:
               call hyb_to_SH(this%LFr_Mloc(:,:,nR), nl_lm%LFrLM(:,nR), l_R(nR))
               call hyb_to_SH(this%LFt_Mloc(:,:,nR), nl_lm%LFtLM(:,nR), l_R(nR))
               call hyb_to_SH(this%LFp_Mloc(:,:,nR), nl_lm%LFpLM(:,nR), l_R(nR))
            end if
         end if

         if ( l_heat .and. (.not. l_bound) ) then
            call hyb_to_qst(this%VSr_Mloc(:,:,nR), this%VSt_Mloc(:,:,nR), &
                 &          this%VSp_Mloc(:,:,nR), nl_lm%VSrLM(:,nR),     &
                 &          nl_lm%VStLM(:,nR), nl_lm%VSpLM(:,nR), l_R(nR))

            if ( l_anel ) then
               call hyb_to_SH(this%ViscHeat_Mloc(:,:,nR), nl_lm%ViscHeatLM(:,nR), &
                    &         l_R(nR))
               if ( l_mag_nl .and. nR>n_r_LCR ) then
                  call hyb_to_SH(this%OhmLoss_Mloc(:,:,nR), nl_lm%OhmLossLM(:,nR), &
                       &         l_R(nR))

               end if
            end if
         end if

         if ( l_chemical_conv .and. (.not. l_bound) ) then
            call hyb_to_qst(this%VXir_Mloc(:,:,nR), this%VXit_Mloc(:,:,nR), &
                 &          this%VXip_Mloc(:,:,nR), nl_lm%VXirLM(:,nR),     &
                 &          nl_lm%VXitLM(:,nR), nl_lm%VXipLM(:,nR), l_R(nR))
         end if

         if ( l_mag_nl ) then
            if ( .not.l_bound .and. nR>n_r_LCR ) then
               call hyb_to_qst(this%VxBr_Mloc(:,:,nR), this%VxBt_Mloc(:,:,nR), &
                    &          this%VxBp_Mloc(:,:,nR), nl_lm%VxBrLM(:,nR),     &
                    &          nl_lm%VxBtLM(:,nR), nl_lm%VxBpLM(:,nR), l_R(nR))
            else
               call hyb_to_sphertor(this%VxBt_Mloc(:,:,nR), this%VxBp_Mloc(:,:,nR), &
                    &               nl_lm%VxBtLM(:,nR), nl_lm%VxBpLM(:,nR), l_R(nR))
            end if
         end if

         if ( lRmsCalc ) then
            call hyb_to_sphertor(this%dpdt_Mloc(:,:,nR), this%dpdp_Mloc(:,:,nR), &
                 &               nl_lm%PFt2LM(:,nR), nl_lm%PFp2LM(:,nR), l_R(nR))
            call hyb_to_sphertor(this%CFt2_Mloc(:,:,nR), this%CFp2_Mloc(:,:,nR), &
                 &               nl_lm%CFt2LM(:,nR), nl_lm%CFp2LM(:,nR), l_R(nR))
            call hyb_to_qst(this%dtVr_Mloc(:,:,nR), this%dtVt_Mloc(:,:,nR), &
                 &          this%dtVp_Mloc(:,:,nR), nl_lm%dtVrLM(:,nR),     &
                 &          nl_lm%dtVtLM(:,nR), nl_lm%dtVpLM(:,nR), l_R(nR))
            if ( l_conv_nl ) then
               call hyb_to_sphertor(this%Advt2_Mloc(:,:,nR), this%Advp2_Mloc(:,:,nR),&
                    &               nl_lm%Advt2LM(:,nR), nl_lm%Advp2LM(:,nR), l_R(nR))
            end if
            if ( l_adv_curl ) then !-- Kinetic pressure : 1/2 d u^2 / dr
               call hyb_to_SH(this%dpkindr_Mloc(:,:,nR), nl_lm%dpkindrLM(:,nR), l_R(nR))
            end if
            if ( l_mag_nl .and. nR>n_r_LCR ) then
               call hyb_to_sphertor(this%LFt2_Mloc(:,:,nR), this%LFp2_Mloc(:,:,nR), &
                    &               nl_lm%LFt2LM(:,nR), nl_lm%LFp2LM(:,nR), l_R(nR))
            end if
         end if

      end do
      call shtns_load_cfg(0)

   end subroutine leg_hyb_to_spec
!-----------------------------------------------------------------------------------
   subroutine transp_Thloc_to_Mloc(this, lRmsCalc)

      class(hybrid_3D_arrays_t) :: this
      logical, intent(in) :: lRmsCalc

      call transpose_m_theta_many(this%Advr_Thloc, this%Advr_Mloc, 1)
      call transpose_m_theta_many(this%Advt_Thloc, this%Advt_Mloc, 1)
      call transpose_m_theta_many(this%Advp_Thloc, this%Advp_Mloc, 1)

      if ( l_heat ) then
         call transpose_m_theta_many(this%VSr_Thloc, this%VSr_Mloc, 1)
         call transpose_m_theta_many(this%VSt_Thloc, this%VSt_Mloc, 1)
         call transpose_m_theta_many(this%VSp_Thloc, this%VSp_Mloc, 1)
      end if

      if ( l_chemical_conv ) then
         call transpose_m_theta_many(this%VXir_Thloc, this%VXir_Mloc, 1)
         call transpose_m_theta_many(this%VXit_Thloc, this%VXit_Mloc, 1)
         call transpose_m_theta_many(this%VXip_Thloc, this%VXip_Mloc, 1)
      end if

      if ( l_anel ) then
         call transpose_m_theta_many(this%ViscHeat_Thloc, this%ViscHeat_Mloc, 1)
         if ( l_mag ) then
            call transpose_m_theta_many(this%OhmLoss_Thloc, this%OhmLoss_Mloc, 1)
         end if
      end if

      if ( l_mag_nl ) then
         call transpose_m_theta_many(this%VxBr_Thloc, this%VxBr_Mloc, 1)
         call transpose_m_theta_many(this%VxBt_Thloc, this%VxBt_Mloc, 1)
         call transpose_m_theta_many(this%VxBp_Thloc, this%VxBp_Mloc, 1)
      end if

      if ( lRmsCalc ) then
         call transpose_m_theta_many(this%dtVr_Thloc, this%dtVr_Mloc, 1)
         call transpose_m_theta_many(this%dtVt_Thloc, this%dtVt_Mloc, 1)
         call transpose_m_theta_many(this%dtVp_Thloc, this%dtVp_Mloc, 1)
         call transpose_m_theta_many(this%PFt2_Thloc, this%PFt2_Mloc, 1)
         call transpose_m_theta_many(this%PFp2_Thloc, this%PFp2_Mloc, 1)
         call transpose_m_theta_many(this%CFt2_Thloc, this%CFt2_Mloc, 1)
         call transpose_m_theta_many(this%CFp2_Thloc, this%CFp2_Mloc, 1)
         call transpose_m_theta_many(this%Advt2_Thloc, this%Advt2_Mloc, 1)
         call transpose_m_theta_many(this%Advp2_Thloc, this%Advp2_Mloc, 1)
         if ( l_adv_curl ) then
            call transpose_m_theta_many(this%dpkindr_Thloc, this%dpkindr_Mloc, 1)
         end if
         if ( l_mag_LF ) then
            call transpose_m_theta_many(this%LFr_Thloc, this%LFr_Mloc, 1)
            call transpose_m_theta_many(this%LFt_Thloc, this%LFt_Mloc, 1)
            call transpose_m_theta_many(this%LFp_Thloc, this%LFp_Mloc, 1)
            call transpose_m_theta_many(this%LFt2_Thloc, this%LFt2_Mloc, 1)
            call transpose_m_theta_many(this%LFp2_Thloc, this%LFp2_Mloc, 1)
         end if
      end if


   end subroutine transp_Thloc_to_Mloc
!-----------------------------------------------------------------------------------
end module hybrid_space_mod
