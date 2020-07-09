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
   use sht, only: scal_to_hyb, scal_to_grad_hyb, torpol_to_hyb,          &
       &            torpol_to_curl_hyb, pol_to_grad_hyb, torpol_to_dphhyb, &
       &            pol_to_curlr_hyb, hyb_to_SH, hyb_to_qst, hyb_to_sphertor
   use mpi_thetap_mod, only: transpose_theta_m_many, transpose_m_theta_many

   implicit none

   private

   complex(cp), target, allocatable :: vel_Mloc(:,:,:,:), gradvel_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: s_Mloc(:,:,:), grads_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: p_Mloc(:,:,:), gradp_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: xi_Mloc(:,:,:), mag_Mloc(:,:,:,:)

   complex(cp), target, allocatable :: NSadv_Mloc(:,:,:,:), heatadv_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: compadv_Mloc(:,:,:,:), emf_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: RMS_Mloc(:,:,:,:), LF_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: LF2_Mloc(:,:,:,:), PF2_Mloc(:,:,:,:)
   complex(cp), target, allocatable :: dtV_Mloc(:,:,:,:), anel_Mloc(:,:,:,:)

   type, public :: hybrid_3D_arrays_t
     
      complex(cp), contiguous, pointer :: vr_Mloc(:,:,:), vt_Mloc(:,:,:), vp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dvrdr_Mloc(:,:,:), dvtdr_Mloc(:,:,:), dvpdr_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: cvr_Mloc(:,:,:), dsdr_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dvrdt_Mloc(:,:,:), dvrdp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dvtdp_Mloc(:,:,:), dvpdp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: br_Mloc(:,:,:), bt_Mloc(:,:,:), bp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: cbr_Mloc(:,:,:), cbt_Mloc(:,:,:), cbp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: cvt_Mloc(:,:,:), cvp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dsdt_Mloc(:,:,:), dsdp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dpdt_Mloc(:,:,:), dpdp_Mloc(:,:,:)

      complex(cp), allocatable :: vel_pThloc(:,:,:,:), gradvel_pThloc(:,:,:,:)
      complex(cp), allocatable :: s_pThloc(:,:,:), grads_pThloc(:,:,:,:)
      complex(cp), allocatable :: p_pThloc(:,:,:), gradp_pThloc(:,:,:,:)
      complex(cp), allocatable :: xi_pThloc(:,:,:), mag_pThloc(:,:,:,:)

      complex(cp), contiguous, pointer :: Advr_Mloc(:,:,:), Advt_Mloc(:,:,:), Advp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: LFr_Mloc(:,:,:), LFt_Mloc(:,:,:), LFp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: VSr_Mloc(:,:,:), VSt_Mloc(:,:,:), VSp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: VXir_Mloc(:,:,:), VXit_Mloc(:,:,:), VXip_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: VxBr_Mloc(:,:,:), VxBt_Mloc(:,:,:), VxBp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: ViscHeat_Mloc(:,:,:), OhmLoss_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dpkindr_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: dtVr_Mloc(:,:,:), dtVt_Mloc(:,:,:), dtVp_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: PFt2_Mloc(:,:,:), PFp2_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: CFt2_Mloc(:,:,:), CFp2_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: Advt2_Mloc(:,:,:), Advp2_Mloc(:,:,:)
      complex(cp), contiguous, pointer :: LFt2_Mloc(:,:,:), LFp2_Mloc(:,:,:)

      complex(cp), allocatable :: NSadv_pThloc(:,:,:,:), heatadv_pThloc(:,:,:,:)
      complex(cp), allocatable :: compadv_pThloc(:,:,:,:), emf_pThloc(:,:,:,:)
      complex(cp), allocatable :: RMS_pThloc(:,:,:,:), LF_pThloc(:,:,:,:)
      complex(cp), allocatable :: LF2_pThloc(:,:,:,:), PF2_pThloc(:,:,:,:)
      complex(cp), allocatable :: dtV_pThloc(:,:,:,:), anel_pThloc(:,:,:,:)
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


      if ( l_adv_curl ) then
         allocate( vel_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,6) )
         bytes_allocated=bytes_allocated+6*n_r_loc*n_m_loc*n_theta_max*&
         &               SIZEOF_DEF_COMPLEX
      else
         allocate( vel_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,4) )
         bytes_allocated=bytes_allocated+4*n_r_loc*n_m_loc*n_theta_max*&
         &               SIZEOF_DEF_COMPLEX
      end if
      this%vr_Mloc(1:,1:,nRstart:) => &
      &          vel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
      this%vt_Mloc(1:,1:,nRstart:) => &
      &          vel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
      this%vp_Mloc(1:,1:,nRstart:) => &
      &          vel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
      this%cvr_Mloc(1:,1:,nRstart:) => &
      &          vel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,4)
      if ( l_adv_curl ) then
         this%cvt_Mloc(1:,1:,nRstart:) => &
         &          vel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,5)
         this%cvp_Mloc(1:,1:,nRstart:) => &
         &          vel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,6)
      end if

      allocate( gradvel_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,7) )
      bytes_allocated=bytes_allocated+7*n_r_loc*n_m_loc*n_theta_max*&
      &               SIZEOF_DEF_COMPLEX
      this%dvrdr_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
      this%dvtdr_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
      this%dvpdr_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
      this%dvrdp_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,4)
      this%dvtdp_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,5)
      this%dvpdp_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,6)
      this%dvrdt_Mloc(1:,1:,nRstart:) => &
      &          gradvel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,7)

      allocate( p_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      allocate( s_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
      bytes_allocated=bytes_allocated+2*n_r_loc*n_m_loc*n_theta_max*&
      &               SIZEOF_DEF_COMPLEX

      if ( l_adv_curl ) then
         allocate( this%vel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop, 6) )
         bytes_allocated=bytes_allocated+6*n_theta_loc*n_r_loc*(n_phi_max/2+1)* &
         &               SIZEOF_DEF_COMPLEX
      else
         allocate( this%vel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop, 4) )
         bytes_allocated=bytes_allocated+4*n_theta_loc*n_r_loc*(n_phi_max/2+1)* &
         &               SIZEOF_DEF_COMPLEX
      end if
      this%vel_pThloc(:,:,:,:)=zero

      allocate( this%p_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( this%s_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop) )
      this%p_pThloc(:,:,:)=zero
      this%s_pThloc(:,:,:)=zero
      bytes_allocated=bytes_allocated+2*n_theta_loc*n_r_loc*(n_phi_max/2+1)* &
      &               SIZEOF_DEF_COMPLEX

      allocate( this%gradvel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop, 7) )
      bytes_allocated=bytes_allocated+7*n_theta_loc*n_r_loc*(n_phi_max/2+1)* &
      &               SIZEOF_DEF_COMPLEX
      this%gradvel_pThloc(:,:,:,:)=zero

      if ( l_chemical_conv ) then
         allocate( xi_Mloc(n_theta_max,n_m_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+n_theta_max*n_m_loc*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         allocate( this%xi_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop) )
         this%xi_pThloc(:,:,:)=zero
         bytes_allocated = bytes_allocated+n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
      end if

      if ( l_mag ) then
         allocate( mag_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,6) )
         bytes_allocated = bytes_allocated+6*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%br_Mloc(1:,1:,nRstart:) => &
         &          mag_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%bt_Mloc(1:,1:,nRstart:) => &
         &          mag_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%bp_Mloc(1:,1:,nRstart:) => &
         &          mag_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
         this%cbr_Mloc(1:,1:,nRstart:) => &
         &          mag_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,4)
         this%cbt_Mloc(1:,1:,nRstart:) => &
         &          mag_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,5)
         this%cbp_Mloc(1:,1:,nRstart:) => &
         &          mag_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,6)

         allocate( this%mag_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,6) )
         bytes_allocated = bytes_allocated+6*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX

         this%mag_pThloc(:,:,:,:)=zero
      end if

      if ( l_viscBcCalc ) then
         allocate( grads_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%dsdr_Mloc(1:,1:,nRstart:) => &
         &          grads_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%dsdt_Mloc(1:,1:,nRstart:) => &
         &          grads_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%dsdp_Mloc(1:,1:,nRstart:) => &
         &          grads_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)

         allocate( this%grads_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
      else
         allocate( grads_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,1) )
         bytes_allocated = bytes_allocated+n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%dsdr_Mloc(1:,1:,nRstart:) => &
         &          grads_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         allocate( this%grads_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,1) )
         bytes_allocated = bytes_allocated+n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
      end if
      this%grads_pThloc(:,:,:,:)=zero

      if ( l_RMS ) then
         allocate( gradp_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,2) )
         bytes_allocated = bytes_allocated+2*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%dpdt_Mloc(1:,1:,nRstart:) => &
         &          gradp_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%dpdp_Mloc(1:,1:,nRstart:) => &
         &          gradp_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)

         allocate( this%gradp_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
         bytes_allocated = bytes_allocated+2*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated+2*n_theta_loc*n_m_max*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%gradp_pThloc(:,:,:,:)=zero
      end if

      allocate( this%NSadv_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
      bytes_allocated = bytes_allocated+3*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
      &                 SIZEOF_DEF_COMPLEX
      this%NSadv_pThloc(:,:,:,:)=zero

      if ( l_heat ) then
         allocate( this%heatadv_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%heatadv_pThloc(:,:,:,:)=zero
      end if 

      if ( l_chemical_conv ) then
         allocate( this%compadv_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%compadv_pThloc(:,:,:,:)=zero
      end if

      if ( l_anel ) then
         if ( l_mag ) then
            allocate( this%anel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
            bytes_allocated = bytes_allocated+2*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate( this%anel_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,1) )
            bytes_allocated = bytes_allocated+n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
            &                 SIZEOF_DEF_COMPLEX
         end if
      end if

      if ( l_mag_nl ) then
         allocate( this%emf_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%emf_pThloc(:,:,:,:)=zero
      end if

      allocate( NSadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3) )
      bytes_allocated = bytes_allocated+3*n_r_loc*n_m_loc*n_theta_max* &
      &                 SIZEOF_DEF_COMPLEX
      this%Advr_Mloc(1:,1:,nRstart:) => &
      &          NSadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
      this%Advt_Mloc(1:,1:,nRstart:) => &
      &          NSadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
      this%Advp_Mloc(1:,1:,nRstart:) => &
      &          NSadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)

      if ( l_heat ) then
         allocate( heatadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%VSr_Mloc(1:,1:,nRstart:) => &
         &          heatadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%VSt_Mloc(1:,1:,nRstart:) => &
         &          heatadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%VSp_Mloc(1:,1:,nRstart:) => &
         &          heatadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
      end if

      if ( l_chemical_conv ) then
         allocate( compadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%VXir_Mloc(1:,1:,nRstart:) => &
         &          compadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%VXit_Mloc(1:,1:,nRstart:) => &
         &          compadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%VXip_Mloc(1:,1:,nRstart:) => &
         &          compadv_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
      end if

      if ( l_anel ) then
         if ( l_mag ) then
            allocate( anel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2) )
            bytes_allocated = bytes_allocated+2*n_r_loc*n_m_loc*n_theta_max* &
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate( anel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1) )
            bytes_allocated = bytes_allocated+n_r_loc*n_m_loc*n_theta_max* &
            &                 SIZEOF_DEF_COMPLEX
         end if
         this%ViscHeat_Mloc(1:,1:,nRstart:) => &
         &          anel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         if ( l_mag ) then
            this%OhmLoss_Mloc(1:,1:,nRstart:) => &
            &          anel_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         end if
      end if

      if ( l_mag_nl ) then
         allocate( emf_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%VxBr_Mloc(1:,1:,nRstart:) => &
         &          emf_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%VxBt_Mloc(1:,1:,nRstart:) => &
         &          emf_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%VxBp_Mloc(1:,1:,nRstart:) => &
         &          emf_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
      end if

      if ( l_RMS ) then
         allocate( this%PF2_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
         bytes_allocated = bytes_allocated+2*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%PF2_pThloc(:,:,:,:)=zero
         if ( l_adv_curl ) then
            allocate( this%RMS_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,5) )
            bytes_allocated = bytes_allocated+5*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate( this%RMS_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,4) )
            bytes_allocated = bytes_allocated+4*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
            &                 SIZEOF_DEF_COMPLEX
         end if
         this%RMS_pThloc(:,:,:,:)=zero
         allocate( this%dtV_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
         &                 SIZEOF_DEF_COMPLEX
         this%dtV_pThloc(:,:,:,:)=zero

         if ( l_mag ) then
            allocate( this%LF_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,3) )
            allocate( this%LF2_pThloc(n_phi_max/2+1,nThetaStart:nThetaStop,nRstart:nRstop,2) )
            bytes_allocated = bytes_allocated+5*n_theta_loc*(n_phi_max/2+1)*n_r_loc*&
            &                 SIZEOF_DEF_COMPLEX
            this%LF_pThloc(:,:,:,:)=zero
            this%LF2_pThloc(:,:,:,:)=zero
         end if

         allocate( PF2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,2) )
         bytes_allocated = bytes_allocated+2*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%PFt2_Mloc(1:,1:,nRstart:) => &
         &          PF2_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%PFp2_Mloc(1:,1:,nRstart:) => &
         &          PF2_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         if ( l_adv_curl ) then
            allocate( RMS_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,5) )
            bytes_allocated = bytes_allocated+5*n_r_loc*n_m_loc*n_theta_max* &
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate( RMS_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,4) )
            bytes_allocated = bytes_allocated+4*n_r_loc*n_m_loc*n_theta_max* &
            &                 SIZEOF_DEF_COMPLEX
         end if
         this%CFt2_Mloc(1:,1:,nRstart:) => &
         &          RMS_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%CFp2_Mloc(1:,1:,nRstart:) => &
         &          RMS_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%Advt2_Mloc(1:,1:,nRstart:) => &
         &          RMS_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
         this%Advp2_Mloc(1:,1:,nRstart:) => &
         &          RMS_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,4)
         if ( l_adv_curl ) then
            this%dpkindr_Mloc(1:,1:,nRstart:) => &
            &          RMS_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,5)
         end if

         allocate( dtV_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,3) )
         bytes_allocated = bytes_allocated+3*n_r_loc*n_m_loc*n_theta_max* &
         &                 SIZEOF_DEF_COMPLEX
         this%dtVr_Mloc(1:,1:,nRstart:) => &
         &          dtV_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
         this%dtVt_Mloc(1:,1:,nRstart:) => &
         &          dtV_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
         this%dtVp_Mloc(1:,1:,nRstart:) => &
         &          dtV_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)

         if ( l_mag ) then
            allocate( LF_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,3) )
            this%LFr_Mloc(1:,1:,nRstart:) => &
            &          LF_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
            this%LFt_Mloc(1:,1:,nRstart:) => &
            &          LF_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
            this%LFp_Mloc(1:,1:,nRstart:) => &
            &          LF_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,3)
            allocate( LF2_Mloc(n_theta_max,n_m_loc,nRstart:nRstop,2) )
            this%LFt2_Mloc(1:,1:,nRstart:) => &
            &          LF2_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,1)
            this%LFp2_Mloc(1:,1:,nRstart:) => &
            &          LF2_Mloc(1:n_theta_max,1:n_m_loc,nRstart:nRstop,2)
            bytes_allocated = bytes_allocated+5*n_r_loc*n_m_loc*n_theta_max* &
            &                 SIZEOF_DEF_COMPLEX
         end if
      end if

   end subroutine initialize
!-----------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation for arrays in hybrid space
      !

      class(hybrid_3D_arrays_t) :: this

      if ( l_RMS ) then
         deallocate( gradp_Mloc )
         deallocate( this%RMS_pThloc, this%dtV_pThloc, this%PF2_pThloc )
         if ( l_mag ) deallocate( this%LF_pThloc, this%LF2_pThloc )
         deallocate( RMS_Mloc, dtV_Mloc, PF2_Mloc )
         if ( l_mag ) deallocate( LF_Mloc, LF2_Mloc )
      end if
      if ( l_mag ) deallocate( this%mag_pThloc, mag_Mloc )
      if ( l_chemical_conv ) then
         deallocate( xi_Mloc, this%xi_pThloc )
      end if

      deallocate( vel_Mloc, gradvel_Mloc )
      deallocate( s_Mloc, p_Mloc, grads_Mloc)
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
               call scal_to_hyb(s(:,nR), s_Mloc(:,:,nR), l_R(nR))
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
            if ( lPressCalc ) call scal_to_hyb(p(:,nR), p_Mloc(:,:,nR), l_R(nR))

            !-- Composition
            if ( l_chemical_conv ) call scal_to_hyb(xi(:,nR), xi_Mloc(:,:,nR),&
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
            call transpose_theta_m_many(s_Mloc, this%s_pThloc, n_r_loc)
            if ( lVisc ) then
               call transpose_theta_m_many(grads_Mloc, this%grads_pThloc, 3*n_r_loc)
            end if
            if ( l_HT .and. (.not. lVisc) ) then
               call transpose_theta_m_many(grads_Mloc, this%grads_pThloc, n_r_loc)
            endif
         end if

         !-- Pressure
         if ( lPressCalc ) call transpose_theta_m_many(p_Mloc, &
                                &                 this%p_pThloc, n_r_loc)

         !-- Composition
         if ( l_chemical_conv ) call transpose_theta_m_many(xi_Mloc, &
                                     &                 this%xi_pThloc, n_r_loc)
         if ( l_adv_curl ) then
            call transpose_theta_m_many(vel_Mloc, this%vel_pThloc, 6*n_r_loc)
            if ( lVisc .or. lPowerCalc .or. lRmsCalc .or. lFluxProfCalc &
            &    .or. lTOCalc .or. lHelCalc .or. lPerpParCalc .or.      &
            &    ( l_frame .and. l_movie_oc .and. l_store_frame) ) then
               call transpose_theta_m_many(gradvel_Mloc, this%gradvel_pThloc, &
                    &                      7*n_r_loc)
            end if
         else
            call transpose_theta_m_many(vel_Mloc, this%vel_pThloc, 4*n_r_loc)
            call transpose_theta_m_many(gradvel_Mloc, this%gradvel_pThloc, &
                 &                      7*n_r_loc)
         end if
      end if

      if ( l_mag .or. l_mag_LF ) then
         call transpose_theta_m_many(mag_Mloc, this%mag_pThloc, 6*n_r_loc)
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

      call transpose_m_theta_many(this%NSadv_pThloc, NSadv_Mloc, 3*n_r_loc)

      if ( l_heat ) then
         call transpose_m_theta_many(this%heatadv_pThloc, heatadv_Mloc, 3*n_r_loc)
      end if

      if ( l_chemical_conv ) then
         call transpose_m_theta_many(this%compadv_pThloc, compadv_Mloc, 3*n_r_loc)
      end if

      if ( l_anel ) then
         if ( l_mag ) then
            call transpose_m_theta_many(this%anel_pThloc, anel_Mloc, 2*n_r_loc)
         else
            call transpose_m_theta_many(this%anel_pThloc, anel_Mloc, n_r_loc)
         end if
      end if

      if ( l_mag_nl ) then
         call transpose_m_theta_many(this%emf_pThloc, emf_Mloc, 3*n_r_loc)
      end if

      if ( lRmsCalc ) then
         call transpose_m_theta_many(this%dtV_pThloc, dtV_Mloc, 3*n_r_loc)
         call transpose_m_theta_many(this%PF2_pThloc, PF2_Mloc, 2*n_r_loc)
         if ( l_adv_curl ) then
            call transpose_m_theta_many(this%RMS_pThloc, RMS_Mloc, 5*n_r_loc)
         else
            call transpose_m_theta_many(this%RMS_pThloc, RMS_Mloc, 4*n_r_loc)
         end if
         if ( l_mag_LF ) then
            call transpose_m_theta_many(this%LF_pThloc, LF_Mloc, 3*n_r_loc)
            call transpose_m_theta_many(this%LF2_pThloc, LF2_Mloc, 2*n_r_loc)
         end if
      end if

   end subroutine transp_Thloc_to_Mloc
!-----------------------------------------------------------------------------------
end module hybrid_space_mod
