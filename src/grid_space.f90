!----------------------------------------------------------------------------
module grid_space
   !
   ! This module is used to compute the nonlinear products in physical space
   ! :math:`(\theta,\phi)`.
   !

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_phi_max, n_theta_max, n_theta_loc, nRstart, &
       &                 nRstop, nThetaStart, nThetaStop, n_r_loc
   use radial_functions, only: or2, orho1, beta, otemp1, visc, r, or3, &
       &                       lambda, or4, or1, alpha0, temp0, opressure0
   use physical_parameters, only: LFfac, n_r_LCR, CorFac, prec_angle,    &
       &                          ThExpNb, ViscHeatFac, oek, po, DissNb, &
       &                          dilution_fac, ra, opr, polind, strat, radratio
   use horizontal_data, only: sinTheta, cosTheta, phi, &
       &                      O_sin_theta_E2, cosn_theta_E2, O_sin_theta
   use parallel_mod, only: get_openmp_blocks
   use constants, only: two, third
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF, &
       &            l_RMS, l_chemical_conv, l_precession,             &
       &            l_centrifuge, l_adv_curl
   use time_schemes, only: type_tscheme
   use communications, only: slice_f, gather_f

   implicit none

   private

   type, public :: grid_space_t
      !----- Nonlinear terms in phi/theta space:
      real(cp), allocatable :: Advr(:,:,:), Advt(:,:,:), Advp(:,:,:)
      real(cp), allocatable :: LFr(:,:,:), LFt(:,:,:), LFp(:,:,:)
      real(cp), allocatable :: PCr(:,:,:), PCt(:,:,:), PCp(:,:,:)
      real(cp), allocatable :: CAr(:,:,:), CAt(:,:,:)
      real(cp), allocatable :: VxBr(:,:,:), VxBt(:,:,:), VxBp(:,:,:)
      real(cp), allocatable :: VSr(:,:,:), VSt(:,:,:), VSp(:,:,:)
      real(cp), allocatable :: VXir(:,:,:), VXit(:,:,:), VXip(:,:,:)
      real(cp), allocatable :: ViscHeat(:,:,:), OhmLoss(:,:,:)

      !---- RMS calculations
      real(cp), allocatable :: Advt2(:,:,:), Advp2(:,:,:), LFt2(:,:,:), LFp2(:,:,:)
      real(cp), allocatable :: CFt2(:,:,:), CFp2(:,:,:), dpdtc(:,:,:), dpdpc(:,:,:)
      real(cp), allocatable :: dtVr(:,:,:), dtVp(:,:,:), dtVt(:,:,:)
      real(cp), allocatable :: dpkindrc(:,:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), allocatable :: vrc(:,:,:), vtc(:,:,:), vpc(:,:,:)
      real(cp), allocatable :: dvrdrc(:,:,:), dvtdrc(:,:,:), dvpdrc(:,:,:)
      real(cp), allocatable :: cvrc(:,:,:), sc(:,:,:), drSc(:,:,:)
      real(cp), allocatable :: dvrdtc(:,:,:), dvrdpc(:,:,:)
      real(cp), allocatable :: dvtdpc(:,:,:), dvpdpc(:,:,:)
      real(cp), allocatable :: brc(:,:,:), btc(:,:,:), bpc(:,:,:)
      real(cp), allocatable :: cbrc(:,:,:), cbtc(:,:,:), cbpc(:,:,:)
      real(cp), allocatable :: pc(:,:,:), xic(:,:,:), cvtc(:,:,:), cvpc(:,:,:)
      real(cp), allocatable :: dsdtc(:,:,:), dsdpc(:,:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: get_nl

   end type grid_space_t

   real(cp), allocatable :: vr_old(:,:,:), vt_old(:,:,:), vp_old(:,:,:)

contains

!----------------------------------------------------------------------------
   subroutine initialize(this)
      !
      ! Memory allocation of arrays in grid space
      !

      class(grid_space_t) :: this

      allocate( this%Advr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),     &
      &         this%Advt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),     &
      &         this%Advp(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%LFr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),      &
      &         this%LFt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),      &
      &         this%LFp(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%VxBr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),     &
      &         this%VxBt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),     &
      &         this%VxBp(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%VSr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),      &
      &         this%VSt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),      &
      &         this%VSp(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%ViscHeat(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &          this%OhmLoss(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      bytes_allocated=bytes_allocated + 14*n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL

      if ( l_precession ) then
         allocate( this%PCr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%PCt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%PCp(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      if ( l_centrifuge ) then
         allocate( this%CAr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%CAt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 2*n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXir(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%VXit(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%VXip(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*n_phi_max*n_theta_loc*n_r_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &         this%vtc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &         this%vpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%dvrdrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &         this%dvtdrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%dvpdrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &         this%cvrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%dvrdtc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &         this%dvrdpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%dvtdpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
      &         this%dvpdpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%brc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),    &
      &         this%btc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),    &
      &         this%bpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      this%btc=1.0e50_cp
      this%bpc=1.0e50_cp
      allocate( this%cbrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),   &
      &         this%cbtc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),   &
      &         this%cbpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%sc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),     &
      &         this%drSc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%pc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      allocate( this%dsdtc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),  &
      &         this%dsdpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
      bytes_allocated=bytes_allocated + 22*n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1,nRstart:nRstop) ) ! TODO allocate this as (1,1,1) and fix the other parts of the code which refer to this variable
      end if

      if ( l_adv_curl ) then
         allocate( this%cvtc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%cvpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated+2*n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         allocate ( this%Advt2(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%Advp2(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         allocate ( this%dtVr(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),  &
         &         this%dtVt(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),   &
         &         this%dtVp(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         allocate ( this%LFt2(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),  &
         &         this%LFp2(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         allocate ( this%CFt2(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop),  &
         &         this%CFp2(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         allocate ( this%dpdtc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop), &
         &         this%dpdpc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 11*n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL

         allocate( vr_old(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         allocate( vp_old(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         allocate( vt_old(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*n_phi_max*n_theta_loc*n_r_loc*&
         &               SIZEOF_DEF_REAL

         this%dtVr(:,:,:)=0.0_cp
         this%dtVt(:,:,:)=0.0_cp
         this%dtVp(:,:,:)=0.0_cp
         vt_old(:,:,:) =0.0_cp
         vr_old(:,:,:) =0.0_cp
         vp_old(:,:,:) =0.0_cp

         if ( l_adv_curl ) then
            allocate ( this%dpkindrc(nThetaStart:nThetaStop,n_phi_max,nRstart:nRstop) )
            bytes_allocated=bytes_allocated + n_phi_max*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
         end if
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation of arrays in grid space
      !

      class(grid_space_t) :: this

      deallocate( this%Advr, this%Advt, this%Advp, this%LFr, this%LFt, this%LFp )
      deallocate( this%VxBr, this%VxBt, this%VxBp, this%VSr, this%VSt, this%VSp )
      if ( l_chemical_conv ) deallocate( this%VXir, this%VXit, this%VXip )
      if ( l_precession ) deallocate( this%PCr, this%PCt, this%PCp )
      if ( l_centrifuge ) deallocate( this%CAr, this%CAt )
      if ( l_adv_curl ) deallocate( this%cvtc, this%cvpc )
      deallocate( this%ViscHeat, this%OhmLoss )

      !----- Fields calculated from these help arrays by legtf:
      deallocate( this%vrc,this%vtc,this%vpc )
      deallocate( this%dvrdrc,this%dvtdrc,this%dvpdrc,this%cvrc )
      deallocate( this%dvrdtc,this%dvrdpc,this%dvtdpc,this%dvpdpc )
      deallocate( this%brc,this%btc,this%bpc,this%cbrc,this%cbtc,this%cbpc )
      deallocate( this%sc,this%drSc, this%pc, this%xic )
      deallocate( this%dsdtc, this%dsdpc )

      !-- RMS Calculations
      if ( l_RMS ) then
         deallocate ( this%Advt2, this%Advp2, this%LFt2, this%LFp2 )
         deallocate ( this%CFt2, this%CFp2, this%dpdtc, this%dpdpc )
         deallocate ( this%dtVr, this%dtVt, this%dtVp )
         if (allocated(vr_old)) deallocate ( vr_old )
         if (allocated(vt_old)) deallocate ( vt_old )
         if (allocated(vp_old)) deallocate ( vp_old )
         if ( l_adv_curl ) deallocate ( this%dpkindrc )
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_nl(this, time, tscheme, nR, nBc, lRmsCalc)
      !
      !  calculates non-linear products in grid-space for radial
      !  level ``nR`` and returns them in arrays
      !
      !  if ``nBc>0`` velocities are zero only the :math:`(\vec{u}\times\vec{B})`
      !  contributions need to be calculated
      !

      class(grid_space_t) :: this

      !-- Input of variables:
      real(cp),            intent(in) :: time    ! instant in time
      class(type_tscheme), intent(in) :: tscheme ! time scheme
      integer,             intent(in) :: nR      ! radial level
      logical,             intent(in) :: lRmsCalc
      integer,             intent(in) :: nBc

      !-- Local variables:
      integer :: nPhi, nPhStart, nPhStop
      real(cp) :: posnalp, O_dt

      if ( l_precession ) posnalp=-two*oek*po*sin(prec_angle)

      !$omp parallel default(shared) private(nPhStart,nPhStop,nPhi)
      nPhStart=1; nPhStop=n_phi_max
      call get_openmp_blocks(nPhStart,nPhStop)

      do nPhi=nPhStart,nPhStop

         if ( l_mag_LF .and. (nBc == 0 .or. lRmsCalc) .and. nR>n_r_LCR ) then
            !------ Get the Lorentz force:
            !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
            this%LFr(:,nPhi,nR)=  LFfac*O_sin_theta_E2(nThetaStart:nThetaStop) * (   &
            &        this%cbtc(:,nPhi,nR)*this%bpc(:,nPhi,nR) -       &
            &        this%cbpc(:,nPhi,nR)*this%btc(:,nPhi,nR) )

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            this%LFt(:,nPhi,nR)=  LFfac*or4(nR)*O_sin_theta_E2(nThetaStart:nThetaStop) * ( &
            &        this%cbpc(:,nPhi,nR)*this%brc(:,nPhi,nR) -             &
            &        this%cbrc(:,nPhi,nR)*this%bpc(:,nPhi,nR) )

            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            this%LFp(:,nPhi,nR)=  LFfac*or4(nR)*O_sin_theta_E2(nThetaStart:nThetaStop) * ( &
            &        this%cbrc(:,nPhi,nR)*this%btc(:,nPhi,nR) -             &
            &        this%cbtc(:,nPhi,nR)*this%brc(:,nPhi,nR) )
         end if      ! Lorentz force required ?

         if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

            if ( l_adv_curl ) then ! Advection is \curl{u} \times u
               this%Advr(:,nPhi,nR)=  - O_sin_theta_E2(nThetaStart:nThetaStop) * (    &
               &        this%cvtc(:,nPhi,nR)*this%vpc(:,nPhi,nR) -     &
               &        this%cvpc(:,nPhi,nR)*this%vtc(:,nPhi,nR) )

               this%Advt(:,nPhi,nR)= -or4(nR)*O_sin_theta_E2(nThetaStart:nThetaStop) * ( &
               &        this%cvpc(:,nPhi,nR)*this%vrc(:,nPhi,nR) -        &
               &        this%cvrc(:,nPhi,nR)*this%vpc(:,nPhi,nR) )

               this%Advp(:,nPhi,nR)= -or4(nR)*O_sin_theta_E2(nThetaStart:nThetaStop) * ( &
               &        this%cvrc(:,nPhi,nR)*this%vtc(:,nPhi,nR) -        &
               &        this%cvtc(:,nPhi,nR)*this%vrc(:,nPhi,nR) )
            else ! Advection is u\grad u
               !------ Get Advection:
               this%Advr(:,nPhi,nR)=          -or2(nR)*orho1(nR) * (  &
               &                                this%vrc(:,nPhi,nR) * &
               &                     (       this%dvrdrc(:,nPhi,nR) - &
               &    ( two*or1(nR)+beta(nR) )*this%vrc(:,nPhi,nR) ) +  &
               &                      O_sin_theta_E2(nThetaStart:nThetaStop) * (     &
               &                                this%vtc(:,nPhi,nR) * &
               &                     (       this%dvrdtc(:,nPhi,nR) - &
               &                  r(nR)*      this%vtc(:,nPhi,nR) ) + &
               &                                this%vpc(:,nPhi,nR) * &
               &                     (       this%dvrdpc(:,nPhi,nR) - &
               &                    r(nR)*      this%vpc(:,nPhi,nR) ) ) )

               this%Advt(:,nPhi,nR)=or4(nR)*O_sin_theta_E2(nThetaStart:nThetaStop)*orho1(nR) * (  &
               &                                         -this%vrc(:,nPhi,nR) * &
               &                                   (   this%dvtdrc(:,nPhi,nR) - &
               &                             beta(nR)*this%vtc(:,nPhi,nR) )   + &
               &                                          this%vtc(:,nPhi,nR) * &
               &                    ( cosn_theta_E2(nThetaStart:nThetaStop)*this%vtc(:,nPhi,nR) + &
               &                                       this%dvpdpc(:,nPhi,nR) + &
               &                                   this%dvrdrc(:,nPhi,nR) )   + &
               &                                          this%vpc(:,nPhi,nR) * &
               &                    ( cosn_theta_E2(nThetaStart:nThetaStop)*this%vpc(:,nPhi,nR) - &
               &                                       this%dvtdpc(:,nPhi,nR) )  )

               this%Advp(:,nPhi,nR)= or4(nR)*O_sin_theta_E2(nThetaStart:nThetaStop)*orho1(nR) * (  &
               &                                          -this%vrc(:,nPhi,nR) * &
               &                                      ( this%dvpdrc(:,nPhi,nR) - &
               &                              beta(nR)*this%vpc(:,nPhi,nR) )   - &
               &                                           this%vtc(:,nPhi,nR) * &
               &                                      ( this%dvtdpc(:,nPhi,nR) + &
               &                                      this%cvrc(:,nPhi,nR) )   - &
               &                     this%vpc(:,nPhi,nR) * this%dvpdpc(:,nPhi,nR) )
            end if

         end if  ! Navier-Stokes nonlinear advection term ?

         if ( l_heat_nl .and. nBc == 0 ) then
            !------ Get V S, the divergence of it is entropy advection:
            this%VSr(:,nPhi,nR)=this%vrc(:,nPhi,nR)*this%sc(:,nPhi,nR)
            this%VSt(:,nPhi,nR)=or2(nR)*this%vtc(:,nPhi,nR)*this%sc(:,nPhi,nR)
            this%VSp(:,nPhi,nR)=or2(nR)*this%vpc(:,nPhi,nR)*this%sc(:,nPhi,nR)
         end if     ! heat equation required ?

         if ( l_chemical_conv .and. nBc == 0 ) then
            this%VXir(:,nPhi,nR)=this%vrc(:,nPhi,nR)*this%xic(:,nPhi,nR)
            this%VXit(:,nPhi,nR)=or2(nR)*this%vtc(:,nPhi,nR)*this%xic(:,nPhi,nR)
            this%VXip(:,nPhi,nR)=or2(nR)*this%vpc(:,nPhi,nR)*this%xic(:,nPhi,nR)
         end if     ! chemical composition equation required ?

         if ( l_precession .and. nBc == 0 ) then
            this%PCr(:,nPhi,nR)=posnalp*O_sin_theta(nThetaStart:nThetaStop)*r(nR)*(                       &
            &                cos(oek*time+phi(nPhi))*this%vpc(:,nPhi,nR)*cosTheta(nThetaStart:nThetaStop)+&
            &                sin(oek*time+phi(nPhi))*this%vtc(:,nPhi,nR) )
            this%PCt(:,nPhi,nR)=   -posnalp*O_sin_theta(nThetaStart:nThetaStop)*or2(nR)*(            &
            &               cos(oek*time+phi(nPhi))*this%vpc(:,nPhi,nR)       + &
            &               sin(oek*time+phi(nPhi))*or1(nR)*this%vrc(:,nPhi,nR) )
            this%PCp(:,nPhi,nR)= posnalp*O_sin_theta(nThetaStart:nThetaStop)*cos(oek*time+phi(nPhi))*&
            &                 or2(nR)*(this%vtc(:,nPhi,nR)-or1(nR)*             &
            &                 this%vrc(:,nPhi,nR)*cosTheta(nThetaStart:nThetaStop))
         end if ! precession term required ?

         if ( l_centrifuge .and. nBc ==0 ) then
            !if ( l_anel ) then
            !   this%CAr(:,nPhi,nR) = dilution_fac*r(nR)*sinTheta(nThetaStart:nThetaStop)*sinTheta(nThetaStart:nThetaStop)* &
            !   &       ( -ra*opr*this%sc(:,nPhi,nR) )
            !   !-- neglect pressure contribution
            !   !& + polind*DissNb*oek*opressure0(nR)*this%pc(:,nPhi,nR) )
            !   this%CAt(:,nPhi,nR) = dilution_fac*r(nR)*sinTheta(nThetaStart:nThetaStop)*cosTheta(nThetaStart:nThetaStop)* &
            !   &       ( -ra*opr*this%sc(:,nPhi,nR) )
            !   !-- neglect pressure contribution
            !   !& + polind*DissNb*oek*opressure0(nR)*this%pc(:,nPhi,nR) )
            !else
            this%CAr(:,nPhi,nR) = -dilution_fac*r(nR)*sinTheta(nThetaStart:nThetaStop)*sinTheta(nThetaStart:nThetaStop)*ra*opr* &
            &                       this%sc(:,nPhi,nR)
            this%CAt(:,nPhi,nR) = -dilution_fac*r(nR)*sinTheta(nThetaStart:nThetaStop)*cosTheta(nThetaStart:nThetaStop)*ra*opr* &
            &                       this%sc(:,nPhi,nR)
            !end if
         end if ! centrifuge

         if ( l_mag_nl ) then

            if ( nBc == 0 .and. nR>n_r_LCR ) then
               !------ Get (V x B) , the curl of this is the dynamo term:
               this%VxBr(:,nPhi,nR)=  orho1(nR)*O_sin_theta_E2(nThetaStart:nThetaStop) * (  &
               &              this%vtc(:,nPhi,nR)*this%bpc(:,nPhi,nR) -      &
               &              this%vpc(:,nPhi,nR)*this%btc(:,nPhi,nR) )

               this%VxBt(:,nPhi,nR)=  orho1(nR)*or4(nR) * (   &
               &       this%vpc(:,nPhi,nR)*this%brc(:,nPhi,nR) - &
               &       this%vrc(:,nPhi,nR)*this%bpc(:,nPhi,nR) )

               this%VxBp(:,nPhi,nR)=   orho1(nR)*or4(nR) * (   &
               &        this%vrc(:,nPhi,nR)*this%btc(:,nPhi,nR) - &
               &        this%vtc(:,nPhi,nR)*this%brc(:,nPhi,nR) )
            else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary
               this%VxBt(:,nPhi,nR)= or4(nR)*orho1(nR)*this%vpc(:,nPhi,nR)*this%brc(:,nPhi,nR)
               this%VxBp(:,nPhi,nR)=-or4(nR)*orho1(nR)*this%vtc(:,nPhi,nR)*this%brc(:,nPhi,nR)
            else if ( nBc == 2 ) then  ! rigid boundary :
               !-- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):
               this%VxBt(:,nPhi,nR)=or4(nR)*orho1(nR)*this%vpc(:,nPhi,nR)*this%brc(:,nPhi,nR)
               this%VxBp(:,nPhi,nR)= 0.0_cp
            end if  ! boundary ?

         end if ! l_mag_nl ?

         if ( l_anel .and. nBc == 0 ) then
            !------ Get viscous heating
            this%ViscHeat(:,nPhi,nR)=      or4(nR)*                  &
            &                     orho1(nR)*otemp1(nR)*visc(nR)*( &
            &     two*(                     this%dvrdrc(:,nPhi,nR) - & ! (1)
            &     (two*or1(nR)+beta(nR))*this%vrc(:,nPhi,nR) )**2  + &
            &     two*( cosn_theta_E2(nThetaStart:nThetaStop)*   this%vtc(:,nPhi,nR) + &
            &                               this%dvpdpc(:,nPhi,nR) + &
            &                               this%dvrdrc(:,nPhi,nR) - & ! (2)
            &     or1(nR)*               this%vrc(:,nPhi,nR) )**2  + &
            &     two*(                     this%dvpdpc(:,nPhi,nR) + &
            &           cosn_theta_E2(nThetaStart:nThetaStop)*   this%vtc(:,nPhi,nR) + & ! (3)
            &     or1(nR)*               this%vrc(:,nPhi,nR) )**2  + &
            &          ( two*               this%dvtdpc(:,nPhi,nR) + &
            &                                 this%cvrc(:,nPhi,nR) - & ! (6)
            &    two*cosn_theta_E2(nThetaStart:nThetaStop)*this%vpc(:,nPhi,nR) )**2  + &
            &                        O_sin_theta_E2(nThetaStart:nThetaStop) * (     &
            &         ( r(nR)*              this%dvtdrc(:,nPhi,nR) - &
            &           (two+beta(nR)*r(nR))*  this%vtc(:,nPhi,nR) + & ! (4)
            &     or1(nR)*            this%dvrdtc(:,nPhi,nR) )**2  + &
            &         ( r(nR)*              this%dvpdrc(:,nPhi,nR) - &
            &           (two+beta(nR)*r(nR))*  this%vpc(:,nPhi,nR) + & ! (5)
            &     or1(nR)*            this%dvrdpc(:,nPhi,nR) )**2 )- &
            &    two*third*(  beta(nR)*        this%vrc(:,nPhi,nR) )**2 )

            if ( l_mag_nl .and. nR>n_r_LCR ) then
               !------ Get ohmic losses
               this%OhmLoss(:,nPhi,nR)= or2(nR)*otemp1(nR)*lambda(nR)*  &
               &    ( or2(nR)*                this%cbrc(:,nPhi,nR)**2 + &
               &      O_sin_theta_E2(nThetaStart:nThetaStop)*   this%cbtc(:,nPhi,nR)**2 + &
               &      O_sin_theta_E2(nThetaStart:nThetaStop)*   this%cbpc(:,nPhi,nR)**2  )
            end if ! if l_mag_nl ?

         end if  ! Viscous heating and Ohmic losses ?

         if ( lRmsCalc ) then
            this%dpdtc(:,nPhi,nR)=this%dpdtc(:,nPhi,nR)*or1(nR)
            this%dpdpc(:,nPhi,nR)=this%dpdpc(:,nPhi,nR)*or1(nR)
            this%CFt2(:,nPhi,nR)=-two*CorFac*cosTheta(nThetaStart:nThetaStop)*this%vpc(:,nPhi,nR)*or1(nR)
            this%CFp2(:,nPhi,nR)= two*CorFac*sinTheta(nThetaStart:nThetaStop)* (or1(nR)*cosTheta(nThetaStart:nThetaStop)*&
            &                             O_sin_theta(nThetaStart:nThetaStop)*this%vtc(:,nPhi,nR) + &
            &                        or2(nR)*sinTheta(nThetaStart:nThetaStop)*this%vrc(:,nPhi,nR) )
            if ( l_conv_nl ) then
               this%Advt2(:,nPhi,nR)=r(nR)*sinTheta(nThetaStart:nThetaStop)*sinTheta(nThetaStart:nThetaStop)*this%Advt(:,nPhi,nR)
               this%Advp2(:,nPhi,nR)=r(nR)*sinTheta(nThetaStart:nThetaStop)*sinTheta(nThetaStart:nThetaStop)*this%Advp(:,nPhi,nR)
            end if
            if ( l_mag_LF .and. nR > n_r_LCR ) then
               this%LFt2(:,nPhi,nR)=r(nR)*sinTheta(nThetaStart:nThetaStop)*sinTheta(nThetaStart:nThetaStop)*this%LFt(:,nPhi,nR)
               this%LFp2(:,nPhi,nR)=r(nR)*sinTheta(nThetaStart:nThetaStop)*sinTheta(nThetaStart:nThetaStop)*this%LFp(:,nPhi,nR)
            end if

            if ( l_adv_curl ) then
               this%dpdtc(:,nPhi,nR)=this%dpdtc(:,nPhi,nR)-or3(nR)*( or2(nR)*   &
               &            this%vrc(:,nPhi,nR)*this%dvrdtc(:,nPhi,nR) -        &
               &            this%vtc(:,nPhi,nR)*(this%dvrdrc(:,nPhi,nR)+        &
               &            this%dvpdpc(:,nPhi,nR)+cosn_theta_E2(nThetaStart:nThetaStop) *     &
               &            this%vtc(:,nPhi,nR))+ this%vpc(:,nPhi,nR)*(         &
               &            this%cvrc(:,nPhi,nR)+this%dvtdpc(:,nPhi,nR)-        &
               &            cosn_theta_E2(nThetaStart:nThetaStop)*this%vpc(:,nPhi,nR)) )
               this%dpdpc(:,nPhi,nR)=this%dpdpc(:,nPhi,nR)- or3(nR)*( or2(nR)*  &
               &            this%vrc(:,nPhi,nR)*this%dvrdpc(:,nPhi,nR) +        &
               &            this%vtc(:,nPhi,nR)*this%dvtdpc(:,nPhi,nR) +        &
               &            this%vpc(:,nPhi,nR)*this%dvpdpc(:,nPhi,nR) )
               if ( l_conv_nl ) then
                  this%Advt2(:,nPhi,nR)=this%Advt2(:,nPhi,nR)-or3(nR)*( or2(nR)*&
                  &            this%vrc(:,nPhi,nR)*this%dvrdtc(:,nPhi,nR) -     &
                  &            this%vtc(:,nPhi,nR)*(this%dvrdrc(:,nPhi,nR)+     &
                  &            this%dvpdpc(:,nPhi,nR)+cosn_theta_E2(nThetaStart:nThetaStop) *  &
                  &            this%vtc(:,nPhi,nR))+ this%vpc(:,nPhi,nR)*(      &
                  &            this%cvrc(:,nPhi,nR)+this%dvtdpc(:,nPhi,nR)-     &
                  &            cosn_theta_E2(nThetaStart:nThetaStop)*this%vpc(:,nPhi,nR)) )
                  this%Advp2(:,nPhi,nR)=this%Advp2(:,nPhi,nR)-or3(nR)*( or2(nR)*&
                  &            this%vrc(:,nPhi,nR)*this%dvrdpc(:,nPhi,nR) +     &
                  &            this%vtc(:,nPhi,nR)*this%dvtdpc(:,nPhi,nR) +     &
                  &            this%vpc(:,nPhi,nR)*this%dvpdpc(:,nPhi,nR) )
               end if

               !- dpkin/dr = 1/2 d (u^2) / dr = ur*dur/dr+ut*dut/dr+up*dup/dr
               this%dpkindrc(:,nPhi,nR)=or4(nR)*this%vrc(:,nPhi,nR)*(         &
               &                         this%dvrdrc(:,nPhi,nR)-           &
               &                         two*or1(nR)*this%vrc(:,nPhi,nR)) +&
               &                         or2(nR)*O_sin_theta_E2(nThetaStart:nThetaStop)*( &
               &                                 this%vtc(:,nPhi,nR)*(     &
               &                         this%dvtdrc(:,nPhi,nR)-           &
               &                         or1(nR)*this%vtc(:,nPhi,nR) ) +   &
               &                                 this%vpc(:,nPhi,nR)*(     &
               &                         this%dvpdrc(:,nPhi,nR)-           &
               &                         or1(nR)*this%vpc(:,nPhi,nR) ) )
            end if
         end if

         if ( l_RMS .and. tscheme%istage == 1 ) then
            O_dt = 1.0_cp/tscheme%dt(1)
            this%dtVr(:,nPhi,nR)=O_dt*or2(nR)*(this%vrc(:,nPhi,nR)-vr_old(:,nPhi,nR))
            this%dtVt(:,nPhi,nR)=O_dt*or1(nR)*(this%vtc(:,nPhi,nR)-vt_old(:,nPhi,nR))
            this%dtVp(:,nPhi,nR)=O_dt*or1(nR)*(this%vpc(:,nPhi,nR)-vp_old(:,nPhi,nR))

            vr_old(:,nPhi,nR)=this%vrc(:,nPhi,nR)
            vt_old(:,nPhi,nR)=this%vtc(:,nPhi,nR)
            vp_old(:,nPhi,nR)=this%vpc(:,nPhi,nR)
         end if

      end do
      !$omp end parallel

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space
