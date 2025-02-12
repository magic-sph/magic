module general_arrays_mod

   implicit none

   private

   type, public, abstract :: general_arrays_t

   end type general_arrays_t

end module general_arrays_mod
!----------------------------------------------------------------------------
module grid_space_arrays_mod
   !
   ! This module is used to compute the nonlinear products in physical space
   ! :math:`(\theta,\phi)`.
   !

   use general_arrays_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_phi_max, nlat_padded
   use radial_functions, only: or2, orho1, beta, otemp1, visc, r, or3, &
       &                       lambda, or4, or1
   use physical_parameters, only: LFfac, n_r_LCR, prec_angle, ViscHeatFac,    &
        &                         oek, po, dilution_fac, ra, opr, OhmLossFac, &
        &                         epsPhase, phaseDiffFac, penaltyFac, tmelt
   use horizontal_data, only: sinTheta, cosTheta, phi, O_sin_theta_E2, &
       &                      cosn_theta_E2, O_sin_theta
   use parallel_mod, only: get_openmp_blocks
   use constants, only: two, third, one
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF, l_adv_curl, &
       &            l_chemical_conv, l_precession, l_centrifuge, l_phase_field

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space:
      real(cp), allocatable :: Advr(:,:), Advt(:,:), Advp(:,:)
      real(cp), allocatable :: LFr(:,:), LFt(:,:), LFp(:,:)
      real(cp), allocatable :: PCr(:,:), PCt(:,:), PCp(:,:)
      real(cp), allocatable :: CAr(:,:), CAt(:,:)
      real(cp), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)
      real(cp), allocatable :: VSr(:,:), VSt(:,:), VSp(:,:)
      real(cp), allocatable :: VXir(:,:), VXit(:,:), VXip(:,:)
      real(cp), allocatable :: heatTerms(:,:), phiTerms(:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), allocatable :: vrc(:,:), vtc(:,:), vpc(:,:)
      real(cp), allocatable :: dvrdrc(:,:), dvtdrc(:,:), dvpdrc(:,:)
      real(cp), allocatable :: cvrc(:,:), sc(:,:), drSc(:,:)
      real(cp), allocatable :: dvrdtc(:,:), dvrdpc(:,:)
      real(cp), allocatable :: dvtdpc(:,:), dvpdpc(:,:)
      real(cp), allocatable :: brc(:,:), btc(:,:), bpc(:,:)
      real(cp), allocatable :: cbrc(:,:), cbtc(:,:), cbpc(:,:)
      real(cp), allocatable :: pc(:,:), xic(:,:), cvtc(:,:), cvpc(:,:)
      real(cp), allocatable :: dsdtc(:,:), dsdpc(:,:), phic(:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: get_nl

   end type grid_space_arrays_t

contains

   subroutine initialize(this)
      !
      ! Memory allocation of arrays in grid space
      !

      class(grid_space_arrays_t) :: this

      allocate( this%Advr(nlat_padded,n_phi_max), this%Advt(nlat_padded,n_phi_max) )
      this%Advr(:,:)=0.0_cp
      this%Advt(:,:)=0.0_cp
      allocate( this%Advp(nlat_padded,n_phi_max), this%LFr(nlat_padded,n_phi_max) )
      this%Advp(:,:)=0.0_cp
      this%LFr(:,:) =0.0_cp
      allocate( this%LFt(nlat_padded,n_phi_max), this%LFp(nlat_padded,n_phi_max) )
      this%LFt(:,:) =0.0_cp
      this%LFp(:,:) =0.0_cp
      allocate( this%VxBr(nlat_padded,n_phi_max), this%VxBt(nlat_padded,n_phi_max) )
      this%VxBr(:,:)=0.0_cp
      this%VxBt(:,:)=0.0_cp
      allocate( this%VxBp(nlat_padded,n_phi_max), this%VSr(nlat_padded,n_phi_max) )
      this%VxBp(:,:)=0.0_cp
      this%VSr(:,:) =0.0_cp
      allocate( this%VSt(nlat_padded,n_phi_max), this%VSp(nlat_padded,n_phi_max) )
      this%VSt(:,:) =0.0_cp
      this%VSp(:,:) =0.0_cp
      allocate( this%heatTerms(nlat_padded,n_phi_max) )
      this%heatTerms(:,:)=0.0_cp
      bytes_allocated=bytes_allocated + 13*n_phi_max*nlat_padded*SIZEOF_DEF_REAL

      if ( l_precession ) then
         allocate( this%PCr(nlat_padded,n_phi_max), this%PCt(nlat_padded,n_phi_max) )
         this%PCr(:,:)=0.0_cp
         this%PCt(:,:)=0.0_cp
         allocate( this%PCp(nlat_padded,n_phi_max) )
         this%PCp(:,:)=0.0_cp
         bytes_allocated=bytes_allocated + 3*n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      end if

      if ( l_centrifuge ) then
         allocate( this%CAr(nlat_padded,n_phi_max), this%CAt(nlat_padded,n_phi_max) )
         this%CAr(:,:)=0.0_cp
         this%CAt(:,:)=0.0_cp
         bytes_allocated=bytes_allocated + 2*n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXir(nlat_padded,n_phi_max), this%VXit(nlat_padded,n_phi_max))
         this%VXir(:,:)=0.0_cp
         this%VXit(:,:)=0.0_cp
         allocate( this%VXip(nlat_padded,n_phi_max) )
         this%VXip(:,:)=0.0_cp
         bytes_allocated=bytes_allocated + 3*n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(nlat_padded,n_phi_max),this%vtc(nlat_padded,n_phi_max) )
      this%vrc(:,:)   =0.0_cp
      this%vtc(:,:)   =0.0_cp
      allocate( this%vpc(nlat_padded,n_phi_max), this%dvrdrc(nlat_padded,n_phi_max) )
      this%vpc(:,:)   =0.0_cp
      this%dvrdrc(:,:)=0.0_cp
      allocate(this%dvtdrc(nlat_padded,n_phi_max), this%dvpdrc(nlat_padded,n_phi_max))
      this%dvtdrc(:,:)=0.0_cp
      this%dvpdrc(:,:)=0.0_cp
      allocate( this%cvrc(nlat_padded,n_phi_max), this%dvrdtc(nlat_padded,n_phi_max) )
      this%cvrc(:,:)  =0.0_cp
      this%dvrdtc(:,:)=0.0_cp
      allocate( this%dvrdpc(nlat_padded,n_phi_max), this%dvtdpc(nlat_padded,n_phi_max) )
      this%dvrdpc(:,:)=0.0_cp
      this%dvtdpc(:,:)=0.0_cp
      allocate( this%dvpdpc(nlat_padded,n_phi_max), this%brc(nlat_padded,n_phi_max) )
      this%dvpdpc(:,:)=0.0_cp
      this%brc(:,:)   =0.0_cp
      allocate( this%btc(nlat_padded,n_phi_max),this%bpc(nlat_padded,n_phi_max) )
      this%btc(:,:)   =1.0e50_cp
      this%bpc(:,:)   =1.0e50_cp
      allocate( this%cbrc(nlat_padded,n_phi_max),this%cbtc(nlat_padded,n_phi_max) )
      this%cbrc(:,:)  =0.0e0_cp
      this%cbtc(:,:)  =0.0e0_cp
      allocate( this%cbpc(nlat_padded,n_phi_max), this%sc(nlat_padded,n_phi_max) )
      this%cbpc(:,:)  =0.0e0_cp
      this%sc(:,:)    =0.0e0_cp
      allocate( this%drSc(nlat_padded,n_phi_max), this%pc(nlat_padded,n_phi_max) )
      this%drSc(:,:)  =0.0e0_cp
      this%pc(:,:)    =0.0e0_cp
      allocate( this%dsdtc(nlat_padded,n_phi_max),this%dsdpc(nlat_padded,n_phi_max) )
      this%dsdtc(:,:) =0.0_cp
      this%dsdpc(:,:) =0.0_cp
      bytes_allocated=bytes_allocated + 22*n_phi_max*nlat_padded*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(nlat_padded,n_phi_max) )
         this%xic(:,:)=0.0_cp
         bytes_allocated=bytes_allocated + n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1) )
      end if

      if ( l_phase_field ) then
         allocate( this%phic(nlat_padded,n_phi_max),this%phiTerms(nlat_padded,n_phi_max) )
         this%phic(:,:)    =0.0_cp
         this%phiTerms(:,:)=0.0_cp
         bytes_allocated=bytes_allocated + 2*n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      else
         allocate( this%phic(1,1), this%phiTerms(1,1) )
      end if

      if ( l_adv_curl ) then
         allocate( this%cvtc(nlat_padded,n_phi_max), this%cvpc(nlat_padded,n_phi_max) )
         this%cvtc(:,:)=0.0_cp
         this%cvpc(:,:)=0.0_cp
         bytes_allocated=bytes_allocated+2*n_phi_max*nlat_padded*SIZEOF_DEF_REAL
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation of arrays in grid space
      !

      class(grid_space_arrays_t) :: this

      deallocate( this%Advr, this%Advt, this%Advp, this%LFr, this%LFt, this%LFp )
      deallocate( this%VxBr, this%VxBt, this%VxBp, this%VSr, this%VSt, this%VSp )
      if ( l_chemical_conv ) deallocate( this%VXir, this%VXit, this%VXip )
      if ( l_precession ) deallocate( this%PCr, this%PCt, this%PCp )
      if ( l_centrifuge ) deallocate( this%CAr, this%CAt )
      if ( l_adv_curl ) deallocate( this%cvtc, this%cvpc )
      if ( l_phase_field ) deallocate( this%phic, this%phiTerms )
      deallocate( this%heatTerms )

      !----- Fields calculated from these help arrays by legtf:
      deallocate( this%vrc,this%vtc,this%vpc )
      deallocate( this%dvrdrc,this%dvtdrc,this%dvpdrc,this%cvrc )
      deallocate( this%dvrdtc,this%dvrdpc,this%dvtdpc,this%dvpdpc )
      deallocate( this%brc,this%btc,this%bpc,this%cbrc,this%cbtc,this%cbpc )
      deallocate( this%sc,this%drSc, this%pc, this%xic )
      deallocate( this%dsdtc, this%dsdpc )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_nl(this, time, nR, nBc, lRmsCalc)
      !
      !  calculates non-linear products in grid-space for radial
      !  level ``nR`` and returns them in arrays
      !
      !  if ``nBc>0`` velocities are zero only the :math:`(\vec{u}\times\vec{B})`
      !  contributions need to be calculated
      !

      class(grid_space_arrays_t) :: this

      !-- Input of variables:
      real(cp),            intent(in) :: time    ! instant in time
      integer,             intent(in) :: nR      ! radial level
      logical,             intent(in) :: lRmsCalc
      integer,             intent(in) :: nBc

      !-- Local variables:
      integer :: nPhi, nPhStart, nPhStop
      real(cp) :: posnalp

      if ( l_precession ) posnalp=-two*oek*po*sin(prec_angle)

      !$omp parallel default(shared) private(nPhStart,nPhStop,nPhi)
      nPhStart=1; nPhStop=n_phi_max
      call get_openmp_blocks(nPhStart,nPhStop)

      do nPhi=nPhStart,nPhStop

         if ( l_mag_LF .and. (nBc == 0 .or. lRmsCalc) .and. nR>n_r_LCR ) then
            !------ Get the Lorentz force:
            !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
            this%LFr(:,nPhi)=  LFfac*O_sin_theta_E2(:) * (   &
            &        this%cbtc(:,nPhi)*this%bpc(:,nPhi) -    &
            &        this%cbpc(:,nPhi)*this%btc(:,nPhi) )

            !---- LFt= 1/(E*Pm) * sin(theta)/r * ( curl(B)_p*B_r - curl(B)_r*B_p )
            this%LFt(:,nPhi)=  LFfac*or4(nR) * (             &
            &        this%cbpc(:,nPhi)*this%brc(:,nPhi) -    &
            &        this%cbrc(:,nPhi)*this%bpc(:,nPhi) )

            !---- LFp= 1/(E*Pm) * sin(theta)/r * ( curl(B)_r*B_t - curl(B)_t*B_r )
            this%LFp(:,nPhi)=  LFfac*or4(nR) * (             &
            &        this%cbrc(:,nPhi)*this%btc(:,nPhi) -    &
            &        this%cbtc(:,nPhi)*this%brc(:,nPhi) )
         end if      ! Lorentz force required ?

         if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

            if ( l_adv_curl ) then ! Advection is \curl{u} \times u
               this%Advr(:,nPhi)=  - O_sin_theta_E2(:) * (    &
               &        this%cvtc(:,nPhi)*this%vpc(:,nPhi) -  &
               &        this%cvpc(:,nPhi)*this%vtc(:,nPhi) )

               this%Advt(:,nPhi)= -or4(nR)* (                    &
               &        this%cvpc(:,nPhi)*this%vrc(:,nPhi) -     &
               &        this%cvrc(:,nPhi)*this%vpc(:,nPhi) )

               this%Advp(:,nPhi)= -or4(nR)* (                    &
               &        this%cvrc(:,nPhi)*this%vtc(:,nPhi) -     &
               &        this%cvtc(:,nPhi)*this%vrc(:,nPhi) )
            else ! Advection is u\grad u
               !------ Get Advection:
               this%Advr(:,nPhi)=          -or2(nR)*orho1(nR) * (  &
               &                                this%vrc(:,nPhi) * &
               &                     (       this%dvrdrc(:,nPhi) - &
               &    ( two*or1(nR)+beta(nR) )*this%vrc(:,nPhi) ) +  &
               &                      O_sin_theta_E2(:) * (        &
               &                                this%vtc(:,nPhi) * &
               &                     (       this%dvrdtc(:,nPhi) - &
               &                  r(nR)*      this%vtc(:,nPhi) ) + &
               &                                this%vpc(:,nPhi) * &
               &                     (       this%dvrdpc(:,nPhi) - &
               &                    r(nR)*      this%vpc(:,nPhi) ) ) )

               this%Advt(:,nPhi)=or4(nR)*orho1(nR) * (                       &
               &                                         -this%vrc(:,nPhi) * &
               &                                   (   this%dvtdrc(:,nPhi) - &
               &                             beta(nR)*this%vtc(:,nPhi) )   + &
               &                                          this%vtc(:,nPhi) * &
               &                    ( cosn_theta_E2(:)*this%vtc(:,nPhi) +    &
               &                                       this%dvpdpc(:,nPhi) + &
               &                                   this%dvrdrc(:,nPhi) )   + &
               &                                          this%vpc(:,nPhi) * &
               &                    ( cosn_theta_E2(:)*this%vpc(:,nPhi) -    &
               &                                       this%dvtdpc(:,nPhi) )  )

               this%Advp(:,nPhi)= or4(nR)*orho1(nR) * (                       &
               &                                          -this%vrc(:,nPhi) * &
               &                                      ( this%dvpdrc(:,nPhi) - &
               &                              beta(nR)*this%vpc(:,nPhi) )   - &
               &                                           this%vtc(:,nPhi) * &
               &                                      ( this%dvtdpc(:,nPhi) + &
               &                                      this%cvrc(:,nPhi) )   - &
               &                     this%vpc(:,nPhi) * this%dvpdpc(:,nPhi) )
            end if

         end if  ! Navier-Stokes nonlinear advection term ?

         !-- Coriolis force in physical space
         !this%Advr(:,nPhi)=this%Advr(:,nPhi) + two*CorFac*r(nR)*this%vpc(:,nPhi)
         !this%Advt(:,nPhi)=this%Advt(:,nPhi) + two*CorFac*cosTheta(:)* &
         !&                 this%vpc(:,nPhi)*or2(nR)
         !this%Advp(:,nPhi)=this%Advp(:,nPhi) - two*CorFac*sinTheta(:)*or1(nR)*    &
         !&                 (or1(nR)*cosTheta(:)*O_sin_theta(:)*this%vtc(:,nPhi) + &
         !&                                 or2(nR)*sinTheta(:)*this%vrc(:,nPhi) )

         if ( l_heat_nl .and. nBc == 0 ) then
            !------ Get V S, the divergence of it is entropy advection:
            this%VSr(:,nPhi)=this%vrc(:,nPhi)*this%sc(:,nPhi)
            this%VSt(:,nPhi)=or2(nR)*this%vtc(:,nPhi)*this%sc(:,nPhi)
            this%VSp(:,nPhi)=or2(nR)*this%vpc(:,nPhi)*this%sc(:,nPhi)
         end if     ! heat equation required ?

         if ( l_chemical_conv .and. nBc == 0 ) then
            this%VXir(:,nPhi)=this%vrc(:,nPhi)*this%xic(:,nPhi)
            this%VXit(:,nPhi)=or2(nR)*this%vtc(:,nPhi)*this%xic(:,nPhi)
            this%VXip(:,nPhi)=or2(nR)*this%vpc(:,nPhi)*this%xic(:,nPhi)
         end if     ! chemical composition equation required ?

         if ( l_phase_field .and. nBc == 0 ) then
            this%Advr(:,nPhi)=this%Advr(:,nPhi)-this%phic(:,nPhi)*this%vrc(:,nPhi)/ &
            &                 epsPhase**2/penaltyFac**2
            this%Advt(:,nPhi)=this%Advt(:,nPhi)-or2(nR)*this%phic(:,nPhi)* &
            &                 this%vtc(:,nPhi)/epsPhase**2/penaltyFac**2
            this%Advp(:,nPhi)=this%Advp(:,nPhi)-or2(nR)*this%phic(:,nPhi)* &
            &                 this%vpc(:,nPhi)/epsPhase**2/penaltyFac**2
            this%phiTerms(:,nPhi)=-one/epsPhase**2* this%phic(:,nPhi)*       &
            &                      (one-this%phic(:,nPhi))*(                 &
            &                      phaseDiffFac*(one-two*this%phic(:,nPhi))+ &
            &                      this%sc(:,nPhi)-tmelt)
         end if ! Nonlinear terms for the phase field equation

         if ( l_precession .and. nBc == 0 ) then
            this%PCr(:,nPhi)=posnalp*O_sin_theta(:)*r(nR)*(                       &
            &                cos(oek*time+phi(nPhi))*this%vpc(:,nPhi)*cosTheta(:)+&
            &                sin(oek*time+phi(nPhi))*this%vtc(:,nPhi) )
            this%PCt(:,nPhi)=   -posnalp*sinTheta(:)*or2(nR)*(               &
            &               cos(oek*time+phi(nPhi))*this%vpc(:,nPhi)       + &
            &               sin(oek*time+phi(nPhi))*or1(nR)*this%vrc(:,nPhi) )
            this%PCp(:,nPhi)= posnalp*sinTheta(:)*cos(oek*time+phi(nPhi))*   &
            &                 or2(nR)*(this%vtc(:,nPhi)-or1(nR)*             &
            &                 this%vrc(:,nPhi)*cosTheta(:))
         end if ! precession term required ?

         if ( l_centrifuge .and. nBc ==0 ) then
            !if ( l_anel ) then
            !   this%CAr(:,nPhi) = dilution_fac*r(nR)*sinTheta(:)**4* &
            !   &       ( -ra*opr*this%sc(:,nPhi) )
            !   !-- neglect pressure contribution
            !   !& + polind*DissNb*oek*opressure0(nR)*this%pc(:,nPhi) )
            !   this%CAt(:,nPhi) = dilution_fac*r(nR)*sinTheta(:)**3*cosTheta(:)* &
            !   &       ( -ra*opr*this%sc(:,nPhi) )
            !   !-- neglect pressure contribution
            !   !& + polind*DissNb*oek*opressure0(nR)*this%pc(:,nPhi) )
            !else
            this%CAr(:,nPhi) = -dilution_fac*r(nR)*sinTheta(:)**4*ra*opr* &
            &                       this%sc(:,nPhi)
            this%CAt(:,nPhi) = -dilution_fac*r(nR)*sinTheta(:)**3*cosTheta(:)*ra*opr* &
            &                       this%sc(:,nPhi)
            !end if
         end if ! centrifuge

         if ( l_mag_nl ) then

            if ( nBc == 0 .and. nR>n_r_LCR ) then
               !------ Get (V x B) , the curl of this is the dynamo term:
               this%VxBr(:,nPhi)=  orho1(nR)*O_sin_theta_E2(:) * (  &
               &              this%vtc(:,nPhi)*this%bpc(:,nPhi) -   &
               &              this%vpc(:,nPhi)*this%btc(:,nPhi) )

               this%VxBt(:,nPhi)=  orho1(nR)*or4(nR) * (   &
               &       this%vpc(:,nPhi)*this%brc(:,nPhi) - &
               &       this%vrc(:,nPhi)*this%bpc(:,nPhi) )

               this%VxBp(:,nPhi)=   orho1(nR)*or4(nR) * (   &
               &        this%vrc(:,nPhi)*this%btc(:,nPhi) - &
               &        this%vtc(:,nPhi)*this%brc(:,nPhi) )
            else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary
               this%VxBt(:,nPhi)= or4(nR)*orho1(nR)*this%vpc(:,nPhi)*this%brc(:,nPhi)
               this%VxBp(:,nPhi)=-or4(nR)*orho1(nR)*this%vtc(:,nPhi)*this%brc(:,nPhi)
            else if ( nBc == 2 ) then  ! rigid boundary :
               !-- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):
               this%VxBt(:,nPhi)=or4(nR)*orho1(nR)*this%vpc(:,nPhi)*this%brc(:,nPhi)
               this%VxBp(:,nPhi)= 0.0_cp
            end if  ! boundary ?

         end if ! l_mag_nl ?

         if ( l_anel .and. nBc == 0 ) then
            !------ Get viscous heating
            this%heatTerms(:,nPhi)=ViscHeatFac*or4(nR)*           &
            &                     orho1(nR)*otemp1(nR)*visc(nR)*( &
            &     two*(                     this%dvrdrc(:,nPhi) - & ! (1)
            &     (two*or1(nR)+beta(nR))*this%vrc(:,nPhi) )**2  + &
            &     two*( cosn_theta_E2(:)*   this%vtc(:,nPhi) +    &
            &                               this%dvpdpc(:,nPhi) + &
            &                               this%dvrdrc(:,nPhi) - & ! (2)
            &     or1(nR)*               this%vrc(:,nPhi) )**2  + &
            &     two*(                     this%dvpdpc(:,nPhi) + &
            &           cosn_theta_E2(:)*   this%vtc(:,nPhi) +    & ! (3)
            &     or1(nR)*               this%vrc(:,nPhi) )**2  + &
            &          ( two*               this%dvtdpc(:,nPhi) + &
            &                                 this%cvrc(:,nPhi) - & ! (6)
            &    two*cosn_theta_E2(:)*this%vpc(:,nPhi) )**2  +    &
            &                        O_sin_theta_E2(:) * (        &
            &         ( r(nR)*              this%dvtdrc(:,nPhi) - &
            &           (two+beta(nR)*r(nR))*  this%vtc(:,nPhi) + & ! (4)
            &     or1(nR)*            this%dvrdtc(:,nPhi) )**2  + &
            &         ( r(nR)*              this%dvpdrc(:,nPhi) - &
            &           (two+beta(nR)*r(nR))*  this%vpc(:,nPhi) + & ! (5)
            &     or1(nR)*            this%dvrdpc(:,nPhi) )**2 )- &
            &    two*third*(  beta(nR)*        this%vrc(:,nPhi) )**2 )

            if ( l_mag_nl .and. nR>n_r_LCR ) then
               !------ Get ohmic losses
               this%heatTerms(:,nPhi)=this%heatTerms(:,nPhi)+        &
               &       OhmLossFac*   or2(nR)*otemp1(nR)*lambda(nR)*  &
               &    ( or2(nR)*             this%cbrc(:,nPhi)**2 +    &
               &      O_sin_theta_E2(:)*   this%cbtc(:,nPhi)**2 +    &
               &      O_sin_theta_E2(:)*   this%cbpc(:,nPhi)**2  )
            end if ! if l_mag_nl ?

         end if  ! Viscous heating and Ohmic losses ?

      end do
      !$omp end parallel

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
