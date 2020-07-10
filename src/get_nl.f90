module general_arrays_mod

   implicit none

   private

   type, public, abstract :: general_arrays_t

   end type general_arrays_t

end module general_arrays_mod
!----------------------------------------------------------------------------
module grid_space_arrays_mod

   use general_arrays_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: nrp, n_phi_max, n_theta_max, n_theta_loc, nRstart, &
       &                 nRstop, nThetaStart, nThetaStop, n_r_loc
   use radial_functions, only: or2, orho1, beta, otemp1, visc, r, or3, &
       &                       lambda, or4, or1, alpha0, temp0, opressure0
   use physical_parameters, only: LFfac, n_r_LCR, CorFac, prec_angle,    &
       &                          ThExpNb, ViscHeatFac, oek, po, DissNb, &
       &                          dilution_fac, ra, opr, polind, strat, radratio
   use blocking, only: nfs, sizeThetaB
   use horizontal_data, only: osn2, cosn2, sinTheta, cosTheta, osn1, phi, &
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

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space:
      real(cp), allocatable :: Advr(:,:), Advt(:,:), Advp(:,:)
      real(cp), allocatable :: LFr(:,:), LFt(:,:), LFp(:,:)
      real(cp), allocatable :: PCr(:,:), PCt(:,:), PCp(:,:)
      real(cp), allocatable :: CAr(:,:), CAt(:,:)
      real(cp), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)
      real(cp), allocatable :: VSr(:,:), VSt(:,:), VSp(:,:)
      real(cp), allocatable :: VXir(:,:), VXit(:,:), VXip(:,:)
      real(cp), allocatable :: ViscHeat(:,:), OhmLoss(:,:)

      !---- RMS calculations
      real(cp), allocatable :: Advt2(:,:), Advp2(:,:), LFt2(:,:), LFp2(:,:)
      real(cp), allocatable :: CFt2(:,:), CFp2(:,:), dpdtc(:,:), dpdpc(:,:)
      real(cp), allocatable :: dtVr(:,:), dtVp(:,:), dtVt(:,:)
      real(cp), allocatable :: dpkindrc(:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), allocatable :: vrc(:,:), vtc(:,:), vpc(:,:)
      real(cp), allocatable :: dvrdrc(:,:), dvtdrc(:,:), dvpdrc(:,:)
      real(cp), allocatable :: cvrc(:,:), sc(:,:), drSc(:,:)
      real(cp), allocatable :: dvrdtc(:,:), dvrdpc(:,:)
      real(cp), allocatable :: dvtdpc(:,:), dvpdpc(:,:)
      real(cp), allocatable :: brc(:,:), btc(:,:), bpc(:,:)
      real(cp), allocatable :: cbrc(:,:), cbtc(:,:), cbpc(:,:)
      real(cp), allocatable :: pc(:,:), xic(:,:), cvtc(:,:), cvpc(:,:)
      real(cp), allocatable :: dsdtc(:,:), dsdpc(:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: output
      procedure :: output_nl_input
      procedure :: get_nl

   end type grid_space_arrays_t

   real(cp), allocatable :: vr_old(:,:,:), vt_old(:,:,:), vp_old(:,:,:)

contains

!----------------------------------------------------------------------------
   subroutine initialize(this)

      class(grid_space_arrays_t) :: this

      allocate( this%Advr(nrp,nThetaStart:nThetaStop),     &
      &         this%Advt(nrp,nThetaStart:nThetaStop),     &
      &         this%Advp(nrp,nThetaStart:nThetaStop) )
      allocate( this%LFr(nrp,nThetaStart:nThetaStop),      &
      &         this%LFt(nrp,nThetaStart:nThetaStop),      &
      &         this%LFp(nrp,nThetaStart:nThetaStop) )
      allocate( this%VxBr(nrp,nThetaStart:nThetaStop),     &
      &         this%VxBt(nrp,nThetaStart:nThetaStop),     &
      &         this%VxBp(nrp,nThetaStart:nThetaStop) )
      allocate( this%VSr(nrp,nThetaStart:nThetaStop),      &
      &         this%VSt(nrp,nThetaStart:nThetaStop),      &
      &         this%VSp(nrp,nThetaStart:nThetaStop) )
      allocate( this%ViscHeat(nrp,nThetaStart:nThetaStop), &
      &          this%OhmLoss(nrp,nThetaStart:nThetaStop) )
      bytes_allocated=bytes_allocated + 14*nrp*n_theta_loc*SIZEOF_DEF_REAL

      if ( l_precession ) then
         allocate( this%PCr(nrp,nThetaStart:nThetaStop), &
         &         this%PCt(nrp,nThetaStart:nThetaStop), &
         &         this%PCp(nrp,nThetaStart:nThetaStop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_loc*SIZEOF_DEF_REAL
      end if

      if ( l_centrifuge ) then
         allocate( this%CAr(nrp,nThetaStart:nThetaStop), &
         &         this%CAt(nrp,nThetaStart:nThetaStop) )
         bytes_allocated=bytes_allocated + 2*nrp*n_theta_loc*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXir(nrp,nThetaStart:nThetaStop), &
         &         this%VXit(nrp,nThetaStart:nThetaStop), &
         &         this%VXip(nrp,nThetaStart:nThetaStop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_loc*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(nrp,nThetaStart:nThetaStop), &
      &         this%vtc(nrp,nThetaStart:nThetaStop), &
      &         this%vpc(nrp,nThetaStart:nThetaStop) )
      allocate( this%dvrdrc(nrp,nThetaStart:nThetaStop), &
      &         this%dvtdrc(nrp,nThetaStart:nThetaStop) )
      allocate( this%dvpdrc(nrp,nThetaStart:nThetaStop), &
      &         this%cvrc(nrp,nThetaStart:nThetaStop) )
      allocate( this%dvrdtc(nrp,nThetaStart:nThetaStop), &
      &         this%dvrdpc(nrp,nThetaStart:nThetaStop) )
      allocate( this%dvtdpc(nrp,nThetaStart:nThetaStop), &
      &         this%dvpdpc(nrp,nThetaStart:nThetaStop) )
      allocate( this%brc(nrp,nThetaStart:nThetaStop),    &
      &         this%btc(nrp,nThetaStart:nThetaStop),    &
      &         this%bpc(nrp,nThetaStart:nThetaStop) )
      this%btc=1.0e50_cp
      this%bpc=1.0e50_cp
      allocate( this%cbrc(nrp,nThetaStart:nThetaStop),   &
      &         this%cbtc(nrp,nThetaStart:nThetaStop),   &
      &         this%cbpc(nrp,nThetaStart:nThetaStop) )
      allocate( this%sc(nrp,nThetaStart:nThetaStop),     &
      &         this%drSc(nrp,nThetaStart:nThetaStop) )
      allocate( this%pc(nrp,nThetaStart:nThetaStop) )
      allocate( this%dsdtc(nrp,nThetaStart:nThetaStop),  &
      &         this%dsdpc(nrp,nThetaStart:nThetaStop) )
      bytes_allocated=bytes_allocated + 22*nrp*n_theta_loc*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(nrp,nThetaStart:nThetaStop) )
         bytes_allocated=bytes_allocated + nrp*n_theta_loc*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1) )
      end if

      if ( l_adv_curl ) then
         allocate( this%cvtc(nrp,nThetaStart:nThetaStop), &
         &         this%cvpc(nrp,nThetaStart:nThetaStop) )
         bytes_allocated=bytes_allocated+2*nrp*n_theta_loc*SIZEOF_DEF_REAL
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         allocate ( this%Advt2(nrp,nThetaStart:nThetaStop), &
         &         this%Advp2(nrp,nThetaStart:nThetaStop) )
         allocate ( this%dtVr(nrp,nThetaStart:nThetaStop),  &
         &         this%dtVt(nrp,nThetaStart:nThetaStop),   &
         &         this%dtVp(nrp,nThetaStart:nThetaStop) )
         allocate ( this%LFt2(nrp,nThetaStart:nThetaStop),  &
         &         this%LFp2(nrp,nThetaStart:nThetaStop) )
         allocate ( this%CFt2(nrp,nThetaStart:nThetaStop),  &
         &         this%CFp2(nrp,nThetaStart:nThetaStop) )
         allocate ( this%dpdtc(nrp,nThetaStart:nThetaStop), &
         &         this%dpdpc(nrp,nThetaStart:nThetaStop) )
         bytes_allocated=bytes_allocated + 11*nrp*n_theta_loc*SIZEOF_DEF_REAL

         allocate( vr_old(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         allocate( vp_old(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         allocate( vt_old(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_loc*n_r_loc*&
         &               SIZEOF_DEF_REAL

         this%dtVr(:,:)=0.0_cp
         this%dtVt(:,:)=0.0_cp
         this%dtVp(:,:)=0.0_cp
         vt_old(:,:,:) =0.0_cp
         vr_old(:,:,:) =0.0_cp
         vp_old(:,:,:) =0.0_cp

         if ( l_adv_curl ) then
            allocate ( this%dpkindrc(nrp, nThetaStart:nThetaStop) )
            bytes_allocated=bytes_allocated + nrp*n_theta_loc*SIZEOF_DEF_REAL
         end if
      end if

   end subroutine initialize
   
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_space_arrays_t) :: this

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
   subroutine output(this)

      class(grid_space_arrays_t) :: this

      write(*,"(A,3ES20.12)") "Advr,Advt,Advp = ",sum(this%Advr), &
                                   sum(this%Advt),sum(this%Advp)

   end subroutine output
!----------------------------------------------------------------------------
   subroutine output_nl_input(this)

      class(grid_space_arrays_t) :: this

      write(*,"(A,6ES20.12)") "vr,vt,vp = ",sum(this%vrc),sum(this%vtc), &
                                            sum(this%vpc)

   end subroutine output_nl_input
!----------------------------------------------------------------------------
   subroutine get_nl(this, time, tscheme, nR, nBc, lRmsCalc)
      !
      !  calculates non-linear products in grid-space for radial
      !  level nR and returns them in arrays wnlr1-3, snlr1-3, bnlr1-3
      !
      !  if nBc >0 velocities are zero only the (vxB)
      !  contributions to bnlr2-3 need to be calculated
      !
      !  vr...sr: (input) velocity, magnetic field comp. and derivs, entropy
      !                   on grid points
      !  nR: (input) radial level
      !  i1: (input) range of points in theta for which calculation is done
      !

      class(grid_space_arrays_t) :: this

      !-- Input of variables:
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: nR
      logical,             intent(in) :: lRmsCalc
      integer,             intent(in) :: nBc

      !-- Local variables:
      integer :: n_th, nThStart, nThStop
      real(cp) :: cnt, rsnt, snt, posnalp, O_dt

      !$omp parallel default(shared) private(nThStart,nThStop) &
      !$omp private(n_th,cnt,snt,rsnt,posnalp)
      nThStart=nThetaStart; nThStop=nThetaStop
      call get_openmp_blocks(nThStart,nThStop)

      do n_th=nThStart,nThStop

         if ( l_mag_LF .and. (nBc == 0 .or. lRmsCalc) .and. nR>n_r_LCR ) then
            !------ Get the Lorentz force:
            !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
            this%LFr(:,n_th)=  LFfac*O_sin_theta_E2(n_th) * (   &
            &        this%cbtc(:,n_th)*this%bpc(:,n_th) -       &
            &        this%cbpc(:,n_th)*this%btc(:,n_th) )

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            this%LFt(:,n_th)=  LFfac*or4(nR)*O_sin_theta_E2(n_th) * ( &
            &        this%cbpc(:,n_th)*this%brc(:,n_th) -             &
            &        this%cbrc(:,n_th)*this%bpc(:,n_th) )

            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            this%LFp(:,n_th)=  LFfac*or4(nR)*O_sin_theta_E2(n_th) * ( &
            &        this%cbrc(:,n_th)*this%btc(:,n_th) -             &
            &        this%cbtc(:,n_th)*this%brc(:,n_th) )
         end if      ! Lorentz force required ?

         if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

            if ( l_adv_curl ) then ! Advection is \curl{u} \times u
               this%Advr(:,n_th)=  - O_sin_theta_E2(n_th) * (    &
               &        this%cvtc(:,n_th)*this%vpc(:,n_th) -     &
               &        this%cvpc(:,n_th)*this%vtc(:,n_th) )

               this%Advt(:,n_th)= -or4(nR)*O_sin_theta_E2(n_th) * ( &
               &        this%cvpc(:,n_th)*this%vrc(:,n_th) -        &
               &        this%cvrc(:,n_th)*this%vpc(:,n_th) )

               this%Advp(:,n_th)= -or4(nR)*O_sin_theta_E2(n_th) * ( &
               &        this%cvrc(:,n_th)*this%vtc(:,n_th) -        &
               &        this%cvtc(:,n_th)*this%vrc(:,n_th) )
            else ! Advection is u\grad u
               !------ Get Advection:
               this%Advr(:,n_th)=          -or2(nR)*orho1(nR) * (  &
               &                                this%vrc(:,n_th) * &
               &                     (       this%dvrdrc(:,n_th) - &
               &    ( two*or1(nR)+beta(nR) )*this%vrc(:,n_th) ) +  &
               &                      O_sin_theta_E2(n_th) * (     &
               &                                this%vtc(:,n_th) * &
               &                     (       this%dvrdtc(:,n_th) - &
               &                  r(nR)*      this%vtc(:,n_th) ) + &
               &                                this%vpc(:,n_th) * &
               &                     (       this%dvrdpc(:,n_th) - &
               &                    r(nR)*      this%vpc(:,n_th) ) ) )

               this%Advt(:,n_th)=or4(nR)*O_sin_theta_E2(n_th)*orho1(nR) * (  &
               &                                         -this%vrc(:,n_th) * &
               &                                   (   this%dvtdrc(:,n_th) - &
               &                             beta(nR)*this%vtc(:,n_th) )   + &
               &                                          this%vtc(:,n_th) * &
               &                    ( cosn_theta_E2(n_th)*this%vtc(:,n_th) + &
               &                                       this%dvpdpc(:,n_th) + &
               &                                   this%dvrdrc(:,n_th) )   + &
               &                                          this%vpc(:,n_th) * &
               &                    ( cosn_theta_E2(n_th)*this%vpc(:,n_th) - &
               &                                       this%dvtdpc(:,n_th) )  )

               this%Advp(:,n_th)= or4(nR)*O_sin_theta_E2(n_th)*orho1(nR) * (  &
               &                                          -this%vrc(:,n_th) * &
               &                                      ( this%dvpdrc(:,n_th) - &
               &                              beta(nR)*this%vpc(:,n_th) )   - &
               &                                           this%vtc(:,n_th) * &
               &                                      ( this%dvtdpc(:,n_th) + &
               &                                      this%cvrc(:,n_th) )   - &
               &                     this%vpc(:,n_th) * this%dvpdpc(:,n_th) )
            end if

         end if  ! Navier-Stokes nonlinear advection term ?

         if ( l_heat_nl .and. nBc == 0 ) then
            !------ Get V S, the divergence of it is entropy advection:
            this%VSr(:,n_th)=this%vrc(:,n_th)*this%sc(:,n_th)
            this%VSt(:,n_th)=or2(nR)*this%vtc(:,n_th)*this%sc(:,n_th)
            this%VSp(:,n_th)=or2(nR)*this%vpc(:,n_th)*this%sc(:,n_th)
         end if     ! heat equation required ?

         if ( l_chemical_conv .and. nBc == 0 ) then
            this%VXir(:,n_th)=this%vrc(:,n_th)*this%xic(:,n_th)
            this%VXit(:,n_th)=or2(nR)*this%vtc(:,n_th)*this%xic(:,n_th)
            this%VXip(:,n_th)=or2(nR)*this%vpc(:,n_th)*this%xic(:,n_th)
         end if     ! chemical composition equation required ?

         if ( l_precession .and. nBc == 0 ) then
            posnalp=-two*oek*po*sin(prec_angle)*O_sin_theta(n_th)
            cnt=cosTheta(n_th)
            this%PCr(:,n_th)=posnalp*r(nR)*(cos(oek*time+phi(:))*      &
            &                this%vpc(:,n_th)*cnt+sin(oek*time+phi(:))*&
            &                this%vtc(:,n_th))
            this%PCt(:,n_th)=   -posnalp*or2(nR)*(                    &
            &               cos(oek*time+phi(:))*this%vpc(:,n_th)     &
            &      +sin(oek*time+phi(:))*or1(nR)*this%vrc(:,n_th) )
            this%PCp(:,n_th)= posnalp*cos(oek*time+phi(:))*or2(nR)*(   &
            &                 this%vtc(:,n_th)-or1(nR)*this%vrc(:,n_th)*cnt)
         end if ! precession term required ?

         if ( l_centrifuge .and. nBc ==0 ) then
            snt=sinTheta(n_th)
            cnt=cosTheta(n_th)
            rsnt=r(nR)*snt
            !if ( l_anel ) then
            !   this%CAr(:,n_th) = dilution_fac*rsnt*snt* &
            !   &       ( -ra*opr*this%sc(:,n_th) )
            !   !-- neglect pressure contribution
            !   !& + polind*DissNb*oek*opressure0(nR)*this%pc(:,n_th) )
            !   this%CAt(:,n_th) = dilution_fac*rsnt*cnt* &
            !   &       ( -ra*opr*this%sc(:,n_th) )
            !   !-- neglect pressure contribution
            !   !& + polind*DissNb*oek*opressure0(nR)*this%pc(:,n_th) )
            !else
            this%CAr(:,n_th) = -dilution_fac*rsnt*snt*ra*opr* &
            &                       this%sc(:,n_th)
            this%CAt(:,n_th) = -dilution_fac*rsnt*cnt*ra*opr* &
            &                       this%sc(:,n_th)
            !end if
         end if ! centrifuge

         if ( l_mag_nl ) then

            if ( nBc == 0 .and. nR>n_r_LCR ) then
               !------ Get (V x B) , the curl of this is the dynamo term:
               this%VxBr(:,n_th)=  orho1(nR)*O_sin_theta_E2(n_th) * (  &
               &              this%vtc(:,n_th)*this%bpc(:,n_th) -      &
               &              this%vpc(:,n_th)*this%btc(:,n_th) )

               this%VxBt(:,n_th)=  orho1(nR)*or4(nR) * (   &
               &       this%vpc(:,n_th)*this%brc(:,n_th) - &
               &       this%vrc(:,n_th)*this%bpc(:,n_th) )

               this%VxBp(:,n_th)=   orho1(nR)*or4(nR) * (   &
               &        this%vrc(:,n_th)*this%btc(:,n_th) - &
               &        this%vtc(:,n_th)*this%brc(:,n_th) )
            else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary
               this%VxBt(:,n_th)= or4(nR)*orho1(nR)*this%vpc(:,n_th)*this%brc(:,n_th)
               this%VxBp(:,n_th)=-or4(nR)*orho1(nR)*this%vtc(:,n_th)*this%brc(:,n_th)
            else if ( nBc == 2 ) then  ! rigid boundary :
               !-- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):
               this%VxBt(:,n_th)=or4(nR)*orho1(nR)*this%vpc(:,n_th)*this%brc(:,n_th)
               this%VxBp(:,n_th)= 0.0_cp
            end if  ! boundary ?

         end if ! l_mag_nl ?

         if ( l_anel .and. nBc == 0 ) then
            !------ Get viscous heating
            this%ViscHeat(:,n_th)=      or4(nR)*                  &
            &                     orho1(nR)*otemp1(nR)*visc(nR)*( &
            &     two*(                     this%dvrdrc(:,n_th) - & ! (1)
            &     (two*or1(nR)+beta(nR))*this%vrc(:,n_th) )**2  + &
            &     two*( cosn_theta_E2(n_th)*   this%vtc(:,n_th) + &
            &                               this%dvpdpc(:,n_th) + &
            &                               this%dvrdrc(:,n_th) - & ! (2)
            &     or1(nR)*               this%vrc(:,n_th) )**2  + &
            &     two*(                     this%dvpdpc(:,n_th) + &
            &           cosn_theta_E2(n_th)*   this%vtc(:,n_th) + & ! (3)
            &     or1(nR)*               this%vrc(:,n_th) )**2  + &
            &          ( two*               this%dvtdpc(:,n_th) + &
            &                                 this%cvrc(:,n_th) - & ! (6)
            &    two*cosn_theta_E2(n_th)*this%vpc(:,n_th) )**2  + &
            &                        O_sin_theta_E2(n_th) * (     &
            &         ( r(nR)*              this%dvtdrc(:,n_th) - &
            &           (two+beta(nR)*r(nR))*  this%vtc(:,n_th) + & ! (4)
            &     or1(nR)*            this%dvrdtc(:,n_th) )**2  + &
            &         ( r(nR)*              this%dvpdrc(:,n_th) - &
            &           (two+beta(nR)*r(nR))*  this%vpc(:,n_th) + & ! (5)
            &     or1(nR)*            this%dvrdpc(:,n_th) )**2 )- &
            &    two*third*(  beta(nR)*        this%vrc(:,n_th) )**2 )

            if ( l_mag_nl .and. nR>n_r_LCR ) then
               !------ Get ohmic losses
               this%OhmLoss(:,n_th)= or2(nR)*otemp1(nR)*lambda(nR)*  &
               &    ( or2(nR)*                this%cbrc(:,n_th)**2 + &
               &      O_sin_theta_E2(n_th)*   this%cbtc(:,n_th)**2 + &
               &      O_sin_theta_E2(n_th)*   this%cbpc(:,n_th)**2  )
            end if ! if l_mag_nl ?

         end if  ! Viscous heating and Ohmic losses ?

         if ( lRmsCalc ) then
            snt=sinTheta(n_th)
            cnt=cosTheta(n_th)
            rsnt=r(nR)*snt
            this%dpdtc(:,n_th)=this%dpdtc(:,n_th)*or1(nR)
            this%dpdpc(:,n_th)=this%dpdpc(:,n_th)*or1(nR)
            this%CFt2(:,n_th)=-two*CorFac*cnt*this%vpc(:,n_th)*or1(nR)
            this%CFp2(:,n_th)= two*CorFac*snt* (cnt*this%vtc(:,n_th)/rsnt +   &
            &                     or2(nR)*snt*this%vrc(:,n_th) )
            if ( l_conv_nl ) then
               this%Advt2(:,n_th)=rsnt*snt*this%Advt(:,n_th)
               this%Advp2(:,n_th)=rsnt*snt*this%Advp(:,n_th)
            end if
            if ( l_mag_LF .and. nR > n_r_LCR ) then
               this%LFt2(:,n_th)=rsnt*snt*this%LFt(:,n_th)
               this%LFp2(:,n_th)=rsnt*snt*this%LFp(:,n_th)
            end if

            if ( l_adv_curl ) then
               this%dpdtc(:,n_th)=this%dpdtc(:,n_th)-or3(nR)*( or2(nR)*   &
               &            this%vrc(:,n_th)*this%dvrdtc(:,n_th) -        &
               &            this%vtc(:,n_th)*(this%dvrdrc(:,n_th)+        &
               &            this%dvpdpc(:,n_th)+cosn_theta_E2(n_th) *     &
               &            this%vtc(:,n_th))+ this%vpc(:,n_th)*(         &
               &            this%cvrc(:,n_th)+this%dvtdpc(:,n_th)-        &
               &            cosn_theta_E2(n_th)*this%vpc(:,n_th)) )
               this%dpdpc(:,n_th)=this%dpdpc(:,n_th)- or3(nR)*( or2(nR)*  &
               &            this%vrc(:,n_th)*this%dvrdpc(:,n_th) +        &
               &            this%vtc(:,n_th)*this%dvtdpc(:,n_th) +        &
               &            this%vpc(:,n_th)*this%dvpdpc(:,n_th) )
               if ( l_conv_nl ) then
                  this%Advt2(:,n_th)=this%Advt2(:,n_th)-or3(nR)*( or2(nR)*&
                  &            this%vrc(:,n_th)*this%dvrdtc(:,n_th) -     &
                  &            this%vtc(:,n_th)*(this%dvrdrc(:,n_th)+     &
                  &            this%dvpdpc(:,n_th)+cosn_theta_E2(n_th) *  &
                  &            this%vtc(:,n_th))+ this%vpc(:,n_th)*(      &
                  &            this%cvrc(:,n_th)+this%dvtdpc(:,n_th)-     &
                  &            cosn_theta_E2(n_th)*this%vpc(:,n_th)) )
                  this%Advp2(:,n_th)=this%Advp2(:,n_th)-or3(nR)*( or2(nR)*&
                  &            this%vrc(:,n_th)*this%dvrdpc(:,n_th) +     &
                  &            this%vtc(:,n_th)*this%dvtdpc(:,n_th) +     &
                  &            this%vpc(:,n_th)*this%dvpdpc(:,n_th) )
               end if

               !- dpkin/dr = 1/2 d (u^2) / dr = ur*dur/dr+ut*dut/dr+up*dup/dr
               this%dpkindrc(:,n_th)=or4(nR)*this%vrc(:,n_th)*(         &
               &                         this%dvrdrc(:,n_th)-           &
               &                         two*or1(nR)*this%vrc(:,n_th)) +&
               &                         or2(nR)*O_sin_theta_E2(n_th)*( &
               &                                 this%vtc(:,n_th)*(     &
               &                         this%dvtdrc(:,n_th)-           &
               &                         or1(nR)*this%vtc(:,n_th) ) +   &
               &                                 this%vpc(:,n_th)*(     &
               &                         this%dvpdrc(:,n_th)-           &
               &                         or1(nR)*this%vpc(:,n_th) ) )
            end if
         end if

         if ( l_RMS .and. tscheme%istage == 1 ) then
            O_dt = 1.0_cp/tscheme%dt(1)
            this%dtVr(:,n_th)=O_dt*or2(nR)*(this%vrc(:,n_th)-vr_old(:,n_th,nR))
            this%dtVt(:,n_th)=O_dt*or1(nR)*(this%vtc(:,n_th)-vt_old(:,n_th,nR))
            this%dtVp(:,n_th)=O_dt*or1(nR)*(this%vpc(:,n_th)-vp_old(:,n_th,nR))

            vr_old(:,n_th,nR)=this%vrc(:,n_th)
            vt_old(:,n_th,nR)=this%vtc(:,n_th)
            vp_old(:,n_th,nR)=this%vpc(:,n_th)
         end if

      end do
      !$omp end parallel

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
