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
   use horizontal_data, only: osn2, cosn2, sinTheta, cosTheta, osn1, phi
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
   real(cp), allocatable :: vr_old_dist(:,:,:), vt_old_dist(:,:,:), vp_old_dist(:,:,:)

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

         allocate( vr_old_dist(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         allocate( vp_old_dist(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         allocate( vt_old_dist(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_loc*n_r_loc*&
         &               SIZEOF_DEF_REAL

         this%dtVr(:,:)=0.0_cp
         this%dtVt(:,:)=0.0_cp
         this%dtVp(:,:)=0.0_cp
         vt_old_dist(:,:,:) =0.0_cp
         vr_old_dist(:,:,:) =0.0_cp
         vp_old_dist(:,:,:) =0.0_cp

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
         if (allocated(vr_old_dist)) deallocate ( vr_old_dist )
         if (allocated(vt_old_dist)) deallocate ( vt_old_dist )
         if (allocated(vp_old_dist)) deallocate ( vp_old_dist )
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
      integer :: n_th, nThetaNHS, n_phi, nThStart, nThStop
      real(cp) :: or4sn2, csn2, cnt, rsnt, snt, posnalp, O_dt

      !$omp parallel default(shared) private(nThStart,nThStop) &
      !$omp private(n_th,nThetaNHS,n_phi) &
      !$omp private(or4sn2,csn2,cnt,snt,rsnt,posnalp)
      nThStart=nThetaStart; nThStop=nThetaStop
      call get_openmp_blocks(nThStart,nThStop)


      if ( l_mag_LF .and. (nBc == 0 .or. lRmsCalc) .and. nR>n_r_LCR ) then
         !------ Get the Lorentz force:
         do n_th=nThStart,nThStop

            nThetaNHS=(n_th+1)/2
            or4sn2   =or4(nR)*osn2(nThetaNHS)

            do n_phi=1,n_phi_max
               !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
               this%LFr(n_phi,n_th)=  LFfac*osn2(nThetaNHS) * (        &
               &        this%cbtc(n_phi,n_th)*this%bpc(n_phi,n_th) - &
               &        this%cbpc(n_phi,n_th)*this%btc(n_phi,n_th) )
            end do

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            do n_phi=1,n_phi_max
               this%LFt(n_phi,n_th)=           LFfac*or4sn2 * (        &
               &        this%cbpc(n_phi,n_th)*this%brc(n_phi,n_th) - &
               &        this%cbrc(n_phi,n_th)*this%bpc(n_phi,n_th) )
            end do

            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            do n_phi=1,n_phi_max
               this%LFp(n_phi,n_th)=           LFfac*or4sn2 * (        &
               &        this%cbrc(n_phi,n_th)*this%btc(n_phi,n_th) - &
               &        this%cbtc(n_phi,n_th)*this%brc(n_phi,n_th) )
            end do

         end do   ! theta loop
      end if      ! Lorentz force required ?

      if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

         if ( l_adv_curl ) then ! Advection is \curl{u} \times u

            do n_th=nThStart,nThStop ! loop over theta points in block

               nThetaNHS=(n_th+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)

               do n_phi=1,n_phi_max
                  this%Advr(n_phi,n_th)=       - osn2(nThetaNHS) * (    &
                  &        this%cvtc(n_phi,n_th)*this%vpc(n_phi,n_th) - &
                  &        this%cvpc(n_phi,n_th)*this%vtc(n_phi,n_th) )

                  this%Advt(n_phi,n_th)=        -        or4sn2 * (     &
                  &        this%cvpc(n_phi,n_th)*this%vrc(n_phi,n_th) - &
                  &        this%cvrc(n_phi,n_th)*this%vpc(n_phi,n_th) )

                  this%Advp(n_phi,n_th)=        -        or4sn2 * (     &
                  &        this%cvrc(n_phi,n_th)*this%vtc(n_phi,n_th) - &
                  &        this%cvtc(n_phi,n_th)*this%vrc(n_phi,n_th) )
               end do

            end do

         else ! Advection is u\grad u

            !------ Get Advection:
            do n_th=nThStart,nThStop ! loop over theta points in block
               nThetaNHS=(n_th+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               csn2     =cosn2(nThetaNHS)
               if ( mod(n_th,2) == 0 ) csn2=-csn2 ! South, odd function in theta

               do n_phi=1,n_phi_max
                  this%Advr(n_phi,n_th)=          -or2(nR)*orho1(nR) * (  &
                  &                                this%vrc(n_phi,n_th) * &
                  &                     (       this%dvrdrc(n_phi,n_th) - &
                  &    ( two*or1(nR)+beta(nR) )*this%vrc(n_phi,n_th) ) +  &
                  &                               osn2(nThetaNHS) * (     &
                  &                                this%vtc(n_phi,n_th) * &
                  &                     (       this%dvrdtc(n_phi,n_th) - &
                  &                  r(nR)*      this%vtc(n_phi,n_th) ) + &
                  &                                this%vpc(n_phi,n_th) * &
                  &                     (       this%dvrdpc(n_phi,n_th) - &
                  &                    r(nR)*      this%vpc(n_phi,n_th) ) ) )
               end do

               do n_phi=1,n_phi_max
                  this%Advt(n_phi,n_th)=         or4sn2*orho1(nR) * (  &
                  &                            -this%vrc(n_phi,n_th) * &
                  &                      (   this%dvtdrc(n_phi,n_th) - &
                  &                beta(nR)*this%vtc(n_phi,n_th) )   + &
                  &                             this%vtc(n_phi,n_th) * &
                  &                      ( csn2*this%vtc(n_phi,n_th) + &
                  &                          this%dvpdpc(n_phi,n_th) + &
                  &                      this%dvrdrc(n_phi,n_th) )   + &
                  &                             this%vpc(n_phi,n_th) * &
                  &                      ( csn2*this%vpc(n_phi,n_th) - &
                  &                          this%dvtdpc(n_phi,n_th) )  )
               end do

               do n_phi=1,n_phi_max
                  this%Advp(n_phi,n_th)=         or4sn2*orho1(nR) * (  &
                  &                            -this%vrc(n_phi,n_th) * &
                  &                        ( this%dvpdrc(n_phi,n_th) - &
                  &                beta(nR)*this%vpc(n_phi,n_th) )   - &
                  &                             this%vtc(n_phi,n_th) * &
                  &                        ( this%dvtdpc(n_phi,n_th) + &
                  &                        this%cvrc(n_phi,n_th) )   - &
                  &       this%vpc(n_phi,n_th) * this%dvpdpc(n_phi,n_th) )
               end do

            end do ! theta loop
         end if

      end if  ! Navier-Stokes nonlinear advection term ?

      if ( l_heat_nl .and. nBc == 0 ) then
         !------ Get V S, the divergence of it is entropy advection:
         do n_th=nThStart,nThStop
            do n_phi=1,n_phi_max     ! calculate v*s components
               this%VSr(n_phi,n_th)= &
               &    this%vrc(n_phi,n_th)*this%sc(n_phi,n_th)
               this%VSt(n_phi,n_th)= &
               &    or2(nR)*this%vtc(n_phi,n_th)*this%sc(n_phi,n_th)
               this%VSp(n_phi,n_th)= &
               &    or2(nR)*this%vpc(n_phi,n_th)*this%sc(n_phi,n_th)
            end do
         end do  ! theta loop
      end if     ! heat equation required ?

      if ( l_chemical_conv .and. nBc == 0 ) then
         !------ Get V S, the divergence of it is the advection of chem comp:
         do n_th=nThStart,nThStop
            do n_phi=1,n_phi_max     ! calculate v*s components
               this%VXir(n_phi,n_th)= &
               &    this%vrc(n_phi,n_th)*this%xic(n_phi,n_th)
               this%VXit(n_phi,n_th)= &
               &    or2(nR)*this%vtc(n_phi,n_th)*this%xic(n_phi,n_th)
               this%VXip(n_phi,n_th)= &
               &    or2(nR)*this%vpc(n_phi,n_th)*this%xic(n_phi,n_th)
            end do
         end do  ! theta loop
      end if     ! chemical composition equation required ?

      if ( l_precession .and. nBc == 0 ) then

         do n_th=nThStart,nThStop
            nThetaNHS=(n_th+1)/2
            posnalp=-two*oek*po*sin(prec_angle)*osn1(nThetaNHS)
            cnt=cosTheta(n_th)
            do n_phi=1,n_phi_max
               this%PCr(n_phi,n_th)=posnalp*r(nR)*(cos(oek*time+phi(n_phi))* &
               &                                  this%vpc(n_phi,n_th)*cnt  &
               &            +sin(oek*time+phi(n_phi))*this%vtc(n_phi,n_th))
               this%PCt(n_phi,n_th)=   -posnalp*or2(nR)*(                   &
               &               cos(oek*time+phi(n_phi))*this%vpc(n_phi,n_th) &
               &      +sin(oek*time+phi(n_phi))*or1(nR)*this%vrc(n_phi,n_th) )
               this%PCp(n_phi,n_th)= posnalp*cos(oek*time+phi(n_phi))*       &
               &              or2(nR)*(      this%vtc(n_phi,n_th)-          &
               &                     or1(nR)*this%vrc(n_phi,n_th)*cnt)
            end do
         end do ! theta loop
      end if ! precession term required ?

      if ( l_centrifuge .and. nBc ==0 ) then
         do n_th=nThStart,nThStop
            nThetaNHS=(n_th+1)/2
            snt=sinTheta(n_th)
            cnt=cosTheta(n_th)
            rsnt=r(nR)*snt
            do n_phi=1,n_phi_max
               if ( l_anel ) then
                  this%CAr(n_phi,n_th) = dilution_fac*rsnt*snt* &
                  &       ( -ra*opr*this%sc(n_phi,n_th) )
                  !-- neglect pressure contribution
                  !& + polind*DissNb*oek*opressure0(nR)*this%pc(n_phi,n_th) )
                  this%CAt(n_phi,n_th) = dilution_fac*rsnt*cnt* &
                  &       ( -ra*opr*this%sc(n_phi,n_th) )
                  !-- neglect pressure contribution
                  !& + polind*DissNb*oek*opressure0(nR)*this%pc(n_phi,n_th) )
               else
                  this%CAr(n_phi,n_th) = -dilution_fac*rsnt*snt*ra*opr* &
                  &                       this%sc(n_phi,n_th)
                  this%CAt(n_phi,n_th) = -dilution_fac*rsnt*cnt*ra*opr* &
                  &                       this%sc(n_phi,n_th)
               end if
            end do ! phi loop
         end do ! theta loop
      end if ! centrifuge

      if ( l_mag_nl ) then

         if ( nBc == 0 .and. nR>n_r_LCR ) then

            !------ Get (V x B) , the curl of this is the dynamo term:
            do n_th=nThStart,nThStop
               nThetaNHS=(n_th+1)/2

               do n_phi=1,n_phi_max
                  this%VxBr(n_phi,n_th)=  orho1(nR)*osn2(nThetaNHS) * (      &
                  &              this%vtc(n_phi,n_th)*this%bpc(n_phi,n_th) - &
                  &              this%vpc(n_phi,n_th)*this%btc(n_phi,n_th) )
               end do

               do n_phi=1,n_phi_max
                  this%VxBt(n_phi,n_th)=  orho1(nR)*or4(nR) * (       &
                  &       this%vpc(n_phi,n_th)*this%brc(n_phi,n_th) - &
                  &       this%vrc(n_phi,n_th)*this%bpc(n_phi,n_th) )
               end do

               do n_phi=1,n_phi_max
                  this%VxBp(n_phi,n_th)=   orho1(nR)*or4(nR) * (       &
                  &        this%vrc(n_phi,n_th)*this%btc(n_phi,n_th) - &
                  &        this%vtc(n_phi,n_th)*this%brc(n_phi,n_th) )
               end do
            end do   ! theta loop

         else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary

            do n_th=nThStart,nThStop
               do n_phi=1,n_phi_max
                  this%VxBt(n_phi,n_th)=  or4(nR) * orho1(nR) * &
                  &        this%vpc(n_phi,n_th)*this%brc(n_phi,n_th)
                  this%VxBp(n_phi,n_th)= -or4(nR) * orho1(nR) * &
                  &        this%vtc(n_phi,n_th)*this%brc(n_phi,n_th)
               end do
            end do

         else if ( nBc == 2 ) then  ! rigid boundary :

            !----- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):

            do n_th=nThStart,nThStop
               do n_phi=1,n_phi_max
                  this%VxBt(n_phi,n_th)= or4(nR) * orho1(nR) * &
                  &    this%vpc(n_phi,n_th)*this%brc(n_phi,n_th)

                  this%VxBp(n_phi,n_th)= 0.0_cp
               end do
            end do

         end if  ! boundary ?

      end if ! l_mag_nl ?

      if ( l_anel .and. nBc == 0 ) then
         !------ Get viscous heating
         do n_th=nThStart,nThStop ! loop over theta points in block
            nThetaNHS=(n_th+1)/2
            csn2     =cosn2(nThetaNHS)
            if ( mod(n_th,2) == 0 ) csn2=-csn2 ! South, odd function in theta

            do n_phi=1,n_phi_max
               this%ViscHeat(n_phi,n_th)=      or4(nR)*                  &
               &                     orho1(nR)*otemp1(nR)*visc(nR)*(     &
               &     two*(                     this%dvrdrc(n_phi,n_th) - & ! (1)
               &     (two*or1(nR)+beta(nR))*this%vrc(n_phi,n_th) )**2  + &
               &     two*( csn2*                  this%vtc(n_phi,n_th) + &
               &                               this%dvpdpc(n_phi,n_th) + &
               &                               this%dvrdrc(n_phi,n_th) - & ! (2)
               &     or1(nR)*               this%vrc(n_phi,n_th) )**2  + &
               &     two*(                     this%dvpdpc(n_phi,n_th) + &
               &           csn2*                  this%vtc(n_phi,n_th) + & ! (3)
               &     or1(nR)*               this%vrc(n_phi,n_th) )**2  + &
               &          ( two*               this%dvtdpc(n_phi,n_th) + &
               &                                 this%cvrc(n_phi,n_th) - & ! (6)
               &      two*csn2*             this%vpc(n_phi,n_th) )**2  + &
               &                                 osn2(nThetaNHS) * (     &
               &         ( r(nR)*              this%dvtdrc(n_phi,n_th) - &
               &           (two+beta(nR)*r(nR))*  this%vtc(n_phi,n_th) + & ! (4)
               &     or1(nR)*            this%dvrdtc(n_phi,n_th) )**2  + &
               &         ( r(nR)*              this%dvpdrc(n_phi,n_th) - &
               &           (two+beta(nR)*r(nR))*  this%vpc(n_phi,n_th) + & ! (5)
               &     or1(nR)*            this%dvrdpc(n_phi,n_th) )**2 )- &
               &    two*third*(  beta(nR)*        this%vrc(n_phi,n_th) )**2 )
            end do
         end do ! theta loop

         if ( l_mag_nl .and. nR>n_r_LCR ) then
            !------ Get ohmic losses
            do n_th=nThStart,nThStop ! loop over theta points in block
               nThetaNHS=(n_th+1)/2
               do n_phi=1,n_phi_max
                  this%OhmLoss(n_phi,n_th)= or2(nR)*otemp1(nR)*lambda(nR)*  &
                  &    ( or2(nR)*                this%cbrc(n_phi,n_th)**2 + &
                  &      osn2(nThetaNHS)*        this%cbtc(n_phi,n_th)**2 + &
                  &      osn2(nThetaNHS)*        this%cbpc(n_phi,n_th)**2  )
               end do
            end do ! theta loop

         end if ! if l_mag_nl ?

      end if  ! Viscous heating and Ohmic losses ?

      if ( lRmsCalc ) then
         do n_th=nThStart,nThStop ! loop over theta points in block
            snt=sinTheta(n_th)
            cnt=cosTheta(n_th)
            rsnt=r(nR)*snt
            do n_phi=1,n_phi_max
               this%dpdtc(n_phi,n_th)=this%dpdtc(n_phi,n_th)*or1(nR)
               this%dpdpc(n_phi,n_th)=this%dpdpc(n_phi,n_th)*or1(nR)
               this%CFt2(n_phi,n_th)=-two*CorFac*cnt*this%vpc(n_phi,n_th)*or1(nR)
               this%CFp2(n_phi,n_th)= two*CorFac*snt* (                &
               &                     cnt*this%vtc(n_phi,n_th)/rsnt +   &
               &                     or2(nR)*snt*this%vrc(n_phi,n_th) )
               if ( l_conv_nl ) then
                  this%Advt2(n_phi,n_th)=rsnt*snt*this%Advt(n_phi,n_th)
                  this%Advp2(n_phi,n_th)=rsnt*snt*this%Advp(n_phi,n_th)
               end if
               if ( l_mag_LF .and. nR > n_r_LCR ) then
                  this%LFt2(n_phi,n_th)=rsnt*snt*this%LFt(n_phi,n_th)
                  this%LFp2(n_phi,n_th)=rsnt*snt*this%LFp(n_phi,n_th)
               end if
            end do
         end do

         if ( l_adv_curl ) then
            do n_th=nThStart,nThStop ! loop over theta points in block
               nThetaNHS=(n_th+1)/2
               csn2     =cosn2(nThetaNHS)
               if ( mod(n_th,2) == 0 ) csn2=-csn2 ! South, odd function in theta
               do n_phi=1,n_phi_max
                  this%dpdtc(n_phi,n_th)=this%dpdtc(n_phi,n_th)-              &
                  &                      or3(nR)*( or2(nR)*                   &
                  &            this%vrc(n_phi,n_th)*this%dvrdtc(n_phi,n_th) - &
                  &            this%vtc(n_phi,n_th)*(this%dvrdrc(n_phi,n_th)+ &
                  &            this%dvpdpc(n_phi,n_th)+csn2 *                 &
                  &            this%vtc(n_phi,n_th))+ this%vpc(n_phi,n_th)*(  &
                  &            this%cvrc(n_phi,n_th)+this%dvtdpc(n_phi,n_th)- &
                  &            csn2*this%vpc(n_phi,n_th)) )
                  this%dpdpc(n_phi,n_th)=this%dpdpc(n_phi,n_th)-              &
                  &                         or3(nR)*( or2(nR)*                &
                  &            this%vrc(n_phi,n_th)*this%dvrdpc(n_phi,n_th) + &
                  &            this%vtc(n_phi,n_th)*this%dvtdpc(n_phi,n_th) + &
                  &            this%vpc(n_phi,n_th)*this%dvpdpc(n_phi,n_th) )
                  if ( l_conv_nl ) then
                     this%Advt2(n_phi,n_th)=this%Advt2(n_phi,n_th)-              &
                     &                      or3(nR)*( or2(nR)*                   &
                     &            this%vrc(n_phi,n_th)*this%dvrdtc(n_phi,n_th) - &
                     &            this%vtc(n_phi,n_th)*(this%dvrdrc(n_phi,n_th)+ &
                     &            this%dvpdpc(n_phi,n_th)+csn2 *                 &
                     &            this%vtc(n_phi,n_th))+ this%vpc(n_phi,n_th)*(  &
                     &            this%cvrc(n_phi,n_th)+this%dvtdpc(n_phi,n_th)- &
                     &            csn2*this%vpc(n_phi,n_th)) )
                     this%Advp2(n_phi,n_th)=this%Advp2(n_phi,n_th)-              &
                     &                      or3(nR)*( or2(nR)*                   &
                     &            this%vrc(n_phi,n_th)*this%dvrdpc(n_phi,n_th) + &
                     &            this%vtc(n_phi,n_th)*this%dvtdpc(n_phi,n_th) + &
                     &            this%vpc(n_phi,n_th)*this%dvpdpc(n_phi,n_th) )
                  end if

                  !- dpkin/dr = 1/2 d (u^2) / dr = ur*dur/dr+ut*dut/dr+up*dup/dr
                  this%dpkindrc(n_phi,n_th)=or4(nR)*this%vrc(n_phi,n_th)*(     &
                  &                         this%dvrdrc(n_phi,n_th)-           &
                  &                         two*or1(nR)*this%vrc(n_phi,n_th)) +&
                  &                         or2(nR)*osn2(nThetaNHS)*(          &
                  &                                 this%vtc(n_phi,n_th)*(     &
                  &                         this%dvtdrc(n_phi,n_th)-           &
                  &                         or1(nR)*this%vtc(n_phi,n_th) ) +   &
                  &                                 this%vpc(n_phi,n_th)*(     &
                  &                         this%dvpdrc(n_phi,n_th)-           &
                  &                         or1(nR)*this%vpc(n_phi,n_th) ) )
               end do
            end do
         end if
      end if

      if ( l_RMS .and. tscheme%istage == 1 ) then
         O_dt = 1.0_cp/tscheme%dt(1)
         do n_th=nThStart,nThStop
            do n_phi=1,n_phi_max
               this%dtVr(n_phi,n_th)=O_dt*or2(nR)*(this%vrc(n_phi,n_th)- &
               &                             vr_old_dist(n_phi,n_th,nR))
               this%dtVt(n_phi,n_th)=O_dt*or1(nR)*(this%vtc(n_phi,n_th)- &
               &                             vt_old_dist(n_phi,n_th,nR))
               this%dtVp(n_phi,n_th)=O_dt*or1(nR)*(this%vpc(n_phi,n_th)- &
               &                             vp_old_dist(n_phi,n_th,nR))

               vr_old_dist(n_phi,n_th,nR)=this%vrc(n_phi,n_th)
               vt_old_dist(n_phi,n_th,nR)=this%vtc(n_phi,n_th)
               vp_old_dist(n_phi,n_th,nR)=this%vpc(n_phi,n_th)
            end do
         end do
      end if

      !$omp end parallel

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod

!----------------------------------------------------------------------------
module grid_space_arrays_3d_mod

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
   use horizontal_data, only: osn2, cosn2, sinTheta, cosTheta, osn1, phi
   use parallel_mod, only: get_openmp_blocks
   use constants, only: two, third
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF, &
       &            l_RMS, l_chemical_conv, l_precession, l_mag,      &
       &            l_centrifuge, l_adv_curl, l_viscBcCalc, l_heat
   use time_schemes, only: type_tscheme
   use communications, only: slice_f, gather_f

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_3D_arrays_t
      !----- Nonlinear terms in phi/theta space:
      real(cp), allocatable :: PCr(:,:,:), PCt(:,:,:), PCp(:,:,:)
      real(cp), allocatable :: CAr(:,:,:), CAt(:,:,:)
      real(cp), pointer :: NSadv(:,:,:,:), heatadv(:,:,:,:), LF(:,:,:,:)
      real(cp), pointer :: compadv(:,:,:,:), anel(:,:,:,:), emf(:,:,:,:)
      real(cp), pointer :: Advr(:,:,:), Advt(:,:,:), Advp(:,:,:)
      real(cp), pointer :: LFr(:,:,:), LFt(:,:,:), LFp(:,:,:)
      real(cp), pointer :: VxBr(:,:,:), VxBt(:,:,:), VxBp(:,:,:)
      real(cp), pointer :: VSr(:,:,:), VSt(:,:,:), VSp(:,:,:)
      real(cp), pointer :: VXir(:,:,:), VXit(:,:,:), VXip(:,:,:)
      real(cp), pointer :: ViscHeat(:,:,:), OhmLoss(:,:,:)

      !---- RMS calculations
      real(cp), pointer :: RMS(:,:,:,:), dtV(:,:,:,:), LF2(:,:,:,:)
      real(cp), pointer :: Advt2(:,:,:), Advp2(:,:,:), LFt2(:,:,:), LFp2(:,:,:)
      real(cp), pointer :: CFt2(:,:,:), CFp2(:,:,:)
      real(cp), pointer :: dtVr(:,:,:), dtVp(:,:,:), dtVt(:,:,:)
      real(cp), pointer :: dpkindrc(:,:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), pointer :: vel(:,:,:,:), gradvel(:,:,:,:), xic(:,:,:)
      real(cp), pointer :: mag(:,:,:,:), grads(:,:,:,:), gradp(:,:,:,:)
      real(cp), pointer :: pc(:,:,:), sc(:,:,:), dpdtc(:,:,:), dpdpc(:,:,:)
      real(cp), pointer :: vrc(:,:,:), vtc(:,:,:), vpc(:,:,:)
      real(cp), pointer :: dvrdrc(:,:,:), dvtdrc(:,:,:), dvpdrc(:,:,:)
      real(cp), pointer :: cvrc(:,:,:), drSc(:,:,:)
      real(cp), pointer :: dvrdtc(:,:,:), dvrdpc(:,:,:)
      real(cp), pointer :: dvtdpc(:,:,:), dvpdpc(:,:,:)
      real(cp), pointer :: brc(:,:,:), btc(:,:,:), bpc(:,:,:)
      real(cp), pointer :: cbrc(:,:,:), cbtc(:,:,:), cbpc(:,:,:)
      real(cp), pointer :: cvtc(:,:,:), cvpc(:,:,:)
      real(cp), pointer :: dsdtc(:,:,:), dsdpc(:,:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: output
      procedure :: output_nl_input
      procedure :: get_nl

   end type grid_3D_arrays_t

   real(cp), allocatable :: vr_old(:,:,:), vt_old(:,:,:), vp_old(:,:,:)
   real(cp), allocatable :: vr_old_dist(:,:,:), vt_old_dist(:,:,:), vp_old_dist(:,:,:)

contains

!----------------------------------------------------------------------------
   subroutine initialize(this)

      class(grid_3D_arrays_t) :: this

      allocate( this%NSadv(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
      this%Advr(1:,nThetaStart:,nRstart:) => &
      &         this%NSadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%Advt(1:,nThetaStart:,nRstart:) => &
      &         this%NSadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
      this%Advp(1:,nThetaStart:,nRstart:) => &
      &         this%NSadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
      bytes_allocated=bytes_allocated+3*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL

      allocate( this%anel(nrp,nThetaStart:nThetaStop,nRstart:nRstop,2) )
      this%ViscHeat(1:,nThetaStart:,nRstart:) => &
      &         this%anel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%OhmLoss(1:,nThetaStart:,nRstart:) => &
      &         this%anel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
      bytes_allocated=bytes_allocated+2*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL

      if ( l_mag .or. l_mag_LF ) then
         allocate( this%LF(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%LFr(1:,nThetaStart:,nRstart:) => &
         &         this%LF(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%LFt(1:,nThetaStart:,nRstart:) => &
         &         this%LF(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%LFp(1:,nThetaStart:,nRstart:) => &
         &         this%LF(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)

         allocate( this%emf(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%VxBr(1:,nThetaStart:,nRstart:) => &
         &         this%emf(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%VxBt(1:,nThetaStart:,nRstart:) => &
         &         this%emf(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%VxBp(1:,nThetaStart:,nRstart:) => &
         &         this%emf(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
         bytes_allocated=bytes_allocated+6*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      if ( l_heat ) then
         allocate( this%heatadv(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%VSr(1:,nThetaStart:,nRstart:) => &
         &         this%heatadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%VSt(1:,nThetaStart:,nRstart:) => &
         &         this%heatadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%VSp(1:,nThetaStart:,nRstart:) => &
         &         this%heatadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
         bytes_allocated=bytes_allocated+3*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      if ( l_precession ) then
         allocate( this%PCr(nrp,nThetaStart:nThetaStop,nRstart:nRstop), &
         &         this%PCt(nrp,nThetaStart:nThetaStop,nRstart:nRstop), &
         &         this%PCp(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      if ( l_centrifuge ) then
         allocate( this%CAr(nrp,nThetaStart:nThetaStop,nRstart:nRstop), &
         &         this%CAt(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 2*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate( this%compadv(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         this%VXir(1:,nThetaStart:,nRstart:) => &
         &         this%compadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%VXit(1:,nThetaStart:,nRstart:) => &
         &         this%compadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%VXip(1:,nThetaStart:,nRstart:) => &
         &         this%compadv(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
         bytes_allocated=bytes_allocated+3*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      if ( l_adv_curl ) then
         allocate( this%vel(nrp,nThetaStart:nThetaStop,nRstart:nRstop,6) )
         bytes_allocated=bytes_allocated+6*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      else
         allocate( this%vel(nrp,nThetaStart:nThetaStop,nRstart:nRstop,4) )
         bytes_allocated=bytes_allocated+4*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if
      this%vrc(1:,nThetaStart:,nRstart:) => &
      &          this%vel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%vtc(1:,nThetaStart:,nRstart:) => &
      &          this%vel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
      this%vpc(1:,nThetaStart:,nRstart:) => &
      &          this%vel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
      this%cvrc(1:,nThetaStart:,nRstart:) => &
      &          this%vel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,4)
      if ( l_adv_curl ) then
         this%cvtc(1:,nThetaStart:,nRstart:) => &
         &          this%vel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,5)
         this%cvpc(1:,nThetaStart:,nRstart:) => &
         &          this%vel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,6)
      end if

      allocate( this%gradvel(nrp,nThetaStart:nThetaStop,nRstart:nRstop,7) )
      bytes_allocated=bytes_allocated+7*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      this%dvrdrc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
      this%dvtdrc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
      this%dvpdrc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
      this%dvrdpc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,4)
      this%dvtdpc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,5)
      this%dvpdpc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,6)
      this%dvrdtc(1:,nThetaStart:,nRstart:) => &
      &          this%gradvel(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,7)

      if ( l_mag .or. l_mag_LF ) then
         allocate( this%mag(nrp,nThetaStart:nThetaStop,nRstart:nRstop,6) )
         bytes_allocated=bytes_allocated+6*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
         this%brc(1:,nThetaStart:,nRstart:) => &
         &          this%mag(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%btc(1:,nThetaStart:,nRstart:) => &
         &          this%mag(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%bpc(1:,nThetaStart:,nRstart:) => &
         &          this%mag(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
         this%cbrc(1:,nThetaStart:,nRstart:) => &
         &          this%mag(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,4)
         this%cbtc(1:,nThetaStart:,nRstart:) => &
         &          this%mag(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,5)
         this%cbpc(1:,nThetaStart:,nRstart:) => &
         &          this%mag(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,6)
      end if

      allocate( this%sc(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( this%pc(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
      bytes_allocated=bytes_allocated+2*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         bytes_allocated=bytes_allocated+nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1,nRstart:nRstop) )
      end if

      if ( l_viscBcCalc ) then
         allocate( this%grads(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated=bytes_allocated+3*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      else 
         allocate( this%grads(nrp,nThetaStart:nThetaStop,nRstart:nRstop,1) )
         bytes_allocated=bytes_allocated+nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
      end if
      this%drSc(1:,nThetaStart:,nRstart:) => &
      &          this%grads(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
      if ( l_viscBcCalc ) then
         this%dsdtc(1:,nThetaStart:,nRstart:) => &
         &          this%grads(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%dsdpc(1:,nThetaStart:,nRstart:) => &
         &          this%grads(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         allocate( this%gradp(nrp,nThetaStart:nThetaStop,nRstart:nRstop,2) )
         bytes_allocated=bytes_allocated+2*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
         this%dpdtc(1:,nThetaStart:,nRstart:) => &
         &          this%gradp(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%dpdpc(1:,nThetaStart:,nRstart:) => &
         &          this%gradp(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)

         if ( l_adv_curl ) then
            allocate( this%RMS(nrp,nThetaStart:nThetaStop,nRstart:nRstop,5) )
            bytes_allocated=bytes_allocated+5*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
         else
            allocate( this%RMS(nrp,nThetaStart:nThetaStop,nRstart:nRstop,4) )
            bytes_allocated=bytes_allocated+4*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
         end if
         this%CFt2(1:,nThetaStart:,nRstart:) => &
         &         this%RMS(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%CFp2(1:,nThetaStart:,nRstart:) => &
         &         this%RMS(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%Advt2(1:,nThetaStart:,nRstart:) => &
         &         this%RMS(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)
         this%Advp2(1:,nThetaStart:,nRstart:) => &
         &         this%RMS(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,4)
         if ( l_adv_curl ) then
            this%dpkindrc(1:,nThetaStart:,nRstart:) => &
            &         this%RMS(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,5)
         end if

         allocate( this%dtV(nrp,nThetaStart:nThetaStop,nRstart:nRstop,3) )
         bytes_allocated=bytes_allocated+3*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
         this%dtVr(1:,nThetaStart:,nRstart:) => &
         &         this%dtV(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
         this%dtVt(1:,nThetaStart:,nRstart:) => &
         &         this%dtV(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         this%dtVp(1:,nThetaStart:,nRstart:) => &
         &         this%dtV(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,3)

         if ( l_mag ) then
            allocate( this%LF2(nrp,nThetaStart:nThetaStop,nRstart:nRstop,2) )
            bytes_allocated=bytes_allocated+2*nrp*n_theta_loc*n_r_loc*SIZEOF_DEF_REAL
            this%LFt2(1:,nThetaStart:,nRstart:) => &
            &         this%LF2(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,1)
            this%LFp2(1:,nThetaStart:,nRstart:) => &
            &         this%LF2(1:nrp,nThetaStart:nThetaStop,nRstart:nRstop,2)
         end if

         allocate( vr_old_dist(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         allocate( vp_old_dist(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         allocate( vt_old_dist(nrp,nThetaStart:nThetaStop,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_loc*n_r_loc**&
         &               SIZEOF_DEF_REAL

         this%dtVr(:,:,:)=0.0_cp
         this%dtVt(:,:,:)=0.0_cp
         this%dtVp(:,:,:)=0.0_cp
         vt_old_dist(:,:,:)=0.0_cp
         vr_old_dist(:,:,:)=0.0_cp
         vp_old_dist(:,:,:)=0.0_cp

      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_3D_arrays_t) :: this

      deallocate( this%NSadv, this%anel )
      if ( l_heat ) deallocate( this%heatadv )
      if ( l_mag .or. l_mag_LF ) deallocate( this%LF, this%emf )
      if ( l_chemical_conv ) deallocate( this%compadv )
      if ( l_precession ) deallocate( this%PCr, this%PCt, this%PCp )
      if ( l_centrifuge ) deallocate( this%CAr, this%CAt )

      !----- Fields calculated from these help arrays by legtf:
      deallocate( this%vel, this%gradvel, this%sc, this%grads, this%pc )
      if ( l_chemical_conv ) deallocate( this%xic )
      if ( l_mag .or. l_mag_LF ) deallocate ( this%mag)

      !-- RMS Calculations
      if ( l_RMS ) then
         deallocate ( this%dtV, this%RMS, this%LF2, this%gradp )
         if (allocated(vr_old)) deallocate ( vr_old )
         if (allocated(vt_old)) deallocate ( vt_old )
         if (allocated(vp_old)) deallocate ( vp_old )
         if (allocated(vr_old_dist)) deallocate ( vr_old_dist )
         if (allocated(vt_old_dist)) deallocate ( vt_old_dist )
         if (allocated(vp_old_dist)) deallocate ( vp_old_dist )
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine output(this)

      class(grid_3D_arrays_t) :: this

      write(*,"(A,3ES20.12)") "Advr,Advt,Advp = ",sum(this%Advr), &
                                   sum(this%Advt),sum(this%Advp)

   end subroutine output
!----------------------------------------------------------------------------
   subroutine output_nl_input(this)

      class(grid_3D_arrays_t) :: this

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

      class(grid_3D_arrays_t) :: this

      !-- Input of variables:
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: nR
      logical,             intent(in) :: lRmsCalc
      integer,             intent(in) :: nBc

      !-- Local variables:
      integer :: n_th, nThetaNHS, n_phi, nThStart, nThStop
      real(cp) :: or4sn2, csn2, cnt, rsnt, snt, posnalp, O_dt

      !$omp parallel default(shared) private(nThStart,nThStop) &
      !$omp private(n_th,nThetaNHS,n_phi) &
      !$omp private(or4sn2,csn2,cnt,snt,rsnt,posnalp)
      nThStart=nThetaStart; nThStop=nThetaStop
      call get_openmp_blocks(nThStart,nThStop)

      if ( l_mag_LF .and. (nBc == 0 .or. lRmsCalc) .and. nR>n_r_LCR ) then
         !------ Get the Lorentz force:
         do n_th=nThStart,nThStop

            nThetaNHS=(n_th+1)/2
            or4sn2   =or4(nR)*osn2(nThetaNHS)

            do n_phi=1,n_phi_max
               !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
               this%LFr(n_phi,n_th,nR)=  LFfac*osn2(nThetaNHS) * (         &
               &        this%cbtc(n_phi,n_th,nR)*this%bpc(n_phi,n_th,nR) - &
               &        this%cbpc(n_phi,n_th,nR)*this%btc(n_phi,n_th,nR) )
            end do

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            do n_phi=1,n_phi_max
               this%LFt(n_phi,n_th,nR)=           LFfac*or4sn2 * (         &
               &        this%cbpc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR) - &
               &        this%cbrc(n_phi,n_th,nR)*this%bpc(n_phi,n_th,nR) )
            end do

            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            do n_phi=1,n_phi_max
               this%LFp(n_phi,n_th,nR)=           LFfac*or4sn2 * (         &
               &        this%cbrc(n_phi,n_th,nR)*this%btc(n_phi,n_th,nR) - &
               &        this%cbtc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR) )
            end do

         end do   ! theta loop
      end if      ! Lorentz force required ?

      if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

         if ( l_adv_curl ) then ! Advection is \curl{u} \times u

            do n_th=nThStart,nThStop ! loop over theta points in block

               nThetaNHS=(n_th+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)

               do n_phi=1,n_phi_max
                  this%Advr(n_phi,n_th,nR)=       - osn2(nThetaNHS) * (       &
                  &        this%cvtc(n_phi,n_th,nR)*this%vpc(n_phi,n_th,nR) - &
                  &        this%cvpc(n_phi,n_th,nR)*this%vtc(n_phi,n_th,nR) )

                  this%Advt(n_phi,n_th,nR)=        -        or4sn2 * (        &
                  &        this%cvpc(n_phi,n_th,nR)*this%vrc(n_phi,n_th,nR) - &
                  &        this%cvrc(n_phi,n_th,nR)*this%vpc(n_phi,n_th,nR) )

                  this%Advp(n_phi,n_th,nR)=        -        or4sn2 * (        &
                  &        this%cvrc(n_phi,n_th,nR)*this%vtc(n_phi,n_th,nR) - &
                  &        this%cvtc(n_phi,n_th,nR)*this%vrc(n_phi,n_th,nR) )
               end do

            end do

         else ! Advection is u\grad u

            !------ Get Advection:
            do n_th=nThStart,nThStop ! loop over theta points in block
               nThetaNHS=(n_th+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               csn2     =cosn2(nThetaNHS)
               if ( mod(n_th,2) == 0 ) csn2=-csn2 ! South, odd function in theta

               do n_phi=1,n_phi_max
                  this%Advr(n_phi,n_th,nR)=          -or2(nR)*orho1(nR) * (  &
                  &                                this%vrc(n_phi,n_th,nR) * &
                  &                     (       this%dvrdrc(n_phi,n_th,nR) - &
                  &    ( two*or1(nR)+beta(nR) )*this%vrc(n_phi,n_th,nR) ) +  &
                  &                               osn2(nThetaNHS) * (        &
                  &                                this%vtc(n_phi,n_th,nR) * &
                  &                     (       this%dvrdtc(n_phi,n_th,nR) - &
                  &                  r(nR)*      this%vtc(n_phi,n_th,nR) ) + &
                  &                                this%vpc(n_phi,n_th,nR) * &
                  &                     (       this%dvrdpc(n_phi,n_th,nR) - &
                  &                    r(nR)*      this%vpc(n_phi,n_th,nR) ) ) )
               end do

               do n_phi=1,n_phi_max
                  this%Advt(n_phi,n_th,nR)=         or4sn2*orho1(nR) * (  &
                  &                            -this%vrc(n_phi,n_th,nR) * &
                  &                      (   this%dvtdrc(n_phi,n_th,nR) - &
                  &                beta(nR)*this%vtc(n_phi,n_th,nR) )   + &
                  &                             this%vtc(n_phi,n_th,nR) * &
                  &                      ( csn2*this%vtc(n_phi,n_th,nR) + &
                  &                          this%dvpdpc(n_phi,n_th,nR) + &
                  &                      this%dvrdrc(n_phi,n_th,nR) )   + &
                  &                             this%vpc(n_phi,n_th,nR) * &
                  &                      ( csn2*this%vpc(n_phi,n_th,nR) - &
                  &                          this%dvtdpc(n_phi,n_th,nR) )  )
               end do

               do n_phi=1,n_phi_max
                  this%Advp(n_phi,n_th,nR)=         or4sn2*orho1(nR) * (  &
                  &                            -this%vrc(n_phi,n_th,nR) * &
                  &                        ( this%dvpdrc(n_phi,n_th,nR) - &
                  &                beta(nR)*this%vpc(n_phi,n_th,nR) )   - &
                  &                             this%vtc(n_phi,n_th,nR) * &
                  &                        ( this%dvtdpc(n_phi,n_th,nR) + &
                  &                        this%cvrc(n_phi,n_th,nR) )   - &
                  &       this%vpc(n_phi,n_th,nR) * this%dvpdpc(n_phi,n_th,nR) )
               end do

            end do ! theta loop
         end if

      end if  ! Navier-Stokes nonlinear advection term ?

      if ( l_heat_nl .and. nBc == 0 ) then
         !------ Get V S, the divergence of it is entropy advection:
         do n_th=nThStart,nThStop
            do n_phi=1,n_phi_max     ! calculate v*s components
               this%VSr(n_phi,n_th,nR)= &
               &    this%vrc(n_phi,n_th,nR)*this%sc(n_phi,n_th,nR)
               this%VSt(n_phi,n_th,nR)= &
               &    or2(nR)*this%vtc(n_phi,n_th,nR)*this%sc(n_phi,n_th,nR)
               this%VSp(n_phi,n_th,nR)= &
               &    or2(nR)*this%vpc(n_phi,n_th,nR)*this%sc(n_phi,n_th,nR)
            end do
         end do  ! theta loop
      end if     ! heat equation required ?

      if ( l_chemical_conv .and. nBc == 0 ) then
         !------ Get V S, the divergence of it is the advection of chem comp:
         do n_th=nThStart,nThStop
            do n_phi=1,n_phi_max     ! calculate v*s components
               this%VXir(n_phi,n_th,nR)= &
               &    this%vrc(n_phi,n_th,nR)*this%xic(n_phi,n_th,nR)
               this%VXit(n_phi,n_th,nR)= &
               &    or2(nR)*this%vtc(n_phi,n_th,nR)*this%xic(n_phi,n_th,nR)
               this%VXip(n_phi,n_th,nR)= &
               &    or2(nR)*this%vpc(n_phi,n_th,nR)*this%xic(n_phi,n_th,nR)
            end do
         end do  ! theta loop
      end if     ! chemical composition equation required ?

      if ( l_precession .and. nBc == 0 ) then

         do n_th=nThStart,nThStop
            nThetaNHS=(n_th+1)/2
            posnalp=-two*oek*po*sin(prec_angle)*osn1(nThetaNHS)
            cnt=cosTheta(n_th)
            do n_phi=1,n_phi_max
               this%PCr(n_phi,n_th,nR)=posnalp*r(nR)*(cos(oek*time+phi(n_phi))* &
               &                                  this%vpc(n_phi,n_th,nR)*cnt   &
               &            +sin(oek*time+phi(n_phi))*this%vtc(n_phi,n_th,nR))
               this%PCt(n_phi,n_th,nR)=   -posnalp*or2(nR)*(                    &
               &               cos(oek*time+phi(n_phi))*this%vpc(n_phi,n_th,nR) &
               &      +sin(oek*time+phi(n_phi))*or1(nR)*this%vrc(n_phi,n_th,nR) )
               this%PCp(n_phi,n_th,nR)= posnalp*cos(oek*time+phi(n_phi))*       &
               &              or2(nR)*(      this%vtc(n_phi,n_th,nR)-           &
               &                     or1(nR)*this%vrc(n_phi,n_th,nR)*cnt)
            end do
         end do ! theta loop
      end if ! precession term required ?

      if ( l_centrifuge .and. nBc ==0 ) then
         do n_th=nThStart,nThStop
            nThetaNHS=(n_th+1)/2
            snt=sinTheta(n_th)
            cnt=cosTheta(n_th)
            rsnt=r(nR)*snt
            do n_phi=1,n_phi_max
               if ( l_anel ) then
                  this%CAr(n_phi,n_th,nR) = dilution_fac*rsnt*snt* &
                  &       ( -ra*opr*this%sc(n_phi,n_th,nR) )
                  !-- neglect pressure contribution
                  !& + polind*DissNb*oek*opressure0(nR)*this%pc(n_phi,n_th,nR) )
                  this%CAt(n_phi,n_th,nR) = dilution_fac*rsnt*cnt* &
                  &       ( -ra*opr*this%sc(n_phi,n_th,nR) )
                  !-- neglect pressure contribution
                  !& + polind*DissNb*oek*opressure0(nR)*this%pc(n_phi,n_th,nR) )
               else
                  this%CAr(n_phi,n_th,nR) = -dilution_fac*rsnt*snt*ra*opr* &
                  &                       this%sc(n_phi,n_th,nR)
                  this%CAt(n_phi,n_th,nR) = -dilution_fac*rsnt*cnt*ra*opr* &
                  &                       this%sc(n_phi,n_th,nR)
               end if
            end do ! phi loop
         end do ! theta loop
      end if ! centrifuge

      if ( l_mag_nl ) then

         if ( nBc == 0 .and. nR>n_r_LCR ) then

            !------ Get (V x B) , the curl of this is the dynamo term:
            do n_th=nThStart,nThStop
               nThetaNHS=(n_th+1)/2

               do n_phi=1,n_phi_max
                  this%VxBr(n_phi,n_th,nR)=  orho1(nR)*osn2(nThetaNHS) * (         &
                  &              this%vtc(n_phi,n_th,nR)*this%bpc(n_phi,n_th,nR) - &
                  &              this%vpc(n_phi,n_th,nR)*this%btc(n_phi,n_th,nR) )
               end do

               do n_phi=1,n_phi_max
                  this%VxBt(n_phi,n_th,nR)=  orho1(nR)*or4(nR) * (          &
                  &       this%vpc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR) - &
                  &       this%vrc(n_phi,n_th,nR)*this%bpc(n_phi,n_th,nR) )
               end do

               do n_phi=1,n_phi_max
                  this%VxBp(n_phi,n_th,nR)=   orho1(nR)*or4(nR) * (          &
                  &        this%vrc(n_phi,n_th,nR)*this%btc(n_phi,n_th,nR) - &
                  &        this%vtc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR) )
               end do
            end do   ! theta loop

         else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary

            do n_th=nThStart,nThStop
               do n_phi=1,n_phi_max
                  this%VxBt(n_phi,n_th,nR)=  or4(nR) * orho1(nR) * &
                  &        this%vpc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR)
                  this%VxBp(n_phi,n_th,nR)= -or4(nR) * orho1(nR) * &
                  &        this%vtc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR)
               end do
            end do

         else if ( nBc == 2 ) then  ! rigid boundary :

            !----- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):

            do n_th=nThStart,nThStop
               do n_phi=1,n_phi_max
                  this%VxBt(n_phi,n_th,nR)= or4(nR) * orho1(nR) * &
                  &    this%vpc(n_phi,n_th,nR)*this%brc(n_phi,n_th,nR)

                  this%VxBp(n_phi,n_th,nR)= 0.0_cp
               end do
            end do

         end if  ! boundary ?

      end if ! l_mag_nl ?

      if ( l_anel .and. nBc == 0 ) then
         !------ Get viscous heating
         do n_th=nThStart,nThStop ! loop over theta points in block
            nThetaNHS=(n_th+1)/2
            csn2     =cosn2(nThetaNHS)
            if ( mod(n_th,2) == 0 ) csn2=-csn2 ! South, odd function in theta

            do n_phi=1,n_phi_max
               this%ViscHeat(n_phi,n_th,nR)=      or4(nR)*                  &
               &                     orho1(nR)*otemp1(nR)*visc(nR)*(        &
               &     two*(                     this%dvrdrc(n_phi,n_th,nR) - & ! (1)
               &     (two*or1(nR)+beta(nR))*this%vrc(n_phi,n_th,nR) )**2  + &
               &     two*( csn2*                  this%vtc(n_phi,n_th,nR) + &
               &                               this%dvpdpc(n_phi,n_th,nR) + &
               &                               this%dvrdrc(n_phi,n_th,nR) - & ! (2)
               &     or1(nR)*               this%vrc(n_phi,n_th,nR) )**2  + &
               &     two*(                     this%dvpdpc(n_phi,n_th,nR) + &
               &           csn2*                  this%vtc(n_phi,n_th,nR) + & ! (3)
               &     or1(nR)*               this%vrc(n_phi,n_th,nR) )**2  + &
               &          ( two*               this%dvtdpc(n_phi,n_th,nR) + &
               &                                 this%cvrc(n_phi,n_th,nR) - & ! (6)
               &      two*csn2*             this%vpc(n_phi,n_th,nR) )**2  + &
               &                                 osn2(nThetaNHS) * (        &
               &         ( r(nR)*              this%dvtdrc(n_phi,n_th,nR) - &
               &           (two+beta(nR)*r(nR))*  this%vtc(n_phi,n_th,nR) + & ! (4)
               &     or1(nR)*            this%dvrdtc(n_phi,n_th,nR) )**2  + &
               &         ( r(nR)*              this%dvpdrc(n_phi,n_th,nR) - &
               &           (two+beta(nR)*r(nR))*  this%vpc(n_phi,n_th,nR) + & ! (5)
               &     or1(nR)*            this%dvrdpc(n_phi,n_th,nR) )**2 )- &
               &    two*third*(  beta(nR)*        this%vrc(n_phi,n_th,nR) )**2 )
            end do
         end do ! theta loop

         if ( l_mag_nl .and. nR>n_r_LCR ) then
            !------ Get ohmic losses
            do n_th=nThStart,nThStop ! loop over theta points in block
               nThetaNHS=(n_th+1)/2
               do n_phi=1,n_phi_max
                  this%OhmLoss(n_phi,n_th,nR)= or2(nR)*otemp1(nR)*lambda(nR)*  &
                  &    ( or2(nR)*                this%cbrc(n_phi,n_th,nR)**2 + &
                  &      osn2(nThetaNHS)*        this%cbtc(n_phi,n_th,nR)**2 + &
                  &      osn2(nThetaNHS)*        this%cbpc(n_phi,n_th,nR)**2  )
               end do
            end do ! theta loop

         end if ! if l_mag_nl ?

      end if  ! Viscous heating and Ohmic losses ?

      if ( lRmsCalc ) then
         do n_th=nThStart,nThStop ! loop over theta points in block
            snt=sinTheta(n_th)
            cnt=cosTheta(n_th)
            rsnt=r(nR)*snt
            do n_phi=1,n_phi_max
               this%dpdtc(n_phi,n_th,nR)=this%dpdtc(n_phi,n_th,nR)*or1(nR)
               this%dpdpc(n_phi,n_th,nR)=this%dpdpc(n_phi,n_th,nR)*or1(nR)
               this%CFt2(n_phi,n_th,nR)=-two*CorFac*cnt*this%vpc(n_phi,n_th,nR)* &
               &                         or1(nR)
               this%CFp2(n_phi,n_th,nR)= two*CorFac*snt* (                &
               &                     cnt*this%vtc(n_phi,n_th,nR)/rsnt +   &
               &                     or2(nR)*snt*this%vrc(n_phi,n_th,nR) )
               if ( l_conv_nl ) then
                  this%Advt2(n_phi,n_th,nR)=rsnt*snt*this%Advt(n_phi,n_th,nR)
                  this%Advp2(n_phi,n_th,nR)=rsnt*snt*this%Advp(n_phi,n_th,nR)
               end if
               if ( l_mag_LF .and. nR > n_r_LCR ) then
                  this%LFt2(n_phi,n_th,nR)=rsnt*snt*this%LFt(n_phi,n_th,nR)
                  this%LFp2(n_phi,n_th,nR)=rsnt*snt*this%LFp(n_phi,n_th,nR)
               end if
            end do
         end do

         if ( l_adv_curl ) then
            do n_th=nThStart,nThStop ! loop over theta points in block
               nThetaNHS=(n_th+1)/2
               csn2     =cosn2(nThetaNHS)
               if ( mod(n_th,2) == 0 ) csn2=-csn2 ! South, odd function in theta
               do n_phi=1,n_phi_max
                  this%dpdtc(n_phi,n_th,nR)=this%dpdtc(n_phi,n_th,nR)-              &
                  &                      or3(nR)*( or2(nR)*                         &
                  &            this%vrc(n_phi,n_th,nR)*this%dvrdtc(n_phi,n_th,nR) - &
                  &            this%vtc(n_phi,n_th,nR)*(this%dvrdrc(n_phi,n_th,nR)+ &
                  &            this%dvpdpc(n_phi,n_th,nR)+csn2 *                    &
                  &            this%vtc(n_phi,n_th,nR))+ this%vpc(n_phi,n_th,nR)*(  &
                  &            this%cvrc(n_phi,n_th,nR)+this%dvtdpc(n_phi,n_th,nR)- &
                  &            csn2*this%vpc(n_phi,n_th,nR)) )
                  this%dpdpc(n_phi,n_th,nR)=this%dpdpc(n_phi,n_th,nR)-              &
                  &                         or3(nR)*( or2(nR)*                      &
                  &            this%vrc(n_phi,n_th,nR)*this%dvrdpc(n_phi,n_th,nR) + &
                  &            this%vtc(n_phi,n_th,nR)*this%dvtdpc(n_phi,n_th,nR) + &
                  &            this%vpc(n_phi,n_th,nR)*this%dvpdpc(n_phi,n_th,nR) )
                  if ( l_conv_nl ) then
                     this%Advt2(n_phi,n_th,nR)=this%Advt2(n_phi,n_th,nR)-              &
                     &                   or3(nR)*( or2(nR)*                         &
                     &         this%vrc(n_phi,n_th,nR)*this%dvrdtc(n_phi,n_th,nR) - &
                     &         this%vtc(n_phi,n_th,nR)*(this%dvrdrc(n_phi,n_th,nR)+ &
                     &         this%dvpdpc(n_phi,n_th,nR)+csn2 *                    &
                     &         this%vtc(n_phi,n_th,nR))+ this%vpc(n_phi,n_th,nR)*(  &
                     &         this%cvrc(n_phi,n_th,nR)+this%dvtdpc(n_phi,n_th,nR)- &
                     &         csn2*this%vpc(n_phi,n_th,nR)) )
                     this%Advp2(n_phi,n_th,nR)=this%Advp2(n_phi,n_th,nR)-              &
                     &                   or3(nR)*( or2(nR)*                         &
                     &         this%vrc(n_phi,n_th,nR)*this%dvrdpc(n_phi,n_th,nR) + &
                     &         this%vtc(n_phi,n_th,nR)*this%dvtdpc(n_phi,n_th,nR) + &
                     &         this%vpc(n_phi,n_th,nR)*this%dvpdpc(n_phi,n_th,nR) )
                  end if

                  !- dpkin/dr = 1/2 d (u^2) / dr = ur*dur/dr+ut*dut/dr+up*dup/dr
                  this%dpkindrc(n_phi,n_th,nR)=or4(nR)*this%vrc(n_phi,n_th,nR)*(  &
                  &                         this%dvrdrc(n_phi,n_th,nR)-           &
                  &                         two*or1(nR)*this%vrc(n_phi,n_th,nR)) +&
                  &                         or2(nR)*osn2(nThetaNHS)*(             &
                  &                                 this%vtc(n_phi,n_th,nR)*(     &
                  &                         this%dvtdrc(n_phi,n_th,nR)-           &
                  &                         or1(nR)*this%vtc(n_phi,n_th,nR) ) +   &
                  &                                 this%vpc(n_phi,n_th,nR)*(     &
                  &                         this%dvpdrc(n_phi,n_th,nR)-           &
                  &                         or1(nR)*this%vpc(n_phi,n_th,nR) ) )
               end do
            end do
         end if
      end if

      if ( l_RMS .and. tscheme%istage == 1 ) then
         O_dt = 1.0_cp/tscheme%dt(1)
         do n_th=nThStart,nThStop
            do n_phi=1,n_phi_max
               this%dtVr(n_phi,n_th,nR)=O_dt*or2(nR)*(this%vrc(n_phi,n_th,nR)- &
               &                             vr_old_dist(n_phi,n_th,nR))
               this%dtVt(n_phi,n_th,nR)=O_dt*or1(nR)*(this%vtc(n_phi,n_th,nR)- &
               &                             vt_old_dist(n_phi,n_th,nR))
               this%dtVp(n_phi,n_th,nR)=O_dt*or1(nR)*(this%vpc(n_phi,n_th,nR)- &
               &                             vp_old_dist(n_phi,n_th,nR))

               vr_old_dist(n_phi,n_th,nR)=this%vrc(n_phi,n_th,nR)
               vt_old_dist(n_phi,n_th,nR)=this%vtc(n_phi,n_th,nR)
               vp_old_dist(n_phi,n_th,nR)=this%vpc(n_phi,n_th,nR)
            end do
         end do
      end if

      !$omp end parallel

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_3d_mod
