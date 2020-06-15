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
   use truncation, only: nrp, n_phi_max, n_theta_max
   use radial_data, only: nRstart, nRstop
   use radial_functions, only: or2, orho1, beta, otemp1, visc, r, or3, &
       &                       lambda, or4, or1, alpha0, temp0, opressure0
   use physical_parameters, only: LFfac, n_r_LCR, CorFac, prec_angle,    &
        &                         ThExpNb, ViscHeatFac, oek, po, DissNb, &
        &                         dilution_fac, ra, opr, polind, strat, radratio
   use blocking, only: nfs, sizeThetaB
   use horizontal_data, only: osn2, cosn2, sinTheta, cosTheta, osn1, phi
   use parallel_mod, only: get_openmp_blocks
   use constants, only: two, third
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF, &
       &            l_RMS, l_chemical_conv, l_precession,             &
       &            l_centrifuge, l_adv_curl
   use time_schemes, only: type_tscheme

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
#ifdef WITH_SHTNS
      procedure :: get_nl_shtns
#endif

   end type grid_space_arrays_t

   real(cp), allocatable :: vr_old(:,:,:), vt_old(:,:,:), vp_old(:,:,:)

contains

   subroutine initialize(this)

      class(grid_space_arrays_t) :: this

      allocate( this%Advr(nrp,nfs), this%Advt(nrp,nfs), this%Advp(nrp,nfs) )
      allocate( this%LFr(nrp,nfs), this%LFt(nrp,nfs), this%LFp(nrp,nfs) )
      allocate( this%VxBr(nrp,nfs), this%VxBt(nrp,nfs), this%VxBp(nrp,nfs) )
      allocate( this%VSr(nrp,nfs), this%VSt(nrp,nfs), this%VSp(nrp,nfs) )
      allocate( this%ViscHeat(nrp,nfs), this%OhmLoss(nrp,nfs) )
      bytes_allocated=bytes_allocated + 14*nrp*nfs*SIZEOF_DEF_REAL

      if ( l_precession ) then
         allocate( this%PCr(nrp,nfs), this%PCt(nrp,nfs), this%PCp(nrp,nfs) )
         bytes_allocated=bytes_allocated + 3*nrp*nfs*SIZEOF_DEF_REAL
      end if

      if ( l_centrifuge ) then
         allocate( this%CAr(nrp,nfs), this%CAt(nrp,nfs) )
         bytes_allocated=bytes_allocated + 2*nrp*nfs*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXir(nrp,nfs), this%VXit(nrp,nfs), this%VXip(nrp,nfs) )
         bytes_allocated=bytes_allocated + 3*nrp*nfs*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(nrp,nfs),this%vtc(nrp,nfs),this%vpc(nrp,nfs) )
      allocate( this%dvrdrc(nrp,nfs),this%dvtdrc(nrp,nfs) )
      allocate( this%dvpdrc(nrp,nfs),this%cvrc(nrp,nfs) )
      allocate( this%dvrdtc(nrp,nfs),this%dvrdpc(nrp,nfs) )
      allocate( this%dvtdpc(nrp,nfs),this%dvpdpc(nrp,nfs) )
      allocate( this%brc(nrp,nfs),this%btc(nrp,nfs),this%bpc(nrp,nfs) )
      this%btc=1.0e50_cp
      this%bpc=1.0e50_cp
      allocate( this%cbrc(nrp,nfs),this%cbtc(nrp,nfs),this%cbpc(nrp,nfs) )
      allocate( this%sc(nrp,nfs),this%drSc(nrp,nfs) )
      allocate( this%pc(nrp,nfs) )
      allocate( this%dsdtc(nrp,nfs),this%dsdpc(nrp,nfs) )
      bytes_allocated=bytes_allocated + 22*nrp*nfs*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(nrp,nfs) )
         bytes_allocated=bytes_allocated + nrp*nfs*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1) )
      end if

      if ( l_adv_curl ) then
         allocate( this%cvtc(nrp,nfs), this%cvpc(nrp,nfs) )
         bytes_allocated=bytes_allocated+2*nrp*nfs*SIZEOF_DEF_REAL
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         allocate ( this%Advt2(nrp,nfs), this%Advp2(nrp,nfs) )
         allocate ( this%dtVr(nrp,nfs), this%dtVt(nrp,nfs), this%dtVp(nrp,nfs) )
         allocate ( this%LFt2(nrp,nfs), this%LFp2(nrp,nfs) )
         allocate ( this%CFt2(nrp,nfs), this%CFp2(nrp,nfs) )
         allocate ( this%dpdtc(nrp,nfs), this%dpdpc(nrp,nfs) )
         bytes_allocated=bytes_allocated + 11*nrp*nfs*SIZEOF_DEF_REAL

         allocate( vt_old(nrp,n_theta_max,nRstart:nRstop) )
         allocate( vp_old(nrp,n_theta_max,nRstart:nRstop) )
         allocate( vr_old(nrp,n_theta_max,nRstart:nRstop) )
         bytes_allocated=bytes_allocated + 3*nrp*n_theta_max*(nRstop-nRstart+1)*&
         &               SIZEOF_DEF_REAL

         this%dtVr(:,:)=0.0_cp
         this%dtVt(:,:)=0.0_cp
         this%dtVp(:,:)=0.0_cp
         vt_old(:,:,:) =0.0_cp
         vr_old(:,:,:) =0.0_cp
         vp_old(:,:,:) =0.0_cp

         if ( l_adv_curl ) then
            allocate ( this%dpkindrc(nrp, nfs) )
            bytes_allocated=bytes_allocated + nrp*nfs*SIZEOF_DEF_REAL
         end if
      end if
      !write(*,"(A,I15,A)") "grid_space_arrays: allocated ",bytes_allocated,"B."

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
         deallocate ( vr_old, vt_old, vp_old )
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
#ifdef WITH_SHTNS
   subroutine get_nl_shtns(this, time, tscheme, nR, nBc, lRmsCalc)
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
      nThStart=1; nThStop=n_theta_max
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
               &                             vr_old(n_phi,n_th,nR))
               this%dtVt(n_phi,n_th)=O_dt*or1(nR)*(this%vtc(n_phi,n_th)- &
               &                             vt_old(n_phi,n_th,nR))
               this%dtVp(n_phi,n_th)=O_dt*or1(nR)*(this%vpc(n_phi,n_th)- &
               &                             vp_old(n_phi,n_th,nR))

               vr_old(n_phi,n_th,nR)=this%vrc(n_phi,n_th)
               vt_old(n_phi,n_th,nR)=this%vtc(n_phi,n_th)
               vp_old(n_phi,n_th,nR)=this%vpc(n_phi,n_th)
            end do
         end do
      end if

      !$omp end parallel


   end subroutine get_nl_shtns
#endif
!----------------------------------------------------------------------------
   subroutine get_nl(this,time,tscheme,nR,nBc,nThetaStart,lRmsCalc)
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
      integer,             intent(in) :: nBc
      integer,             intent(in) :: nThetaStart
      logical,             intent(in) :: lRmsCalc

      !-- Local variables:
      integer :: nTheta
      integer :: nThetaLast,nThetaB,nThetaNHS
      integer :: nPhi
      real(cp) :: or2sn2,or4sn2,csn2,snt,cnt,rsnt,posnalp, O_dt

      nThetaLast=nThetaStart-1

      if ( l_mag_LF .and. (nBc == 0 .or. lRmsCalc) .and. nR>n_r_LCR ) then
         !------ Get the Lorentz force:
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB

            nTheta   =nTheta+1
            nThetaNHS=(nTheta+1)/2
            or4sn2   =or4(nR)*osn2(nThetaNHS)

            do nPhi=1,n_phi_max
               !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
               this%LFr(nPhi,nThetaB)=  LFfac*osn2(nThetaNHS) * (        &
               &        this%cbtc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) - &
               &        this%cbpc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) )
            end do
            this%LFr(n_phi_max+1,nThetaB)=0.0_cp
            this%LFr(n_phi_max+2,nThetaB)=0.0_cp

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            do nPhi=1,n_phi_max
               this%LFt(nPhi,nThetaB)=           LFfac*or4sn2 * (        &
               &        this%cbpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) - &
               &        this%cbrc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) )
            end do
            this%LFt(n_phi_max+1,nThetaB)=0.0_cp
            this%LFt(n_phi_max+2,nThetaB)=0.0_cp
            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            do nPhi=1,n_phi_max
               this%LFp(nPhi,nThetaB)=           LFfac*or4sn2 * (        &
               &        this%cbrc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) - &
               &        this%cbtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) )
            end do
            this%LFp(n_phi_max+1,nThetaB)=0.0_cp
            this%LFp(n_phi_max+2,nThetaB)=0.0_cp

         end do   ! theta loop
      end if      ! Lorentz force required ?

      if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

         if ( l_adv_curl ) then ! Advection is u \omega

            !------ Get Advection:
            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB ! loop over theta points in block
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               do nPhi=1,n_phi_max
                  this%Advr(nPhi,nThetaB)=       - osn2(nThetaNHS) * (      &
                  &        this%cvtc(nPhi,nThetaB)*this%vpc(nPhi,nThetaB) - &
                  &        this%cvpc(nPhi,nThetaB)*this%vtc(nPhi,nThetaB) )

                  this%Advt(nPhi,nThetaB)=        -        or4sn2 * (       &
                  &        this%cvpc(nPhi,nThetaB)*this%vrc(nPhi,nThetaB) - &
                  &        this%cvrc(nPhi,nThetaB)*this%vpc(nPhi,nThetaB) )

                  this%Advp(nPhi,nThetaB)=        -        or4sn2 * (       &
                  &        this%cvrc(nPhi,nThetaB)*this%vtc(nPhi,nThetaB) - &
                  &        this%cvtc(nPhi,nThetaB)*this%vrc(nPhi,nThetaB) )
               end do
               this%Advr(n_phi_max+1,nThetaB)=0.0_cp
               this%Advr(n_phi_max+2,nThetaB)=0.0_cp
               this%Advt(n_phi_max+1,nThetaB)=0.0_cp
               this%Advt(n_phi_max+2,nThetaB)=0.0_cp
               this%Advp(n_phi_max+1,nThetaB)=0.0_cp
               this%Advp(n_phi_max+2,nThetaB)=0.0_cp
            end do

         else ! Advection is u \grad u

            !------ Get Advection:
            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB ! loop over theta points in block
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               csn2     =cosn2(nThetaNHS)
               if ( mod(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

               do nPhi=1,n_phi_max
                  this%Advr(nPhi,nThetaB)=          -or2(nR)*orho1(nR) * (  &
                  &                                this%vrc(nPhi,nThetaB) * &
                  &                     (       this%dvrdrc(nPhi,nThetaB) - &
                  &    ( two*or1(nR)+beta(nR) )*this%vrc(nPhi,nThetaB) ) +  &
                  &                               osn2(nThetaNHS) * (       &
                  &                                this%vtc(nPhi,nThetaB) * &
                  &                     (       this%dvrdtc(nPhi,nThetaB) - &
                  &                  r(nR)*      this%vtc(nPhi,nThetaB) ) + &
                  &                                this%vpc(nPhi,nThetaB) * &
                  &                     (       this%dvrdpc(nPhi,nThetaB) - &
                  &                    r(nR)*      this%vpc(nPhi,nThetaB) ) ) )
               end do
               this%Advr(n_phi_max+1,nThetaB)=0.0_cp
               this%Advr(n_phi_max+2,nThetaB)=0.0_cp
               do nPhi=1,n_phi_max
                  this%Advt(nPhi,nThetaB)=         or4sn2*orho1(nR) * (  &
                  &                            -this%vrc(nPhi,nThetaB) * &
                  &                      (   this%dvtdrc(nPhi,nThetaB) - &
                  &                beta(nR)*this%vtc(nPhi,nThetaB) )   + &
                  &                             this%vtc(nPhi,nThetaB) * &
                  &                      ( csn2*this%vtc(nPhi,nThetaB) + &
                  &                          this%dvpdpc(nPhi,nThetaB) + &
                  &                      this%dvrdrc(nPhi,nThetaB) )   + &
                  &                             this%vpc(nPhi,nThetaB) * &
                  &                      ( csn2*this%vpc(nPhi,nThetaB) - &
                  &                          this%dvtdpc(nPhi,nThetaB) )  )
               end do
               this%Advt(n_phi_max+1,nThetaB)=0.0_cp
               this%Advt(n_phi_max+2,nThetaB)=0.0_cp
               do nPhi=1,n_phi_max
                  this%Advp(nPhi,nThetaB)=         or4sn2*orho1(nR) * (  &
                  &                            -this%vrc(nPhi,nThetaB) * &
                  &                        ( this%dvpdrc(nPhi,nThetaB) - &
                  &                beta(nR)*this%vpc(nPhi,nThetaB) )   - &
                  &                             this%vtc(nPhi,nThetaB) * &
                  &                        ( this%dvtdpc(nPhi,nThetaB) + &
                  &                        this%cvrc(nPhi,nThetaB) )   - &
                  &       this%vpc(nPhi,nThetaB) * this%dvpdpc(nPhi,nThetaB) )
               end do
               this%Advp(n_phi_max+1,nThetaB)=0.0_cp
               this%Advp(n_phi_max+2,nThetaB)=0.0_cp
            end do ! theta loop

         end if ! Curl form or non curl form

      end if  ! Navier-Stokes nonlinear advection term ?

      if ( l_heat_nl .and. nBc == 0 ) then
         !------ Get V S, the divergence of the is entropy advection:
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB
            nTheta   =nTheta+1
            nThetaNHS=(nTheta+1)/2
            or2sn2=or2(nR)*osn2(nThetaNHS)
            do nPhi=1,n_phi_max     ! calculate v*s components
               this%VSr(nPhi,nThetaB)= &
               &    this%vrc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
               this%VSt(nPhi,nThetaB)= &
               &    or2sn2*this%vtc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
               this%VSp(nPhi,nThetaB)= &
               &    or2sn2*this%vpc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
            end do
            this%VSr(n_phi_max+1,nThetaB)=0.0_cp
            this%VSr(n_phi_max+2,nThetaB)=0.0_cp
            this%VSt(n_phi_max+1,nThetaB)=0.0_cp
            this%VSt(n_phi_max+2,nThetaB)=0.0_cp
            this%VSp(n_phi_max+1,nThetaB)=0.0_cp
            this%VSp(n_phi_max+2,nThetaB)=0.0_cp
         end do  ! theta loop
      end if     ! heat equation required ?

      if ( l_chemical_conv .and. nBc == 0 ) then
         !------ Get V Xi, the divergence of the is advection of chemical comp:
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB
            nTheta   =nTheta+1
            nThetaNHS=(nTheta+1)/2
            or2sn2=or2(nR)*osn2(nThetaNHS)
            do nPhi=1,n_phi_max     ! calculate v*s components
               this%VXir(nPhi,nThetaB)= &
               &    this%vrc(nPhi,nThetaB)*this%xic(nPhi,nThetaB)
               this%VXit(nPhi,nThetaB)= &
               &    or2sn2*this%vtc(nPhi,nThetaB)*this%xic(nPhi,nThetaB)
               this%VXip(nPhi,nThetaB)= &
               &    or2sn2*this%vpc(nPhi,nThetaB)*this%xic(nPhi,nThetaB)
            end do
            this%VXir(n_phi_max+1,nThetaB)=0.0_cp
            this%VXir(n_phi_max+2,nThetaB)=0.0_cp
            this%VXit(n_phi_max+1,nThetaB)=0.0_cp
            this%VXit(n_phi_max+2,nThetaB)=0.0_cp
            this%VXip(n_phi_max+1,nThetaB)=0.0_cp
            this%VXip(n_phi_max+2,nThetaB)=0.0_cp
         end do  ! theta loop
      end if     ! chemical composition equation required ?

      if ( l_precession .and. nBc == 0 ) then
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB
            nTheta=nTheta+1
            nThetaNHS=(nTheta+1)/2
            posnalp=-two*oek*po*sin(prec_angle)*osn1(nThetaNHS)
            cnt=cosTheta(nTheta)
            do nPhi=1,n_phi_max
               this%PCr(nPhi,nThetaB)=posnalp*r(nR)*(cos(oek*time+phi(nPhi))* &
               &                                  this%vpc(nPhi,nThetaB)*cnt  &
               &            +sin(oek*time+phi(nPhi))*this%vtc(nPhi,nThetaB))
               this%PCt(nPhi,nThetaB)= -posnalp*or2(nR)*(                     &
               &               cos(oek*time+phi(nPhi))*this%vpc(nPhi,nThetaB) &
               &      +sin(oek*time+phi(nPhi))*or1(nR)*this%vrc(nPhi,nThetaB) )
               this%PCp(nPhi,nThetaB)= posnalp*cos(oek*time+phi(nPhi))*       &
               &              or2(nR)*(      this%vtc(nPhi,nThetaB)-          &
               &                     or1(nR)*this%vrc(nPhi,nThetaB)*cnt)
            end do
            this%PCr(n_phi_max+1,nThetaB)=0.0_cp
            this%PCr(n_phi_max+2,nThetaB)=0.0_cp
            this%PCt(n_phi_max+1,nThetaB)=0.0_cp
            this%PCt(n_phi_max+2,nThetaB)=0.0_cp
            this%PCp(n_phi_max+1,nThetaB)=0.0_cp
            this%PCp(n_phi_max+2,nThetaB)=0.0_cp
         end do ! theta loop
      end if ! precession term required ?

      if ( l_centrifuge .and. nBc ==0 ) then
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB
            nTheta=nTheta+1
            nThetaNHS=(nTheta+1)/2
            snt=sinTheta(nTheta)
            cnt=cosTheta(nTheta)
            rsnt=r(nR)*snt
            do nPhi=1,n_phi_max
               if ( l_anel ) then
                  this%CAr(nPhi,nThetaB) = dilution_fac*rsnt*snt* &
                       &  ( -ra*opr*this%sc(nPhi,nThetaB) )
                       !-- neglect pressure contribution
                       !& + polind*DissNb*oek*opressure0(nR)*this%pc(nPhi,nThetaB) )
                  this%CAt(nPhi,nThetaB) = dilution_fac*rsnt*cnt* &
                       &  ( -ra*opr*this%sc(nPhi,nThetaB) )
                       !-- neglect pressure contribution
                       !& + polind*DissNb*oek*opressure0(nR)*this%pc(nPhi,nThetaB) )
               else
                  this%CAr(nPhi,nThetaB) = -dilution_fac*rsnt*snt*ra*opr*this%sc(nPhi,nThetaB)
                  this%CAt(nPhi,nThetaB) = -dilution_fac*rsnt*cnt*ra*opr*this%sc(nPhi,nThetaB)
               end if
            end do
            this%CAr(n_phi_max+1,nThetaB)=0.0_cp
            this%CAr(n_phi_max+2,nThetaB)=0.0_cp
            this%CAt(n_phi_max+1,nThetaB)=0.0_cp
            this%CAt(n_phi_max+2,nThetaB)=0.0_cp
         end do
      end if ! centrifuge

      if ( l_mag_nl ) then

         if ( nBc == 0 .and. nR>n_r_LCR ) then

            !------ Get (V x B) , the curl of this is the dynamo term:
            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               or4sn2=or4(nR)*osn2(nThetaNHS)

               do nPhi=1,n_phi_max
                  this%VxBr(nPhi,nThetaB)=  orho1(nR)*osn2(nThetaNHS) * (        &
                  &              this%vtc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) - &
                  &              this%vpc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) )
               end do
               this%VxBr(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBr(n_phi_max+2,nThetaB)=0.0_cp

               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)=  orho1(nR)*or4sn2 * (        &
                  &     this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) - &
                  &     this%vrc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) )
               end do
               this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBt(n_phi_max+2,nThetaB)=0.0_cp

               do nPhi=1,n_phi_max
                  this%VxBp(nPhi,nThetaB)=   orho1(nR)*or4sn2 * (        &
                  &      this%vrc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) - &
                  &      this%vtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) )
               end do
               this%VxBp(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBp(n_phi_max+2,nThetaB)=0.0_cp
            end do   ! theta loop

         else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary

            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)=  or4sn2 * orho1(nR) * &
                  &    this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
                  this%VxBp(nPhi,nThetaB)= -or4sn2 * orho1(nR) * &
                  &    this%vtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
               end do
               this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBt(n_phi_max+2,nThetaB)=0.0_cp
               this%VxBp(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBp(n_phi_max+2,nThetaB)=0.0_cp
            end do

         else if ( nBc == 2 ) then  ! rigid boundary :

            !----- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):

            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)= or4sn2 * orho1(nR) * &
                  &    this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
                  this%VxBp(nPhi,nThetaB)= 0.0_cp
               end do
               this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBt(n_phi_max+2,nThetaB)=0.0_cp
               this%VxBp(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBp(n_phi_max+2,nThetaB)=0.0_cp
            end do

         end if  ! boundary ?

      end if ! l_mag_nl ?

      if ( l_anel .and. nBc == 0 ) then
         !------ Get viscous heating
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB ! loop over theta points in block
            nTheta   =nTheta+1
            nThetaNHS=(nTheta+1)/2
            csn2     =cosn2(nThetaNHS)
            if ( mod(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta

            do nPhi=1,n_phi_max
               this%ViscHeat(nPhi,nThetaB)=      or4(nR)*                  &
               &                     orho1(nR)*otemp1(nR)*visc(nR)*(       &
               &     two*(                     this%dvrdrc(nPhi,nThetaB) - & ! (1)
               &     (two*or1(nR)+beta(nR))*this%vrc(nphi,nThetaB) )**2  + &
               &     two*( csn2*                  this%vtc(nPhi,nThetaB) + &
               &                               this%dvpdpc(nphi,nThetaB) + &
               &                               this%dvrdrc(nPhi,nThetaB) - & ! (2)
               &     or1(nR)*               this%vrc(nPhi,nThetaB) )**2  + &
               &     two*(                     this%dvpdpc(nphi,nThetaB) + &
               &           csn2*                  this%vtc(nPhi,nThetaB) + & ! (3)
               &     or1(nR)*               this%vrc(nPhi,nThetaB) )**2  + &
               &          ( two*               this%dvtdpc(nPhi,nThetaB) + &
               &                                 this%cvrc(nPhi,nThetaB) - & ! (6)
               &      two*csn2*             this%vpc(nPhi,nThetaB) )**2  + &
               &                                 osn2(nThetaNHS) * (       &
               &         ( r(nR)*              this%dvtdrc(nPhi,nThetaB) - &
               &           (two+beta(nR)*r(nR))*  this%vtc(nPhi,nThetaB) + & ! (4)
               &     or1(nR)*            this%dvrdtc(nPhi,nThetaB) )**2  + &
               &         ( r(nR)*              this%dvpdrc(nPhi,nThetaB) - &
               &           (two+beta(nR)*r(nR))*  this%vpc(nPhi,nThetaB) + & ! (5)
               &     or1(nR)*            this%dvrdpc(nPhi,nThetaB) )**2 )- &
               &    two*third*(  beta(nR)*        this%vrc(nPhi,nThetaB) )**2 )
            end do
            this%ViscHeat(n_phi_max+1,nThetaB)=0.0_cp
            this%ViscHeat(n_phi_max+2,nThetaB)=0.0_cp
         end do ! theta loop

         if ( l_mag_nl .and. nR>n_r_LCR ) then
            !------ Get ohmic losses
            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB ! loop over theta points in block
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               do nPhi=1,n_phi_max
                  this%OhmLoss(nPhi,nThetaB)= or2(nR)*otemp1(nR)*lambda(nR)*  &
                  &    ( or2(nR)*                this%cbrc(nPhi,nThetaB)**2 + &
                  &      osn2(nThetaNHS)*        this%cbtc(nPhi,nThetaB)**2 + &
                  &      osn2(nThetaNHS)*        this%cbpc(nPhi,nThetaB)**2  )
               end do
               this%OhmLoss(n_phi_max+1,nThetaB)=0.0_cp
               this%OhmLoss(n_phi_max+2,nThetaB)=0.0_cp
            end do ! theta loop

         end if ! if l_mag_nl ?

      end if  ! Viscous heating and Ohmic losses ?

      if ( lRmsCalc ) then
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB ! loop over theta points in block
            nTheta   =nTheta+1
            nThetaNHS=(nTheta+1)/2
            snt=sinTheta(nTheta)
            cnt=cosTheta(nTheta)
            rsnt=r(nR)*snt
            do nPhi=1,n_phi_max
               this%dpdtc(nPhi,nThetaB)=this%dpdtc(nPhi,nThetaB)*or1(nR)* &
               &                        osn2(nThetaNHS)
               this%dpdpc(nPhi,nThetaB)=this%dpdpc(nPhi,nThetaB)*or1(nR)* &
               &                        osn2(nThetaNHS)
               this%CFt2(nPhi,nThetaB)=-2*CorFac *cnt*this%vpc(nPhi,nThetaB)* &
               &                       or1(nR)*osn2(nThetaNHS)
               this%CFp2(nPhi,nThetaB)=2*CorFac * (                      &
               &                     cnt*this%vtc(nPhi,nThetaB)/rsnt +   &
               &                     or2(nR)*snt*this%vrc(nPhi,nThetaB) )/snt
               if ( l_conv_nl ) then
                  this%Advt2(nPhi,nThetaB)=r(nR)*this%Advt(nPhi,nThetaB)
                  this%Advp2(nPhi,nThetaB)=r(nR)*this%Advp(nPhi,nThetaB)
               end if
               if ( l_mag_LF .and. nR > n_r_LCR ) then
                  this%LFt2(nPhi,nThetaB)=r(nR)*this%LFt(nPhi,nThetaB)
                  this%LFp2(nPhi,nThetaB)=r(nR)*this%LFp(nPhi,nThetaB)
               end if
               this%dpdtc(n_phi_max+1,nThetaB)=0.0_cp
               this%dpdtc(n_phi_max+2,nThetaB)=0.0_cp
               this%dpdpc(n_phi_max+1,nThetaB)=0.0_cp
               this%dpdpc(n_phi_max+2,nThetaB)=0.0_cp
               this%CFt2(n_phi_max+1,nThetaB)=0.0_cp
               this%CFt2(n_phi_max+2,nThetaB)=0.0_cp
               this%CFp2(n_phi_max+1,nThetaB)=0.0_cp
               this%CFp2(n_phi_max+2,nThetaB)=0.0_cp
               if ( l_conv_nl ) then
                  this%Advt2(n_phi_max+1,nThetaB)=0.0_cp
                  this%Advt2(n_phi_max+2,nThetaB)=0.0_cp
                  this%Advp2(n_phi_max+1,nThetaB)=0.0_cp
                  this%Advp2(n_phi_max+2,nThetaB)=0.0_cp
               end if
               if ( l_mag_nl .and. nR > n_r_LCR ) then
                  this%LFt2(n_phi_max+1,nThetaB)=0.0_cp
                  this%LFt2(n_phi_max+2,nThetaB)=0.0_cp
                  this%LFp2(n_phi_max+1,nThetaB)=0.0_cp
                  this%LFp2(n_phi_max+2,nThetaB)=0.0_cp
               end if
            end do
         end do

         if ( l_adv_curl ) then
            nTheta=nThetaLast
            do nThetaB=1,sizeThetaB ! loop over theta points in block
               nTheta   =nTheta+1
               nThetaNHS=(nTheta+1)/2
               snt=sinTheta(nTheta)
               csn2     =cosn2(nThetaNHS)
               if ( mod(nTheta,2) == 0 ) csn2=-csn2 ! South, odd function in theta
               do nPhi=1,n_phi_max
                  this%dpdtc(nPhi,nThetaB)=this%dpdtc(nPhi,nThetaB)-              &
                  &            osn2(nThetaNHS)*or3(nR)*( or2(nR)*                 &
                  &            this%vrc(nPhi,nThetaB)*this%dvrdtc(nPhi,nThetaB) - &
                  &            this%vtc(nPhi,nThetaB)*(this%dvrdrc(nPhi,nThetaB)+ &
                  &            this%dvpdpc(nPhi,nThetaB)+csn2 *                   &
                  &            this%vtc(nPhi,nThetaB))+ this%vpc(nPhi,nThetaB)*(  &
                  &            this%cvrc(nPhi,nThetaB)+this%dvtdpc(nPhi,nThetaB)- &
                  &            csn2*this%vpc(nPhi,nThetaB)) )
                  this%dpdpc(nPhi,nThetaB)=this%dpdpc(nPhi,nThetaB)-              &
                  &            osn2(nThetaNHS)*or3(nR)*( or2(nR)*                 &
                  &            this%vrc(nPhi,nThetaB)*this%dvrdpc(nPhi,nThetaB) + &
                  &            this%vtc(nPhi,nThetaB)*this%dvtdpc(nPhi,nThetaB) + &
                  &            this%vpc(nPhi,nThetaB)*this%dvpdpc(nPhi,nThetaB) )
                  if ( l_conv_nl ) then
                     this%Advt2(nPhi,nThetaB)=this%Advt2(nPhi,nThetaB)-          &
                     &        osn2(nThetaNHS)*or3(nR)*( or2(nR)*                 &
                     &        this%vrc(nPhi,nThetaB)*this%dvrdtc(nPhi,nThetaB) - &
                     &        this%vtc(nPhi,nThetaB)*(this%dvrdrc(nPhi,nThetaB)+ &
                     &          this%dvpdpc(nPhi,nThetaB)+csn2 *                 &
                     &        this%vtc(nPhi,nThetaB))+ this%vpc(nPhi,nThetaB)*(  &
                     &        this%cvrc(nPhi,nThetaB)+this%dvtdpc(nPhi,nThetaB)- &
                     &            csn2*this%vpc(nPhi,nThetaB)) )
                     this%Advp2(nPhi,nThetaB)=this%Advp2(nPhi,nThetaB)-          &
                     &        osn2(nThetaNHS)*or3(nR)*( or2(nR)*                 &
                     &        this%vrc(nPhi,nThetaB)*this%dvrdpc(nPhi,nThetaB) + &
                     &        this%vtc(nPhi,nThetaB)*this%dvtdpc(nPhi,nThetaB) + &
                     &        this%vpc(nPhi,nThetaB)*this%dvpdpc(nPhi,nThetaB) )
                  end if

                  !- dpkin/dr = 1/2 d (u^2) / dr = ur*dur/dr+ut*dut/dr+up*dup/dr
                  this%dpkindrc(nPhi,nThetaB)=or4(nR)*this%vrc(nPhi,nThetaB)*(   &
                  &                         this%dvrdrc(nPhi,nThetaB)-           &
                  &                         two*or1(nR)*this%vrc(nPhi,nThetaB)) +&
                  &                         or2(nR)*osn2(nThetaNHS)*(            &
                  &                                 this%vtc(nPhi,nThetaB)*(     &
                  &                         this%dvtdrc(nPhi,nThetaB)-           &
                  &                         or1(nR)*this%vtc(nPhi,nThetaB) ) +   &
                  &                                 this%vpc(nPhi,nThetaB)*(     &
                  &                         this%dvpdrc(nPhi,nThetaB)-           &
                  &                         or1(nR)*this%vpc(nPhi,nThetaB) ) )
               end do
               this%dpkindrc(n_phi_max+1,nThetaB)=0.0_cp
               this%dpkindrc(n_phi_max+2,nThetaB)=0.0_cp
            end do
         end if
      end if

      if ( l_RMS .and. tscheme%istage == 1 ) then
         O_dt = 1.0_cp/tscheme%dt(1)
         nTheta=nThetaLast
         do nThetaB=1,sizeThetaB ! loop over theta points in block
            nTheta   =nTheta+1
            nThetaNHS=(nTheta+1)/2
            do nPhi=1,n_phi_max
               this%dtVr(nPhi,nThetaB)=O_dt*or2(nR)*(this%vrc(nPhi,nThetaB)- &
               &                             vr_old(nPhi,nTheta,nR))
               this%dtVt(nPhi,nThetaB)=O_dt*or1(nR)*(this%vtc(nPhi,nThetaB)- &
               &                       vt_old(nPhi,nTheta,nR))*osn2(nThetaNHS)
               this%dtVp(nPhi,nThetaB)=O_dt*or1(nR)*(this%vpc(nPhi,nThetaB)- &
               &                       vp_old(nPhi,nTheta,nR))*osn2(nThetaNHS)

               vr_old(nPhi,nTheta,nR)=this%vrc(nPhi,nThetaB)
               vt_old(nPhi,nTheta,nR)=this%vtc(nPhi,nThetaB)
               vp_old(nPhi,nTheta,nR)=this%vpc(nPhi,nThetaB)
            end do
            this%dtVr(n_phi_max+1,nThetaB)=0.0_cp
            this%dtVr(n_phi_max+2,nThetaB)=0.0_cp
            this%dtVt(n_phi_max+1,nThetaB)=0.0_cp
            this%dtVt(n_phi_max+2,nThetaB)=0.0_cp
            this%dtVp(n_phi_max+1,nThetaB)=0.0_cp
            this%dtVp(n_phi_max+2,nThetaB)=0.0_cp
         end do
      end if

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
