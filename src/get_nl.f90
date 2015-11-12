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
   use truncation, only: nrp, n_phi_max
   use radial_functions, only: or2, orho1, beta, otemp1, visc, r, &
                               lambda, or4, or1
   use physical_parameters, only: LFfac, n_r_LCR, CorFac
   use blocking, only: nfs, sizeThetaB
   use horizontal_data, only: osn2, cosn2, sinTheta, cosTheta
   use constants, only: two, third
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space: 
      real(cp), allocatable :: Advr(:,:), Advt(:,:), Advp(:,:)
      real(cp), allocatable :: LFr(:,:),  LFt(:,:),  LFp(:,:)
      real(cp), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)
      real(cp), allocatable :: VSr(:,:),  VSt(:,:),  VSp(:,:)
      real(cp), allocatable :: ViscHeat(:,:), OhmLoss(:,:)

      !---- RMS calculations
      real(cp), allocatable :: Advt2(:,:), Advp2(:,:)
      real(cp), allocatable :: LFt2(:,:), LFp2(:,:)
      real(cp), allocatable :: CFt2(:,:), CFp2(:,:)
      real(cp), allocatable :: p1(:,:), p2(:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), pointer :: vrc(:,:), vtc(:,:), vpc(:,:)
      real(cp), pointer :: dvrdrc(:,:), dvtdrc(:,:), dvpdrc(:,:)
      real(cp), pointer :: cvrc(:,:), sc(:,:), drSc(:,:)
      real(cp), pointer :: dvrdtc(:,:), dvrdpc(:,:)
      real(cp), pointer :: dvtdpc(:,:), dvpdpc(:,:)
      real(cp), pointer :: brc(:,:), btc(:,:), bpc(:,:)
      real(cp), pointer :: cbrc(:,:), cbtc(:,:), cbpc(:,:)
      real(cp), pointer :: pc(:,:)
      real(cp), pointer :: dsdtc(:,:), dsdpc(:,:)

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

contains

   subroutine initialize(this)

      class(grid_space_arrays_t) :: this
      integer :: size_in_bytes

      allocate( this%Advr(nrp,nfs) )
      allocate( this%Advt(nrp,nfs) )
      allocate( this%Advp(nrp,nfs) )
      allocate( this%LFr(nrp,nfs) )
      allocate( this%LFt(nrp,nfs) )
      allocate( this%LFp(nrp,nfs) )
      allocate( this%VxBr(nrp,nfs) )
      allocate( this%VxBt(nrp,nfs) )
      allocate( this%VxBp(nrp,nfs) )
      allocate( this%VSr(nrp,nfs) )
      allocate( this%VSt(nrp,nfs) )
      allocate( this%VSp(nrp,nfs) )
      allocate( this%ViscHeat(nrp,nfs) )
      allocate( this%OhmLoss(nrp,nfs) )
      size_in_bytes=14*nrp*nfs*SIZEOF_DEF_REAL

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
      size_in_bytes=size_in_bytes + 21*nrp*nfs*SIZEOF_DEF_REAL

      !-- RMS Calculations
      allocate ( this%Advt2(nrp,nfs) )
      allocate ( this%Advp2(nrp,nfs) )
      allocate ( this%LFt2(nrp,nfs) )
      allocate ( this%LFp2(nrp,nfs) )
      allocate ( this%CFt2(nrp,nfs) )
      allocate ( this%CFp2(nrp,nfs) )
      allocate ( this%p1(nrp,nfs) )
      allocate ( this%p2(nrp,nfs) )
      size_in_bytes=size_in_bytes + 8*nrp*nfs*SIZEOF_DEF_REAL
      !write(*,"(A,I15,A)") "grid_space_arrays: allocated ",size_in_bytes,"B."

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_space_arrays_t) :: this
      integer :: size_in_bytes

      deallocate( this%Advr )
      deallocate( this%Advt )
      deallocate( this%Advp )
      deallocate( this%LFr )
      deallocate( this%LFt )
      deallocate( this%LFp )
      deallocate( this%VxBr )
      deallocate( this%VxBt )
      deallocate( this%VxBp )
      deallocate( this%VSr )
      deallocate( this%VSt )
      deallocate( this%VSp )
      deallocate( this%ViscHeat )
      deallocate( this%OhmLoss )
      size_in_bytes=14*nrp*nfs*SIZEOF_DEF_REAL

      !----- Fields calculated from these help arrays by legtf:
      deallocate( this%vrc,this%vtc,this%vpc )
      deallocate( this%dvrdrc,this%dvtdrc )
      deallocate( this%dvpdrc,this%cvrc )
      deallocate( this%dvrdtc,this%dvrdpc )
      deallocate( this%dvtdpc,this%dvpdpc )
      deallocate( this%brc,this%btc,this%bpc )
      deallocate( this%cbrc,this%cbtc,this%cbpc )
      deallocate( this%sc,this%drSc )
      deallocate( this%pc )
      deallocate( this%dsdtc, this%dsdpc )
      size_in_bytes=size_in_bytes + 21*nrp*nfs*SIZEOF_DEF_REAL

      !-- RMS Calculations
      deallocate ( this%Advt2 )
      deallocate ( this%Advp2 )
      deallocate ( this%LFt2 )
      deallocate ( this%LFp2 )
      deallocate ( this%CFt2 )
      deallocate ( this%CFp2 )
      deallocate ( this%p1 )
      deallocate ( this%p2 )
      size_in_bytes=size_in_bytes + 8*nrp*nfs*SIZEOF_DEF_REAL
      write(*,"(A,I15,A)") "grid_space_arrays: deallocated ",size_in_bytes,"B."

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
   subroutine get_nl_shtns(this, nR, nBc)
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
      integer, intent(in) :: nR
      integer, intent(in) :: nBc

      !-- Local variables:
      integer :: nThetaB, nThetaNHS
      integer :: nPhi
      real(cp) :: or2sn2, or4sn2, csn2


      if ( l_mag_LF .and. nBc == 0 .and. nR>n_r_LCR ) then
         !------ Get the Lorentz force:
         !$OMP PARALLEL DO default(none) &
         !$OMP& private(nThetaB, nPhi, nThetaNHS, or4sn2) &
         !$OMP& shared(n_phi_max, sizeThetaB) &
         !$OMP& shared(or4, nR, osn2, this, LFfac)
         do nThetaB=1,sizeThetaB

            nThetaNHS=(nThetaB+1)/2
            or4sn2   =or4(nR)*osn2(nThetaNHS)

            do nPhi=1,n_phi_max
               !---- LFr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
               this%LFr(nPhi,nThetaB)=  LFfac*osn2(nThetaNHS) * (        &
                        this%cbtc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) - &
                        this%cbpc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) )
            end do
            ! this%LFr(n_phi_max+1,nThetaB)=0.0_cp
            ! this%LFr(n_phi_max+2,nThetaB)=0.0_cp

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            do nPhi=1,n_phi_max
               this%LFt(nPhi,nThetaB)=           LFfac*or4sn2 * (        &
                        this%cbpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) - &
                        this%cbrc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) )
            end do
            ! this%LFt(n_phi_max+1,nThetaB)=0.0_cp
            ! this%LFt(n_phi_max+2,nThetaB)=0.0_cp
            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            do nPhi=1,n_phi_max
               this%LFp(nPhi,nThetaB)=           LFfac*or4sn2 * (        &
                        this%cbrc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) - &
                        this%cbtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) )
            end do
            ! this%LFp(n_phi_max+1,nThetaB)=0.0_cp
            ! this%LFp(n_phi_max+2,nThetaB)=0.0_cp

         end do   ! theta loop
         !$OMP END PARALLEL DO
      end if      ! Lorentz force required ?

      if ( l_conv_nl .and. nBc == 0 ) then

         !------ Get Advection:
         !$OMP PARALLEL DO default(none) &
         !$OMP& private(nThetaB, nPhi, nThetaNHS, or4sn2, csn2) &
         !$OMP& shared(this, or2, beta, orho1, r, sizeThetaB, or4, nR, osn2) &
         !$OMP& shared(n_phi_max, cosn2, or1)
         do nThetaB=1,sizeThetaB ! loop over theta points in block
            nThetaNHS=(nThetaB+1)/2
            or4sn2   =or4(nR)*osn2(nThetaNHS)
            csn2     =cosn2(nThetaNHS)
            if ( mod(nThetaB,2) == 0 ) csn2=-csn2 ! South, odd function in theta

            do nPhi=1,n_phi_max
               this%Advr(nPhi,nThetaB)=          -or2(nR)*orho1(nR) * (  &
                                                this%vrc(nPhi,nThetaB) * &
                                     (       this%dvrdrc(nPhi,nThetaB) - &
                    ( two*or1(nR)+beta(nR) )*this%vrc(nPhi,nThetaB) ) +  &
                                               osn2(nThetaNHS) * (       &
                                                this%vtc(nPhi,nThetaB) * &
                                     (       this%dvrdtc(nPhi,nThetaB) - &
                                  r(nR)*      this%vtc(nPhi,nThetaB) ) + &
                                                this%vpc(nPhi,nThetaB) * &
                                     (       this%dvrdpc(nPhi,nThetaB) - &
                                    r(nR)*      this%vpc(nPhi,nThetaB) ) ) )
            end do
            ! this%Advr(n_phi_max+1,nThetaB)=0.0_cp
            ! this%Advr(n_phi_max+2,nThetaB)=0.0_cp
            do nPhi=1,n_phi_max
               this%Advt(nPhi,nThetaB)=         or4sn2*orho1(nR) * (  &
                                            -this%vrc(nPhi,nThetaB) * &
                                      (   this%dvtdrc(nPhi,nThetaB) - &
                                beta(nR)*this%vtc(nPhi,nThetaB) )   + &
                                             this%vtc(nPhi,nThetaB) * &
                                      ( csn2*this%vtc(nPhi,nThetaB) + &
                                          this%dvpdpc(nPhi,nThetaB) + &
                                      this%dvrdrc(nPhi,nThetaB) )   + &
                                             this%vpc(nPhi,nThetaB) * &
                                      ( csn2*this%vpc(nPhi,nThetaB) - &
                                          this%dvtdpc(nPhi,nThetaB) )  )
            end do
            ! this%Advt(n_phi_max+1,nThetaB)=0.0_cp
            ! this%Advt(n_phi_max+2,nThetaB)=0.0_cp
            do nPhi=1,n_phi_max
               this%Advp(nPhi,nThetaB)=         or4sn2*orho1(nR) * (  &
                                            -this%vrc(nPhi,nThetaB) * &
                                        ( this%dvpdrc(nPhi,nThetaB) - &
                                beta(nR)*this%vpc(nPhi,nThetaB) )   - &
                                             this%vtc(nPhi,nThetaB) * &
                                        ( this%dvtdpc(nPhi,nThetaB) + &
                                        this%cvrc(nPhi,nThetaB) )   - &
                       this%vpc(nPhi,nThetaB) * this%dvpdpc(nPhi,nThetaB) )
            end do
            ! this%Advp(n_phi_max+1,nThetaB)=0.0_cp
            ! this%Advp(n_phi_max+2,nThetaB)=0.0_cp
         end do ! theta loop
         !$OMP END PARALLEL DO

      end if  ! Navier-Stokes nonlinear advection term ?

      if ( l_heat_nl .and. nBc == 0 ) then
         !------ Get V S, the divergence of the is entropy advection:
         !$OMP PARALLEL DO default(none) &
         !$OMP& private(nThetaB, nPhi, nThetaNHS, or2sn2) &
         !$OMP& shared(this, or2, osn2, sizeThetaB, nR, n_phi_max)
         do nThetaB=1,sizeThetaB
            nThetaNHS=(nThetaB+1)/2
            or2sn2=or2(nR)*osn2(nThetaNHS)
            do nPhi=1,n_phi_max     ! calculate v*s components
               this%VSr(nPhi,nThetaB)= &
                    this%vrc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
               this%VSt(nPhi,nThetaB)= &
                    or2sn2*this%vtc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
               this%VSp(nPhi,nThetaB)= &
                    or2sn2*this%vpc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
            end do
            ! this%VSr(n_phi_max+1,nThetaB)=0.0_cp
            ! this%VSr(n_phi_max+2,nThetaB)=0.0_cp
            ! this%VSt(n_phi_max+1,nThetaB)=0.0_cp
            ! this%VSt(n_phi_max+2,nThetaB)=0.0_cp
            ! this%VSp(n_phi_max+1,nThetaB)=0.0_cp
            ! this%VSp(n_phi_max+2,nThetaB)=0.0_cp
         end do  ! theta loop
         !$OMP END PARALLEL DO
      end if     ! heat equation required ?

      if ( l_mag_nl ) then

         if ( nBc == 0 .and. nR>n_r_LCR ) then

            !------ Get (V x B) , the curl of this is the dynamo term:
            !$OMP PARALLEL DO default(none) &
            !$OMP& private(nThetaB, nPhi, nThetaNHS, or4sn2) &
            !$OMP& shared(this, or4, osn2, nR, sizeThetaB, n_phi_max) &
            !$OMP& shared(orho1)
            do nThetaB=1,sizeThetaB
               nThetaNHS=(nThetaB+1)/2
               or4sn2=or4(nR)*osn2(nThetaNHS)

               do nPhi=1,n_phi_max
                  this%VxBr(nPhi,nThetaB)=  orho1(nR)*osn2(nThetaNHS) * (        &
                                 this%vtc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) - &
                                 this%vpc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) )
               end do
               ! this%VxBr(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBr(n_phi_max+2,nThetaB)=0.0_cp

               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)=  orho1(nR)*or4sn2 * (        &
                        this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) - &
                        this%vrc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) )
               end do
               ! this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBt(n_phi_max+2,nThetaB)=0.0_cp

               do nPhi=1,n_phi_max
                  this%VxBp(nPhi,nThetaB)=   orho1(nR)*or4sn2 * (        &
                         this%vrc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) - &
                         this%vtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) )
               end do
               ! this%VxBp(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBp(n_phi_max+2,nThetaB)=0.0_cp
            end do   ! theta loop
            !$OMP END PARALLEL DO

         else if ( nBc == 1 .or. nR<=n_r_LCR ) then ! stress free boundary

            !$OMP PARALLEL DO default(none) &
            !$OMP& private(nThetaB, nPhi, nThetaNHS, or4sn2)            &
            !$OMP& shared(this, sizeThetaB, n_phi_max, or4, osn2)       &
            !$OMP& shared(nR, orho1)
            do nThetaB=1,sizeThetaB
               nThetaNHS=(nThetaB+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)=  or4sn2 * orho1(nR) * &
                       this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
                  this%VxBp(nPhi,nThetaB)= -or4sn2 * orho1(nR) * &
                       this%vtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
               end do
               ! this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBt(n_phi_max+2,nThetaB)=0.0_cp
               ! this%VxBp(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBp(n_phi_max+2,nThetaB)=0.0_cp
            end do
            !$OMP END PARALLEL DO

         else if ( nBc == 2 ) then  ! rigid boundary :

            !----- Only vp /= 0 at boundary allowed (rotation of boundaries about z-axis):
            !$OMP PARALLEL DO default(none) &
            !$OMP& private(nThetaB, nPhi, nThetaNHS, or4sn2) &
            !$OMP& shared(this, sizeThetaB, or4, orho1, n_phi_max, nR, osn2)
            do nThetaB=1,sizeThetaB
               nThetaNHS=(nThetaB+1)/2
               or4sn2   =or4(nR)*osn2(nThetaNHS)
               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)= or4sn2 * orho1(nR) * &
                       this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
                  this%VxBp(nPhi,nThetaB)= 0.0_cp
               end do
               ! this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBt(n_phi_max+2,nThetaB)=0.0_cp
               ! this%VxBp(n_phi_max+1,nThetaB)=0.0_cp
               ! this%VxBp(n_phi_max+2,nThetaB)=0.0_cp
            end do
            !$OMP END PARALLEL DO

         end if  ! boundary ?

      end if ! l_mag_nl ?

      if ( l_anel .and. nBc == 0 ) then
         !------ Get viscous heating
         !$OMP PARALLEL DO default(none) &
         !$OMP& private(nThetaB, nPhi, nThetaNHS, csn2) &
         !$OMP& shared(this, nR, sizeThetaB, n_phi_max, cosn2, otemp1) &
         !$OMP& shared(beta, orho1, or4, osn2, visc, r, or1)
         do nThetaB=1,sizeThetaB ! loop over theta points in block
            nThetaNHS=(nThetaB+1)/2
            csn2     =cosn2(nThetaNHS)
            if ( mod(nThetaB,2) == 0 ) csn2=-csn2 ! South, odd function in theta

            do nPhi=1,n_phi_max
               this%ViscHeat(nPhi,nThetaB)=      or4(nR)*                  &
                                     orho1(nR)*otemp1(nR)*visc(nR)*(       &
                     two*(                     this%dvrdrc(nPhi,nThetaB) - & ! (1)
                     (two*or1(nR)+beta(nR))*this%vrc(nphi,nThetaB) )**2  + &
                     two*( csn2*                  this%vtc(nPhi,nThetaB) + &
                                               this%dvpdpc(nphi,nThetaB) + &
                                               this%dvrdrc(nPhi,nThetaB) - & ! (2)
                     or1(nR)*               this%vrc(nPhi,nThetaB) )**2  + &
                     two*(                     this%dvpdpc(nphi,nThetaB) + &
                           csn2*                  this%vtc(nPhi,nThetaB) + & ! (3)
                     or1(nR)*               this%vrc(nPhi,nThetaB) )**2  + &
                          ( two*               this%dvtdpc(nPhi,nThetaB) + &
                                                 this%cvrc(nPhi,nThetaB) - & ! (6)
                      two*csn2*             this%vpc(nPhi,nThetaB) )**2  + &
                                                 osn2(nThetaNHS) * (       &
                         ( r(nR)*              this%dvtdrc(nPhi,nThetaB) - &
                           (two+beta(nR)*r(nR))*  this%vtc(nPhi,nThetaB) + & ! (4)
                     or1(nR)*            this%dvrdtc(nPhi,nThetaB) )**2  + &
                         ( r(nR)*              this%dvpdrc(nPhi,nThetaB) - &
                           (two+beta(nR)*r(nR))*  this%vpc(nPhi,nThetaB) + & ! (5)
                     or1(nR)*            this%dvrdpc(nPhi,nThetaB) )**2 )- &
                    two*third*(  beta(nR)*        this%vrc(nPhi,nThetaB) )**2 )
            end do
            ! this%ViscHeat(n_phi_max+1,nThetaB)=0.0_cp
            ! this%ViscHeat(n_phi_max+2,nThetaB)=0.0_cp
         end do ! theta loop
         !$OMP END PARALLEL DO

         if ( l_mag_nl .and. nR>n_r_LCR ) then
            !------ Get ohmic losses
            !$OMP PARALLEL DO default(none) &
            !$OMP& private(nThetaB, nPhi, nThetaNHS) &
            !$OMP& shared(this, nR, sizeThetaB, n_phi_max) &
            !$OMP& shared(or2, otemp1, lambda, osn2)
            do nThetaB=1,sizeThetaB ! loop over theta points in block
               nThetaNHS=(nThetaB+1)/2
               do nPhi=1,n_phi_max
                  this%OhmLoss(nPhi,nThetaB)= or2(nR)*otemp1(nR)*lambda(nR)*  &
                       ( or2(nR)*                this%cbrc(nPhi,nThetaB)**2 + &
                         osn2(nThetaNHS)*        this%cbtc(nPhi,nThetaB)**2 + &
                         osn2(nThetaNHS)*        this%cbpc(nPhi,nThetaB)**2  )
               end do
               ! this%OhmLoss(n_phi_max+1,nThetaB)=0.0_cp
               ! this%OhmLoss(n_phi_max+2,nThetaB)=0.0_cp
            end do ! theta loop
            !$OMP END PARALLEL DO

         end if ! if l_mag_nl ?

      end if  ! Viscous heating and Ohmic losses ?

   end subroutine get_nl_shtns
#endif
!----------------------------------------------------------------------------
   subroutine get_nl(this,nR,nBc,nThetaStart,lRmsCalc)
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
      integer, intent(in) :: nR
      integer, intent(in) :: nBc
      integer, intent(in) :: nThetaStart
      logical, intent(in) :: lRmsCalc

      !-- Local variables:
      integer :: nTheta
      integer :: nThetaLast,nThetaB,nThetaNHS
      integer :: nPhi
      real(cp) :: or2sn2,or4sn2,csn2,snt,cnt,rsnt

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
                        this%cbtc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) - &
                        this%cbpc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) )
            end do
            this%LFr(n_phi_max+1,nThetaB)=0.0_cp
            this%LFr(n_phi_max+2,nThetaB)=0.0_cp

            !---- LFt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
            do nPhi=1,n_phi_max
               this%LFt(nPhi,nThetaB)=           LFfac*or4sn2 * (        &
                        this%cbpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) - &
                        this%cbrc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) )
            end do
            this%LFt(n_phi_max+1,nThetaB)=0.0_cp
            this%LFt(n_phi_max+2,nThetaB)=0.0_cp
            !---- LFp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
            do nPhi=1,n_phi_max
               this%LFp(nPhi,nThetaB)=           LFfac*or4sn2 * (        &
                        this%cbrc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) - &
                        this%cbtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) )
            end do
            this%LFp(n_phi_max+1,nThetaB)=0.0_cp
            this%LFp(n_phi_max+2,nThetaB)=0.0_cp

         end do   ! theta loop
      end if      ! Lorentz force required ?

      if ( l_conv_nl .and. (nBc == 0 .or. lRmsCalc) ) then

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
                                                this%vrc(nPhi,nThetaB) * &
                                     (       this%dvrdrc(nPhi,nThetaB) - &
                    ( two*or1(nR)+beta(nR) )*this%vrc(nPhi,nThetaB) ) +  &
                                               osn2(nThetaNHS) * (       &
                                                this%vtc(nPhi,nThetaB) * &
                                     (       this%dvrdtc(nPhi,nThetaB) - &
                                  r(nR)*      this%vtc(nPhi,nThetaB) ) + &
                                                this%vpc(nPhi,nThetaB) * &
                                     (       this%dvrdpc(nPhi,nThetaB) - &
                                    r(nR)*      this%vpc(nPhi,nThetaB) ) ) )
            end do
            this%Advr(n_phi_max+1,nThetaB)=0.0_cp
            this%Advr(n_phi_max+2,nThetaB)=0.0_cp
            do nPhi=1,n_phi_max
               this%Advt(nPhi,nThetaB)=         or4sn2*orho1(nR) * (  &
                                            -this%vrc(nPhi,nThetaB) * &
                                      (   this%dvtdrc(nPhi,nThetaB) - &
                                beta(nR)*this%vtc(nPhi,nThetaB) )   + &
                                             this%vtc(nPhi,nThetaB) * &
                                      ( csn2*this%vtc(nPhi,nThetaB) + &
                                          this%dvpdpc(nPhi,nThetaB) + &
                                      this%dvrdrc(nPhi,nThetaB) )   + &
                                             this%vpc(nPhi,nThetaB) * &
                                      ( csn2*this%vpc(nPhi,nThetaB) - &
                                          this%dvtdpc(nPhi,nThetaB) )  )
            end do
            this%Advt(n_phi_max+1,nThetaB)=0.0_cp
            this%Advt(n_phi_max+2,nThetaB)=0.0_cp
            do nPhi=1,n_phi_max
               this%Advp(nPhi,nThetaB)=         or4sn2*orho1(nR) * (  &
                                            -this%vrc(nPhi,nThetaB) * &
                                        ( this%dvpdrc(nPhi,nThetaB) - &
                                beta(nR)*this%vpc(nPhi,nThetaB) )   - &
                                             this%vtc(nPhi,nThetaB) * &
                                        ( this%dvtdpc(nPhi,nThetaB) + &
                                        this%cvrc(nPhi,nThetaB) )   - &
                       this%vpc(nPhi,nThetaB) * this%dvpdpc(nPhi,nThetaB) )
            end do
            this%Advp(n_phi_max+1,nThetaB)=0.0_cp
            this%Advp(n_phi_max+2,nThetaB)=0.0_cp
         end do ! theta loop

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
                    this%vrc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
               this%VSt(nPhi,nThetaB)= &
                    or2sn2*this%vtc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
               this%VSp(nPhi,nThetaB)= &
                    or2sn2*this%vpc(nPhi,nThetaB)*this%sc(nPhi,nThetaB)
            end do
            this%VSr(n_phi_max+1,nThetaB)=0.0_cp
            this%VSr(n_phi_max+2,nThetaB)=0.0_cp
            this%VSt(n_phi_max+1,nThetaB)=0.0_cp
            this%VSt(n_phi_max+2,nThetaB)=0.0_cp
            this%VSp(n_phi_max+1,nThetaB)=0.0_cp
            this%VSp(n_phi_max+2,nThetaB)=0.0_cp
         end do  ! theta loop
      end if     ! heat equation required ?

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
                                 this%vtc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) - &
                                 this%vpc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) )
               end do
               this%VxBr(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBr(n_phi_max+2,nThetaB)=0.0_cp

               do nPhi=1,n_phi_max
                  this%VxBt(nPhi,nThetaB)=  orho1(nR)*or4sn2 * (        &
                        this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) - &
                        this%vrc(nPhi,nThetaB)*this%bpc(nPhi,nThetaB) )
               end do
               this%VxBt(n_phi_max+1,nThetaB)=0.0_cp
               this%VxBt(n_phi_max+2,nThetaB)=0.0_cp

               do nPhi=1,n_phi_max
                  this%VxBp(nPhi,nThetaB)=   orho1(nR)*or4sn2 * (        &
                         this%vrc(nPhi,nThetaB)*this%btc(nPhi,nThetaB) - &
                         this%vtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB) )
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
                       this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
                  this%VxBp(nPhi,nThetaB)= -or4sn2 * orho1(nR) * &
                       this%vtc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
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
                       this%vpc(nPhi,nThetaB)*this%brc(nPhi,nThetaB)
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
                                     orho1(nR)*otemp1(nR)*visc(nR)*(       &
                     two*(                     this%dvrdrc(nPhi,nThetaB) - & ! (1)
                     (two*or1(nR)+beta(nR))*this%vrc(nphi,nThetaB) )**2  + &
                     two*( csn2*                  this%vtc(nPhi,nThetaB) + &
                                               this%dvpdpc(nphi,nThetaB) + &
                                               this%dvrdrc(nPhi,nThetaB) - & ! (2)
                     or1(nR)*               this%vrc(nPhi,nThetaB) )**2  + &
                     two*(                     this%dvpdpc(nphi,nThetaB) + &
                           csn2*                  this%vtc(nPhi,nThetaB) + & ! (3)
                     or1(nR)*               this%vrc(nPhi,nThetaB) )**2  + &
                          ( two*               this%dvtdpc(nPhi,nThetaB) + &
                                                 this%cvrc(nPhi,nThetaB) - & ! (6)
                      two*csn2*             this%vpc(nPhi,nThetaB) )**2  + &
                                                 osn2(nThetaNHS) * (       &
                         ( r(nR)*              this%dvtdrc(nPhi,nThetaB) - &
                           (two+beta(nR)*r(nR))*  this%vtc(nPhi,nThetaB) + & ! (4)
                     or1(nR)*            this%dvrdtc(nPhi,nThetaB) )**2  + &
                         ( r(nR)*              this%dvpdrc(nPhi,nThetaB) - &
                           (two+beta(nR)*r(nR))*  this%vpc(nPhi,nThetaB) + & ! (5)
                     or1(nR)*            this%dvrdpc(nPhi,nThetaB) )**2 )- &
                    two*third*(  beta(nR)*        this%vrc(nPhi,nThetaB) )**2 )
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
                       ( or2(nR)*                this%cbrc(nPhi,nThetaB)**2 + &
                         osn2(nThetaNHS)*        this%cbtc(nPhi,nThetaB)**2 + &
                         osn2(nThetaNHS)*        this%cbpc(nPhi,nThetaB)**2  )
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
            snt=sinTheta(nTheta)
            cnt=cosTheta(nTheta)
            rsnt=r(nR)*snt
            do nPhi=1,n_phi_max
               this%p1(nPhi,nThetaB)=this%pc(nPhi,nThetaB)/snt
               this%p2(nPhi,nThetaB)=cnt*this%p1(nPhi,nThetaB)
               this%CFt2(nPhi,nThetaB)=-2*CorFac *cnt*this%vpc(nPhi,nThetaB)/rsnt
               this%CFp2(nPhi,nThetaB)=2*CorFac * (                      &
                                     cnt*this%vtc(nPhi,nThetaB)/rsnt +   &
                                     or2(nR)*snt*this%vrc(nPhi,nThetaB) )
               if ( l_conv_nl ) then
                  this%Advt2(nPhi,nThetaB)=rsnt*this%Advt(nPhi,nThetaB)
                  this%Advp2(nPhi,nThetaB)=rsnt*this%Advp(nPhi,nThetaB)
               end if
               if ( l_mag_LF .and. nR > n_r_LCR ) then
                  this%LFt2(nPhi,nThetaB)=rsnt*this%LFt(nPhi,nThetaB)
                  this%LFp2(nPhi,nThetaB)=rsnt*this%LFp(nPhi,nThetaB)
               end if
               this%p1(n_phi_max+1,nThetaB)=0.0_cp
               this%p1(n_phi_max+2,nThetaB)=0.0_cp
               this%p2(n_phi_max+1,nThetaB)=0.0_cp
               this%p2(n_phi_max+2,nThetaB)=0.0_cp
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
      end if

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
