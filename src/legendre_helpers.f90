module leg_helper_mod

   use precision_mod
   use truncation, only: lm_max,l_max
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: or2
   use torsional_oscillations, only: ddzASL
   use Grenoble, only: lGrenoble, b0, db0, ddb0
   use blocking, only: lm2l, lm2m, lm2
   use horizontal_data, only: dLh
   use logic, only: l_conv, l_mag_kin, l_heat, l_mag, l_movie_oc,    &
       &            l_mag_LF, l_fluxProfs
   use fields, only: s_Rloc,ds_Rloc, z_Rloc,dz_Rloc, p_Rloc,dp_Rloc, &
       &             b_Rloc,db_Rloc,ddb_Rloc, aj_Rloc,dj_Rloc,       &
       &             w_Rloc,dw_Rloc,ddw_Rloc, omega_ic,omega_ma
   use const, only: zero, one, two
#ifdef WITH_SHTNS
   use shtns
#endif

   implicit none

   private

   type, public :: leg_helper_t
      ! Help arrays for Legendre transform calculated in legPrepG:
      ! Parallelizatio note: these are the R-distributed versions
      ! of the field scalars.
      complex(cp), allocatable :: dLhw(:), dLhdw(:), dLhz(:), dLhb(:), dLhj(:)
      complex(cp), allocatable :: vhG(:), vhC(:), dvhdrG(:), dvhdrC(:)
      complex(cp), allocatable :: bhG(:), bhC(:), cbhG(:), cbhC(:)
      !----- R-distributed versions of scalar fields (see c_fields.f):
      complex(cp), allocatable :: sR(:), dsR(:), preR(:), dpR(:)
      real(cp), allocatable :: zAS(:), dzAS(:), ddzAS(:) ! used in TO
      real(cp) :: omegaIC,omegaMA
      complex(cp), allocatable :: bCMB(:)
#ifdef WITH_SHTNS
      real(kind=8), allocatable :: shtns_s(:)
      real(kind=8), allocatable :: shtns_p(:)
      real(kind=8), allocatable :: shtns_drs(:)

      real(kind=8), allocatable :: shtns_dsdt(:)
      real(kind=8), allocatable :: shtns_dsdp(:)

      real(kind=8), allocatable :: shtns_dpdt(:)
      real(kind=8), allocatable :: shtns_dpdp(:)

      real(kind=8), allocatable :: shtns_vr(:)
      real(kind=8), allocatable :: shtns_vt(:)
      real(kind=8), allocatable :: shtns_vp(:)

      real(kind=8), allocatable :: shtns_dvrdr(:)
      real(kind=8), allocatable :: shtns_dvtdr(:)
      real(kind=8), allocatable :: shtns_dvpdr(:)

      real(kind=8), allocatable :: shtns_dvrdt(:)
      real(kind=8), allocatable :: shtns_dvrdp(:)

      real(kind=8), allocatable :: shtns_cvr(:)

      real(kind=8), allocatable :: shtns_br(:)
      real(kind=8), allocatable :: shtns_bt(:)
      real(kind=8), allocatable :: shtns_bp(:)

      real(kind=8), allocatable :: shtns_cbr(:)
      real(kind=8), allocatable :: shtns_cbt(:)
      real(kind=8), allocatable :: shtns_cbp(:)
#endif
   contains

      procedure :: initialize
      procedure :: legPrepG

   end type leg_helper_t

   public :: legPrep, legPrep_IC

contains

   subroutine initialize(this,lm_max,lm_maxMag,l_max)
#ifdef WITH_SHTNS
      use truncation, only : n_phi_max, n_theta_max
#endif

      class(leg_helper_t) :: this
      integer,intent(in) :: lm_max,lm_maxMag,l_max

      allocate( this%dLhw(lm_max) )
      allocate( this%dLhdw(lm_max) )
      allocate( this%dLhz(lm_max) )
      allocate( this%dLhb(lm_max) )
      allocate( this%dLhj(lm_max) )
      allocate( this%vhG(lm_max) )
      allocate( this%vhC(lm_max) )
      allocate( this%dvhdrG(lm_max) )
      allocate( this%dvhdrC(lm_max) )
      allocate( this%bhG(lm_maxMag) )
      allocate( this%bhC(lm_maxMag) )
      allocate( this%cbhG(lm_maxMag) )
      allocate( this%cbhC(lm_maxMag) )
      !----- R-distributed versions of scalar fields (see c_fields.f):
      allocate( this%sR(lm_max),this%dsR(lm_max) )
      allocate( this%preR(lm_max),this%dpR(lm_max) )
      allocate( this%zAS(l_max+1),this%dzAS(l_max+1),this%ddzAS(l_max+1) ) ! used in TO

      allocate( this%bCMB(lm_maxMag) )

#ifdef WITH_SHTNS
      allocate(this%shtns_s(n_phi_max * n_theta_max))
      allocate(this%shtns_drs(n_phi_max * n_theta_max))
      allocate(this%shtns_p(n_phi_max * n_theta_max))

      allocate(this%shtns_dsdt(n_phi_max * n_theta_max))
      allocate(this%shtns_dsdp(n_phi_max * n_theta_max))

      allocate(this%shtns_dpdt(n_phi_max * n_theta_max))
      allocate(this%shtns_dpdp(n_phi_max * n_theta_max))

      allocate(this%shtns_vr(n_phi_max * n_theta_max))
      allocate(this%shtns_vt(n_phi_max * n_theta_max))
      allocate(this%shtns_vp(n_phi_max * n_theta_max))

      allocate(this%shtns_dvrdr(n_phi_max * n_theta_max))
      allocate(this%shtns_dvtdr(n_phi_max * n_theta_max))
      allocate(this%shtns_dvpdr(n_phi_max * n_theta_max))

      allocate(this%shtns_dvrdt(n_phi_max * n_theta_max))
      allocate(this%shtns_dvrdp(n_phi_max * n_theta_max))

      allocate(this%shtns_cvr(n_phi_max * n_theta_max))

      allocate(this%shtns_br(n_phi_max * n_theta_max))
      allocate(this%shtns_bt(n_phi_max * n_theta_max))
      allocate(this%shtns_bp(n_phi_max * n_theta_max))

      allocate(this%shtns_cbr(n_phi_max * n_theta_max))
      allocate(this%shtns_cbt(n_phi_max * n_theta_max))
      allocate(this%shtns_cbp(n_phi_max * n_theta_max))

      this%shtns_bt=1.0d50
      this%shtns_bp=1.0d50
#endif
   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine legPrepG(this,nR,nBc,lDeriv,lRmsCalc,l_frame, &
        &              lTOnext,lTOnext2,lTOcalc, lFluxProfCalc)
      !
      !  Purpose of this subroutine is to prepare Legendre transforms     
      !  from (r,l,m) space to (r,theta,m) space by calculating           
      !  auxiliary arrays dpdw,dpddw, ....... dLhj which contain          
      !  combinations of harmonic coeffs multiplied with (l,m)-dependend  
      !  factors as well as the radial dependence:
      !
      !     * nBc  =0 standard inner radial grid point                       
      !     * nBc  =1 radial velocity zero, spatial derivs not needed        
      !     * nBc  =2 all velocity comp. zero, spatial derivs not needed     
      !
      !  lDeriv=.true. field derivatives required                        

      class(leg_helper_t) :: this

      !-- Input of variables
      integer, intent(in) :: nR          ! radial level
      integer, intent(in) :: nBc         ! boundary condition
      logical, intent(in) :: lDeriv      ! get also field derivatives !
      logical, intent(in) :: lRmsCalc    ! Rms force balance ?
      logical, intent(in) :: l_frame     ! Calculate movie frame?
      logical, intent(in) :: lTOnext     ! for TO output
      logical, intent(in) :: lTOnext2
      logical, intent(in) :: lTOcalc
      logical, intent(in) :: lFluxProfCalc

      !-- Input of scalar fields in LM-distributed space
      !   These have to be collected and stored in the
      !   R-distributed output variables

      !-- Local variables:
      integer :: lm,l,m
      complex(cp) :: dbd


      if ( nR == n_r_icb ) this%omegaIC=omega_ic
      if ( nR == n_r_cmb ) this%omegaMA=omega_ma
      if ( l_conv .or. l_mag_kin ) then

         if ( l_heat ) then
#ifdef WITH_SHTNS
            call scal_to_spat(s_Rloc(:, nR), this%shtns_s)
            if (lFluxProfCalc) then
               call scal_to_spat(p_Rloc(:, nR), this%shtns_p)
           end if
           call scal_to_grad_spat(s_Rloc(:, nR), this%shtns_dsdt, &
                                  this%shtns_dsdp)
           ! if (lViscBcCalc .or. l_HT) then
           call scal_to_spat(ds_Rloc(:, nR), this%shtns_drs)
           ! endif
#else
            do lm=1,lm_max
               this%sR(lm) =s_Rloc(lm,nR)   ! used for plotting and Rms
               this%dsR(lm)=ds_Rloc(lm,nR)  ! used for plotting and Rms
            end do
#endif
#ifdef WITH_SHTNS
            call scal_to_spat(s_Rloc(:, nR), this%shtns_s)
#endif
         end if
         if ( lTOnext .or. lTOnext2 .or. lTOCalc ) then
            do lm=1,lm_max
               l=lm2l(lm)
               m=lm2m(lm)
               if ( l <= l_max .and. m == 0 ) then
                  this%zAS(l+1)  =real(z_Rloc(lm,nR))   ! used in TO
                  this%dzAS(l+1) =real(dz_Rloc(lm,nR))  ! used in TO (anelastic)
                  this%ddzAS(l+1)=ddzASL(l+1,nR)        ! used in TO
               end if
            end do
         end if
         if ( lRmsCalc .or. l_fluxProfs ) then
            this%preR(1)=zero
            this%dpR(1)=zero
            do lm=2,lm_max
               this%preR(lm)=p_Rloc(lm,nR)    ! used for Rms in get_td (anelastic)
               this%dpR(lm)=dp_Rloc(lm,nR)  ! used for Rms in get_td
            end do
         end if
         if ( l_mag .and. l_frame .and. l_movie_oc .and. nR == n_r_cmb ) then
            this%bCMB(1)=zero ! used in s_store_movie_frame.f
            do lm=2,lm_max
               this%bCMB(lm)=b_Rloc(lm,nR)  ! used for movie output of surface field
            end do
         end if

         if ( nBc /= 2 ) then ! nBc=2 is flag for fixed boundary
#ifdef WITH_SHTNS
            call torpol_to_spat(w_Rloc(:, nR), dw_Rloc(:, nR),  z_Rloc(:, nR), nR, &
                this%shtns_vr, this%shtns_vt, this%shtns_vp)
#endif
            this%dLhw(1)=zero
            this%vhG(1) =zero
            this%vhC(1) =zero
            do lm=2,lm_max
               this%dLhw(lm)=dLh(lm)*w_Rloc(lm,nR)
               this%vhG(lm) =dw_Rloc(lm,nR) - &
                    cmplx(-aimag(z_Rloc(lm,nR)),real(z_Rloc(lm,nR)),kind=cp)
               this%vhC(lm) =dw_Rloc(lm,nR) + &
                    cmplx(-aimag(z_Rloc(lm,nR)),real(z_Rloc(lm,nR)),kind=cp)
            end do
         end if

         if ( lDeriv ) then
#ifdef WITH_SHTNS
           call pol_to_curlr_spat(z_Rloc(:, nR), this%shtns_cvr)
           call torpol_to_spat(dw_Rloc(:, nR), ddw_Rloc(:, nR), dz_Rloc(:, nR), nR, &
               this%shtns_dvrdr, this%shtns_dvtdr, this%shtns_dvpdr)
           call pol_to_grad_spat(w_Rloc(:, nR), &
               this%shtns_dvrdt, this%shtns_dvrdp)
#else
            this%dLhdw(1) =zero
            this%dLhz(1)  =zero
            this%dvhdrG(1)=zero
            this%dvhdrC(1)=zero
            do lm=2,lm_max
               this%dLhz(lm)  =dLh(lm)*z_Rloc(lm,nR)
               this%dLhdw(lm) =dLh(lm)*dw_Rloc(lm,nR)
               this%dvhdrG(lm)=ddw_Rloc(lm,nR) - &
                    cmplx(-aimag(dz_Rloc(lm,nR)),real(dz_Rloc(lm,nR)),kind=cp)
               this%dvhdrC(lm)=ddw_Rloc(lm,nR) + &
                    cmplx(-aimag(dz_Rloc(lm,nR)),real(dz_Rloc(lm,nR)),kind=cp)
            end do
#endif
         end if

      end if

      if ( l_mag .or. l_mag_LF ) then

         !PRINT*,"aj: ",SUM(ABS(aj(:,nR))),SUM(ABS(dLh))
         !PRINT*,"dj: ",SUM(ABS(dj(:,nR)))
#ifdef WITH_SHTNS
         call torpol_to_spat(b_Rloc(:, nR), db_Rloc(:, nR),  aj_Rloc(:, nR), nR, &
             this%shtns_br, this%shtns_bt, this%shtns_bp)
#else
         this%dLhb(1)=zero
         this%bhG(1) =zero
         this%bhC(1) =zero
         do lm=2,lm_max
            this%dLhb(lm)=dLh(lm)*b_Rloc(lm,nR)
            this%bhG(lm) =db_Rloc(lm,nR) - &
                 cmplx(-aimag(aj_Rloc(lm,nR)),real(aj_Rloc(lm,nR)),kind=cp)
            this%bhC(lm) =db_Rloc(lm,nR) + &
                 cmplx(-aimag(aj_Rloc(lm,nR)),real(aj_Rloc(lm,nR)),kind=cp)
         end do
#endif
#ifndef WITH_SHTNS
         if ( lGrenoble ) then ! Add dipole imposed by inner core
            lm=lm2(1,0)
            this%dLhb(lm)=this%dLhb(lm)+dLh(lm)*b0(nR)
            this%bhG(lm) =this%bhG(lm)+db0(nR)
            this%bhC(lm) =this%bhC(lm)+db0(nR)
         end if
#endif
         if ( lDeriv ) then
#ifdef WITH_SHTNS
            call torpol_to_curl_spat(b_Rloc(:, nR), db_Rloc(:, nR),     &
                ddb_Rloc(:, nR), aj_Rloc(:, nR), dj_Rloc(:, nR), nR,    &
                this%shtns_cbr, this%shtns_cbt, this%shtns_cbp)
#else
            this%dLhj(1)=zero
            this%cbhG(1)=zero
            this%cbhC(1)=zero
            do lm=2,lm_max
               this%dLhj(lm)=dLh(lm)*aj_Rloc(lm,nR)
               dbd     =or2(nR)*this%dLhb(lm)-ddb_Rloc(lm,nR)
               this%cbhG(lm)=dj_Rloc(lm,nR)-cmplx(-aimag(dbd),real(dbd),kind=cp)
               this%cbhC(lm)=dj_Rloc(lm,nR)+cmplx(-aimag(dbd),real(dbd),kind=cp)
            end do
            if ( lGrenoble ) then ! Add dipole imposed by inner core
               lm=lm2(1,0)
               this%cbhG(lm)=this%cbhG(lm)+cmplx(0.0_cp,ddb0(nR),kind=cp)
               this%cbhC(lm)=this%cbhC(lm)-cmplx(0.0_cp,ddb0(nR),kind=cp)
            end if
#endif
         end if

      end if   ! magnetic terms required ?

   end subroutine legPrepG
!------------------------------------------------------------------------------
   subroutine legPrep(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc, &
                      r,lDeriv,lHor,dLhw,vhG,vhC,dLhz,cvhG,&
                      cvhC)
      !
      !  Purpose of this subroutine is to prepare Legendre transforms     
      !  from (r,l,m) space to (r,theta,m) space by calculating           
      !  auxiliary arrays w, dw, ddw, ....... which contain               
      !  combinations of harmonic coeffs multiplied with (l,m)-dependend  
      !  factors as well as possible the radial dependencies.
      !
      !  lDeriv=.true. field derivatives required  for curl of field     
      !
    
      !-- Input variables:
      integer,     intent(in) :: lm_max
      complex(cp), intent(in) :: w(lm_max)
      complex(cp), intent(in) :: dw(lm_max)
      complex(cp), intent(in) :: ddw(lm_max)
      complex(cp), intent(in) :: z(lm_max)
      complex(cp), intent(in) :: dz(lm_max)
      real(cp),    intent(in) :: dLh(lm_max)
      integer,     intent(in) :: l_max,minc
      real(cp),    intent(in) :: r
      logical,     intent(in) :: lHor
      logical,     intent(in) :: lDeriv

      !-- Output variable:
      complex(cp), intent(out) :: dLhw(*),dLhz(*)
      complex(cp), intent(out) :: vhG(*),vhC(*)
      complex(cp), intent(out) :: cvhG(*),cvhC(*)

      !-- Local variables:
      integer :: lm,l,m
      real(cp) :: Or_e2
      complex(cp) :: help


      lm=0
      do m=0,l_max,minc
         do l=m,l_max
            lm=lm+1
            dLhw(lm)=dLh(lm)*w(lm)
            if ( lHor ) then
               vhG(lm) =dw(lm)-cmplx(-aimag(z(lm)),real(z(lm)),kind=cp)
               vhC(lm) =dw(lm)+cmplx(-aimag(z(lm)),real(z(lm)),kind=cp)
            end if
         end do
      end do
      dLhw(1)=zero
      if ( lHor ) then
         vhG(1) =zero
         vhC(1) =zero
      end if

      if ( lDeriv ) then
         Or_e2=one/r**2
         lm=0
         do m=0,l_max,minc
            do l=m,l_max
               lm=lm+1
               dLhz(lm)=dLh(lm)*z(lm)
               help=dLh(lm)*Or_e2*w(lm)-ddw(lm)
               cvhG(lm)=dz(lm)-cmplx(-aimag(help),real(help),kind=cp)
               cvhC(lm)=dz(lm)+cmplx(-aimag(help),real(help),kind=cp)
            end do
         end do
         dLhz(1)=zero
         cvhG(1)=zero
         cvhc(1)=zero
      end if

   end subroutine legPrep
!------------------------------------------------------------------------------
   subroutine legPrep_IC(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc, &
                         r,r_ICB,lDeriv,lHor,lCondIC,dLhw,vhG,&
                         vhC,dLhz,cvhG,cvhC)
      !
      !  Purpose of this subroutine is to prepare Legendre transforms     
      !  from (r,l,m) space to (r,theta,m) space by calculating           
      !  auxiliary arrays dLhw,vhG, ......... which contain               
      !  combinations of harmonic coeffs multiplied with (l,m)-dependend  
      !  factors as well as possible the radial dependencies.             
      !
      !  lDeriv=.true. field derivatives required  for curl of field     
      !
      !  Note: This routine is used for the inner core magnetic field     
      !  which has a special radial function ansatz. It can also be       
      !  used to prepare the calculation of a field in an insulating      
      !  inner core for lCondIC=.false.. For this the w has to be the     
      !  outer core poloidal field and nR is the grid point for the ICB.  
      !  In any case legTF can be used for the following Legendre         
      !  transform and fftJW for the Fourier transform.                   
      !

      !-- Input variables:
      integer,     intent(in) :: lm_max
      complex(cp), intent(in) :: w(lm_max)
      complex(cp), intent(in) :: dw(lm_max)
      complex(cp), intent(in) :: ddw(lm_max)
      complex(cp), intent(in) :: z(lm_max)
      complex(cp), intent(in) :: dz(lm_max)
      real(cp),    intent(in) :: dLh(lm_max)
      integer ,    intent(in) :: l_max,minc
      real(cp),    intent(in) :: r,r_ICB
      logical,     intent(in) :: lHor
      logical,     intent(in) :: lDeriv
      logical,     intent(in) :: lCondIC

      !-- Output variables:
      complex(cp), intent(out) :: dLhw(lm_max),dLhz(lm_max)
      complex(cp), intent(out) :: vhG(lm_max),vhC(lm_max)
      complex(cp), intent(out) :: cvhG(lm_max),cvhC(lm_max)

      !-- Local variables:
      integer :: lm,l,m
      complex(cp) :: help1,help2

      real(cp) :: rRatio,rDep(0:l_max),rDep2(0:l_max)


      rRatio  =r/r_ICB
      rDep(0) =rRatio
      rDep2(0)=one/r_ICB ! rDep2=rDep/r
      do l=1,l_max
         rDep(l) =rDep(l-1)*rRatio
         rDep2(l)=rDep2(l-1)*rRatio
      end do

      lm=0
      do m=0,l_max,minc
         do l=m,l_max
            lm=lm+1
            dLhw(lm)=rDep(l)*dLh(lm)*w(lm)
            if ( lHor ) then
               if ( lCondIC ) then
                  help1=rDep2(l)*((l+1)*w(lm)+r*dw(lm))
                  help2=rDep(l)*z(lm)
                  vhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=cp )
                  vhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=cp )
               else
                  vhG(lm)=rDep2(l)*(l+1)*w(lm) ! Only poloidal
                  vhC(lm)=vhG(lm)
               end if
            end if
         end do
      end do
      dLhw(1)=zero
      if ( lHor ) then
         vhG(1) =zero
         vhC(1) =zero
      end if

      if ( lDeriv ) then
         lm=0
         do m=0,l_max,minc
            do l=m,l_max
               lm=lm+1
               if ( lCondIC ) then
                  dLhz(lm)=rDep(l)*dLh(lm)*z(lm)
                  help1=rDep(l)*( (l+1)*z(lm)/r+dz(lm) )
                  help2=rDep(l)*(-two*(l+1)/r*dw(lm)-ddw(lm))
                  cvhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=cp)
                  cvhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=cp)
               else
                  dLhz(lm)=zero
                  cvhG(lm)=zero
                  cvhC(lm)=zero
               end if
            end do
         end do
         dLhz(1)=zero
         cvhG(1)=zero
         cvhc(1)=zero
      end if

   end subroutine legPrep_IC
!------------------------------------------------------------------------------
end module leg_helper_mod
