module leg_helper_mod

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
   use const, only: zero

   implicit none

   private
   
   type, public :: leg_helper_t
      !----- Help arrays for Legendre transform calculated in legPrepG:
      !      Parallelizatio note: these are the R-distributed versions
      !      of the field scalars.
      complex(kind=8), allocatable :: dLhw(:), dLhdw(:), dLhz(:), dLhb(:), dLhj(:)
      complex(kind=8), allocatable :: vhG(:), vhC(:), dvhdrG(:), dvhdrC(:)
      complex(kind=8), allocatable :: bhG(:), bhC(:), cbhG(:), cbhC(:)
      !----- R-distributed versions of scalar fields (see c_fields.f):
      complex(kind=8), allocatable :: sR(:), dsR(:), preR(:), dpR(:)
      real(kind=8), allocatable :: zAS(:), dzAS(:), ddzAS(:) ! used in TO
      real(kind=8) :: omegaIC,omegaMA
      complex(kind=8), allocatable :: bCMB(:)
 
   contains
 
      procedure :: initialize
      procedure :: legPrepG
 
   end type leg_helper_t

   public :: legPrep, legPrep_IC

contains

   subroutine initialize(this,lm_max,lm_maxMag,l_max)

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
      allocate( this%zAS(l_max+1),this%dzAS(l_max+1),&
           & this%ddzAS(l_max+1) ) ! used in TO

      allocate( this%bCMB(lm_maxMag) )

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine legPrepG(this,nR,nBc,lDeriv,lRmsCalc,l_frame, &
        &              lTOnext,lTOnext2,lTOcalc)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to prepare Legendre transforms     |
      !  |  from (r,l,m) space to (r,theta,m) space by calculating           |
      !  |  auxiliary arrays dpdw,dpddw, ....... dLhj which contain          |
      !  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
      !  |  factors as well as the radial dependence.                        |
      !  |    nBc  =0 standard inner radial grid point                       |
      !  |    nBc  =1 radial velocity zero, spatial derivs not needed        |
      !  |    nBc  =2 all velocity comp. zero, spatial derivs not needed     |
      !  |   lDeriv=.true. field derivatives required                        |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      class(leg_helper_t) :: this
    
      !-- Input of variables
      integer, intent(in) :: nR              ! radial level
      integer, intent(in) :: nBc             ! boundary condition
      logical, intent(in) :: lDeriv          ! get also field derivatives !
      logical, intent(in) :: lRmsCalc        ! Rms force balance ?
      logical, intent(in) :: l_frame         ! Calculate movie frame?
      logical, intent(in) :: lTOnext         ! for TO output
      logical, intent(in) :: lTOnext2
      logical, intent(in) :: lTOcalc
    
      !-- Input of scalar fields in LM-distributed space
      !   These have to be collected and stored in the
      !   R-distributed output variables
    
      !-- Local variables:
      integer :: lm,l,m
      complex(kind=8) :: dbd
    
    
      if ( nR == n_r_icb ) this%omegaIC=omega_ic
      if ( nR == n_r_cmb ) this%omegaMA=omega_ma
      if ( l_conv .or. l_mag_kin ) then
    
         if ( l_heat ) then
            do lm=1,lm_max
               this%sR(lm) =s_Rloc(lm,nR)   ! used for plotting and Rms
               this%dsR(lm)=ds_Rloc(lm,nR)  ! used for plotting and Rms
            end do
         end if
         if ( lTOnext .or. lTOnext2 .or. lTOCalc ) then
            do lm=1,lm_max
               l=lm2l(lm)
               m=lm2m(lm)
               if ( l <= l_max .and. m == 0 ) then
                  this%zAS(l+1)  =real(z_Rloc(lm,nR))   ! used in TO
                  this%dzAS(l+1) =real(dz_Rloc(lm,nR))  ! used in TO (anelastic)
                  this%ddzAS(l+1)=ddzASL(l+1,nR)    ! used in TO
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
            this%dLhw(1)=zero
            this%vhG(1) =zero
            this%vhC(1) =zero
            do lm=2,lm_max
               this%dLhw(lm)=dLh(lm)*w_Rloc(lm,nR)
               this%vhG(lm) =dw_Rloc(lm,nR) - &
                    cmplx(-aimag(z_Rloc(lm,nR)),real(z_Rloc(lm,nR)),kind=kind(0d0))
               this%vhC(lm) =dw_Rloc(lm,nR) + &
                    cmplx(-aimag(z_Rloc(lm,nR)),real(z_Rloc(lm,nR)),kind=kind(0d0))
            end do
         end if
    
         if ( lDeriv ) then
            this%dLhdw(1) =zero
            this%dLhz(1)  =zero
            this%dvhdrG(1)=zero
            this%dvhdrC(1)=zero
            do lm=2,lm_max
               this%dLhz(lm)  =dLh(lm)*z_Rloc(lm,nR)
               this%dLhdw(lm) =dLh(lm)*dw_Rloc(lm,nR)
               this%dvhdrG(lm)=ddw_Rloc(lm,nR) - &
                    cmplx(-aimag(dz_Rloc(lm,nR)),real(dz_Rloc(lm,nR)),kind=kind(0d0))
               this%dvhdrC(lm)=ddw_Rloc(lm,nR) + &
                    cmplx(-aimag(dz_Rloc(lm,nR)),real(dz_Rloc(lm,nR)),kind=kind(0d0))
            end do
         end if
    
      end if
    
      if ( l_mag .or. l_mag_LF ) then
    
         !PRINT*,"aj: ",SUM(ABS(aj(:,nR))),SUM(ABS(dLh))
         !PRINT*,"dj: ",SUM(ABS(dj(:,nR)))
         this%dLhb(1)=zero
         this%bhG(1) =zero
         this%bhC(1) =zero
         do lm=2,lm_max
            this%dLhb(lm)=dLh(lm)*b_Rloc(lm,nR)
            this%bhG(lm) =db_Rloc(lm,nR) - &
                 cmplx(-aimag(aj_Rloc(lm,nR)),real(aj_Rloc(lm,nR)),kind=kind(0d0))
            this%bhC(lm) =db_Rloc(lm,nR) + &
                 cmplx(-aimag(aj_Rloc(lm,nR)),real(aj_Rloc(lm,nR)),kind=kind(0d0))
         end do
         if ( lGrenoble ) then ! Add dipole imposed by inner core
            lm=lm2(1,0)
            this%dLhb(lm)=this%dLhb(lm)+dLh(lm)*b0(nR)
            this%bhG(lm) =this%bhG(lm)+db0(nR)
            this%bhC(lm) =this%bhC(lm)+db0(nR)
         end if
         if ( lDeriv ) then
            this%dLhj(1)=zero
            this%cbhG(1)=zero
            this%cbhC(1)=zero
            do lm=2,lm_max
               this%dLhj(lm)=dLh(lm)*aj_Rloc(lm,nR)
               dbd     =or2(nR)*this%dLhb(lm)-ddb_Rloc(lm,nR)
               this%cbhG(lm)=dj_Rloc(lm,nR)-cmplx(-aimag(dbd),real(dbd), &
                                                        kind=kind(0d0))
               this%cbhC(lm)=dj_Rloc(lm,nR)+cmplx(-aimag(dbd),real(dbd), &
                                                        kind=kind(0d0))
            end do
            if ( lGrenoble ) then ! Add dipole imposed by inner core
               lm=lm2(1,0)
               this%cbhG(lm)=this%cbhG(lm)+cmplx(0.D0,ddb0(nR), &
                                                             kind=kind(0d0))
               this%cbhC(lm)=this%cbhC(lm)-cmplx(0.D0,ddb0(nR), &
                                                             kind=kind(0d0))
            end if
         end if
    
      end if   ! magnetic terms required ?
    
   end subroutine legPrepG
!------------------------------------------------------------------------------
   subroutine legPrep(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc, &
                      r,lDeriv,lHor,dLhw,vhG,vhC,dLhz,cvhG,&
                      cvhC)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to prepare Legendre transforms     |
      !  |  from (r,l,m) space to (r,theta,m) space by calculating           |
      !  |  auxiliary arrays w, dw, ddw, ....... which contain               |
      !  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
      !  |  factors as well as possible the radial dependencies.             |
      !  |   lHor=.true. horizontal components required                      |
      !  |   lDeriv=.true. field derivatives required  for curl of field     |
      !  +-------------------------------------------------------------------+
    
      implicit none
    
      !-- Input variables:
      integer,         intent(in) :: lm_max
      complex(kind=8), intent(in) :: w(lm_max)
      complex(kind=8), intent(in) :: dw(lm_max)
      complex(kind=8), intent(in) :: ddw(lm_max)
      complex(kind=8), intent(in) :: z(lm_max)
      complex(kind=8), intent(in) :: dz(lm_max)
      real(kind=8),    intent(in) :: dLh(lm_max)
      integer,         intent(in) :: l_max,minc
      real(kind=8),    intent(in) :: r
      logical,         intent(in) :: lHor
      logical,         intent(in) :: lDeriv
    
      !-- Output variable:
      complex(kind=8), intent(out) :: dLhw(*),dLhz(*)
      complex(kind=8), intent(out) :: vhG(*),vhC(*)
      complex(kind=8), intent(out) :: cvhG(*),cvhC(*)
    
      !-- Local variables:
      integer :: lm,l,m
      real(kind=8) :: Or_e2
      complex(kind=8) :: help
    
    
      lm=0
      do m=0,l_max,minc
         do l=m,l_max
            lm=lm+1
            dLhw(lm)=dLh(lm)*w(lm)
            if ( lHor ) then
               vhG(lm) =dw(lm)-cmplx(-aimag(z(lm)),real(z(lm)),kind=kind(0d0))
               vhC(lm) =dw(lm)+cmplx(-aimag(z(lm)),real(z(lm)),kind=kind(0d0))
            end if
         end do
      end do
      dLhw(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
      if ( lHor ) then
         vhG(1) =cmplx(0.D0,0.D0,kind=kind(0d0))
         vhC(1) =cmplx(0.D0,0.D0,kind=kind(0d0))
      end if
    
      if ( lDeriv ) then
         Or_e2=1.D0/r**2
         lm=0
         do m=0,l_max,minc
            do l=m,l_max
               lm=lm+1
               dLhz(lm)=dLh(lm)*z(lm)
               help=dLh(lm)*Or_e2*w(lm)-ddw(lm)
               cvhG(lm)=dz(lm)-cmplx(-aimag(help),real(help),kind=kind(0d0))
               cvhC(lm)=dz(lm)+cmplx(-aimag(help),real(help),kind=kind(0d0))
            end do
         end do
         dLhz(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
         cvhG(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
         cvhc(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
      end if
    
   end subroutine legPrep
!------------------------------------------------------------------------------
   subroutine legPrep_IC(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc, &
                         r,r_ICB,lDeriv,lHor,lCondIC,dLhw,vhG,&
                         vhC,dLhz,cvhG,cvhC)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to prepare Legendre transforms     |
      !  |  from (r,l,m) space to (r,theta,m) space by calculating           |
      !  |  auxiliary arrays dLhw,vhG, ......... which contain               |
      !  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
      !  |  factors as well as possible the radial dependencies.             |
      !  |   lHor=.true. horizontal components required                      |
      !  |   lDeriv=.true. field derivatives required  for curl of field     |
      !  |  Note: This routine is used for the inner core magnetic field     |
      !  |  which has a special radial function ansatz. It can also be       |
      !  |  used to prepare the calculation of a field in an insulating      |
      !  |  inner core for lCondIC=.false.. For this the w has to be the     |
      !  |  outer core poloidal field and nR is the grid point for the ICB.  |
      !  |  In any case legTF can be used for the following Legendre         |
      !  |  transform and fftJW for the Fourier transform.                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,         intent(in) :: lm_max
      complex(kind=8), intent(in) :: w(lm_max)
      complex(kind=8), intent(in) :: dw(lm_max)
      complex(kind=8), intent(in) :: ddw(lm_max)
      complex(kind=8), intent(in) :: z(lm_max)
      complex(kind=8), intent(in) :: dz(lm_max)
      real(kind=8),    intent(in) :: dLh(lm_max)
      integer ,        intent(in) :: l_max,minc
      real(kind=8),    intent(in) :: r,r_ICB
      logical,         intent(in) :: lHor
      logical,         intent(in) :: lDeriv
      logical,         intent(in) :: lCondIC

      !-- Output variables:
      complex(kind=8), intent(out) :: dLhw(lm_max),dLhz(lm_max)
      complex(kind=8), intent(out) :: vhG(lm_max),vhC(lm_max)
      complex(kind=8), intent(out) :: cvhG(lm_max),cvhC(lm_max)
       
      !-- Local variables:
      integer :: lm,l,m
      complex(kind=8) :: help1,help2
              
      real(kind=8) :: rRatio,rDep(0:l_max-1),rDep2(0:l_max-1)


      rRatio  =r/r_ICB
      rDep(0) =rRatio
      rDep2(0)=1.D0/r_ICB ! rDep2=rDep/r
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
                  vhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=kind(0d0) )
                  vhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=kind(0d0) )
               else
                  vhG(lm)=rDep2(l)*(l+1)*w(lm) ! Only poloidal
                  vhC(lm)=vhG(lm)
               end if
            end if
         end do
      end do
      dLhw(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
      if ( lHor ) then
         vhG(1) =cmplx(0.D0,0.D0,kind=kind(0d0))
         vhC(1) =cmplx(0.D0,0.D0,kind=kind(0d0))
      end if

      if ( lDeriv ) then
         lm=0
         do m=0,l_max,minc
            do l=m,l_max
               lm=lm+1
               if ( lCondIC ) then
                  dLhz(lm)=rDep(l)*dLh(lm)*z(lm)
                  help1=rDep(l)*( (l+1)*z(lm)/r+dz(lm) )
                  help2=rDep(l)*(-2.D0*(l+1)/r*dw(lm)-ddw(lm))
                  cvhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=kind(0d0))
                  cvhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=kind(0d0))
               else
                  dLhz(lm)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  cvhG(lm)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  cvhC(lm)=cmplx(0.D0,0.D0,kind=kind(0d0))
               end if
            end do
         end do
         dLhz(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
         cvhG(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
         cvhc(1)=cmplx(0.D0,0.D0,kind=kind(0d0))
      end if

   end subroutine legPrep_IC
!------------------------------------------------------------------------------
end module leg_helper_mod
