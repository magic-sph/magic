!$Id$
module radial_functions

   use truncation, only: n_r_max, n_cheb_max, n_r_ic_max
   use matrices, only: s0Mat,s0Pivot
#ifdef WITH_MKL_LU
   use lapack95, only: getrf,getrs
#else
   use algebra, only: sgesl,sgefa
#endif
   use const, only: sq4pi, one, two, three, four, half
   use physical_parameters
   use num_param, only: alpha
   use logic, only: l_mag, l_cond_ic, l_heat, l_anelastic_liquid, &
                    l_isothermal, l_anel, l_newmap
   use init_costf ! Everything is needed
   use chebyshev_polynoms_mod ! Everything is needed
   use cosine_transform, only: costf1
   use radial_der, only: get_dr
 
   implicit none

   private
 
   !-- arrays depending on r:
   real(cp), public, allocatable :: r(:)
   real(cp), public, allocatable :: r_ic(:)
   real(cp), public, allocatable :: O_r_ic(:)
   real(cp), public, allocatable :: O_r_ic2(:)
   real(cp), public, allocatable :: or1(:),or2(:),or3(:),or4(:)
   real(cp), public, allocatable :: otemp1(:),rho0(:),temp0(:)
   real(cp), public, allocatable :: dtemp0(:),d2temp0(:)
   real(cp), public, allocatable :: dentropy0(:)
   real(cp), public, allocatable :: orho1(:),orho2(:)
   real(cp), public, allocatable :: beta(:), dbeta(:)
   real(cp), public, allocatable :: drx(:),ddrx(:),dddrx(:)
   real(cp), public :: dr_fac,dr_fac_ic,alpha1,alpha2
   real(cp), public :: topcond, botcond
   real(cp), public :: r_cmb,r_icb,r_surface
 
   !-- arrays for buoyancy, depend on Ra and Pr:
   real(cp), public, allocatable :: rgrav(:),agrav(:)
 
   !-- chebychev polynomials, derivatives and integral:
   real(cp), public :: cheb_norm                   ! cheb normalisation 
   real(cp), public, allocatable :: cheb(:,:)       ! Chebychev polynomials
   real(cp), public, allocatable :: dcheb(:,:)      ! first radial derivative
   real(cp), public, allocatable :: d2cheb(:,:)     ! second radial derivative
   real(cp), public, allocatable :: d3cheb(:,:)     ! third radial derivative
   real(cp), public, allocatable :: cheb_int(:)     ! array for cheb integrals !
   integer, public :: nDi_costf1
   integer, public, allocatable :: i_costf_init(:)   ! info for transform
   integer, public :: nDd_costf1
   real(cp), public, allocatable ::  d_costf_init(:)   ! info for tranform
 
   !-- same for inner core:
   real(cp), public :: cheb_norm_ic
   real(cp), public, allocatable :: cheb_ic(:,:)   
   real(cp), public, allocatable :: dcheb_ic(:,:)  
   real(cp), public, allocatable :: d2cheb_ic(:,:) 
   real(cp), public, allocatable :: cheb_int_ic(:)        
   integer, public :: nDi_costf1_ic
 
   integer, public, allocatable :: i_costf1_ic_init(:)
   integer, public :: nDd_costf1_ic
 
   real(cp), public, allocatable ::  d_costf1_ic_init(:)
   integer, public :: nDi_costf2_ic
 
   integer, public, allocatable :: i_costf2_ic_init(:)
   integer, public :: nDd_costf2_ic
 
   real(cp), public, allocatable ::  d_costf2_ic_init(:)
 
   !-- Radius functions for cut-back grid without boundaries:
   !-- (and for the nonlinear mapping)
   integer, public, allocatable :: i_costf_initC(:)   ! info for transform
   real(cp), public, allocatable ::  d_costf_initC(:)   ! info for tranform
   real(cp), public, allocatable ::  rC(:),dr_facC(:)
   real(cp), public :: alph1,alph2
   integer, public :: n_r_maxC,n_cheb_maxC,nCut
 
   real(cp), public, allocatable :: lambda(:),dLlambda(:),jVarCon(:)
   real(cp), public, allocatable :: sigma(:)
   real(cp), public, allocatable :: kappa(:),dLkappa(:)
   real(cp), public, allocatable :: visc(:),dLvisc(:)
   real(cp), public, allocatable :: divKtemp0(:)
   real(cp), public, allocatable :: epscProf(:)

   public :: initialize_radial_functions, radial, transportProperties

contains

   subroutine initialize_radial_functions

      nDi_costf1=2*n_r_max+2
      nDd_costf1=2*n_r_max+5

      nDi_costf1_ic=2*n_r_ic_max+2
      nDd_costf1_ic=2*n_r_ic_max+5
      nDi_costf2_ic=2*n_r_ic_max
      nDd_costf2_ic=2*n_r_ic_max+n_r_ic_max/2+5

      ! allocate the arrays
      allocate( r(n_r_max) )
      allocate( r_ic(n_r_ic_max) )
      allocate( O_r_ic(n_r_ic_max) )
      allocate( O_r_ic2(n_r_ic_max) )
      allocate( or1(n_r_max),or2(n_r_max),or3(n_r_max),or4(n_r_max) )
      allocate( otemp1(n_r_max),rho0(n_r_max),temp0(n_r_max) )
      allocate( dtemp0(n_r_max),d2temp0(n_r_max),dentropy0(n_r_max) )
      allocate( orho1(n_r_max),orho2(n_r_max) )
      allocate( beta(n_r_max), dbeta(n_r_max) )
      allocate( drx(n_r_max),ddrx(n_r_max),dddrx(n_r_max) )
      allocate( rgrav(n_r_max),agrav(n_r_max) )

      allocate( cheb(n_r_max,n_r_max) )     ! Chebychev polynomials
      allocate( dcheb(n_r_max,n_r_max) )    ! first radial derivative
      allocate( d2cheb(n_r_max,n_r_max) )   ! second radial derivative
      allocate( d3cheb(n_r_max,n_r_max) )   ! third radial derivative
      allocate( cheb_int(n_r_max) )         ! array for cheb integrals !
      allocate( i_costf_init(nDi_costf1) )  ! info for transform
      allocate( d_costf_init(nDd_costf1) ) ! info for tranform

      allocate( cheb_ic(n_r_ic_max,n_r_ic_max) )
      allocate( dcheb_ic(n_r_ic_max,n_r_ic_max) )
      allocate( d2cheb_ic(n_r_ic_max,n_r_ic_max) )
      allocate( cheb_int_ic(n_r_ic_max) )
      allocate( i_costf1_ic_init(nDi_costf1_ic) )
      allocate( d_costf1_ic_init(nDd_costf1_ic) )
      allocate( i_costf2_ic_init(nDi_costf2_ic) )
      allocate( d_costf2_ic_init(nDd_costf2_ic) )

      allocate( lambda(n_r_max),dLlambda(n_r_max),jVarCon(n_r_max) )
      allocate( sigma(n_r_max) )
      allocate( kappa(n_r_max),dLkappa(n_r_max) )
      allocate( visc(n_r_max),dLvisc(n_r_max) )
      allocate( epscProf(n_r_max),divKtemp0(n_r_max) )

      allocate( i_costf_initC(nDi_costf1) )   ! info for transform
      allocate( d_costf_initC(nDd_costf1) )   ! info for tranform
      allocate( rC(n_r_max),dr_facC(n_r_max) )

   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine radial
      !--------------------------------------------------------------------
      !  Calculates everything needed for radial functions, transfroms etc.
      !--------------------------------------------------------------------

      !-- Local variables:
      integer :: n_r,n_cheb,n_cheb_int
      integer :: n_r_ic_tot,k
      integer :: n_const(1)

      !integer :: n_r_start
      real(cp) :: fac_int
      real(cp) :: r_cheb(n_r_max)
      real(cp) :: r_cheb_ic(2*n_r_ic_max-1),r_ic_2(2*n_r_ic_max-1)
      real(cp) :: alphaT(n_r_max), drho0(n_r_max)
      real(cp) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

      real(cp) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9 ! polynomial fit for density
      real(cp) :: temptop,gravtop,rhotop
      real(cp) :: hcomp,CompNb,GrunNb
      real(cp) :: dtemp0cond(n_r_max),dtemp0ad(n_r_max),hcond(n_r_max)
      real(cp) :: func(n_r_max)

      real(cp) :: rrOcmb,rStrat
      real(cp) :: gravFit(n_r_max),rhoFit(n_r_max) ! normalised to rcmb
      real(cp) :: w1(n_r_max),w2(n_r_max)

      integer :: stop_signal
      integer :: filehandle,nCheb

      !-- Radial grid point:
      !   radratio is aspect ratio
      !   radratio = (inner core r) / (CMB r) = r_icb/r_cmb
      r_cmb=one/(one-radratio)
      r_icb=r_cmb-one
      r_surface=2.8209_cp    ! in units of (r_cmb-r_icb)

      cheb_norm=sqrt(two/real(n_r_max-1,kind=cp))
      dr_fac=two/(r_cmb-r_icb)

      if ( l_newmap ) then
         alpha1         =alph1
         alpha2         =alph2
         paraK=atan(alpha1*(1+alpha2))/atan(alpha1*(1-alpha2))
         paraX0=(paraK-1)/(paraK+1)
         lambd=atan(alpha1*(1-alpha2))/(1-paraX0)
      else
         alpha1         =0.0_cp
         alpha2         =0.0_cp
      end if

      !-- Start with outer core:
      !   cheb_grid calculates the n_r_max gridpoints, these
      !   are the extrema of a Cheb pylonomial of degree n_r_max-1,
      !   r_cheb are the grid_points in the Cheb intervall [-1,1]
      !   and r are these points mapped to the intervall [r_icb,r_cmb]:
      call cheb_grid(r_icb,r_cmb,n_r_max-1,r,r_cheb, &
                           alpha1,alpha2,paraX0,lambd)
#if 0
      do n_r=1,n_r_max
         write(*,"(I3,2ES20.12)") n_r,r_cheb(n_r),r(n_r)
      end do
#endif

      if ( l_newmap ) then
         do n_r=1,n_r_max
            drx(n_r) =                          (2*alpha1) /      &
                 ((1+alpha1**2*(2*r(n_r)-r_icb-r_cmb-alpha2)**2)* &
                 lambd)
            ddrx(n_r) = -(8*alpha1**3*(2*r(n_r)-r_icb-r_cmb-alpha2)) / &
                 ((1+alpha1**2*(-2*r(n_r)+r_icb+r_cmb+alpha2)**2)**2*  & 
                 lambd)
            dddrx(n_r) =               (16*alpha1**3*(-1+3*alpha1**2* &
                                (-2*r(n_r)+r_icb+r_cmb+alpha2)**2)) / &
                 ((1+alpha1**2*(-2*r(n_r)+r_icb+r_cmb+alpha2)**2)**3* &
                 lambd)
         end do
      else
         do n_r=1,n_r_max
            drx(n_r)=two/(r_cmb-r_icb)
            ddrx(n_r)=0
            dddrx(n_r)=0
         end do
      end if

      !-- Calculate chebs and their derivatives up to degree n_r_max-1
      !   on the n_r radial grid points r:
      call get_chebs(n_r_max,r_icb,r_cmb,r_cheb,n_r_max,       &
                     cheb,dcheb,d2cheb,d3cheb,n_r_max,n_r_max, &
                     drx,ddrx,dddrx)

#if 0
      open(NEWUNIT=filehandle,file="r_cheb.dat")
      do n_r=1,n_r_max
         write(filehandle,"(2ES20.12)") r_cheb(n_r),r(n_r)
      end do
      close(filehandle)
#endif

      !-- Initialize fast cos transform for chebs:
      call init_costf1(n_r_max,i_costf_init,nDi_costf1, &
                       d_costf_init,nDd_costf1)

      or1=one/r        ! 1/r
      or2=or1*or1       ! 1/r**2
      or3=or1*or2       ! 1/r**3
      or4=or2*or2       ! 1/r**4

      !-- Fit to an interior model
      if ( index(interior_model,'JUP') /= 0 ) then
         a0=-122.36071577_cp
         a1= 440.86067831_cp
         a2=-644.82401602_cp
         a3= 491.00495215_cp
         a4=-201.96655354_cp
         a5=  37.38863965_cp
         a6=  -4.60312999_cp
         a7=   4.46020423_cp

         do n_r=1,n_r_max
            rrOcmb = r(n_r)/r_cmb*r_cut_model
            rhoFit(n_r) = a7 + a6*rrOcmb   + a5*rrOcmb**2 + a4*rrOcmb**3 &
                             + a3*rrOcmb**4+ a2*rrOcmb**5 + a1*rrOcmb**6 &
                             + a0*rrOcmb**7
            gravFit(n_r)=four*rrOcmb - three*rrOcmb**2
            temp0(n_r)=rhoFit(n_r)**(one/polind)
         end do

         ! To normalise to the outer radius
         temptop=temp0(1)
         rhotop=rhoFit(1)
         gravtop=gravFit(1)

         temp0  =temp0/temptop
         gravFit=gravFit/gravtop

         ! Derivative of the temperature needed to get alpha_T
         call get_dr(temp0,dtemp0,n_r_max,n_cheb_max,w1, &
                     w2,i_costf_init,d_costf_init,drx)

         alphaT=-dtemp0/(gravFit*temp0)

         ! Dissipation number needed in the dissipation numbers
         DissNb=alphaT(1)
         ViscHeatFac=DissNb*pr/raScaled
         if (l_mag) then
            OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
         end if
         ! Adiabatic: buoyancy term is linked to the temperature gradient

         !       dT
         !      ---- =  -Di * alpha_T * T * grav
         !       dr
         rgrav=-BuoFac*dtemp0/DissNb
         rho0=rhoFit/rhotop

         call get_dr(rho0,drho0,n_r_max,n_cheb_max,w1, &
                     w2,i_costf_init,d_costf_init,drx)
         beta=drho0/rho0
         call get_dr(beta,dbeta,n_r_max,n_cheb_max,w1,     &
                     w2,i_costf_init,d_costf_init,drx)
         call get_dr(dtemp0,d2temp0,n_r_max,n_cheb_max,w1, &
                  w2,i_costf_init,d_costf_init,drx)
         dentropy0=0.0_cp
         
      else if ( index(interior_model,'SAT') /= 0 ) then

         ! the shell can't be thicker than eta=0.15, because the fit doesn't work
         ! below that (in Nadine's profile, that's where the IC is, anyway)
         a0=  7791.6205_cp
         a1=-38964.7491_cp
         a2= 82576.2667_cp
         a3=-96511.4441_cp
         a4= 67847.2219_cp
         a5=-29393.1585_cp
         a6=  7745.12023_cp
         a7= -1177.98473_cp
         a8=    86.0013409_cp
         a9=     1.11379407_cp

         do n_r=1,n_r_max
            rrOcmb = r(n_r)/r_cmb*r_cut_model
            rhoFit(n_r) = a9 + a8*rrOcmb   + a7*rrOcmb**2 + a6*rrOcmb**3 &
                             + a5*rrOcmb**4+ a4*rrOcmb**5 + a3*rrOcmb**6 &
                             + a2*rrOcmb**7+ a1*rrOcmb**8 + a0*rrOcmb**9
            gravFit(n_r)=four*rrOcmb - three*rrOcmb**2
            temp0(n_r)=rhoFit(n_r)**(one/polind)
         end do

         ! To normalise to the outer radius
         temptop=temp0(1)
         rhotop=rhoFit(1)
         gravtop=gravFit(1)

         temp0  =temp0/temptop
         gravFit=gravFit/gravtop

         ! Derivative of the temperature needed to get alpha_T
         call get_dr(temp0,dtemp0,n_r_max,n_cheb_max,w1, &
                     w2,i_costf_init,d_costf_init,drx)

         alphaT=-dtemp0/(gravFit*temp0)

         ! Inverse of the Froude number needed in the dissipation numbers
         DissNb=alphaT(1)
         ViscHeatFac=DissNb*pr/raScaled
         if (l_mag) then
            OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
         end if
         ! Adiabatic: buoyancy term is linked to the temperature gradient

         !       dT
         !      ---- =  Di * alpha_T * T * grav
         !       dr

         ! N.B. rgrav is not gravity but the whole RHS !!!
         rgrav=-BuoFac*dtemp0/DissNb
         rho0=rhoFit/rhotop

         call get_dr(rho0,drho0,n_r_max,n_cheb_max,w1, &
                     w2,i_costf_init,d_costf_init,drx)
         beta=drho0/rho0
         call get_dr(beta,dbeta,n_r_max,n_cheb_max,w1,     &
                     w2,i_costf_init,d_costf_init,drx)
         call get_dr(dtemp0,d2temp0,n_r_max,n_cheb_max,w1, &
                  w2,i_costf_init,d_costf_init,drx)
         dentropy0=0.0_cp

      else if ( index(interior_model,'SUN') /= 0 ) then

         a7=  113.63001006_cp
         a6= -691.6084317_cp
         a5= 1615.06990369_cp
         a4=-1570.0073169_cp
         a3=  -24.81006594_cp
         a2= 1336.03589943_cp
         a1=-1038.72509351_cp
         a0=  260.41507794_cp

         do n_r=1,n_r_max
            rrOcmb = r(n_r)/r_cmb*r_cut_model
            rhoFit(n_r) = a7 + a6*rrOcmb   + a5*rrOcmb**2 + a4*rrOcmb**3 &
                             + a3*rrOcmb**4+ a2*rrOcmb**5 + a1*rrOcmb**6 &
                             + a0*rrOcmb**7
            gravFit(n_r)=four*rrOcmb - three*rrOcmb**2
            temp0(n_r)=rhoFit(n_r)**(one/polind)
         end do

         ! To normalise to the outer radius
         temptop=temp0(1)
         rhotop=rhoFit(1)
         gravtop=gravFit(1)

         temp0  =temp0/temptop
         gravFit=gravFit/gravtop

         ! Derivative of the temperature needed to get alpha_T
         call get_dr(temp0,dtemp0,n_r_max,n_cheb_max,w1, &
                     w2,i_costf_init,d_costf_init,drx)

         alphaT=-dtemp0/(gravFit*temp0)

         ! Dissipation number needed in the dissipation numbers
         DissNb=alphaT(1)
         ViscHeatFac=DissNb*pr/raScaled
         if (l_mag) then
            OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
         end if
         ! Adiabatic: buoyancy term is linked to the temperature gradient

         !       dT
         !      ---- =  -Di * alpha_T * T * grav
         !       dr
         rgrav=-BuoFac*dtemp0/DissNb
         rho0=rhoFit/rhotop

         call get_dr(rho0,drho0,n_r_max,n_cheb_max,w1, &
                     w2,i_costf_init,d_costf_init,drx)
         beta=drho0/rho0
         call get_dr(beta,dbeta,n_r_max,n_cheb_max,w1,     &
                     w2,i_costf_init,d_costf_init,drx)
         call get_dr(dtemp0,d2temp0,n_r_max,n_cheb_max,w1, &
                  w2,i_costf_init,d_costf_init,drx)
         dentropy0=0.0_cp

      else if ( index(interior_model,'EARTH') /= 0 ) then
         DissNb=0.3929_cp ! Di = \alpha_O g d / c_p
         CompNb=0.0566_cp ! Co = \alpha_O T_O
         GrunNb=1.5_cp ! Gruneisen paramater
         hcomp =2.2_cp*r_cmb

         alphaT=(one+0.6_cp*r**2/hcomp**2)/(one+0.6_cp/2.2_cp**2)
         rgrav =(r-0.6_cp*r**3/hcomp**2)/(r_cmb*(one-0.6_cp/2.2_cp**2))

         !dentropy0 = -half*(ampStrat+one)*(one-tanh(slopeStrat*(r-rStrat)))+ &
         !            & ampStrat

         !! d ln(temp0) / dr
         !dtemp0=epsS*dentropy0-DissNb*alphaT*rgrav

         !call getBackground(dtemp0,0.0_cp,temp0)
         !temp0=exp(temp0) ! this was ln(T_0)
         !dtemp0=dtemp0*temp0

         !drho0=-CompNb*epsS*alphaT*temp0*dentropy0-DissNb/GrunNb*alphaT*rgrav
         !call getBackground(drho0,0.0_cp,rho0)
         !rho0=exp(rho0) ! this was ln(rho_0)
         !beta=drho0

         hcond = (one-0.4469_cp*(r/r_cmb)**2)/(one-0.4469_cp)
         hcond = hcond/hcond(1)
         temp0 = (one+GrunNb*(r_icb**2-r**2)/hcomp**2)
         temp0 = temp0/temp0(1)
         dtemp0cond=-cmbHflux/(r**2*hcond)
          
         do k=1,10 ! 10 iterations is enough to converge
            dtemp0ad=-DissNb*alphaT*rgrav*temp0-epsS*temp0(n_r_max)
            n_const=minloc(abs(dtemp0ad-dtemp0cond))
            rStrat=r(n_const(1))
            func=half*(tanh(slopeStrat*(r-rStrat))+one)

            if ( rStrat<r_cmb ) then
               dtemp0=func*dtemp0cond+(one-func)*dtemp0ad
            else
               dtemp0=dtemp0ad
            end if

            call getBackground(dtemp0,one,temp0)
         end do

         dentropy0=dtemp0/temp0/epsS+DissNb*alphaT*rgrav/epsS
         drho0=-CompNb*epsS*alphaT*temp0*dentropy0-DissNb*alphaT*rgrav/GrunNb
         call getBackground(drho0,0.0_cp,rho0)
         rho0=exp(rho0) ! this was ln(rho_0)
         beta=drho0

         ! The final stuff is always required
         call get_dr(beta,dbeta,n_r_max,n_cheb_max,w1,     &
                &    w2,i_costf_init,d_costf_init,drx)
         call get_dr(dtemp0,d2temp0,n_r_max,n_cheb_max,w1, &
                &    w2,i_costf_init,d_costf_init,drx)

         ViscHeatFac=DissNb*pr/raScaled
         if (l_mag) then
            OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
         end if

         ! N.B. rgrav is not gravity but alpha * grav
         rgrav = BuoFac*alphaT*rgrav

      else  !-- Usual polytropic reference state
         ! g(r) = g0 + g1*r/ro + g2*(ro/r)**2
         ! Default values: g0=0, g1=1, g2=0
         ! An easy way to change gravity
         rgrav=BuoFac*(g0+g1*r/r_cmb+g2*(r_cmb/r)**2)
         dentropy0=0.0_cp

         if (l_anel) then
            if (l_isothermal) then
               DissNb=strat /( g0+half*g1*(one+radratio) +g2/radratio )
               ViscHeatFac=0.0_cp
               temp0=one
               rho0=exp(-DissNb*(g0*(r-r_cmb) + &
                    g1/(two*r_cmb)*(r**2-r_cmb**2) - &
                    g2*(r_cmb**2/r-r_cmb)))

               beta =-DissNb*rgrav/BuoFac
               dbeta=-DissNb*(g1/r_cmb-two*g2*r_cmb**2*or3)
               dtemp0=0.0_cp
               d2temp0=0.0_cp
            else
               DissNb=( exp(strat/polind)-one )/ &
                 ( g0+half*g1*(one+radratio) +g2/radratio )
               ViscHeatFac=DissNb*pr/raScaled
               temp0=-DissNb*( g0*r+half*g1*r**2/r_cmb-g2*r_cmb**2/r ) + &
                      one + DissNb*r_cmb*(g0+half*g1-g2)
               rho0=temp0**polind

               !-- Computation of beta= dln rho0 /dr and dbeta=dbeta/dr
               beta=-polind*DissNb*rgrav/temp0/BuoFac
               dbeta=-polind*DissNb/temp0**2 * &
                    ((g1/r_cmb-two*g2*r_cmb**2*or3)* &
                    temp0  + DissNb*rgrav**2/BuoFac**2)
               dtemp0=-DissNb*rgrav/BuoFac
               d2temp0=-DissNb*(g1/r_cmb-two*g2*r_cmb**2*or3)
            end if
            if (l_mag) then
               OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
            end if
         end if
      end if

      agrav=alpha*rgrav

      if ( .not. l_heat ) then
         rgrav=0.0_cp
         agrav=0.0_cp
      end if

      !-- Get additional functions of r:
      if ( l_anel ) then
         orho1=one/rho0
         orho2=orho1*orho1
         otemp1=one/temp0
      else
         rho0     =one
         temp0    =one
         otemp1   =one
         orho1    =one
         orho2    =one
         beta     =0.0_cp
         dbeta    =0.0_cp
         dtemp0   =0.0_cp
         d2temp0  =0.0_cp
         dentropy0=0.0_cp
      end if

      !-- Factors for cheb integrals:
      cheb_int(1)=one   ! Integration constant chosen !
      do n_cheb=3,n_r_max,2
         cheb_int(n_cheb)  =-one/real(n_cheb*(n_cheb-2),kind=cp)
         cheb_int(n_cheb-1)= 0.0_cp
      end do


      !-- Proceed with inner core:

      if ( n_r_ic_max > 0 ) then

         n_r_ic_tot=2*n_r_ic_max-1

         !----- cheb_grid calculates the n_r_ic_tot gridpoints,
         !      these are the extrema of a Cheb of degree n_r_ic_tot-1.
         call cheb_grid(-r_icb,r_icb,n_r_ic_tot-1, &
                         r_ic_2,r_cheb_ic,0.0_cp,0.0_cp,0.0_cp,0.0_cp)

         !----- Store first n_r_ic_max points of r_ic_2 to r_ic:
         do n_r=1,n_r_ic_max-1
            r_ic(n_r)   =r_ic_2(n_r)
            O_r_ic(n_r) =one/r_ic(n_r)
            O_r_ic2(n_r)=O_r_ic(n_r)*O_r_ic(n_r)
         end do
         n_r=n_r_ic_max
         r_ic(n_r)   =0.0_cp
         O_r_ic(n_r) =0.0_cp
         O_r_ic2(n_r)=0.0_cp

         !-- Get no of point on graphical output grid:
         !   No is set to -1 to indicate that point is not on graphical output grid.

      end if

      if ( n_r_ic_max > 0 .and. l_cond_ic ) then

         dr_fac_ic=two/(two*r_icb)
         cheb_norm_ic=sqrt(two/real(n_r_ic_max-1,kind=cp))

         !----- Calculate the even Chebs and their derivative:
         !      n_r_ic_max even chebs up to degree 2*n_r_ic_max-2
         !      at the n_r_ic_max first points in r_ic [r_icb,0].
         !      NOTE invers order in r_ic!
         call get_chebs_even(n_r_ic_max,-r_icb,r_icb,r_cheb_ic, &
                                   n_r_ic_max,cheb_ic,dcheb_ic, &
                                d2cheb_ic,n_r_ic_max,n_r_ic_max)

         !----- Initialize transforms:
         call init_costf1(n_r_ic_max,i_costf1_ic_init,nDi_costf1_ic, &
                          d_costf1_ic_init,nDd_costf1_ic)
         call init_costf2(n_r_ic_max-1,i_costf2_ic_init,nDi_costf2_ic, &
                          d_costf2_ic_init,nDd_costf2_ic)


         !----- Factors for cheb integrals, only even contribution:
         fac_int=one/dr_fac_ic   ! thats 1 for the outer core
         cheb_int_ic(1)=fac_int   ! Integration constant chosen !
         do n_cheb=2,n_r_ic_max
            n_cheb_int=2*n_cheb-1
            cheb_int_ic(n_cheb)=-fac_int / real(n_cheb_int*(n_cheb_int-2),kind=cp)
         end do

      end if

   end subroutine radial
!------------------------------------------------------------------------------
   subroutine transportProperties

      integer :: n_r

      real(cp) :: a,b,c,s1,s2,r0
      real(cp) :: dsigma0
      real(cp) :: dvisc(n_r_max), dkappa(n_r_max), dsigma(n_r_max)
      real(cp) :: condBot(n_r_max), condTop(n_r_max)
      real(cp) :: func(n_r_max), kcond(n_r_max)
      real(cp) :: a0,a1,a2,a3,a4,a5
      real(cp) :: kappatop,rrOcmb
      real(cp) :: w1(n_r_max),w2(n_r_max)

      !-- Variable conductivity:

      if ( imagcon == -10 ) then
         nVarCond=1
         lambda  =r**5.0_cp
         sigma   =one/lambda
         dLlambda=5.0_cp/r
      else if ( l_mag ) then
          if ( nVarCond == 0 ) then
             lambda  =one
             sigma   =one
             dLlambda=0.0_cp
          else if ( nVarCond == 1 ) then
             b =log(three)/con_FuncWidth
             r0=con_radratio*r_cmb
             s1=tanh(b*(r0-r_cmb))
             s2=tanh(b*(r0-r_icb))
             a =(-one+con_LambdaOut)/(s1-s2)
             c =(s1-s2*con_LambdaOut)/(s1-s2)
             sigma   = a*tanh(b*(r0-r))+c
             dsigma  =-a*b/cosh(b*(r0-r))
             lambda  =one/sigma
             dLlambda=-dsigma/sigma
          else if ( nVarCond == 2 ) then

             r0=con_radratio*r_cmb
             !------ Use grid point closest to r0:
             do n_r=1,n_r_max
                if ( r(n_r) < r0 )then
                   r0=r(n_r)
                   EXIT
                end if
             end do
             dsigma0=(con_LambdaMatch-1)*con_DecRate /(r0-r_icb)
             do n_r=1,n_r_max
                if ( r(n_r) < r0 ) then
                   sigma(n_r)   = one+(con_LambdaMatch-1)* &
                       ( (r(n_r)-r_icb)/(r0-r_icb) )**con_DecRate
                   dsigma(n_r)  =  dsigma0 * &
                       ( (r(n_r)-r_icb)/(r0-r_icb) )**(con_DecRate-1)
                else
                   sigma(n_r)  =con_LambdaMatch * &
                       exp(dsigma0/con_LambdaMatch*(r(n_r)-r0))
                   dsigma(n_r) = dsigma0* &
                       exp(dsigma0/con_LambdaMatch*(r(n_r)-r0))
                end if
                lambda(n_r)  = one/sigma(n_r)
                dLlambda(n_r)=-dsigma(n_r)/sigma(n_r)
             end do
          else if ( nVarCond == 3 ) then ! Magnetic diff propto 1/rho
             lambda=rho0(n_r_max)/rho0
             sigma=one/lambda
             call get_dr(lambda,dsigma,n_r_max,n_cheb_max, &
                         w1,w2,i_costf_init,d_costf_init,drx)
             dLlambda=dsigma/lambda
          else if ( nVarCond == 4 ) then ! Profile
             lambda=(rho0/rho0(n_r_max))**difExp
             sigma=one/lambda
             call get_dr(lambda,dsigma,n_r_max,n_cheb_max, &
                         w1,w2,i_costf_init,d_costf_init,drx)
             dLlambda=dsigma/lambda
          end if
      end if

      !-- Variable thermal diffusivity
      if ( l_heat ) then
         if ( nVarDiff == 0 ) then
            kappa  =one
            dLkappa=0.0_cp
         else if ( nVarDiff == 1 ) then ! Constant conductivity
            ! kappa(n_r)=one/rho0(n_r) Denise's version
            kappa=rho0(n_r_max)/rho0
            call get_dr(kappa,dkappa,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 2 ) then ! Profile
            kappa=(rho0/rho0(n_r_max))**difExp
            call get_dr(kappa,dkappa,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 3 ) then ! polynomial fit to a model
            if ( radratio < 0.19 ) then
               write(*,*) '! NOTE: with this polynomial fit     '
               write(*,*) '! for variable thermal conductivity  '
               write(*,*) '! considering radratio < 0.2 may lead'
               write(*,*) '! to strange profiles'
               stop
            end if
            a0 = -0.32839722_cp
            a1 =  one
            a2 = -1.16153274_cp
            a3 =  0.63741485_cp
            a4 = -0.15812944_cp
            a5 =  0.01034262_cp
            do n_r=1,n_r_max
               rrOcmb = r(n_r)/r_cmb*r_cut_model
               kappa(n_r)= a5 + a4*rrOcmb    + a3*rrOcmb**2 &
                              + a2*rrOcmb**3 + a1*rrOcmb**4 &
                                             + a0*rrOcmb**5

            ENDDO
            kappatop=kappa(1) ! normalise by the value at the top
            kappa=kappa/kappatop
            call get_dr(kappa,dkappa,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 4) then ! Earth case
            !condTop=r_cmb**2*dtemp0(1)*or2/dtemp0
            !do n_r=2,n_r_max
            !  if ( r(n_r-1)>rStrat .and. r(n_r)<=rStrat ) then
            !     if ( r(n_r-1)-rStrat < rStrat-r(n_r) ) then
            !        n_const=n_r-1
            !     else
            !        n_const=n_r
            !     end if
            !  end if
            !end do
            !condBot=(rho0/rho0(n_const))*condTop(n_const)
            !func=half*(tanh(slopeStrat*(r-rStrat))+one)
            !kcond=condTop*func+condBot*(1-func)
            !kcond=kcond/kcond(n_r_max)
            !kappa=kcond/rho0
            !call get_dr(kappa,dkappa,n_r_max,n_cheb_max, &
            !            w1,w2,i_costf_init,d_costf_init,drx)
            !dLkappa=dkappa/kappa

            ! Alternative scenario
            kcond=(one-0.4469_cp*(r/r_cmb)**2)/(one-0.4469_cp)
            kcond=kcond/kcond(1)
            kappa=kcond/rho0
            call get_dr(kappa,dkappa,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
         end if
      end if

      !-- Eps profiles
      !-- The remaining division by rho will happen in s_updateS.F90
      if ( nVarEps == 0 ) then
         ! eps is constant
         epscProf=otemp1
      else if ( nVarEps == 1 ) then
         ! rho*eps in the RHS
         epscProf=rho0*otemp1
      end if

      !-- Variable viscosity
      if ( nVarVisc == 0 ) then ! default: constant kinematic viscosity
         visc  =one
         dLvisc=0.0_cp
      else if ( nVarVisc == 1 ) then ! Constant dynamic viscosity
         visc=rho0(n_r_max)/rho0
         call get_dr(visc,dvisc,n_r_max,n_cheb_max, &
                     w1,w2,i_costf_init,d_costf_init,drx)
         dLvisc=dvisc/visc
      else if ( nVarVisc == 2 ) then ! Profile
         visc=(rho0/rho0(n_r_max))**difExp
         call get_dr(visc,dvisc,n_r_max,n_cheb_max, &
                     w1,w2,i_costf_init,d_costf_init,drx)
         dLvisc=dvisc/visc
      end if

      if ( l_anelastic_liquid ) then
         divKtemp0=rho0*kappa*(d2temp0+(beta+dLkappa+two*or1)*dtemp0)*sq4pi
      else
         divKtemp0=0.0_cp
      end if

   end subroutine transportProperties
!------------------------------------------------------------------------------
   subroutine getBackground(input,boundaryVal,output)

      !-- Input variables:
      real(cp), intent(in) :: input(n_r_max)
      real(cp), intent(in) :: boundaryVal

      !-- Output variables:
      real(cp), intent(out) :: output(n_r_max)

      !-- Local variables:
      real(cp) :: rhs(n_r_max)
      real(cp) :: tmp(n_r_max)
      integer :: n_cheb,n_r,info


      do n_cheb=1,n_r_max
         do n_r=2,n_r_max
            s0Mat(n_r,n_cheb)=cheb_norm*dcheb(n_cheb,n_r)
         end do
      end do

      !-- boundary conditions
      do n_cheb=1,n_cheb_max
         s0Mat(1,n_cheb)=cheb_norm
         s0Mat(n_r_max,n_cheb)=0.0_cp
      end do

      !-- fill with zeros
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            s0Mat(1,n_cheb)=0.0_cp
         end do
      end if

      !-- renormalize
      do n_r=1,n_r_max
         s0Mat(n_r,1)      =half*s0Mat(n_r,1)
         s0Mat(n_r,n_r_max)=half*s0Mat(n_r,n_r_max)
      end do

#ifdef WITH_MKL_LU
      call getrf(s0Mat,s0Pivot,info)
#else
      call sgefa(s0Mat,n_r_max,n_r_max,s0Pivot,info)
#endif
      if ( info /= 0 ) then
         write(*,*) '! Singular Matrix in getBackground!'
         stop '20'
      end if

      do n_r=2,n_r_max
         rhs(n_r)=input(n_r)
      end do
      rhs(1)=boundaryVal

      !-- Solve for s0:
#ifdef WITH_MKL_LU
      call getrs(s0Mat,s0Pivot,rhs)
#else
      call sgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
#endif

      !-- Copy result to s0:
      do n_r=1,n_r_max
         output(n_r)=rhs(n_r)
      end do

      !-- Set cheb-modes > n_cheb_max to zero:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            output(n_cheb)=0.0_cp
         end do
      end if

      !-- Transform to radial space:
      call costf1(output,tmp,i_costf_init,d_costf_init)

   end subroutine getBackground
!------------------------------------------------------------------------------
end module radial_functions
