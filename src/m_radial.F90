!$Id$
!*******************************************************************
!  Common block containing radial functions, all THREADPRIVATE
!  These arrays and numbers are calculated by routine s_radial.f.
!*******************************************************************
MODULE radial_functions
  use truncation
  use radial_data
  implicit none

  !-- arrays depending on r:
  REAL(kind=8),allocatable :: r(:)
  REAL(kind=8),ALLOCATABLE :: r_ic(:)
  REAL(kind=8),ALLOCATABLE :: O_r_ic(:)
  REAL(kind=8),ALLOCATABLE :: O_r_ic2(:)
  REAL(kind=8),ALLOCATABLE :: or1(:),or2(:),or3(:),or4(:)
  REAL(kind=8),ALLOCATABLE :: otemp1(:),rho0(:),temp0(:)
  REAL(kind=8),ALLOCATABLE :: orho1(:),orho2(:)
  REAL(kind=8),ALLOCATABLE :: beta(:), dbeta(:)
  REAL(kind=8),ALLOCATABLE :: drx(:),ddrx(:),dddrx(:)
  REAL(kind=8) :: dr_fac,dr_fac_ic,alpha1,alpha2
  REAL(kind=8) :: topcond, botcond
  REAL(kind=8) :: r_cmb,r_icb,r_surface

  !-- arrays for buoyancy, depend on Ra and Pr:
  REAL(kind=8),allocatable :: rgrav(:),agrav(:)
  !$OMP THREADPRIVATE( r,r_ic,or1,or2,or3,or4,O_r_ic,O_r_ic2 )            
  !$OMP THREADPRIVATE( rgrav,agrav )             
  !$OMP THREADPRIVATE( drx,ddrx,dddrx,alpha1,alpha2 )
  !$OMP THREADPRIVATE( dr_fac,dr_fac_ic,r_cmb,r_icb,r_surface )
  !$OMP THREADPRIVATE( topcond, botcond )
  !$OMP THREADPRIVATE( otemp1,temp0,rho0,orho1,orho2,beta,dbeta )


  !-- chebychev polynomials, derivatives and integral:
  REAL(kind=8) :: cheb_norm                   ! cheb normalisation 
  REAL(kind=8),ALLOCATABLE :: cheb(:,:)       ! Chebychev polynomials
  REAL(kind=8),ALLOCATABLE :: dcheb(:,:)      ! first radial derivative
  REAL(kind=8),ALLOCATABLE :: d2cheb(:,:)     ! second radial derivative
  REAL(kind=8),ALLOCATABLE :: d3cheb(:,:)     ! third radial derivative
  REAL(kind=8),ALLOCATABLE :: cheb_int(:)           ! array for cheb integrals !
  INTEGER nDi_costf1
  !PARAMETER (nDi_costf1=2*n_r_max+2)
  INTEGER,ALLOCATABLE :: i_costf_init(:)   ! info for transform
  INTEGER :: nDd_costf1
  !PARAMETER (nDd_costf1=2*n_r_max+5)
  REAL(kind=8),ALLOCATABLE ::  d_costf_init(:)   ! info for tranform
  !$OMP THREADPRIVATE( cheb_norm,cheb,dcheb,d2cheb,d3cheb,cheb_int )
  !$OMP THREADPRIVATE( d_costf_init,i_costf_init )

  !-- same for inner core:
  REAL(kind=8) :: cheb_norm_ic
  REAL(kind=8),ALLOCATABLE :: cheb_ic(:,:)   
  REAL(kind=8),ALLOCATABLE :: dcheb_ic(:,:)  
  REAL(kind=8),ALLOCATABLE :: d2cheb_ic(:,:) 
  REAL(kind=8),ALLOCATABLE :: cheb_int_ic(:)        
  INTEGER nDi_costf1_ic

  INTEGER,allocatable :: i_costf1_ic_init(:)
  INTEGER :: nDd_costf1_ic

  REAL(kind=8),ALLOCATABLE ::  d_costf1_ic_init(:)
  INTEGER :: nDi_costf2_ic

  INTEGER,ALLOCATABLE :: i_costf2_ic_init(:)
  INTEGER :: nDd_costf2_ic

  REAL(kind=8),ALLOCATABLE ::  d_costf2_ic_init(:)
  !$OMP THREADPRIVATE( cheb_norm_ic,cheb_ic,dcheb_ic,d2cheb_ic )
  !$OMP THREADPRIVATE( cheb_int_ic )
  !$OMP THREADPRIVATE( d_costf1_ic_init,d_costf2_ic_init )
  !$OMP THREADPRIVATE( i_costf1_ic_init,i_costf2_ic_init )

  !-- Radius functions for cut-back grid without boundaries:
  !-- (and for the nonlinear mapping)
  INTEGER,allocatable :: i_costf_initC(:)   ! info for transform
  REAL(kind=8),ALLOCATABLE ::  d_costf_initC(:)   ! info for tranform
  REAL(kind=8),ALLOCATABLE ::  rC(:),dr_facC(:)
  REAL(kind=8) :: alph1,alph2
  INTEGER :: n_r_maxC,n_cheb_maxC,nCut

  PUBLIC :: initialize_radial_functions, radial
  REAL(kind=8),ALLOCATABLE :: lambda(:),dLlambda(:),jVarCon(:)
  REAL(kind=8),ALLOCATABLE :: sigma(:)
  REAL(kind=8),ALLOCATABLE :: kappa(:),dLkappa(:)
  REAL(kind=8),ALLOCATABLE :: visc(:),dLvisc(:)
  REAL(kind=8),ALLOCATABLE :: epscProf(:)

  !------------------------------------------------------------------------------
CONTAINS
  SUBROUTINE initialize_radial_functions

    nDi_costf1=2*n_r_max+2
    nDd_costf1=2*n_r_max+5

    nDi_costf1_ic=2*n_r_ic_max+2
    nDd_costf1_ic=2*n_r_ic_max+5
    nDi_costf2_ic=2*n_r_ic_max
    nDd_costf2_ic=2*n_r_ic_max+n_r_ic_max/2+5

    !$OMP PARALLEL
    ! allocate the arrays
    ALLOCATE( r(n_r_max) )
    ALLOCATE( r_ic(n_r_ic_max) )
    ALLOCATE( O_r_ic(n_r_ic_max) )
    ALLOCATE( O_r_ic2(n_r_ic_max) )
    ALLOCATE( or1(n_r_max),or2(n_r_max),or3(n_r_max),or4(n_r_max) )
    ALLOCATE( otemp1(n_r_max),rho0(n_r_max),temp0(n_r_max) )
    ALLOCATE( orho1(n_r_max),orho2(n_r_max) )
    ALLOCATE( beta(n_r_max), dbeta(n_r_max) )
    ALLOCATE( drx(n_r_max),ddrx(n_r_max),dddrx(n_r_max) )
    ALLOCATE( rgrav(n_r_max),agrav(n_r_max) )

    ALLOCATE( cheb(n_r_max,n_r_max) )     ! Chebychev polynomials
    ALLOCATE( dcheb(n_r_max,n_r_max) )    ! first radial derivative
    ALLOCATE( d2cheb(n_r_max,n_r_max) )   ! second radial derivative
    ALLOCATE( d3cheb(n_r_max,n_r_max) )   ! third radial derivative
    ALLOCATE( cheb_int(n_r_max) )         ! array for cheb integrals !
    ALLOCATE( i_costf_init(nDi_costf1) )  ! info for transform
    ALLOCATE( d_costf_init(nDd_costf1) ) ! info for tranform

    ALLOCATE( cheb_ic(n_r_ic_max,n_r_ic_max) )
    ALLOCATE( dcheb_ic(n_r_ic_max,n_r_ic_max) )
    ALLOCATE( d2cheb_ic(n_r_ic_max,n_r_ic_max) )
    ALLOCATE( cheb_int_ic(n_r_ic_max) )
    ALLOCATE( i_costf1_ic_init(nDi_costf1_ic) )
    ALLOCATE( d_costf1_ic_init(nDd_costf1_ic) )
    ALLOCATE( i_costf2_ic_init(nDi_costf2_ic) )
    ALLOCATE( d_costf2_ic_init(nDd_costf2_ic) )
    !$OMP END PARALLEL

    ALLOCATE( lambda(n_r_max),dLlambda(n_r_max),jVarCon(n_r_max) )
    ALLOCATE( sigma(n_r_max) )
    ALLOCATE( kappa(n_r_max),dLkappa(n_r_max) )
    ALLOCATE( visc(n_r_max),dLvisc(n_r_max) )
    ALLOCATE( epscProf(n_r_max) )

    ALLOCATE( i_costf_initC(nDi_costf1) )   ! info for transform
    ALLOCATE( d_costf_initC(nDd_costf1) )   ! info for tranform
    ALLOCATE( rC(n_r_max),dr_facC(n_r_max) )

  END SUBROUTINE initialize_radial_functions
  !********************************************************************
  SUBROUTINE radial
    !********************************************************************

    !--------------------------------------------------------------------
    !  Calculates everything needed for radial functions, transfroms etc.
    !--------------------------------------------------------------------
    USE physical_parameters
    USE num_param
    USE logic
    USE output_data
    USE usefull, only: check_dim

    IMPLICIT NONE

    !-- LOCAL VARIABLES:
    INTEGER :: n_r,n_cheb,n_cheb_int
    INTEGER :: n_r_ic_tot
    INTEGER,PARAMETER :: n_r_ic_tot_local=200
    !INTEGER :: n_r_start
    REAL(kind=8) :: fac_int
    REAL(kind=8) :: r_cheb(n_r_max)
    REAL(kind=8) :: r_cheb_ic(n_r_ic_tot_local)
    REAL(kind=8) :: r_ic_2(n_r_ic_tot_local)
    REAL(kind=8) :: ofr ! inverse of the Froude number
    REAL(kind=8) :: dtemp0(n_r_max),alphaT(n_r_max)
    REAL(kind=8) :: drho0(n_r_max)
    REAL(kind=8) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

    REAL(kind=8) :: a0,a1,a2,a3,a4,a5,a6,a7 ! polynomial fit for density
    REAL(kind=8) :: temptop,gravtop,rhotop
    REAL(kind=8) :: rrOcmb,gravFit(n_r_max),rhoFit(n_r_max) ! normalised to rcmb
    REAL(kind=8) :: w1(n_r_max),w2(n_r_max)

    INTEGER :: stop_signal

    !-- end of declaration
    !-------------------------------------------------------------------

    !-- Polynomial fit for a Jupiter model

    IF ( l_interior_model) THEN
       a0=-122.36071577d0
       a1= 440.86067831d0
       a2=-644.82401602d0
       a3= 491.00495215d0
       a4=-201.96655354d0
       a5=  37.38863965d0
       a6=  -4.60312999d0
       a7=   4.46020423d0
    END IF

    !-- Radial grid point:
    !   radratio is aspect ratio
    !   radratio = (inner core r) / (CMB r) = r_icb/r_cmb
    r_cmb=1.D0/(1.D0-radratio)
    r_icb=r_cmb-1.D0
    r_surface=2.8209D0    ! in units of (r_cmb-r_icb)

    cheb_norm=DSQRT(2.D0/DBLE(n_r_max-1))
    dr_fac=2.D0/(r_cmb-r_icb)

    if (l_newmap) then
       alpha1         =alph1
       alpha2         =alph2
       paraK=DATAN(alpha1*(1+alpha2))/DATAN(alpha1*(1-alpha2))
       paraX0=(paraK-1)/(paraK+1)
       lambd=DATAN(alpha1*(1-alpha2))/(1-paraX0)
    else
       alpha1         =0.D0
       alpha2         =0.D0
    end if

    !-- Start with outer core:
    !   cheb_x_map_e calculates the n_r_max gridpoints, these
    !   are the extrema of a Cheb pylonomial of degree n_r_max-1,
    !   r_cheb are the grid_points in the Cheb intervall [-1,1]
    !   and r are these points mapped to the intervall [r_icb,r_cmb]:
    CALL cheb_x_map_e_xr(r_icb,r_cmb,n_r_max-1,r,r_cheb, &
         !                         alpha1,alpha2,paraK,paraX0,lambd)
         alpha1,alpha2,paraX0,lambd)

    if (l_newmap) then
       do n_r=1,n_r_max
          drx(n_r) =                          (2*alpha1) / &
               ((1+alpha1**2*(2*r(n_r)-r_icb-r_cmb-alpha2)**2)* &
               lambd)
          ddrx(n_r) = -(8*alpha1**3*(2*r(n_r)-r_icb-r_cmb-alpha2)) / &
               ((1+alpha1**2*(-2*r(n_r)+r_icb+r_cmb+alpha2)**2)**2* &
               lambd)
          dddrx(n_r) =          (16*alpha1**3*(-1+3*alpha1**2* &
               (-2*r(n_r)+r_icb+r_cmb+alpha2)**2)) / &
               ((1+alpha1**2*(-2*r(n_r)+r_icb+r_cmb+alpha2)**2)**3* &
               lambd)
       end do
    else
       do n_r=1,n_r_max
          drx(n_r)=2.D0/(r_cmb-r_icb)
          ddrx(n_r)=0
          dddrx(n_r)=0
       end do
    end if

    !-- Calculate chebs and their derivatives up to degree n_r_max-1
    !   on the n_r radial grid points r:
    CALL get_chebs(n_r_max,r_icb,r_cmb,r_cheb,n_r_max, &
             cheb,dcheb,d2cheb,d3cheb,n_r_max,n_r_max, &
                                       drx,ddrx,dddrx)

    !-- Initialize fast cos transform for chebs:
    CALL init_costf1(n_r_max,i_costf_init,nDi_costf1, &
                     d_costf_init,nDd_costf1)

    !-- Fit to an interior model
    IF ( l_interior_model ) THEN

       DO n_r=1,n_r_max
          rrOcmb = r(n_r)/r_cmb*r_cut_model
          rhoFit(n_r) = a7 + a6*rrOcmb   + a5*rrOcmb**2 + a4*rrOcmb**3 &
                           + a3*rrOcmb**4+ a2*rrOcmb**5 + a1*rrOcmb**6 &
                           + a0*rrOcmb**7
          gravFit(n_r)=4.d0*rrOcmb - 3.d0*rrOcmb**2
          temp0(n_r)=rhoFit(n_r)**(1.d0/polind)
       END DO

       ! To normalise to the outer radius
       temptop=temp0(1)
       rhotop=rhoFit(1)
       gravtop=gravFit(1)

       temp0  =temp0/temptop
       gravFit=gravFit/gravtop

       ! Derivative of the temperature needed to get alpha_T
       CALL get_dr(temp0,dtemp0,1,1,1,n_r_max,n_cheb_max,w1, &
                   w2,i_costf_init,d_costf_init,drx)

       alphaT=-dtemp0/(gravFit*temp0)

       ! Inverse of the Froude number needed in the dissipation numbers
       ofr=alphaT(1)
       ViscHeatFac=ofr*pr/raScaled
       IF (l_mag) THEN
          OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
       END IF

       !-- Usual polytropic reference state
    ELSE
       IF (l_anel) THEN
          ofr=( DEXP(strat/polind)-1.D0 )/ &
               ( g0+0.5D0*g1*(1.D0+radratio) +g2/radratio )
          IF (l_isothermal) THEN
             ofr=strat /( g0+0.5D0*g1*(1.D0+radratio) +g2/radratio )
             ViscHeatFac=0.D0
          ELSE
             ViscHeatFac=ofr*pr/raScaled
          END IF
          IF (l_mag) THEN
             OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
          END IF
       END IF
    END IF

    !-- Get additional functions of r:
    DO n_r=1,n_r_max
       IF ( l_heat ) THEN
          IF ( l_interior_model ) THEN
             ! Adiabatic: buoyancy term is linked to the temperature gradient

             !       dT
             ! Fr * ---- =  alpha_T * T * grav
             !       dr

             ! N.B. rgrav is not gravity but the whole RHS !!!

             rgrav(n_r)=-BuoFac * dtemp0(n_r)/ofr
          ELSE
             ! g(r) = g0 + g1*r/ro + g2*(ro/r)**2
             ! Default values: g0=0, g1=1, g2=0
             ! An easy way to change gravity
             rgrav(n_r)=BuoFac*( g0 + &
                  g1*r(n_r)/r_cmb + &
                  g2*(r_cmb/r(n_r))**2 )
          END IF
          agrav(n_r)=alpha*rgrav(n_r)
       ELSE
          rgrav(n_r)=0.D0
          agrav(n_r)=0.D0
       END IF
       IF ( l_anel ) THEN
          IF ( l_isothermal ) THEN
             temp0(n_r)=1.D0
             rho0(n_r)=DEXP(-ofr*(g0*(r(n_r)-r_cmb) + &
                  g1/(2.d0*r_cmb)*(r(n_r)**2-r_cmb**2) - &
                  g2*(r_cmb**2/r(n_r)-r_cmb)))
          ELSE IF ( l_interior_model ) THEN
             rho0(n_r)=rhoFit(n_r)/rhotop
             ! temp0(n_r)=rho0(n_r)**(1.d0/polind)
          ELSE
             temp0(n_r)=-ofr*( g0*r(n_r) + &
                  0.5D0*g1*r(n_r)**2/r_cmb - &
                  g2*r_cmb**2/r(n_r) ) + &
                  1.D0 + ofr*r_cmb*(g0+0.5D0*g1-g2)
             rho0(n_r)=temp0(n_r)**polind
          END IF
          orho1(n_r)=1.D0/rho0(n_r)
          orho2(n_r)=orho1(n_r)*orho1(n_r)
          otemp1(n_r)=1.D0/temp0(n_r)
       ELSE
          rho0(n_r)=1.D0
          temp0(n_r)=1.D0
          otemp1(n_r)=1.D0
          orho1(n_r)=1.D0
          orho2(n_r)=1.D0
       END IF
       or1(n_r)=1.D0/r(n_r)             ! 1/r
       or2(n_r)=or1(n_r)*or1(n_r)       ! 1/r**2
       or3(n_r)=or1(n_r)*or2(n_r)       ! 1/r**3
       or4(n_r)=or2(n_r)*or2(n_r)       ! 1/r**4
    END DO

    !-- Computation of beta= dln rho0 /dr and dbeta=dbeta/dr
    IF (l_anel) THEN
       ! better than the derivative the log(rho0)
       ! especially for large stratification Nrho=5
       IF ( l_isothermal ) THEN
          beta =-ofr*rgrav/BuoFac
          dbeta=-ofr*(g1/r_cmb-2.D0*g2*r_cmb**2*or3)
       ELSE IF ( l_interior_model ) THEN
          CALL get_dr(rho0,drho0,1,1,1,n_r_max,n_cheb_max,w1, &
                      w2,i_costf_init,d_costf_init,drx)
          beta=drho0/rho0
          CALL get_dr(beta,dbeta,1,1,1,n_r_max,n_cheb_max,w1, &
                      w2,i_costf_init,d_costf_init,drx)
       ELSE
         beta=-polind*ofr*rgrav/temp0/BuoFac
         dbeta=-polind*ofr/temp0**2 * &
              ((g1/r_cmb-2.D0*g2*r_cmb**2*or3)* &
              temp0  + ofr*rgrav**2/BuoFac**2)
       END IF

    ELSE
       beta =0.D0
       dbeta=0.D0
    END IF

    !-- Factors for cheb integrals:
    cheb_int(1)=1.d0   ! Integration constant chosen !
    DO n_cheb=3,n_r_max,2
       cheb_int(n_cheb)  =-1.D0/DBLE(n_cheb*(n_cheb-2))
       cheb_int(n_cheb-1)= 0.D0
    END DO


    !-- Proceed with inner core:

    IF ( n_r_ic_max > 0 ) THEN

       n_r_ic_tot=2*n_r_ic_max-1
       CALL check_dim(n_r_ic_tot,n_r_ic_tot_local, &
            'n_r_ic_tot_local','rad_func',stop_signal)
       IF ( stop_signal == 1 ) STOP


       !----- cheb_x_map_e calculates the n_r_ic_tot gridpoints,
       !      these are the extrema of a Cheb of degree n_r_ic_tot-1.
       CALL cheb_x_map_e(-r_icb,r_icb,n_r_ic_tot-1, &
            r_ic_2,r_cheb_ic)

       !----- Store first n_r_ic_max points of r_ic_2 to r_ic:
       DO n_r=1,n_r_ic_max-1
          r_ic(n_r)   =r_ic_2(n_r)
          O_r_ic(n_r) =1.D0/r_ic(n_r)
          O_r_ic2(n_r)=O_r_ic(n_r)*O_r_ic(n_r)
       END DO
       n_r=n_r_ic_max
       r_ic(n_r)   =0.D0
       O_r_ic(n_r) =0.D0
       O_r_ic2(n_r)=0.D0

       !-- Get no of point on graphical output grid:
       !   No is set to -1 to indicate that point is not on graphical output grid.

    END IF

    if ( n_r_ic_max > 0 .AND. l_cond_ic ) then

       dr_fac_ic=2.D0/(2.D0*r_icb)
       cheb_norm_ic=DSQRT(2.D0/DBLE(n_r_ic_max-1))

       !----- Calculate the even Chebs and their derivative:
       !      n_r_ic_max even chebs up to degree 2*n_r_ic_max-2
       !      at the n_r_ic_max first points in r_ic [r_icb,0].
       !      NOTE invers order in r_ic!
       CALL get_chebs_even(n_r_ic_max,-r_icb,r_icb,r_cheb_ic, &
            n_r_ic_max,cheb_ic,dcheb_ic, &
            d2cheb_ic,n_r_ic_max,n_r_ic_max)

       !----- Initialize transforms:
       CALL init_costf1(n_r_ic_max, &
            i_costf1_ic_init,nDi_costf1_ic, &
            d_costf1_ic_init,nDd_costf1_ic)
       CALL init_costf2(n_r_ic_max-1, &
            i_costf2_ic_init,nDi_costf2_ic, &
            d_costf2_ic_init,nDd_costf2_ic)


       !----- Factors for cheb integrals, only even contribution:
       fac_int=1.D0/dr_fac_ic   ! thats 1 for the outer core
       cheb_int_ic(1)=fac_int   ! Integration constant chosen !
       DO n_cheb=2,n_r_ic_max
          n_cheb_int=2*n_cheb-1
          cheb_int_ic(n_cheb)=-fac_int / &
               dble(n_cheb_int*(n_cheb_int-2))
       END DO

    END IF

    RETURN
  END SUBROUTINE radial
  !--------------------------------------------------------------------
  SUBROUTINE transportProperties

    USE physical_parameters
    USE logic, ONLY: l_mag,l_heat
    USE init_fields, ONLY: imagcon

    IMPLICIT NONE

    INTEGER :: n_r

    REAL(kind=8) :: a,b,c,s1,s2,r0
    REAL(kind=8) :: dsigma0
    REAL(kind=8) :: dsigma(n_r_max)
    REAL(kind=8) :: dvisc(n_r_max),dkappa(n_r_max)
    REAL(kind=8) :: a0,a1,a2,a3,a4,a5
    REAL(kind=8) :: kappatop,rrOcmb
    REAL(kind=8) :: w1(n_r_max),w2(n_r_max)

!-- Variable conductivity:
    IF ( imagcon == -10 ) THEN
        nVarCond=1
        lambda  =r**5.D0
        sigma   =1.D0/lambda
        dLlambda=5.D0/r
    ELSE IF ( l_mag ) THEN
        IF ( nVarCond == 0 ) THEN
            lambda  =1.D0
            sigma   =1.D0
            dLlambda=0.D0
        ELSE IF ( nVarCond == 1 ) THEN
            b =DLOG(3.D0)/con_FuncWidth
            r0=con_radratio*r_cmb
            s1=TANH(b*(r0-r_cmb))
            s2=TANH(b*(r0-r_icb))
            a =(-1.D0+con_LambdaOut)/(s1-s2)
            c =(s1-s2*con_LambdaOut)/(s1-s2)
            sigma   = a*TANH(b*(r0-r))+c
            dsigma  =-a*b/COSH(b*(r0-r))
            lambda  =1.D0/sigma
            dLlambda=-dsigma/sigma
        ELSE IF ( nVarCond == 2 ) THEN

            r0=con_radratio*r_cmb
        !------ Use grid point closest to r0:
            DO n_r=1,n_r_max
                IF ( r(n_r) < r0 )THEN
                    r0=r(n_r)
                    EXIT
                END IF
            END DO
            dsigma0=(con_LambdaMatch-1)*con_DecRate /(r0-r_icb)
            DO n_r=1,n_r_max
                IF ( r(n_r) < r0 ) THEN
                    sigma(n_r)   = 1.D0+(con_LambdaMatch-1)* &
                        ( (r(n_r)-r_icb)/(r0-r_icb) )**con_DecRate
                    dsigma(n_r)  =  dsigma0 * &
                        ( (r(n_r)-r_icb)/(r0-r_icb) )**(con_DecRate-1)
                ELSE
                    sigma(n_r)  =con_LambdaMatch * &
                        DEXP(dsigma0/con_LambdaMatch*(r(n_r)-r0))
                    dsigma(n_r) = dsigma0* &
                        DEXP(dsigma0/con_LambdaMatch*(r(n_r)-r0))
                END IF
                lambda(n_r)  = 1.D0/sigma(n_r)
                dLlambda(n_r)=-dsigma(n_r)/sigma(n_r)
            END DO
        ELSE IF ( nVarCond == 3 ) THEN ! Magnetic diff propto 1/rho
            lambda=rho0(n_r_max)/rho0
            sigma=1.d0/lambda
            CALL get_dr(lambda,dsigma,1,1,1,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLlambda=dsigma/lambda
        ELSE IF ( nVarCond == 4 ) THEN ! Profile
            lambda=(rho0/rho0(n_r_max))**difExp
            sigma=1.d0/lambda
            CALL get_dr(lambda,dsigma,1,1,1,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLlambda=dsigma/lambda
        END IF
    END IF

!-- Variable thermal diffusivity
    IF ( l_heat ) THEN
        IF ( nVarDiff == 0 ) THEN
            kappa  =1.D0
            dLkappa=0.D0
        ELSE IF ( nVarDiff == 1 ) THEN ! Constant conductivity
            ! kappa(n_r)=1.d0/rho0(n_r) Denise's version
            kappa=rho0(n_r_max)/rho0
            CALL get_dr(kappa,dkappa,1,1,1,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
        ELSE IF ( nVarDiff == 2 ) THEN ! Profile
            kappa=(rho0/rho0(n_r_max))**difExp
            CALL get_dr(kappa,dkappa,1,1,1,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
        ELSE IF ( nVarDiff == 3 ) THEN ! polynomial fit to a model
            IF ( radratio < 0.19 ) THEN
                WRITE(*,*) '! NOTE: with this polynomial fit     '
                WRITE(*,*) '! for variable thermal conductivity  '
                WRITE(*,*) '! considering radratio < 0.2 may lead'
                WRITE(*,*) '! to strange profiles'
                STOP
            END IF
            a0 = -0.32839722d0
            a1 =  1.d0
            a2 = -1.16153274d0
            a3 =  0.63741485d0
            a4 = -0.15812944d0
            a5 =  0.01034262d0
            DO n_r=1,n_r_max
                rrOcmb = r(n_r)/r_cmb*r_cut_model
                kappa(n_r)= a5 + a4*rrOcmb    + a3*rrOcmb**2 &
                               + a2*rrOcmb**3 + a1*rrOcmb**4 &
                                              + a0*rrOcmb**5

            ENDDO
            kappatop=kappa(1) ! normalise by the value at the top
            kappa=kappa/kappatop
            CALL get_dr(kappa,dkappa,1,1,1,n_r_max,n_cheb_max, &
                        w1,w2,i_costf_init,d_costf_init,drx)
            dLkappa=dkappa/kappa
        END IF
    END IF

!-- Eps profiles
    IF ( nVarEps == 0 ) THEN
        ! eps is constant
        epscProf=orho1*otemp1
    ELSE IF ( nVarEps == 1 ) THEN
        ! rho*eps is constant
        epscProf=otemp1
    END IF

!-- Variable viscosity
    IF ( nVarVisc == 0 ) THEN ! default: constant kinematic viscosity
        visc  =1.d0
        dLvisc=0.d0
    ELSE IF ( nVarVisc == 1 ) THEN ! Constant dynamic viscosity
        visc=rho0(n_r_max)/rho0
        CALL get_dr(visc,dvisc,1,1,1,n_r_max,n_cheb_max, &
                    w1,w2,i_costf_init,d_costf_init,drx)
        dLvisc=dvisc/visc
    ELSE IF ( nVarVisc == 2 ) THEN ! Profile
        visc=(rho0/rho0(n_r_max))**difExp
        CALL get_dr(visc,dvisc,1,1,1,n_r_max,n_cheb_max, &
                    w1,w2,i_costf_init,d_costf_init,drx)
        dLvisc=dvisc/visc
    END IF
  END SUBROUTINE transportProperties
END MODULE radial_functions
