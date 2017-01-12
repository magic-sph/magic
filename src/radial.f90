module radial_functions
   !
   !  This module initiates all the radial functions (transport properties, density,
   !  temperature, cheb transforms, etc.)
   !

   use truncation, only: n_r_max, n_cheb_max, n_r_ic_max, fd_ratio, &
       &                 fd_stretch, fd_order
   use algebra, only: sgesl,sgefa
   use constants, only: sq4pi, one, two, three, four, half
   use physical_parameters
   use logic, only: l_mag, l_cond_ic, l_heat, l_anelastic_liquid,  &
       &            l_isothermal, l_anel, l_newmap, l_non_adia,    &
       &            l_TP_form, l_temperature_diff, l_single_matrix,&
       &            l_finite_diff
   use chebyshev_polynoms_mod ! Everything is needed
   use cosine_transform_odd
   use cosine_transform_even
   use radial_scheme, only: type_rscheme
   use chebyshev, only: type_cheb_odd
   use f_differences, only: type_fd
   !-- remove useless calls to ddr and dddr
   use radial_der, only: get_dr, get_ddr, get_dddr, get_dr
   use mem_alloc, only: bytes_allocated
   use useful, only: logWrite
   use parallel_mod, only: rank
   use output_data, only: tag
   use finite_differences, only: get_FD_grid, type_stencil

   !-- to be removed
   use integration, only: rInt_R
 
   implicit none

   private
 
   !-- arrays depending on r:
   real(cp), public, allocatable :: r(:)         ! radii
   real(cp), public, allocatable :: r_ic(:)      ! IC radii
   real(cp), public, allocatable :: O_r_ic(:)    ! Inverse of IC radii
   real(cp), public, allocatable :: O_r_ic2(:)   ! Inverse of square of IC radii
   real(cp), public, allocatable :: or1(:)       ! :math:`1/r`
   real(cp), public, allocatable :: or2(:)       ! :math:`1/r^2`
   real(cp), public, allocatable :: or3(:)       ! :math:`1/r^3`
   real(cp), public, allocatable :: or4(:)       ! :math:`1/r^4`
   real(cp), public, allocatable :: otemp1(:)    ! Inverse of background temperature
   real(cp), public, allocatable :: rho0(:)      ! Inverse of background density
   real(cp), public, allocatable :: temp0(:)     ! Background temperature
   real(cp), public, allocatable :: dLtemp0(:)   ! Inverse of temperature scale height
   real(cp), public, allocatable :: ddLtemp0(:)  ! :math:`d/dr(1/T dT/dr)` 
   real(cp), private, allocatable :: d2temp0(:)  ! Second rad. derivative of background temperature
   real(cp), public, allocatable :: dentropy0(:) ! Radial gradient of background entropy
   real(cp), public, allocatable :: orho1(:)     ! :math:`1/\tilde{\rho}`
   real(cp), public, allocatable :: orho2(:)     ! :math:`1/\tilde{\rho}^2`
   real(cp), public, allocatable :: beta(:)      ! Inverse of density scale height drho0/rho0
   real(cp), public, allocatable :: dbeta(:)     ! Radial gradient of beta

   real(cp), public, allocatable :: alpha0(:)    ! Thermal expansion coefficient
   real(cp), public, allocatable :: dLalpha0(:)  ! :math:`1/\alpha d\alpha/dr`
   real(cp), public, allocatable :: ddLalpha0(:) ! :math:`d/dr(1/alpha d\alpha/dr)`
   real(cp), public, allocatable :: ogrun(:)     ! :math:`1/\Gamma`

   real(cp), public, allocatable :: drx(:)       ! First derivative of non-linear mapping (see Bayliss and Turkel, 1990)
   real(cp), public, allocatable :: ddrx(:)      ! Second derivative of non-linear mapping
   real(cp), public, allocatable :: dddrx(:)     ! Third derivative of non-linear mapping
   real(cp), public :: dr_fac                    ! :math:`2/d`, where :math:`d=r_o-r_i`
   real(cp), public :: dr_fac_ic                 ! For IC: :math:`2/(2 r_i)`
   real(cp), public :: alpha1   ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alpha2   ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
   real(cp), public :: r_cmb                     ! OC radius
   real(cp), public :: r_icb                     ! IC radius
   real(cp), public :: r_surface                 ! Surface radius for extrapolation
 
   !-- arrays for buoyancy, depend on Ra and Pr:
   real(cp), public, allocatable :: rgrav(:)     ! Buoyancy term `dtemp0/Di`
 
   !-- chebychev polynomials, derivatives and integral:
   real(cp), public :: cheb_norm                    ! Chebyshev normalisation 
   real(cp), public, allocatable :: cheb(:,:)       ! Chebyshev polynomials
   real(cp), public, allocatable :: dcheb(:,:)      ! First radial derivative
   real(cp), public, allocatable :: d2cheb(:,:)     ! Second radial derivative
   real(cp), public, allocatable :: d3cheb(:,:)     ! Third radial derivative
   real(cp), public, allocatable :: cheb_int(:)     ! Array for cheb integrals
   integer, public :: nDi_costf1                     ! Radii for transform
   integer, public :: nDd_costf1                     ! Radii for transform
   type(costf_odd_t), public :: chebt_oc
   type(costf_odd_t), public :: chebt_ic
   type(costf_even_t), public :: chebt_ic_even

   real(cp), public, allocatable :: eye(:,:)    ! Identity matrix
   real(cp), public, allocatable :: drMat(:,:)  ! First radial derivative
   real(cp), public, allocatable :: d2rMat(:,:) ! Second radial derivative
   real(cp), public, allocatable :: d3rMat(:,:) ! Third radial derivative
   type(type_stencil), public :: stencil_oc

   class(type_rscheme), public, pointer :: rscheme_oc
 
   !-- same for inner core:
   real(cp), public :: cheb_norm_ic                      ! Chebyshev normalisation for IC
   real(cp), public, allocatable :: cheb_ic(:,:)         ! Chebyshev polynomials for IC
   real(cp), public, allocatable :: dcheb_ic(:,:)        ! First radial derivative of cheb_ic
   real(cp), public, allocatable :: d2cheb_ic(:,:)       ! Second radial derivative cheb_ic
   real(cp), public, allocatable :: cheb_int_ic(:)       ! Array for integrals of cheb for IC
   integer, public :: nDi_costf1_ic                      ! Radii for transform
 
   integer, public :: nDd_costf1_ic                      ! Radii for transform
 
   integer, public :: nDi_costf2_ic                      ! Radii for transform
 
   integer, public :: nDd_costf2_ic                      ! Radii for transform
 
   !-- Radius functions for cut-back grid without boundaries:
   !-- (and for the nonlinear mapping)
   real(cp), public :: alph1       ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alph2       ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
 
   real(cp), public, allocatable :: lambda(:)     ! Array of magnetic diffusivity
   real(cp), public, allocatable :: dLlambda(:)   ! Derivative of magnetic diffusivity
   real(cp), public, allocatable :: jVarCon(:)    ! Analytical solution for toroidal field potential aj (see init_fields.f90)
   real(cp), public, allocatable :: sigma(:)      ! Electrical conductivity
   real(cp), public, allocatable :: kappa(:)      ! Thermal diffusivity
   real(cp), public, allocatable :: dLkappa(:)    ! Derivative of thermal diffusivity
   real(cp), public, allocatable :: visc(:)       ! Kinematic viscosity
   real(cp), public, allocatable :: dLvisc(:)     ! Derivative of kinematic viscosity
   real(cp), public, allocatable :: divKtemp0(:)  ! Term for liquid anelastic approximation
   real(cp), public, allocatable :: epscProf(:)   ! Sources in heat equations

   public :: initialize_radial_functions, radial, transportProperties, &
   &         finalize_radial_functions

contains

   subroutine initialize_radial_functions
      !
      ! Initial memory allocation
      !

      integer :: n_in

      ! allocate the arrays
      allocate( r(n_r_max) )
      allocate( r_ic(n_r_ic_max) )
      allocate( O_r_ic(n_r_ic_max) )
      allocate( O_r_ic2(n_r_ic_max) )
      allocate( or1(n_r_max),or2(n_r_max),or3(n_r_max),or4(n_r_max) )
      allocate( otemp1(n_r_max),rho0(n_r_max),temp0(n_r_max) )
      allocate( dLtemp0(n_r_max),d2temp0(n_r_max),dentropy0(n_r_max) )
      allocate( ddLtemp0(n_r_max) )
      allocate( orho1(n_r_max),orho2(n_r_max) )
      allocate( beta(n_r_max), dbeta(n_r_max) )
      allocate( alpha0(n_r_max), dLalpha0(n_r_max), ddLalpha0(n_r_max) )
      allocate( drx(n_r_max),ddrx(n_r_max),dddrx(n_r_max) )
      allocate( rgrav(n_r_max), ogrun(n_r_max) )
      bytes_allocated = bytes_allocated + &
                        (24*n_r_max+3*n_r_ic_max)*SIZEOF_DEF_REAL

      allocate( lambda(n_r_max),dLlambda(n_r_max),jVarCon(n_r_max) )
      allocate( sigma(n_r_max) )
      allocate( kappa(n_r_max),dLkappa(n_r_max) )
      allocate( visc(n_r_max),dLvisc(n_r_max) )
      allocate( epscProf(n_r_max),divKtemp0(n_r_max) )
      bytes_allocated = bytes_allocated + 10*n_r_max*SIZEOF_DEF_REAL


      if ( .not. l_finite_diff ) then

         allocate( cheb(n_r_max,n_r_max) )     ! Chebychev polynomials
         allocate( dcheb(n_r_max,n_r_max) )    ! first radial derivative
         allocate( d2cheb(n_r_max,n_r_max) )   ! second radial derivative
         allocate( d3cheb(n_r_max,n_r_max) )   ! third radial derivative
         allocate( cheb_int(n_r_max) )         ! array for cheb integrals !
         bytes_allocated = bytes_allocated + &
                           (4*n_r_max*n_r_max+n_r_max)*SIZEOF_DEF_REAL

         nDi_costf1=2*n_r_max+2
         nDd_costf1=2*n_r_max+5
         call chebt_oc%initialize(n_r_max,nDi_costf1,nDd_costf1)

         allocate ( type_cheb_odd :: rscheme_oc )

         nDi_costf1_ic=2*n_r_ic_max+2
         nDd_costf1_ic=2*n_r_ic_max+5
         nDi_costf2_ic=2*n_r_ic_max
         nDd_costf2_ic=2*n_r_ic_max+n_r_ic_max/2+5

         allocate( cheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( dcheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( d2cheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( cheb_int_ic(n_r_ic_max) )
         bytes_allocated = bytes_allocated + &
                           (3*n_r_ic_max*n_r_ic_max+n_r_ic_max)*SIZEOF_DEF_REAL

         call chebt_ic%initialize(n_r_ic_max,nDi_costf1_ic,nDd_costf1_ic)

         n_in = n_cheb_max

      else

         allocate ( type_fd :: rscheme_oc )

         !-- To be removed
         allocate( eye(n_r_max,n_r_max) )      ! identity matrix
         allocate( drMat(n_r_max,n_r_max) )    ! first radial derivative
         allocate( d2rMat(n_r_max,n_r_max) )   ! second radial derivative
         allocate( d3rMat(n_r_max,n_r_max) )   ! second radial derivative
         bytes_allocated = bytes_allocated + &
                           (4*n_r_max*n_r_max)*SIZEOF_DEF_REAL
         call stencil_oc%initialize(n_r_max,fd_order)
         !--

         n_in = fd_order

      end if
      call rscheme_oc%initialize(n_r_max,n_in)


   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine finalize_radial_functions

      deallocate( r, r_ic, O_r_ic, O_r_ic2, or1, or2, or3, or4 )
      deallocate( otemp1, rho0, temp0, dLtemp0, d2temp0, dentropy0 )
      deallocate( ddLtemp0, orho1, orho2, beta, dbeta, alpha0 )
      deallocate( ddLalpha0, dLalpha0, drx, ddrx, dddrx, rgrav, ogrun )
      deallocate( lambda, dLlambda, jVarCon, sigma, kappa, dLkappa )
      deallocate( visc, dLvisc, epscProf, divKtemp0 )

      if ( .not. l_finite_diff ) then
         deallocate( cheb, dcheb, d2cheb, d3cheb, cheb_int )
         deallocate( cheb_ic, dcheb_ic, d2cheb_ic, cheb_int_ic )
         call chebt_oc%finalize()
         call chebt_ic%finalize()
         if ( n_r_ic_max > 0 .and. l_cond_ic ) call chebt_ic_even%finalize()
      else
         deallocate( eye, drMat)
         call stencil_oc%finalize()
      end if

      call rscheme_oc%finalize()

   end subroutine finalize_radial_functions
!------------------------------------------------------------------------------
   subroutine radial
      !
      !  Calculates everything needed for radial functions, transforms etc.
      !

      !-- Local variables:
      integer :: n_r,n_cheb,n_cheb_int
      integer :: n_r_ic_tot, k, i
      integer :: n_const(1)

      !integer :: n_r_start
      real(cp) :: fac_int
      real(cp) :: r_cheb(n_r_max)
      real(cp) :: r_cheb_ic(2*n_r_ic_max-1),r_ic_2(2*n_r_ic_max-1)
      real(cp) :: drho0(n_r_max),dtemp0(n_r_max)
      real(cp) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

      real(cp) :: hcomp,fac
      real(cp) :: dtemp0cond(n_r_max),dtemp0ad(n_r_max),hcond(n_r_max)
      real(cp) :: func(n_r_max)

      real(cp), allocatable :: coeffDens(:), coeffTemp(:), coeffAlpha(:)
      real(cp) :: w1(n_r_max), w2(n_r_max), rrOcmb(n_r_max)
      character(len=80) :: message
      character(len=76) :: fileName
      integer :: fileHandle

      real(cp) :: ratio1, ratio2

     ! To be removed:
      real(cp) :: f(n_r_max), df(n_r_max), ddf(n_r_max), dddf(n_r_max)
      complex(cp) :: f_c(1,n_r_max), df_c(1,n_r_max), ddf_c(1,n_r_max), dddf_c(1,n_r_max)
      integer :: j

      !-- Radial grid point:
      !   radratio is aspect ratio
      !   radratio = (inner core r) / (CMB r) = r_icb/r_cmb
      r_cmb=one/(one-radratio)
      r_icb=r_cmb-one
      r_surface=2.8209_cp    ! in units of (r_cmb-r_icb)


      if ( .not. l_finite_diff ) then

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
         !   r_cheb are the grid_points in the Cheb interval [-1,1]
         !   and r are these points mapped to the interval [r_icb,r_cmb]:
         call cheb_grid(r_icb,r_cmb,n_r_max-1,r,r_cheb, &
              &               alpha1,alpha2,paraX0,lambd)

         if ( l_newmap ) then
            do n_r=1,n_r_max
               drx(n_r) =                          (two*alpha1) /      &
                    ((one+alpha1**2*(two*r(n_r)-r_icb-r_cmb-alpha2)**2)* &
                    lambd)
               ddrx(n_r) = -(8.0_cp*alpha1**3*(two*r(n_r)-r_icb-r_cmb-alpha2)) / &
                    ((one+alpha1**2*(-two*r(n_r)+r_icb+r_cmb+alpha2)**2)**2*  & 
                    lambd)
               dddrx(n_r) =        (16.0_cp*alpha1**3*(-one+three*alpha1**2* &
                                     (-two*r(n_r)+r_icb+r_cmb+alpha2)**2)) / &
                    ((one+alpha1**2*(-two*r(n_r)+r_icb+r_cmb+alpha2)**2)**3* &
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
              &         cheb,dcheb,d2cheb,d3cheb,n_r_max,n_r_max, &
              &         drx,ddrx,dddrx)
      else
         call get_FD_grid(fd_stretch, fd_ratio, r_icb, r_cmb, r)

         call stencil_oc%get_FD_coeffs(r)

         call stencil_oc%nullify_epsilon()

         call stencil_oc%get_FD_matder(n_r_max,eye,drMat,d2rMat,d3rMat)


      end if

      if ( .not. l_finite_diff ) then
         ratio1=alph1
         ratio2=alph2
      else
         ratio1=fd_stretch
         ratio2=fd_ratio
      end if
      call rscheme_oc%get_grid(n_r_max, r_icb, r_cmb, ratio1, ratio2, r)

      call rscheme_oc%get_der_mat(n_r_max)

      if ( l_finite_diff ) then
         !-- To be removed
         do n_r=1,n_r_max
            f_c(1,n_r)=sin(r(n_r))
         end do

         call get_dr(f_c, df_c, 1, 1, 1, n_r_max, rscheme_oc)
         call get_ddr(f_c, df_c, ddf_c, 1, 1, 1, n_r_max, rscheme_oc)
         call get_dddr(f_c, df_c, ddf_c, dddf_c, 1, 1, 1, n_r_max, rscheme_oc)

         open(newunit=fileHandle, file='test.txt', status='unknown')

         do n_r=1,n_r_max
         write(fileHandle, '(5ES20.12)') r(n_r), real(f_c(1,n_r)), &
         &                               real(df_c(1,n_r)), real(ddf_c(1,n_r)), &
         &                               real(dddf_c(1,n_r))
         end do

         f(:)=real(f_c(1,:))

         print*, rInt_R(f,r,rscheme_oc)
         
         do i=1,n_r_max
            df(i)=0.0_cp
            ddf(i)=0.0_cp
            dddf(i)=0.0_cp
            do j=1,n_r_max
               df(i) = df(i) + drMat(i,j)*f(j)
               ddf(i) = ddf(i) + d2rMat(i,j)*f(j)
               dddf(i) = dddf(i) + d3rMat(i,j)*f(j)
            end do
         end do

         do n_r=1,n_r_max
            ! print*, r(n_r),real(df_c(1,n_r)),df(n_r)
            ! print*, r(n_r),real(ddf_c(1,n_r)),ddf(n_r)
            print*, r(n_r),real(dddf_c(1,n_r)),dddf(n_r)
         end do

         close(fileHandle)
         ! stop
      end if

      if ( rank == 0 ) then
         fileName = 'radius.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         do n_r=1,n_r_max
            write(fileHandle,'(I4, ES16.8)') n_r, r(n_r)
         end do
         close(fileHandle)
      end if
      !stop

      or1=one/r         ! 1/r
      or2=or1*or1       ! 1/r**2
      or3=or1*or2       ! 1/r**3
      or4=or2*or2       ! 1/r**4

      !-- Get entropy gradient
      call getEntropyGradient ! By default this is zero

      !-- Fit to an interior model
      if ( index(interior_model,'JUP') /= 0 ) then

         if ( l_non_adia ) then
            rrOcmb(:) = r(:)*r_cut_model/r_cmb
            rgrav(:)=  83.4792166296_cp*rrOcmb(:)+7.0970761715_cp*rrOcmb(:)**2 &
            &        -112.0000517766_cp*rrOcmb(:)**3+47.3447404648_cp*rrOcmb(:)**4

            allocate ( coeffAlpha(10), coeffTemp(10) )
            coeffAlpha = [ -12.9483344953_cp, 12.7631620079_cp, -60.0717008192_cp, &
               &           1.41916870466_cp, 755.055391736_cp, -1938.08838168_cp,  &
               &           952.893688457_cp, 2544.71502695_cp, -3703.20551213_cp,  &
               &           1440.95591192_cp]
            coeffTemp = [ 1.24462655e+05_cp, -2.85767595e+06_cp, 3.04794799e+07_cp, &
               &         -1.72807386e+08_cp, 5.83323621e+08_cp, -1.23322830e+09_cp, &
               &          1.64950647e+09_cp, -1.35626053e+09_cp, 6.25695974e+08_cp, &
               &          -1.23976043e+08_cp ]
            ! Temperature is only required to estimate ThExpNb
            alpha0(:)=0.0_cp
            temp0(:) =0.0_cp
            do i=1,10
               alpha0(:) = alpha0(:)+coeffAlpha(i)*rrOcmb(:)**(i-1)
               temp0(:)  = temp0(:) +coeffTemp(i) *rrOcmb(:)**(i-1)
            end do
            alpha0(:)=exp(alpha0(:)) ! Polynomial fit was on ln(alpha0)
            DissNb   =alpha0(1)*rgrav(1)*(rrOcmb(1)-rrOcmb(n_r_max))*6.9894e7 &
            &         /1.5e4_cp ! 1.5e4 is cp 6.9894e7 is R_J
            ThExpNb  =alpha0(1)*temp0(1)
            alpha0(:)=alpha0(:)/alpha0(1)
            rgrav(:) =rgrav(:)/rgrav(1)

            ! d ln(temp0) / dr
            dtemp0(:)=epsS*dentropy0(:)-DissNb*alpha0(:)*rgrav(:)
            call getBackground(dtemp0,0.0_cp,temp0)
            temp0=exp(temp0) ! this was ln(T_0)
            dtemp0=dtemp0*temp0

            !-- Radial profile for the Grüneisen parameter (from French et al.)
            ogrun(:) = one/(0.57_cp-0.17_cp*tanh(50.0_cp*(rrOcmb(:)-0.88_cp)))
            GrunNb = one/ogrun(1)
            ogrun(:) = ogrun(:)/ogrun(1)

            drho0=-ThExpNb*epsS*alpha0*temp0*dentropy0-DissNb/GrunNb*ogrun*alpha0*rgrav
            call getBackground(drho0,0.0_cp,rho0)
            rho0=exp(rho0) ! this was ln(rho_0)
            beta=drho0

            ! The final stuff is always required
            call get_dr(beta,dbeta,n_r_max,rscheme_oc)
            call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
            call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
            dLalpha0=dLalpha0/alpha0 ! d log (alpha) / dr
            call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
            dLtemp0 = dtemp0/temp0
            ddLtemp0 =-(dtemp0/temp0)**2+d2temp0/temp0

            !-- Multiply the gravity by alpha0 and temp0
            rgrav(:)=rgrav(:)*alpha0(:)*temp0(:)

            !nVarDiff = 5
            deallocate(coeffAlpha, coeffTemp)
         else
            allocate( coeffDens(8), coeffTemp(10) )
            coeffDens = [4.46020423_cp, -4.60312999_cp, 37.38863965_cp,       &
               &         -201.96655354_cp, 491.00495215_cp, -644.82401602_cp, &
               &         440.86067831_cp, -122.36071577_cp] 

            coeffTemp = [0.999735638_cp, 0.0111053831_cp, 2.70889691_cp,  &
               &         -83.5604443_cp, 573.151526_cp, -1959.41844_cp,   &
               &         3774.39367_cp, -4159.56327_cp, 2447.75300_cp,    &
               &         -596.464198_cp]

            call polynomialBackground(coeffDens,coeffTemp)
            deallocate( coeffDens, coeffTemp)
         end if

      else if ( index(interior_model,'SAT') /= 0 ) then

         ! the shell can't be thicker than eta=0.15, because the fit doesn't 
         ! work below that (in Nadine's profile, that's where the IC is anyway)
         ! also r_cut_model maximum is 0.999, because rho is negative beyond
         ! that

         allocate( coeffDens(4), coeffTemp(9) )
         coeffDens = [-0.33233543_cp, 0.90904075_cp, -0.9265371_cp, &
            &         0.34973134_cp ]

         coeffTemp = [1.00294605_cp,-0.44357815_cp,13.9295826_cp,  &
            &         -137.051347_cp,521.181670_cp,-1044.41528_cp, &
            &         1166.04926_cp,-683.198387_cp, 162.962632_cp ]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'SUN') /= 0 ) then

         ! rho is negative beyond r_cut_model=0.9965
         ! radratio should be 0.7 (size of the Sun's CZ)

         allocate( coeffDens(6), coeffTemp(4) )
         coeffDens = [-24.83750402_cp, 231.79029994_cp, -681.72774358_cp, &
            &         918.30741266_cp,-594.30093367_cp, 150.76802942_cp ]

         coeffTemp = [5.53715416_cp, -8.10611274_cp, 1.7350452_cp, &
            &         0.83470843_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'GLIESE229B') /= 0 ) then
         ! Use also nVarDiff=2 with difExp=0.52

         allocate( coeffDens(8), coeffTemp(5) )
         coeffDens = [0.99879163_cp,0.15074601_cp,-4.20328423_cp,   &
            &         6.43542034_cp,-12.67297113_cp,21.68593078_cp, &
            &         -17.74832309_cp,5.35405134_cp]

         coeffTemp = [0.99784506_cp,0.16540448_cp,-3.44594354_cp,   &
            &         3.68189750_cp,-1.39046384_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'COROT3B') /= 0 ) then
         ! Use also nVarDiff=2 with difExp=0.62

         allocate( coeffDens(7), coeffTemp(7) )
         coeffDens = [1.00035987_cp,-0.01294658_cp,-2.78586315_cp,  &
            &         0.70289860_cp,2.59463562_cp,-1.65868190_cp,   &
            &         0.15984718_cp]

         coeffTemp = [1.00299303_cp,-0.33722671_cp,1.71340063_cp,     & 
            &         -12.50287121_cp,21.52708693_cp,-14.91959338_cp, &
            &         3.52970611_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'KOI889B') /= 0 ) then
         ! Use also nVarDiff=2 with difExp=0.68

         allocate( coeffDens(6), coeffTemp(6) )
         coeffDens = [1.01038678_cp,-0.17615484_cp,-1.50567127_cp,  &
            &         -1.65738032_cp,4.20394427_cp,-1.87394994_cp]

         coeffTemp = [1.02100249_cp,-0.60750867_cp,3.23371939_cp,   &
            &         -12.80774142_cp,15.37629271_cp,-6.19288785_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'EARTH') /= 0 ) then
         DissNb =0.3929_cp ! Di = \alpha_O g d / c_p
         ThExpNb=0.0566_cp ! Co = \alpha_O T_O
         GrunNb =1.5_cp ! Gruneisen paramater
         hcomp  =2.2_cp*r_cmb

         alpha0=(one+0.6_cp*r**2/hcomp**2)/(one+0.6_cp/2.2_cp**2)
         rgrav =(r-0.6_cp*r**3/hcomp**2)/(r_cmb*(one-0.6_cp/2.2_cp**2))

         !dentropy0 = -half*(ampStrat+one)*(one-tanh(slopeStrat*(r-rStrat)))+ &
         !            & ampStrat

         !! d ln(temp0) / dr
         !dtemp0=epsS*dentropy0-DissNb*alpha0*rgrav

         !call getBackground(dtemp0,0.0_cp,temp0)
         !temp0=exp(temp0) ! this was ln(T_0)
         !dtemp0=dtemp0*temp0

         !drho0=-ThExpNb*epsS*alpha0*temp0*dentropy0-DissNb/GrunNb*alpha0*rgrav
         !call getBackground(drho0,0.0_cp,rho0)
         !rho0=exp(rho0) ! this was ln(rho_0)
         !beta=drho0

         hcond = (one-0.4469_cp*(r/r_cmb)**2)/(one-0.4469_cp)
         hcond = hcond/hcond(1)
         temp0 = (one+GrunNb*(r_icb**2-r**2)/hcomp**2)
         temp0 = temp0/temp0(1)
         dtemp0cond=-cmbHflux/(r**2*hcond)
          
         do k=1,10 ! 10 iterations is enough to converge
            dtemp0ad=-DissNb*alpha0*rgrav*temp0-epsS*temp0(n_r_max)
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

         dentropy0=dtemp0/temp0/epsS+DissNb*alpha0*rgrav/epsS
         !drho0=-ThExpNb*epsS*alpha0*temp0*dentropy0-DissNb*alpha0*rgrav/GrunNb
         drho0=-DissNb*alpha0*rgrav/GrunNb
         call getBackground(drho0,0.0_cp,rho0)
         rho0=exp(rho0) ! this was ln(rho_0)
         beta=drho0

         ! The final stuff is always required
         call get_dr(beta,dbeta,n_r_max,rscheme_oc)
         call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
         call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
         dLalpha0=dLalpha0/alpha0 ! d log (alpha) / dr
         call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
         dLtemp0 = dtemp0/temp0
         call get_dr(dLtemp0,ddLtemp0,n_r_max,rscheme_oc)

         ! N.B. rgrav is not gravity but alpha * grav
         rgrav = alpha0*rgrav

         l_non_adia = .true.

      else  !-- Usual polytropic reference state
         ! g(r) = g0 + g1*r/ro + g2*(ro/r)**2
         ! Default values: g0=0, g1=1, g2=0
         ! An easy way to change gravity
         rgrav(:)=g0+g1*r(:)/r_cmb+g2*(r_cmb/r)**2

         if ( l_anel ) then

            if ( l_non_adia ) then

               dtemp0(:)=-DissNb*rgrav(:)
               !-- Use d2temp0 as a work array
               d2temp0(:)=-epsS*dentropy0(:)
               call getBackground(dtemp0,1.0_cp,temp0,d2temp0)
               dtemp0(:)=epsS*temp0(:)*dentropy0(:)-DissNb*rgrav(:)

               drho0=-ThExpNb*epsS*dentropy0-DissNb/GrunNb*rgrav/temp0
               call getBackground(drho0,0.0_cp,rho0)
               rho0=exp(rho0) ! this was ln(rho_0)
               beta=drho0

               ! The final stuff is always required
               call get_dr(beta,dbeta,n_r_max,rscheme_oc)
               call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)

               dLtemp0 = dtemp0/temp0
               ddLtemp0 =-(dtemp0/temp0)**2+d2temp0/temp0

            else !-- Adiabatic reference state

               if ( l_isothermal ) then ! Gruneisen is zero in this limit
                  fac      =strat /( g0+half*g1*(one+radratio) +g2/radratio )
                  DissNb   =0.0_cp
                  GrunNb   =0.0_cp
                  temp0    =one
                  rho0     =exp(-fac*(g0*(r-r_cmb) +      &
                            g1/(two*r_cmb)*(r**2-r_cmb**2) - &
                            g2*(r_cmb**2/r-r_cmb)))

                  beta     =-fac*rgrav
                  dbeta    =-fac*(g1/r_cmb-two*g2*r_cmb**2*or3)
                  d2temp0  =0.0_cp
                  dLtemp0  =0.0_cp
                  ddLtemp0 =0.0_cp
               else
                  if ( strat == 0.0_cp .and. DissNb /= 0.0_cp ) then
                     strat = polind* log(( g0+half*g1*(one+radratio)+g2/radratio )* &
                                         DissNb+1)
                  else
                     DissNb=( exp(strat/polind)-one )/ &
                            ( g0+half*g1*(one+radratio) +g2/radratio )
                  end if
                  GrunNb   =one/polind
                  temp0    =-DissNb*( g0*r+half*g1*r**2/r_cmb-g2*r_cmb**2/r ) + &
                            one + DissNb*r_cmb*(g0+half*g1-g2)
                  rho0     =temp0**polind

                  !-- Computation of beta= dln rho0 /dr and dbeta=dbeta/dr
                  beta     =-polind*DissNb*rgrav/temp0
                  dbeta    =-polind*DissNb/temp0**2 *         &
                            ((g1/r_cmb-two*g2*r_cmb**2*or3)*  &
                            temp0  + DissNb*rgrav**2)
                  dtemp0   =-DissNb*rgrav
                  d2temp0  =-DissNb*(g1/r_cmb-two*g2*r_cmb**2*or3)
                  dLtemp0  =dtemp0/temp0
                  ddLtemp0 =-(dtemp0/temp0)**2+d2temp0/temp0
               end if
            end if

            !-- Thermal expansion coefficient (1/T for an ideal gas)
            alpha0   =one/temp0
            ogrun    =one
            dLalpha0 =-dLtemp0
            ddLalpha0=-ddLtemp0

         end if
      end if

      if ( l_anel ) then
         call logWrite('')
         call logWrite('!      This is an anelastic model')
         if ( l_TP_form ) then
            call logWrite('! You use temperature and pressure as thermodynamic variables')
         else
            call logWrite('! You use entropy and pressure as thermodynamic variables')
         end if
         if ( l_temperature_diff ) then
            call logWrite('! You use temperature diffusion')
         else
            call logWrite('! You use entropy diffusion')
         end if
         if ( l_single_matrix ) then
            call logWrite('! You solve one big matrix (3*n_r_max,3*n_r_max)')
         else
            call logWrite('! You solve two small matrices')
         end if
         call logWrite('! The key parameters are the following')
         write(message,'(''!      DissNb ='',ES16.6)') DissNb
         call logWrite(message)
         write(message,'(''!      ThExpNb='',ES16.6)') ThExpNb
         call logWrite(message)
         write(message,'(''!      GrunNb ='',ES16.6)') GrunNb
         call logWrite(message)
         write(message,'(''!      N_rho  ='',ES16.6)') strat
         call logWrite(message)
         write(message,'(''!      pol_ind='',ES16.6)') polind
         call logWrite(message)
         call logWrite('')
      end if

      !-- Get additional functions of r:
      if ( l_anel ) then
         orho1      =one/rho0
         orho2      =orho1*orho1
         otemp1     =one/temp0
         ViscHeatFac=DissNb*pr/raScaled
         if (l_mag) then
            OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
         else
            OhmLossFac=0.0_cp
         end if
      else
         rho0     =one
         temp0    =one
         otemp1   =one
         orho1    =one
         orho2    =one
         alpha0   =one
         ogrun    =one
         beta     =0.0_cp
         dbeta    =0.0_cp
         dLalpha0 =0.0_cp
         ddLalpha0=0.0_cp
         dLtemp0  =0.0_cp
         ddLtemp0 =0.0_cp
         d2temp0  =0.0_cp
         ViscHeatFac=0.0_cp
         OhmLossFac =0.0_cp
      end if

      !-- Factors for cheb integrals:
      if ( .not. l_finite_diff ) then
         cheb_int(1)=one   ! Integration constant chosen !
         do n_cheb=3,n_r_max,2
            cheb_int(n_cheb)  =-one/real(n_cheb*(n_cheb-2),kind=cp)
            cheb_int(n_cheb-1)= 0.0_cp
         end do
      end if

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
         call chebt_ic_even%initialize(n_r_ic_max-1,nDi_costf2_ic,nDd_costf2_ic)

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
      !
      ! Calculates the transport properties: electrical conductivity,
      ! kinematic viscosity and thermal conductivity.
      !

      integer :: n_r

      real(cp) :: a,b,c,s1,s2,r0
      real(cp) :: dsigma0
      real(cp) :: dvisc(n_r_max), dkappa(n_r_max), dsigma(n_r_max)
      !real(cp) :: condBot(n_r_max), condTop(n_r_max)
      !real(cp) :: func(n_r_max)
      real(cp) :: kcond(n_r_max)
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
                   exit
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
             call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
             dLlambda=dsigma/lambda
          else if ( nVarCond == 4 ) then ! Profile
             lambda=(rho0/rho0(n_r_max))**difExp
             sigma=one/lambda
             call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
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
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 2 ) then ! Profile
            kappa=(rho0/rho0(n_r_max))**difExp
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 3 ) then ! polynomial fit to a model
            if ( radratio < 0.19_cp ) then
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
               &              + a2*rrOcmb**3 + a1*rrOcmb**4 &
               &                             + a0*rrOcmb**5

            end do
            kappatop=kappa(1) ! normalise by the value at the top
            kappa=kappa/kappatop
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
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
            !call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            !dLkappa=dkappa/kappa

            ! Alternative scenario
            kcond=(one-0.4469_cp*(r/r_cmb)**2)/(one-0.4469_cp)
            kcond=kcond/kcond(1)
            kappa=kcond/rho0
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 5 ) then ! Thermal background equilibrium
            kcond(:)=one/(r(:)*r(:)*dLtemp0(:)*temp0(:))
            kcond=kcond/kcond(1)
            kappa=kcond/rho0
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         end if
      end if

      !-- Eps profiles
      !-- The remaining division by rho will happen in updateS.f90
      if ( nVarEps == 0 ) then
         ! eps is constant
         if ( l_anelastic_liquid .or. l_TP_form ) then
            epscProf(:)=one
         else
            epscProf(:)=otemp1(:)
         end if
      else if ( nVarEps == 1 ) then
         ! rho*eps in the RHS
         if ( l_anelastic_liquid .or. l_TP_form ) then
            epscProf(:)=rho0(:)
         else
            epscProf(:)=rho0(:)*otemp1(:)
         end if
      end if

      !-- Variable viscosity
      if ( nVarVisc == 0 ) then ! default: constant kinematic viscosity
         visc  =one
         dLvisc=0.0_cp
      else if ( nVarVisc == 1 ) then ! Constant dynamic viscosity
         visc=rho0(n_r_max)/rho0
         call get_dr(visc,dvisc,n_r_max,rscheme_oc)
         dLvisc=dvisc/visc
      else if ( nVarVisc == 2 ) then ! Profile
         visc=(rho0/rho0(n_r_max))**difExp
         call get_dr(visc,dvisc,n_r_max,rscheme_oc)
         dLvisc=dvisc/visc
      end if

      if ( l_anelastic_liquid .or. l_non_adia ) then
         divKtemp0=rho0*kappa*(d2temp0+(beta+dLkappa+two*or1)*temp0*dLtemp0)*sq4pi
      else
         divKtemp0=0.0_cp
      end if

   end subroutine transportProperties
!------------------------------------------------------------------------------
   subroutine getEntropyGradient
      !
      ! This subroutine allows to calculate the background entropy gradient
      ! in case stable stratification is required
      !

      integer :: n_r

      if ( nVarEntropyGrad == 0 ) then ! Default: isentropic
         dEntropy0(:)=0.0_cp
         l_non_adia = .false.
      else if ( nVarEntropyGrad == 1 ) then ! Takehiro
         if ( rStrat <= r_icb ) then
            dentropy0(:) = ampStrat
         else
            dentropy0(:) = -half*(ampStrat+one)*(one-tanh(slopeStrat*(r(:)-rStrat)))&
            &              + ampStrat
         end if
         l_non_adia = .true.
      else if ( nVarEntropyGrad == 2 ) then ! Flat + linear
         if ( rStrat <= r_icb ) then
            dentropy0(:) = ampStrat
         else
            do n_r=1,n_r_max
               if ( r(n_r) <= rStrat ) then
                  dentropy0(n_r)=-one
               else
                  dentropy0(n_r)=(ampStrat+one)*(r(n_r)-r_cmb)/(r_cmb-rStrat) + &
                  &              ampStrat
               end if
            end do
         end if
         l_non_adia = .true.
      else if ( nVarEntropyGrad == 3 ) then ! SSL
         if ( rStrat <= r_icb ) then
            dentropy0(:) = -one
         else
            dentropy0(:) = 0.25_cp*(ampStrat+one)*(one+tanh(slopeStrat*(r(:)-rStrat)))&
            &              *(one-tanh(slopeStrat*(r(:)-rStrat-thickStrat)))           &
            &              - one
         end if
         l_non_adia = .true.
      end if

   end subroutine getEntropyGradient
!------------------------------------------------------------------------------
   subroutine getBackground(input,boundaryVal,output,coeff)
      !
      ! Linear solver of the form: df/dx +coeff*f= input with f(1)=boundaryVal
      ! 

      !-- Input variables:
      real(cp), intent(in) :: input(n_r_max)
      real(cp), intent(in) :: boundaryVal
      real(cp), optional, intent(in) :: coeff(n_r_max)

      !-- Output variables:
      real(cp), intent(out) :: output(n_r_max)

      !-- Local variables:
      real(cp) :: rhs(n_r_max)
      integer :: n_r,info,n_r_out, n_cheb
      real(cp) :: workMat_fac(n_r_max, 2)

      real(cp), allocatable :: workMat(:,:)
      integer, allocatable :: workPivot(:)

      allocate( workMat(n_r_max,n_r_max) )
      allocate( workPivot(n_r_max) )

      do n_r_out=1,n_r_max
         do n_r=2,n_r_max
            if ( present(coeff) ) then
               workMat(n_r,n_r_out)=rscheme_oc%rnorm*(              &
               &                      rscheme_oc%drMat(n_r,n_r_out)+&
               &            coeff(n_r)*rscheme_oc%rMat(n_r,n_r_out) )
            else
               workMat(n_r,n_r_out)=rscheme_oc%rnorm*               &
               &                    rscheme_oc%drMat(n_r,n_r_out)
            end if
         end do
      end do

      !-- boundary conditions
      do n_r_out=1,rscheme_oc%n_max
         workMat(1,n_r_out)=rscheme_oc%rnorm*rscheme_oc%rMat(1,n_r_out)
         if ( rscheme_oc%version == 'cheb' ) workMat(n_r_max,n_r_out)=0.0_cp
      end do

      !-- fill with zeros
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            workMat(1,n_r_out)=0.0_cp
         end do
      end if

      !-- renormalize
      do n_r=1,n_r_max
        workMat(n_r,1)      =rscheme_oc%boundary_fac*workMat(n_r,1)
        workMat(n_r,n_r_max)=rscheme_oc%boundary_fac*workMat(n_r,n_r_max)
      end do

      do n_r=1,n_r_max
         workMat_fac(n_r,1)=one/maxval(abs(workMat(n_r,:)))
      end do
      do n_r=1,n_r_max
         workMat(n_r,:)=workMat(n_r,:)*workMat_fac(n_r,1)
      end do
      do n_r=1,n_r_max
         workMat_fac(n_r,2)=one/maxval(abs(workMat(:,n_r)))
      end do
      do n_r=1,n_r_max
         workMat(:,n_r)=workMat(:,n_r)*workMat_fac(n_r,2)
      end do


      call sgefa(workMat,n_r_max,n_r_max,workPivot,info)

      if ( info /= 0 ) then
         write(*,*) '! Singular Matrix in getBackground!'
         stop '20'
      end if

      do n_r=2,n_r_max
         rhs(n_r)=input(n_r)
      end do
      rhs(1)=boundaryVal

      do n_r=1,n_r_max
         rhs(n_r)=rhs(n_r)*workMat_fac(n_r,1)
      end do

      !-- Solve for s0:
      call sgesl(workMat,n_r_max,n_r_max,workPivot,rhs)

      !-- Copy result to s0:
      do n_r=1,n_r_max
         output(n_r)=rhs(n_r)*workMat_fac(n_r,2)
      end do

      !-- Set cheb-modes > n_cheb_max to zero:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            output(n_r_out)=0.0_cp
         end do
      end if

      !-- Transform to radial space:
      call rscheme_oc%costf1(output)

      deallocate( workMat, workPivot )

   end subroutine getBackground
!------------------------------------------------------------------------------
   subroutine polynomialBackground(coeffDens,coeffTemp)
      !
      ! This subroutine allows to calculate a reference state based on an input
      ! polynomial function.
      !

      !-- Input variables
      real(cp), intent(in) :: coeffDens(:)
      real(cp), intent(in) :: coeffTemp(:)

      !-- Local variables
      real(cp) :: rrOcmb(n_r_max),gravFit(n_r_max)
      real(cp) :: drho0(n_r_max),dtemp0(n_r_max)
      real(cp) :: w1(n_r_max),w2(n_r_max)

      integer ::  nDens,nTemp,i

      nDens = size(coeffDens)
      nTemp = size(coeffTemp)
      rrOcmb(:) = r(:)*r_cut_model/r_cmb
      gravFit(:)=four*rrOcmb(:)-three*rrOcmb(:)**2

      ! Set to zero initially
      rho0(:) =0.0_cp
      temp0(:)=0.0_cp

      do i=1,nDens
         rho0(:) = rho0(:)+coeffDens(i)*rrOcmb(:)**(i-1)
      end do

      do i=1,nTemp
         temp0(:) = temp0(:)+coeffTemp(i)*rrOcmb(:)**(i-1)
      end do

      ! Normalise to the outer radius
      temp0(:)   = temp0(:)/temp0(1)
      rho0(:)    = rho0(:)/rho0(1)
      gravFit(:) = gravFit(:)/gravFit(1)

      ! Derivative of the temperature needed to get alpha_T
      call get_dr(temp0,dtemp0,n_r_max,rscheme_oc)

      alpha0(:)=-dtemp0(:)/(gravFit(:)*temp0(:))

      ! Dissipation number
      DissNb=alpha0(1)
      alpha0(:)=alpha0(:)/alpha0(1)

      ogrun(:)=alpha0(:)*temp0(:)

      ! Adiabatic: buoyancy term is linked to the temperature gradient

      !       dT
      !      ---- =  -Di * alpha_T * T * grav
      !       dr
      rgrav(:)=-dtemp0(:)/DissNb

      call get_dr(rho0,drho0,n_r_max,rscheme_oc)
      beta(:)=drho0(:)/rho0(:)
      call get_dr(beta,dbeta,n_r_max,rscheme_oc)
      call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
      call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
      dLalpha0(:)=dLalpha0(:)/alpha0(:) ! d log (alpha) / dr
      call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
      dLtemp0(:)=dtemp0(:)/temp0(:)
      call get_dr(dLtemp0,ddLtemp0,n_r_max,rscheme_oc)
      dentropy0(:)=0.0_cp

   end subroutine polynomialBackground
!------------------------------------------------------------------------------
end module radial_functions
