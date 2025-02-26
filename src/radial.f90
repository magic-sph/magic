module radial_functions
   !
   !  This module initiates all the radial functions (transport properties, density,
   !  temperature, cheb transforms, etc.)
   !

   use iso_fortran_env, only: output_unit
   use truncation, only: n_r_max, n_cheb_max, n_r_ic_max, fd_ratio, rcut_l, &
       &                 fd_stretch, fd_order, fd_order_bound, l_max
   use algebra, only: prepare_mat, solve_mat
   use constants, only: sq4pi, one, two, three, four, half, pi, zero
   use physical_parameters
   use logic, only: l_mag, l_cond_ic, l_heat, l_anelastic_liquid,  &
       &            l_isothermal, l_anel, l_non_adia, l_centrifuge,&
       &            l_temperature_diff, l_single_matrix, l_var_l,  &
       &            l_finite_diff, l_newmap, l_full_sphere,        &
       &            l_chemical_conv, l_drag, l_tidal
   use radial_data, only: nRstart, nRstop
   use chebyshev_polynoms_mod ! Everything is needed
   use cosine_transform_odd
   use cosine_transform_even
   use radial_scheme, only: type_rscheme
   use chebyshev, only: type_cheb_odd
   use finite_differences, only: type_fd
   use radial_der, only: get_dr
   use mem_alloc, only: bytes_allocated
   use useful, only: logWrite, abortRun
   use parallel_mod, only: rank
   use output_data, only: tag
   use num_param, only: alph1, alph2

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
   real(cp), public, allocatable :: dxicond(:)   ! Radial gradient of chemical composition
   real(cp), public, allocatable :: orho1(:)     ! :math:`1/\tilde{\rho}`
   real(cp), public, allocatable :: orho2(:)     ! :math:`1/\tilde{\rho}^2`
   real(cp), public, allocatable :: beta(:)      ! Inverse of density scale height drho0/rho0
   real(cp), public, allocatable :: dbeta(:)     ! Radial gradient of beta
   real(cp), public, allocatable :: ddbeta(:)    ! 2nd derivative of beta

   real(cp), public, allocatable :: alpha0(:)    ! Thermal expansion coefficient
   real(cp), public, allocatable :: dLalpha0(:)  ! :math:`1/\alpha d\alpha/dr`
   real(cp), public, allocatable :: ddLalpha0(:) ! :math:`d/dr(1/alpha d\alpha/dr)`
   real(cp), public, allocatable :: ogrun(:)     ! :math:`1/\Gamma`
   real(cp), public, allocatable :: drag0(:)     ! Radiative drag neutrinos

   real(cp), public :: dr_fac_ic                 ! For IC: :math:`2/(2 r_i)`
   real(cp), public :: alpha1   ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alpha2   ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
   real(cp), public :: r_cmb                     ! OC radius
   real(cp), public :: r_icb                     ! IC radius
   real(cp), public :: r_surface                 ! Surface radius for extrapolation in units of (r_cmb-r_icb)

   !-- arrays for buoyancy, depend on Ra and Pr:
   real(cp), public, allocatable :: rgrav(:)     ! Buoyancy term `dtemp0/Di`

   !-- chebychev polynomials, derivatives and integral:
   real(cp), public, allocatable :: cheb_int(:)     ! Array for cheb integrals
   integer, public :: nDd_costf1                    ! Radii for transform
   type(costf_odd_t), public :: chebt_ic
   type(costf_even_t), public :: chebt_ic_even

   !-- Radial scheme
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

   real(cp), public, allocatable :: lambda(:)     ! Array of magnetic diffusivity
   real(cp), public, allocatable :: dLlambda(:)   ! Derivative of magnetic diffusivity
   real(cp), public, allocatable :: jVarCon(:)    ! Analytical solution for toroidal field potential aj (see init_fields.f90)
   real(cp), public, allocatable :: sigma(:)      ! Electrical conductivity
   real(cp), public, allocatable :: kappa(:)      ! Thermal diffusivity
   real(cp), public, allocatable :: dLkappa(:)    ! Derivative of thermal diffusivity
   real(cp), public, allocatable :: visc(:)       ! Kinematic viscosity
   real(cp), public, allocatable :: dLvisc(:)     ! Derivative of kinematic viscosity
   real(cp), public, allocatable :: ddLvisc(:)    ! 2nd derivative of kinematic viscosity
   real(cp), public, allocatable :: divKtemp0(:)  ! Term for liquid anelastic approximation
   real(cp), public, allocatable :: epscProf(:)   ! Sources in heat equations
   real(cp), public, allocatable :: dr_top_ic(:)  ! Derivative in real space for r=r_i

   integer, public, allocatable :: l_R(:) ! Variable degree with radius

   public :: initialize_radial_functions, radial, transportProperties, &
   &         finalize_radial_functions

contains

   subroutine initialize_radial_functions
      !
      ! Initial memory allocation
      !

      integer :: n_in, n_in_2

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
      allocate( beta(n_r_max), dbeta(n_r_max), ddbeta(n_r_max) )
      allocate( alpha0(n_r_max), dLalpha0(n_r_max), ddLalpha0(n_r_max) )
      allocate( rgrav(n_r_max), ogrun(n_r_max) )
      allocate(drag0(n_r_max))
      bytes_allocated = bytes_allocated+(23*n_r_max+3*n_r_ic_max)*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( dxicond(n_r_max) )
         dxicond(:)=0.0_cp
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
      end if

      allocate( lambda(n_r_max),dLlambda(n_r_max),jVarCon(n_r_max) )
      allocate( sigma(n_r_max),kappa(n_r_max),dLkappa(n_r_max) )
      allocate( visc(n_r_max),dLvisc(n_r_max),ddLvisc(n_r_max) )
      allocate( epscProf(n_r_max),divKtemp0(n_r_max) )
      bytes_allocated = bytes_allocated + 11*n_r_max*SIZEOF_DEF_REAL

      !allocate ( l_R(nRstart:nRstop) )
      !bytes_allocated = bytes_allocated +(nRstop-nRstart+1)*SIZEOF_INTEGER
      allocate ( l_R(1:n_r_max) )
      bytes_allocated = bytes_allocated +n_r_max*SIZEOF_INTEGER

      if ( .not. l_full_sphere ) then
         nDi_costf1_ic=2*n_r_ic_max+2
         nDd_costf1_ic=2*n_r_ic_max+5
         nDi_costf2_ic=2*n_r_ic_max
         nDd_costf2_ic=2*n_r_ic_max+n_r_ic_max/2+5

         allocate( cheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( dcheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( d2cheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( cheb_int_ic(n_r_ic_max) )
         bytes_allocated = bytes_allocated + &
         &                 (3*n_r_ic_max*n_r_ic_max+n_r_ic_max)*SIZEOF_DEF_REAL

         call chebt_ic%initialize(n_r_ic_max,nDi_costf1_ic,nDd_costf1_ic)

         allocate ( dr_top_ic(n_r_ic_max) )
         bytes_allocated = bytes_allocated+n_r_ic_max*SIZEOF_DEF_REAL
      end if

      if ( .not. l_finite_diff ) then

         allocate( cheb_int(n_r_max) )         ! array for cheb integrals !
         bytes_allocated = bytes_allocated + n_r_max*SIZEOF_DEF_REAL

         allocate ( type_cheb_odd :: rscheme_oc )

         n_in = n_cheb_max
         if ( l_newmap ) then
            n_in_2 = 1
         else
            n_in_2 = 0
         end if

      else

         allocate ( type_fd :: rscheme_oc )

         n_in   = fd_order
         n_in_2 = fd_order_bound

      end if
      call rscheme_oc%initialize(n_r_max,n_in,n_in_2)

   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine finalize_radial_functions
      !
      ! Memory deallocation of radial functions
      !

      deallocate( l_R )
      deallocate( r, r_ic, O_r_ic, O_r_ic2, or1, or2, or3, or4 )
      deallocate( otemp1, rho0, temp0, dLtemp0, d2temp0, dentropy0 )
      deallocate( ddLtemp0, orho1, orho2, beta, dbeta, ddbeta, alpha0 )
      deallocate( ddLalpha0, dLalpha0, rgrav, ogrun )
      deallocate( lambda, dLlambda, jVarCon, sigma, kappa, dLkappa )
      deallocate( visc, dLvisc, ddLvisc, epscProf, divKtemp0 )
      deallocate(drag0)

      if ( l_chemical_conv ) deallocate(dxicond)

      if ( .not. l_full_sphere ) then
         deallocate( dr_top_ic )
         deallocate( cheb_ic, dcheb_ic, d2cheb_ic, cheb_int_ic )
         call chebt_ic%finalize()
         if ( n_r_ic_max > 0 .and. l_cond_ic ) call chebt_ic_even%finalize()
      end if

      if ( .not. l_finite_diff ) deallocate( cheb_int )

      call rscheme_oc%finalize()

   end subroutine finalize_radial_functions
!------------------------------------------------------------------------------
   subroutine radial()
      !
      !  Calculates everything needed for radial functions, transforms etc.
      !

      !-- Local variables:
      integer :: n_r,n_cheb,n_cheb_int
      integer :: n_r_ic_tot, k, i
      integer :: n_const(1)

      !integer :: n_r_start
      real(cp) :: fac_int
      real(cp) :: r_cheb_ic(2*n_r_ic_max-1),r_ic_2(2*n_r_ic_max-1)
      real(cp) :: drho0(n_r_max),dtemp0(n_r_max)

      real(cp) :: hcomp,fac
      real(cp) :: dtemp0cond(n_r_max),dtemp0ad(n_r_max),hcond(n_r_max)
      real(cp) :: func(n_r_max)

      real(cp), allocatable :: coeffDens(:), coeffTemp(:), coeffAlpha(:), coeffDentropy(:)
      real(cp), allocatable :: coeffGrav(:), coeffGrun(:)
      real(cp) :: rrOcmb(n_r_max)
      real(cp) :: dr_top(2*n_r_ic_max-1)
      character(len=80) :: message
      character(len=76) :: fileName
      integer :: fileHandle
      real(cp) :: ratio1, ratio2, diff, coeff

      !-- Radial grid point:
      !   radratio is aspect ratio
      !   radratio = (inner core r) / (CMB r) = r_icb/r_cmb
      if (l_tidal) then  !ARS
         r_cmb=one
         r_icb=radratio
      else
         r_cmb=one/(one-radratio)
         r_icb=r_cmb-one
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

      if ( rank == 0 ) then
         fileName = 'radius.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         do n_r=1,n_r_max
            write(fileHandle,'(I4, ES16.8)') n_r, r(n_r)
         end do
         close(fileHandle)
      end if

      if ( l_full_sphere ) then
         or1(:n_r_max-1)=one/r(:n_r_max-1)       ! 1/r
         or2(:n_r_max-1)=or1(:n_r_max-1)*or1(:n_r_max-1)  ! 1/r**2
         or3(:n_r_max-1)=or1(:n_r_max-1)*or2(:n_r_max-1)  ! 1/r**3
         or4(:n_r_max-1)=or2(:n_r_max-1)*or2(:n_r_max-1)  ! 1/r**4
         or1(n_r_max)=0.0_cp
         or2(n_r_max)=0.0_cp
         or3(n_r_max)=0.0_cp
         or4(n_r_max)=0.0_cp
      else
         or1(:)=one/r(:)       ! 1/r
         or2(:)=or1(:)*or1(:)  ! 1/r**2
         or3(:)=or1(:)*or2(:)  ! 1/r**3
         or4(:)=or2(:)*or2(:)  ! 1/r**4
      end if

      !-- Determine the max. degree for each radial level
      if ( l_var_l ) then ! Nat's form from Marti et al. (2014)
         !l_R(:) = int(one+(l_max-one)*sqrt(r(nRstart:nRstop)/r_cmb))
         l_R(:) = int(one+(l_max-one)*sqrt(r(:)/r_cmb))
         l_R(:) = int(one+l_max*sqrt(r(:)/r_cmb/rcut_l))
         do n_r=1,n_r_max
            if ( l_R(n_r) > l_max ) l_R(n_r)=l_max
         end do
         call logWrite('! Spherical harmonic degree varies with radius')
         call logWrite('! It increases between l_min and l_max following:')
         write(message,'(''!   l_min ='',i5, '' at r ='', f7.4)') minval(l_R(:)), &
         &     r(minloc(l_R(:),dim=1))
         call logWrite(message)
         write(message,'(''!   l_max ='',i5, '' at r ='', f7.4)') maxval(l_R(:)), &
         &     r(maxloc(l_R(:),dim=1))
         call logWrite(message)
         call logWrite('')
      else ! Default is constant l
         l_R(:) = l_max
      end if

      !-- Get entropy gradient
      call getEntropyGradient() ! By default this is zero
      drag0(:) = zero
      !-- Fit to an interior model
      if ( index(interior_model,'JUP') /= 0 ) then

         if ( l_non_adia ) then
            rrOcmb(:) = r(:)*r_cut_model/r_cmb
            rgrav(:)=  83.4792166296_cp*rrOcmb(:)+7.0970761715_cp*rrOcmb(:)**2 &
            &        -112.0000517766_cp*rrOcmb(:)**3+47.3447404648_cp*rrOcmb(:)**4

            allocate ( coeffAlpha(10), coeffTemp(10) )
            coeffAlpha = [ -12.9483344953_cp, 12.7631620079_cp, -60.0717008192_cp, &
            &              1.41916870466_cp, 755.055391736_cp, -1938.08838168_cp,  &
            &              952.893688457_cp, 2544.71502695_cp, -3703.20551213_cp,  &
            &              1440.95591192_cp]
            coeffTemp = [ 1.24462655e+05_cp, -2.85767595e+06_cp, 3.04794799e+07_cp,&
            &            -1.72807386e+08_cp, 5.83323621e+08_cp, -1.23322830e+09_cp,&
            &             1.64950647e+09_cp, -1.35626053e+09_cp, 6.25695974e+08_cp,&
            &             -1.23976043e+08_cp ]
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

            drho0=-ThExpNb*epsS*alpha0*temp0*dentropy0-DissNb/GrunNb*ogrun* &
            &      alpha0*rgrav
            call getBackground(drho0,0.0_cp,rho0)
            rho0=exp(rho0) ! this was ln(rho_0)
            beta=drho0

            ! The final stuff is always required
            call get_dr(beta,dbeta,n_r_max,rscheme_oc)
            call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
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
            &            -201.96655354_cp, 491.00495215_cp, -644.82401602_cp, &
            &            440.86067831_cp, -122.36071577_cp]

            coeffTemp = [0.999735638_cp, 0.0111053831_cp, 2.70889691_cp,  &
            &            -83.5604443_cp, 573.151526_cp, -1959.41844_cp,   &
            &            3774.39367_cp, -4159.56327_cp, 2447.75300_cp,    &
            &            -596.464198_cp]

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
         &            0.34973134_cp ]

         coeffTemp = [1.00294605_cp,-0.44357815_cp,13.9295826_cp,  &
         &            -137.051347_cp,521.181670_cp,-1044.41528_cp, &
         &            1166.04926_cp,-683.198387_cp, 162.962632_cp ]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'SUN') /= 0 ) then

         ! rho is negative beyond r_cut_model=0.9965
         ! radratio should be 0.7 (size of the Sun's CZ)
         ! This is Model S by JCD (e.g.: JCD+ Science 1996)

         allocate( coeffDens(6), coeffTemp(4) )
         coeffDens = [-24.83750402_cp, 231.79029994_cp, -681.72774358_cp, &
         &            918.30741266_cp,-594.30093367_cp, 150.76802942_cp ]

         coeffTemp = [5.53715416_cp, -8.10611274_cp, 1.7350452_cp, &
         &            0.83470843_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'GLIESE229B') /= 0 ) then
         ! Use also nVarDiff=2 with difExp=0.52

         allocate( coeffDens(8), coeffTemp(5) )
         coeffDens = [0.99879163_cp,0.15074601_cp,-4.20328423_cp,   &
         &            6.43542034_cp,-12.67297113_cp,21.68593078_cp, &
         &            -17.74832309_cp,5.35405134_cp]

         coeffTemp = [0.99784506_cp,0.16540448_cp,-3.44594354_cp,   &
         &            3.68189750_cp,-1.39046384_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'COROT3B') /= 0 ) then
         ! Use also nVarDiff=2 with difExp=0.62

         allocate( coeffDens(7), coeffTemp(7) )
         coeffDens = [1.00035987_cp,-0.01294658_cp,-2.78586315_cp,  &
         &            0.70289860_cp,2.59463562_cp,-1.65868190_cp,   &
         &            0.15984718_cp]

         coeffTemp = [1.00299303_cp,-0.33722671_cp,1.71340063_cp,     &
         &            -12.50287121_cp,21.52708693_cp,-14.91959338_cp, &
         &            3.52970611_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'KOI889B') /= 0 ) then
         ! Use also nVarDiff=2 with difExp=0.68

         allocate( coeffDens(6), coeffTemp(6) )
         coeffDens = [1.01038678_cp,-0.17615484_cp,-1.50567127_cp,  &
         &            -1.65738032_cp,4.20394427_cp,-1.87394994_cp]

         coeffTemp = [1.02100249_cp,-0.60750867_cp,3.23371939_cp,   &
         &            -12.80774142_cp,15.37629271_cp,-6.19288785_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'PNS_0V2S') /= 0 ) then
         allocate( coeffDens(6), coeffTemp(6) , coeffGrav(6) )
         coeffDens = [ 59.62268063_cp, -217.70124346_cp, &
              375.50063468_cp, -373.12779855_cp, 203.17403157_cp, -46.46897763_cp]

         coeffTemp = [-1.97868071_cp,   47.94900753_cp, -138.94186903_cp,  &
              176.44670046_cp, -109.98867268_cp,   27.51435739_cp]

         coeffGrav = [4.923996388662253_cp, -1.293518887100434_cp, &
              &     -26.037072977949265_cp, 55.59002270949802_cp,  &
              &     -45.9222941604744_cp, 13.739417490581724_cp]

         call polynomialBackground(coeffDens,coeffTemp,coeffGrav)
         deallocate( coeffDens, coeffTemp, coeffGrav )

      else if ( index(interior_model,'PNS_1S') /= 0 ) then
         allocate( coeffDens(6), coeffTemp(6) )
         coeffDens = [ 40.16045646_cp, -109.48920696_cp,  186.02134783_cp, &
              -204.74056538_cp, 105.44630174_cp,  -16.39774328_cp]

         coeffTemp = [ -20.1583557_cp, 160.49531633_cp, -412.24886308_cp, &
              520.77920394_cp, -334.07999285_cp, 86.21287296_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'PNS_2S') /= 0 ) then
         allocate( coeffDens(6), coeffTemp(6) , coeffGrav(6) )
         coeffDens = [14.5649560668_cp, -26.1462235672_cp, 52.8432513624_cp, &
              -82.837369536_cp, 58.2910747189_cp, -15.7212352085_cp]

         coeffTemp = [-6.33078724295_cp, 74.0516826175_cp, -202.268812325_cp, &
              258.6129266_cp, -161.782859987_cp, 38.7217973339_cp]

         coeffGrav = [0.0961974297099_cp, 3.06791065527_cp, 3.91870769838_cp, &
              -14.5167546463_cp, 11.5773447456_cp, -3.14354794653_cp]

         call polynomialBackground(coeffDens,coeffTemp,coeffGrav)
         deallocate( coeffDens, coeffTemp, coeffGrav )

      else if ( index(interior_model,'PNS_5S') /= 0 ) then
         allocate( coeffDens(6), coeffTemp(6) )
         coeffDens = [5.56892773949_cp, -0.532496149403_cp, -1.63704759716_cp, &
              -10.02819404_cp, 13.8052508126_cp, -6.16775910685_cp]

         coeffTemp = [4.32220656771_cp, -0.729605496043_cp, 0.860028766401_cp, &
              -9.29530062905_cp, 9.28443065708_cp, -3.44528845237_cp]

         call polynomialBackground(coeffDens,coeffTemp)
         deallocate( coeffDens, coeffTemp)

      else if ( index(interior_model,'PNS_SZ_0V2S') /= 0 ) then

         rrOcmb(:) = r(:)*r_cut_model/r_cmb
         allocate( coeffDens(11), coeffTemp(6) , coeffGrav(6), coeffDentropy(7))

         coeffTemp=[196.28690966001918_cp,-1105.906145023404_cp, &
              2554.7606628747303_cp, -2981.0616231930344_cp,     &
              1746.1572768996048_cp,-409.23926817272013_cp]

         coeffGrav=[-54.17356002706597_cp, 360.21272192378154_cp,&
              -876.2038899394643_cp, 1032.0971147446267_cp,      &
              -598.5680136354003_cp, 137.63554371201047_cp]

         coeffDens=[-2309723.2541369465_cp, 23358205.896083623_cp,&
              -103179719.55008087_cp,259280343.79339626_cp,       &
              -401991097.61232275_cp,383870769.0656396_cp,        &
              -199782584.99709758_cp,18574457.14277596_cp,        &
              41318238.790227205_cp,-23288018.044316825_cp,       &
              4149129.7702482087_cp]

         coeffDentropy=[-1296.7861983167516_cp,8524.65622292761_cp,&
              -23141.841394462575_cp,33187.22479826618_cp,           &
              -26488.875853067482_cp,11143.388161137767_cp,          &
              -1926.7650976608754_cp]

!         coeffAlpha=[-71.38015335545124_cp, 446.3489846730089_cp, &
!                     -1111.749355758775_cp, 1383.514366395739_cp, &
!                     -858.303326659732_cp, 212.57006335992187_cp]

         rho0(:) =0.0_cp
         temp0(:)=0.0_cp
         rgrav(:)=0.0_cp
         dentropy0(:)=0.0_cp
!         alpha0(:)=0.0_cp

         do i=1,11
            rho0(:) = rho0(:)+coeffDens(i)*rrOcmb(:)**(i-1)
         end do
         do i=1,7
            dentropy0(:)= dentropy0(:)+coeffDentropy(i)*rrOcmb(:)**(i-1)
         end do
         do i=1,6
            rgrav(:)  = rgrav(:) +coeffGrav(i)*rrOcmb(:)**(i-1)
            temp0(:)  = temp0(:) +coeffTemp(i)*rrOcmb(:)**(i-1)
!            alpha0(:) = alpha0(:)+coeffAlpha(i)*rr0cmb(:)**(i-1)
         end do

         ! Dissipation number
         !        DissNb   =alpha0(1)*rgrav(1)*(rrOcmb(1)-rrOcmb(n_r_max))*6.9894e7 &
         !        &         /1.5e4_cp ! 1.5e4 is cp 6.9894e7 is R_J
         !ThExpNb  =alpha0(1)*temp0(1) ! equal to 1 cuz normalisation
         dentropy0(:) = dentropy0(:)/dentropy0(1)
         temp0(:) =temp0(:)/temp0(1)

         ! Adiabatic: buoyancy term is linked to the temperature gradient
         !       dT
         !      ---- =  -Di * alpha_T * T * grav
         !       dr
         !         rgrav(:)=-dtemp0(:)/DissNb
         ! Non adiabactic: CAREFUL It's d ln(temp0)
         ! dtemp0(:)=epsS*dentropy0(:)-DissNb*alpha0(:)*rgrav(:)
         ! epsS = 1 ?

         call get_dr(temp0,dtemp0,n_r_max,rscheme_oc)
         dLtemp0(:)=dtemp0(:)/temp0(:)
!         dentropy0(:) = dtemp0(:)/temp0(:)/epsS+DissNb*alpha0(:)*rgrav(:)/epsS
!         dtemp0(:)+DissNb*alpha0(:)*rgrav(:) !To compare with fit from entropy gradient
         alpha0(:)=(dentropy0(:) - dLtemp0(:))/rgrav(:)
         ThExpNB = 1.19425932_cp
         DissNb   =alpha0(1) ! or take the predicted value 1.35789
         alpha0(:)=alpha0(:)/alpha0(1)

         !-- Multiply the gravity by alpha0 and temp0
         rgrav(:)=rgrav(:)*alpha0(:)*temp0(:)

         call get_dr(rho0,drho0,n_r_max,rscheme_oc)
         beta(:)=drho0(:)/rho0(:)
         call get_dr(beta,dbeta,n_r_max,rscheme_oc)
         call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
         call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
         call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
         dLalpha0(:)=dLalpha0(:)/alpha0(:) ! d log (alpha) / dr
         call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
         call get_dr(dLtemp0,ddLtemp0,n_r_max,rscheme_oc)
         !         dentropy0(:)=0.0_cp

         !- \Gamma = 1/\rho/c_v ( \partial p/\partial T)_\rho = (d\ln T/d \ln \rho)_s
         ogrun(:)=beta(:)/dLtemp0(:)
         GrunNb  =one/ogrun(1)
         strat   =log(rho0(n_r_max)/rho0(1))
         polind  =ogrun(1)

         l_non_adia = .true.

         !         call polynomialBackground(coeffDens,coeffTemp,coeffGrav)
         deallocate( coeffDens, coeffTemp, coeffGrav, coeffDentropy)

      else if (index(interior_model,'HMNS_17MS_SZ') /= 0 ) then

         rrOcmb(:) = r(:)*r_cut_model/r_cmb
         allocate( coeffDens(14), coeffTemp(15) , coeffGrav(14), coeffDentropy(14),coeffAlpha(14),coeffGrun(14))

         coeffDens =[8.35611339e+00_cp,8.74463291e+01_cp,3.93867939e+02_cp,1.00096732e+03_cp,&
                     1.58160389e+03_cp,1.61734228e+03_cp,1.08343137e+03_cp,4.69310868e+02_cp,&
                     1.22509898e+02_cp,1.55908307e+01_cp,9.05893611e-01_cp,5.62581360e-01_cp,&
                     -3.44411392e+00_cp,1.28395168e-05_cp]
         coeffTemp=[-3.84479404e+05_cp,  3.20226839e+06_cp, -1.18068483e+07_cp,  2.51872234e+07_cp,&
              -3.39009704e+07_cp,  2.89707086e+07_cp, -1.40487741e+07_cp,  1.15447232e+06_cp,&
              3.47092643e+06_cp, -2.66927071e+06_cp,  1.03261814e+06_cp, -2.38022047e+05_cp,&
              3.24775159e+04_cp, -2.40872491e+03_cp,  7.99707136e+01_cp] -[-3.98974982e-04_cp,&
              -2.34290259e-03_cp,  1.71269197e-02_cp,  9.32427496e-03_cp, &
              2.52106711e-02_cp, -1.50569007e-02_cp,  4.12280355e-02_cp,  4.60373797e-03_cp, &
              -1.30971847e-03_cp,  2.29658000e-03_cp, -3.06276721e-04_cp,  3.28790455e-04_cp,&
              2.58963555e-05_cp,  3.72774002e-06_cp, -1.92877536e-08_cp]
       coeffGrav=[5.12277059e+02_cp,5.22698735e+03_cp,2.34743426e+04_cp,6.09337549e+04_cp,&
              1.00879708e+05_cp,1.10774223e+05_cp,8.11629844e+04_cp,3.87830547e+04_cp,&
              1.13586630e+04_cp,1.76008849e+03_cp,8.42637402e+01_cp,-5.93291222e+00_cp,&
              -1.71807798e+00_cp,1.75079064e-03_cp]+[ 9.02172133e-08_cp,  1.76800404e-06_cp, &
              4.31695698e-05_cp,  4.84447810e-07_cp, &
              -1.31286375e-04_cp,  2.88932337e-04_cp, -7.39409006e-06_cp,  3.26496956e-05_cp,&
              -3.11257245e-05_cp, -1.27490648e-06_cp, -4.76843383e-08_cp,  1.15867227e-09_cp,&
              4.83339146e-09_cp, -2.29841645e-12_cp]
         coeffDentropy=[-2.71186941e+03_cp,-2.72715819e+04_cp,-1.21184650e+05_cp,-3.12318398e+05_cp,&
              -5.14840606e+05_cp,-5.64265219e+05_cp,-4.13689488e+05_cp,-1.98671893e+05_cp,&
              -5.91741985e+04_cp,-9.74141046e+03_cp,-6.80402916e+02_cp,-1.83548419e+01_cp,&
              -7.74080965e+00_cp,6.84110279e-03_cp]+[ 3.63495019e-06_cp, -3.64270709e-05_cp,&
              2.61578753e-04_cp,  1.51584973e-05_cp, &
              -3.68003093e-04_cp, -4.80020070e-04_cp,  3.40328552e-06_cp,  2.39534245e-04_cp,&
              -8.25224561e-07_cp, -1.77185939e-06_cp,  3.64710786e-07_cp, -4.35362324e-09_cp,&
              1.99811190e-09_cp,  4.68442125e-12_cp]
         coeffAlpha=[9.99440317e+02_cp,9.17953447e+03_cp,3.68937854e+04_cp,8.51100015e+04_cp,&
              1.24178202e+05_cp,1.18971545e+05_cp,7.51727200e+04_cp,3.05910833e+04_cp,&
              7.55168185e+03_cp,9.89230714e+02_cp,4.41378488e+01_cp,-9.10433411e-01_cp, &
              7.34535366e-01_cp,4.59734984e-04_cp]
         coeffGrun=[-2.49531149e+02_cp,-2.14574152e+03_cp,-7.98060427e+03_cp,-1.67920613e+04_cp,&
              -2.19272303e+04_cp,-1.82978917e+04_cp,-9.62716909e+03_cp,-2.97363833e+03_cp,&
              -4.22358839e+02_cp,8.67206864e+00_cp,5.15346545e+00_cp,-1.16955754e+00_cp,&
              7.24613039e-02_cp,-8.39110481e-05_cp] + [-5.62233993e-08_cp, -2.93029325e-06_cp, &
              1.60124728e-07_cp,  3.16922415e-05_cp, -3.04100031e-05_cp, -4.20354772e-05_cp, &
              2.62865979e-06_cp, -2.31091053e-06_cp, -4.94001711e-07_cp, -1.41527678e-10_cp, &
              -9.29925470e-10_cp, -4.53345739e-09_cp, -1.09581649e-11_cp, -1.89740937e-14_cp]

         rho0(:) =0.0_cp
         temp0(:)=0.0_cp
         rgrav(:)=0.0_cp
         dentropy0(:)=0.0_cp
         alpha0(:)=0.0_cp
         ogrun(:) =0.0_cp
         do i=1,14
            rho0(:) = rho0(:)+coeffDens(i)*log(rrOcmb(:))**(14-i)
            rgrav(:) = rgrav(:)+coeffGrav(i)*log(rrOcmb(:))**(14-i)
            dentropy0(:) = dentropy0(:)+coeffDentropy(i)*log(rrOcmb(:))**(14-i)
            alpha0(:) = alpha0(:)+coeffAlpha(i)*log(rrOcmb(:))**(14-i)
            ogrun(:) = ogrun(:)+coeffGrun(i)*log(rrOcmb(:))**(14-i)
         end do
         rho0(:) =exp(rho0(:))
         rgrav(:)=exp(rgrav(:))
         dentropy0(:)=exp(dentropy0(:))
         alpha0(:)=exp(alpha0(:))
         ogrun(:) =exp(ogrun(:))
         do i=1,15
            temp0(:)= temp0(:)+coeffTemp(i)*(rrOcmb(:))**(15-i)
         end do

         dentropy0(:) = dentropy0(:)/dentropy0(1)
!         temp0(:) =temp0(:)/temp0(1)
         rho0(:) =rho0(:)/rho0(1)
         rgrav(:) =rgrav(:)/rgrav(1)
         ogrun(:) = ogrun(:)/ogrun(1)
         alpha0(:) = alpha0(:)/alpha0(1)

         if (ra*dentropy0(1)<0) then
            dentropy0(:) = -dentropy0(:)
         end if

         rgrav(:)=rgrav(:)*alpha0(:)*temp0(:)
         call get_dr(temp0,dtemp0,n_r_max,rscheme_oc)
         dLtemp0(:)=dtemp0(:)/temp0(:)
         call get_dr(rho0,drho0,n_r_max,rscheme_oc)
         beta(:)=drho0(:)/rho0(:)
         call get_dr(beta,dbeta,n_r_max,rscheme_oc)
         call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
         call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
         call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
         dLalpha0(:)=dLalpha0(:)/alpha0(:) ! d log (alpha) / dr
         call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
         call get_dr(dLtemp0,ddLtemp0,n_r_max,rscheme_oc)

         DissNb=2.125815_cp*(1-radratio)/0.77472611
         ThExpNB= 2.1687894433_cp*(1-radratio)/0.77472611
         GrunNb =  0.508170_cp
         polind = 1/GrunNb
         strat   =log(rho0(n_r_max)/rho0(1))
         l_non_adia = .true.

         if (l_drag) then
            drag0(:) = (6000.0_cp*(0.5390537718711352_cp*temp0)**6)*1.7342490808689306_cp !damping rate in viscous units
         end if

         deallocate( coeffDens, coeffTemp, coeffGrav, coeffDentropy,coeffAlpha,coeffGrun)

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
         call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
         call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
         call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
         dLalpha0=dLalpha0/alpha0 ! d log (alpha) / dr
         call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
         dLtemp0 = dtemp0/temp0
         call get_dr(dLtemp0,ddLtemp0,n_r_max,rscheme_oc)

         ! N.B. rgrav is not gravity but alpha * grav
         rgrav = alpha0*rgrav

         !-- ogrun
         ogrun(:)=one/GrunNb

         l_non_adia = .true.

      else if (index(interior_model,'MESA_5M_ZAMS') /= 0) then

         l_non_adia = .true.
         rrOcmb = r*r_cut_model/r_cmb

         allocate( coeffAlpha(12), coeffTemp(12), coeffGrav(12), coeffGrun(12),&
         &         coeffDens(12) )

         coeffAlpha=[1.07280013e+00_cp, -7.35271452e-02_cp,  1.35373501e+00_cp,&
                 &  -2.70033254e+01_cp,  1.34218482e+02_cp, -1.85203678e+02_cp,&
                 &  -4.03705669e+02_cp,  1.65171816e+03_cp, -1.91235980e+03_cp,&
                 &   4.33699991e+02_cp,  6.74102545e+02_cp, -3.67415688e+02_cp]

         coeffTemp=[1.71233163e+01_cp,  2.16376708e-01_cp, -2.33617981e+01_cp, &
              &     2.95986199e+02_cp, -3.15614535e+03_cp,  1.80265801e+04_cp, &
              &    -5.99840212e+04_cp,  1.23788205e+05_cp, -1.61112426e+05_cp, &
              &     1.28833197e+05_cp, -5.78439543e+04_cp,  1.11724125e+04_cp]

         coeffGrav=[1.31049693e+02_cp,  1.04893208e+06_cp,  1.93822208e+06_cp, &
              &    -5.15044761e+07_cp,  4.41798276e+08_cp, -2.66511461e+09_cp, &
              &     1.02063273e+10_cp, -2.43712757e+10_cp,  3.64094408e+10_cp, &
              &    -3.32103850e+10_cp,  1.69572548e+10_cp, -3.72150208e+09_cp]

         coeffGrun=[1.57026325e+00_cp, -4.46247381e-02_cp,  3.03734351e-01_cp, &
              &    -8.28253239e+00_cp, -2.63234309e+01_cp,  6.05835596e+02_cp, &
              &    -2.89779916e+03_cp,  6.90494749e+03_cp, -9.29010895e+03_cp, &
              &     7.04974613e+03_cp, -2.73663642e+03_cp,  3.98015351e+02_cp]

         coeffDens=[3.05336249e+00_cp, -2.52357616e-01_cp, -3.47786718e-01_cp, &
              &    -3.82246141e+02_cp,  4.44024160e+03_cp, -2.80734664e+04_cp, &
              &     1.02942741e+05_cp, -2.30659612e+05_cp,  3.22111291e+05_cp, &
                   -2.74248402e+05_cp,  1.30456564e+05_cp, -2.66069272e+04_cp]

         alpha0(:)=0.0_cp
         temp0(:) =0.0_cp
         rgrav(:) =0.0_cp
         ogrun(:) =0.0_cp
         dtemp0(:)=0.0_cp
         rho0(:)  =0.0_cp
         drho0(:) =0.0_cp

         do i=1,12
            alpha0(:) = alpha0(:)+coeffAlpha(i)*rrOcmb(:)**(i-1)
            temp0(:)  = temp0(:) +coeffTemp(i) *rrOcmb(:)**(i-1)
            rgrav(:)  = rgrav(:) +coeffGrav(i) *rrOcmb(:)**(i-1)
            ogrun(:)  = ogrun(:) +coeffGrun(i) *rrOcmb(:)**(i-1)
            rho0(:)   = rho0(:)  +coeffDens(i) *rrOcmb(:)**(i-1)
         end do

         do i=2,12 ! d (lnT)/ dr
            dtemp0(:) = dtemp0(:)+coeffTemp(i)*(i-1)*rrOcmb(:)**(i-2)
            drho0(:)  = drho0(:) +coeffDens(i)*(i-1)*rrOcmb(:)**(i-2)
         end do

         alpha0(:)=alpha0(:)/alpha0(1)
         rgrav(:) =rgrav(:)/rgrav(1)

         !Normalise
         call getBackground(dtemp0,0.0_cp,temp0)
         call getBackground(drho0,0.0_cp,rho0)
         temp0 = exp(temp0) ! Fit was on ln(temp0)
         rho0  = exp(rho0)  ! Fit was on ln(rho0)
         drho0 = drho0 * rho0 ! Fits were on ln(rho0) and ln(temp0)
         dtemp0= dtemp0* temp0

         beta = drho0

         !-- Final stuff
         call get_dr(beta,dbeta,n_r_max,rscheme_oc)
         call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
         call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
         call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
         dLalpha0=dLalpha0/alpha0 ! d log (alpha) / dr
         call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
         dLtemp0 = dtemp0/temp0
         ddLtemp0 =-(dtemp0/temp0)**2+d2temp0/temp0

         !-- Multiply the gravity by alpha0 and temp0
         rgrav(:)=rgrav(:)*alpha0(:)*temp0(:)

      else  !-- Usual polytropic reference state
         ! g(r) = g0 + g1*r/ro + g2*(ro/r)**2
         ! Default values: g0=0, g1=1, g2=0
         ! An easy way to change gravity
         rgrav(:)=g0+g1*r(:)/r_cmb+g2*r_cmb**2*or2(:)

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
               call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
               call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)

               dLtemp0 = dtemp0/temp0
               ddLtemp0 =-(dtemp0/temp0)**2+d2temp0/temp0

               ogrun(:) = one/GrunNb

            else !-- Adiabatic reference state

               if ( l_isothermal ) then ! Gruneisen is zero in this limit
                  if ( l_full_sphere ) then
                     fac = strat/(half*g1*(one+radratio))/(r_cmb-r_icb)
                  else
                     fac=strat / ( g0+half*g1*(one+radratio)+g2/radratio ) &
                     &         / (r_cmb-r_icb)
                  end if
                  DissNb     =0.0_cp
                  GrunNb     =0.0_cp
                  ogrun(:)   =0.0_cp
                  temp0(:)   =one
                  rho0(:)    =exp(-fac*(g0*r(:)+half*g1*r(:)**2/r_cmb-g2*r_cmb**2/&
                  &               r(:))+fac*r_cmb*(g0+half*g1-g2))

                  beta(:)    =-fac*rgrav(:)
                  dbeta(:)   =-fac*(g1/r_cmb-two*g2*r_cmb**2*or3(:))
                  ddbeta(:)  =6.0_cp*fac*g2*r_cmb**2*or4(:)
                  d2temp0(:) =0.0_cp
                  dLtemp0(:) =0.0_cp
                  ddLtemp0(:)=0.0_cp
               else
                  if ( strat == 0.0_cp .and. DissNb /= 0.0_cp ) then
                     if ( l_full_sphere ) then
                        strat = polind* log( (half*g1*(one+radratio))*  &
                        &                    (r_cmb-r_icb)*DissNb+1 )
                     else
                        strat = polind* log( (g0+half*g1*(one+radratio)+g2/radratio)*&
                        &                    (r_cmb-r_icb)*DissNb+1 )
                     end if
                  else
                     if ( l_full_sphere ) then
                        DissNb=( exp(strat/polind)-one )/(r_cmb-r_icb)/ &
                        &      ( half*g1*(one+radratio) )
                     else
                        DissNb=( exp(strat/polind)-one )/(r_cmb-r_icb)/ &
                        &      ( g0+half*g1*(one+radratio) +g2/radratio )
                     end if
                  end if
                  GrunNb  =one/polind
                  ogrun(:)=one/GrunNb
                  temp0(:)=-DissNb*( g0*r(:)+half*g1*r(:)**2/r_cmb-   &
                  &         g2*r_cmb**2*or1(:) ) + one + DissNb*r_cmb* &
                  &         (g0+half*g1-g2)
                  rho0(:) =temp0**polind

                  !if (l_centrifuge) opressure0(:) = temp0**(-polind-1)

                  !-- Computation of beta= dln rho0 /dr and dbeta=dbeta/dr
                  beta(:)     =-polind*DissNb*rgrav(:)/temp0(:)
                  dbeta(:)    =-polind*DissNb/temp0(:)**2 *         &
                  &            ((g1/r_cmb-two*g2*r_cmb**2*or3(:))*  &
                  &            temp0(:)  + DissNb*rgrav(:)**2)
                  ddbeta(:)   =-polind*DissNb/temp0(:)**3 * (   temp0(:)*          &
                  &            ( 6.0_cp*g2*r_cmb**2*or4(:)*temp0(:)+               &
                  &              DissNb*rgrav(:)*(g1/r_cmb-two*g2*r_cmb**2*or3(:)) &
                  &            )+two*DissNb*rgrav(:)*(DissNb*rgrav(:)**2+          &
                  &             (g1/r_cmb-two*g2*r_cmb**2*or3(:))*temp0(:) ) )
                  dtemp0(:)   =-DissNb*rgrav(:)
                  d2temp0(:)  =-DissNb*(g1/r_cmb-two*g2*r_cmb**2*or3(:))
                  dLtemp0(:)  =dtemp0(:)/temp0(:)
                  ddLtemp0(:) =-(dtemp0(:)/temp0(:))**2+d2temp0(:)/temp0(:)
               end if

            end if

            !-- Thermal expansion coefficient (1/T for an ideal gas)
            alpha0(:)   =one/temp0(:)
            dLalpha0(:) =-dLtemp0(:)
            ddLalpha0(:)=-ddLtemp0(:)

         end if
      end if

      if ( l_anel ) then
         call logWrite('')
         call logWrite('!      This is an anelastic model')
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
         orho1(:)   =one/rho0(:)
         orho2(:)   =orho1(:)*orho1(:)
         otemp1(:)  =one/temp0(:)
         ViscHeatFac=DissNb*pr/raScaled
         if (l_mag) then
            OhmLossFac=ViscHeatFac/(ekScaled*prmag**2)
         else
            OhmLossFac=0.0_cp
         end if
      else ! Boussinesq
         rho0(:)     =one
         temp0(:)    =one
         !opressure0(:)=one
         otemp1(:)   =one
         orho1(:)    =one
         orho2(:)    =one
         alpha0(:)   =one
         ogrun(:)    =0.0_cp
         beta(:)     =0.0_cp
         dbeta(:)    =0.0_cp
         ddbeta(:)   =0.0_cp
         dLalpha0(:) =0.0_cp
         ddLalpha0(:)=0.0_cp
         dLtemp0(:)  =0.0_cp
         ddLtemp0(:) =0.0_cp
         d2temp0(:)  =0.0_cp
         ViscHeatFac =0.0_cp
         OhmLossFac  =0.0_cp
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
         call cheb_grid(-r_icb,r_icb,n_r_ic_tot-1,r_ic_2,r_cheb_ic,  &
              &         0.0_cp,0.0_cp,0.0_cp,0.0_cp,.false.)

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

      end if

      if ( n_r_ic_max > 0 .and. l_cond_ic ) then

         dr_fac_ic=two/(two*r_icb)
         cheb_norm_ic=sqrt(two/real(n_r_ic_max-1,kind=cp))

         !----- Calculate the even Chebs and their derivative:
         !      n_r_ic_max even chebs up to degree 2*n_r_ic_max-2
         !      at the n_r_ic_max first points in r_ic [r_icb,0].
         !      NOTE invers order in r_ic!
         call get_chebs_even(n_r_ic_max,-r_icb,r_icb,r_cheb_ic, &
              &                    n_r_ic_max,cheb_ic,dcheb_ic, &
              &                 d2cheb_ic,n_r_ic_max,n_r_ic_max)

         !----- Initialize transforms:
         call chebt_ic_even%initialize(n_r_ic_max-1,nDi_costf2_ic,nDd_costf2_ic)

         !----- Factors for cheb integrals, only even contribution:
         fac_int=one/dr_fac_ic   ! thats 1 for the outer core
         cheb_int_ic(1)=fac_int   ! Integration constant chosen !
         do n_cheb=2,n_r_ic_max
            n_cheb_int=2*n_cheb-1
            cheb_int_ic(n_cheb)=-fac_int / real(n_cheb_int*(n_cheb_int-2),kind=cp)
         end do

         dr_top(1)=(two*(2*n_r_ic_max-2)*(2*n_r_ic_max-2)+one)/6.0_cp
         do n_r=2,2*n_r_ic_max-1
            diff = two* sin( half*(n_r-1)*pi/(2*n_r_ic_max-2) ) * &
            &           sin( half*(n_r-1)*pi/(2*n_r_ic_max-2) )
            if ( mod(n_r,2) == 0 ) then
               coeff=-one
            else
               coeff=one
            end if
            dr_top(n_r)=two*coeff/diff
         end do
         dr_top(2*n_r_ic_max-1) = half*dr_top(2*n_r_ic_max-1)

         !-- This array can be used to compute the derivatives of the
         !-- inner core quantities at radius r_i in real space. This is
         !-- based on Chebyshev differentiation matrices in real space
         !-- and make use of symmetry properties.
         do n_r=1,n_r_ic_max-1
            dr_top_ic(n_r) =dr_top(n_r)+dr_top(2*n_r_ic_max-n_r)
         end do
         dr_top_ic(n_r_ic_max) =dr_top(n_r_ic_max)

         !-- Cheb factor of the form 2/(b-a)=2/(2*ri)=1/ri
         dr_top_ic(:) = dr_top_ic(:)*dr_fac_ic

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
      real(cp) :: a0,a1,a2,a3,a4,a5,a6,a7
      real(cp) :: kappatop,rrOcmb, ampVisc, ampKap, slopeVisc, slopeKap

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
          else if ( nVarCond == 5 ) then
             ! 'PNS magnetic diffusivity profile'
             ! this formula assume degenerate protons
             lambda = temp0**2/rho0**(3./2.)
             sigma = one/lambda
             call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
             dLlambda=dsigma/lambda
          else if ( nVarCond == 6 ) then
             ! PNS magnetic diffusivity profile
             ! degenerate relativistic electrons
             ! non degenerate protons
             ! Thompson & Duncan, ApJ (1993)
             ! This formula neglect the Ye dependence
             lambda = rho0**(-1./3.)
             sigma = one/lambda
             call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
             dLlambda=dsigma/lambda
          else if ( nVarCond == 7 ) then
             ! PNS magnetic diffusivity profile
             ! Thompson & Duncan, ApJ (1993)
             ! lambda \propto (\rho Y_e)**(-1/3)
             if ( index(interior_model, 'PNS_5S') /= 0 ) then
                !! Y_e polynomial fit
                a0 = 3.6070973312516563
                a1 = 0.024205264127999757
                a2 = -1.623849893276523
                a3 = 1.9655839620692974
                a4 = -2.9551512201779615
                do n_r=1,n_r_max
                   rrOcmb = r(n_r)/r_cmb*r_cut_model
                   lambda(n_r)= (rho0(n_r)*(a0 + a1*rrOcmb + a2*rrOcmb**2 + a3*rrOcmb**3 + a4*rrOcmb**4))**(-1./3.)
                end do
             else
                call abortRun('Conductivity profile not defined at this time')
             end if
             lambda = lambda/lambda(1)
             sigma = one/lambda
             call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
             dLlambda=dsigma/lambda
          else if ( nVarCond == 8 ) then ! fit PNS magnetic diffusivity
             if ( index(interior_model, 'PNS_2S') /= 0 ) then
                a0 = -29.3468638288
                a1 = 335.929657508
                a2 = -1617.65787705
                a3 = 4286.501036
                a4 = -6746.07848196
                a5 = 6308.95415428
                a6 = -3247.58459547
                a7 = 710.280067712
             else
                call abortRun('Conductivity profile not defined')
             end if
             do n_r=1,n_r_max
                rrOcmb = r(n_r)/r_cmb*r_cut_model
                lambda(n_r)= a0 + a1*rrOcmb    + a2*rrOcmb**2 &
                     + a3*rrOcmb**3 + a4*rrOcmb**4 &
                     + a5*rrOcmb**5 + a6*rrOcmb**6 &
                     + a7*rrOcmb**7
             end do
             lambda = lambda/lambda(1)
             sigma = one/lambda
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
               write(output_unit,*) '! NOTE: with this polynomial fit     '
               write(output_unit,*) '! for variable thermal conductivity  '
               write(output_unit,*) '! considering radratio < 0.2 may lead'
               write(output_unit,*) '! to strange profiles'
               call abortRun('Stop the run in radial.f90')
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
         else if ( nVarDiff == 6 ) then ! Jump in the stratified layer
            ampKap = 10.0_cp
            slopeKap = 100.0_cp
            if ( rStrat <= r_icb ) then
               kappa(:) = one
            else
               kappa(:)=(-half*(ampKap-one)*tanh(slopeKap*(r(:)-rStrat))+      &
                    &          half*(ampKap+one))*(half*(ampKap-one)*tanh(slopeKap* &
                    &          (r(:)-rStrat-thickStrat))+half*(ampKap+one))/ampKap
            end if
            dLkappa(:)=ampKap*(slopeKap*(-half*ampKap + half)*(-tanh(slopeKap*(r &
            &         - rStrat))**2 + 1)*(half*ampKap + (half*ampKap - half)*    &
            &         tanh(slopeKap*(r - rStrat - thickStrat)) + half)/ampKap +  &
            &         slopeKap*(half*ampKap - half)*(-tanh(slopeKap*(r - rStrat -&
            &         thickStrat))**2 + 1)*(half*ampKap + (-half*ampKap + half)* &
            &         tanh(slopeKap*(r - rStrat)) + half)/ampKap)/((half*ampKap +&
            &         (-half*ampKap + half)*tanh(slopeKap*(r - rStrat)) + half)* &
            &         (half*ampKap + (half*ampKap - half)*tanh(slopeKap*(r -     &
            &         rStrat - thickStrat)) + half))
          else if ( nVarDiff == 7 ) then ! Bottom stratified
             ampKap = 10.0_cp
             slopeKap = 30.0_cp
             if ( rStrat <= r_icb ) then
                kappa(:) = one
             else
                kappa(:)=(half*(ampKap-one)*tanh(slopeKap* &
                &       (r(:)-rStrat))+half*(ampKap+one))/ampKap
             end if
             dLkappa(:)=slopeKap*(half*ampKap-half)*(-tanh(slopeKap*(r(:)-rStrat))**2&
             &         +one)/(half*ampKap+(half*ampKap-half)*tanh(slopeKap*(r(:)-    &
             &         rStrat))+half)
         else if (nVarDiff == 8) then ! fit PNS thermal diffusivity
            ! warning: reversed order for coefficients a_i
            ! compared to nVarDiff==3 above
            if ( index(interior_model,'PNS_0V2S') /= 0 ) then
               a0 = -15.4045660979
               a1 = 167.416944961
               a2 = -773.523043238
               a3 = 1971.51233444
               a4 = -2995.52988304
               a5 = 2717.20296952
               a6 = -1364.75385894
               a7 = 294.078979841
            else if ( index(interior_model, 'PNS_1S') /= 0 ) then
               a0 = -1892.77098523
               a1 = 17646.3340418
               a2 = -70253.2335883
               a3 = 154843.933552
               a4 = -204092.633188
               a5 = 160897.639922
               a6 = -70265.4062536
               a7 = 13117.1344588
            else if ( index(interior_model, 'PNS_2S') /= 0 ) then
               a0 = -95.9735911526
               a1 = 1087.08663925
               a2 = -5203.4624924
               a3 = 13653.9408117
               a4 = -21218.8320365
               a5 = 19540.2406991
               a6 = -9879.84517647
               a7 = 2117.83593177
            else if ( index(interior_model, 'PNS_5S') /= 0 ) then
               a0 = 0.0556428940038
               a1 = 1.57277204681
               a2 = -17.4234841566
               a3 = 92.3174440711
               a4 = -256.540091374
               a5 = 387.399928018
               a6 = -299.886696812
               a7 = 93.4823970453
            else
               call abortRun('Stop the run in radial.f90. Thermal diffusivity profile is not defined !')
            endif
            do n_r=1,n_r_max
               rrOcmb = r(n_r)/r_cmb*r_cut_model
               kappa(n_r)= a0 + a1*rrOcmb    + a2*rrOcmb**2 &
                    + a3*rrOcmb**3 + a4*rrOcmb**4 &
                    + a5*rrOcmb**5 + a6*rrOcmb**6 &
                    + a7*rrOcmb**7
            end do
            kappatop=kappa(1) ! normalise by the value at the top
            kappa=kappa/kappatop
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 9 ) then
            ! PNS thermal diffusivity
            ! approximation given by Thompson & Duncan (93), Eq. 7
            kappa = otemp1*rho0**(-2.0/3.0)
            kappatop=kappa(1) ! normalise by the value at the top
            kappa=kappa/kappatop
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         else if ( nVarDiff == 10 ) then
            ! PNS thermal diffusivity
            ! approximation given by Thompson & Duncan (93), Eq. 7
            kappa = rho0**(-4.0/3.0)
            kappatop=kappa(1) ! normalise by the value at the top
            kappa=kappa/kappatop
            call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
            dLkappa=dkappa/kappa
         end if
      end if

      !-- Eps profiles
      !-- The remaining division by rho will happen in updateS.f90
      if ( nVarEps == 0 ) then
         ! eps is constant
         if ( l_anelastic_liquid ) then
            epscProf(:)=one
         else
            epscProf(:)=otemp1(:)
         end if
      else if ( nVarEps == 1 ) then
         ! rho*eps in the RHS
         if ( l_anelastic_liquid ) then
            epscProf(:)=rho0(:)
         else
            epscProf(:)=rho0(:)*otemp1(:)
         end if
      else if ( nVarEps == 2 ) then
         ! rho*temp*eps in the RHS
         if ( l_anelastic_liquid ) then
            epscProf(:)=rho0(:)*temp0(:)
         else
            epscProf(:)=rho0(:)
         end if
      else if ( nVarEps == 3 ) then
         ! eps*rho**2*temp**(-3)*exp(-Bn/T) in the RHS
         if ( l_anelastic_liquid ) then
            epscProf(:)=rho0(:)**2*temp0(:)**(-3)*exp(-Bn/temp0(:))
         else
            epscProf(:)=rho0(:)**2*temp0(:)**(-4)*exp(-Bn/temp0(:))
         end if
      end if

      !-- Variable viscosity
      if ( nVarVisc == 0 ) then ! default: constant kinematic viscosity
         visc(:)   =one
         dLvisc(:) =0.0_cp
         ddLvisc(:)=0.0_cp
      else if ( nVarVisc == 1 ) then ! Constant dynamic viscosity
         visc(:)=rho0(n_r_max)/rho0(:)
         call get_dr(visc,dvisc,n_r_max,rscheme_oc)
         dLvisc(:)=dvisc(:)/visc(:)
         call get_dr(dLvisc,ddLvisc,n_r_max,rscheme_oc)
      else if ( nVarVisc == 2 ) then ! Profile
         visc(:)=(rho0/rho0(n_r_max))**difExp
         call get_dr(visc,dvisc,n_r_max,rscheme_oc)
         dLvisc(:)=dvisc(:)/visc(:)
         call get_dr(dLvisc,ddLvisc,n_r_max,rscheme_oc)
      else if ( nVarVisc == 3 ) then ! Jump in the stratified layer
         ampVisc = 10.0_cp
         slopeVisc = 100.0_cp
         if ( rStrat <= r_icb ) then
            visc(:) = one
         else
            visc(:)=(-half*(ampVisc-one)*tanh(slopeVisc*(r(:)-rStrat))+       &
            &         half*(ampVisc+one))*(half*(ampVisc-one)*tanh(slopeVisc* &
            &         (r(:)-rStrat-thickStrat))+half*(ampVisc+one))/ampVisc
         end if
         dLvisc(:)=ampVisc*(slopeVisc*(-half*ampVisc + half)*(-tanh(slopeVisc*(r &
         &         - rStrat))**2 + 1)*(half*ampVisc + (half*ampVisc - half)*     &
         &         tanh(slopeVisc*(r - rStrat - thickStrat)) + half)/ampVisc +   &
         &         slopeVisc*(half*ampVisc - half)*(-tanh(slopeVisc*(r - rStrat -&
         &         thickStrat))**2 + 1)*(half*ampVisc + (-half*ampVisc + half)*  &
         &         tanh(slopeVisc*(r - rStrat)) + half)/ampVisc)/((half*ampVisc +&
         &         (-half*ampVisc + half)*tanh(slopeVisc*(r - rStrat)) + half)*  &
         &         (half*ampVisc + (half*ampVisc - half)*tanh(slopeVisc*(r -     &
         &         rStrat - thickStrat)) + half))

         ddLvisc(:)=-ampVisc*slopeVisc*(-half*ampVisc + half)*(slopeVisc*(-half* &
         &          ampVisc + half)*(-tanh(slopeVisc*(r(:) - rStrat))**2 + 1)*   &
         &          (half*ampVisc + (half*ampVisc - half)*tanh(slopeVisc*(r(:) - &
         &          rStrat - thickStrat)) + half)/ampVisc + slopeVisc*(half*     &
         &          ampVisc - half)*(-tanh(slopeVisc*(r(:) - rStrat - thickStrat)&
         &          )**2 + 1)*(half*ampVisc + (-half*ampVisc + half)*tanh(       &
         &          slopeVisc*(r(:) - rStrat)) + half)/ampVisc)*(-tanh(slopeVisc*&
         &          (r(:) - rStrat))**2 + 1)/((half*ampVisc + (-half*ampVisc +   &
         &          half)*tanh(slopeVisc*(r(:) - rStrat)) + half)**2*(half*      &
         &          ampVisc + (half*ampVisc - half)*tanh(slopeVisc*(r(:) - rStrat&
         &          - thickStrat)) + half)) - ampVisc*slopeVisc*(half*ampVisc -  &
         &          half)*(slopeVisc*(-half*ampVisc + half)*(-tanh(slopeVisc*(   &
         &          r(:) - rStrat))**2 + 1)*(half*ampVisc + (half*ampVisc - half)&
         &          *tanh(slopeVisc*(r(:) - rStrat - thickStrat)) + half)/ampVisc&
         &          + slopeVisc*(half*ampVisc - half)*(-tanh(slopeVisc*(r(:) -   &
         &          rStrat - thickStrat))**2 + 1)*(half*ampVisc + (-half*ampVisc &
         &          + half)*tanh(slopeVisc*(r(:) - rStrat)) + half)/ampVisc)*(   &
         &          -tanh(slopeVisc*(r(:) - rStrat - thickStrat))**2 + 1)/((half*&
         &          ampVisc + (-half*ampVisc + half)*tanh(slopeVisc*(r(:) -      &
         &          rStrat)) + half)*(half*ampVisc + (half*ampVisc - half)*tanh( &
         &          slopeVisc*(r(:) - rStrat - thickStrat)) + half)**2)+ampVisc*(&
         &          2*slopeVisc**2*(-half*ampVisc + half)*(half*ampVisc - half)*(&
         &          -tanh(slopeVisc*(r(:) - rStrat))**2 + 1)*(-tanh(slopeVisc*(  &
         &          r(:) - rStrat - thickStrat))**2 + 1)/ampVisc - 2*slopeVisc**2&
         &          *(-half*ampVisc + half)*(-tanh(slopeVisc*(r(:) - rStrat))**2 &
         &          + 1)*(half*ampVisc + (half*ampVisc - half)*tanh(slopeVisc*(  &
         &          r(:) - rStrat - thickStrat)) + half)*tanh(slopeVisc*(r(:) -  &
         &          rStrat))/ampVisc - 2*slopeVisc**2*(half*ampVisc - half)*(    &
         &          -tanh(slopeVisc*(r(:) - rStrat - thickStrat))**2 + 1)*(half* &
         &          ampVisc + (-half*ampVisc + half)*tanh(slopeVisc*(r(:) -      &
         &          rStrat)) + half)*tanh(slopeVisc*(r(:) - rStrat - thickStrat))&
         &          /ampVisc)/((half*ampVisc + (-half*ampVisc + half)*tanh(      &
         &          slopeVisc*(r(:) - rStrat)) + half)*(half*ampVisc + (half*    &
         &          ampVisc - half)*tanh(slopeVisc*(r(:) - rStrat - thickStrat)) &
         &          + half))
      else if ( nVarVisc == 4 ) then ! Bottom stratified
         ampVisc = 10.0_cp
         slopeVisc = 30.0_cp
         if ( rStrat <= r_icb ) then
            visc(:) = one
         else
            visc(:)=(half*(ampVisc-one)*tanh(slopeVisc* &
            &       (r(:)-rStrat))+half*(ampVisc+one))/ampVisc
         end if
         dLvisc(:)=slopeVisc*(half*ampVisc-half)*(-tanh(slopeVisc*(r(:)-rStrat))**2&
         &         +one)/(half*ampVisc+(half*ampVisc-half)*tanh(slopeVisc*(r(:)-   &
         &         rStrat))+half)

         ddLvisc(:)=-slopeVisc**2*(half*ampVisc-half)**2*(-tanh(slopeVisc*(r(:)- &
         &          rStrat))**2+one)**2/(half*ampVisc+(half*ampVisc-half)*       &
         &          tanh(slopeVisc*(r(:)-rStrat))+half)**2-2*slopeVisc**2*(half* &
         &          ampVisc-half)*(-tanh(slopeVisc*(r(:)-rStrat))**2+one)*tanh(  &
         &          slopeVisc*(r(:)-rStrat))/(half*ampVisc+(half*ampVisc-half)*  &
         &          tanh(slopeVisc*(r(:)-rStrat))+half)
!My version of nVarVisc==3
!      else if ( nVarVisc == 3) then ! Profile nVarCond = 1
!         b =log(three)/con_FuncWidth
!         r0=con_radratio*r_cmb
!         s1=tanh(b*(r0-r_cmb))!
!         s2=tanh(b*(r0-r_icb))
!         a =(-one+con_LambdaOut)/(s1-s2)
!         c =(s1-s2*con_LambdaOut)/(s1-s2)
!         sigma   = a*tanh(b*(r0-r))+c
!         dsigma  =-a*b/cosh(b*(r0-r))
!         visc  =one/(a*tanh(b*(r0-r))+c)
!         dLvisc=(a*b/(cosh(b*(r0-r))))/(a*tanh(b*(r0-r))+c)
      else if ( nVarVisc == 5 ) then ! neutrino viscosity
         ! Neutrino viscosity
         ! Guilet et al, MNRAS 447, 3992-4003 (2015)
         ! Eq. (10)
         visc = temp0**2*rho0**(-2)
         visc = visc/visc(1)
         call get_dr(visc,dvisc,n_r_max,rscheme_oc)
         dLvisc(:)=dvisc(:)/visc(:)
         call get_dr(dLvisc,ddLvisc,n_r_max,rscheme_oc)
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

     integer :: n_r, i
     real(cp), allocatable :: coeffdEntropy0(:)
     real(cp) :: N_profile(n_r_max)

      if ( nVarEntropyGrad == 0 ) then ! Default: isentropic
         dEntropy0(:)=0.0_cp
         l_non_adia = .false.
      else if ( nVarEntropyGrad == 1  .or. rStrat >= r_cmb) then ! Takehiro
         if ( rStrat <= r_icb ) then
            dentropy0(:) = ampStrat
         else
            dentropy0(:) = -half*(ampStrat+one)*(one-tanh(slopeStrat*(r(:)-rStrat)))&
            &              + ampStrat
         end if
         l_non_adia = .true.
      else if ( nVarEntropyGrad == 2  .or. rStrat >= r_cmb) then ! Flat + linear
         if ( rStrat <= r_icb ) then
            dentropy0(:) = ampStrat
         else
            do n_r=1,n_r_max
               if ( r(n_r) <= rStrat  .or. rStrat >= r_cmb) then
                  dentropy0(n_r)=-one
               else
                  dentropy0(n_r)=(ampStrat+one)*(r(n_r)-r_cmb)/(r_cmb-rStrat) + &
                  &              ampStrat
               end if
            end do
         end if
         l_non_adia = .true.
      else if ( nVarEntropyGrad == 3 ) then ! SSL
         if ( rStrat <= r_icb  .or. rStrat >= r_cmb) then
            dentropy0(:) = -one
         else
            dentropy0(:) = 0.25_cp*(ampStrat+one)*(one+tanh(slopeStrat*(r(:)-rStrat)))&
            &              *(one-tanh(slopeStrat*(r(:)-rStrat-thickStrat)))           &
            &              - one
         end if
         l_non_adia = .true.
      else if ( nVarEntropyGrad == 4 ) then ! modified Takehiro
         if ( rStrat <= r_icb .or. rStrat >= r_cmb ) then
            dentropy0(:) = (r(:)**3-r_cmb**3)/(r_cmb**3 -r_icb**3)*(r_icb/r(:))**2
         else
            dentropy0(:) = half*(-ampStrat+(r(:)**3-r_cmb**3)/(r_cmb**3 -r_icb**3)* &
            &              (r_icb/r(:))**2)*(one-tanh(slopeStrat*(r(:)-rStrat))) +  &
            &              ampStrat
         end if
         l_non_adia = .true.
      else if ( nVarEntropyGrad == 5 ) then ! uniform volumic heat without strat
         dentropy0(:) = (r(:)**3-r_cmb**3)/(r_cmb**3 -r_icb**3)*(r_icb/r(:))**2
         l_non_adia = .true.
      else if (nVarEntropyGrad == 30) then !For g mode resonance with constant BV freq
         dentropy0(:) = -one/(g0+g1*r(:)/r_cmb+g2*r_cmb**2*or2(:))
         l_non_adia = .true.
      else if (nVarEntropyGrad == 31) then !NS freq model
         allocate(coeffdEntropy0(21))
         coeffdEntropy0(:)=[-8.21312802e-10_cp, -3.84090223e-08_cp, -7.39055209e-07_cp, -6.93324577e-06_cp,&
              & -2.03920934e-05_cp, 2.13283982e-04_cp, 2.23488684e-03_cp, 1.96336517e-03_cp, -1.05522774e-01_cp,&
              &-9.89930211e-01_cp, -4.93105883e+00_cp, -1.59444033e+01_cp, -3.48515802e+01_cp, -5.02358057e+01_cp, &
              & -4.18018822e+01_cp, -6.40566603e+00_cp, 2.80473923e+01_cp, 3.48563797e+01_cp, &
              & 2.06455828e+01_cp, 7.85773117e+00_cp, 3.27009642e-01_cp]
         N_profile(:)=0.0_cp
         rStrat=0.8338751096437745_cp
         do n_r=1,n_r_max
            if (r(n_r)>rStrat) then
               N_profile(n_r)=(3.27412712_cp*log(r(n_r)/r_cmb))
            else if (r(n_r)>0.003_cp) then
               do i=1,21
                  N_profile(n_r) = N_profile(n_r)+coeffdEntropy0(i)*(log(r(n_r)/r_cmb))**(21-i)
               end do
            else
               N_profile(n_r)=-1e50_cp
            end if
         end do
         dentropy0(:)=-(exp(N_profile(:)))**2/(g0+g1*r(:)/r_cmb+g2*r_cmb**2*or2(:))!to keep BV freq profile cst
         if (l_full_sphere) then
            dentropy0(n_r_max)=zero
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
      integer :: n_r,info,n_r_out
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

      call prepare_mat(workMat,n_r_max,n_r_max,workPivot,info)

      if ( info /= 0 ) then
         call abortRun('! Singular Matrix in getBackground!')
      end if

      do n_r=2,n_r_max
         rhs(n_r)=input(n_r)
      end do
      rhs(1)=boundaryVal

      do n_r=1,n_r_max
         rhs(n_r)=rhs(n_r)*workMat_fac(n_r,1)
      end do

      !-- Solve for s0:
      call solve_mat(workMat,n_r_max,n_r_max,workPivot,rhs)

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
   subroutine polynomialBackground(coeffDens,coeffTemp,coeffGrav)
      !
      ! This subroutine allows to calculate a reference state based on an input
      ! polynomial function.
      !

      !-- Input variables
      real(cp),           intent(in) :: coeffDens(:)
      real(cp),           intent(in) :: coeffTemp(:)
      real(cp), optional, intent(in) :: coeffGrav(:)

      !-- Local variables
      real(cp) :: rrOcmb(n_r_max),gravFit(n_r_max)
      real(cp) :: drho0(n_r_max),dtemp0(n_r_max)
      integer :: nGrav,nDens,nTemp,i

      rrOcmb(:) = r(:)*r_cut_model/r_cmb

      !-- Assemble gravity profile
      if ( present(coeffGrav) ) then
         nGrav=size(coeffGrav)
         gravFit(:)=0.0_cp
         do i=1,nGrav
            gravFit(:) = gravFit(:)+coeffGrav(i)*rrOcmb(:)**(i-1)
         end do
      else
         gravFit(:)=four*rrOcmb(:)-three*rrOcmb(:)**2
      end if

      ! Set to zero initially
      nDens = size(coeffDens)
      nTemp = size(coeffTemp)
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
      DissNb   =alpha0(1)
      alpha0(:)=alpha0(:)/alpha0(1)

      ! Adiabatic: buoyancy term is linked to the temperature gradient

      !       dT
      !      ---- =  -Di * alpha_T * T * grav
      !       dr
      rgrav(:)=-dtemp0(:)/DissNb

      call get_dr(rho0,drho0,n_r_max,rscheme_oc)
      beta(:)=drho0(:)/rho0(:)
      call get_dr(beta,dbeta,n_r_max,rscheme_oc)
      call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)
      call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
      call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
      dLalpha0(:)=dLalpha0(:)/alpha0(:) ! d log (alpha) / dr
      call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
      dLtemp0(:)=dtemp0(:)/temp0(:)
      call get_dr(dLtemp0,ddLtemp0,n_r_max,rscheme_oc)
      dentropy0(:)=0.0_cp

      !- \Gamma = 1/\rho/c_v ( \partial p/\partial T)_\rho = (d\ln T/d \ln \rho)_s
      ogrun(:)=beta(:)/dLtemp0(:)
      GrunNb  =one/ogrun(1)
      strat   =log(rho0(n_r_max)/rho0(1))
      polind  =ogrun(1)

   end subroutine polynomialBackground
!------------------------------------------------------------------------------
end module radial_functions
