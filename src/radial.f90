module radial_functions
   !
   !  This module initiates all the radial functions (transport properties, density,
   !  temperature, cheb transforms, etc.)
   !

   use iso_fortran_env, only: output_unit
   use truncation, only: n_r_max, n_cheb_max, n_r_ic_max, fd_ratio, rcut_l, &
       &                 fd_stretch, fd_order, fd_order_bound, l_max, n_cheb_ic_max
   use algebra, only: prepare_mat, solve_mat
   use constants, only: sq4pi, one, two, three, four, half, pi
   use physical_parameters
   use logic, only: l_mag, l_cond_ic, l_heat, l_anelastic_liquid,  &
       &            l_isothermal, l_anel, l_non_adia, l_centrifuge,&
       &            l_temperature_diff, l_single_matrix, l_var_l,  &
       &            l_finite_diff, l_newmap, l_full_sphere,        &
       &            l_chemical_conv
   use radial_data, only: nRstart, nRstop
   use chebyshev_polynoms_mod ! Everything is needed
   use cosine_transform_odd
   use cosine_transform_even
   use radial_scheme, only: type_rscheme
   use chebyshev, only: type_cheb_odd
   use finite_differences, only: type_fd
   use radial_der, only: get_dr
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use useful, only: logWrite, abortRun
   use parallel_mod, only: rank
   use output_data, only: tag
   use num_param, only: alph1, alph2
   use special, only: l_curr, fac_loop

   implicit none

   private

   !-- arrays depending on r:
#ifdef WITH_OMP_GPU
   !$omp declare target (r, r_ic, O_r_ic, O_r_ic2, or1, or2, or3, or4, &
   !$omp&                otemp1, rho0, temp0, dLtemp0, dentropy0, ddLtemp0, &
   !$omp&                orho1, orho2, beta, dbeta, ddbeta, alpha0, dLalpha0, ddLalpha0, &
   !$omp&                rgrav, ogrun, &
   !$omp&                lambda, dLlambda, jVarCon, sigma, kappa, dLkappa, &
   !$omp&                visc, dLvisc, ddLvisc, epscProf, divKtemp0, l_R, &
   !$omp&                cheb_ic, dcheb_ic, d2cheb_ic, cheb_int_ic, dr_top_ic, cheb_int, &
   !$omp&                dxicond) !-- Note: Compiler does not accept for rscheme_oc
#endif
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
   integer :: nDd_costf1_ic                      ! Radii for transform
   integer :: nDi_costf2_ic                      ! Radii for transform
   integer :: nDd_costf2_ic                      ! Radii for transform

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

   subroutine initialize_radial_functions()
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
      bytes_allocated = bytes_allocated+(22*n_r_max+3*n_r_ic_max)*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: r, r_ic, O_r_ic, O_r_ic2, or1, or2, or3, or4, &
      !$omp&                             otemp1, rho0, temp0, dLtemp0, dentropy0, ddLtemp0, &
      !$omp&                             orho1, orho2, beta, dbeta, ddbeta, alpha0, dLalpha0, ddLalpha0, &
      !$omp&                             rgrav, ogrun)
      gpu_bytes_allocated = gpu_bytes_allocated+(21*n_r_max+3*n_r_ic_max)*SIZEOF_DEF_REAL !-- d2temp0 is not on GPU
#endif

      if ( l_chemical_conv ) then
         allocate( dxicond(n_r_max) )
         dxicond(:)=0.0_cp
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc:dxicond)
         gpu_bytes_allocated = gpu_bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
      end if

      allocate( lambda(n_r_max),dLlambda(n_r_max),jVarCon(n_r_max) )
      allocate( sigma(n_r_max),kappa(n_r_max),dLkappa(n_r_max) )
      allocate( visc(n_r_max),dLvisc(n_r_max),ddLvisc(n_r_max) )
      allocate( epscProf(n_r_max),divKtemp0(n_r_max) )
      bytes_allocated = bytes_allocated + 11*n_r_max*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: lambda, dLlambda, jVarCon, sigma, kappa, dLkappa, &
      !$omp&                             visc, dLvisc, ddLvisc, epscProf, divKtemp0)
      gpu_bytes_allocated = gpu_bytes_allocated + 11*n_r_max*SIZEOF_DEF_REAL
#endif

      !allocate ( l_R(nRstart:nRstop) )
      !bytes_allocated = bytes_allocated +(nRstop-nRstart+1)*SIZEOF_INTEGER
      allocate ( l_R(1:n_r_max) )
      bytes_allocated = bytes_allocated +n_r_max*SIZEOF_INTEGER
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: l_R)
      gpu_bytes_allocated = gpu_bytes_allocated +n_r_max*SIZEOF_INTEGER
#endif

      if ( .not. l_full_sphere ) then
         nDd_costf1_ic=2*n_r_ic_max+5
         nDi_costf2_ic=2*n_r_ic_max
         nDd_costf2_ic=2*n_r_ic_max+n_r_ic_max/2+5

         allocate( cheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( dcheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( d2cheb_ic(n_r_ic_max,n_r_ic_max) )
         allocate( cheb_int_ic(n_r_ic_max) )
         bytes_allocated = bytes_allocated + &
         &                 (3*n_r_ic_max*n_r_ic_max+n_r_ic_max)*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: cheb_ic, dcheb_ic, d2cheb_ic, cheb_int_ic)
         gpu_bytes_allocated = gpu_bytes_allocated + &
         &                 (3*n_r_ic_max*n_r_ic_max+n_r_ic_max)*SIZEOF_DEF_REAL
#endif

         call chebt_ic%initialize(n_r_ic_max,n_cheb_ic_max,nDd_costf1_ic)

         allocate ( dr_top_ic(n_r_ic_max) )
         bytes_allocated = bytes_allocated+n_r_ic_max*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: dr_top_ic)
         gpu_bytes_allocated = gpu_bytes_allocated+n_r_ic_max*SIZEOF_DEF_REAL
#endif
      end if

      if ( .not. l_finite_diff ) then

         allocate( cheb_int(n_r_max) )         ! array for cheb integrals !
         bytes_allocated = bytes_allocated + n_r_max*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: cheb_int)
         gpu_bytes_allocated = gpu_bytes_allocated + n_r_max*SIZEOF_DEF_REAL
#endif

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

#ifdef WITH_OMP_GPU
      call rscheme_oc%initialize(n_r_max,n_in,n_in_2,.true.)
      !$omp target enter data map(alloc: rscheme_oc)
#else
      call rscheme_oc%initialize(n_r_max,n_in,n_in_2)
#endif

   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine finalize_radial_functions()
      !
      ! Memory deallocation of radial functions
      !

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: l_R, r, r_ic, O_r_ic, O_r_ic2, or1, or2, or3, or4, &
      !$omp&                             otemp1, rho0, temp0, dLtemp0, dentropy0, &
      !$omp&                             ddLtemp0, orho1, orho2, beta, dbeta, ddbeta, alpha0, &
      !$omp&                             ddLalpha0, dLalpha0, rgrav, ogrun, &
      !$omp&                             lambda, dLlambda, jVarCon, sigma, kappa, dLkappa, &
      !$omp&                             visc, dLvisc, ddLvisc, epscProf, divKtemp0)

      if ( l_curr ) then
         !$omp target exit data map(delete : fac_loop)
      end if
#endif

      deallocate( l_R )
      deallocate( r, r_ic, O_r_ic, O_r_ic2, or1, or2, or3, or4 )
      deallocate( otemp1, rho0, temp0, dLtemp0, d2temp0, dentropy0 )
      deallocate( ddLtemp0, orho1, orho2, beta, dbeta, ddbeta, alpha0 )
      deallocate( ddLalpha0, dLalpha0, rgrav, ogrun )
      deallocate( lambda, dLlambda, jVarCon, sigma, kappa, dLkappa )
      deallocate( visc, dLvisc, ddLvisc, epscProf, divKtemp0 )

      if ( l_chemical_conv ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dxicond)
#endif
         deallocate(dxicond)
      end if

      if ( .not. l_full_sphere ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: dr_top_ic, cheb_ic, dcheb_ic, d2cheb_ic, cheb_int_ic)
#endif
         deallocate( dr_top_ic )
         deallocate( cheb_ic, dcheb_ic, d2cheb_ic, cheb_int_ic )
         call chebt_ic%finalize()
         if ( n_r_ic_max > 0 .and. l_cond_ic ) call chebt_ic_even%finalize()
      end if

      if ( .not. l_finite_diff ) then
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: cheb_int)
#endif
         deallocate( cheb_int )
      end if

#ifdef WITH_OMP_GPU
      !$omp target exit data map(release : rscheme_oc)
      call rscheme_oc%finalize(.true.)
#else
      call rscheme_oc%finalize()
#endif

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

      real(cp), allocatable :: coeffDens(:), coeffTemp(:), coeffAlpha(:)
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
      r_cmb=one/(one-radratio)
      r_icb=r_cmb-one

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

         ! the shell can't be thinner than eta=0.38947 and r_cut_model <= 0.95, because
         ! the fits don't work beyond that. The Preising et al. 2023 model
         ! includes a region of strong ! Helium concentration at the
         ! bottom which is probably stably stratified and has
         ! been cut out from this fit.

         if (l_non_adia) then
            rrOcmb(:) = r(:)*r_cut_model/r_cmb
            allocate(coeffTemp(8),coeffDens(6),coeffGrav(5),coeffAlpha(8))

            coeffTemp= [ -77880.376328818_cp, 1089237.928095523_cp,    &
            &           -5718555.197998112_cp, 16027297.768664116_cp,  &
            &           -26159709.956658602_cp, 24913613.036177561_cp, &
            &           -12819667.849508610_cp, 2746322.598433222_cp ]
            coeffDens= [ 8.813043492_cp, -51.578717073_cp, 145.416402036_cp, &
            &           -211.845381533_cp, 152.131339032_cp, -42.964602602_cp]
            coeffGrav= [ 115.392383256_cp, -483.892696422_cp, 935.515330783_cp,&
            &           -826.605750144_cp, 270.984339687_cp ]
            coeffAlpha= [ -871.964023446_cp, 10634.774082757_cp, -55010.331729823_cp,   &
            &            154178.756043949_cp, -252573.849597133_cp, 241773.428704681_cp,&
            &            -125275.544905827_cp, 27136.257157590_cp ]

            alpha0(:)=0.0_cp
            temp0(:) =0.0_cp
            rgrav(:) =0.0_cp
            rho0(:)  =0.0_cp

            do i=1,8
               alpha0(:) = alpha0(:)+coeffAlpha(i)*rrOcmb(:)**(i-1)
               temp0(:)  = temp0(:) +coeffTemp(i) *rrOcmb(:)**(i-1)
            end do

            do i=1,5
               rgrav(:)  = rgrav(:) +coeffGrav(i) *rrOcmb(:)**(i-1)
            end do

            do i=1,6
               rho0(:)   = rho0(:) +coeffDens(i) *rrOcmb(:)**(i-1)
            end do

            alpha0(:)=exp(alpha0(:)) ! Polynomial fit was on ln(alpha0)

            DissNb   =alpha0(1)*rgrav(1)*(rrOcmb(1)-rrOcmb(n_r_max))*5.8232e7_cp &
            &         /1.73e4_cp ! 1.73e4 is cp 5.8232e7_cp is R_S

            ThExpNb  =alpha0(1)*temp0(1)
            alpha0(:)=alpha0(:)/alpha0(1)
            rgrav(:) =rgrav(:)/rgrav(1)
            temp0(:) =temp0(:)/temp0(1)
            rho0(:)  =rho0(:)/rho0(1)

            strat   =log(rho0(n_r_max)/rho0(1)) ! Just for printing

            !-- Radial profile for the Grüneisen parameter (from Preising et al. 2023)
            ogrun(:) = ( -0.25_cp * (0.7_cp - 0.14_cp)*(             &
            &             1.0_cp+tanh(50.0_cp*(rrOcmb(:)-0.65_cp)))  &
            &           *(1.0_cp-tanh(50.0_cp*(rrOcmb(:)-0.8_cp)))) + 0.7_cp
            GrunNb = ogrun(1)
            ogrun(:) = one/ogrun(:)
            polind   = ogrun(1) ! Just for the print
            ogrun(:) = ogrun(:)/ogrun(1)

            ! The final stuff is always required
            call get_dr(rho0,drho0,n_r_max,rscheme_oc)
            beta(:) = drho0(:)/rho0(:)
            call get_dr(beta,dbeta,n_r_max,rscheme_oc)
            call get_dr(dbeta,ddbeta,n_r_max,rscheme_oc)

            call get_dr(temp0,dtemp0,n_r_max,rscheme_oc)
            call get_dr(dtemp0,d2temp0,n_r_max,rscheme_oc)
            call get_dr(alpha0,dLalpha0,n_r_max,rscheme_oc)
            dLalpha0=dLalpha0/alpha0 ! d log (alpha) / dr
            call get_dr(dLalpha0,ddLalpha0,n_r_max,rscheme_oc)
            dLtemp0 = dtemp0/temp0
            ddLtemp0 =-(dtemp0/temp0)**2+d2temp0/temp0

            !-- Multiply the gravity by alpha0 and temp0
            rgrav(:)=rgrav(:)*alpha0(:)*temp0(:)

            deallocate(coeffTemp,coeffDens,coeffGrav,coeffAlpha)
         else
            allocate(coeffTemp(8),coeffDens(6))
            coeffTemp= [ -77880.376328818_cp, 1089237.928095523_cp,      &
            &            -5718555.197998112_cp, 16027297.768664116_cp,   &
            &            -26159709.956658602_cp, 24913613.036177561_cp,  &
            &            -12819667.849508610_cp, 2746322.598433222_cp ]
            coeffDens= [ 8.813043492_cp, -51.578717073_cp, 145.416402036_cp, &
            &           -211.845381533_cp, 152.131339032_cp, -42.964602602_cp]

            call polynomialBackground(coeffDens,coeffTemp)
            deallocate( coeffDens, coeffTemp)
         end if

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
         rrOcmb(:) = r(:)*r_cut_model/r_cmb

         allocate( coeffAlpha(12), coeffTemp(12), coeffGrav(12), coeffGrun(12),&
         &         coeffDens(12) )

         coeffAlpha=[1.07280013e+00_cp, -7.35271452e-02_cp,  1.35373501e+00_cp,&
         &          -2.70033254e+01_cp,  1.34218482e+02_cp, -1.85203678e+02_cp,&
         &          -4.03705669e+02_cp,  1.65171816e+03_cp, -1.91235980e+03_cp,&
         &           4.33699991e+02_cp,  6.74102545e+02_cp, -3.67415688e+02_cp]

         coeffTemp=[1.71233163e+01_cp,  2.16376708e-01_cp, -2.33617981e+01_cp, &
         &          2.95986199e+02_cp, -3.15614535e+03_cp,  1.80265801e+04_cp, &
         &         -5.99840212e+04_cp,  1.23788205e+05_cp, -1.61112426e+05_cp, &
         &          1.28833197e+05_cp, -5.78439543e+04_cp,  1.11724125e+04_cp]

         coeffGrav=[1.31049693e+02_cp,  1.04893208e+06_cp,  1.93822208e+06_cp, &
         &         -5.15044761e+07_cp,  4.41798276e+08_cp, -2.66511461e+09_cp, &
         &          1.02063273e+10_cp, -2.43712757e+10_cp,  3.64094408e+10_cp, &
         &         -3.32103850e+10_cp,  1.69572548e+10_cp, -3.72150208e+09_cp]

         coeffGrun=[1.57026325e+00_cp, -4.46247381e-02_cp,  3.03734351e-01_cp, &
         &         -8.28253239e+00_cp, -2.63234309e+01_cp,  6.05835596e+02_cp, &
         &         -2.89779916e+03_cp,  6.90494749e+03_cp, -9.29010895e+03_cp, &
         &          7.04974613e+03_cp, -2.73663642e+03_cp,  3.98015351e+02_cp]

         coeffDens=[3.05336249e+00_cp, -2.52357616e-01_cp, -3.47786718e-01_cp, &
         &         -3.82246141e+02_cp,  4.44024160e+03_cp, -2.80734664e+04_cp, &
         &          1.02942741e+05_cp, -2.30659612e+05_cp,  3.22111291e+05_cp, &
         &         -2.74248402e+05_cp,  1.30456564e+05_cp, -2.66069272e+04_cp]

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
         r_ic(1:n_r_ic_max-1)   =r_ic_2(1:n_r_ic_max-1)
         O_r_ic(1:n_r_ic_max-1) =one/r_ic(1:n_r_ic_max-1)
         O_r_ic2(1:n_r_ic_max-1)=O_r_ic(1:n_r_ic_max-1)*O_r_ic(1:n_r_ic_max-1)
         r_ic(n_r_ic_max)   =0.0_cp
         O_r_ic(n_r_ic_max) =0.0_cp
         O_r_ic2(n_r_ic_max)=0.0_cp

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
   subroutine transportProperties()
      !
      ! Calculates the transport properties: electrical conductivity,
      ! kinematic viscosity and thermal conductivity.
      !
      real(cp) :: a,b,c,s1,s2,r0,a0,a1,a2,a3,a4,a5
      real(cp) :: dsigma0, ampVisc, ampKap, slopeVisc, slopeKap
      real(cp) :: dvisc(n_r_max), dkappa(n_r_max), dsigma(n_r_max)
      real(cp) :: rrOcmb(n_r_max), kcond(n_r_max)
      !real(cp) :: condBot(n_r_max), condTop(n_r_max), func(n_r_max)

      !-- Variable conductivity:
      if ( imagcon == -10 ) then
         nVarCond=1
         lambda(:)  =r(:)**5.0_cp
         sigma(:)   =one/lambda(:)
         dLlambda(:)=5.0_cp/r(:)
      else if ( l_mag ) then
         select case(nVarCond)

            case(1)
               b =log(three)/con_FuncWidth
               r0=con_radratio*r_cmb
               s1=tanh(b*(r0-r_cmb))
               s2=tanh(b*(r0-r_icb))
               a =(-one+con_LambdaOut)/(s1-s2)
               c =(s1-s2*con_LambdaOut)/(s1-s2)
               sigma(:)   = a*tanh(b*(r0-r(:)))+c
               dsigma(:)  =-a*b/cosh(b*(r0-r(:)))
               lambda(:)  =one/sigma(:)
               dLlambda(:)=-dsigma(:)/sigma(:)

            case(2) ! Two-branches solution from Gomez-Perez
               r0=con_radratio*r_cmb
               !------ Find the grid point closest to r0:
               r0=r(minloc(abs(r(:)-r0),dim=1))
               dsigma0=(con_LambdaMatch-1)*con_DecRate/(r0-r_icb)
               where ( r < r0 )
                  sigma=one+(con_LambdaMatch-1)*((r-r_icb)/(r0-r_icb))**con_DecRate
                  dsigma=dsigma0*((r-r_icb)/(r0-r_icb))**(con_DecRate-1)
               else where
                  sigma=con_LambdaMatch*exp(dsigma0/con_LambdaMatch*(r-r0))
                  dsigma=dsigma0*exp(dsigma0/con_LambdaMatch*(r-r0))
               end where
               lambda(:)  =one/sigma(:)
               dLlambda(:)=-dsigma(:)/sigma(:)

            case(3) ! Magnetic diff propto 1/rho
               lambda(:)=rho0(n_r_max)/rho0(:)
               sigma(:)=one/lambda(:)
               call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
               dLlambda(:)=dsigma(:)/lambda(:)

            case(4) ! Profile related to a power-law of density
               lambda=(rho0(:)/rho0(n_r_max))**difExp
               sigma(:)=one/lambda(:)
               call get_dr(lambda,dsigma,n_r_max,rscheme_oc)
               dLlambda(:)=dsigma(:)/lambda(:)

            case default ! Constant diffusivity
               lambda(:)  =one
               sigma(:)   =one
               dLlambda(:)=0.0_cp

         end select
      end if

      !-- Variable thermal diffusivity
      if ( l_heat ) then

         select case(nVarDiff)

            case(1) ! Constant conductivity
               ! kappa(:)=one/rho0(:) Denise Tortorella's version
               kappa(:)=rho0(n_r_max)/rho0(:)
               call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
               dLkappa(:)=dkappa(:)/kappa(:)

            case(2) ! Related to a power-law of the density profile
               kappa(:)=(rho0(:)/rho0(n_r_max))**difExp
               call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
               dLkappa(:)=dkappa(:)/kappa(:)

            case(3) ! Polynomial fit to a model
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
               rrOcmb(:) = r(:)/r_cmb*r_cut_model
               kappa(:)= a5 + a4*rrOcmb(:) + a3*rrOcmb(:)**2 + a2*rrOcmb(:)**3 &
               &            + a1*rrOcmb(:)**4 + a0*rrOcmb(:)**5
               kappa(:)=kappa(:)/kappa(1) ! normalize by the value at the outer boundary
               call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
               dLkappa(:)=dkappa(:)/kappa(:)

            case(4) ! Earth case
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
               kcond(:)=(one-0.4469_cp*(r(:)/r_cmb)**2)/(one-0.4469_cp)
               kcond(:)=kcond(:)/kcond(1)
               kappa(:)=kcond(:)/rho0
               call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
               dLkappa(:)=dkappa(:)/kappa(:)

            case(5) ! Thermal background equilibrium
               kcond(:)=one/(r(:)*r(:)*dLtemp0(:)*temp0(:))
               kcond(:)=kcond(:)/kcond(1)
               kappa(:)=kcond(:)/rho0(:)
               call get_dr(kappa,dkappa,n_r_max,rscheme_oc)
               dLkappa(:)=dkappa(:)/kappa(:)

            case(6) ! Jump in the stratified layer
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

            case(7) ! Bottom stratified
               ampKap = 10.0_cp
               slopeKap = 30.0_cp
               if ( rStrat <= r_icb ) then
                  kappa(:)=one
               else
                  kappa(:)=(half*(ampKap-one)*tanh(slopeKap* &
                  &       (r(:)-rStrat))+half*(ampKap+one))/ampKap
               end if
               dLkappa(:)=slopeKap*(half*ampKap-half)*(-tanh(slopeKap*(r(:)-rStrat))**2&
               &         +one)/(half*ampKap+(half*ampKap-half)*tanh(slopeKap*(r(:)-    &
               &         rStrat))+half)

            case default
               kappa(:)  =one
               dLkappa(:)=0.0_cp

         end select

      end if

      !-- Eps profiles
      !-- The remaining division by rho will happen in updateS.f90
      select case(nVarEps)

         case(1) ! rho*eps in the RHS
            if ( l_anelastic_liquid ) then
               epscProf(:)=rho0(:)
            else
               epscProf(:)=rho0(:)*otemp1(:)
            end if

         case(2) ! rho*temp*eps in the RHS
            if ( l_anelastic_liquid ) then
               epscProf(:)=rho0(:)*temp0(:)
            else
               epscProf(:)=rho0(:)
            end if
         case(3) ! eps*rho**2*temp**(-3)*exp(-Bn/T) in the RHS
            if ( l_anelastic_liquid ) then
               epscProf(:)=rho0(:)**2*temp0(:)**(-3)*exp(-Bn/temp0(:))
            else
               epscProf(:)=rho0(:)**2*temp0(:)**(-4)*exp(-Bn/temp0(:))
            end if

         case default ! eps is constant
            if ( l_anelastic_liquid ) then
               epscProf(:)=one
            else
               epscProf(:)=otemp1(:)
            end if

      end select

      !-- Variable viscosity
      select case(nVarVisc)

         case(1) ! Constant dynamic viscosity
            visc(:)=rho0(n_r_max)/rho0(:)
            call get_dr(visc,dvisc,n_r_max,rscheme_oc)
            dLvisc(:)=dvisc(:)/visc(:)
            call get_dr(dLvisc,ddLvisc,n_r_max,rscheme_oc)

         case(2) ! Profile based on a power-law of density
            visc(:)=(rho0/rho0(n_r_max))**difExp
            call get_dr(visc,dvisc,n_r_max,rscheme_oc)
            dLvisc(:)=dvisc(:)/visc(:)
            call get_dr(dLvisc,ddLvisc,n_r_max,rscheme_oc)

         case(3) ! Jump in the stratified layer
            ampVisc  =10.0_cp
            slopeVisc=100.0_cp
            if ( rStrat <= r_icb ) then
               visc(:)=one
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

         case(4) ! Bottom stratified
            ampVisc  =10.0_cp
            slopeVisc=30.0_cp
            if ( rStrat <= r_icb ) then
               visc(:)=one
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

         case default ! default: constant kinematic viscosity
            visc(:)   =one
            dLvisc(:) =0.0_cp
            ddLvisc(:)=0.0_cp

      end select

      if ( l_anelastic_liquid .or. l_non_adia ) then
         divKtemp0=rho0*kappa*(d2temp0+(beta+dLkappa+two*or1)*temp0*dLtemp0)*sq4pi
      else
         divKtemp0=0.0_cp
      end if

   end subroutine transportProperties
!------------------------------------------------------------------------------
   subroutine getEntropyGradient()
      !
      ! This subroutine allows to calculate the background entropy gradient
      ! in case stable stratification is required
      !

      integer :: k

      select case(nVarEntropyGrad)

         case(1) ! Takehiro
            if ( rStrat <= r_icb .or. rStrat >= r_cmb ) then
               dentropy0(:) = ampStrat
            else
               dentropy0(:) = -half*(ampStrat+one)*(one-tanh(slopeStrat*(r(:)-rStrat)))&
               &              + ampStrat
            end if
            l_non_adia = .true.

         case(2) ! Flat + linear
            if ( rStrat <= r_icb .or. rStrat >= r_cmb ) then
               dentropy0(:) = ampStrat
            else
               where ( r <= rStrat )
                  dentropy0 = -one
               else where
                  dentropy0 = (ampStrat+one)*(r-r_cmb)/(r_cmb-rStrat) + ampStrat
               end where
            end if
            l_non_adia = .true.

         case(3) ! SSL
            if ( rStrat <= r_icb  .or. rStrat >= r_cmb) then
               dentropy0(:) = -one
            else
               dentropy0(:) = 0.25_cp*(ampStrat+one)*(one+tanh(slopeStrat*(r(:)-rStrat)))&
               &              *(one-tanh(slopeStrat*(r(:)-rStrat-thickStrat)))           &
               &              - one
            end if
            l_non_adia = .true.

         case(4) ! modified Takehiro
            if ( rStrat <= r_icb .or. rStrat >= r_cmb ) then
               dentropy0(:) = (r(:)**3-r_cmb**3)/(r_cmb**3-r_icb**3)*(r_icb/r(:))**2
            else
               dentropy0(:) = half*(-ampStrat+(r(:)**3-r_cmb**3)/(r_cmb**3-r_icb**3)* &
               &              (r_icb/r(:))**2)*(one-tanh(slopeStrat*(r(:)-rStrat))) + &
               &              ampStrat
            end if
            l_non_adia = .true.

         case(5) ! uniform volumic heat without strat
            dentropy0(:) = (r(:)**3-r_cmb**3)/(r_cmb**3-r_icb**3)*(r_icb/r(:))**2
            l_non_adia = .true.

         case(6) ! Chemical heating + linear stratification
            if ( rStrat <= r_icb .or. rStrat >= r_cmb ) then
               dentropy0(:) = (r(:)**3-r_cmb**3)/(r_cmb**3-r_icb**3)*(r_icb/r(:))**2
            else
               dentropy0(:) = (r(:)**3-r_cmb**3)/(r_cmb**3-r_icb**3)*(r_icb/r(:))**2
               where ( r >= rStrat )
                  dentropy0 = dentropy0 + ampStrat*(r-rStrat)/(r_cmb-rStrat)
               end where
            end if
            l_non_adia = .true.

         case(7) ! N stratified layers, given by arrays
            dentropy0(:) = -one

            ! Check if bottom SSL exists
            if ( ampStrat_arr(1) > 0.0_cp ) then
               dentropy0(:) = dentropy0(:) + (ampStrat_arr(1)+one) *         &
               &             (one - half * (one-tanh(slopeStrat_arr(1) *  &
               &                               (rStrat_arr(1) - r(:)))) )
            end if

            ! Subsequent ones
            do k=2,nSSLmax
               dentropy0(:) = dentropy0(:) + ( 0.25_cp * (ampStrat_arr(k)+one)* &
               &                 (one+tanh(slopeStrat_arr(k) * (r(:)-rStrat_arr(k)))) *  &
               &                 (one - tanh(slopeStrat_arr(k) * (r(:) -   &
               &                 rStrat_arr(k) - thickStrat_arr(k))) ) )
            end do

            l_non_adia = .true.

         case default ! Default: isentropic
            dentropy0(:)=0.0_cp
            l_non_adia = .false.

      end select

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
     workMat(:,1)      =rscheme_oc%boundary_fac*workMat(:,1)
     workMat(:,n_r_max)=rscheme_oc%boundary_fac*workMat(:,n_r_max)

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

      rhs(2:n_r_max)=input(2:n_r_max)
      rhs(1)=boundaryVal

      rhs(:)=rhs(:)*workMat_fac(:,1)

      !-- Solve for s0:
      call solve_mat(workMat,n_r_max,n_r_max,workPivot,rhs)

      !-- Copy result to s0:
      output(:)=rhs(:)*workMat_fac(:,2)

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
