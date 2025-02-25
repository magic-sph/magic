module horizontal_data
   !
   !  Module containing functions depending on longitude
   !  and latitude plus help arrays depending on degree and order
   !

   use truncation, only: l_max, n_theta_max, n_phi_max, nlat_padded, &
       &                 lm_max, n_m_max, minc, m_min, m_max, l_axi
   use radial_functions, only: r_cmb
   use physical_parameters, only: ek
   use num_param, only: difeta, difnu, difkap, ldif, ldifexp, difchem
   use blocking, only: lm2l, lm2m
   use logic, only: l_non_rot, l_scramble_theta
   use plms_theta, only: plm_theta
   use fft
   use constants, only: pi, zero, one, two, half
   use precision_mod
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   !-- Arrays depending on theta (colatitude):
   integer, public, allocatable :: n_theta_cal2ord(:)
   integer, public, allocatable :: n_theta_ord2cal(:)
   real(cp), public, allocatable :: theta_ord(:)       ! Gauss points (unscrambled)
   real(cp), public, allocatable :: sinTheta_E2(:)     ! :math:`\sin^2\theta`
   real(cp), public, allocatable :: O_sin_theta(:)     ! :math:`1/\sin\theta`
   real(cp), public, allocatable :: O_sin_theta_E2(:)  ! :math:`1/\sin^2\theta`
   real(cp), public, allocatable :: cosn_theta_E2(:)   ! :math:`\cos\theta/\sin^2\theta`
   real(cp), public, allocatable :: sinTheta(:)        ! :math:`\sin\theta`
   real(cp), public, allocatable :: cosTheta(:)        ! :math:`\cos\theta`

   !-- Phi (longitude)
   real(cp), public, allocatable :: phi(:)

   !-- Legendres:
   real(cp), public, allocatable :: gauss(:)
   real(cp), public, allocatable :: dPl0Eq(:)

   !-- Arrays depending on l and m:
   complex(cp), public, allocatable :: dPhi(:) ! :math:`\mathrm{i} m`
   real(cp), public, allocatable :: dLh(:)     ! :math:`\ell(\ell+1)`
   real(cp), public, allocatable :: dTheta1S(:),dTheta1A(:)
   real(cp), public, allocatable :: dTheta2S(:),dTheta2A(:)
   real(cp), public, allocatable :: dTheta3S(:),dTheta3A(:)
   real(cp), public, allocatable :: dTheta4S(:),dTheta4A(:)
   real(cp), public, allocatable :: hdif_B(:),hdif_V(:),hdif_S(:),hdif_Xi(:)
   
   public :: initialize_horizontal_data, horizontal, finalize_horizontal_data, &
   &         gauleg

contains

   subroutine initialize_horizontal_data()
      !
      ! Memory allocation of horizontal functions
      !

      allocate( n_theta_cal2ord(n_theta_max) )
      allocate( n_theta_ord2cal(n_theta_max) )
      allocate( theta_ord(n_theta_max) )
      allocate( sinTheta_E2(nlat_padded) )
      allocate( O_sin_theta(nlat_padded) )
      allocate( O_sin_theta_E2(nlat_padded) )
      allocate( sinTheta(nlat_padded) )
      allocate( cosn_theta_E2(nlat_padded) )
      allocate( cosTheta(nlat_padded) )
      bytes_allocated = bytes_allocated+2*n_theta_max*SIZEOF_INTEGER+&
      &                 (n_theta_max+6*nlat_padded)*SIZEOF_DEF_REAL
      O_sin_theta(:)   =0.0_cp
      O_sin_theta_E2(:)=0.0_cp
      sinTheta(:)      =0.0_cp
      cosn_theta_E2(:) =0.0_cp
      cosTheta(:)      =0.0_cp
      sinTheta_E2(:)   =0.0_cp

      !-- Phi (longitude)
      allocate( phi(n_phi_max) )
      bytes_allocated = bytes_allocated+n_phi_max*SIZEOF_DEF_REAL

      !-- Legendres:
      allocate( gauss(n_theta_max) )
      allocate( dPl0Eq(l_max+1) )
      bytes_allocated = bytes_allocated+(n_theta_max+l_max+1)*SIZEOF_DEF_REAL

      !-- Arrays depending on l and m:
      allocate( dPhi(lm_max), dLh(lm_max) )
      allocate( dTheta1S(lm_max),dTheta1A(lm_max) )
      allocate( dTheta2S(lm_max),dTheta2A(lm_max) )
      allocate( dTheta3S(lm_max),dTheta3A(lm_max) )
      allocate( dTheta4S(lm_max),dTheta4A(lm_max) )
      allocate( hdif_B(0:l_max),hdif_V(0:l_max),hdif_S(0:l_max))
      allocate( hdif_Xi(0:l_max) )
      bytes_allocated = bytes_allocated+(10*lm_max+4*(l_max+1))*SIZEOF_DEF_REAL

   end subroutine initialize_horizontal_data
!------------------------------------------------------------------------------
   subroutine finalize_horizontal_data
      !
      ! Memory deallocation of horizontal functions
      !

      deallocate( cosn_theta_E2, sinTheta, cosTheta, theta_ord, n_theta_cal2ord )
      deallocate( sinTheta_E2, O_sin_theta, O_sin_theta_E2, phi )
      deallocate( gauss, dPl0Eq, n_theta_ord2cal )
      deallocate( dPhi, dLh, dTheta1S, dTheta1A )
      deallocate( dTheta2S, dTheta2A, dTheta3S, dTheta3A, dTheta4S, dTheta4A )
      deallocate( hdif_B, hdif_V, hdif_S, hdif_Xi )

      if ( .not. l_axi ) call finalize_fft()

   end subroutine finalize_horizontal_data
!------------------------------------------------------------------------------
   subroutine horizontal()
      !
      !  Calculates functions of :math:`\theta` and :math:`\phi`, for example
      !  the Legendre functions,  and functions of degree :math:`\ell` and
      !  order :math:`m` of the legendres.
      !

      !-- Local variables:
      integer :: norm,n_theta,n_phi
      integer :: l,m,lm
      real(cp) :: clm(0:l_max+1,0:l_max+1), Pl0Eq(l_max+1), tmp_gauss(n_theta_max)
      real(cp) :: colat, fac

      norm=2 ! norm chosen so that a surface integral over
             ! any ylm**2 is 1.

      !-- Calculate grid points and weights for the
      !   Gauss-Legendre integration of the plms:
      call gauleg(-one,one,theta_ord,tmp_gauss,n_theta_max)

      !-- Legendre polynomials and cos(theta) derivative:
      do n_theta=1,n_theta_max/2  ! Loop over colat in NHS

         colat=theta_ord(n_theta)

         ! Get dP for all degrees and order m=0 at the equator only
         ! Usefull to estimate the flow velocity at the equator
         call plm_theta(half*pi,l_max,0,0,minc,Pl0Eq,dPl0Eq,l_max+1,norm)

         if ( l_scramble_theta ) then
            O_sin_theta(2*n_theta-1)   =one/sin(colat)
            O_sin_theta(2*n_theta  )   =one/sin(colat)
            O_sin_theta_E2(2*n_theta-1)=one/(sin(colat)*sin(colat))
            O_sin_theta_E2(2*n_theta  )=one/(sin(colat)*sin(colat))
            sinTheta(2*n_theta-1)      =sin(colat)
            sinTheta(2*n_theta  )      =sin(colat)
            sinTheta_E2(2*n_theta-1)   =sin(colat)**2
            sinTheta_E2(2*n_theta  )   =sin(colat)**2
            cosTheta(2*n_theta-1)      =cos(colat)
            cosTheta(2*n_theta  )      =-cos(colat)
            cosn_theta_E2(2*n_theta-1) =cos(colat)/sin(colat)/sin(colat)
            cosn_theta_E2(2*n_theta)   =-cos(colat)/sin(colat)/sin(colat)
         else
            O_sin_theta(n_theta)                 =one/sin(colat)
            O_sin_theta(n_theta_max-n_theta+1)   =one/sin(colat)
            O_sin_theta_E2(n_theta)              =one/(sin(colat)*sin(colat))
            O_sin_theta_E2(n_theta_max-n_theta+1)=one/(sin(colat)*sin(colat))
            sinTheta(n_theta)                    =sin(colat)
            sinTheta(n_theta_max-n_theta+1)      =sin(colat)
            sinTheta_E2(n_theta)                 =sin(colat)**2
            sinTheta_E2(n_theta_max-n_theta+1)   =sin(colat)**2
            cosTheta(n_theta)                    =cos(colat)
            cosTheta(n_theta_max-n_theta+1)      =-cos(colat)
            cosn_theta_E2(n_theta)               =cos(colat)/sin(colat)/sin(colat)
            cosn_theta_E2(n_theta_max-n_theta+1) =-cos(colat)/sin(colat)/sin(colat)
         end if
      end do

      !-- Resort thetas in the alternating north/south order they
      !   are used for the calculations:
      if ( l_scramble_theta ) then
         do n_theta=1,n_theta_max/2
            n_theta_ord2cal(n_theta)              =2*n_theta-1
            n_theta_ord2cal(n_theta_max-n_theta+1)=2*n_theta
            n_theta_cal2ord(2*n_theta-1)=n_theta
            n_theta_cal2ord(2*n_theta)  =n_theta_max-n_theta+1
            gauss(2*n_theta-1)          =tmp_gauss(n_theta)
            gauss(2*n_theta)            =tmp_gauss(n_theta_max-n_theta+1)
         end do
      else
         do n_theta=1,n_theta_max
            n_theta_ord2cal(n_theta)=n_theta
            n_theta_cal2ord(n_theta)=n_theta
            gauss(n_theta)          =tmp_gauss(n_theta)
         end do
      end if

      !----- Same for longitude output grid:
      fac=two*pi/real(n_phi_max*minc,cp)
      do n_phi=1,n_phi_max
         phi(n_phi)=fac*real(n_phi-1,cp)
      end do

      !-- Initialize fast fourier transform for phis:
      if ( .not. l_axi ) call init_fft(n_phi_max)

      !-- Build arrays depending on degree l and order m
      !   and hyperdiffusion factors:
      do m=m_min,m_max,minc  ! Build auxiliary array clm
         do l=m,l_max+1
            clm(l,m)=sqrt( real((l+m)*(l-m),cp) / real((2*l-1)*(2*l+1),cp) )
         end do
      end do

      !---- Operators for derivatives:
      do lm=1,lm_max
         l=lm2l(lm)
         m=lm2m(lm)

         !-- Phi derivative:
         dPhi(lm)=cmplx(0.0_cp,real(m,cp),cp)
         !-- Negative horizontal Laplacian *r^2
         dLh(lm)     =real(l*(l+1),cp)                 ! = qll1
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 )
         dTheta1S(lm)=real(l+1,cp)        *clm(l,m)    ! = qcl1
         dTheta1A(lm)=real(l,cp)          *clm(l+1,m)  ! = qcl
         !-- Operator ( sin(thetaR) * d/d theta )
         dTheta2S(lm)=real(l-1,cp)        *clm(l,m)    ! = qclm1
         dTheta2A(lm)=real(l+2,cp)        *clm(l+1,m)  ! = qcl2
         !-- Operator ( sin(theta) * d/d theta + cos(theta) dLh )
         dTheta3S(lm)=real((l-1)*(l+1),cp)*clm(l,m)    ! = q0l1lm1(lm)
         dTheta3A(lm)=real(l*(l+2),cp)    *clm(l+1,m)  ! = q0cll2(lm)
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 ) * dLh
         dTheta4S(lm)=dTheta1S(lm)*real((l-1)*l,cp)
         dTheta4A(lm)=dTheta1A(lm)*real((l+1)*(l+2),cp)
      end do ! lm

      !-- Hyperdiffusion
      do l=0,l_max

         !-- Hyperdiffusion
         hdif_B(l) =one
         hdif_V(l) =one
         hdif_S(l) =one
         hdif_Xi(l)=one
         if ( ldifexp > 0 ) then

            if ( ldif >= 0 .and. l > ldif ) then

               !-- Kuang and Bloxham type:
               !                 hdif_B(l) =one+difeta *real(l+1-ldif,cp)**ldifexp
               !                 hdif_V(l) =one+ difnu *real(l+1-ldif,cp)**ldifexp
               !                 hdif_S(l) =one+difkap *real(l+1-ldif,cp)**ldifexp
               !                 hdif_Xi(l)=one+difchem*real(l+1-ldif,cp)**ldifexp

               !-- Old type:
               hdif_B(l) = one + difeta  * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp
               hdif_V(l) = one + difnu   * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp
               hdif_S(l) = one + difkap  * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp
               hdif_Xi(l)= one + difchem * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp

            else if ( ldif < 0 ) then

               !-- Grote and Busse type:
               hdif_B(l) =(one+difeta *real(l,cp)**ldifexp ) / &
               &          (one+difeta *real(-ldif,cp)**ldifexp )
               hdif_V(l) =(one+difnu  *real(l,cp)**ldifexp ) / &
               &          (one+difnu  *real(-ldif,cp)**ldifexp )
               hdif_S(l) =(one+difkap *real(l,cp)**ldifexp ) / &
               &          (one+difkap *real(-ldif,cp)**ldifexp )
               hdif_Xi(l)=(one+difchem*real(l,cp)**ldifexp ) / &
               &          (one+difchem*real(-ldif,cp)**ldifexp )

            end if
         else if ( ldifexp <-1) then !hyperdiffusion of the large scale velocity for m=0 modes
            if (ldif >=0 .and. l <ldif) then
               hdif_V(l) = one + difnu*(real(l,cp))**(ldifexp+1) !
            end if
         end if
      end do

   end subroutine horizontal
!------------------------------------------------------------------------------
   subroutine gauleg(sinThMin,sinThMax,theta_ord,gauss,n_th_max)
      !
      ! Subroutine is based on a NR code.
      ! Calculates N zeros of legendre polynomial P(l=N) in
      ! the interval [sinThMin,sinThMax].
      ! Zeros are returned in radiants theta_ord(i)
      ! The respective weights for Gauss-integration are given in gauss(i).
      !

      !-- Input variables:
      real(cp), intent(in) :: sinThMin ! lower bound in radiants
      real(cp), intent(in) :: sinThMax ! upper bound in radiants
      integer,  intent(in) :: n_th_max ! desired maximum degree

      !-- Output variables:
      real(cp), intent(out) :: theta_ord(n_th_max) ! zeros cos(theta)
      real(cp), intent(out) :: gauss(n_th_max)     ! associated Gauss-Legendre weights

      !-- Local variables:
      integer :: m,i,j
      real(cp) :: sinThMean,sinThDiff,p1,p2,p3,pp,z,z1
      real(cp), parameter :: eps = 10.0_cp*epsilon(one)

      ! use symmetry
      m=(n_th_max+1)/2

      !-- Map on symmetric interval:
      sinThMean=half*(sinThMax+sinThMin)
      sinThDiff=half*(sinThMax-sinThMin)

      do i=1,m
         !----- Initial guess for zeros:
         z  = cos( pi*( (real(i,cp)-0.25_cp)/(real(n_th_max,cp)+half)) )
         z1 = z+10.0_cp*eps

         do while( abs(z-z1) > eps)
            !----- Use recurrence to calulate P(l=n_th_max,z=cos(theta))
            p1=one
            p2=0.0_cp
            ! do loop over degree !
            do j=1,n_th_max
               p3=p2
               p2=p1
               p1=( real(2*j-1,cp)*z*p2-real(j-1,cp)*p3 )/real(j,cp)
            end do

            !----- Newton method to refine zero: pp is derivative !
            pp=real(n_th_max,cp)*(z*p1-p2)/(z*z-one)
            z1=z
            z=z1-p1/pp
         end do

         !----- Another zero found
         theta_ord(i)           =acos(sinThMean+sinThDiff*z)
         theta_ord(n_th_max+1-i)=acos(sinThMean-sinThDiff*z)
         gauss(i)               =two*sinThDiff/((one-z*z)*pp*pp)
         gauss(n_th_max+1-i)    =gauss(i)

      end do

   end subroutine gauleg
!------------------------------------------------------------------------------
end module horizontal_data
