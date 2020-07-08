module horizontal_data
   !
   !  Module containing functions depending on longitude
   !  and latitude plus help arrays depending on degree and order
   !

   use truncation
   use radial_functions, only: r_cmb
   use physical_parameters, only: ek
   use num_param, only: difeta, difnu, difkap, ldif, ldifexp, difchem
   use logic, only: l_non_rot, l_RMS
   use plms_theta, only: plm_theta
   use fft
   use constants, only: pi, zero, one, two, half
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use LMmapping, only: map_dist_st, map_glbl_st
   use communications, only: slice_Flm_real

   implicit none

   private

   !-- Arrays depending on theta (colatitude):
   integer, public, allocatable :: n_theta_cal2ord(:)
   real(cp), public, allocatable :: theta(:)
   real(cp), public, allocatable :: theta_ord(:)
   real(cp), public, allocatable :: sn2(:)
   real(cp), public, allocatable :: osn2(:)
   real(cp), public, allocatable :: cosn2(:)
   real(cp), public, allocatable :: osn1(:)
   real(cp), public, allocatable :: O_sin_theta(:)
   real(cp), public, allocatable :: O_sin_theta_E2(:)
   real(cp), public, allocatable :: cosn_theta_E2(:)
   real(cp), public, allocatable :: sinTheta(:)
   real(cp), public, allocatable :: cosTheta(:)

   !-- Phi (longitude)
   real(cp), public, allocatable :: phi(:)

   !-- Legendres:
   real(cp), public, allocatable :: Plm(:,:)
   real(cp), public, allocatable :: Plm_loc(:,:)
   real(cp), public, allocatable :: dPlm_loc(:,:)
   real(cp), public, allocatable :: wPlm(:,:)
   real(cp), public, allocatable :: wdPlm(:,:)
   real(cp), public, allocatable :: dPlm(:,:)
   real(cp), public, allocatable :: gauss(:)
   real(cp), public, allocatable :: dPl0Eq(:)

   !-- Arrays depending on l and m:
   complex(cp), public, allocatable :: dPhi(:)
   real(cp), public, allocatable :: dLh(:)
   real(cp), public, allocatable :: dTheta1S(:),dTheta1A(:)
   real(cp), public, allocatable :: D_mc2m(:)
   real(cp), public, allocatable :: hdif_B(:),hdif_V(:),hdif_S(:),hdif_Xi(:)

   complex(cp), public, allocatable :: dPhi_loc(:)
   real(cp), public, allocatable :: dLh_loc(:)
   real(cp), public, allocatable :: dTheta1S_loc(:),dTheta1A_loc(:)
   real(cp), public, allocatable :: dTheta2S_loc(:),dTheta2A_loc(:)
   real(cp), public, allocatable :: dTheta3S_loc(:),dTheta3A_loc(:)
   real(cp), public, allocatable :: dTheta4S_loc(:),dTheta4A_loc(:)
   
   !-- Limiting l for a given m, used in legtf
   integer, public, allocatable :: lStart(:),lStop(:)
   integer, public, allocatable :: lStartP(:),lStopP(:)
   logical, public, allocatable :: lmOdd(:),lmOddP(:)

   public :: initialize_horizontal_data, horizontal, finalize_horizontal_data

contains

   subroutine initialize_horizontal_data

      allocate( n_theta_cal2ord(n_theta_max) )
      allocate( theta(n_theta_max) )
      allocate( theta_ord(n_theta_max) )
      allocate( sn2(n_theta_max/2) )
      allocate( osn2(n_theta_max/2) )
      allocate( cosn2(n_theta_max/2) )
      allocate( osn1(n_theta_max/2) )
      allocate( O_sin_theta(n_theta_max) )
      allocate( O_sin_theta_E2(n_theta_max) )
      allocate( sinTheta(n_theta_max) )
      allocate( cosn_theta_E2(n_theta_max) )
      allocate( cosTheta(n_theta_max) )
      bytes_allocated = bytes_allocated+n_theta_max*SIZEOF_INTEGER+&
      &                 9*n_theta_max*SIZEOF_DEF_REAL

      !-- Phi (longitude)
      allocate( phi(n_phi_max) )
      bytes_allocated = bytes_allocated+n_phi_max*SIZEOF_DEF_REAL

      !-- Legendres:
      allocate( gauss(n_theta_max) )
      allocate( dPl0Eq(l_max+1) )
      bytes_allocated = bytes_allocated+(n_theta_max+l_max+1)*SIZEOF_DEF_REAL

#ifndef WITH_SHTNS
      allocate( Plm(lm_max,n_theta_max/2) )
      allocate( Plm_loc(n_lm_loc,n_theta_max/2) )
      allocate( dPlm_loc(n_lm_loc,n_theta_max/2) )
      allocate( wPlm(lmP_max,n_theta_max/2) )
      allocate( dPlm(lm_max,n_theta_max/2) )
      bytes_allocated = bytes_allocated+(lm_max*n_theta_max+ &
      &                 lmP_max*n_theta_max/2)*SIZEOF_DEF_REAL

      if ( l_RMS ) then
         allocate( wdPlm(lmP_max,n_theta_max/2) )
         bytes_allocated = bytes_allocated*lmP_max*n_theta_max/2*SIZEOF_DEF_REAL
      end if
#endif

      !-- Arrays depending on l and m:
      allocate( dPhi(lm_max), dLh(lm_max) )
      allocate( dTheta1S(lm_max),dTheta1A(lm_max) )
      allocate( D_mc2m(n_m_max) )
      allocate( hdif_B(0:l_max),hdif_V(0:l_max),hdif_S(0:l_max) )
      allocate( hdif_Xi(0:l_max) )
      bytes_allocated = bytes_allocated+(4*lm_max+n_m_max+4*(l_max+1))* &
      &                 SIZEOF_DEF_REAL
      
      !-- Distributed arrays depending on l and m:)
      !>@ TODO : decide which one we keep 
      allocate( dPhi_loc(n_lm_loc), dLh_loc(n_lm_loc) )
      allocate( dTheta1S_loc(n_lm_loc),dTheta1A_loc(n_lm_loc) )
      allocate( dTheta2S_loc(n_lm_loc),dTheta2A_loc(n_lm_loc) )
      allocate( dTheta3S_loc(n_lm_loc),dTheta3A_loc(n_lm_loc) )
      allocate( dTheta4S_loc(n_lm_loc),dTheta4A_loc(n_lm_loc) )
      bytes_allocated = bytes_allocated+10*n_lm_loc*SIZEOF_DEF_REAL

      !-- Limiting l for a given m, used in legtf
      allocate( lStart(n_m_max),lStop(n_m_max) )
      allocate( lStartP(n_m_max),lStopP(n_m_max) )
      allocate( lmOdd(n_m_max),lmOddP(n_m_max) )
      bytes_allocated = bytes_allocated+6*n_m_max*SIZEOF_INTEGER

   end subroutine initialize_horizontal_data
!------------------------------------------------------------------------------
   subroutine finalize_horizontal_data

      deallocate( cosn_theta_E2, sinTheta, cosTheta, theta, theta_ord, n_theta_cal2ord )
      deallocate( sn2, osn2, cosn2, osn1, O_sin_theta, O_sin_theta_E2, phi )
      deallocate( gauss, dPl0Eq )
#ifndef WITH_SHTNS
      deallocate( Plm, Plm_loc, dPlm_loc, wPlm, dPlm )
      if ( l_RMS ) deallocate( wdPlm )
#endif
      deallocate( dPhi, dLh, dTheta1S, dTheta1A )
      deallocate( D_mc2m, hdif_B, hdif_V, hdif_S, hdif_Xi )
      deallocate( lStart, lStop, lStartP, lStopP, lmOdd, lmOddP )

      ! Deallocate the distributed fields 
      deallocate( dLh_loc, dPhi_loc, dTheta1S_loc, dTheta1A_loc )
      deallocate( dTheta2S_loc, dTheta2A_loc, dTheta3S_loc, dTheta3A_loc )
      deallocate( dTheta4S_loc, dTheta4A_loc )
      
      if ( .not. l_axi ) call finalize_fft()

   end subroutine finalize_horizontal_data
!------------------------------------------------------------------------------
   subroutine horizontal
      !
      !  Calculates functions of theta and phi, for exmample the
      !  Legendre functions, and functions of degree l and order m
      !  of the legendres.
      !

      !-- Local variables:
      integer :: norm,n_theta,n_phi
      integer :: l,m,lm,mc
      real(cp) :: ampnu!,q0
      real(cp) :: clm(0:l_max+1,0:l_max+1)
#ifndef WITH_SHTNS
      integer :: lmP
      real(cp) :: plma(lmP_max)
      real(cp) :: dtheta_plma(lmP_max)
#endif
      real(cp) :: colat
      real(cp) :: fac
      real(cp) :: Pl0Eq(l_max+1)

      norm=2 ! norm chosen so that a surface integral over
             ! any ylm**2 is 1.

      !-- Calculate grid points and weights for the
      !   Gauss-Legendre integration of the plms:
      call gauleg(-one,one,theta_ord,gauss,n_theta_max)

      !-- Legendre polynomials and cos(theta) derivative:
      !   Note: the following functions are only stored for the northern hemisphere
      !         southern hemisphere values differ from the northern hemisphere
      !         values by a sign that depends on the symmetry of the function.
      !         The only asymmetric function (sign=-1) stored here is cosn2 !
      do n_theta=1,n_theta_max/2  ! Loop over colat in NHS

         colat=theta_ord(n_theta)

#ifndef WITH_SHTNS
         !----- plmtheta calculates plms and their derivatives
         !      up to degree and order l_max+1 and m_max at
         !      the points cos(theta_ord(n_theta)):
         call plm_theta(colat,l_max+1,m_max,minc, &
              &         plma,dtheta_plma,lmP_max,norm)
         do lmP=1,lmP_max
            l=map_glbl_st%lmP2l(lmP)
            if ( l <= l_max ) then
               lm=map_glbl_st%lmP2lm(lmP)
               Plm(lm,n_theta) =plma(lmP)
               dPlm(lm,n_theta)=dtheta_plma(lmP)
            end if
            wPlm(lmP,n_theta) =two*pi*gauss(n_theta)*plma(lmP)
            if ( l_RMS ) wdPlm(lmP,n_theta)=two*pi*gauss(n_theta)*dtheta_plma(lmP)
         end do

         call slice_Flm_real(Plm(:,n_theta), Plm_loc(:,n_theta))
         call slice_Flm_real(dPlm(:,n_theta), dPlm_loc(:,n_theta))
#endif

         ! Get dP for all degrees and order m=0 at the equator only
         ! Usefull to estimate the flow velocity at the equator
         call plm_theta(half*pi,l_max,0,minc,Pl0Eq,dPl0Eq,l_max+1,norm)

         !-- More functions stored to obscure the code:
         sn2(n_theta)               =sin(colat)**2
         osn1(n_theta)              =one/sin(colat)
         osn2(n_theta)              =osn1(n_theta)*osn1(n_theta)
         cosn2(n_theta)             =cos(colat)*osn2(n_theta)
         O_sin_theta(2*n_theta-1)   =one/sin(colat)
         O_sin_theta(2*n_theta  )   =one/sin(colat)
         O_sin_theta_E2(2*n_theta-1)=one/(sin(colat)*sin(colat))
         O_sin_theta_E2(2*n_theta  )=one/(sin(colat)*sin(colat))
         sinTheta(2*n_theta-1)      =sin(colat)
         sinTheta(2*n_theta  )      =sin(colat)
         cosTheta(2*n_theta-1)      =cos(colat)
         cosTheta(2*n_theta  )      =-cos(colat)
         cosn_theta_E2(2*n_theta-1) =cos(colat)/sin(colat)/sin(colat)
         cosn_theta_E2(2*n_theta)   =-cos(colat)/sin(colat)/sin(colat)
      end do


      !-- Resort thetas in the alternating north/south order they
      !   are used for the calculations:
      do n_theta=1,n_theta_max/2
         n_theta_cal2ord(2*n_theta-1)=n_theta
         n_theta_cal2ord(2*n_theta)  =n_theta_max-n_theta+1
         theta(2*n_theta-1)          =theta_ord(n_theta)
         theta(2*n_theta)            =theta_ord(n_theta_max-n_theta+1)
      end do


      !----- Same for longitude output grid:
      fac=two*pi/real(n_phi_max*minc,cp)
      do n_phi=1,n_phi_max
         phi(n_phi)=fac*real(n_phi-1,cp)
      end do


      !-- Initialize fast fourier transform for phis:
      if ( .not. l_axi ) call init_fft(n_phi_max)

      !-- Build arrays depending on degree l and order m
      !   and hyperdiffusion factors:
      do m=0,m_max,minc  ! Build auxiliary array clm
         do l=m,l_max+1
            clm(l,m)=sqrt( real((l+m)*(l-m),cp) / real((2*l-1)*(2*l+1),cp) )
         end do
      end do

      do lm=1,lm_max
         l=map_glbl_st%lm2l(lm)
         m=map_glbl_st%lm2m(lm)

         !---- Operators for derivatives:
         !-- Phi derivative:
         dPhi(lm)=cmplx(0.0_cp,real(m,cp),cp)
         !-- Negative horizontal Laplacian *r^2
         dLh(lm)     =real(l*(l+1),cp)                 ! = qll1
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 )
         dTheta1S(lm)=real(l+1,cp)        *clm(l,m)    ! = qcl1
         dTheta1A(lm)=real(l,cp)          *clm(l+1,m)  ! = qcl

      end do ! lm

      !>@ TODO: probably keep only the following once everything is finished
      do lm=1,n_lm_loc
         l=map_dist_st%lm2l(lm)
         m=map_dist_st%lm2m(lm)

         !---- Operators for derivatives:
         !-- Phi derivative:
         dPhi_loc(lm)=cmplx(0.0_cp,real(m,cp),cp)
         !-- Negative horizontal Laplacian *r^2
         dLh_loc(lm)     =real(l*(l+1),cp)                 ! = qll1
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 )
         dTheta1S_loc(lm)=real(l+1,cp)        *clm(l,m)    ! = qcl1
         dTheta1A_loc(lm)=real(l,cp)          *clm(l+1,m)  ! = qcl
         !-- Operator ( sin(thetaR) * d/d theta )
         dTheta2S_loc(lm)=real(l-1,cp)        *clm(l,m)    ! = qclm1
         dTheta2A_loc(lm)=real(l+2,cp)        *clm(l+1,m)  ! = qcl2
         !-- Operator ( sin(theta) * d/d theta + cos(theta) dLh )
         dTheta3S_loc(lm)=real((l-1)*(l+1),cp)*clm(l,m)    ! = q0l1lm1(lm)
         dTheta3A_loc(lm)=real(l*(l+2),cp)    *clm(l+1,m)  ! = q0cll2(lm)
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 ) * dLh
         dTheta4S_loc(lm)=dTheta1S_loc(lm)*real((l-1)*l,cp)
         dTheta4A_loc(lm)=dTheta1A_loc(lm)*real((l+1)*(l+2),cp)
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
            !                 hdif_B(l)=
            !     *                   one+difeta*real(l+1-ldif,cp)**ldifexp
            !                 hdif_V(l)=
            !     *                   one+ difnu*real(l+1-ldif,cp)**ldifexp
            !                 hdif_S(l)=
            !     &                   one+difkap*real(l+1-ldif,cp)**ldifexp

            !-- Old type:
               hdif_B(l)= one + difeta * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp
               hdif_V(l)= one + difnu * ( real(l+1-ldif,cp) / &
               &                            real(l_max+1-ldif,cp) )**ldifexp
               hdif_S(l)= one + difkap * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp
               hdif_Xi(l)= one + difchem * ( real(l+1-ldif,cp) / &
               &                             real(l_max+1-ldif,cp) )**ldifexp

             else if ( ldif < 0 ) then

             !-- Grote and Busse type:
                hdif_B(l)= (one+difeta*real(l,cp)**ldifexp ) / &
                &           (one+difeta*real(-ldif,cp)**ldifexp )
                hdif_V(l)= (one+difnu*real(l,cp)**ldifexp ) / &
                &           (one+difnu*real(-ldif,cp)**ldifexp )
                hdif_S(l)= (one+difkap*real(l,cp)**ldifexp ) / &
                &           (one+difkap*real(-ldif,cp)**ldifexp )
                hdif_Xi(l)=(one+difchem*real(l,cp)**ldifexp ) / &
                &           (one+difchem*real(-ldif,cp)**ldifexp )

             end if

         else

            if ( l == l_max .and. .not. l_non_rot ) then
            !  Chose ampnu so that the viscous force is at least as
            !  strong as the viscous force for l=l_max:
            !  We can turn this argument around and state that
            !  for example for Ek=1e-4 l_max should be 221.
               ampnu=(r_cmb**2/real(l_max*(l_max+1),cp))*(two/ek)
               ampnu=max(one,ampnu)
               hdif_V(l)=ampnu*hdif_V(l)
            end if

         end if

      end do

      !-- Build auxiliary index arrays for Legendre transform:
      !   lStartP, lStopP give start and end positions in lmP-block.
      !   lStart, lStop give start and end positions in lm-block.
      lStartP(1)=1
      lStopP(1) =l_max+2
      lStart(1) =1
      lStop(1)  =l_max+1
      D_mc2m(1)=0
      if ( mod(l_max,2) == 0 ) then
         lmOdd(1) =.true.
         lmOddP(1)=.false.
      else
         lmOdd(1) =.false.
         lmOddP(1)=.true.
      end if
      do mc=2,n_m_max
         m=(mc-1)*minc
         D_mc2m(mc) =real(m,cp)
         lStartP(mc)=lStopP(mc-1)+1
         lStopP(mc) =lStartP(mc) +l_max-m+1
         lStart(mc) =lStop(mc-1) +1
         lStop(mc)  =lStart(mc)  +l_max-m
         if ( mod(lStop(mc)-lStart(mc),2) == 0 ) then
            lmOdd(mc) =.true.
            lmOddP(mc)=.false.
         else
            lmOdd(mc) =.false.
            lmOddP(mc)=.true.
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
