module horizontal_data
   !
   !  Module containing functions depending on longitude 
   !  and latitude plus help arrays depending on degree and order
   !

   use truncation, only: l_max, lmP_max, n_theta_max, n_phi_max, &
                         lm_max, n_m_max, minc, m_max
   use radial_functions, only: r_cmb
   use physical_parameters, only: ek
   use num_param, only: difeta, difnu, difkap, ldif, ldifexp
   use blocking, only: lmP2l, lmP2lm, lm2l, lm2m
   use logic, only: l_non_rot
   use plms_theta, only: plm_theta, plm_thetaAS
   use fft
   use constants, only: pi, zero, one, two, half
   use precision_mod
   use mem_alloc, only: bytes_allocated
 
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
   real(cp), public, allocatable :: sinTheta(:)
   real(cp), public, allocatable :: cosTheta(:)
 
   !-- Phi (longitude)
   real(cp), public, allocatable :: phi(:)
 
   !-- Legendres:
   real(cp), public, allocatable :: Plm(:,:)
   real(cp), public, allocatable :: wPlm(:,:)
   real(cp), public, allocatable :: dPlm(:,:)
   real(cp), public, allocatable :: gauss(:)
   real(cp), public, allocatable :: dPl0Eq(:)
 
   !-- Arrays depending on l and m:
   complex(cp), public, allocatable :: dPhi(:)
   complex(cp), public, allocatable :: dPhi0(:)
   complex(cp), public, allocatable :: dPhi02(:)
   real(cp), public, allocatable :: dLh(:)
   real(cp), public, allocatable :: dTheta1S(:),dTheta1A(:)
   real(cp), public, allocatable :: dTheta2S(:),dTheta2A(:)
   real(cp), public, allocatable :: dTheta3S(:),dTheta3A(:)
   real(cp), public, allocatable :: dTheta4S(:),dTheta4A(:)
   real(cp), public, allocatable :: D_m(:),D_l(:),D_lP1(:)
   real(cp), public, allocatable :: D_mc2m(:)
   real(cp), public, allocatable :: hdif_B(:),hdif_V(:),hdif_S(:)
 
   !-- Limiting l for a given m, used in legtf
   integer, public, allocatable :: lStart(:),lStop(:)
   integer, public, allocatable :: lStartP(:),lStopP(:)
   logical, public, allocatable :: lmOdd(:),lmOddP(:)

   public :: initialize_horizontal_data, horizontal

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
      allocate( cosTheta(n_theta_max) )
      bytes_allocated = bytes_allocated+n_theta_max*SIZEOF_INTEGER+&
                        8*n_theta_max*SIZEOF_DEF_REAL

      !-- Phi (longitude)
      allocate( phi(n_phi_max) )
      bytes_allocated = bytes_allocated+n_phi_max*SIZEOF_DEF_REAL

      !-- Legendres:
      allocate( Plm(lm_max,n_theta_max/2) )
      allocate( wPlm(lmP_max,n_theta_max/2) )
      allocate( dPlm(lm_max,n_theta_max/2) )
      allocate( gauss(n_theta_max) )
      allocate( dPl0Eq(l_max+1) )
      bytes_allocated = bytes_allocated+(lm_max*n_theta_max+ &
                        lmP_max*n_theta_max/2+n_theta_max+l_max+1)*SIZEOF_DEF_REAL

      !-- Arrays depending on l and m:
      allocate( dPhi(lm_max) )
      allocate( dPhi0(lm_max) )
      allocate( dPhi02(lm_max) )
      allocate( dLh(lm_max) )
      allocate( dTheta1S(lm_max),dTheta1A(lm_max) )
      allocate( dTheta2S(lm_max),dTheta2A(lm_max) )
      allocate( dTheta3S(lm_max),dTheta3A(lm_max) )
      allocate( dTheta4S(lm_max),dTheta4A(lm_max) )
      allocate( D_m(lm_max),D_l(lm_max),D_lP1(lm_max) )
      allocate( D_mc2m(n_m_max) )
      allocate( hdif_B(lm_max),hdif_V(lm_max),hdif_S(lm_max) )
      bytes_allocated = bytes_allocated+(18*lm_max+n_m_max)*SIZEOF_DEF_REAL

      !-- Limiting l for a given m, used in legtf
      allocate( lStart(n_m_max),lStop(n_m_max) )
      allocate( lStartP(n_m_max),lStopP(n_m_max) )
      allocate( lmOdd(n_m_max),lmOddP(n_m_max) )
      bytes_allocated = bytes_allocated+6*n_m_max*SIZEOF_INTEGER

   end subroutine initialize_horizontal_data
!------------------------------------------------------------------------------
   subroutine horizontal
      !
      !  Calculates functions of theta and phi, for exmample the
      !  Legendre functions, and functions of degree l and order m
      !  of the legendres.
      !

      !-- Local variables:
      integer :: norm,n_theta,n_phi
      integer :: l,m,lm,lmP,mc
      real(cp) :: ampnu!,q0
      real(cp) :: clm(0:l_max+1,0:l_max+1)
      real(cp) :: plma(lmP_max)
      real(cp) :: dtheta_plma(lmP_max)
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

          !----- plmtheta calculates plms and their derivatives
          !      up to degree and order l_max+1 and m_max at
          !      the points cos(theta_ord(n_theta)):
          call plm_theta(colat,l_max+1,m_max,minc, &
                         plma,dtheta_plma,lmP_max,norm)
          do lmP=1,lmP_max
              l=lmP2l(lmP)
              if ( l <= l_max ) then
                  lm=lmP2lm(lmP)
                  Plm(lm,n_theta) =plma(lmP)
                  dPlm(lm,n_theta)=dtheta_plma(lmP)
              end if
              wPlm(lmP,n_theta)=two*pi*gauss(n_theta)*plma(lmP)
          end do

          ! Get dP for all degrees and order m=0 at the equator only
          ! Usefull to estimate the flow velocity at the equator
          call plm_thetaAS(half*pi,l_max,Pl0Eq,dPl0Eq,l_max+1,norm)

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
      call init_fft(n_phi_max)

      !-- Build arrays depending on degree l and order m
      !   and hyperdiffusion factors:
      do m=0,m_max,minc  ! Build auxiliary array clm
          do l=m,l_max+1
              clm(l,m)=sqrt( real((l+m)*(l-m),cp) / real((2*l-1)*(2*l+1),cp) )
          end do
      end do

      do lm=1,lm_max
          l=lm2l(lm)
          m=lm2m(lm)

          !-- Help arrays:
          D_l(lm)  =real(l,cp)
          D_lP1(lm)=real(l+1,cp)
          D_m(lm)  =real(m,cp)

          !---- Operators for derivatives:

          !-- Phi derivate:
          dPhi(lm)=cmplx(0.0_cp,real(m,cp),cp)
          if ( l < l_max ) then
              dPhi0(lm)    =cmplx(0.0_cp,real(m,cp),cp)
              dPhi02(lm)   =dPhi0(lm)*dPhi0(lm)
          else
              dPhi0(lm)    =zero
              dPhi02(lm)   =zero
          end if
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

          !-- Hyperdiffusion
          hdif_B(lm)=one
          hdif_V(lm)=one
          hdif_S(lm)=one
          if ( ldifexp > 0 ) then

              if ( ldif >= 0 .and. l > ldif ) then

              !-- Kuang and Bloxham type:
              !                 hdif_B(lm)=
              !     *                   one+difeta*real(l+1-ldif,cp)**ldifexp
              !                 hdif_V(lm)=
              !     *                   one+ difnu*real(l+1-ldif,cp)**ldifexp
              !                 hdif_S(lm)=
              !     &                   one+difkap*real(l+1-ldif,cp)**ldifexp

              !-- Old type:
                  hdif_B(lm)= one + difeta * ( real(l+1-ldif,cp) / &
                                                real(l_max+1-ldif,cp) )**ldifexp
                  hdif_V(lm)= one + difnu * ( real(l+1-ldif,cp) / &
                                               real(l_max+1-ldif,cp) )**ldifexp
                  hdif_S(lm)= one + difkap * ( real(l+1-ldif,cp) / &
                                                real(l_max+1-ldif,cp) )**ldifexp

              else if ( ldif < 0 ) then

              !-- Grote and Busse type:
                  hdif_B(lm)= (one+difeta*real(l,cp)**ldifexp ) / &
                              (one+difeta*real(-ldif,cp)**ldifexp )
                  hdif_V(lm)= (one+difnu*real(l,cp)**ldifexp ) / &
                              (one+difnu*real(-ldif,cp)**ldifexp )
                  hdif_S(lm)= (one+difkap*real(l,cp)**ldifexp ) / &
                              (one+difkap*real(-ldif,cp)**ldifexp )
                               
              end if

          else

              if ( l == l_max .and. .not. l_non_rot ) then
              !  Chose ampnu so that the viscous force is at least as
              !  strong as the viscous force for l=l_max:
              !  We can turn this argument around and state that
              !  for example for Ek=1e-4 l_max should be 221.
                  ampnu=(r_cmb**2/real(l_max*(l_max+1),cp))*(two/ek)
                  ampnu=max(one,ampnu)
                  hdif_V(lm)=ampnu*hdif_V(lm)
              end if

          end if

      end do ! lm


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

#if 0
      write(*,"(A,I6,A,I10,A,I10,A)")                                           &
           & "cache info of first element, all have dimension 1:lm_max=",lm_max,&
           & " = ",lm_max*16,"B (C), ",lm_max*8,"B (R)"
      call print_cache_info_dcmplx("C: dPhi0"//C_NULL_CHAR,dPhi0(1))
      call print_cache_info_dcmplx("C: dPhi"//C_NULL_CHAR,dPhi(1))
      call print_cache_info_dreal("R: dTheta1A"//C_NULL_CHAR,dTheta1A(1))
      call print_cache_info_dreal("R: dTheta1S"//C_NULL_CHAR,dTheta1S(1))
      call print_cache_info_dreal("R: dTheta2A"//C_NULL_CHAR,dTheta2A(1))
      call print_cache_info_dreal("R: dTheta2S"//C_NULL_CHAR,dTheta2S(1))
      call print_cache_info_dreal("R: dTheta3A"//C_NULL_CHAR,dTheta3A(1))
      call print_cache_info_dreal("R: dTheta3S"//C_NULL_CHAR,dTheta3S(1))
      call print_cache_info_dreal("R: dTheta4A"//C_NULL_CHAR,dTheta4A(1))
      call print_cache_info_dreal("R: dTheta4S"//C_NULL_CHAR,dTheta4S(1))
#endif

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
         theta_ord(i)              =acos(sinThMean+sinThDiff*z)
         theta_ord(n_th_max+1-i)=acos(sinThMean-sinThDiff*z)
         gauss(i)                  =two*sinThDiff/((one-z*z)*pp*pp)
         gauss(n_th_max+1-i)    =gauss(i)
    
      end do
     
   end subroutine gauleg
!------------------------------------------------------------------------------
end module horizontal_data
