!$Id$
module horizontal_data
   !------------------------------------------------------------------
   !  Module containing functions depending on longitude 
   !  and latitude plus help arrays depending on degree and order
   !------------------------------------------------------------------

   use truncation, only: l_max, lmP_max, n_theta_max, n_phi_max, &
                         lm_max, n_m_max, minc, m_max
   use radial_functions, only: r_cmb
   use physical_parameters, only: ek
   use num_param, only: difeta, difnu, difkap, ldif, ldifexp
   use blocking, only: lmP2l, lmP2lm, lm2l, lm2m
   use logic, only: l_non_rot
   use plms_theta, only: plm_theta, plm_thetaAS
#if (FFTLIB==JW)
   use fft_JW, only: init_fft
#elif (FFTLIB==MKL)
   use fft_MKL, only: init_fft
#endif
   use const, only: pi
 
   implicit none

   private
 
   !-- Arrays depending on theta (colatitude):
   integer, public, allocatable :: n_theta_cal2ord(:)
   real(kind=8), public, allocatable :: theta(:)
   real(kind=8), public, allocatable :: theta_ord(:)
   real(kind=8), public, allocatable :: sn2(:)
   real(kind=8), public, allocatable :: osn2(:)
   real(kind=8), public, allocatable :: cosn2(:)
   real(kind=8), public, allocatable :: osn1(:)
   real(kind=8), public, allocatable :: O_sin_theta(:)
   real(kind=8), public, allocatable :: O_sin_theta_E2(:)
   real(kind=8), public, allocatable :: sinTheta(:)
   real(kind=8), public, allocatable :: cosTheta(:)
 
   !-- Phi (longitude)
   real(kind=8), public, allocatable :: phi(:)
 
   !-- Legendres:
   real(kind=8), public, allocatable :: Plm(:,:)
   real(kind=8), public, allocatable :: wPlm(:,:)
   real(kind=8), public, allocatable :: dPlm(:,:)
   real(kind=8), public, allocatable :: gauss(:)
   real(kind=8), public, allocatable :: dPl0Eq(:)
 
   !-- Arrays depending on l and m:
   complex(kind=8), public, allocatable :: dPhi(:)
   complex(kind=8), public, allocatable :: dPhi0(:)
   complex(kind=8), public, allocatable :: dPhi02(:)
   real(kind=8), public, allocatable :: dLh(:)
   real(kind=8), public, allocatable :: dTheta1S(:),dTheta1A(:)
   real(kind=8), public, allocatable :: dTheta2S(:),dTheta2A(:)
   real(kind=8), public, allocatable :: dTheta3S(:),dTheta3A(:)
   real(kind=8), public, allocatable :: dTheta4S(:),dTheta4A(:)
   real(kind=8), public, allocatable :: D_m(:),D_l(:),D_lP1(:)
   real(kind=8), public, allocatable :: D_mc2m(:)
   real(kind=8), public, allocatable :: hdif_B(:),hdif_V(:),hdif_S(:)
 
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

      !-- Phi (longitude)
      allocate( phi(n_phi_max) )

      !-- Legendres:
      allocate( Plm(lm_max,n_theta_max/2) )
      allocate( wPlm(lmP_max,n_theta_max/2) )
      allocate( dPlm(lm_max,n_theta_max/2) )
      allocate( gauss(n_theta_max) )
      allocate( dPl0Eq(l_max+1) )

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

      !-- Limiting l for a given m, used in legtf
      allocate( lStart(n_m_max),lStop(n_m_max) )
      allocate( lStartP(n_m_max),lStopP(n_m_max) )
      allocate( lmOdd(n_m_max),lmOddP(n_m_max) )

   end subroutine initialize_horizontal_data
!------------------------------------------------------------------------------
   subroutine horizontal
      !----------------------------------------------------------------
      !  Calculates functions of theta and phi, for exmample the
      !  Legendre functions, and functions of degree l and order m
      !  of the legendres.
      !----------------------------------------------------------------

      !-- Local variables:
      integer :: norm,n_theta,n_phi
      integer :: l,m,lm,lmP,mc
      real(kind=8) :: ampnu!,q0
      real(kind=8) :: clm(0:l_max+1,0:l_max+1)
      real(kind=8) :: plma(lmP_max)
      real(kind=8) :: dtheta_plma(lmP_max)
      real(kind=8) :: colat
      real(kind=8) :: fac
      real(kind=8) :: Pl0Eq(l_max+1)

      norm=2 ! norm chosen so that a surface integral over
             ! any ylm**2 is 1.

      !-- Calculate grid points and weights for the
      !   Gauss-Legendre integration of the plms:
      call gauleg(-1.d0,1.d0,theta_ord,gauss,n_theta_max)

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
              wPlm(lmP,n_theta)=2.d0*pi*gauss(n_theta)*plma(lmP)
          end do

          ! Get dP for all degrees and order m=0 at the equator only
          ! Usefull to estimate the flow velocity at the equator
          call plm_thetaAS(pi/2.d0,l_max,Pl0Eq,dPl0Eq,l_max+1,norm)

          !-- More functions stored to obscure the code:
          sn2(n_theta)               =dsin(colat)**2
          osn1(n_theta)              =1.D0/dsin(colat)
          osn2(n_theta)              =osn1(n_theta)*osn1(n_theta)
          cosn2(n_theta)             =dcos(colat)*osn2(n_theta)
          O_sin_theta(2*n_theta-1)   =1.D0/dsin(colat)
          O_sin_theta(2*n_theta  )   =1.D0/dsin(colat)
          O_sin_theta_E2(2*n_theta-1)=1.D0/(dsin(colat)*dsin(colat))
          O_sin_theta_E2(2*n_theta  )=1.D0/(dsin(colat)*dsin(colat))
          sinTheta(2*n_theta-1)      =dsin(colat)
          sinTheta(2*n_theta  )      =dsin(colat)
          cosTheta(2*n_theta-1)      =dcos(colat)
          cosTheta(2*n_theta  )      =-dcos(colat)
                    
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
      fac=2.d0*pi/dble(n_phi_max*minc)
      do n_phi=1,n_phi_max
          phi(n_phi)=fac*dble(n_phi-1)
      end do


      !-- Initialize fast fourier transform for phis:
      call init_fft(n_phi_max)

      !-- Build arrays depending on degree l and order m
      !   and hyperdiffusion factors:
      do m=0,m_max,minc  ! Build auxiliary array clm
          do l=m,l_max+1
              clm(l,m)=dsqrt( dble((l+m)*(l-m)) / dble((2*l-1)*(2*l+1)) )
          end do
      end do

      do lm=1,lm_max
          l=lm2l(lm)
          m=lm2m(lm)

          !-- Help arrays:
          D_l(lm)  =dble(l)
          D_lP1(lm)=dble(l+1)
          D_m(lm)  =dble(m)

          !---- Operators for derivatives:

          !-- Phi derivate:
          dPhi(lm)=cmplx(0.D0,dble(m),kind=kind(0d0))
          if ( l < l_max ) then
              dPhi0(lm)    =cmplx(0.D0,dble(m),kind=kind(0d0))
              dPhi02(lm)   =dPhi0(lm)*dPhi0(lm)
          else
              dPhi0(lm)    =cmplx(0.D0,0.D0,kind=kind(0d0))
              dPhi02(lm)   =cmplx(0.D0,0.D0,kind=kind(0d0))
          end if
          !-- Negative horizontal Laplacian *r^2
          dLh(lm)     =dble(l*(l+1))                 ! = qll1
          !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 )
          dTheta1S(lm)=dble(l+1)        *clm(l,m)    ! = qcl1
          dTheta1A(lm)=dble(l)          *clm(l+1,m)  ! = qcl
          !-- Operator ( sin(thetaR) * d/d theta )
          dTheta2S(lm)=dble(l-1)        *clm(l,m)    ! = qclm1
          dTheta2A(lm)=dble(l+2)        *clm(l+1,m)  ! = qcl2
          !-- Operator ( sin(theta) * d/d theta + cos(theta) dLh )
          dTheta3S(lm)=dble((l-1)*(l+1))*clm(l,m)    ! = q0l1lm1(lm)
          dTheta3A(lm)=dble(l*(l+2))    *clm(l+1,m)  ! = q0cll2(lm)
          !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 ) * dLh
          dTheta4S(lm)=dTheta1S(lm)*dble((l-1)*l)
          dTheta4A(lm)=dTheta1A(lm)*dble((l+1)*(l+2))

          !-- Hyperdiffusion
          hdif_B(lm)=1.D0
          hdif_V(lm)=1.D0
          hdif_S(lm)=1.D0
          if ( ldifexp > 0 ) then

              if ( ldif >= 0 .and. l > ldif ) then

              !-- Kuang and Bloxham type:
              !                 hdif_B(lm)=
              !     *                   1.D0+difeta*dble(l+1-ldif)**ldifexp
              !                 hdif_V(lm)=
              !     *                   1.D0+ difnu*dble(l+1-ldif)**ldifexp
              !                 hdif_S(lm)=
              !     &                   1.D0+difkap*dble(l+1-ldif)**ldifexp

              !-- Old type:
                  hdif_B(lm)= 1.D0 + difeta * ( dble(l+1-ldif) / &
                                                dble(l_max+1-ldif) )**ldifexp
                  hdif_V(lm)= 1.D0 + difnu * ( dble(l+1-ldif) / &
                                               dble(l_max+1-ldif) )**ldifexp
                  hdif_S(lm)= 1.D0 + difkap * ( dble(l+1-ldif) / &
                                                dble(l_max+1-ldif) )**ldifexp

              else if ( ldif < 0 ) then

              !-- Grote and Busse type:
                  hdif_B(lm)= (1.D0+difeta*dble(l)**ldifexp ) / &
                              (1.D0+difeta*dble(-ldif)**ldifexp )
                  hdif_V(lm)= (1.D0+difnu*dble(l)**ldifexp ) / &
                              (1.D0+difnu*dble(-ldif)**ldifexp )
                  hdif_S(lm)= (1.D0+difkap*dble(l)**ldifexp ) / &
                              (1.D0+difkap*dble(-ldif)**ldifexp )
                               
              end if

          else

              if ( l == l_max .and. .not. l_non_rot ) then
              !  Chose ampnu so that the viscous force is at least as
              !  strong as the viscous force for l=l_max:
              !  We can turn this argument around and state that
              !  for example for Ek=1e-4 l_max should be 221.
                  ampnu=(r_cmb**2/dble(l_max*(l_max+1)))*(2.D0/ek)
                  ampnu=Dmax1(1.D0,ampnu)
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
          D_mc2m(mc) =dble(m)
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
   subroutine gauleg(sinThMin,sinThMax,theta_ord,gauss,n_theta_max)
      !----------------------------------------------------------------------
      ! Subroutine is based on a NR code.
      ! Calculates N zeros of legendre polynomial P(l=N) in
      ! the interval [sinThMin,sinThMax].
      ! Zeros are returned in radiants theta_ord(i)
      ! The respective weights for Gauss-integration are given in gauss(i).
      !----------------------------------------------------------------------

      !-- Input variables:
      real(kind=8), intent(in) :: sinThMin,sinThMax ! lower/upper bound in radiants
      integer,      intent(in) :: n_theta_max  ! desired maximum degree
    
      !-- Output variables:
      real(kind=8), intent(out) :: theta_ord(n_theta_max) ! zeros cos(theta)
      real(kind=8), intent(out) :: gauss(n_theta_max) ! associated Gauss-Legendre weights
    
      !-- Local variables:
      integer                :: m,i,j
      real(kind=8)           :: sinThMean,sinThDiff,p1,p2,p3,pp,z,z1
      real(kind=8), parameter :: eps = 10.0D0*EPSILON(1.0D0)
    
      m=(n_theta_max+1)/2  ! use symmetry
    
      !-- Map on symmetric interval:
      sinThMean=0.5D0*(sinThMax+sinThMin)
      sinThDiff=0.5D0*(sinThMax-sinThMin)
    
      do i=1,m
         !----- Initial guess for zeros:
         z  = dcos( pi*( (dble(i)-0.25D0)/(dble(n_theta_max)+0.5D0)) )
         z1 = z+10*eps
     
         do while( dabs(z-z1) > eps)
            !----- Use recurrence to calulate P(l=n_theta_max,z=cos(theta))
            p1=1.D0
            p2=0.D0
            do j=1,n_theta_max   ! do loop over degree !
               p3=p2
               p2=p1
               p1=( dble(2*j-1)*z*p2-dbLe(j-1)*p3 )/dble(j)
            end do
      
            !----- Newton method to refine zero: pp is derivative !
            pp=dble(n_theta_max)*(z*p1-p2)/(z*z-1.D0)
            z1=z
            z=z1-p1/pp
         end do
     
         !----- Another zero found
         theta_ord(i)              =dacos(sinThMean+sinThDiff*z)
         theta_ord(n_theta_max+1-i)=dacos(sinThMean-sinThDiff*z)
         gauss(i)                  =2.D0*sinThDiff/((1.D0-z*z)*pp*pp)
         gauss(n_theta_max+1-i)    =gauss(i)
    
      end do
     
   end subroutine gauleg
!------------------------------------------------------------------------------
end module horizontal_data
