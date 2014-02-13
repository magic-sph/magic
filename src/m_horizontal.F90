!$Id$
!********************************************************************
!  Module containing functions depending on longitude 
!  and latitude plus help arrays depending on degree and order
!********************************************************************

MODULE horizontal_data
  USE truncation
  use omp_lib
  !USE cutils,ONLY: print_cache_info_dcmplx,print_cache_info_dreal,print_cache_info_integer

  IMPLICIT NONE

  !-- Arrays depending on theta (colatitude):
  INTEGER,ALLOCATABLE :: n_theta_cal2ord(:)
  REAL(kind=8),ALLOCATABLE :: theta(:)
  REAL(kind=8),ALLOCATABLE :: theta_ord(:)
  REAL(kind=8),ALLOCATABLE :: sn2(:)
  REAL(kind=8),ALLOCATABLE :: osn2(:)
  REAL(kind=8),ALLOCATABLE :: cosn2(:)
  REAL(kind=8),ALLOCATABLE :: osn1(:)
  REAL(kind=8),ALLOCATABLE :: O_sin_theta(:)
  REAL(kind=8),ALLOCATABLE :: O_sin_theta_E2(:)
  REAL(kind=8),ALLOCATABLE :: sinTheta(:)
  REAL(kind=8),ALLOCATABLE :: cosTheta(:)

  !-- Phi (longitude)
  REAL(kind=8),ALLOCATABLE :: phi(:)

  !-- Legendres:
  REAL(kind=8),ALLOCATABLE :: Plm(:,:)
  REAL(kind=8),ALLOCATABLE :: wPlm(:,:)
  REAL(kind=8),ALLOCATABLE :: dPlm(:,:)
  REAL(kind=8),ALLOCATABLE :: gauss(:)
  REAL(kind=8),ALLOCATABLE :: dPl0Eq(:)

  !-- Arrays depending on l and m:
  COMPLEX(kind=8),ALLOCATABLE :: dPhi(:)
  COMPLEX(kind=8),ALLOCATABLE :: dPhi0(:)
  COMPLEX(kind=8),ALLOCATABLE :: dPhi02(:)
  REAL(kind=8),ALLOCATABLE :: dLh(:)
  REAL(kind=8),ALLOCATABLE :: dTheta1S(:),dTheta1A(:)
  REAL(kind=8),ALLOCATABLE :: dTheta2S(:),dTheta2A(:)
  REAL(kind=8),ALLOCATABLE :: dTheta3S(:),dTheta3A(:)
  REAL(kind=8),ALLOCATABLE :: dTheta4S(:),dTheta4A(:)
  REAL(kind=8),ALLOCATABLE :: D_m(:),D_l(:),D_lP1(:)
  REAL(kind=8),ALLOCATABLE :: D_mc2m(:)
  REAL(kind=8),ALLOCATABLE :: hdif_B(:),hdif_V(:),hdif_S(:)

  !-- Limiting l for a given m, used in legtf
  INTEGER,ALLOCATABLE :: lStart(:),lStop(:)
  INTEGER,ALLOCATABLE :: lStartP(:),lStopP(:)
  LOGICAL,ALLOCATABLE :: lmOdd(:),lmOddP(:)

  !----------------------------------------------------------------------
CONTAINS
  SUBROUTINE initialize_horizontal_data

    ALLOCATE( n_theta_cal2ord(n_theta_max) )
    ALLOCATE( theta(n_theta_max) )
    ALLOCATE( theta_ord(n_theta_max) )
    ALLOCATE( sn2(n_theta_max/2) )
    ALLOCATE( osn2(n_theta_max/2) )
    ALLOCATE( cosn2(n_theta_max/2) )
    ALLOCATE( osn1(n_theta_max/2) )
    ALLOCATE( O_sin_theta(n_theta_max) )
    ALLOCATE( O_sin_theta_E2(n_theta_max) )
    ALLOCATE( sinTheta(n_theta_max) )
    ALLOCATE( cosTheta(n_theta_max) )

    !-- Phi (longitude)
    ALLOCATE( phi(n_phi_max) )

    !-- Legendres:
    ALLOCATE( Plm(lm_max,n_theta_max/2) )
    ALLOCATE( wPlm(lmP_max,n_theta_max/2) )
    ALLOCATE( dPlm(lm_max,n_theta_max/2) )
    ALLOCATE( gauss(n_theta_max) )
    ALLOCATE( dPl0Eq(l_max+1) )

    !-- Arrays depending on l and m:
    ALLOCATE( dPhi(lm_max) )
    ALLOCATE( dPhi0(lm_max) )
    ALLOCATE( dPhi02(lm_max) )
    ALLOCATE( dLh(lm_max) )
    ALLOCATE( dTheta1S(lm_max),dTheta1A(lm_max) )
    ALLOCATE( dTheta2S(lm_max),dTheta2A(lm_max) )
    ALLOCATE( dTheta3S(lm_max),dTheta3A(lm_max) )
    ALLOCATE( dTheta4S(lm_max),dTheta4A(lm_max) )
    ALLOCATE( D_m(lm_max),D_l(lm_max),D_lP1(lm_max) )
    ALLOCATE( D_mc2m(n_m_max) )
    ALLOCATE( hdif_B(lm_max),hdif_V(lm_max),hdif_S(lm_max) )

    !-- Limiting l for a given m, used in legtf
    ALLOCATE( lStart(n_m_max),lStop(n_m_max) )
    ALLOCATE( lStartP(n_m_max),lStopP(n_m_max) )
    ALLOCATE( lmOdd(n_m_max),lmOddP(n_m_max) )

  END SUBROUTINE initialize_horizontal_data
!***************************************************************
  SUBROUTINE horizontal
!----------------------------------------------------------------
!  Calculates functions of theta and phi, for exmample the
!  Legendre functions, and functions of degree l and order m
!  of the legendres.
!----------------------------------------------------------------

  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE logic
  USE output_data
  USE plms_theta
#if (FFTLIB==JW)
  USE fft_JW, only: init_fft
#elif (FFTLIB==MKL)
  USE fft_MKL, only: init_fft
#endif

!-- LOCAL VARIABLES:
    INTEGER :: norm,n_theta,n_phi
    INTEGER :: l,m,lm,lmP,mc
    REAL(kind=8) :: pi2
    REAL(kind=8) :: ampnu!,q0
    REAL(kind=8) :: clm(0:l_max+1,0:l_max+1)
    REAL(kind=8) :: plma(lmP_max)
    REAL(kind=8) :: dtheta_plma(lmP_max)
    REAL(kind=8) :: colat
    REAL(kind=8) :: fac
    REAL(kind=8) :: Pl0Eq(l_max+1)

     
!-- End of declaration
!-------------------------------------------------------------------------


    pi2=8.D0*DATAN(1.D0)
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
            
    DO n_theta=1,n_theta_max/2  ! Loop over colat in NHS

        colat=theta_ord(n_theta)

    !----- plmtheta calculates plms and their derivatives
    !      up to degree and order l_max+1 and m_max at
    !      the points cos(theta_ord(n_theta)):
        CALL plm_theta(colat,l_max+1,m_max,minc, &
                       plma,dtheta_plma,lmP_max,norm)
        DO lmP=1,lmP_max
            l=lmP2l(lmP)
            IF ( l <= l_max ) THEN
                lm=lmP2lm(lmP)
                Plm(lm,n_theta) =plma(lmP)
                dPlm(lm,n_theta)=dtheta_plma(lmP)
            END IF
            wPlm(lmP,n_theta)=pi2*gauss(n_theta)*plma(lmP)
        END DO

        ! Get dP for all degrees and order m=0 at the equator only
        ! Usefull to estimate the flow velocity at the equator
        CALL plm_thetaAS(pi2/4.d0,l_max,Pl0Eq,dPl0Eq,l_max+1,norm)

    !----- More functions stored to obscure the code:
        sn2(n_theta)               =DSIN(colat)**2
        osn1(n_theta)              =1.D0/DSIN(colat)
        osn2(n_theta)              =osn1(n_theta)*osn1(n_theta)
        cosn2(n_theta)             =DCOS(colat)*osn2(n_theta)
        O_sin_theta(2*n_theta-1)   =1.D0/DSIN(colat)
        O_sin_theta(2*n_theta  )   =1.D0/DSIN(colat)
        O_sin_theta_E2(2*n_theta-1)=1.D0/(DSIN(colat)*DSIN(colat))
        O_sin_theta_E2(2*n_theta  )=1.D0/(DSIN(colat)*DSIN(colat))
        sinTheta(2*n_theta-1)      =DSIN(colat)
        sinTheta(2*n_theta  )      =DSIN(colat)
        cosTheta(2*n_theta-1)      =DCOS(colat)
        cosTheta(2*n_theta  )      =-DCOS(colat)
                  
    END DO


!-- Resort thetas in the alternating north/south order they
!   are used for the calculations:
    DO n_theta=1,n_theta_max/2
        n_theta_cal2ord(2*n_theta-1)=n_theta
        n_theta_cal2ord(2*n_theta)  =n_theta_max-n_theta+1
        theta(2*n_theta-1)          =theta_ord(n_theta)
        theta(2*n_theta)            =theta_ord(n_theta_max-n_theta+1)
    END DO
       

!----- Same for longitude output grid:
    fac=pi2/DBLE(n_phi_max*minc)
    DO n_phi=1,n_phi_max
        phi(n_phi)=fac*DBLE(n_phi-1)
    END DO


!-- Initialize fast fourier transform for phis:
    !CALL init_fft(n_phi_max,i_fft_init,nDi_fft, &
    !              d_fft_init,nDd_fft)
    CALL init_fft(n_phi_max)

!-- Build arrays depending on degree l and order m
!   and hyperdiffusion factors:

    DO m=0,m_max,minc  ! Build auxiliary array clm
        DO l=m,l_max+1
            clm(l,m)=DSQRT( DBLE((l+m)*(l-m)) / &
                            DBLE((2*l-1)*(2*l+1)) )
        END DO
    END DO

    DO lm=1,lm_max
        l=lm2l(lm)
        m=lm2m(lm)

    !---- Help arrays:
        D_l(lm)  =DBLE(l)
        D_lP1(lm)=DBLE(l+1)
        D_m(lm)  =DBLE(m)

    !---- Operators for derivatives:

    !------- Phi derivate:
        dPhi(lm)=CMPLX(0.D0,DBLE(m),KIND=KIND(0d0))
        IF ( l < l_max ) THEN
            dPhi0(lm)    =CMPLX(0.D0,DBLE(m),KIND=KIND(0d0))
            dPhi02(lm)   =dPhi0(lm)*dPhi0(lm)
        ELSE
            dPhi0(lm)    =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dPhi02(lm)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END IF
    !------- Negative horizontal Laplacian *r^2
        dLh(lm)     =DBLE(l*(l+1))                 ! = qll1
    !------- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 )
        dTheta1S(lm)=DBLE(l+1)        *clm(l,m)    ! = qcl1
        dTheta1A(lm)=DBLE(l)          *clm(l+1,m)  ! = qcl
    !------- Operator ( sin(thetaR) * d/d theta )
        dTheta2S(lm)=DBLE(l-1)        *clm(l,m)    ! = qclm1
        dTheta2A(lm)=DBLE(l+2)        *clm(l+1,m)  ! = qcl2
    !------- Operator ( sin(theta) * d/d theta + cos(theta) dLh )
        dTheta3S(lm)=DBLE((l-1)*(l+1))*clm(l,m)    ! = q0l1lm1(lm)
        dTheta3A(lm)=DBLE(l*(l+2))    *clm(l+1,m)  ! = q0cll2(lm)
    !------- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 ) * dLh
        dTheta4S(lm)=dTheta1S(lm)*DBLE((l-1)*l)
        dTheta4A(lm)=dTheta1A(lm)*DBLE((l+1)*(l+2))

    !--- Hyperdiffusion
        hdif_B(lm)=1.D0
        hdif_V(lm)=1.D0
        hdif_S(lm)=1.D0
        IF ( ldifexp > 0 ) THEN

            IF ( ldif >= 0 .AND. l > ldif ) THEN

            !-- Kuang and Bloxham type:
            !                 hdif_B(lm)=
            !     *                   1.D0+difeta*DBLE(l+1-ldif)**ldifexp
            !                 hdif_V(lm)=
            !     *                   1.D0+ difnu*DBLE(l+1-ldif)**ldifexp
            !                 hdif_S(lm)=
            !     &                   1.D0+difkap*DBLE(l+1-ldif)**ldifexp

            !-- Old type:
                hdif_B(lm)= 1.D0 + difeta * ( &
                                             DBLE(l+1-ldif) / &
                                             DBLE(l_max+1-ldif) )**ldifexp
                hdif_V(lm)= 1.D0 + difnu * ( &
                                             DBLE(l+1-ldif) / &
                                             DBLE(l_max+1-ldif) )**ldifexp
                hdif_S(lm)= 1.D0 + difkap * ( &
                                             DBLE(l+1-ldif) / &
                                             DBLE(l_max+1-ldif) )**ldifexp

            ELSE IF ( ldif < 0 ) THEN

            !-- Grote and Busse type:
                hdif_B(lm)= &
                            (1.D0+difeta*DBLE(l)**ldifexp ) / &
                            (1.D0+difeta*DBLE(-ldif)**ldifexp )
                hdif_V(lm) = &
                            (1.D0+difnu*DBLE(l)**ldifexp ) / &
                            (1.D0+difnu*DBLE(-ldif)**ldifexp )
                hdif_S(lm)= &
                            (1.D0+difkap*DBLE(l)**ldifexp ) / &
                            (1.D0+difkap*DBLE(-ldif)**ldifexp )
                             
            END IF

        ELSE

            IF ( l == l_max .AND. .NOT. l_non_rot ) THEN
            !  Chose ampnu so that the viscous force is at least as
            !  strong as the viscous force for l=l_max:
            !  We can turn this argument around and state that
            !  for example for Ek=1e-4 l_max should be 221.
                ampnu=(r_cmb**2/DBLE(l_max*(l_max+1)))*(2.D0/ek)
                ampnu=DMAX1(1.D0,ampnu)
                hdif_V(lm)=ampnu*hdif_V(lm)
            END IF

        END IF

    END DO ! lm


!-- Build auxiliary index arrays for Legendre transform:
!   lStartP, lStopP give start and end positions in lmP-block.
!   lStart, lStop give start and end positions in lm-block.
    lStartP(1)=1
    lStopP(1) =l_max+2
    lStart(1) =1
    lStop(1)  =l_max+1
    D_mc2m(1)=0
    IF ( MOD(l_max,2) == 0 ) THEN
        lmOdd(1) =.TRUE.
        lmOddP(1)=.FALSE.
    ELSE
        lmOdd(1) =.FALSE.
        lmOddP(1)=.TRUE.
    ENDIF
    DO mc=2,n_m_max
        m=(mc-1)*minc
        D_mc2m(mc) =DBLE(m)
        lStartP(mc)=lStopP(mc-1)+1
        lStopP(mc) =lStartP(mc) +l_max-m+1
        lStart(mc) =lStop(mc-1) +1
        lStop(mc)  =lStart(mc)  +l_max-m
        IF ( MOD(lStop(mc)-lStart(mc),2) == 0 ) THEN
            lmOdd(mc) =.TRUE.
            lmOddP(mc)=.FALSE.
        ELSE
            lmOdd(mc) =.FALSE.
            lmOddP(mc)=.TRUE.
        END IF
    END DO

#if 0
    WRITE(*,"(A,I6,A,I10,A,I10,A)") "cache info of first element, all have dimension 1:lm_max=",lm_max,&
         &" = ",lm_max*16,"B (C), ",lm_max*8,"B (R)"
    CALL print_cache_info_dcmplx("C: dPhi0"//C_NULL_CHAR,dPhi0(1))
    CALL print_cache_info_dcmplx("C: dPhi"//C_NULL_CHAR,dPhi(1))
    CALL print_cache_info_dreal("R: dTheta1A"//C_NULL_CHAR,dTheta1A(1))
    CALL print_cache_info_dreal("R: dTheta1S"//C_NULL_CHAR,dTheta1S(1))
    CALL print_cache_info_dreal("R: dTheta2A"//C_NULL_CHAR,dTheta2A(1))
    CALL print_cache_info_dreal("R: dTheta2S"//C_NULL_CHAR,dTheta2S(1))
    CALL print_cache_info_dreal("R: dTheta3A"//C_NULL_CHAR,dTheta3A(1))
    CALL print_cache_info_dreal("R: dTheta3S"//C_NULL_CHAR,dTheta3S(1))
    CALL print_cache_info_dreal("R: dTheta4A"//C_NULL_CHAR,dTheta4A(1))
    CALL print_cache_info_dreal("R: dTheta4S"//C_NULL_CHAR,dTheta4S(1))
#endif

    RETURN
    end SUBROUTINE horizontal

!---------------------------------------------------------------------------
END MODULE horizontal_data
