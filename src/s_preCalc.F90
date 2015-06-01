!$Id$
!***********************************************************************
SUBROUTINE preCalc
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/15/02  by JW. --------------!

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to initialize the calc values,     |
  !  |  arrays, constants that are used all over the code.               |
  !  |  The stuff is stored in the common blocks.                        |
  !  |  MPI: This is called by every processors.                         |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters,ONLY:nVarEps,pr,prmag,ra,rascaled,ek,ekscaled, &
       & opr,opm,o_sr,radratio,sigma_ratio,CorFac,LFfac,BuoFac,PolInd,   &
       & nVarCond,nVarDiff,nVarVisc,rho_ratio_ic,rho_ratio_ma,epsc,epsc0,&
       & ktops,kbots,interior_model,r_LCR,n_r_LCR
  USE num_param
  USE init_fields
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE const
  USE integration, ONLY: rInt_R

  IMPLICIT NONE

  !---- Local:
  REAL(kind=8) :: sq4pi,c1,help
  REAL(kind=8) :: delmin,sr_top,si_top,sr_bot,si_bot
  INTEGER :: n,n_r,l,m,l_bot,m_bot,l_top,m_top
  CHARACTER(len=76) :: fileName
  CHARACTER(len=80) :: message

  REAL(kind=8) :: mom(n_r_max)

  REAL(kind=8) :: facIH

  !-- end of declaration
  !---------------------------------------------------------------


  !-- Determine scales depending on n_tScale,n_lScale :
  IF ( n_tScale == 0 ) THEN
     !----- Viscous time scale:
     tScale=1.D0
  ELSE IF ( n_tScale == 1 ) THEN
     !----- Magnetic time scale:
     tScale=1.D0/prmag
  ELSE IF ( n_tScale == 2 ) THEN
     !----- Thermal time scale:
     tScale=1.D0/pr
  END IF
  IF ( n_lScale == 0 ) THEN
     !----- Outer Core:
     lScale=1.D0
  ELSE IF ( n_lScale == 1 ) THEN
     !----- Total Core:
     lScale=(1.D0-radratio)
  END IF

  !---- Scale according to scdIFf:
  vScale  =lScale/tScale
  pScale  =tScale**2/lScale**2
  eScale  =vScale*vScale/enscale
  raScaled=ra/lScale**3
  ekScaled=ek*lScale**2

  IF ( l_cond_ic ) O_sr=1.D0/sigma_ratio

  opr=1.D0/pr
  IF ( l_mag ) THEN
     opm=1.D0/prmag
  ELSE
     opm=0.D0
  END IF

  ! Note: CorFac is the factor in front of the Coriolis force. In the scaling
  !       used here (viscous time scale) this is simply the inverse Ekman num:
  !          CorFac=1/Ek
  !       In the non-rotating case there is no Coriolis force and
  !          CorFac=0
  !       LFfac is the factor the Lorentz-force in the Navier-Stokes
  !       equation is multiplied with. This depends on the magnetic
  !       scale chosen. In the standart rotating case the B-scale
  !       is sqrt(\rho \lambda \mu \Omega) where \rho is core density,
  !       \lambda is magnetic diffusivity, \mu is magnetic permeability
  !       and \Omega is rotation rate. In this scaling the square of the
  !       dimensionless magnetic field is identical with the Elasser number:
  !          Els=B^2 / (\rho\lambda\mu\Omega)
  !       The factor to be used in front of the LF is then
  !          LFfac=1/(Ek Pm).
  !       This can also be used to bring kinetic and magnetic energy to
  !       the same scale since Ek=Em/(Ek Pm).
  !       If the system is non rotating we chose the magnetic scale
  !       sqrt(\rho \mu) \nu / L where \nu is the viscous diffusivity
  !       and L=ro-ri is the chosen length scale. Then LFfac=1.
  !       Note that for kinematic dynamos the magnetic field strength
  !       is arbitrary. We nevertheless still use the same LFfac.

  IF ( l_non_rot ) THEN
     CorFac=0.d0
     IF ( l_mag .OR. l_mag_LF ) THEN
        LFfac=1D0
     ELSE
        LFfac=0d0
     END IF
  ELSE
     CorFac=1.D0/ekScaled
     IF ( l_mag .OR. l_mag_LF ) THEN
        LFfac=1.D0/(ekScaled*prmag)
     ELSE
        LFfac=0d0
     END IF
  END IF

  ! Note: BuoFac is the factor used in front of the buoyancy force.
  !       In the scaling used here its
  !          BuoFac=Ra/Pr
  !       where Ra= ( g0 \alpha \delta T L^3)/(\nu \kappa)
  !       with g0 the CMB gravity, \alpha the thermal expansivity,
  !       \delta T the temperature scale, and \kappa the thermal
  !       diffusivity
  BuoFac=raScaled/pr
  IF ( index(interior_model,'JUP') /= 0 ) THEN
     polind=1.d0/0.45d0
  ELSE IF ( index(interior_model,'SAT') /= 0 ) THEN
     polind=1.d0/0.5d0
  ELSE IF ( index(interior_model,'SUN') /= 0 ) THEN
     polind=1.d0/0.6d0
  END IF

  dtStart=dtStart/tScale
  dtMax  =dtMax/tScale
  IF ( .NOT. l_non_rot ) dtMax  =MIN(dtMax,intfac*ekScaled)
  dtMin  =dtMax/1.0D6

  !-- Calculate radial functions for all threads (chebs,r,.....):
  CALL radial
  IF ( ( l_plotmap ).and.(rank.eq.0) ) THEN
     fileName='rNM.'//TAG
     OPEN(99,FILE=fileName,STATUS='UNKNOWN')
     DO n_r=1,n_r_max
        WRITE(99,'(I4,4D12.4)') n_r, &
             r(n_r)-r_icb,drx(n_r),ddrx(n_r),dddrx(n_r)
     END DO
     CLOSE(99)
  END IF

  CALL transportProperties

  IF ( ( l_anel ).and.( rank.eq.0 ) ) THEN
     ! Write the equilibrium setup in anel.TAG
     fileName='anel.'//TAG
     OPEN(99,FILE=fileName,STATUS='UNKNOWN')
     WRITE(99,'(8a15)') 'radius', 'temp0', 'rho0', 'beta', &
         &       'dbeta', 'grav', 'ds0/dr', 'div(k grad T)'
     DO n_r=1,n_r_max
        WRITE(99,'(8e15.7)') r(n_r),1./otemp1(n_r),        &
         &   rho0(n_r),beta(n_r),dbeta(n_r),               &
         &   rgrav(n_r)/BuoFac,dentropy0(n_r),             &
         &   divKtemp0(n_r)
     END DO
     CLOSE(99)
  END IF

  !-- Write radial profiles
  IF ( l_mag .AND. nVarCond > 0 ) THEN
     fileName='varCond.'//TAG
     OPEN(99,FILE=fileName,STATUS='UNKNOWN')
     WRITE(99,'(4a15)') 'radius', 'sigma', 'lambda', 'dLlambda'
     DO n_r=n_r_max,1,-1
        WRITE(99,'(4e15.7)') r(n_r),sigma(n_r),lambda(n_r), &
             dLlambda(n_r)
     END DO
     CLOSE(99)
  END IF

  IF ( ( l_heat .AND. nVarDiff > 0  .OR. nVarVisc > 0).and.( rank.eq.0 ) ) THEN
     fileName='varDiff.'//TAG
     OPEN(99,FILE=fileName,STATUS='UNKNOWN')
     WRITE(99,'(5a15)') 'radius', 'conductivity', 'kappa', &
          'dLkappa', 'Prandtl'
     DO n_r=n_r_max,1,-1
        WRITE(99,'(5D15.7)') r(n_r),kappa(n_r)*rho0(n_r), &
             kappa(n_r),dLkappa(n_r), &
             pr*visc(n_r)/kappa(n_r)
     END DO
     CLOSE(99)
  END IF

  IF ( ( nVarVisc > 0 ).and.(rank.eq.0) ) THEN
     fileName='varVisc.'//TAG
     OPEN(99,FILE=fileName,STATUS='UNKNOWN')
     WRITE(99,'(7a15)') 'radius', 'dynVisc', 'kinVisc', &
          'dLvisc', 'Ekman', 'Prandtl', 'Pm'
     IF ( l_mag ) THEN
        DO n_r=n_r_max,1,-1
           WRITE(99,'(7D15.7)') r(n_r),visc(n_r)*rho0(n_r), &
                visc(n_r),dLvisc(n_r), &
                ek*visc(n_r), &
                pr*visc(n_r)/kappa(n_r), &
                prmag*visc(n_r)/lambda(n_r)
        END DO
     ELSE
        DO n_r=n_r_max,1,-1
           WRITE(99,'(7D15.7)') r(n_r),visc(n_r)*rho0(n_r), &
                visc(n_r),dLvisc(n_r), &
                ek*visc(n_r), &
                pr*visc(n_r)/kappa(n_r), &
                prmag
        END DO
     END IF
     CLOSE(99)
  END IF

  l_LCR=.FALSE.
  n_r_LCR=0
  DO n_r=1,n_r_max
     IF ( r_LCR<=r(n_r)/r_CMB ) THEN
         l_LCR=.TRUE.
         n_r_LCR=n_r
     END IF
  END DO
  IF ( n_r_LCR==1 ) THEN
     l_LCR=.FALSE.
     n_r_LCR=0
  END IF

  !-- Compute some constants:
  !zero  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
  !pi    =4.D0*DATAN(1.D0)
  vol_ic=4.D0*pi/3.D0*r_icb**3             ! Inner core volume
  vol_oc=4.D0*pi/3.D0*(r_cmb**3-r_icb**3)  ! Outer core volume
  surf_cmb=4.D0*pi*r_cmb**2                ! Surface of CMB

  !-- Initialize everything that has to do with the horizontal representation
  !   on all threads:
  CALL horizontal

  !-- Computation of the average density (usefull to compute Re and Rm)
  IF ( l_anel ) THEN
     DO n_r=1,n_r_max
        mom(n_r)=r(n_r)**2 * rho0(n_r)
     END DO
     mass=4.D0*pi/vol_oc* &
          rInt_R(mom,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
  ELSE
     mass=1.D0
  END IF

  !-- Calculate auxiliary arrays containing effective Courant grid intervals:
  c1=1.D0/DBLE(l_max*(l_max+1))
  delxh2(1)      =c1*r_cmb**2
  delxh2(n_r_max)=c1*r_icb**2
  delxr2(1)      =(r(1)-r(2))**2
  delxr2(n_r_max)=(r(n_r_max-1)-r(n_r_max))**2
  do n_r=2,n_r_max-1
     delxh2(n_r)=c1*r(n_r)**2
     delmin=min((r(n_r-1)-r(n_r)),(r(n_r)-r(n_r+1)))
     delxr2(n_r)=delmin*delmin
  END do

  !-- Constants used for rotating core or mantle:
  y10_norm=0.5D0*DSQRT(3.D0/pi)  ! y10=y10_norm * cos(theta)
  y11_norm=0.5D0*DSQRT(1.5D0/pi) ! y11=y11_norm * sin(theta)

  !----- Proportionality factor of (l=1,m=0) toroidal velocity potential
  !      and inner core rotation rate:
  c_z10_omega_ic=y10_norm/(r(n_r_max)*r(n_r_max))

  !----- Proportionality factor of (l=1,m=0) toroidal velocity potential
  !      and mantle rotation rate:
  c_z10_omega_ma=y10_norm/(r(1)*r(1))

  !----- Inner-core normalized moment of inertia:
  c_moi_ic=8.D0*pi/15.D0*r_icb**5*rho_ratio_ic*rho0(n_r_max)

  !----- Outer-core normalized moment of inertia:
  ! _moi_oc=8.D0*pi/15.D0*(r_cmb**5-r_icb**5) ! rho=cst
  DO n_r=1,n_r_max
     mom(n_r)=r(n_r)**4 * rho0(n_r)
  END DO
  c_moi_oc=8.D0*pi/3.D0* &
       rInt_R(mom,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)

  !----- Mantle normalized moment of inertia:
  c_moi_ma=8.D0*pi/15.D0*(r_surface**5-r_cmb**5)*rho_ratio_ma

  !----- IC normalised moment of inertia / r_icb**4 * 3/(8 pi)
  c_dt_z10_ic=0.2D0*r_icb*rho_ratio_ic*rho0(n_r_max)

  !----- Mantle normalised moment of inertia / r_cmb**4 * 3/(8 pi)
  c_dt_z10_ma=0.2D0*r_cmb*rho_ratio_ma * &
       ( (r_surface/r_cmb)**5 - 1.D0 )

  !----- Proportionality factor for ic lorentz_torque as used in
  !      ic torque-equation (z10):
  c_lorentz_ic=0.25D0*DSQRT(3.D0/pi)/(r(n_r_max)*r(n_r_max))

  !----- Proportionality factor for mantle lorentz_torque as used in
  !      mantle torque-equation (z10):
  c_lorentz_ma=0.25D0*DSQRT(3.D0/pi)/(r(1)*r(1))

  !-- Set thermal boundary conditions for fixed temp. on both boundaries:
  !----- Extract tops and bots
  IF ( l_heat ) THEN

     !-- Renormalisation of heat flow coeffs and heat source.
     !      This has been copied from the Christensen and
     !      means the we use a different normalisation of
     !      input tops,bots and epsc0 than in the code.
     !      Spherical harmonics of the input tops, bots and epsc0
     !      are normalised so that the integral of Y*Y(c.c.)
     !      (c.c. stands for conjugate complex)
     !      over a spherical surface of radius 1 gives 4 Pi.
     !      In the code we use totally normalized spherical harmonic,
     !      i.e. the integral of Y*Y(c.c.) over a spherical surface
     !      of radius 1 is unity
     sq4pi=DSQRT(16.D0*DATAN(1.D0))
     epsc=epsc0*sq4pi

     DO m=0,m_max,minc
        DO l=m,l_max
           bots(l,m)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           tops(l,m)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           DO n=1,n_s_bounds
              l_bot =INT(s_bot(4*n-3))
              m_bot =INT(s_bot(4*n-2))
              sr_bot=s_bot(4*n-1)
              si_bot=s_bot(4*n)
              l_top =INT(s_top(4*n-3))
              m_top =INT(s_top(4*n-2))
              sr_top=s_top(4*n-1)
              si_top=s_top(4*n)
              IF ( l_bot == l .AND. m_bot == m .AND. &
                   CMPLX(sr_bot,si_bot,KIND=KIND(0d0)) /= (0.D0,0.D0) &
                   ) THEN
                 IF ( m == 0 ) si_bot=0.D0
                 bots(l,m)=sq4pi*CMPLX(sr_bot,si_bot,KIND=KIND(0d0))
                 IF ( kbots == 2 ) bots(l,m)=bots(l,m)*lScale
              END IF
              IF ( l_top == l .AND. m_top == m .AND. &
                   CMPLX(sr_top,si_top,KIND=KIND(0d0)) /= (0.D0,0.D0) &
                   ) THEN
                 IF ( m == 0 ) si_top=0.D0
                 tops(l,m)=sq4pi*CMPLX(sr_top,si_top,KIND=KIND(0d0))
                 IF ( ktops == 2 ) tops(l,m)=tops(l,m)*lScale
              END IF
           END DO
        END DO
     END DO

     IF ( nVarEps==0 ) THEN
        facIH=vol_oc
     ELSE IF ( nVarEps==1 ) THEN
        facIH=mass*vol_oc
     END IF

     IF ( ktops == 1 .AND. kbots == 1 ) THEN

        tops(0,0)=-r_icb**2/(r_icb**2+r_cmb**2)*sq4pi
        bots(0,0)= r_cmb**2/(r_icb**2+r_cmb**2)*sq4pi

     ELSE IF ( ktops == 2 .AND. kbots == 2 ) THEN

        IF ( REAL(bots(0,0)) > 0.D0 ) THEN
           WRITE(*,*)
           WRITE(*,*) '! NOTE: you have supplied'
           WRITE(*,*) '! s_bot(l=0,m=0)>0 which '
           WRITE(*,*) '! means there is a heat '
           WRITE(*,*) '! flux into the inner core.'
           WRITE(*,*) '! This is unrealistic!'
           WRITE(*,*) '! Use s_bot(l=0,m=0)<0 !'
           STOP
        END IF

        !--- |epsc0|=1 signifies that the heat production rate is used as
        !    temperature scale! I make sure here that the heat flux through
        !    inner and outer boundary, bots and tops, balance the sources.
        IF ( DABS(epsc0) == 1.D0 ) THEN

           !--- Make sure that all the heat comes out at CMB
           IF ( tops(0,0) == 0.D0 ) THEN

              !--- Compensate by flux from ICB:
              !    all over the core :
              IF ( epsc0 >= 0 ) THEN
                 WRITE(*,*) '! NOTE: when the flux through the '
                 WRITE(*,*) '! outer boundary is zero, sinks in'
                 WRITE(*,*) '! the outer core need to balance  '
                 WRITE(*,*) '! the flux from the ICB. Thus we  '
                 WRITE(*,*) '! need epsc<0 !                   '
                 STOP
              END IF
              bots(0,0)=epsc*pr*facIH/(4.d0*pi*r_icb**2 * &
                   rho0(n_r_max)*temp0(n_r_max))
              CALL logWrite( &
                   '! CMB heat flux set to balance volume sources!')

           ELSE IF ( tops(0,0) /= 0.D0 .AND. &
                bots(0,0) == 0.D0 ) THEN

              !--- Correct tops to balance inner sources/sinks:
              IF ( epsc0 <= 0 .AND. bots(0,0) == 0.D0  ) THEN
                 WRITE(*,*) '! NOTE: when the flux through the '
                 WRITE(*,*) '! ICB is zero we need sources in  '
                 WRITE(*,*) '! the outer core which means      '
                 WRITE(*,*) '! epsc0>0.                        '
                 STOP
              END IF
              IF ( DABS(REAL(tops(0,0))) == sq4pi ) &
                   CALL logWrite( &
                   '! You intend to use the CMB flux as buoy. scale??')
              tops(0,0)=-facIH*epsc*pr/(4.D0*pi*r_cmb**2) + &
                   radratio**2*rho0(n_r_max)*temp0(n_r_max)*bots(0,0)
              IF ( tops(0,0) /= help ) CALL logWrite( &
                   '!!!! WARNING: CMB heat flux corrected !!!!')

           ELSE IF ( tops(0,0) /= 0.D0 .AND. &
                bots(0,0) /= 0.D0 ) THEN

              !--- Correct tops to balance inner sources/sinks:
              IF ( epsc0 <= 0 .AND. bots(0,0) == 0.D0  ) THEN
                 WRITE(*,*) '! NOTE: when the flux through the '
                 WRITE(*,*) '! ICB is zero we need sources in  '
                 WRITE(*,*) '! the outer core which means      '
                 WRITE(*,*) '! epsc0>0.                        '
                 STOP
              END IF
              help=4.D0*pi/pr/facIH *            &
                   (r_icb**2*REAL(bots(0,0))*rho0(n_r_max)*temp0(n_r_max) - &
                   r_cmb**2*REAL(tops(0,0)))
              IF ( help /= epsc ) THEN
                 WRITE(*,*) '! NOTE: when flux BC through the '
                 WRITE(*,*) '! ICB and CMB are used the sources '
                 WRITE(*,*) '! have to balance the total flux.'
                 STOP
              END IF

           END IF

        ELSE IF ( epsc0 == 0.D0 .AND.    &
             ( tops(0,0) /= 0.D0 .OR. &
             bots(0,0) /= 0.D0 ) ) THEN
           !--- Correct epsc0 to balance the difference between
           !    flux through the inner and outer boundary:
           epsc=4.D0*pi/pr/facIH *          &
                (r_icb**2*REAL(bots(0,0))*rho0(n_r_max)*temp0(n_r_max) - &
                r_cmb**2*REAL(tops(0,0)))
           CALL logWrite( &
                '! Sources introduced to balance surface heat flux!')
           WRITE(message,'(''!      epsc0='',D16.6)') epsc/sq4pi
           CALL logWrite(message)
        ELSE IF ( epsc0 /= 0.D0 .AND.     &
             tops(0,0) /= 0.D0 .AND. &
             bots(0,0) /= 0.D0 ) THEN
           help=4.D0*pi/pr/facIH *          &
                (r_icb**2*REAL(bots(0,0))*rho0(n_r_max)*temp0(n_r_max) - &
                r_cmb**2*REAL(tops(0,0)))
           IF ( help /= epsc ) THEN
              WRITE(*,*) '! NOTE: when flux BC through the '
              WRITE(*,*) '! ICB and/or CMB is used the sources '
              WRITE(*,*) '! have to balance it.'
              STOP
           END IF
        END IF

     END IF

     IF ( l_anelastic_liquid ) THEN
        bots(0,0)=0.D0
        tops(0,0)=0.D0
        epsc=0.D0
        !epsc=4.D0*pi/pr/epsS/vol_oc *          &
        !     (r_icb**2*dtemp0(n_r_max)*rho0(n_r_max)*kappa(n_r_max) - &
        !      r_cmb**2*dtemp0(1)*rho0(1)*kappa(1))*sq4pi
     END IF
     IF ( ktops == 1 ) THEN
        WRITE(message,                                 &
             &  '(''! Constant temp. at CMB T='',D16.6)') &
             &   REAL(tops(0,0))/sq4pi
        CALL logWrite(message)
     ELSE IF ( ktops == 2 ) THEN
        help=surf_cmb*REAL(tops(0,0))/sq4pi
        WRITE(message,                                         &
             &  '(''! Const. total CMB buoy. flux    ='',D16.6)') &
             &  help
        CALL logWrite(message)
     END IF
     IF ( kbots == 1 ) THEN
        WRITE(message,                                 &
             &  '(''! Constant temp. at ICB T='',D16.6)') &
             &  REAL(bots(0,0))/sq4pi
        CALL logWrite(message)
     ELSE IF ( kbots == 2 ) THEN
        help=surf_cmb*radratio**2*rho0(n_r_max)*temp0(n_r_max)*REAL(bots(0,0))/sq4pi
        WRITE(message,                                         &
             &  '(''! Const. total ICB buoy. flux    ='',D16.6)') &
             &  help
        CALL logWrite(message)
     END IF
     help=facIH*pr*epsc/sq4pi
     WRITE(message,'(''! Total vol. buoy. source ='',D16.6)') &
          &  help
     CALL logWrite(message)

  END IF

  RETURN
end SUBROUTINE preCalc


!---------------------------------------------------------------------------
!-- end of subroutine preCalc
!---------------------------------------------------------------------------
