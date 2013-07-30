!$Id$
!***********************************************************************
    SUBROUTINE initS(s,lmStart,lmStop)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to initialize the entropy field    |
!  |  according to the input control parameters.                       |
!  |  For init_s1 < 100: random noise initialized                      |
!  |                  the noise spectrum decays as l ^ (init_s1-1)     |
!  |                  with peak amplitude amp_s1  for l=1              |
!  |      init_s1 >=100: a specific harmonic mode initialized          |
!  |              with amplitude amp_s1.                               |
!  |              init_s1 is interpreted as number llmm                |
!  |              where ll: harmonic degree, mm: harmonic order.       |
!  |      init_s2 >100 : a second harmonic mode initialized            |
!  |              with amplitude amp_s2.                               |
!  |              init_s2 is again interpreted as number llmm          |
!  |              where ll: harmonic degree, mm: harmonic order.       |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE init_fields
    USE blocking
    USE horizontal_data
    USE const
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
    USE usefull, ONLY: random
    USE algebra, ONLY: sgesl,sgefa

    IMPLICIT NONE

    INTEGER :: lmStart,lmStop

!-- output:
    COMPLEX(kind=8) :: s(lm_max,n_r_max)

!-- local:
    INTEGER :: n_r,lm,l,m,lm00,lmMin
    REAL(kind=8) :: x,rr,c_r,c_i,s_r,s_i
    REAL(kind=8) :: ra1,ra2
    REAL(kind=8) :: s0(n_r_max),s1(n_r_max)

    INTEGER :: nTheta,n,nThetaStart,nThetaB,nPhi,nS
    REAL(kind=8) :: xL,yL,zL,rH,angleL,s00,s00P
    REAL(kind=8) :: mata(n_impS_max,n_impS_max)
    REAL(kind=8) :: amp(n_impS_max)
    INTEGER :: pivot(n_impS_max)
    REAL(kind=8) :: xS(n_impS_max),yS(n_impS_max)
    REAL(kind=8) :: zS(n_impS_max),sFac(n_impS_max)
    REAL(kind=8) :: sCMB(nrp,nfs)
    COMPLEX(kind=8) :: sLM(lmP_max)
    INTEGER :: info,i,j

!-- end of declaration
!----------------------------------------------------------------------
     
!--  Calculate conductive solution for entropy; (l=0,m=0)-mode
!    (returned in array s0)

    lm00=lm2(0,0)
    lmMin=MAX(lmStart,2)

    IF ( .NOT. l_start_file ) THEN

        IF ( lmStart <= lm00 .AND. lmStop >= lm00 ) THEN
            CALL s_cond(s0)

        !--- Initialize (l=0,m=0)-mode with s0:
            open(unit=999, file='scond.dat')
            DO n_r=1,n_r_max
                s(1,n_r)=s0(n_r)
                write(999,*) r(n_r), s0(n_r)/SQRT(4.*DACOS(-1.d0))
            END DO
            close(999)
        END IF

    END IF

!-- Radial dependence of perturbation in s1:
    DO n_r=1,n_r_max
        x=2.d0*r(n_r)-r_cmb-r_icb
        s1(n_r)=1.d0-3.d0*x**2+3.d0*x**4-x**6
    END DO

    IF ( init_s1 < 100 .AND. init_s1 > 0 ) THEN

    !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
         
        rr=random(1.d0)
        DO lm=lmMin,lmStop
            ra1=(-1.d0+2.d0*random(0.d0))*amp_s1/D_l(lm)**(init_s1-1)
            ra2=(-1.d0+2.d0*random(0.d0))*amp_s1/D_l(lm)**(init_s1-1)
            DO n_r=1,n_r_max
                c_r=ra1*s1(n_r)
                c_i=ra2*s1(n_r)
                IF ( lm2m(lm) > 0 ) THEN  ! non axisymmetric modes
                    s(lm,n_r)=s(lm,n_r)+CMPLX(c_r,c_i,KIND=KIND(0d0))
                ELSE
                    s(lm,n_r)=s(lm,n_r)+CMPLX(c_r,0.d0,KIND=KIND(0d0))
                END IF
            END DO
        END DO
         
    ELSE  IF ( init_s1 >= 100 ) THEN

    !-- Initialize one or two modes specifically

    !----- Initialize first mode:
        l=init_s1/100
        if ( l.gt.99 ) l=init_s1/1000
        m=mod(init_s1,100)
        if ( l.gt.99 ) m=mod(init_s1,1000)
        IF ( mod(m,minc) /= 0 ) THEN
            WRITE(*,*) &
                '! Wave number of mode for entropy initialisation'
            WRITE(*,*) &
                '! not compatible with phi-symmetry:',m
            STOP
        END IF
        IF ( l > l_max .OR. l < m ) THEN
            WRITE(*,*) '! Degree of mode for entropy initialisation'
            WRITE(*,*) '! > l_max or < m !',l
            STOP
        END IF
        lm=lm2(l,m)

        IF ( lmMin <= lm .AND. lmStop >= lm ) THEN
            DO n_r=1,n_r_max
                c_r=s1(n_r)*amp_s1
                s(lm,n_r)=s(lm,n_r)+CMPLX(c_r,0.d0,KIND=KIND(0d0))
            END DO

            WRITE(*,'(/'' ! Entropy initialized at mode:'', &
              &  '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_s1
        END IF

    !----- Initialize second mode:
        IF ( init_s2 > 99 ) THEN
            m=mod(init_s2,100)
            IF ( mod(m,minc) /= 0 ) THEN
                WRITE(*,*) &
                    '! Wave number of mode for entropy initialisation'
                WRITE(*,*) &
                    '! not compatible with phi-symmetry:',m
                STOP
            END IF
            l=init_s2/100
            IF ( l > l_max .OR. l < m ) THEN
                WRITE(*,*) &
                    '! Degree of mode for entropy initialisation'
                WRITE(*,*) &
                    '! > l_max or < m !',l
                STOP
            END IF

            lm=lm2(l,m)
            IF ( lmMin <= lm .AND. lmStop >= lm ) THEN
                s_r=amp_s2
                s_i=0.d0
                IF ( amp_s2 < 0.d0 .AND. m /= 0 ) THEN
                !-------- Sin(phi)-mode initialized for amp_s2<0
                    s_r=0.d0
                    s_i=amp_s2
                END IF
                DO n_r=1,n_r_max
                    c_r=s1(n_r)*s_r
                    c_i=s1(n_r)*s_i
                    s(lm,n_r)=s(lm,n_r)+CMPLX(c_r,c_i,KIND=KIND(0d0))
                END DO
                WRITE(6,'('' ! Second mode:'', &
                  &  '' l='',i3,'' m='',i3,'' Ampl='',f8.5/)') l,m,amp_s2
            END IF

        END IF

    END IF

    IF ( lmStart > lm00 .OR. impS == 0 ) RETURN

!--- Now care for the prescribed boundary condition:

    IF ( minc /= 1 ) THEN
        WRITE(*,*) '! impS doesnt work for minc.NE.1'
        STOP
    END IF

    IF ( ABS(impS) == 1 ) THEN
        n_impS=2
        peakS(2)=-peakS(1)
        thetaS(2)=pi-thetaS(1)
        phiS(2)  =pi+phiS(1)
        IF ( phiS(2) > 2*pi ) phiS(2)=phiS(2)-2*pi
        widthS(2)=widthS(1)
    END IF

!--- Determine the peak value vector in (xS,yS,zS) space.
!       Then get the proportionality factors for the linear dependence
!       of the mean (l=0,m=0) contribution on the total peak amplitude
!       amp:
    DO nS=1,n_impS

        xS(nS)=DSIN(thetaS(nS))*DCOS(phiS(nS))
        yS(nS)=DSIN(thetaS(nS))*DSIN(phiS(nS))
        zS(nS)=DCOS(thetaS(nS))

        nTheta=0
        DO n=1,nThetaBs ! loop over the theta blocks

            nThetaStart=(n-1)*sizeThetaB+1
            DO nThetaB=1,sizeThetaB
                nTheta=nTheta+1
                DO nPhi=1,n_phi_max
                    xL=sinTheta(nTheta)*DCOS(phi(nPhi))
                    yL=sinTheta(nTheta)*DSIN(phi(nPhi))
                    zL=cosTheta(nTheta)
                    rH=DSQRT((xS(nS)-xL)**2 + &
                       (yS(nS)-yL)**2+(zS(nS)-zL)**2)
                !------ Opening angleL with peak value vector:
                    angleL=2.D0*DABS(DASIN(rH/2))
                    IF ( angleL <= widthS(nS) ) THEN
                        sCMB(nPhi,nThetaB)= &
                                   (DCOS(angleL/widthS(nS)*pi)+1)/2
                    ELSE
                        sCMB(nPhi,nThetaB)=0.D0
                    END IF
                END DO
            END DO
        !------ Transform to spherical hamonic space for each theta block
            CALL fft_thetab(sCMB,-1)
            CALL legTF1(nThetaStart,sLM,sCMB)

        END DO ! Loop over theta blocks

    !--- sFac describes the linear dependence of the (l=0,m=0) mode
    !    on the amplitude peakS, SQRT(4*pi) is a normalisation factor
    !    according to the spherical harmonic function form chosen here.
        sFac(nS)=REAL(sLM(lm00))/DSQRT(4*pi)

    END DO ! Loop over peak

!--- Value due to prescribed (l=0,m=0) contribution
    s00P=REAL(tops(0,0))/DSQRT(4*pi)
    IF ( s00P == 0.D0 .AND. impS < 0 ) THEN
        WRITE(*,*) '! No relative amplitudes possible!'
        WRITE(*,*) '! for impS<0 because the mean value!'
        WRITE(*,*) '! is zero! Refince s_top?'
        STOP
    END IF
    IF ( impS > 0 ) s00P=1.D0

!--- Determine the true amplitudes amp for the peaks by solving linear system:
!    These amplitudes guarantee that the peak as an ampliture peakS
!    above or below the mean (l=0,m=0)
    IF ( n_impS == 1 ) THEN
        amp(1)=peakS(1)/(s00P*(1.D0-sFac(1)))
    ELSE
        DO j=1,n_impS
            amp(j)=-peakS(j)/s00P
            DO i=1,n_impS
                IF ( i == j ) THEN
                    mata(i,j)=sFac(i)-1
                ELSE
                    mata(i,j)=sFac(i)
                END IF
            END DO
        END DO
        CALL sgefa(mata,n_impS_max,n_impS,pivot,info)
        CALL sgesl(mata,n_impS_max,n_impS,pivot,amp)
    END IF
    s00=0.D0
    DO nS=1,n_impS
        s00=s00+sFac(nS)*amp(nS)
    END DO

!--- Now get the total thing so that the mean (l=0,m=0) due
!    to the peaks is zero. The (l=0,m=0) contribution is
!    determined (prescribed) by other means.
    nTheta=0
    DO n=1,nThetaBs ! loop over the theta blocks

        nThetaStart=(n-1)*sizeThetaB+1
        DO nThetaB=1,sizeThetaB
            nTheta=nTheta+1
            DO nPhi=1,n_phi_max
                xL=sinTheta(nTheta)*DCOS(phi(nPhi))
                yL=sinTheta(nTheta)*DSIN(phi(nPhi))
                zL=cosTheta(nTheta)
                sCMB(nPhi,nThetaB)=-s00
                DO nS=1,n_impS
                    rH=DSQRT((xS(nS)-xL)**2 + (yS(nS)-yL)**2+(zS(nS)-zL)**2)
                !------ Opening angle with peak value vector:
                    angleL=2.D0*DABS(DASIN(rH/2))
                    IF ( angleL <= widthS(nS) )                    &
                           sCMB(nPhi,nThetaB)=sCMB(nPhi,nThetaB) + &
                           amp(nS)*(DCOS(angleL/widthS(nS)*pi)+1)/2
                END DO
            END DO
        END DO
    !------ Transform to spherical hamonic space for each theta block
        CALL fft_thetab(sCMB,-1)
        CALL legTF1(nThetaStart,sLM,sCMB)

    END DO ! Loop over theta blocks


!--- Finally store the boundary condition and care for
!    the fact that peakS provides the relative amplitudes
!    in comparison to the (l=0,m=0) contribution when impS<0:
!    Note that the (l=0,m=0) has to be determined by other means
!    for example by setting: s_top= 0 0 -1 0
    DO m=0,l_max,minc
        DO l=m,l_max
            lm=lmP2(l,m)
            IF ( l <= l_max .AND. l > 0 ) tops(l,m)=tops(l,m)+sLM(lm)
        END DO
    END DO


    RETURN
    end SUBROUTINE initS

!---------------------------------------------------------------------------
