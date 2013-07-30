!$Id$
!***********************************************************************
    SUBROUTINE initV(w,z,omega_ic,omega_ma,lmStart,lmStop)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to initialize the velocity field   |
!  |  So far it is only rudimentary and will be expanded later.        |
!  |  Because s is neede for dwdt init_s has to be called before.      |
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
    USE logic
    USE const
    use parallel_mod
    USE usefull, only: random
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
    IMPLICIT NONE

!-- input:
! include 'truncation.f'
! include 'c_blocking.f'
! include 'c_phys_param.f'
! include 'c_num_param.f'
! include 'c_radial.f'
! include 'c_horizontal.f'
! include 'c_init_fields.f'
! include 'c_const.f'
! include 'c_logic.f'
    INTEGER :: lmStart,lmStop

!--output:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    REAL(kind=8) :: omega_ic,omega_ma

!-- local:
    INTEGER :: lm,lmP,lmPS,lmPA,l,m,n,lmStart00
    INTEGER :: nR,nTheta,nThetaB,nThetaStart,nPhi
    REAL(kind=8) :: ra1,ra2,c_r,c_i
    REAL(kind=8) :: amp_r,rExp,rr
    REAL(kind=8) :: rDep(n_r_max)
    INTEGER :: l1m0

    REAL(kind=8) :: ss,ome(nrp,nfs)
    COMPLEX(kind=8) :: omeLM(lmP_max)

!-- end of declaration
!----------------------------------------------------------------------
            
    lmStart00=MAX(lmStart,1)
    l1m0=lm2(1,0)

!-- Initialize rotation according to
!   given inner core and mantel rotation rate:
    IF ( init_v1 == 1 .AND. &
       ( omega_ic1 /= 0.D0 .OR. omega_ma1 /= 0.D0 ) ) THEN

    !-- Approximating the Stewardson solution:
        DO nR=1,n_r_max

            nTheta=0
            DO n=1,nThetaBs ! loop over the theta blocks

                nThetaStart=(n-1)*sizeThetaB+1
                DO nThetaB=1,sizeThetaB
                    nTheta=nTheta+1
                    ss=r(nR)*sinTheta(nTheta)
                !------------ start with constructing rotation rate ome:
                    DO nPhi=1,n_phi_max
                        IF ( ss <= r_icb ) THEN
                            ome(nPhi,nThetaB)=omega_ma1+0.5D0*omega_ic1
                        ELSE
                            ome(nPhi,nThetaB)=omega_ma1
                        END IF
                    END DO
                END DO
            !------------ Transform to spherical hamonic space for each theta block
                CALL fft_thetab(ome,-1)
                CALL legTF1(nThetaStart,omeLM,ome)
            END DO ! End of loop over theta blocks

        !------------ ome now in spherical harmonic space,
        !             apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
        !             additional application of r**2/(l*(l+1)) then yields
        !             the axisymmetric toriodal flow contribution:
            DO lm=lmStart00,lmStop
                l   =lm2l(lm)
                m   =lm2m(lm)
                lmP =lm2lmP(lm)
                lmPS=lmP2lmPS(lmP)   ! l-1
                lmPA=lmP2lmPA(lmP)   ! l+1
                IF ( l > m ) THEN
                    z(lm,nR)=z(lm,nR) + &
                                r(nR)**2/dLh(lm) * ( &
                              dTheta1S(lm)*omeLM(lmPS) - &
                              dTheta1A(lm)*omeLM(lmPA) )
                ELSE IF ( l == m ) THEN
                    z(lm,nR)=z(lm,nR) - &
                                r(nR)**2/dLh(lm) * &
                              dTheta1A(lm)*omeLM(lmPA)
                END IF
            END DO

        END DO ! close loop over radial grid points

    ELSE IF ( init_v1 > 1 ) THEN

    !--- Add random noise toroidal field of all (l,m) modes exept (l=0,m=0):
    !    It decays likes l**(init_v1-1)
    !    Amplitude is chosen so that the (1,0) term resembles amp_v1 *
    !    the 'solid body' rotation set by inner core and mantle rotation.

        rr=random(1.d0)
        rExp=4.
        IF ( omega_ic1 /= 0 ) THEN
            amp_r=amp_v1*omega_ic1*r_ICB**(rExp+1.)/y10_norm
        ELSE
            amp_r=amp_v1*r_ICB**(rExp+1.)/y10_norm
        END IF
        DO nR=1,n_r_max
            rDep(nR)=amp_r/r(nR)**(rExp-1.)
            DO lm=lmStart00,lmStop
                m=lm2m(lm)
                ra1=(-1.d0+2.d0*random(0.d0))/D_l(lm)**(init_v1-1)
                ra2=(-1.d0+2.d0*random(0.d0))/D_l(lm)**(init_v1-1)
                c_r=ra1*rDep(nR)
                c_i=ra2*rDep(nR)
                IF ( m == 0 ) THEN  ! non axisymmetric modes
                    z(lm,nR)=z(lm,nR)+CMPLX(c_r,0.D0,KIND=KIND(0d0))
                ELSE
                    z(lm,nR)=z(lm,nR)+CMPLX(c_r,c_i,KIND=KIND(0d0))
                END IF
            END DO
        END DO

    ELSE IF ( init_v1 < -1 ) THEN

    !--- Add random noise poloidal field of all (l,m) modes exept (l=0,m=0):
    !    It decays likes l**(init_v1-1)
    !    Amplitude is chosen to be comparable to amp * inner core roation speed
    !    at inner core boundary...
        rr=random(1.d0)
        IF ( omega_ic1 /= 0.D0 ) THEN
            amp_r=amp_v1*omega_ic1*r_icb*r_icb/(y10_norm*PI)
        ELSE
            amp_r=amp_v1
        ENDIF
        DO nR=1,n_r_max
            rDep(nR)=-amp_r*DSIN( (r(nR)-r_ICB)*PI )
            DO lm=lmStart00,lmStop
                m=lm2m(lm)
                ra1=(-1.d0+2.d0*random(0.d0))/D_l(lm)**(-init_v1-1)
                ra2=(-1.d0+2.d0*random(0.d0))/D_l(lm)**(-init_v1-1)
                c_r=ra1*rDep(nR)
                c_i=ra2*rDep(nR)
                IF ( m > 0 ) THEN  ! no axisymmetric modes
                    w(lm,nR)=w(lm,nR)+CMPLX(c_r,c_i,KIND=KIND(0d0))
                    z(lm,nR)=z(lm,nR)+CMPLX(c_r,c_i,KIND=KIND(0d0))
                END IF
            END DO
        END DO

    END IF

!----- Caring for IC and mantle rotation rates if this
!      has not been done alread in read_start_file.f:
    IF ( lmStart == 1 ) THEN ! Selects one processor !
        IF ( .NOT. l_start_file ) THEN
            WRITE(*,*) '! NO STARTFILE READ, SETTING Z10!'
            IF ( l_SRIC .OR. &
            l_rot_ic .AND. omega_ic1 /= 0.d0 ) THEN
                omega_ic=omega_ic1*DCOS(omegaOsz_ic1*tShift_ic1) + &
                         omega_ic2*DCOS(omegaOsz_ic2*tShift_ic2)
                WRITE(*,*)
                WRITE(*,*) '! I use prescribed inner core rotation rate:'
                WRITE(*,*) '! omega_ic=',omega_ic
                z(l1m0,n_r_icb)=CMPLX(omega_ic/c_z10_omega_ic,KIND=KIND(0d0))
            ELSE IF ( l_rot_ic .AND. omega_ic1 == 0.D0 ) THEN
                omega_ic=c_z10_omega_ic*REAL(z(l1m0,n_r_icb))
            ELSE
                omega_ic=0.D0
            END IF
            IF ( l_SRMA .OR. &
            l_rot_ma .AND. omega_ma1 /= 0.d0 ) THEN
                omega_ma=omega_ma1*DCOS(omegaOsz_ma1*tShift_ma1) + &
                         omega_ma2*DCOS(omegaOsz_ma2*tShift_ma2)
                WRITE(*,*)
                WRITE(*,*) '! I use prescribed mantle rotation rate:'
                WRITE(*,*) '! omega_ma=',omega_ma
                z(l1m0,n_r_cmb)=CMPLX(omega_ma/c_z10_omega_ma,KIND=KIND(0d0))
            ELSE IF ( l_rot_ma .AND. omega_ma1 == 0.D0 ) THEN
                omega_ma=c_z10_omega_ma*REAL(z(l1m0,n_r_cmb))
            ELSE
                omega_ma=0.D0
            END IF
        ELSE
            IF ( nRotIc == 2 ) omega_ic=omega_ic1
            IF ( nRotMa == 2 ) omega_ma=omega_ma1
        END IF
    END IF ! lmStart=1 ?


    RETURN
    end SUBROUTINE initV

!--------------------------------------------------------------------
