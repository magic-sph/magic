!$Id$
!***********************************************************************
    SUBROUTINE write_rot(time,dt,eKinIC,ekinMA,w,z,dz,b, &
                                      omega_ic,omega_ma, &
                    lorentz_torque_ic,lorentz_torque_ma)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to write information on the        |
!  |  outputfile file_rot.                                             |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data
    USE const
    USE integration, ONLY: rInt

    IMPLICIT NONE

!-- Input of variables:
    REAL(kind=8) :: omega_ic,omega_ma
    REAL(kind=8) :: lorentz_torque_ma,lorentz_torque_ic
    REAL(kind=8) :: time,dt

!-- Input of scalar fields of toroidal flow and poloidal mag. field:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    COMPLEX(kind=8) :: dz(lm_max,n_r_max)
    COMPLEX(kind=8) :: b(lm_maxMag,n_r_maxMag)

!-- Output into rot_file
    REAL(kind=8) :: eKinIC,eKinOC,eKinMA

!-- Local variables:
    INTEGER :: n_r1,n_r2,n_r3,nR
    INTEGER :: l1m0,l1m1
    REAL(kind=8) :: viscous_torque_ic,viscous_torque_ma
    REAL(kind=8) :: AMz,eKinAMz
    REAL(kind=8) :: angular_moment_oc(3)
    REAL(kind=8) :: angular_moment_ic(3)
    REAL(kind=8) :: angular_moment_ma(3)
    COMPLEX(kind=8) :: z10(n_r_max),z11(n_r_max)
    CHARACTER(len=80) :: filename

    REAL(kind=8) :: powerLor,powerVis

    REAL(kind=8) :: AMzLast,eKinAMzLast
    SAVE AMzLast,eKinAMzLast

!-- end of declaration
!-----------------------------------------------------------------------


!-- Calculating viscous torques:
    IF ( l_rot_ic .AND. kbotv == 2 ) THEN
        CALL get_viscous_torque(viscous_torque_ic, &
                                z(2,n_r_max),dz(2,n_r_max),r_icb)
    ELSE
        viscous_torque_ic=0.d0
    END IF
    IF ( l_rot_ma .AND. ktopv == 2 ) THEN
        CALL get_viscous_torque(viscous_torque_ma, &
                                z(2,1),dz(2,1),r_cmb)
    ELSE
        viscous_torque_ma=0.d0
    END IF

    IF ( l_SRIC ) THEN
        powerLor=lorentz_torque_ic*omega_IC
        powerVis=viscous_torque_ic*omega_IC
        OPEN(n_SRIC_file,file=SRIC_file,status="unknown", &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1p,2x,d20.12,4d17.6)') &
                                      time*tScale, &
                                  omega_ic/tScale, &
         (powerLor+powerVis)*vScale*vScale/tScale, &
                    powerVis*vScale*vScale/tScale, &
                    powerLor*vScale*vScale/tScale
        CLOSE(n_SRIC_file)
    END IF
    IF ( l_SRMA ) THEN
        powerLor=lorentz_torque_ma*omega_ma
        powerVis=viscous_torque_ma*omega_ma
        OPEN(n_SRMA_file,file=SRMA_file,status="unknown", &
             POSITION='APPEND')
        WRITE(n_SRMA_file,'(1p,2x,d20.12,4d17.6)') &
                                      time*tScale, &
                                  omega_ma/tScale, &
         (powerLor+powerVis)*vScale*vScale/tScale, &
                    powerVis*vScale*vScale/tScale, &
                    powerLor*vScale*vScale/tScale
        CLOSE(n_SRMA_file)
    END IF


    IF ( l_drift ) THEN
        filename='driftVD.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        n_r1=INT(1.D0/3.D0*(n_r_max-1))
        n_r2=INT(2.D0/3.D0*(n_r_max-1))
        n_r3=n_r_max-1
        WRITE(n_SRIC_file,'(1P,2X,D20.12,24D12.4)') &
                                              time, &
                        z(lm2(minc,minc)    ,n_r1), &
                        z(lm2(2*minc,2*minc),n_r1), &
                        z(lm2(3*minc,3*minc),n_r1), &
                        z(lm2(4*minc,4*minc),n_r1), &
                        z(lm2(minc,minc)    ,n_r2), &
                        z(lm2(2*minc,2*minc),n_r2), &
                        z(lm2(3*minc,3*minc),n_r2), &
                        z(lm2(4*minc,4*minc),n_r2), &
                        z(lm2(minc,minc)    ,n_r3), &
                        z(lm2(2*minc,2*minc),n_r3), &
                        z(lm2(3*minc,3*minc),n_r3), &
                        z(lm2(4*minc,4*minc),n_r3)
        CLOSE(n_SRIC_file)
        filename='driftVQ.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,24D12.4)') &
                                              time, &
                      z(lm2(minc+1,minc)    ,n_r1), &
                      z(lm2(2*minc+1,2*minc),n_r1), &
                      z(lm2(3*minc+1,3*minc),n_r1), &
                      z(lm2(4*minc+1,4*minc),n_r1), &
                      z(lm2(minc+1,minc)    ,n_r2), &
                      z(lm2(2*minc+1,2*minc),n_r2), &
                      z(lm2(3*minc+1,3*minc),n_r2), &
                      z(lm2(4*minc+1,4*minc),n_r2), &
                      z(lm2(minc+1,minc)    ,n_r3), &
                      z(lm2(2*minc+1,2*minc),n_r3), &
                      z(lm2(3*minc+1,3*minc),n_r3), &
                      z(lm2(4*minc+1,4*minc),n_r3)
        CLOSE(n_SRIC_file)
        IF ( l_mag .OR. l_mag_LF ) THEN
           filename='driftBD.'//tag
           OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
                POSITION='APPEND')
           n_r1=n_r_CMB
           n_r2=n_r_ICB
           n_r3=n_r_ic_max/2
           WRITE(n_SRIC_file,'(1P,2X,D20.12,16D12.4)') &
                                                 time, &
                         b(lm2(minc+1,minc)    ,n_r1), &
                         b(lm2(2*minc+1,2*minc),n_r1), &
                         b(lm2(3*minc+1,3*minc),n_r1), &
                         b(lm2(4*minc+1,4*minc),n_r1), &
                         b(lm2(minc+1,minc)    ,n_r2), &
                         b(lm2(2*minc+1,2*minc),n_r2), &
                         b(lm2(3*minc+1,3*minc),n_r2), &
                         b(lm2(4*minc+1,4*minc),n_r2)
           CLOSE(n_SRIC_file)
           filename='driftBQ.'//tag
           OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
                POSITION='APPEND')
           WRITE(n_SRIC_file,'(1P,2X,D20.12,16D12.4)') &
                                                 time, &
                           b(lm2(minc,minc)    ,n_r1), &
                           b(lm2(2*minc,2*minc),n_r1), &
                           b(lm2(3*minc,3*minc),n_r1), &
                           b(lm2(4*minc,4*minc),n_r1), &
                           b(lm2(minc,minc)    ,n_r2), &
                           b(lm2(2*minc,2*minc),n_r2), &
                           b(lm2(3*minc,3*minc),n_r2), &
                           b(lm2(4*minc,4*minc),n_r2)
           CLOSE(n_SRIC_file)
        END IF ! l_mag
    END IF

    IF ( .NOT. l_SRIC .AND. ( l_rot_ic .OR. l_rot_ma ) ) THEN
        IF ( l_save_out ) THEN
            OPEN(n_rot_file,file=rot_file,status='UNKNOWN', &
                 POSITION='APPEND')
        END IF
        WRITE(n_rot_file,'(1P,2X,D20.12,6D14.6)') &
                                     time*tScale, &
                                 omega_ic/tScale, &
              lScale**2*vScale*lorentz_torque_ic, &
              lScale**2*vScale*viscous_torque_ic, &
                                 omega_ma/tScale, &
              lScale**2*vScale*lorentz_torque_ma, &
             -lScale**2*vScale*viscous_torque_ma
        IF ( l_save_out ) CLOSE(n_rot_file)
    END IF

    IF ( l_AM ) THEN
        l1m0=lm2(1,0)
        l1m1=lm2(1,1)
        IF ( l1m1 > 0 ) THEN
           DO nR=1,n_r_max
              z10(nR)=z(l1m0,nR)
              z11(nR)=z(l1m1,nR)
           END DO
        ELSE
           DO nR=1,n_r_max
              z10(nR)=z(l1m0,nR)
              z11(nR)=CMPLX(0d0,0d0,KIND=KIND(0d0))
           END DO
        END IF
        CALL get_angular_moment(z10,z11,omega_ic,omega_ma, &
                                        angular_moment_oc, &
                      angular_moment_ic,angular_moment_ma)
        IF ( l_save_out ) THEN
            OPEN(n_angular_file,FILE=angular_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
        END IF
        AMz=angular_moment_oc(3) + &
            angular_moment_ic(3)+angular_moment_ma(3)
        eKinAMz=0.5d0*(angular_moment_oc(3)**2/c_moi_oc + &
                       angular_moment_ic(3)**2/c_moi_ic + &
                       angular_moment_ma(3)**2/c_moi_ma )
        eKinIC=0.5d0*angular_moment_ic(3)**2/c_moi_ic
        eKinOC=0.5d0*angular_moment_oc(3)**2/c_moi_oc
        eKinMA=0.5d0*angular_moment_ma(3)**2/c_moi_ma
        IF ( AMzLast /= 0.D0 ) THEN
            WRITE(n_angular_file,'(1p,2x,d20.12,5d14.6,7d20.12)') &
                                                     time*tScale, &
                                               angular_moment_oc, &
                                            angular_moment_ic(3), &
                                            angular_moment_ma(3), &
                            AMz,(AMz-AMzLast)/AMzLast/dt,eKinAMz, &
                            (eKinAMz-eKinAMzLast)/eKinAMzLast/dt, &
                                            eKinIC,eKinOC,eKinMA
        END IF
        IF ( l_save_out ) CLOSE(n_angular_file)
        AMzLast=AMz
        eKinAMzLast=eKinAMz
    END IF

    IF ( l_iner ) THEN
        n_r1=INT(1.D0/2.D0*(n_r_max-1))
        filename='inerP.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,21D12.4)') &
                                              time, &
                           REAL(w(lm2(1,1),n_r1)), &
                           REAL(w(lm2(2,1),n_r1)), &
                           REAL(w(lm2(2,2),n_r1)), &
                           REAL(w(lm2(3,1),n_r1)), &
                           REAL(w(lm2(3,2),n_r1)), &
                           REAL(w(lm2(3,3),n_r1)), &
                           REAL(w(lm2(4,1),n_r1)), &
                           REAL(w(lm2(4,2),n_r1)), &
                           REAL(w(lm2(4,3),n_r1)), &
                           REAL(w(lm2(4,4),n_r1)), &
                           REAL(w(lm2(5,1),n_r1)), &
                           REAL(w(lm2(5,2),n_r1)), &
                           REAL(w(lm2(5,3),n_r1)), &
                           REAL(w(lm2(5,4),n_r1)), &
                           REAL(w(lm2(5,5),n_r1)), &
                           REAL(w(lm2(6,1),n_r1)), &
                           REAL(w(lm2(6,2),n_r1)), &
                           REAL(w(lm2(6,3),n_r1)), &
                           REAL(w(lm2(6,4),n_r1)), &
                           REAL(w(lm2(6,5),n_r1)), &
                           REAL(w(lm2(6,6),n_r1))
        CLOSE(n_SRIC_file)

        filename='inerT.'//tag
        OPEN(n_SRIC_file,FILE=filename,STATUS='UNKNOWN', &
             POSITION='APPEND')
        WRITE(n_SRIC_file,'(1P,2X,D20.12,21D12.4)') &
                                              time, &
                           REAL(z(lm2(1,1),n_r1)), &
                           REAL(z(lm2(2,1),n_r1)), &
                           REAL(z(lm2(2,2),n_r1)), &
                           REAL(z(lm2(3,1),n_r1)), &
                           REAL(z(lm2(3,2),n_r1)), &
                           REAL(z(lm2(3,3),n_r1)), &
                           REAL(z(lm2(4,1),n_r1)), &
                           REAL(z(lm2(4,2),n_r1)), &
                           REAL(z(lm2(4,3),n_r1)), &
                           REAL(z(lm2(4,4),n_r1)), &
                           REAL(z(lm2(5,1),n_r1)), &
                           REAL(z(lm2(5,2),n_r1)), &
                           REAL(z(lm2(5,3),n_r1)), &
                           REAL(z(lm2(5,4),n_r1)), &
                           REAL(z(lm2(5,5),n_r1)), &
                           REAL(z(lm2(6,1),n_r1)), &
                           REAL(z(lm2(6,2),n_r1)), &
                           REAL(z(lm2(6,3),n_r1)), &
                           REAL(z(lm2(6,4),n_r1)), &
                           REAL(z(lm2(6,5),n_r1)), &
                           REAL(z(lm2(6,6),n_r1))
        CLOSE(n_SRIC_file)

    END IF


    RETURN
    end SUBROUTINE write_rot
!-----------------------------------------------------------------------
