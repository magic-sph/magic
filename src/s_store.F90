!$Id$
!********************************************************************
    SUBROUTINE store(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic)
!********************************************************************

!--------------------------------------------------------------------

! *** store results on disc file (restart file)
!   In addition to the magnetic field and velocity potentials
!   we store the time derivative terms
!     djdt(lm,nR),dbdt(lm,nR), ......

!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE init_fields
    USE blocking
    USE logic
    USE fieldsLast
    USE output_data

    IMPLICIT NONE

!-- Input of variables:
    REAL(kind=8),INTENT(IN) :: time,dt,dtNew

!-- Input of scalar fields to be stored:
    COMPLEX(kind=8),INTENT(IN) :: w(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: z(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: p(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: s(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: b(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: aj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: b_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: aj_ic(lm_maxMag,n_r_ic_maxMag)

!-- end of declaration
!---------------------------------------------------------------------

!-- Write parameters:
    IF ( .NOT. l_heat ) THEN
        inform=11
    ELSE
        inform=12
    END IF

    WRITE(n_rst_file)          time*tScale,dt*tScale, &
                             ra,pr,prmag,ek,radratio, &
    inform,n_r_max,n_theta_max,n_phi_tot,minc,nalias, &
                               n_r_ic_max,sigma_ratio

!-- Write velocity, pressure, entropy:
    !WRITE(*,"(A,I4,A,2ES20.13)") "w(3,",1,") = ",w(3,1)

    IF ( l_heat ) THEN
        WRITE(n_rst_file) w,z,p,s
        WRITE(n_rst_file) dsdtLast,dwdtLast,dzdtLast,dpdtLast
    ELSE
        WRITE(n_rst_file) w,z,p
        WRITE(n_rst_file) dwdtLast,dzdtLast,dpdtLast
    END IF

!-- Write magnetic field:
    IF ( l_mag ) WRITE(n_rst_file) b,aj,dbdtLast,djdtLast

!-- Write IC magnetic field:
    IF ( l_mag .AND. l_cond_ic ) &
        WRITE(n_rst_file) b_ic,aj_ic,dbdt_icLast,djdt_icLast

!-- Store Lorentz-torques and rotation rates:
    WRITE(n_rst_file) lorentz_torque_icLast, &
                      lorentz_torque_maLast, &
          omega_ic1,omegaOsz_ic1,tOmega_ic1, &
          omega_ic2,omegaOsz_ic2,tOmega_ic2, &
          omega_ma1,omegaOsz_ma1,tOmega_ma1, &
          omega_ma2,omegaOsz_ma2,tOmega_ma2, &
                                       dtNew


    RETURN
    end SUBROUTINE store

!----------------------------------------------------------------------
