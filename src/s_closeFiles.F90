!$Id$
!***********************************************************************
    SUBROUTINE closeFiles
!***********************************************************************

!    !------------ This is release 2 level 10  --------------!
!    !------------ Created on 2/5/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  | Defines names and unit for output files and opens then.           |
!  | MPI: called only by the processor responsible for output !        |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE logic
    USE output_data
    USE parallel_mod

    IMPLICIT NONE

!-- Input: read inputfiles first!
! include 'truncation.f'
! include 'c_logic.f'
! include 'c_output.f'
! include 'c_parallel.f'

!-- Local:
    INTEGER :: n

!-- end of declaration
!-----------------------------------------------------------------------


!-- Close various output files:
    IF ( .NOT. l_save_out ) THEN

        IF ( iAmProc == logProc ) THEN
            CLOSE(n_log_file)
            CLOSE(n_e_kin_file)
            CLOSE(n_par_file)
            IF ( l_AM ) THEN
                CLOSE(n_angular_file)
            END IF
            IF ( l_anel ) THEN
                CLOSE(n_u_square_file)
            END IF
            IF ( l_RMS .OR. l_RMStest ) THEN
                CLOSE(n_dtvrms_file)
                CLOSE(n_dtvasrms_file)
            END IF
            IF ( l_mag ) then
                CLOSE(n_e_mag_oc_file)
                CLOSE(n_e_mag_ic_file)
                CLOSE(n_dipole_file)
                IF ( l_RMS .OR. l_RMStest) THEN
                    CLOSE(n_dtbrms_file)
                    CLOSE(n_dtdrms_file)
                END IF
                IF ( l_cmb_field ) close(n_cmb_file)
                IF ( l_cmb_field .AND. l_movie ) close(n_cmbMov_file)
                IF ( l_r_field ) THEN
                    DO n=1,n_coeff_r_max
                        CLOSE(n_v_r_file(n))
                        CLOSE(n_b_r_file(n))
                    END DO
                ENDIF
                IF ( l_dt_cmb_field ) CLOSE(n_dt_cmb_file)
            END IF
            IF ( l_rot_ic .OR. l_rot_ma .AND.      &
                 .NOT. l_SRIC .AND. .NOT. l_SRMA ) &
                CLOSE(n_rot_file)
            CLOSE(n_misc_file)
            CLOSE(n_power_file)
        END IF
        IF ( iAmProc == movProc .AND. l_movie ) THEN
            IF ( l_movie ) THEN
                DO n=1,n_movies
                    CLOSE(n_movie_file(n))
                END DO
            END IF
        END IF
    END IF


    RETURN
    end SUBROUTINE closeFiles

!-------------------------------------------------------------------
