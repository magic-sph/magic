!$Id$
!***********************************************************************
    SUBROUTINE openFiles
!***********************************************************************

!    !------------ This is release 2 level 10  --------------!
!    !------------ Created on 2/5/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  | Defines names and unit for output files and opens then.           |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE init_fields
    USE logic
    USE output_data
    USE parallel_mod
    USE charmanip, only: length_to_blank

    IMPLICIT NONE

!-- Output written into c_output.f

!-- Local:
    CHARACTER(len=72) :: string
    INTEGER :: n
    INTEGER ::  length

!-- end of declaration
!-----------------------------------------------------------------------


!-- Define output unit numbers:
    n_log_file         =7
    n_start_file       =8
    n_rst_file         =9
    n_graph_file       =10
    n_e_kin_file       =11
    n_kin_spec_file    =12
    n_e_mag_ic_file    =13
    n_e_mag_oc_file    =14
    n_mag_spec_file    =15
    n_u2_spec_file     =121
    n_lp_file          =16
    n_rot_file         =17
    n_dipole_file      =18
    n_signal_file      =19
    n_cmb_file         =20
    n_misc_file        =21
    n_cmbMov_file      =24
    n_SRIC_file        =25
    n_SRMA_file        =26
    n_dt_cmb_file      =27
    n_power_file       =28
    n_u_square_file    =29
    n_par_file         =300
    n_angular_file     =301
    n_dtvrms_file      =302
    n_dtvasrms_file    =303
    n_dtbrms_file      =304
    n_dtdrms_file      =305
    DO n=1,n_coeff_r_max
        n_v_r_file(n)=30+n
    END DO
    DO n=1,n_coeff_r_max
        n_v_r_mov_file(n)=40+n
    END DO
    DO n=1,n_coeff_r_max
        n_b_r_file(n)=50+n
    END DO
    DO n=1,n_coeff_r_max
        n_b_r_mov_file(n)=60+n
    END DO
    DO n=1,n_movies_max
        n_movie_file(n)=70+n
    END DO
    nLF=n_log_file

!-- Define output file names:
    e_kin_file='e_kin.'//tag
    log_file='log.'//tag
    par_file='par.'//tag
    IF ( l_mag ) THEN
        e_mag_ic_file='e_mag_ic.'//tag
        e_mag_oc_file='e_mag_oc.'//tag
        dipole_file='dipole.'//tag
        IF ( l_RMS .OR. l_RMStest) THEN
            dtbrms_file='dtBrms.'//tag
            dtdrms_file='dtDrms.'//tag
        END IF
    END IF
    IF ( l_AM ) THEN
        angular_file='AM.'//tag
    END IF
    IF ( l_anel ) THEN
        u_square_file='u_square.'//tag
    END IF
    IF ( l_RMS .OR. l_RMStest) THEN
        dtvrms_file='dtVrms.'//tag
        dtvasrms_file='dtVAsRms.'//tag
    END IF
    IF ( l_rot_ic .OR. l_rot_ma ) THEN
        rot_file='rot.'//tag
    END IF
    IF ( l_cmb_field ) THEN
        cmb_file   ='B_coeff_cmb.'//tag
        cmbMov_file='B_coeff_cmbMov.'//tag
    END IF
    IF ( l_dt_cmb_field ) THEN
        dt_cmb_file   ='B_coeff_dt_cmb.'//tag
    END IF
    IF ( l_r_field ) THEN
        DO n=1,n_coeff_r_max
            WRITE(string,'(''V_coeff_r'',i1,''.'')') n
            length=length_to_blank(string)
            v_r_file(n)=string(1:length)//tag
            WRITE(string,'(''B_coeff_r'',i1,''.'')') n
            length=length_to_blank(string)
            B_r_file(n)=string(1:length)//tag
        END DO
    END IF
    misc_file ='misc.'//tag
    SRIC_file ='SRIC.'//tag
    SRMA_file ='SRMA.'//tag
    power_file='power.'//tag

!-- Open various output files that will be used throughout the run:
    IF ( .NOT. l_save_out ) THEN

        IF ( iAmProc == logProc ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN')
            OPEN(n_e_kin_file,FILE=e_kin_file,STATUS='NEW')
            OPEN(n_par_file,FILE=par_file,STATUS='NEW')
            IF ( l_RMS .OR. l_RMStest) THEN
                OPEN(n_dtvrms_file,FILE=dtvrms_file,STATUS='NEW')
                OPEN(n_dtvasrms_file,FILE=dtvasrms_file,STATUS='NEW')
            END IF
            IF ( l_anel ) THEN
                OPEN(n_u_square_file,FILE=u_square_file,STATUS='NEW')
            END IF
            IF ( l_AM ) THEN
                OPEN(n_angular_file,FILE=angular_file,STATUS='NEW')
            END IF
            IF ( l_mag ) THEN
                OPEN(n_e_mag_oc_file,FILE=e_mag_oc_file,STATUS='NEW')
                OPEN(n_e_mag_ic_file,FILE=e_mag_ic_file,STATUS='NEW')
                OPEN(n_dipole_file,file=dipole_file,status='new')
                IF ( l_RMS .OR. l_RMStest) THEN
                    OPEN(n_dtbrms_file,FILE=dtbrms_file,STATUS='NEW')
                    OPEN(n_dtdrms_file,FILE=dtdrms_file,STATUS='NEW')
                END IF
                IF ( l_cmb_field ) THEN
                    OPEN(n_cmb_file,file=cmb_file, &
                         STATUS='NEW',FORM='UNFORMATTED')
                    IF ( l_movie ) THEN
                        OPEN(n_cmbMov_file,FILE=cmbMov_file, &
                             STATUS='NEW',FORM='UNFORMATTED')
                    END IF
                END IF
                IF ( l_dt_cmb_field )                    &
                    OPEN(n_dt_cmb_file,FILE=dt_cmb_file, &
                         STATUS='NEW',FORM='UNFORMATTED')
                IF ( l_r_field ) THEN
                    DO n=1,n_coeff_r_max
                        OPEN(n_v_r_file(n),FILE=v_r_file(n), &
                             STATUS='NEW',FORM='UNFORMATTED')
                        OPEN(n_b_r_file(n),FILE=b_r_file(n), &
                             STATUS='NEW',FORM='UNFORMATTED')
                    END DO
                ENDIF
            END IF
            IF ( .NOT. l_SRIC .AND. .NOT. l_SRMA ) THEN
                IF( l_rot_ic .OR. l_rot_ma ) &
                     OPEN(n_rot_file,FILE=rot_file,STATUS="NEW")
            END IF
            OPEN(n_misc_file,FILE=misc_file,STATUS="NEW")
            IF ( l_power ) &
                 OPEN(n_power_file,FILE=power_file,STATUS="NEW")
        END IF

    !          IF ( iAmProc.EQ.TOProc .AND. l_TO ) THEN
    !--------  TO files are openen in outTO !
    !          END IF

    !          IF ( iAmProc.EQ.movProc .AND. l_movie ) THEN
    !-------- open movie files ?
    !          END IF

    END IF


    RETURN
    end SUBROUTINE openFiles

!-------------------------------------------------------------------
