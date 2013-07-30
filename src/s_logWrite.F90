!$Id$
!**************************************************************
    SUBROUTINE logWrite(message)
!**************************************************************

    USE truncation
    USE logic
    USE output_data
    USE parallel_mod

    IMPLICIT NONE

!--- Input:
    CHARACTER(len=*) :: message

!--- Local:

!--------------------------------------------------------------

    IF ( iAmProc == logProc ) THEN
        IF ( l_save_out ) THEN
            OPEN(n_log_file,file=log_file,status='UNKNOWN', &
                 POSITION='APPEND')
        END IF
        WRITE(n_log_file,*)
        WRITE(n_log_file,*) message
        WRITE(*,*)
        WRITE(*,*)          message
        IF ( l_save_out ) CLOSE(n_log_file)
    END IF

    RETURN
    end SUBROUTINE logWrite

!--------------------------------------------------------------
