!**************************************************************
    SUBROUTINE logWrite(message)
!**************************************************************
  USE logic,ONLY: l_save_out
  USE parallel_mod,ONLY: rank
  USE output_data,ONLY: n_log_file,log_file
  IMPLICIT NONE

!--- Input:
    CHARACTER(len=*) :: message

!--------------------------------------------------------------
    IF ( rank.eq.0 ) THEN
        IF ( l_save_out ) THEN
            OPEN(n_log_file,file=log_file,status='UNKNOWN', &
                 POSITION='APPEND')
        END IF
        !WRITE(n_log_file,*)
        WRITE(n_log_file,*) trim(message)
        !WRITE(*,*)
        WRITE(*,*)          TRIM(message)
        IF ( l_save_out ) CLOSE(n_log_file)
    END IF

    RETURN
    end SUBROUTINE logWrite

!--------------------------------------------------------------

