!$Id$
!**************************************************************
    SUBROUTINE safeOpen(nf,name)
!**************************************************************

    USE logic

    IMPLICIT NONE

    INTEGER :: nf
    CHARACTER(len=*) :: name

!--------------------------------------------------------------

    IF ( l_save_out ) THEN
        OPEN(nf,file=name,status='unknown',POSITION='APPEND')
    END IF

    RETURN
    end SUBROUTINE safeOpen

!--------------------------------------------------------------

!**************************************************************
    SUBROUTINE safeClose(nf)
!**************************************************************

    USE logic

    IMPLICIT NONE

    INTEGER :: nf

!--------------------------------------------------------------

    IF ( l_save_out ) THEN
        CLOSE(nf)
    END IF

    RETURN
    end SUBROUTINE safeClose

!--------------------------------------------------------------

