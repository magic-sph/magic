!$Id$
!***********************************************************************
    SUBROUTINE checkMovie
!***********************************************************************

!    !------------ This is release 2 level 10  --------------!
!    !------------ Created on 2/5/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Checks whether movie is desired and opens the respective files.  |
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

!-- input/output
! include 'truncation.f'
! include 'c_logic.f'
! include 'c_output.f'
! include 'c_parallel.f'

!-- Local
    INTEGER :: n

!----------------------------------------------------------------------


    IF ( .NOT. l_movie ) THEN

    ! MPI: this should be known to all procs !
        l_movie_oc=.false.
        l_movie_ic=.false.
        l_HTmovie =.false.
        l_dtBmovie=.false.
        l_store_frame=.false.

    ELSE

    ! MPI: this should be known to all procs !
        CALL get_movie_type

        IF ( iAmProc == movProc ) THEN
        !----- Open movie files on first processor only:
            IF ( .NOT. l_save_out ) THEN
                DO n=1,n_movies
                    OPEN(n_movie_file(n),file=movie_file(n), &
                         STATUS='NEW',FORM='UNFORMATTED')
                END DO
            END IF
        END IF
        IF ( l_dtBmovie .AND. ldtBMem == 0 ) then
            CALL logWrite('! You required dtB caculation.')
            CALL logWrite('! Please set ldtBMem=1 in truncation.f')
            CALL logWrite('! This is needed to reserve memory.')
            STOP
        END IF

    END IF

    RETURN
    end SUBROUTINE checkMovie

!**********************************************************************
