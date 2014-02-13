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

  USE truncation
  USE logic
  USE output_data
  USE parallel_mod,only: rank
  USE movie_data,ONLY: n_movies,n_movie_file
  IMPLICIT NONE

  !-- Local
  INTEGER :: n

  !----------------------------------------------------------------------


  IF ( .NOT. l_movie ) THEN

     l_movie_oc=.false.
     l_movie_ic=.false.
     l_HTmovie =.false.
     l_dtBmovie=.false.
     l_store_frame=.false.

  ELSE

     CALL get_movie_type

     IF ( rank.eq.0 ) THEN
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
