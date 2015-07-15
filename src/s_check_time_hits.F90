!$Id$
!***********************************************************************
    SUBROUTINE check_time_hits(l_new_dt,time,dt,dt_new)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Checks whether a certain dt is required to hit a                 |
!  |  specific output-time.                                            |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE num_param
    USE logic
    USE output_data
    USE parallel_mod,only:rank
    IMPLICIT NONE

!-- input:
! include 'truncation.f'
! include 'c_output.f'
! include 'c_num_param.f'
! include 'c_logic.f'

!-- output: ev. modified dt
    LOGICAL :: l_new_dt ! signfies change of dt !
    REAL(kind=8) :: time,dt,dt_new
     
!-- local:
    INTEGER :: n_dt_hit
    INTEGER, PARAMETER :: n_dt_hit_max=10
    REAL(kind=8) ::  dt_hit(n_dt_hit_max) ! dt for different hit times
    INTEGER :: n                    ! counter
    REAL(kind=8) ::  time_new             ! Next time step


!-- end of declaration
!----------------------------------------------------------------------

    time_new=time+dt
    l_new_dt=.FALSE.

    n_dt_hit=7

    DO n=1,n_dt_hit
        dt_hit(n)=0.D0
    END DO

    DO n=1,n_time_hits
        IF ( t_rst(n) > time .AND. &
             t_rst(n) < time_new ) &
             dt_hit(1)=t_rst(n)-time
        IF ( t_graph(n) > time .AND. &
             t_graph(n) < time_new ) &
             dt_hit(2)=t_graph(n)-time
        IF ( t_log(n) > time .AND. &
             t_log(n) < time_new ) &
             dt_hit(3)=t_log(n)-time
        IF ( t_spec(n) > time .AND. &
             t_spec(n) < time_new ) &
             dt_hit(4)=t_spec(n)-time
        IF ( t_cmb(n) > time .AND. &
             t_cmb(n) < time_new ) &
             dt_hit(5)=t_cmb(n)-time
        IF ( t_movie(n) > time .AND. &
             t_movie(n) < time_new ) &
             dt_hit(6)=t_movie(n)-time
        IF ( t_TO(n) > time .AND. &
             t_TO(n) < time_new ) &
             dt_hit(7)=t_TO(n)-time
        IF ( t_TOmovie(n) > time .AND. &
             t_TOmovie(n) < time_new ) &
             dt_hit(7)=t_TOmovie(n)-time
    END DO

    DO n=1,n_dt_hit
        IF ( dt_hit(n) /= 0.D0 .AND. &
             dt_hit(n) < dt_new ) THEN
            l_new_dt=.TRUE.
            dt_new=dt_hit(n)
        END IF
    END DO

    IF ( l_new_dt ) THEN
        IF ( dt_new < dtMin ) dt_new=dtMin
        time_new=time+dt_new
        WRITE(*, &
      &     '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2d16.6)') &
      &     time_new*tScale,time*tScale
        IF (rank == 0) THEN
           IF ( l_save_out ) THEN
              OPEN(n_log_file,file=log_file,status='unknown', &
                   POSITION='APPEND')
              WRITE(n_log_file, &
                   &     '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2d16.6)') &
                   &     time_new*tScale,time*tScale
              CLOSE(n_log_file)
           ELSE
              WRITE(n_log_file, &
                   &    '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2d16.6)') &
                   &    time_new*tScale,time*tScale
           END IF
        END IF
    END IF

    RETURN
    end SUBROUTINE check_time_hits
              
!----------------------------------------------------------------------
