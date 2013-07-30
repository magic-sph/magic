!$Id$
!***********************************************************************
    SUBROUTINE preCalcTimes(time,n_time_step)
!***********************************************************************

!    !------------ This is release 2 level 10  --------------!
!    !------------ Created on 2/5/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Precalc. after time, time and dthas been read from startfile.    |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE num_param
    USE init_fields
    USE logic
    USE output_data

    IMPLICIT NONE

    REAL(kind=8) ::  time
    INTEGER :: n_time_step

!-- Local:
    LOGICAL :: l_time
    INTEGER :: n

!-----------------------------------------------------------------------


!----- Set time step:
    IF ( l_reset_t ) THEN
        time=0.D0
        n_time_step=0
    END IF

    tmagcon=tmagcon+time

!-- Get output times:
    l_time_hits=.false.
    CALL get_hit_times(t_graph,n_time_hits,n_t_graph,l_time, &
                        t_graph_start,t_graph_stop,dt_graph, &
                  n_graphs,n_graph_step,'graph',time,tScale)
    l_time_hits=l_time_hits.OR.l_time
    CALL get_hit_times(t_rst,n_time_hits,n_t_rst,l_time, &
                          t_rst_start,t_rst_stop,dt_rst, &
                    n_rsts,n_rst_step,'rst',time,tScale)
    l_time_hits=l_time_hits.OR.l_time
    CALL get_hit_times(t_log,n_time_hits,n_t_log,l_time, &
                          t_log_start,t_log_stop,dt_log, &
                    n_logs,n_log_step,'log',time,tScale)
    l_time_hits=l_time_hits.OR.l_time
    CALL get_hit_times(t_spec,n_time_hits,n_t_spec,l_time, &
                         t_spec_start,t_spec_stop,dt_spec, &
                   n_specs,n_spec_step,'spec',time,tScale)
    l_time_hits=l_time_hits.OR.l_time
    IF ( l_cmb_field ) then
        l_cmb_field=.false.
        CALL get_hit_times(t_cmb,n_time_hits,n_t_cmb,l_time, &
                              t_cmb_start,t_cmb_stop,dt_cmb, &
                        n_cmbs,n_cmb_step,'cmb',time,tScale)
        if ( n_cmbs > 0 .OR. n_cmb_step > 0 .OR. &
             l_time ) l_cmb_field= .TRUE. 
        l_time_hits=l_time_hits.OR.l_time
    END IF
    l_dt_cmb_field=l_dt_cmb_field .AND. l_cmb_field
    IF ( l_movie ) then
        CALL get_hit_times(t_movie,n_time_hits,n_t_movie,l_time, &
                            t_movie_start,t_movie_stop,dt_movie, &
                n_movie_frames,n_movie_step,'movie',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
    END IF
    IF ( l_TO ) then
        IF ( n_TOs == 0 .AND. n_t_TO == 0 ) &
             n_TO_step=MAX(3,n_TO_step)
        CALL get_hit_times(t_TO,n_time_hits,n_t_TO,l_time, &
                               t_TO_start,t_TO_stop,dt_TO, &
                         n_TOs,n_TO_step,'TO',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
        IF ( n_TOZs == 0 .AND. n_t_TOZ == 0 ) n_TOZs=3
        CALL get_hit_times(t_TOZ,n_time_hits,n_t_TOZ,l_time, &
                              t_TOZ_start,t_TOZ_stop,dt_TOZ, &
                        n_TOZs,n_TOZ_step,'TOZ',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
    END IF
    IF ( l_TOmovie ) then
        CALL get_hit_times(t_TOmovie,n_time_hits,n_t_TOmovie,l_time, &
                          t_TOmovie_start,t_TOmovie_stop,dt_TOmovie, &
              n_TOmovie_frames,n_TOmovie_step,'TOmovie',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
    END IF

    IF ( l_storePot ) THEN
        l_storeVpot   =.TRUE.
        n_Vpot_step   =n_pot_step
        n_Vpots       =n_pots
        t_Vpot_start  =t_pot_start
        t_Vpot_stop   =t_pot_stop
        dt_Vpot       =dt_pot
        l_storeBpot   =.TRUE.
        n_Bpot_step   =n_pot_step
        n_Bpots       =n_pots
        t_Bpot_start  =t_pot_start
        t_Bpot_stop   =t_pot_stop
        dt_Bpot       =dt_pot
        l_storeTpot   =.TRUE.
        n_Tpot_step   =n_pot_step
        n_Tpots       =n_pots
        t_Tpot_start  =t_pot_start
        t_Tpot_stop   =t_pot_stop
        dt_Tpot       =dt_pot
        DO n=1,n_time_hits
            t_Bpot(n)=t_pot(n)
            t_Vpot(n)=t_pot(n)
            t_Tpot(n)=t_pot(n)
        END DO
    END IF

    IF ( l_storeBpot ) then
        CALL get_hit_times(t_Bpot,n_time_hits,n_t_Bpot,l_time, &
                             t_Bpot_start,t_Bpot_stop,dt_Bpot, &
                       n_Bpots,n_Bpot_step,'Bpot',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
    END IF
    IF ( l_storeVpot ) then
        CALL get_hit_times(t_Vpot,n_time_hits,n_t_Vpot,l_time, &
                             t_Vpot_start,t_Vpot_stop,dt_Vpot, &
                       n_Vpots,n_Vpot_step,'Vpot',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
    END IF
    IF ( l_storeTpot ) then
        CALL get_hit_times(t_Tpot,n_time_hits,n_t_Tpot,l_time, &
                             t_Tpot_start,t_Tpot_stop,dt_Tpot, &
                       n_Tpots,n_Tpot_step,'Tpot',time,tScale)
        l_time_hits=l_time_hits.OR.l_time
    END IF


    RETURN
    end SUBROUTINE preCalcTimes

!-------------------------------------------------------------------------------
