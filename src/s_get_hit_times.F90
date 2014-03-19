!$Id$
!***********************************************************************
    subroutine get_hit_times(t,n_t_max,n_t,l_t,t_start,t_stop,dt, &
                             n_tot,n_step,string,time,tScale)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  This subroutine checks whether any specific times t(*) are given |
!  |  on input. If so, it returns their number n_r and sets l_t        |
!  |  to true. If not, t(*) may also be defined by giving a time step  |
!  |  dt or a number n_tot of desired output times and t_stop>t_start. |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    implicit none

!-- Input/Output
    integer :: n_t_max    ! Dimension of t(*)
    real(kind=8) ::  t(n_t_max) ! Times for output
    real(kind=8) ::  t_start    ! Starting time for output
    real(kind=8) ::  t_stop     ! Stop time for output
    real(kind=8) ::  dt         ! Time step for output
    real(kind=8) ::  time       ! Time of start file
    real(kind=8) ::  tScale
    integer :: n_t        ! No. of output times
    logical :: l_t        ! =.true. if output times are defined
    integer :: n_tot      ! No. of output (times) if no times defined
    integer :: n_step     ! Ouput step in no. of time steps
    character(len=*) :: string

!-- Local
    integer :: n         ! Counter

!-- end of declaration
!-----------------------------------------------------------------------

    t_start=t_start/tScale
    t_stop =t_stop/tScale
    dt     =dt/tScale

!-- Check whether any time is given explicitely:
    l_t=.false.
    n_t=0
    do n=1,n_t_max
        if ( t(n) >= 0.d0 ) then
            t(n)=t(n)/tScale
            l_t=.true.
            n_t=n_t+1
        end if
    end do

!-- Check times should be constructed:
    if ( t_start < time ) t_start=time
    if ( .NOT. l_t .AND. ( dt > 0.d0 .OR. &
       ( n_tot > 0 .AND. t_stop > t_start ) ) ) then

        if ( n_tot > 0 .AND. dt > 0.d0 ) then
            n_t  =n_tot
            n_tot=0
        else if ( dt > 0.d0 ) then
            if ( t_stop > t_start ) then
                n_t=IDINT((t_stop-t_start)/dt)+1
            else
                n_t=n_t_max
            end if
        else if ( n_tot > 0 ) then
            n_t=n_tot
            n_tot=0
            dt=(t_stop-t_start)/dble(n_t-1)
        end if
        if ( n_t > n_t_max ) then
            write(*,*) '! Sorry, maximum no. of times for'
            write(*,*) '! output ',string
            write(*,*) '! is:',n_t_max
            write(*,*) '! Increase n_time_hits in c_output.f!'
            stop
        end if

        l_t=.true.
        if ( t_start == time ) then
            n_t=n_t-1
            t(1)=t_start+dt
        else
            t(1)=t_start
        end if
        do n=2,n_t
            t(n)=t(n-1)+dt
        end do

    end if


    if ( n_tot /= 0 .AND. n_step /= 0 ) then
        write(*,*)
        write(*,*) '! You have to either provide the total'
        write(*,*) '! number or the step for output:'
        WRITE(*,'(A,2(A,I10))') string, "n_tot = ",n_tot,", n_step = ",n_step
        write(*,*) '! I set the step width to zero!'
        n_step=0
    end if

    if ( l_t ) then
        t_start=t(1)
        t_stop =t(n_t)
        dt     =t(2)-t(1)
    end if
            
    return
    end subroutine get_hit_times
!------------------------------------------------------------------------------
