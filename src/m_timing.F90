!$Id$
module timing
  
  IMPLICIT NONE

contains
  !******************************************************************
  SUBROUTINE wallTime(time)
    !******************************************************************
    !   This routine returns the wallclock time in four
    !   integer arguments.
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: time
    INTEGER :: count,countRate,countMax
    INTEGER(kind=8) :: mSeconds

    !------------------------------------------------------------------

    !--- SYSTEM_CLOCK is a Fortran90 subroutine that
    !    returns the wallclock time with a precision
    !    determined by countRate (counts per second).
    !    For example, countRate is set to 100 on the
    !    SP2 under AIX so that the precision is
    !    10 milli seconds. The count is reset when
    !    countMax is reached. On the SP2 under AIX
    !    countMax is set to 24 hours, and the reset
    !    is performed at midnight. Thus, the count
    !    can be used to get the real wallclock time.

    CALL SYSTEM_CLOCK(count,countRate,countMax)

    mSeconds=IDINT(1.D3*DBLE(count)/DBLE(countRate))
    CALL ms2time(mSeconds,time)

  END SUBROUTINE wallTime

  !------------------------------------------------------------------
  SUBROUTINE get_resetTime(resetTime)
    !******************************************************************
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: resetTime
    INTEGER :: count,countRate,countMax
    INTEGER(kind=8) :: mSeconds

    !------------------------------------------------------------------

    CALL SYSTEM_CLOCK(count,countRate,countMax)

    mSeconds=IDINT(1.D3*DBLE(countMax)/DBLE(countRate))
    CALL ms2time(mSeconds,resetTime)

  END SUBROUTINE get_resetTime

  !------------------------------------------------------------------

  SUBROUTINE ms2time(ms,time)
    !******************************************************************
    !  Transforms accumulated milliseconds ms into an four-element
    !  integer arrays time(4) containing the time in
    !  hours=time(1), minutes=time(2), seconds=time(3),
    !  and milliseconds=time(4).
    !------------------------------------------------------------------

    INTEGER(kind=8),parameter :: msecSecond=1000
    INTEGER(kind=8),parameter :: msecMinute=60000
    INTEGER(kind=8),parameter :: msecHour=3600000

    INTEGER(kind=8) :: ms,mSeconds
    INTEGER, dimension(4) :: time
    INTEGER :: seconds,minutes,hours

    !------------------------------------------------------------------

    mSeconds=mS
    hours   =INT(mSeconds/msecHour,4)
    mSeconds=mSeconds-hours*msecHour
    minutes =INT(mSeconds/msecMinute,4)
    mSeconds=mSeconds-minutes*msecMinute
    seconds =INT(mSeconds/msecSecond,4)
    mSeconds=mSeconds-seconds*msecSecond
    time(1) =hours
    time(2) =minutes
    time(3) =seconds
    time(4) =INT(mSeconds,4)

  END SUBROUTINE ms2time

  !------------------------------------------------------------------
  INTEGER(kind=8) FUNCTION time2ms(time)
    !******************************************************************
    !  Transforms a four-element integer arrays time(4) containing the
    !  time in hours=time(1), minutes=time(2), seconds=time(3),
    !  and milliseconds=time(4) into accumulated milliseconds.
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: time
    INTEGER(kind=8) :: msecMinute
    INTEGER(kind=8) :: msecHour
    INTEGER(kind=8) :: msecSecond

    !------------------------------------------------------------------

    msecHour  =3600000
    msecMinute=60000
    msecSecond=1000

    time2ms=time(1)*msecHour+time(2)*msecMinute+time(3)*msecSecond+time(4)

  END FUNCTION time2ms

  !------------------------------------------------------------------
  SUBROUTINE subTime(timeStart,timeStop,timeD)
    !  Returns time passed between timeStop and timeStart.
    !  Note timeStop has to be younger than timeStart, otherwise
    !  24 hours are added. This is necessary on systems like the IBM
    !  where the time counter as reset every day at midnight.
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: timeD,timeStart,timeStop
    INTEGER(kind=8) :: msDiff,msStart,msStop

    INTEGER :: count,countRate,countMax
    INTEGER(kind=8) :: msReset

    INTEGER(kind=8), PARAMETER :: msDay=1000*24*3600

    !------------------------------------------------------------------

    CALL SYSTEM_CLOCK(count,countRate,countMax)
    msReset=IDINT(1.D3*DBLE(countMax)/DBLE(countRate))
    msStart=time2ms(timeStart)
    msStop =time2ms(timeStop)

    IF ( msStop < msStart ) THEN
       IF ( msStop > msStart-1 ) THEN
          ! Some weird inprecision:
          msStop=msStop+1
       ELSE
          IF ( msDay < msReset ) THEN
             IF ( msStop > msStart-msDay ) THEN
                ! Thats a day reset:
                msStop=msStop+msDay
             ELSE IF ( msStop > msStart-msReset ) THEN
                ! Thats the overrun reset:
                msStop=msStop+msReset
             ELSE
                WRITE(*,*) 'Timing problem!'
                WRITE(*,*) 'Negative time in ms:',msStop-msStart
                STOP
             END IF
          ELSE IF ( msStop > msStart-msReset ) THEN
             ! Thats the overrun reset:
             msStop=msStop+msReset
          ELSE
             WRITE(*,*) 'Timing problem!'
             WRITE(*,*) 'Negative time in ms:',msStop-msStart
             STOP
          END IF
       END IF
    END IF
    msDiff=msStop-msStart
    CALL ms2time(msDiff,timeD)

  END SUBROUTINE subTime

  !------------------------------------------------------------------
  LOGICAL FUNCTION lTimeLimit(time,timeMax)
    !  False when time exeeds timeMax
    !------------------------------------------------------------------

    integer, dimension(4) :: time,timeMax

    !------------------------------------------------------------------

    lTimeLimit=.FALSE.
    IF ( time(1) > timeMax(1) ) THEN
       lTimeLimit=.TRUE.
    ELSE IF ( time(1) == timeMax(1) ) THEN
       IF ( time(2) > timeMax(2) ) THEN
          lTimeLimit=.TRUE.
       ELSE IF ( time(2) == timeMax(2) ) THEN
          IF ( time(3) > timeMax(3) ) THEN
             lTimeLimit=.TRUE.
          ELSE IF ( time(3) == timeMax(3) ) THEN
             IF ( time(4) >= timeMax(4) ) &
                  lTimeLimit= .TRUE. 
          END IF
       END IF
    END IF

  end FUNCTION lTimeLimit

  !------------------------------------------------------------------
  SUBROUTINE addTime(time1,time2)
    !******************************************************************
    !  Adds time2 to time1
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: time1,time2
    INTEGER(kind=8) :: t1ms,t2ms

    !------------------------------------------------------------------

    t1ms=time2ms(time1)
    t2ms=time2ms(time2)
    t1ms=t1ms+t2ms
    CALL ms2time(t1ms,time1)

  end SUBROUTINE addTime

  !------------------------------------------------------------------

  SUBROUTINE meanTime(time,n)
    !  Gets mean time
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: time
    INTEGER :: n
    INTEGER(kind=8) :: tms

    !------------------------------------------------------------------

    tms=time2ms(time)
    IF ( n==0 ) THEN
       WRITE(*,*) '! Error in meanTime !'
       WRITE(*,*) '! Cant build mean for zero steps!'
       tms=0
    ELSE
       tms=tms/n
    END IF
    CALL ms2time(tms,time)

    RETURN
  end SUBROUTINE meanTime

  !------------------------------------------------------------------

  LOGICAL FUNCTION lNegTime(time1,time2)
    !  Negative passed time? Means we have passed midnight.
    !  The wallclock time is reset to zero on some
    !  computers at midnight.
    !------------------------------------------------------------------

    INTEGER, dimension(4) :: time1,time2
    INTEGER(kind=8) :: tms1,tms2

    INTEGER :: count,countRate,countMax

    !------------------------------------------------------------------

    CALL SYSTEM_CLOCK(count,countRate,countMax)

    tms1=time2ms(time1)
    tms2=time2ms(time2)
    IF ( countRate*tms2 < countRate*tms1-1 ) THEN
       lNegTime=.TRUE.
    ELSE
       lNegTime=.FALSE.
    END IF

  end FUNCTION lNegTime
  !------------------------------------------------------------------
  SUBROUTINE writeTime(nOut,text,time)
    !  Returns time passed between timeStop and timeStart.
    !  Note timeStop has to be younger than timeStart, otherwise
    !  24 hours are added. This is necessary on systems like the IBM
    !  where the time counter as reset every day at midnight.
    !------------------------------------------------------------------

    !--- Input:
    INTEGER :: nOut
    CHARACTER(len=*) :: text
    INTEGER, dimension(4) :: time
    INTEGER :: timeH

    !--- Local
    INTEGER :: hoursDay,days

    !------------------------------------------------------------------

    hoursDay=24
    days=time(1)/hoursDay
    timeH=time(1)-days*hoursDay
    WRITE(nOut,                                &
         &       '(/,1x,A," ",i2,"d : ",i2,"h : ",i2, &
         &       "m : ",i2,"s : ",i3,"ms",/)')        &
         &       text,days,timeH,time(2),time(3),time(4)

  end SUBROUTINE writeTime
  !------------------------------------------------------------------
end module timing

