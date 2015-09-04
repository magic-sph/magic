!$Id$
module timing

#ifdef WITH_MPI
   use mpi
#endif
   use precision_mod, only: cp, lip
   use parallel_mod, only: rank
 
   implicit none

   private

   integer, parameter :: msecSecond=1000
   integer, parameter :: msecMinute=60000
   integer, parameter :: msecHour  =3600000

   public :: wallTime, ms2time, time2ms, subTime, lTimeLimit, &
             addTime, lNegTime, writeTime, meanTime

contains
   subroutine wallTime(time)
      !------------------------------------------------------------------
      !   This routine returns the wallclock time in four
      !   integer arguments.
      !------------------------------------------------------------------
  
      !-- Output variables
      integer, intent(out) :: time(4)

      !-- Local variables
      integer(lip) :: mSeconds
      real(cp) :: dbl_seconds
#ifndef WITH_MPI
      integer :: count, countRate, countMax
#endif
  
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
  
#ifdef WITH_MPI
      dbl_seconds = MPI_Wtime()
#else
      call SYSTEM_CLOCK(count,countRate,countMax)
      dbl_seconds=real(count,kind=cp)/real(countRate,kind=cp)
#endif
  
      mSeconds = int(1000.0_cp*dbl_seconds,kind=lip)
      call ms2time(mSeconds,time)
 
   end subroutine wallTime
!----------------------------------------------------------------------------
   subroutine ms2time(ms,time)
      !------------------------------------------------------------------
      !  Transforms accumulated milliseconds ms into an four-element
      !  integer arrays time(4) containing the time in
      !  hours=time(1), minutes=time(2), seconds=time(3),
      !  and milliseconds=time(4).
      !------------------------------------------------------------------
  
      !-- Input variables:
      integer(lip), intent(in) :: ms
 
      !-- Output variables:
      integer, intent(out) :: time(4)
 
      !-- Local variables:
      integer(lip) :: mSeconds
      integer :: seconds,minutes,hours
  
      mSeconds=ms
      hours   =int(mSeconds/msecHour)
      mSeconds=mSeconds-hours*msecHour
      minutes =int(mSeconds/msecMinute)
      mSeconds=mSeconds-minutes*msecMinute
      seconds =int(mSeconds/msecSecond)
      mSeconds=mSeconds-seconds*msecSecond
      time(1) =hours
      time(2) =minutes
      time(3) =seconds
      time(4) =int(mSeconds)
 
   end subroutine ms2time
!----------------------------------------------------------------------------
   integer(lip) function time2ms(time)
      !------------------------------------------------------------------
      !  Transforms a four-element integer arrays time(4) containing the
      !  time in hours=time(1), minutes=time(2), seconds=time(3),
      !  and milliseconds=time(4) into accumulated milliseconds.
      !------------------------------------------------------------------
  
      !-- Input variable:
      integer, intent(in) :: time(4)

      time2ms=time(1)*msecHour+time(2)*msecMinute+time(3)*msecSecond+time(4)
  
   end function time2ms
!----------------------------------------------------------------------------
   subroutine subTime(timeStart,timeStop,timeD)
      !------------------------------------------------------------------
      !  Returns time passed between timeStop and timeStart.
      !  Note timeStop has to be younger than timeStart, otherwise
      !  24 hours are added. This is necessary on systems like the IBM
      !  where the time counter as reset every day at midnight.
      !------------------------------------------------------------------
  
      !-- Input variables:
      integer, intent(in) :: timeStart(4), timeStop(4)
 
      !-- Output variable:
      integer, intent(out) :: timeD(4)
  
      !-- Local variables:
      integer(lip) :: msDiff,msStart,msStop
  
      msStart=time2ms(timeStart)
      msStop =time2ms(timeStop)
  
      msDiff=msStop-msStart
      call ms2time(msDiff,timeD)
 
   end subroutine subTime
!----------------------------------------------------------------------------
   logical function lTimeLimit(time,timeMax)
      !------------------------------------------------------------------
      !  True when time exeeds timeMax
      !------------------------------------------------------------------
  
      !-- Input variables
      integer, intent(in) :: time(4),timeMax(4)
  
      lTimeLimit=.false.
      if ( time(1) > timeMax(1) ) then
         lTimeLimit=.true.
      else if ( time(1) == timeMax(1) ) then
         if ( time(2) > timeMax(2) ) then
            lTimeLimit=.true.
         else if ( time(2) == timeMax(2) ) then
            if ( time(3) > timeMax(3) ) then
               lTimeLimit=.true.
            else if ( time(3) == timeMax(3) ) then
               if ( time(4) >= timeMax(4) ) lTimeLimit= .true. 
            end if
         end if
      end if
  
   end function lTimeLimit
!----------------------------------------------------------------------------
   subroutine addTime(time1,time2)
  
      !-- Input variables
      integer :: time1(4),time2(4)
 
      !-- Local variables:
      integer(lip) :: t1ms,t2ms
  
      t1ms=time2ms(time1)
      t2ms=time2ms(time2)
      t1ms=t1ms+t2ms
      call ms2time(t1ms,time1)
 
   end subroutine addTime
!----------------------------------------------------------------------------
   subroutine meanTime(time,n)
 
      !-- Input variables
      integer :: time(4)
      integer, intent(in) :: n
 
      !-- Local variables
      integer(lip) :: tms
  
      tms=time2ms(time)
      if ( n==0 ) then
         if ( rank == 0 ) then
            write(*,*) '! Error in meanTime !'
            write(*,*) '! Cant build mean for zero steps!'
         end if
         tms=0
      else
         tms=tms/n
      end if
      call ms2time(tms,time)
  
   end subroutine meanTime
!----------------------------------------------------------------------------
   logical function lNegTime(time1,time2)
      !------------------------------------------------------------------
      !  Negative passed time? Means we have passed midnight.
      !  The wallclock time is reset to zero on some
      !  computers at midnight.
      !------------------------------------------------------------------
  
      !-- Input variables
      integer, intent(in) :: time1(4),time2(4)
 
      !-- Local variables
      integer(lip) :: tms1,tms2
  
      tms1=time2ms(time1)
      tms2=time2ms(time2)
      if ( tms2 < tms1-1 ) then
         lNegTime=.true.
      else
         lNegTime=.false.
      end if
 
   end function lNegTime
!----------------------------------------------------------------------------
   subroutine writeTime(nOut,text,time)
      !------------------------------------------------------------------
      !  Returns time passed between timeStop and timeStart.
      !  Note timeStop has to be younger than timeStart, otherwise
      !  24 hours are added. This is necessary on systems like the IBM
      !  where the time counter are reset every day at midnight.
      !------------------------------------------------------------------
  
      !-- Input variables:
      integer,          intent(in) :: nOut
      character(len=*), intent(in) :: text
      integer,          intent(in) :: time(4)
  
      !-- Local variables:
      integer :: timeH
      integer :: hoursDay,days
  
      hoursDay=24
      days=time(1)/hoursDay
      timeH=time(1)-days*hoursDay
      if (nOut == 6) then
         write(*,'(/,1X,2A,I4,A,I2,A,I2,A,I2,A,I3,A)')         &
              &     text," ",days,"d : ",timeH,"h : ",time(2), &
              &       "m : ",time(3),"s : ",time(4),"ms"
      else
         write(nOut,'(/,1x,A,A,i4,A,i2,A,i2,A,i2,A,i3,A)')     &
              &     text," ",days,"d : ",timeH,"h : ",time(2), &
              &       "m : ",time(3),"s : ",time(4),"ms"
      end if
 
   end subroutine writeTime
!----------------------------------------------------------------------------
end module timing

