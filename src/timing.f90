module timing
   !
   ! Useful functions for time-stepping
   !

   use iso_fortran_env, only: output_unit
#ifdef WITH_MPI
   use mpimod
#endif
   use precision_mod
   use parallel_mod

   implicit none

   private

   type, public :: timer_type
      integer :: n_counts
      real(cp) :: tStart
      real(cp) :: tTot
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: start_count
      procedure :: stop_count
   end type timer_type

   public :: formatTime

contains

   subroutine initialize(this)

      class(timer_type) :: this

      this%n_counts=0
      this%tStart = 0.0_cp
      this%tTot = 0.0_cp

   end subroutine initialize
!-----------------------------------------------------------------------------
   subroutine finalize(this, message, n_log_file)

      class(timer_type) :: this

      character(len=*), intent(in) :: message
      integer,          intent(in) :: n_log_file

      !-- Local variables
      integer :: n, n_out

      if ( this%n_counts > 0 ) this%tTot=this%tTot/real(this%n_counts,cp)

#ifdef WITH_MPI
      call MPI_AllReduce(MPI_IN_PLACE, this%tTot, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)
      this%tTot=this%tTot/real(n_procs,cp)
#endif

      if ( rank == 0 ) then
         if ( this%n_counts > 0 ) then
            do n=1,2
               if ( n == 1 ) then
                  n_out=n_log_file
               else if ( n == 2 ) then
                  n_out=output_unit
               end if
               call formatTime(n_out, message, this%tTot)
            end do
         end if
      end if

   end subroutine finalize
!-----------------------------------------------------------------------------
   subroutine start_count(this)

      class(timer_type) :: this

#ifndef WITH_MPI
      !-- Local variables
      integer :: counts, count_rate, count_max
#endif

#ifdef WITH_MPI
      this%tStart = MPI_Wtime()
#else
      call system_clock(counts, count_rate, count_max)
      this%tStart = real(counts,kind=cp) / real(count_rate,kind=cp)
#endif

   end subroutine start_count
!-----------------------------------------------------------------------------
   subroutine stop_count(this, l_increment)

      class(timer_type) :: this
      logical, optional, intent(in) :: l_increment

      logical :: l_count
      real(cp) :: tStop
#ifndef WITH_MPI
      !-- Local variables
      integer :: counts, count_rate, count_max
#endif

      if ( present(l_increment) ) then
         l_count=l_increment
      else
         l_count=.true.
      end if

#ifdef WITH_MPI
      tStop = MPI_Wtime()
#else
      call system_clock(counts, count_rate, count_max)
      tStop = real(counts,kind=cp) / real(count_rate,kind=cp)
#endif

      if ( tStop > this%tStart ) then
         this%tTot=this%tTot+(tStop-this%tStart)
         if ( l_count ) this%n_counts=this%n_counts+1
      end if

   end subroutine stop_count
!-----------------------------------------------------------------------------
   subroutine formatTime(n_out, message, time_in_sec)

      !-- Input variable
      integer,          intent(in) :: n_out
      character(len=*), intent(in) :: message
      real(cp),         intent(in) :: time_in_sec

      !-- Local variable
      integer :: n_time(5)
      real(cp) :: remaining_time

      n_time(1) = int(time_in_sec/86400_cp)
      remaining_time = time_in_sec-n_time(1)*86400_cp
      n_time(2) = int(remaining_time/3600_cp)
      remaining_time = remaining_time-n_time(2)*3600_cp
      n_time(3) = int(remaining_time/60_cp)
      remaining_time = remaining_time-n_time(3)*60_cp
      n_time(4) = int(remaining_time)
      remaining_time = remaining_time-n_time(4)
      n_time(5) = int(remaining_time/1e-3_cp)

      if ( time_in_sec > 0.1_cp .and. time_in_sec < 60.0_cp ) then
         write(n_out,'(1x,A,A,f6.3,A/,1x)') message, ' ', time_in_sec, " seconds"
      else if ( time_in_sec >= 60.0_cp .and. time_in_sec < 3.6e3_cp ) then
         write(n_out,'(1x,A,A,I2,A,I2,A,I3,A/,1x)') message, ' ', n_time(3), &
         &     " m ", n_time(4), " s ", n_time(5), " ms"
      else if ( time_in_sec >= 3.6e3_cp .and. time_in_sec < 8.64e4_cp ) then
         write(n_out,'(1x,A,A,I2,A,I2,A,I3,A,I3,A/,1x)') message, ' ',    &
         &     n_time(2), " h ", n_time(3), " m ", n_time(4), " s ",      &
         &     n_time(5), " ms"
      else if ( time_in_sec >= 8.64e4_cp ) then
         write(n_out,'(1x,A,A,I2,A,I2,A,I3,A,I3,A,I3,A/,1x)') message, ' ', &
         &     n_time(1), " d ", n_time(2), " h ", n_time(3), " m ",        &
         &     n_time(4), " s ", n_time(5), " ms"
      else
         write(n_out,'(1x,A,A,es10.3,A/,1x)') message, ' ', time_in_sec, " seconds"
      end if

   end subroutine formatTime
!-----------------------------------------------------------------------------
end module timing

