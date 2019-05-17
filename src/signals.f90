module signals_mod

   use parallel_mod
   use precision_mod
   use timing, only: time2ms
   use charmanip, only: capitalize
   use output_data, only: tag

   implicit none

   private

   integer :: n_sig_file
   real(cp) :: tsig

   public :: initialize_signals, check_signals

contains

   subroutine initialize_signals

      character(len=144) :: file_name

      tsig = 0.0_cp
      file_name = 'signal.'//tag
      open(newunit=n_sig_file, file=file_name, status='unknown')
      write(n_sig_file,'(A3)') 'NOT'
      close(n_sig_file)

   end subroutine initialize_signals
!------------------------------------------------------------------------------
   subroutine read_signal_file(signals)

      !-- Outputs signals
      integer, intent(inout) :: signals(5)

      !-- Local variables:
      character(len=255) :: message
      character(len=76) :: SIG


      signals(:) = 0
      if ( rank == 0 ) then
         !----- Signalling via file signal:
         message='signal'//'.'//tag
         open(newunit=n_sig_file, file=trim(message), status='old')
         read(n_sig_file,*) SIG
         close(n_sig_file)
         if ( len(trim(SIG)) > 0 ) then ! Non blank string ?

            call capitalize(SIG)

            if ( index(SIG,'END')/=0 ) signals(1)=1  !n_stop_signal=1

            if ( index(SIG,'GRA')/=0 ) then
               signals(2)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if

            if ( index(SIG,'RST')/=0 ) then
               signals(3)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if

            if ( index(SIG,'SPE')/=0 ) then
               signals(4)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if

            if ( index(SIG,'POT')/=0 ) then
               signals(5)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if
         end if
      end if

#ifdef WITH_MPI
      call MPI_Bcast(signals,5,MPI_Integer,0,MPI_COMM_WORLD,ierr)
#endif

   end subroutine read_signal_file
!------------------------------------------------------------------------------
   subroutine check_signals(run_time_passed, signals)

      !-- Input variable:
      integer, intent(in) :: run_time_passed(4)

      !-- Output variables
      integer, intent(out) :: signals(5)

      !-- Local variables:
      real(cp) :: run_time
      logical :: l_check_signal

      run_time = time2ms(run_time_passed)*1e-3_cp ! to seconds

      if ( run_time > 0.0_cp ) then
         tsig = tsig+run_time
      end if

      if ( rank == 0 ) then
         if ( tsig > 2.0_cp ) then ! Only check signals every two seconds
            l_check_signal = .true.
            tsig = 0.0_cp
         else
            l_check_signal = .false.
         end if
      end if

#ifdef WITH_MPI
      call MPI_Bcast(l_check_signal,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

      if ( l_check_signal ) then
#ifdef WITH_MPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call read_signal_file(signals)
      else
         signals(:) = 0
      end if

   end subroutine check_signals
!--------------------------------------------------------------------------------
end module signals_mod
