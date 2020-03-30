module signals_mod
   !
   ! This module handles the reading/writing of the signal.TAG file which
   ! allows to communicate with MagIC during its execution
   !


   use parallel_mod
   use precision_mod
   use charmanip, only: capitalize
   use output_data, only: tag

   implicit none

   private

   integer :: n_sig_file ! File handle
   real(cp) :: tsig      ! A timer to check the signal only every 2 seconds

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
      if ( l_master_rank ) then
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
      call MPI_Bcast(signals,5,MPI_Integer,0,comm_r,ierr)
#endif

   end subroutine read_signal_file
!------------------------------------------------------------------------------
   subroutine check_signals(run_time, signals)

      !-- Input variable:
      real(cp), intent(in) :: run_time

      !-- Output variables
      integer, intent(out) :: signals(5)

      !-- Local variables:
      logical :: l_check_signal

      if ( run_time > 0.0_cp ) then
         tsig = tsig+run_time
      end if

      if ( coord_r == 0 ) then
         if ( tsig > 2.0_cp ) then ! Only check signals every two seconds
            l_check_signal = .true.
            tsig = 0.0_cp
         else
            l_check_signal = .false.
         end if
      end if

#ifdef WITH_MPI
      call MPI_Bcast(l_check_signal,1,MPI_LOGICAL,0,comm_r,ierr)
#endif

      if ( l_check_signal ) then
#ifdef WITH_MPI
         call MPI_Barrier(comm_r,ierr)
#endif
         !-- Read the signal file
         call read_signal_file(signals)

#ifdef WITH_MPI
         call MPI_Barrier(comm_r,ierr)
#endif
      else
         signals(:) = 0
      end if

   end subroutine check_signals
!--------------------------------------------------------------------------------
end module signals_mod
