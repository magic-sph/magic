module mem_alloc
   !
   ! This little module is used to estimate the global memory allocation
   ! used in MagIC
   !

   use parallel_mod
   use truncation
   use precision_mod, only: lip, cp
   use output_data, only: tag

   implicit none

   private

   integer(lip), public :: bytes_allocated
   integer :: n_memory_file
   character(len=72) :: memory_file
   integer :: n_ranks_print
   integer, allocatable :: ranks_selected(:)

   public :: initialize_memory_counter, memWrite, finalize_memory_counter

contains

   subroutine initialize_memory_counter

      integer :: iproc

      bytes_allocated = 0 !
      n_ranks_print = min(n_ranks-1, 6)

      if ( n_ranks_print > 0 ) then
         allocate ( ranks_selected(n_ranks_print) )
         if ( n_ranks < 8 ) then
            ranks_selected = [(iproc,iproc=1,n_ranks-1)]
         else
            ranks_selected =[1,2,3,4,n_ranks-2,n_ranks-1]
         end if
      end if

      if ( l_master_rank ) then
         memory_file = 'mem_alloc.'//tag
         open(newunit=n_memory_file, file=memory_file, &
         &    status='unknown', position='append')
      end if

   end subroutine initialize_memory_counter
!------------------------------------------------------------------------------
   subroutine memWrite(origin, bytes_alloc)

      !-- Input variables
      character(len=*), intent(in) :: origin
      integer(lip),     intent(in) :: bytes_alloc

      !-- Local variables
      character(len=48) :: header
      character(len=12) :: st
#ifdef WITH_MPI
      integer :: i, iproc, sr_tag, status(MPI_STATUS_SIZE)
      integer(lip) :: bytes_other_proc
#endif

      header=' !----------------------------------------------'

      if ( l_master_rank ) then
         header(len(header)/2-len(origin)/2-1:len(header)/2+len(origin)/2+1) = &
                   ' '//origin//' '
         !write(n_memory_file, "(A20,A25)") ' !-----------------', origin
         write(n_memory_file, "(A48)") header
         st = human_readable_size(bytes_alloc)
         write(n_memory_file, "(A17,A,I4)") st," allocated on coord_r ", coord_r
      end if

#ifdef WITH_MPI
      sr_tag = 234876 ! arbitrary tag

      if  ( n_ranks_print /= 0 ) then

         do i=1,n_ranks_print
            iproc=ranks_selected(i)
            if ( coord_r == iproc ) then
               call MPI_Send(bytes_alloc,1,MPI_INTEGER8,0,sr_tag+iproc, &
                    &        mpi_comm_world,ierr)
            end if
         end do
         do i=1,n_ranks_print
            iproc=ranks_selected(i)
            if ( l_master_rank ) then
               call MPI_Recv(bytes_other_proc,1,MPI_INTEGER8,iproc,sr_tag+iproc, &
                    &        mpi_comm_world,status,ierr)
               if ( n_ranks > 8 .and. iproc == n_ranks -2 ) then
                  write(n_memory_file, *) "               ..."
               end if
               st = human_readable_size(bytes_other_proc)
               write(n_memory_file, "(A17,A,I4)") st, " allocated on coord_r ", iproc
            end if
         end do

      end if
#endif

      if ( l_master_rank ) write(n_memory_file, *) " "

   end subroutine memWrite
!------------------------------------------------------------------------------
   subroutine finalize_memory_counter

      !-- Local variables
      character(len=12) :: st
#ifdef WITH_MPI
      integer(lip) :: bytes_other_proc
      integer :: i, iproc, sr_tag, status(MPI_STATUS_SIZE)
#endif


      if ( l_master_rank ) then
         write(n_memory_file, *) " "
         write(n_memory_file, *) "!=============================================="
         write(n_memory_file, *) "!           TOTAL MEMORY ALLOCATION           !"
         write(n_memory_file, *) "!=============================================="
         write(n_memory_file, *) " "
         st = human_readable_size(bytes_allocated)
         write(n_memory_file, "(A17,A,I4)") st," allocated on coord_r ", coord_r
      end if

#ifdef WITH_MPI
      sr_tag =62985

      if  ( n_ranks_print /= 0 ) then

         do i=1,n_ranks_print
            iproc=ranks_selected(i)
            if ( coord_r == iproc ) then
               call MPI_Send(bytes_allocated,1,MPI_INTEGER8,0,sr_tag+iproc, &
                    &        mpi_comm_world,ierr)
            end if
         end do
         do i=1,n_ranks_print
            iproc=ranks_selected(i)
            if ( l_master_rank ) then
               call MPI_Recv(bytes_other_proc,1,MPI_INTEGER8,iproc,sr_tag+iproc, &
                    &        mpi_comm_world,status,ierr)
               if ( n_ranks > 8 .and. iproc == n_ranks-2 ) then
                  write(n_memory_file, *) "               ..."
               end if
               st = human_readable_size(bytes_other_proc)
               write(n_memory_file, "(A17,A,I4)") st, " allocated on coord_r ", iproc
            end if
         end do

      end if
#endif

      if ( l_master_rank ) then
         close(n_memory_file)
      end if

   end subroutine finalize_memory_counter
!------------------------------------------------------------------------------
   character(len=12) function human_readable_size(bytes) result(st)

      integer(lip), intent(in) :: bytes

      !-- Local variables
      character(len=1) :: suffix='B'
      character(len=2) :: units(6)
      integer :: i
      real(cp) :: bytes_float

      bytes_float = real(bytes, kind=cp)

      units = ['  ', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi']

      do i=1,6
         if ( bytes_float < 1024.0_cp ) then
            write(st, "(F8.3,A1,A,A)") bytes_float, ' ', units(i), suffix
            exit
         end if
         bytes_float = bytes_float / 1024.0_cp
      end do

   end function human_readable_size
!------------------------------------------------------------------------------
end module mem_alloc
