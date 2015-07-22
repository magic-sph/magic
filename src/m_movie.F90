!$Id$
module movie_data
!********************************************************************
!  Module containing work arrays
!********************************************************************

   use truncation
   use logic, only:  l_store_frame, l_save_out, l_movie, &
                     l_movie_oc, l_movie_ic, l_HTmovie,  &
                     l_dtBmovie, l_store_frame, l_save_out
   use parallel_mod
   use radial_data,only:nRstart,nRstop
 
   implicit none
   
   private
 
   real(kind=8), public :: movieDipColat,movieDipLon
   real(kind=8), public :: movieDipStrength,movieDipStrengthGeo
   real(kind=8), public :: t_movieS(10000)
 
   !-- Info in movie type and were the frames are stored:
   integer, public, parameter :: n_movies_max=30  ! Max no. of different movies
   integer, public, parameter :: n_movie_fields_max=6 ! Max no. of fields per movie
   real(kind=8), public ::  movie_const(n_movies_max)
   character(len=80), public :: movie(n_movies_max)  ! Only for input !
   character(len=72), public :: movie_file(n_movies_max)
 
   logical, public :: lStoreMov(n_movies_max),lICField(n_movies_max)
   integer, public :: n_movies
   integer, public :: n_movie_type(n_movies_max)
   integer, public :: n_movie_surface(n_movies_max)
   integer, public :: n_movie_const(n_movies_max)
   integer, public :: n_movie_fields(n_movies_max)
   integer, public :: n_movie_fields_ic(n_movies_max)
   integer, public :: n_movie_field_type(n_movie_fields_max,n_movies_max)
   integer, public :: n_movie_field_start(n_movie_fields_max,n_movies_max)
   integer, public :: n_movie_field_stop(n_movie_fields_max,n_movies_max)
 
   integer, public :: n_movie_file(n_movies_max)
 
   !-- Work arrays for storing movie frame:
   integer, public :: n_frame_work  
   integer, public :: n_MD
   real(kind=8), public, allocatable :: frames(:)

   public :: initialize_movie_data, finalize_movie_data, &
             movie_gather_frames_to_rank0

contains

   subroutine initialize_movie_data

      integer :: n

      n_MD=lMovieMem*n_movie_work*n_theta_max*n_phi_max*n_r_tot
      n_frame_work=max0(n_MD,1)

      !write(*,"(A,I10,A,F8.3,A)") "Allocating ",n_frame_work*8," bytes = ",&
      !     & (n_frame_work*8)/1024.0/1024.0,"MB for movie frames."
      allocate( frames(n_frame_work) )

      if ( .not. l_movie ) then
         l_movie_oc=.false.
         l_movie_ic=.false.
         l_HTmovie =.false.
         l_dtBmovie=.false.
         l_store_frame=.false.
      else
         call get_movie_type

         if ( rank == 0 ) then
            do n=1,n_movies_max
               n_movie_file(n)=70+n
            end do
            !----- Open movie files on first processor only:
            if ( .not. l_save_out ) then
               do n=1,n_movies
                  open(n_movie_file(n),file=movie_file(n), &
                       status='NEW',form='UNFORMATTED')
               end do
            end if
         end if
         if ( l_dtBmovie .and. ldtBMem == 0 ) then
            call logWrite('! You required dtB caculation.')
            call logWrite('! Please set ldtBMem=1 in truncation.f')
            call logWrite('! This is needed to reserve memory.')
            stop
         end if

      end if

   end subroutine initialize_movie_data
!----------------------------------------------------------------------------
   subroutine finalize_movie_data

      integer :: n

      if ( rank == 0 .and. l_movie ) then
         if ( l_movie ) then
            do n=1,n_movies
               close(n_movie_file(n))
            end do
         end if
      end if

   end subroutine finalize_movie_data
!----------------------------------------------------------------------------
   subroutine movie_gather_frames_to_rank0

      integer :: n_fields,n_surface,n_movie,n_const
      integer :: n_start,n_stop,n_field
      integer :: myTag, status(MPI_STATUS_SIZE)
      integer :: local_start,local_end,irank,sendcount
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(kind=8), allocatable :: field_frames_global(:)
      integer :: max_field_length,field_length
      
      max_field_length=0
      do n_movie=1,n_movies
         n_fields =n_movie_fields(n_movie)
         do n_field=1,n_fields
            n_start = n_movie_field_start(n_field,n_movie)
            n_stop  = n_movie_field_stop(n_field,n_movie)
            field_length=n_stop-n_start+1
            if (field_length > max_field_length) max_field_length=field_length
         end do
      end do
      if (rank == 0) then
         allocate(field_frames_global(max_field_length))
      else
         ! This is only needed for debug runs with boundary check.
         allocate(field_frames_global(1))
      end if

      ! loop over all movies
      do n_movie=1,n_movies
         n_fields =n_movie_fields(n_movie)
         n_surface=n_movie_surface(n_movie)
         n_const  =n_movie_const(n_movie)
         if ( n_surface == -1 ) then ! Earth Surface
            ! theta-phi surface for n_r=1 (CMB)
            ! frames is already existent on rank 0 with all
            ! needed values
            ! do nothing, pass to the next movie
            
            cycle

         else if ( n_surface == 0 ) then ! 3d
            ! 3d, all grid points written to frames
            ! but only n_r=nRstart:nRstop on one rank,
            ! gather needed
            do n_field=1,n_fields
               n_start = n_movie_field_start(n_field,n_movie)
               n_stop  = n_movie_field_stop(n_field,n_movie)
            end do

         else if ( n_surface == 1 ) then ! Surface r=constant
            ! frames is set only for one rank, where n_r=n_const
            ! send to rank 0
            n_start=n_movie_field_start(1,n_movie)
            n_stop =n_movie_field_stop(n_fields,n_movie)
            myTag=7654
            if (rank == 0) then
               if ((nRstart <= n_const) .and. (n_const <= nRstop)) then
                  ! relevant frames already set on rank 0
                  ! do nothing
               else
                  call MPI_Recv(frames(n_start),n_stop-n_start+1,MPI_DOUBLE_PRECISION,&
                       & MPI_ANY_SOURCE,mytag,MPI_COMM_WORLD,status,ierr)
               end if
            else
               if ((nRstart <= n_const) .and. (n_const <= nRstop)) then
                  ! relevant frames are all on this rank  /= 0
                  ! send to rank 0
                  call MPI_Send(frames(n_start),n_stop-n_start+1,MPI_DOUBLE_PRECISION,&
                       & 0,mytag,MPI_COMM_WORLD,ierr)
               end if
            end if

         else if ( n_surface == 2 ) then ! Surface theta=constant
            do n_field=1,n_fields
               n_start = n_movie_field_start(n_field,n_movie)
               n_stop  = n_movie_field_stop(n_field,n_movie)
               field_length = n_stop-n_start+1

               local_start=n_start+(nRstart-1)*n_phi_max
               local_end  =local_start+nr_per_rank*n_phi_max-1
               if (rank == n_procs-1) local_end = local_start+nr_on_last_rank*n_phi_max-1
               if (local_end > n_stop) then
                  write(*,"(A,2I7)") "local_end exceeds n_stop: ",local_end,n_stop
                  stop
               end if
               do irank=0,n_procs-1
                  recvcounts(irank) = nr_per_rank*n_phi_max
                  displs(irank)     = irank*nr_per_rank*n_phi_max
               end do
               recvcounts(n_procs-1) = nr_on_last_rank*n_phi_max
               sendcount=local_end-local_start+1

               call mpi_gatherv(frames(local_start),sendcount,MPI_DOUBLE_PRECISION,&
                    & field_frames_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
                    & 0,MPI_COMM_WORLD,ierr)
               if (rank == 0) then
                  frames(n_start:n_stop)=field_frames_global(1:field_length)
               end if
            end do  ! Do loop over field for one movie
         else if ( iabs(n_surface) == 3 ) then  ! Surface phi=const.
            ! all ranks have a part of the frames array for each movie
            ! we need to gather

            do n_field=1,n_fields
               n_start = n_movie_field_start(n_field,n_movie)
               n_stop  = n_movie_field_stop(n_field,n_movie)
               field_length = n_stop-n_start+1

               local_start=n_start+(nRstart-1)*n_theta_max
               local_end  =local_start+nr_per_rank*n_theta_max-1
               if (rank == n_procs-1) local_end = local_start+ &
                                                  nr_on_last_rank*n_theta_max-1
               if (local_end > n_stop) then
                  write(*,"(A,2I7)") "local_end exceeds n_stop: ",local_end,n_stop
                  stop
               end if
               do irank=0,n_procs-1
                  recvcounts(irank) = nr_per_rank*n_theta_max
                  displs(irank)     = irank*nr_per_rank*n_theta_max
               end do
               recvcounts(n_procs-1) = nr_on_last_rank*n_theta_max
               sendcount=local_end-local_start+1
               call mpi_gatherv(frames(local_start),sendcount,MPI_DOUBLE_PRECISION,&
                    & field_frames_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
                    & 0,MPI_COMM_WORLD,ierr)
               if (rank == 0) then
                  frames(n_start:n_stop)=field_frames_global(1:field_length)
               end if
            end do  ! Do loop over field for one movie

         end if
      end do

      if (rank == 0) deallocate(field_frames_global)

   end subroutine movie_gather_frames_to_rank0
!----------------------------------------------------------------------------
end module movie_data
