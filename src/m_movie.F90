!$Id$
!********************************************************************
!  Module containing work arrays
!********************************************************************

!------------ This is release 1 level 1  --------------!
!------------ Created on 1/17/02  by JW. --------------!

MODULE movie_data
  use truncation
  USE parallel_mod
  USE radial_data,ONLY:nRstart,nRstop
  IMPLICIT NONE

  REAL(kind=8) :: movieDipColat,movieDipLon
  REAL(kind=8) :: movieDipStrength,movieDipStrengthGeo
  REAL(kind=8) :: t_movieS(10000)

  !-- Info in movie type and were the frames are stored:
  INTEGER,PARAMETER :: n_movies_max=30  ! Max no. of different movies
  INTEGER,PARAMETER :: n_movie_fields_max=6 ! Max no. of fields per movie
  REAL(kind=8) ::  movie_const(n_movies_max)
  CHARACTER(len=80) :: movie(n_movies_max)  ! Only for input !
  logical lStoreMov(n_movies_max),lICField(n_movies_max)
  INTEGER :: n_movies
  INTEGER :: n_movie_type(n_movies_max)
  INTEGER :: n_movie_surface(n_movies_max)
  INTEGER :: n_movie_const(n_movies_max)
  INTEGER :: n_movie_fields(n_movies_max)
  INTEGER :: n_movie_fields_ic(n_movies_max)
  INTEGER :: n_movie_field_type(n_movie_fields_max,n_movies_max)
  INTEGER :: n_movie_field_start(n_movie_fields_max,n_movies_max)
  INTEGER :: n_movie_field_stop(n_movie_fields_max,n_movies_max)

  INTEGER :: n_movie_file(n_movies_max)

  !-- Work arrays for storing movie frame:
  INTEGER :: n_frame_work  
  INTEGER :: n_MD
  REAL(kind=8),ALLOCATABLE :: frames(:)
CONTAINS
  SUBROUTINE initialize_movie_data

    n_MD=lMovieMem*n_movie_work*n_theta_max*n_phi_max*n_r_tot
    n_frame_work=MAX0(n_MD,1)

    !WRITE(*,"(A,I10,A,F8.3,A)") "Allocating ",n_frame_work*8," bytes = ",&
    !     & (n_frame_work*8)/1024.0/1024.0,"MB for movie frames."
    ALLOCATE( frames(n_frame_work) )

  END SUBROUTINE initialize_movie_data

  SUBROUTINE movie_gather_frames_to_rank0
    INTEGER :: n_fields,n_surface,n_movie,n_const
    INTEGER :: n_start,n_stop,n_field
    INTEGER :: myTag, status(MPI_STATUS_SIZE)
    INTEGER :: local_start,local_end,irank,sendcount
    INTEGER,DIMENSION(0:n_procs-1) :: recvcounts,displs
    REAL(kind=8), DIMENSION(:),allocatable :: field_frames_global
    INTEGER :: max_field_length,field_length
    
    max_field_length=0
    DO n_movie=1,n_movies
       n_fields =n_movie_fields(n_movie)
       DO n_field=1,n_fields
          n_start = n_movie_field_start(n_field,n_movie)
          n_stop  = n_movie_field_stop(n_field,n_movie)
          field_length=n_stop-n_start+1
          if (field_length > max_field_length) max_field_length=field_length
       END DO
    END DO
    IF (rank == 0) THEN
       ALLOCATE(field_frames_global(max_field_length))
    ELSE
       ! This is only needed for debug runs with boundary check.
       ALLOCATE(field_frames_global(1))
    END IF

    ! loop over all movies
    DO n_movie=1,n_movies
       n_fields =n_movie_fields(n_movie)
       n_surface=n_movie_surface(n_movie)
       n_const  =n_movie_const(n_movie)
       !WRITE(*,"(A,I2,A,3I3)") "Movie nr. ",n_movie,", n_surface,n_fields,n_const = ",&
       !     & n_surface,n_fields,n_const
       IF ( n_surface == -1 ) THEN ! Earth Surface
          ! theta-phi surface for n_r=1 (CMB)
          ! frames is already existent on rank 0 with all
          ! needed values
          ! do nothing, pass to the next movie
          !do n_field=1,n_fields
          !   n_start = n_movie_field_start(n_field,n_movie)
          !   n_stop  = n_movie_field_stop(n_field,n_movie)
             !WRITE(*,"(2(A,I5),A)") "No communication, using frames(",n_start,":",n_stop,")"
          !end do
          
          CYCLE

       ELSE IF ( n_surface == 0 ) THEN ! 3d
          ! 3d, all grid points written to frames
          ! but only n_r=nRstart:nRstop on one rank,
          ! gather needed
          do n_field=1,n_fields
             n_start = n_movie_field_start(n_field,n_movie)
             n_stop  = n_movie_field_stop(n_field,n_movie)
          end do

       ELSE IF ( n_surface == 1 ) THEN ! Surface r=constant
          ! frames is set only for one rank, where n_r=n_const
          ! send to rank 0
          n_start=n_movie_field_start(1,n_movie)
          n_stop =n_movie_field_stop(n_fields,n_movie)
          myTag=7654
          if (rank == 0) then
             IF ((nRstart <= n_const) .AND. (n_const <= nRstop)) THEN
                ! relevant frames already set on rank 0
                ! do nothing
             ELSE
                CALL MPI_Recv(frames(n_start),n_stop-n_start+1,MPI_DOUBLE_PRECISION,&
                     & MPI_ANY_SOURCE,mytag,MPI_COMM_WORLD,status,ierr)
             END IF
          ELSE
             IF ((nRstart <= n_const) .AND. (n_const <= nRstop)) THEN
                ! relevant frames are all on this rank .ne.0
                ! send to rank 0
                CALL MPI_Send(frames(n_start),n_stop-n_start+1,MPI_DOUBLE_PRECISION,&
                     & 0,mytag,MPI_COMM_WORLD,ierr)
             END IF
          END IF

       else if ( n_surface == 2 ) then ! Surface theta=constant
          DO n_field=1,n_fields
             n_start = n_movie_field_start(n_field,n_movie)
             n_stop  = n_movie_field_stop(n_field,n_movie)
             field_length = n_stop-n_start+1

             local_start=n_start+(nRstart-1)*n_phi_max
             local_end  =local_start+nr_per_rank*n_phi_max-1
             IF (rank == n_procs-1) local_end = local_start+nr_on_last_rank*n_phi_max-1
             IF (local_end > n_stop) THEN
                WRITE(*,"(A,2I7)") "local_end exceeds n_stop: ",local_end,n_stop
                STOP
             END IF
             DO irank=0,n_procs-1
                recvcounts(irank) = nr_per_rank*n_phi_max
                displs(irank)     = irank*nr_per_rank*n_phi_max
             END DO
             recvcounts(n_procs-1) = nr_on_last_rank*n_phi_max
             sendcount=local_end-local_start+1

             CALL mpi_gatherv(frames(local_start),sendcount,MPI_DOUBLE_PRECISION,&
                  & field_frames_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
                  & 0,MPI_COMM_WORLD,ierr)
             IF (rank == 0) THEN
                frames(n_start:n_stop)=field_frames_global(1:field_length)
             END IF
          END DO  ! Do loop over field for one movie
       else if ( IABS(n_surface) == 3 ) then  ! Surface phi=const.
          ! all ranks have a part of the frames array for each movie
          ! we need to gather

          DO n_field=1,n_fields
             n_start = n_movie_field_start(n_field,n_movie)
             n_stop  = n_movie_field_stop(n_field,n_movie)
             field_length = n_stop-n_start+1

             local_start=n_start+(nRstart-1)*n_theta_max
             local_end  =local_start+nr_per_rank*n_theta_max-1
             IF (rank == n_procs-1) local_end = local_start+nr_on_last_rank*n_theta_max-1
             IF (local_end > n_stop) THEN
                WRITE(*,"(A,2I7)") "local_end exceeds n_stop: ",local_end,n_stop
                STOP
             END IF
             DO irank=0,n_procs-1
                recvcounts(irank) = nr_per_rank*n_theta_max
                displs(irank)     = irank*nr_per_rank*n_theta_max
                !IF (rank == 0) WRITE(*,"(A,I5,A,I3,A,I5,A)") "Receiving ",recvcounts(irank),&
                !     &" vals from rank ",irank,", writing to frames(",&
                !     & n_start+displs(irank),")"
             END DO
             recvcounts(n_procs-1) = nr_on_last_rank*n_theta_max
             sendcount=local_end-local_start+1
             !WRITE(*,"(2(A,I5),A,ES22.12,2I5)") "Sending frames(",local_start,":",local_end,") = ",&
             !     & SUM(frames(local_start:local_end)),&
             !     &n_start+displs(rank),n_start+displs(rank)+recvcounts(rank)-1

             ! MPI_IN_PLACE would be better, but does not work with Intel MPI 4.1
             CALL mpi_gatherv(frames(local_start),sendcount,MPI_DOUBLE_PRECISION,&
                  & field_frames_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
                  & 0,MPI_COMM_WORLD,ierr)
             IF (rank == 0) THEN
                frames(n_start:n_stop)=field_frames_global(1:field_length)
                !DO irank=0,n_procs-1
                 !  WRITE(*,"(A,I2,2(A,I5),A,ES22.12)") "irank=",irank,&
                 !       &": After gatherv on rank 0, frames(",&
                 !       &n_start+displs(irank),":",n_start+displs(irank)+recvcounts(irank)-1,&
                 !       &") = ",SUM(frames(n_start+displs(irank):n_start+displs(irank)+recvcounts(irank)-1))
                !END DO
             END IF
          END DO  ! Do loop over field for one movie

       END IF
    END DO
    if (rank == 0) DEALLOCATE(field_frames_global)
  END SUBROUTINE movie_gather_frames_to_rank0
END MODULE movie_data
