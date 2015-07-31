!$Id$
#include "perflib_preproc.cpp"
MODULE communications
  use mpi
  USE parallel_mod, ONLY: rank, n_procs, ierr, nr_per_rank, nr_on_last_rank
  USE LMLoop_data, ONLY: llm, ulm
  USE truncation, ONLY: l_max, lm_max, minc, n_r_max, n_r_ic_max
  USE blocking, ONLY: st_map, lo_map, lmStartB, lmStopB
  USE radial_data, ONLY: nRstart, nRstop
  USE logic, ONLY: l_mag, l_conv, l_heat

  implicit none

  private 

  interface get_global_sum
     module procedure get_global_sum_cmplx_2d, get_global_sum_cmplx_1d, &
                      get_global_sum_real_2d
  end interface

  type, public :: lm2r_type
     INTEGER, ALLOCATABLE :: final_wait_array(:), s_request(:), r_request(:)
     COMPLEX(kind=8), pointer :: temp_Rloc(:,:,:), arr_Rloc(:,:,:)
     INTEGER :: count
  end type lm2r_type

  type, public :: gather_type
     INTEGER, ALLOCATABLE :: gather_mpi_type(:)
     INTEGER :: dim2
  end type gather_type

  ! MPI datatypes for the redistribution of the d?dt arrays
  INTEGER, SAVE, ALLOCATABLE :: s_transfer_type(:),s_transfer_type_nr_end(:)
  INTEGER, SAVE, ALLOCATABLE :: r_transfer_type(:), r_transfer_type_nr_end(:)
  INTEGER, SAVE, ALLOCATABLE :: s_transfer_type_cont(:,:)
  INTEGER, SAVE, ALLOCATABLE :: s_transfer_type_nr_end_cont(:,:)
  INTEGER, SAVE, ALLOCATABLE :: r_transfer_type_cont(:,:)
  INTEGER, SAVE, ALLOCATABLE :: r_transfer_type_nr_end_cont(:,:)
  INTEGER :: r_lm_gather_type, r_lm_gather_type_lm_end
  INTEGER, ALLOCATABLE :: s_request(:),r_request(:),final_wait_array(:)
  INTEGER, ALLOCATABLE :: array_of_statuses(:,:)

  public :: gather_from_lo_to_rank0,scatter_from_rank0_to_lo,&
          & gather_all_from_lo_to_rank0
  public :: get_global_sum, myAllGather,&
          & r2lm_redist,r2lo_redist,initialize_communications,&
          & create_lm2r_type!,lo2r_redist,lm2r_redist
  public :: lo2r_redist_start,lo2r_redist_wait

  ! declaration of the types for the redistribution
  !TYPE(lm2r_type),PUBLIC :: lo2r_s, lo2r_ds, lo2r_z, lo2r_dz
  !TYPE(lm2r_type),public :: lo2r_p,lo2r_dp
  !TYPE(lm2r_type), PUBLIC :: lo2r_b, lo2r_db, lo2r_ddb, lo2r_aj, lo2r_dj
  TYPE(lm2r_type), public :: lo2r_w, lo2r_s, lo2r_z, lo2r_p, lo2r_b, lo2r_aj

  TYPE(gather_type), public :: gt_OC,gt_IC,gt_cheb

  COMPLEX(kind=8), ALLOCATABLE :: temp_gather_lo(:)
  COMPLEX(kind=8), ALLOCATABLE :: temp_r2lo(:,:)

contains
  
  SUBROUTINE initialize_communications
    INTEGER :: proc,my_lm_per_rank
    INTEGER(KIND=MPI_ADDRESS_KIND) :: zerolb, extent, sizeof_double_complex
    INTEGER(KIND=MPI_ADDRESS_KIND) :: lb_marker, myextent, true_lb, true_extent
    INTEGER :: base_col_type,temptype
    INTEGER :: blocklengths(3),blocklengths_on_last(3),displs(3),displs_on_last(3)
    integer :: i

    ! first setup the datatype. It is not equal for all ranks. The n_procs-1 rank can
    ! have a smaller datatype.
    ! Due to the different number of radial and lm points for the ranks, we need essentially 
    ! three different datatypes
    ! transfer_type: Standard for the transfers between ranks 0-(n_procs-2)
    ! transfer_type_nr_end: for transfers involving rank (n_procs-1) as receiver
    ! transfer_type_lm_end: for transfers involving rank (n_procs-1) as sender
    ! +----+----+----+----+
    ! |    |    |    |    |
    ! |    |    |    |    |
    ! +----+----+----+----+
    ! |    |    |    |    |
    ! |    |    |    |    |
    ! +----+----+----+----+
    ! nr_per_rank does already exist
    ! lm_per_rank is set here
    ! ATTENTION: for the last rank, the numbers are different and are
    !            stored in nr_on_last_rank and lm_on_last_rank
    ALLOCATE(s_transfer_type(n_procs))
    ALLOCATE(s_transfer_type_nr_end(n_procs))
    ALLOCATE(r_transfer_type(n_procs))
    ALLOCATE(r_transfer_type_nr_end(n_procs))

    ALLOCATE(s_transfer_type_cont(n_procs,3))
    ALLOCATE(s_transfer_type_nr_end_cont(n_procs,3))
    ALLOCATE(r_transfer_type_cont(n_procs,3))
    ALLOCATE(r_transfer_type_nr_end_cont(n_procs,3))

    DO proc=0,n_procs-1
       my_lm_per_rank=lmStopB(proc+1)-lmStartB(proc+1)+1
       !WRITE(*,"(2(A,I4))") "lm_per_rank on rank ", proc," is ",my_lm_per_rank
       CALL MPI_Type_vector(nr_per_rank,my_lm_per_rank,&
            &lm_max,MPI_DOUBLE_COMPLEX,s_transfer_type(proc+1),ierr)
       CALL MPI_Type_commit(s_transfer_type(proc+1),ierr)

       ! The same for the last rank for nR
       CALL MPI_Type_vector(nr_on_last_rank,my_lm_per_rank,&
            &lm_max,MPI_DOUBLE_COMPLEX,s_transfer_type_nr_end(proc+1),ierr)
       CALL MPI_Type_commit(s_transfer_type_nr_end(proc+1),ierr)

       ! we do not need special receive datatypes, as the buffers are contiguous in memory
       ! but for ease of reading, we define the receive datatypes explicitly
       CALL MPI_Type_contiguous(my_lm_per_rank*nr_per_rank,&
            & MPI_DOUBLE_COMPLEX,r_transfer_type(proc+1),ierr)
       CALL MPI_Type_commit(r_transfer_type(proc+1),ierr)
       CALL MPI_Type_contiguous(my_lm_per_rank*nr_on_last_rank,&
            &MPI_DOUBLE_COMPLEX,r_transfer_type_nr_end(proc+1),ierr)
       CALL MPI_Type_commit(r_transfer_type_nr_end(proc+1),ierr)


       ! define the transfer types for the containers
       ! same schema as for the other types
       ! some temporary datatypes, not needed for communication
       ! but only for constructing the final datatypes
       CALL MPI_Type_get_extent(MPI_DOUBLE_COMPLEX,zerolb,sizeof_double_complex,ierr)
       CALL MPI_Type_contiguous(my_lm_per_rank,MPI_DOUBLE_COMPLEX,temptype,ierr)
       zerolb=0
       extent=lm_max*sizeof_double_complex
       CALL MPI_Type_create_resized(temptype,zerolb,extent,base_col_type,ierr)
       !CALL MPI_type_get_extent(base_col_type,lb_marker,myextent,ierr)
       !WRITE(*,"(2(A,I10))") "base_col_type: lb = ",lb_marker,", extent = ",myextent
       blocklengths = (/ nr_per_rank, nr_per_rank, nr_per_rank /)
       displs       = (/ 0,           nr_per_rank, 2*nr_per_rank   /)
       blocklengths_on_last = (/nr_on_last_rank,nr_on_last_rank, nr_on_last_rank/)
       displs_on_last       = (/ 0,             nr_on_last_rank, 2*nr_on_last_rank /)
       DO i=1,3
          CALL MPI_Type_vector(i,nr_per_rank*my_lm_per_rank,n_r_max*my_lm_per_rank,&
               & MPI_DOUBLE_COMPLEX,r_transfer_type_cont(proc+1,i),ierr)
          CALL MPI_Type_commit(r_transfer_type_cont(proc+1,i),ierr)
          CALL MPI_Type_vector(i,nr_on_last_rank*my_lm_per_rank,n_r_max*my_lm_per_rank,&
               & MPI_DOUBLE_COMPLEX,r_transfer_type_nr_end_cont(proc+1,i),ierr)
          CALL MPI_Type_commit(r_transfer_type_nr_end_cont(proc+1,i),ierr)

          CALL MPI_Type_indexed(i,blocklengths(1:i),&
               & displs(1:i),base_col_type,s_transfer_type_cont(proc+1,i),ierr)
          CALL MPI_Type_commit(s_transfer_type_cont(proc+1,i),ierr)
          CALL MPI_Type_indexed(i,blocklengths_on_last(1:i),&
               & displs_on_last(1:i),base_col_type,s_transfer_type_nr_end_cont(proc+1,i),ierr)
          CALL MPI_Type_commit(s_transfer_type_nr_end_cont(proc+1,i),ierr)

#if 0
          IF (i == 3) THEN
             CALL MPI_type_get_extent(r_transfer_type_cont(proc+1,i),lb_marker,myextent,ierr)
             CALL MPI_type_get_true_extent(r_transfer_type_cont(proc+1,i),true_lb,true_extent,ierr)
             WRITE(*,"(2(A,I3),3(A,I10))") "r_transfer_type_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                  & ", extent = ",myextent,", true extent = ",true_extent
             CALL MPI_type_get_extent(s_transfer_type_cont(proc+1,i),lb_marker,myextent,ierr)
             CALL MPI_type_get_true_extent(s_transfer_type_cont(proc+1,i),true_lb,true_extent,ierr)
             WRITE(*,"(2(A,I3),3(A,I10))") "s_transfer_type_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                  & ", extent = ",myextent,", true extent = ",true_extent
             
             CALL MPI_type_get_extent(r_transfer_type_nr_end_cont(proc+1,i),lb_marker,myextent,ierr)
             CALL MPI_type_get_true_extent(r_transfer_type_nr_end_cont(proc+1,i),true_lb,true_extent,ierr)
             WRITE(*,"(2(A,I3),3(A,I10))") "r_transfer_type_nr_end_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                  & ", extent = ",myextent,", true extent = ",true_extent
             CALL MPI_type_get_extent(s_transfer_type_nr_end_cont(proc+1,i),lb_marker,myextent,ierr)
             CALL MPI_type_get_true_extent(s_transfer_type_nr_end_cont(proc+1,i),true_lb,true_extent,ierr)
             WRITE(*,"(2(A,I3),3(A,I10))") "s_transfer_type_nr_end_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                  & ", extent = ",myextent,", true extent = ",true_extent
          END IF

#endif
       END DO
    END DO


    CALL create_gather_type(gt_OC,n_r_max)
    CALL create_gather_type(gt_IC,n_r_ic_max)

    ALLOCATE(s_request(n_procs-1),r_request(n_procs-1))
    ALLOCATE(array_of_statuses(MPI_STATUS_SIZE,2*(n_procs-1)))
    ALLOCATE(final_wait_array(2*(n_procs-1)))

    IF (l_heat) THEN
       CALL create_lm2r_type(lo2r_s,2)
       !CALL create_lm2r_type(lo2r_ds)
    END IF
    IF (l_conv) THEN
       CALL create_lm2r_type(lo2r_z,2)
       !CALL create_lm2r_type(lo2r_dz)
       CALL create_lm2r_type(lo2r_w,3)
       !CALL create_lm2r_type(lo2r_w)
       !CALL create_lm2r_type(lo2r_dw)
       !CALL create_lm2r_type(lo2r_ddw)
       CALL create_lm2r_type(lo2r_p,2)
       !CALL create_lm2r_type(lo2r_dp)
    END IF

    IF (l_mag) THEN
       CALL create_lm2r_type(lo2r_b,3)
       !CALL create_lm2r_type(lo2r_db)
       !CALL create_lm2r_type(lo2r_ddb)
       CALL create_lm2r_type(lo2r_aj,2)
       !CALL create_lm2r_type(lo2r_dj)
    END IF


    ! allocate a temporary array for the gather operations.
    ALLOCATE(temp_r2lo(lm_max,nRstart:nRstop))
    IF (rank == 0) THEN
       ALLOCATE(temp_gather_lo(1:lm_max))
    ELSE
       allocate(temp_gather_lo(1))
    END IF
  END SUBROUTINE initialize_communications


  COMPLEX(kind=8) FUNCTION get_global_sum_cmplx_2d(dwdt_local) RESULT(global_sum)
    COMPLEX(kind=8), INTENT(IN) :: dwdt_local(:,:)
    
    INTEGER :: ierr
    COMPLEX(kind=8) :: local_sum
    
    local_sum = SUM( dwdt_local )
    CALL MPI_Allreduce(local_sum,global_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  END FUNCTION get_global_sum_cmplx_2d

  REAL(kind=8) FUNCTION get_global_sum_real_2d(dwdt_local) RESULT(global_sum)
    real(kind=8), INTENT(IN) :: dwdt_local(:,:)
    
    INTEGER :: ierr
    real(kind=8) :: local_sum
    
    local_sum = SUM( dwdt_local )
    CALL MPI_Allreduce(local_sum,global_sum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  END FUNCTION get_global_sum_real_2d

  COMPLEX(kind=8) FUNCTION get_global_sum_cmplx_1d(arr_local) RESULT(global_sum)
    COMPLEX(kind=8), INTENT(IN) :: arr_local(:)
    
    INTEGER :: lb,ub,ierr,i
    COMPLEX(kind=8) :: local_sum,c,y,t

    lb = LBOUND(arr_local,1)
    ub = UBOUND(arr_local,1)
    
    ! Kahan summation algorithm
    !function KahanSum(input)
    !var sum = 0.0
    !var c = 0.0          //A running compensation for lost low-order bits.
    !for i = 1 to input.length do
    !    y = input[i] - c    //So far, so good: c is zero.
    !    t = sum + y         //Alas, sum is big, y small, so low-order digits of y are lost.
    !    c = (t - sum) - y   //(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
    !    sum = t             //Algebraically, c should always be zero. Beware eagerly optimising compilers!
    !    //Next time around, the lost low part will be added to y in a fresh attempt.
    !return sum
    
    local_sum = 0.0D0
    c = 0.0D0          !A running compensation for lost low-order bits.
    DO i=lb,ub
       y = arr_local(i) - c    !So far, so good: c is zero.
       t = local_sum + y       !Alas, sum is big, y small, so low-order digits of y are lost.
       c = (t - local_sum) - y !(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
       local_sum = t           !Algebraically, c should always be zero. Beware eagerly optimising compilers!
       !Next time around, the lost low part will be added to y in a fresh attempt.
    END DO

    !local_sum = SUM( arr_local )
    CALL MPI_Allreduce(local_sum,global_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  END FUNCTION get_global_sum_cmplx_1d


  SUBROUTINE gather_all_from_lo_to_rank0(self,arr_lo,arr_full)
    type(gather_type) :: self
    COMPLEX(kind=8) :: arr_lo(llm:ulm,self%dim2)
    COMPLEX(kind=8) :: arr_full(1:lm_max,self%dim2)
    
    INTEGER :: ierr,irank,l,m,nR
    !COMPLEX(kind=8) :: temp_lo((1:lm_max,self%dim2)
    COMPLEX(kind=8), allocatable :: temp_lo(:,:)
    INTEGER :: type_size,gather_tag,status(MPI_STATUS_SIZE)

    IF (rank == 0) ALLOCATE(temp_lo(1:lm_max,self%dim2))
    IF (n_procs == 1) THEN
       ! copy the data on rank 0
       DO nR=1,self%dim2
          temp_lo(llm:ulm,nR)=arr_lo(:,nR)
       END DO
    ELSE
       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       gather_tag=1990
       IF (rank == 0) THEN
          DO irank=1,n_procs-1
             CALL MPI_Recv(temp_lo(lmStartB(irank+1),1),1,&
                  &self%gather_mpi_type(irank),irank,gather_tag,&
                  & MPI_COMM_WORLD,status,ierr)
          END DO

          ! copy the data on rank 0
          DO nR=1,self%dim2
             temp_lo(llm:ulm,nR)=arr_lo(:,nR)
          END DO
          !WRITE(*,"(A,I3,A,2ES22.14)") "recving temp_lo(",1+irank*lm_per_rank,") = ",&
          !     & SUM(temp_lo(1+irank*lm_per_rank:1+irank*lm_per_rank+lm_on_last_rank,:))
       ELSE
          ! Now send the data to rank 0
          !WRITE(*,"(A,I5,A,I2)") "Sending ",(ulm-llm+1)*self%dim2," dc from rank ",rank
          !WRITE(*,"(A,2ES22.14)") "sending arr_lo = ", SUM(arr_lo)
          CALL MPI_Send(arr_lo,self%dim2*(ulm-llm+1),MPI_DOUBLE_COMPLEX,0,gather_tag,&
               &MPI_COMM_WORLD,ierr)
       END IF
       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    END IF

    IF (rank == 0) THEN
       ! reorder
       do nR=1,self%dim2
          DO l=0,l_max
             DO m=0,l,minc
                arr_full(st_map%lm2(l,m),nR) = temp_lo(lo_map%lm2(l,m),nR)
             END DO
          END DO
       end do
       DEALLOCATE(temp_lo)
    END IF

  END SUBROUTINE gather_all_from_lo_to_rank0

  SUBROUTINE create_gather_type(self,dim2)
    type(gather_type) :: self
    integer :: dim2

    integer :: proc

    ! Define the datatypes for gather_all_from_lo_to_rank0
    ! the sending array has dimension (llm:ulm,1:dim2)
    ! receiving array has dimension (1:lm_max,1:dim2)

    ALLOCATE(self%gather_mpi_type(0:n_procs-1))
    ! 1. Datatype for the data on one rank 
    DO proc=0,n_procs-1
       CALL MPI_type_vector(dim2,lmStopB(proc+1)-lmStartB(proc+1)+1,&
            & lm_max,MPI_DOUBLE_COMPLEX,&
            & self%gather_mpi_type(proc),ierr)
       CALL MPI_Type_commit(self%gather_mpi_type(proc),ierr)
    END DO
    ! 2. Datatype for the data on the last rank
    !CALL MPI_Type_vector(dim2,lmStopB(n_procs)-lmStartB(n_procs)+1,&
    !     &lm_max,MPI_DOUBLE_COMPLEX,&
    !     & self%gather_mpi_type_end,ierr)
    !CALL MPI_Type_commit(self%gather_mpi_type_end,ierr)
    self%dim2=dim2
  END SUBROUTINE create_gather_type

  SUBROUTINE destroy_gather_type(self)
    type(gather_type) :: self

    integer :: proc

    DO proc=0,n_procs-1
       CALL MPI_Type_free(self%gather_mpi_type(proc),ierr)
    END DO
    DEALLOCATE(self%gather_mpi_type)
  END SUBROUTINE destroy_gather_type


  SUBROUTINE gather_from_lo_to_rank0(arr_lo,arr_full)
    COMPLEX(kind=8) :: arr_lo(llm:ulm)
    COMPLEX(kind=8) :: arr_full(1:lm_max)

    INTEGER :: sendcounts(0:n_procs-1),displs(0:n_procs-1)
    INTEGER :: irank,l,m
    !COMPLEX(kind=8) :: temp_lo(1:lm_max)

    DO irank=0,n_procs-1
       sendcounts(irank) = lmStopB(irank+1)-lmStartB(irank+1)+1
       displs(irank) = lmStartB(irank+1)-1 !irank*lm_per_rank
    END DO
    !sendcounts(n_procs-1) = lm_on_last_rank
    
    CALL MPI_GatherV(arr_lo,sendcounts(rank),MPI_DOUBLE_COMPLEX,&
         &temp_gather_lo,sendcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

    IF (rank == 0) THEN
       ! reorder
       DO l=0,l_max
          DO m=0,l,minc
             arr_full(st_map%lm2(l,m)) = temp_gather_lo(lo_map%lm2(l,m))
          END DO
       END DO
    END IF
    
  END SUBROUTINE gather_from_lo_to_rank0

  subroutine scatter_from_rank0_to_lo(arr_full,arr_lo)
    COMPLEX(kind=8) :: arr_full(1:lm_max)
    COMPLEX(kind=8) :: arr_lo(llm:ulm)

    INTEGER :: sendcounts(0:n_procs-1),displs(0:n_procs-1)
    INTEGER :: irank,l,m
    !COMPLEX(kind=8) :: temp_lo(1:lm_max)

    DO irank=0,n_procs-1
       sendcounts(irank) = lmStopB(irank+1)-lmStartB(irank+1)+1
       displs(irank) = lmStartB(irank+1)-1
    END DO

    if (rank == 0) then
       ! reorder
       DO l=0,l_max
          DO m=0,l,minc
             temp_gather_lo(lo_map%lm2(l,m)) = arr_full(st_map%lm2(l,m))
          END DO
       END DO
    end if

    CALL MPI_ScatterV(temp_gather_lo,sendcounts,displs,MPI_DOUBLE_COMPLEX,&
         & arr_lo,sendcounts(rank),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

  end subroutine scatter_from_rank0_to_lo


  SUBROUTINE create_lm2r_type(self,count)
    type(lm2r_type) :: self
    INTEGER,OPTIONAL,INTENT(IN) :: count

    IF (.NOT.PRESENT(count)) THEN
       self%count=1
    ELSE
       self%count = count
    END IF
    ALLOCATE(self%s_request(n_procs-1))
    ALLOCATE(self%r_request(n_procs-1))
    ALLOCATE(self%final_wait_array(2*(n_procs-1)))
    ALLOCATE(self%temp_Rloc(1:lm_max,nRstart:nRstop,1:self%count))
  END SUBROUTINE create_lm2r_type

  SUBROUTINE destroy_lm2r_type(self)
    type(lm2r_type) :: self

    deallocate(self%temp_Rloc)
    DEALLOCATE(self%s_request)
    DEALLOCATE(self%r_request)
    DEALLOCATE(self%final_wait_array)
  END SUBROUTINE destroy_lm2r_type

  ! --------------------- NONBLOCKING ---------------------
  ! Here comes the nonblocking variant
  SUBROUTINE lm2r_redist_start(self,arr_LMloc,arr_Rloc)
    type(lm2r_type) :: self
    COMPLEX(kind=8),INTENT(IN)  :: arr_LMloc(llm:ulm,n_r_max,*)
    COMPLEX(kind=8),INTENT(OUT) :: arr_Rloc(lm_max,nRstart:nRstop,*)

    ! Local variables
    INTEGER :: send_pe, recv_pe,i,irank
    integer :: transfer_tag=1111

    !PERFON('lm2r_st')

    IF (rank < n_procs-1) THEN
       ! all the ranks from [0,n_procs-2]
       DO irank=0,n_procs-1
          !IF (rank == irank) THEN
          ! just copy
          !   arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
          !ELSE
          ! send_pe: send to this rank
          ! recv_pe: receive from this rank
          send_pe = MODULO(rank+irank,n_procs)
          recv_pe = MODULO(rank-irank+n_procs,n_procs)
          !PRINT*,"send to ",send_pe,",     recv from ",recv_pe
          IF (rank == send_pe) THEN
             !PERFON('loc_copy')
             ! just copy
             DO i=1,self%count
                arr_Rloc(llm:ulm,nRstart:nRstop,i)=arr_LMloc(llm:ulm,nRstart:nRstop,i)
             END DO
             !PERFOFF
          ELSE
             !PERFON('irecv')
             !IF (recv_pe == n_procs-1) THEN
             !   CALL MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),nRstart,1),&
             !        & 1,s_transfer_type_cont(n_procs,self%count),recv_pe,transfer_tag,&
             !        &MPI_COMM_WORLD,self%r_request(irank),ierr)
             !ELSE
                CALL MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),nRstart,1),&
                     & 1,s_transfer_type_cont(recv_pe+1,self%count),recv_pe,transfer_tag,&
                     &MPI_COMM_WORLD,self%r_request(irank),ierr)
             !END IF
             !PERFOFF
             !PERFON('isend')
             IF (send_pe == n_procs-1) THEN
                CALL MPI_Isend(arr_LMloc(llm,1+nr_per_rank*send_pe,1),&
                     & 1,r_transfer_type_nr_end_cont(rank+1,self%count),send_pe,transfer_tag,&
                     &MPI_COMM_WORLD,self%s_request(irank),ierr)
             ELSE
                CALL MPI_Isend(arr_LMloc(llm,1+nr_per_rank*send_pe,1),&
                     & 1,r_transfer_type_cont(rank+1,self%count),send_pe,transfer_tag,&
                     &MPI_COMM_WORLD,self%s_request(irank),ierr)
             END IF
             !PERFOFF
          END IF
       END DO

       i=1
       DO irank=1,n_procs-1
          self%final_wait_array(i)=self%s_request(irank)
          self%final_wait_array(i+1)=self%r_request(irank)
          i = i + 2
       END DO
       !PRINT*,"Waiting for completion of nonblocking communication 1"
       !CALL mpi_waitall(2*n_procs,final_wait_array,array_of_statuses,ierr)
       !PRINT*,"Nonblocking communication 1 is done."
    ELSE
       ! rank  ==  n_procs-1
       ! all receives are with the s_transfer_type_nr_end
       ! all sends are done with r_transfer_type_lm_end
       DO irank=0,n_procs-1
          ! send_pe: send to this rank
          ! recv_pe: receive from this rank
          send_pe = MODULO(rank+irank,n_procs)
          recv_pe = MODULO(rank-irank+n_procs,n_procs)
          !PRINT*,"send to ",send_pe,",     recv from ",recv_pe
          if (rank == send_pe) then
             !PERFON('loc_copy')
             ! just copy
             DO i=1,self%count
                arr_Rloc(llm:ulm,nRstart:nRstop,i)=arr_LMloc(llm:ulm,nRstart:nRstop,i)
             END DO
             !PERFOFF
          ELSE
             !PERFON('irecv')
             CALL MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),nRstart,1),&
                  & 1,s_transfer_type_nr_end_cont(recv_pe+1,self%count),recv_pe,&
                  &transfer_tag,MPI_COMM_WORLD,self%r_request(irank),ierr)
             !PERFOFF
             !PERFON('isend')
             CALL MPI_Isend(arr_LMloc(llm,1+nr_per_rank*send_pe,1),&
                  & 1,r_transfer_type_cont(rank+1,self%count),send_pe,transfer_tag,&
                  &MPI_COMM_WORLD,self%s_request(irank),ierr)
             !PERFOFF
          end if
       END DO
       i=1
       DO irank=1,n_procs-1
          self%final_wait_array(i)=self%s_request(irank)
          self%final_wait_array(i+1)=self%r_request(irank)
          i = i + 2
       END DO
       !PRINT*,"Waiting for completion of nonblocking communication 1"
       !CALL mpi_waitall(2*(n_procs-1),final_wait_array,array_of_statuses,ierr)
       !PRINT*,"Nonblocking communication 1 is done."
    END IF
    !WRITE(*,"(A,I3)") "lm2r_redist_start on n_procs=",n_procs
    !PERFOFF
  END SUBROUTINE lm2r_redist_start

  SUBROUTINE lm2r_redist_wait(self)
    type(lm2r_type) :: self
    integer :: ierr

    integer :: array_of_statuses(MPI_STATUS_SIZE,2*n_procs)

    !PERFON('lm2r_wt')
    !WRITE(*,"(A,I3)") "n_procs = ",n_procs
    !WRITE(*,"(2(A,I3))") "Waiting for ",2*(n_procs-1)," requests,", size(self%final_wait_array)
    CALL mpi_waitall(2*(n_procs-1),self%final_wait_array,array_of_statuses,ierr)
    !PERFOFF
  END SUBROUTINE lm2r_redist_wait

  SUBROUTINE lo2r_redist_start(self,arr_lo,arr_Rloc)
    type(lm2r_type) :: self
    COMPLEX(kind=8), INTENT(IN) :: arr_lo(llm:ulm,1:n_r_max,*)
    COMPLEX(kind=8), TARGET, INTENT(OUT) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)


    PERFON('lo2r_st')
    !CALL lm2r_redist(arr_lo,temp_lo)
    self%arr_Rloc(1:,nRstart:,1:) => arr_Rloc(1:lm_max,nRstart:nRstop,1:self%count)
    !self%arr_Rloc(1:,nRstart:) => arr_Rloc(1:,nRstart:)
    CALL lm2r_redist_start(self,arr_lo,self%temp_Rloc)

    PERFOFF
  END SUBROUTINE lo2r_redist_start

  SUBROUTINE lo2r_redist_wait(self)
    type(lm2r_type) :: self

    ! Local variables
    INTEGER :: nR,l,m,i

    !PERFON("lo2r_wt")
    call lm2r_redist_wait(self)
    ! now in self%temp_Rloc we do have the lo_ordered r-local part
    ! now reorder to the original ordering
    DO i=1,self%count
       DO nR=nRstart,nRstop
          DO l=0,l_max
             DO m=0,l,minc
                self%arr_Rloc(st_map%lm2(l,m),nR,i) = self%temp_Rloc(lo_map%lm2(l,m),nR,i)
             END DO
          END DO
       END DO
    END DO
    !PERFOFF
  END SUBROUTINE lo2r_redist_wait

  SUBROUTINE r2lm_redist(arr_rloc,arr_LMloc)
    COMPLEX(kind=8),INTENT(IN) :: arr_Rloc(lm_max,nRstart:nRstop)
    COMPLEX(kind=8),INTENT(OUT) :: arr_LMloc(llm:ulm,n_r_max)

    ! Local variables
    INTEGER :: send_pe, recv_pe,i,irank
    integer :: transfer_tag=1111
    LOGICAL :: yetComplete(2*(n_procs-1))
    logical :: flag
    INTEGER :: status(MPI_STATUS_SIZE)
    integer :: completeCounter

    !WRITE(*,"(A)") "----------- start r2lm_redist -------------"
    !PERFON('r2lm_dst')
    IF (rank < n_procs-1) THEN
       ! all the ranks from [0,n_procs-2]
       DO irank=0,n_procs-1
          !IF (rank == irank) THEN
          ! just copy
          !   arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
          !ELSE
          ! send_pe: send to this rank
          ! recv_pe: receive from this rank
          send_pe = MODULO(rank+irank,n_procs)
          recv_pe = MODULO(rank-irank+n_procs,n_procs)
          IF (rank == send_pe) THEN
             arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
          ELSE
             CALL MPI_Isend(arr_Rloc(lmStartB(send_pe+1),nRstart),&
                  & 1,s_transfer_type(send_pe+1),send_pe,transfer_tag,&
                  &MPI_COMM_WORLD,s_request(irank),ierr)
             !WRITE(*,"(2(A,I3))") "Sending s_transfer_type(",send_pe+1,") to pe ",send_pe
             IF (recv_pe == n_procs-1) THEN
                CALL MPI_Irecv(arr_LMloc(llm,1+nr_per_rank*recv_pe),&
                     & 1,r_transfer_type_nr_end(rank+1),recv_pe,transfer_tag,&
                     &MPI_COMM_WORLD,r_request(irank),ierr)
                !WRITE(*,"(2(A,I3))") "Receiving r_transfer_type_nr_end(",rank+1,&
                !     &") from pe ",recv_pe
             ELSE
                CALL MPI_Irecv(arr_LMloc(llm,1+nr_per_rank*recv_pe),&
                     & 1,r_transfer_type(rank+1),recv_pe,transfer_tag,&
                     &MPI_COMM_WORLD,r_request(irank),ierr)
                !WRITE(*,"(2(A,I3))") "Receiving r_transfer_type(",rank+1,") from pe ",recv_pe
             END IF
          END IF
       END DO

       i=1
       DO irank=1,n_procs-1
          final_wait_array(i)=s_request(irank)
          final_wait_array(i+1)=r_request(irank)
          i = i + 2
       END DO
       !PRINT*,"Waiting for completion of nonblocking communication 1"
#if 0
       yetComplete=.FALSE.
       completeCounter=0
       i=1
       DO WHILE (completeCounter < 2*(n_procs-1))
          IF (.NOT.yetComplete(i)) CALL mpi_test(final_wait_array(i),flag,status,ierr)
          IF (flag) THEN
             yetComplete(i)=.true.
             !WRITE(*,"(A,I3,A,I10)") "status of request ",i," is ",status(MPI_ERROR)
             completeCounter = completeCounter + 1
          END IF
          i = MODULO(i,2*(n_procs-1))+1
       END DO
#endif
       CALL mpi_waitall(2*(n_procs-1),final_wait_array,array_of_statuses,ierr)
       IF (ierr /= MPI_SUCCESS) WRITE(*,"(A)") "Error with nonblocking comm. 1"
       !PRINT*,"Nonblocking communication 1 is done."
    ELSE
       ! rank  ==  n_procs-1
       ! all receives are with the r_transfer_type_lm_end
       ! all sends are done with s_transfer_type_nr_end
       DO irank=0,n_procs-1
          ! send_pe: send to this rank
          ! recv_pe: receive from this rank
          send_pe = MODULO(rank+irank,n_procs)
          recv_pe = MODULO(rank-irank+n_procs,n_procs)
          !PRINT*,"send to ",send_pe,",     recv from ",recv_pe
          IF (rank == send_pe) THEN
             ! just copy
             arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
          ELSE
             CALL MPI_Irecv(arr_LMloc(llm,1+nr_per_rank*recv_pe),&
                  & 1,r_transfer_type(rank+1),recv_pe,transfer_tag,&
                  &MPI_COMM_WORLD,r_request(irank),ierr)
             !WRITE(*,"(2(A,I3))") "Receiving r_transfer_type(",rank+1,") from pe ",recv_pe
             CALL MPI_Isend(arr_Rloc(lmStartB(send_pe+1),nRstart),&
                  & 1,s_transfer_type_nr_end(send_pe+1),send_pe,&
                  &transfer_tag,MPI_COMM_WORLD,s_request(irank),ierr)
             !WRITE(*,"(2(A,I3))") "Sending s_transfer_type_nr_end(",send_pe+1,") to pe ",send_pe

          END IF
       END DO
       i=1
       DO irank=1,n_procs-1
          final_wait_array(i)=s_request(irank)
          final_wait_array(i+1)=r_request(irank)
          i = i + 2
       END DO
       !PRINT*,"Waiting for completion of nonblocking communication 2"
       CALL mpi_waitall(2*(n_procs-1),final_wait_array,array_of_statuses,ierr)
       IF (ierr /= MPI_SUCCESS) WRITE(*,"(A)") "Error with nonblocking comm. 2"
       !PRINT*,"Nonblocking communication 2 is done."
    end if

    !PERFOFF
    !WRITE(*,"(A)") "----------- end   r2lm_redist -------------"

  END SUBROUTINE r2lm_redist


  SUBROUTINE r2lo_redist(arr_Rloc,arr_lo)
    COMPLEX(kind=8), INTENT(IN) :: arr_Rloc(1:lm_max,nRstart:nRstop)
    COMPLEX(kind=8), INTENT(OUT) :: arr_lo(llm:ulm,1:n_r_max)

    ! Local variables
    INTEGER :: nR,l,m
    ! temporary reordered array
    !COMPLEX(kind=8) :: temp_lo(1:lm_max,nRstart:nRstop)

    ! Just copy the array with permutation
    !PERFON('r2lo_dst')
    DO nR=nRstart,nRstop
       DO l=0,l_max
          DO m=0,l,minc
             temp_r2lo(lo_map%lm2(l,m),nR) = arr_Rloc(st_map%lm2(l,m),nR)
          END DO
       END DO
    END DO

    CALL r2lm_redist(temp_r2lo,arr_lo)
    !PERFOFF
  END SUBROUTINE r2lo_redist


  SUBROUTINE lm2lo_redist(arr_LMloc,arr_lo)
    COMPLEX(kind=8), INTENT(IN) :: arr_LMloc(llm:ulm,1:n_r_max)
    COMPLEX(kind=8), INTENT(OUT) :: arr_lo(llm:ulm,1:n_r_max)

    ! Local variables
    INTEGER :: nR,l,m

    IF (n_procs > 1) THEN
       PRINT*,"lm2lo not yet parallelized"
       stop
    END IF

    DO nR=1,n_r_max
       DO l=0,l_max
          DO m=0,l,minc
             arr_lo(lo_map%lm2(l,m),nR) = arr_LMloc(st_map%lm2(l,m),nR)
          END DO
       END DO
    END DO
    
  END SUBROUTINE lm2lo_redist

  SUBROUTINE lo2lm_redist(arr_lo,arr_LMloc)
    COMPLEX(kind=8), INTENT(IN) :: arr_lo(llm:ulm,1:n_r_max)
    COMPLEX(kind=8), INTENT(OUT) :: arr_LMloc(llm:ulm,1:n_r_max)

    ! Local variables
    INTEGER :: nR,l,m

    IF (n_procs > 1) THEN
       PRINT*,"lo2lm not yet parallelized"
       stop
    END IF

    DO nR=1,n_r_max
       DO l=0,l_max
          DO m=0,l,minc
             arr_LMloc(st_map%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
          END DO
       END DO
    END DO
    
  END SUBROUTINE lo2lm_redist

#ifdef WITH_MPI
#define STANDARD 1001
#define EXTENDED 1002
#define DATATYPE 1003

#define ALLGATHER STANDARD
#if (ALLGATHER==STANDARD)
  subroutine myAllGather(arr,dim1,dim2)

    use blocking
    use parallel_mod
    implicit none

    integer, intent(IN) :: dim1,dim2
    complex(kind=8), intent(INOUT) :: arr(dim1,dim2)

    integer :: lmStart_on_rank,lmStop_on_rank
    integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
    integer :: irank,nR
    !double precision :: local_sum, global_sum, recvd_sum

    !write(*,"(A,ES15.8)") "before: arr = ",sum(real(conjg(arr)*arr))

    lmStart_on_rank = lmStartB(1+rank*nLMBs_per_rank)
    lmStop_on_rank  = lmStopB(MIN((rank+1)*nLMBs_per_rank,nLMBs)) 
    sendcount  = lmStop_on_rank-lmStart_on_rank+1
    do irank=0,n_procs-1
       recvcounts(irank) = lmStopB( MIN((irank+1)*nLMBs_per_rank,nLMBs) ) - lmStartB( 1+irank*nLMBs_per_rank ) + 1
    end do
    DO irank=0,n_procs-1
       !displs(irank) = (nR-1)*dim1 + sum(recvcounts(0:irank-1))
       displs(irank) = SUM(recvcounts(0:irank-1))
    END DO
    !write(*,"(4X,2I4,A,I6)") rank,nR," displs = ",displs(rank)
    !write(*,"(5(A,I4))") "LMBlocks ",1+rank*nLMBs_per_rank,"->",(rank+1)*nLMBs_per_rank,&
    !     &", lm runs from ",lmStart_on_rank," to ",lmStop_on_rank,", recvcounts = ",recvcounts(rank)
    do nR=1,dim2
       !local_sum = sum( real( conjg(arr(lmStart_on_rank:lmStop_on_rank,nR))*&
       !     &           arr(lmStart_on_rank:lmStop_on_rank,nR) ) )
       CALL MPI_AllGatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
            & arr(1,nR),recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
       !recvd_sum = sum( real( &
       !     & conjg(arr( lmStartB(1+(modulo(rank+1,2)*nLMBs_per_rank)):lmStopB(1+(modulo(rank+1,2)+1)*nLMBs_per_rank-1),nR ))*&
       !     &       arr( lmStartB(1+(modulo(rank+1,2)*nLMBs_per_rank)):lmStopB(1+(modulo(rank+1,2)+1)*nLMBs_per_rank-1),nR ) ) )
       !global_sum = sum( real( conjg(arr(:,nR))*arr(:,nR) ) )
       !write(*,"(4X,A,I4,3(A,ES20.13))") "nR = ",nR,": l_sum = ",local_sum,", r_sum = ",recvd_sum,", g_sum = ", global_sum
    end do
    !write(*,"(A)") "---------------------------"
  end subroutine myAllGather

#elif (ALLGATHER==EXTENDED)
  SUBROUTINE myAllGather(arr,edim1,dim2)

    use blocking
    use parallel_mod
    implicit none

    integer, intent(IN) :: edim1,dim2
    complex(kind=8), intent(INOUT) :: arr(edim1, dim2)

    integer :: sendcount,recvcount
    integer :: irank,nR

    recvcount = edim1/n_procs
    do nR=1,dim2
       CALL MPI_AllGather(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
            & arr(1,nR),recvcount,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
    end do
  END SUBROUTINE myAllGather

#elif (ALLGATHER==DATATYPE)
  ! now with a datatype to have one big transfer, not many small ones
  subroutine myAllGather(arr,dim1,dim2)

    use blocking
    use parallel_mod
    implicit none

    integer, intent(IN) :: dim1,dim2
    complex(kind=8), intent(INOUT) :: arr(dim1,dim2)

    integer :: lmStart_on_rank,lmStop_on_rank
    integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
    integer :: irank,nR
    INTEGER :: sendtype, new_sendtype
    INTEGER(KIND=MPI_ADDRESS_KIND) :: lb,extent,extent_dcmplx

    PERFON('mk_dt')
    ! definition of the datatype (will later be pulled out of here)
    ! we assume dim1=lm_max and dim2=n_r_max
    CALL mpi_type_get_extent(MPI_DOUBLE_COMPLEX,lb,extent_dcmplx,ierr)
    CALL mpi_type_vector(dim2,1,dim1,MPI_DOUBLE_COMPLEX,sendtype,ierr)
    lb=0
    extent=extent_dcmplx
    CALL mpi_type_create_resized(sendtype,lb,extent,new_sendtype,ierr)
    CALL mpi_type_commit(new_sendtype,ierr)


    lmStart_on_rank = lmStartB(1+rank*nLMBs_per_rank)
    lmStop_on_rank  = lmStopB(1+(rank+1)*nLMBs_per_rank-1) 
    sendcount  = lmStop_on_rank-lmStart_on_rank+1
    do irank=0,n_procs-1
       recvcounts(irank) = lmStopB( (irank+1)*nLMBs_per_rank ) - lmStartB( 1+irank*nLMBs_per_rank ) + 1
    end do
    DO irank=0,n_procs-1
       !displs(irank) = (nR-1)*dim1 + sum(recvcounts(0:irank-1))
       displs(irank) = SUM(recvcounts(0:irank-1))
    END DO
    PERFOFF
    PERFON('comm')
    CALL MPI_AllGatherV(MPI_IN_PLACE,sendcount,new_sendtype,&
         & arr,recvcounts,displs,new_sendtype,MPI_COMM_WORLD,ierr)
    PERFOFF
  end subroutine myAllGather
#endif
#endif
END MODULE communications
