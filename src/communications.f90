!$Id$
#include "perflib_preproc.cpp"
module communications

#ifdef WITH_MPI
   use mpi
#endif
   use precision_mod, only: cp
   use parallel_mod, only: rank, n_procs, ierr, nr_per_rank, nr_on_last_rank
   use LMLoop_data, only: llm, ulm
   use truncation, only: l_max, lm_max, minc, n_r_max, n_r_ic_max
   use blocking, only: st_map, lo_map, lmStartB, lmStopB
   use radial_data, only: nRstart, nRstop
   use logic, only: l_mag, l_conv, l_heat
 
   implicit none
 
   private 
 
   interface get_global_sum
      module procedure get_global_sum_cmplx_2d, get_global_sum_cmplx_1d, &
                       get_global_sum_real_2d
   end interface
 
   type, public :: lm2r_type
      integer, allocatable :: final_wait_array(:), s_request(:), r_request(:)
      complex(cp), pointer :: temp_Rloc(:,:,:), arr_Rloc(:,:,:)
      integer :: count
   end type lm2r_type
 
   type, public :: gather_type
      integer, allocatable :: gather_mpi_type(:)
      integer :: dim2
   end type gather_type
 
   ! MPI datatypes for the redistribution of the d?dt arrays
   integer, save, allocatable :: s_transfer_type(:),s_transfer_type_nr_end(:)
   integer, save, allocatable :: r_transfer_type(:), r_transfer_type_nr_end(:)
   integer, save, allocatable :: s_transfer_type_cont(:,:)
   integer, save, allocatable :: s_transfer_type_nr_end_cont(:,:)
   integer, save, allocatable :: r_transfer_type_cont(:,:)
   integer, save, allocatable :: r_transfer_type_nr_end_cont(:,:)
   integer :: r_lm_gather_type, r_lm_gather_type_lm_end
   integer, allocatable :: s_request(:),r_request(:),final_wait_array(:)
   integer, allocatable :: array_of_statuses(:,:)
 
   public :: gather_from_lo_to_rank0,scatter_from_rank0_to_lo,&
           & gather_all_from_lo_to_rank0
   public :: get_global_sum, &
           & r2lm_redist,r2lo_redist,initialize_communications,&
           & create_lm2r_type!,lo2r_redist,lm2r_redist
   public :: lo2r_redist_start,lo2r_redist_wait
#ifdef WITH_MPI
   public :: myAllGather
#endif
 
   ! declaration of the types for the redistribution
   !type(lm2r_type),PUBLIC :: lo2r_s, lo2r_ds, lo2r_z, lo2r_dz
   !type(lm2r_type),public :: lo2r_p,lo2r_dp
   !type(lm2r_type), PUBLIC :: lo2r_b, lo2r_db, lo2r_ddb, lo2r_aj, lo2r_dj
   type(lm2r_type), public :: lo2r_w, lo2r_s, lo2r_z, lo2r_p, lo2r_b, lo2r_aj
 
   type(gather_type), public :: gt_OC,gt_IC,gt_cheb
 
   complex(cp), allocatable :: temp_gather_lo(:)
   complex(cp), allocatable :: temp_r2lo(:,:)

contains
  
   subroutine initialize_communications

      integer :: proc,my_lm_per_rank
#ifdef WITH_MPI
      integer(kind=MPI_ADDRESS_KIND) :: zerolb, extent, sizeof_double_complex
      integer(kind=MPI_ADDRESS_KIND) :: lb_marker, myextent, true_lb, true_extent
      integer :: base_col_type,temptype
      integer :: blocklengths(3),blocklengths_on_last(3),displs(3),displs_on_last(3)
      integer :: i

      ! first setup the datatype. It is not equal for all ranks. The n_procs-1 rank can
      ! have a smaller datatype.
      ! Due to the different number of radial and lm points for the ranks, 
      ! we need essentially three different datatypes
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
      allocate(s_transfer_type(n_procs))
      allocate(s_transfer_type_nr_end(n_procs))
      allocate(r_transfer_type(n_procs))
      allocate(r_transfer_type_nr_end(n_procs))

      allocate(s_transfer_type_cont(n_procs,3))
      allocate(s_transfer_type_nr_end_cont(n_procs,3))
      allocate(r_transfer_type_cont(n_procs,3))
      allocate(r_transfer_type_nr_end_cont(n_procs,3))

      do proc=0,n_procs-1
         my_lm_per_rank=lmStopB(proc+1)-lmStartB(proc+1)+1
         !write(*,"(2(A,I4))") "lm_per_rank on rank ", proc," is ",my_lm_per_rank
         call MPI_Type_vector(nr_per_rank,my_lm_per_rank,&
              &lm_max,MPI_DOUBLE_COMPLEX,s_transfer_type(proc+1),ierr)
         call MPI_Type_commit(s_transfer_type(proc+1),ierr)

         ! The same for the last rank for nR
         call MPI_Type_vector(nr_on_last_rank,my_lm_per_rank,&
              &lm_max,MPI_DOUBLE_COMPLEX,s_transfer_type_nr_end(proc+1),ierr)
         call MPI_Type_commit(s_transfer_type_nr_end(proc+1),ierr)

         ! we do not need special receive datatypes, as the buffers are 
         ! contiguous in memory but for ease of reading, we define the 
         ! receive datatypes explicitly
         call MPI_Type_contiguous(my_lm_per_rank*nr_per_rank,&
              & MPI_DOUBLE_COMPLEX,r_transfer_type(proc+1),ierr)
         call MPI_Type_commit(r_transfer_type(proc+1),ierr)
         call MPI_Type_contiguous(my_lm_per_rank*nr_on_last_rank,&
              &MPI_DOUBLE_COMPLEX,r_transfer_type_nr_end(proc+1),ierr)
         call MPI_Type_commit(r_transfer_type_nr_end(proc+1),ierr)


         ! define the transfer types for the containers
         ! same schema as for the other types
         ! some temporary datatypes, not needed for communication
         ! but only for constructing the final datatypes
         call MPI_Type_get_extent(MPI_DOUBLE_COMPLEX,zerolb,sizeof_double_complex,ierr)
         call MPI_Type_contiguous(my_lm_per_rank,MPI_DOUBLE_COMPLEX,temptype,ierr)
         zerolb=0
         extent=lm_max*sizeof_double_complex
         call MPI_Type_create_resized(temptype,zerolb,extent,base_col_type,ierr)
         !call MPI_type_get_extent(base_col_type,lb_marker,myextent,ierr)
         !write(*,"(2(A,I10))") "base_col_type: lb = ",lb_marker,", extent = ",myextent
         blocklengths = (/ nr_per_rank, nr_per_rank, nr_per_rank /)
         displs       = (/ 0,           nr_per_rank, 2*nr_per_rank   /)
         blocklengths_on_last = (/nr_on_last_rank,nr_on_last_rank, nr_on_last_rank/)
         displs_on_last       = (/ 0,             nr_on_last_rank, 2*nr_on_last_rank /)
         do i=1,3
            call MPI_Type_vector(i,nr_per_rank*my_lm_per_rank,n_r_max*my_lm_per_rank,&
                 & MPI_DOUBLE_COMPLEX,r_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_commit(r_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_vector(i,nr_on_last_rank*my_lm_per_rank, &
                 & n_r_max*my_lm_per_rank,MPI_DOUBLE_COMPLEX,      &
                 & r_transfer_type_nr_end_cont(proc+1,i),ierr)
            call MPI_Type_commit(r_transfer_type_nr_end_cont(proc+1,i),ierr)

            call MPI_Type_indexed(i,blocklengths(1:i),&
                 & displs(1:i),base_col_type,s_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_commit(s_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_indexed(i,blocklengths_on_last(1:i),&
                 & displs_on_last(1:i),base_col_type,         &
                 & s_transfer_type_nr_end_cont(proc+1,i),ierr)
            call MPI_Type_commit(s_transfer_type_nr_end_cont(proc+1,i),ierr)

#if 0
            if (i == 3) then
               call MPI_type_get_extent(r_transfer_type_cont(proc+1,i), &
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(r_transfer_type_cont(proc+1,i), &
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                    &
                    & "r_transfer_type_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
               call MPI_type_get_extent(s_transfer_type_cont(proc+1,i), &
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(s_transfer_type_cont(proc+1,i), &
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                    &
                    & "s_transfer_type_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
               
               call MPI_type_get_extent(r_transfer_type_nr_end_cont(proc+1,i), &
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(r_transfer_type_nr_end_cont(proc+1,i),&
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                           &
                    & "r_transfer_type_nr_end_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
               call MPI_type_get_extent(s_transfer_type_nr_end_cont(proc+1,i),&
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(s_transfer_type_nr_end_cont(proc+1,i),&
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                           &
                    & "s_transfer_type_nr_end_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
            end if

#endif
         end do
      end do
#endif


      call create_gather_type(gt_OC,n_r_max)
      call create_gather_type(gt_IC,n_r_ic_max)

#ifdef WITH_MPI
      allocate(s_request(n_procs-1),r_request(n_procs-1))
      allocate(array_of_statuses(MPI_STATUS_SIZE,2*(n_procs-1)))
      allocate(final_wait_array(2*(n_procs-1)))
#endif

      if ( l_heat ) then
         call create_lm2r_type(lo2r_s,2)
         !call create_lm2r_type(lo2r_ds)
      end if
      if ( l_conv ) then
         call create_lm2r_type(lo2r_z,2)
         !call create_lm2r_type(lo2r_dz)
         call create_lm2r_type(lo2r_w,3)
         !call create_lm2r_type(lo2r_w)
         !call create_lm2r_type(lo2r_dw)
         !call create_lm2r_type(lo2r_ddw)
         call create_lm2r_type(lo2r_p,2)
         !call create_lm2r_type(lo2r_dp)
      end if

      if ( l_mag ) then
         call create_lm2r_type(lo2r_b,3)
         !call create_lm2r_type(lo2r_db)
         !call create_lm2r_type(lo2r_ddb)
         call create_lm2r_type(lo2r_aj,2)
         !call create_lm2r_type(lo2r_dj)
      end if


      ! allocate a temporary array for the gather operations.
      allocate(temp_r2lo(lm_max,nRstart:nRstop))
      if ( rank == 0 ) then
         allocate(temp_gather_lo(1:lm_max))
      else
         allocate(temp_gather_lo(1))
      end if

   end subroutine initialize_communications
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_2d(dwdt_local) result(global_sum)

      complex(cp), intent(in) :: dwdt_local(:,:)
      
#ifdef WITH_MPI
      integer :: ierr
      complex(cp) :: local_sum
      
      local_sum = SUM( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DOUBLE_COMPLEX, &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      global_sum= SUM(dwdt_local)
#endif

   end function get_global_sum_cmplx_2d
!-------------------------------------------------------------------------------
   real(cp) function get_global_sum_real_2d(dwdt_local) result(global_sum)

      real(cp), intent(in) :: dwdt_local(:,:)
      
#ifdef WITH_MPI
      integer :: ierr
      real(cp) :: local_sum
      
      local_sum = SUM( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      global_sum= SUM(dwdt_local)
#endif

   end function get_global_sum_real_2d
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_1d(arr_local) result(global_sum)

      complex(cp), intent(in) :: arr_local(:)
      
#ifdef WITH_MPI
      integer :: lb,ub,ierr,i
      complex(cp) :: local_sum,c,y,t

      lb = lbound(arr_local,1)
      ub = ubound(arr_local,1)
      
      ! Kahan summation algorithm
      !function KahanSum(input)
      !var sum = 0.0
      !var c = 0.0          //A running compensation for lost low-order bits.
      !for i = 1 to input.length do
      !    y = input[i] - c    //So far, so good: c is zero.
      !    t = sum + y         //Alas, sum is big, y small, 
      !                        //so low-order digits of y are lost.
      !    c = (t - sum) - y   //(t - sum) recovers the high-order part of y; 
      !                        //subtracting y recovers -(low part of y)
      !    sum = t             //Algebraically, c should always be zero. 
      !                        //Beware eagerly optimising compilers!
      !    //Next time around, the lost low part will be added to y in a fresh attempt.
      !return sum
      
      local_sum = 0.0_cp
      c = 0.0_cp          !A running compensation for lost low-order bits.
      do i=lb,ub
         y = arr_local(i) - c ! So far, so good: c is zero.
         t = local_sum + y    ! Alas, sum is big, y small, 
                              ! so low-order digits of y are lost.
         c = (t - local_sum) - y ! (t - sum) recovers the high-order part of y; 
                                 ! subtracting y recovers -(low part of y)
         local_sum = t           ! Algebraically, c should always be zero. 
                                 ! Beware eagerly optimising compilers
         ! Next time around, the lost low part will be added to y in a fresh attempt.
      end do

      !local_sum = SUM( arr_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DOUBLE_COMPLEX, &
                          MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      global_sum = SUM( arr_local )
#endif

   end function get_global_sum_cmplx_1d
!-------------------------------------------------------------------------------
   subroutine gather_all_from_lo_to_rank0(self,arr_lo,arr_full)

      type(gather_type) :: self
      complex(cp) :: arr_lo(llm:ulm,self%dim2)
      complex(cp) :: arr_full(1:lm_max,self%dim2)
      
      integer :: l,m,nR
#ifdef WITH_MPI
      integer :: ierr,irank
      !complex(cp) :: temp_lo((1:lm_max,self%dim2)
      complex(cp), allocatable :: temp_lo(:,:)
      integer :: type_size,gather_tag,status(MPI_STATUS_SIZE)

      if ( rank == 0 ) allocate(temp_lo(1:lm_max,self%dim2))
      if (n_procs == 1) then
         ! copy the data on rank 0
         do nR=1,self%dim2
            temp_lo(llm:ulm,nR)=arr_lo(:,nR)
         end do
      else
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         gather_tag=1990
         if ( rank == 0 ) then
            do irank=1,n_procs-1
               call MPI_Recv(temp_lo(lmStartB(irank+1),1),1,       &
                    & self%gather_mpi_type(irank),irank,gather_tag,&
                    & MPI_COMM_WORLD,status,ierr)
            end do

            ! copy the data on rank 0
            do nR=1,self%dim2
               temp_lo(llm:ulm,nR)=arr_lo(:,nR)
            end do
            !write(*,"(A,I3,A,2ES22.14)") "recving temp_lo(",1+irank*lm_per_rank,") &
            !     & = ",SUM(temp_lo(1+irank*lm_per_rank:1+irank*lm_per_rank+        &
            !     & lm_on_last_rank,:))
         else
            ! Now send the data to rank 0
            !write(*,"(A,I5,A,I2)") "Sending ",(ulm-llm+1)*self%dim2," &
            !   &    dc from rank ",rank
            !write(*,"(A,2ES22.14)") "sending arr_lo = ", SUM(arr_lo)
            call MPI_Send(arr_lo,self%dim2*(ulm-llm+1),MPI_DOUBLE_COMPLEX, &
                          0,gather_tag,MPI_COMM_WORLD,ierr)
         end if
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end if

      if ( rank == 0 ) then
         ! reorder
         do nR=1,self%dim2
            do l=0,l_max
               do m=0,l,minc
                  arr_full(st_map%lm2(l,m),nR) = temp_lo(lo_map%lm2(l,m),nR)
               end do
            end do
         end do
         deallocate(temp_lo)
      end if
#else
      do nR=1,self%dim2
         do l=0,l_max
            do m=0,l,minc
               arr_full(st_map%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
            end do
         end do
      end do
#endif

   end subroutine gather_all_from_lo_to_rank0
!-------------------------------------------------------------------------------
   subroutine create_gather_type(self,dim2)

      type(gather_type) :: self
      integer :: dim2

      integer :: proc

      ! Define the datatypes for gather_all_from_lo_to_rank0
      ! the sending array has dimension (llm:ulm,1:dim2)
      ! receiving array has dimension (1:lm_max,1:dim2)

#ifdef WITH_MPI
      allocate(self%gather_mpi_type(0:n_procs-1))
      ! 1. Datatype for the data on one rank 
      do proc=0,n_procs-1
         call MPI_type_vector(dim2,lmStopB(proc+1)-lmStartB(proc+1)+1,&
              &               lm_max,MPI_DOUBLE_COMPLEX,              &
              &               self%gather_mpi_type(proc),ierr)
         call MPI_Type_commit(self%gather_mpi_type(proc),ierr)
      end do
#endif
      ! 2. Datatype for the data on the last rank
      !call MPI_Type_vector(dim2,lmStopB(n_procs)-lmStartB(n_procs)+1,&
      !     &lm_max,MPI_DOUBLE_COMPLEX,&
      !     & self%gather_mpi_type_end,ierr)
      !call MPI_Type_commit(self%gather_mpi_type_end,ierr)
      self%dim2=dim2

   end subroutine create_gather_type
!-------------------------------------------------------------------------------
   subroutine destroy_gather_type(self)

      type(gather_type) :: self

      integer :: proc

#ifdef WITH_MPI
      do proc=0,n_procs-1
         call MPI_Type_free(self%gather_mpi_type(proc),ierr)
      end do
#endif
      deallocate(self%gather_mpi_type)

   end subroutine destroy_gather_type
!-------------------------------------------------------------------------------
   subroutine gather_from_lo_to_rank0(arr_lo,arr_full)

      complex(cp) :: arr_lo(llm:ulm)
      complex(cp) :: arr_full(1:lm_max)

      integer :: l,m
#ifdef WITH_MPI
      integer :: sendcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank
      !complex(cp) :: temp_lo(1:lm_max)

      do irank=0,n_procs-1
         sendcounts(irank) = lmStopB(irank+1)-lmStartB(irank+1)+1
         displs(irank) = lmStartB(irank+1)-1 !irank*lm_per_rank
      end do
      !sendcounts(n_procs-1) = lm_on_last_rank
      
      call MPI_GatherV(arr_lo,sendcounts(rank),MPI_DOUBLE_COMPLEX,&
           &           temp_gather_lo,sendcounts,displs,          &
           &           MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

      if ( rank == 0 ) then
         ! reorder
         do l=0,l_max
            do m=0,l,minc
               arr_full(st_map%lm2(l,m)) = temp_gather_lo(lo_map%lm2(l,m))
            end do
         end do
      end if
#else
      do l=0,l_max
         do m=0,l,minc
            arr_full(st_map%lm2(l,m)) = arr_lo(lo_map%lm2(l,m))
         end do
      end do
#endif
    
   end subroutine gather_from_lo_to_rank0
!-------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_lo(arr_full,arr_lo)

      complex(cp) :: arr_full(1:lm_max)
      complex(cp) :: arr_lo(llm:ulm)

      integer :: l,m
#ifdef WITH_MPI
      integer :: sendcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank
      !complex(cp) :: temp_lo(1:lm_max)

      do irank=0,n_procs-1
         sendcounts(irank) = lmStopB(irank+1)-lmStartB(irank+1)+1
         displs(irank) = lmStartB(irank+1)-1
      end do

      if ( rank == 0 ) then
         ! reorder
         do l=0,l_max
            do m=0,l,minc
               temp_gather_lo(lo_map%lm2(l,m)) = arr_full(st_map%lm2(l,m))
            end do
         end do
      end if

      call MPI_ScatterV(temp_gather_lo,sendcounts,displs,MPI_DOUBLE_COMPLEX,&
           &            arr_lo,sendcounts(rank),MPI_DOUBLE_COMPLEX,0,       &
           &            MPI_COMM_WORLD,ierr)
#else
      do l=0,l_max
         do m=0,l,minc
            arr_lo(lo_map%lm2(l,m)) = arr_full(st_map%lm2(l,m))
         end do
      end do
#endif

   end subroutine scatter_from_rank0_to_lo
!-------------------------------------------------------------------------------
   subroutine create_lm2r_type(self,count)

      type(lm2r_type) :: self
      integer, optional, intent(in) :: count

      if (.not. present(count)) then
         self%count=1
      else
         self%count = count
      end if
      allocate(self%s_request(n_procs-1))
      allocate(self%r_request(n_procs-1))
      allocate(self%final_wait_array(2*(n_procs-1)))
      allocate(self%temp_Rloc(1:lm_max,nRstart:nRstop,1:self%count))

   end subroutine create_lm2r_type
!-------------------------------------------------------------------------------
   subroutine destroy_lm2r_type(self)

      type(lm2r_type) :: self

      deallocate(self%temp_Rloc)
      deallocate(self%s_request)
      deallocate(self%r_request)
      deallocate(self%final_wait_array)

   end subroutine destroy_lm2r_type
!-------------------------------------------------------------------------------
  ! --------------------- NONBLOCKING ---------------------
  ! Here comes the nonblocking variant
   subroutine lm2r_redist_start(self,arr_LMloc,arr_Rloc)

      type(lm2r_type) :: self
      complex(cp), intent(in)  :: arr_LMloc(llm:ulm,n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(lm_max,nRstart:nRstop,*)

      integer :: i
#ifdef WITH_MPI
      ! Local variables
      integer :: send_pe,recv_pe,irank
      integer :: transfer_tag=1111

      !PERFON('lm2r_st')

      if ( rank < n_procs-1 ) then
         ! all the ranks from [0,n_procs-2]
         do irank=0,n_procs-1
            !if (rank == irank) then
            ! just copy
            !   arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
            !else
            ! send_pe: send to this rank
            ! recv_pe: receive from this rank
            send_pe = modulo(rank+irank,n_procs)
            recv_pe = modulo(rank-irank+n_procs,n_procs)
            !PRINT*,"send to ",send_pe,",     recv from ",recv_pe
            if (rank == send_pe) then
               !PERFON('loc_copy')
               ! just copy
               do i=1,self%count
                  arr_Rloc(llm:ulm,nRstart:nRstop,i)= &
                       arr_LMloc(llm:ulm,nRstart:nRstop,i)
               end do
               !PERFOFF
            else
               !PERFON('irecv')
               !if (recv_pe == n_procs-1) then
               !   call MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),nRstart,1),   &
               !        &         1,s_transfer_type_cont(n_procs,self%count),&
               !        &         recv_pe,transfer_tag,MPI_COMM_WORLD,       &
               !        &         self%r_request(irank),ierr)
               !else
                  call MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),nRstart,1),     &
                       &         1,s_transfer_type_cont(recv_pe+1,self%count),&
                       &         recv_pe,transfer_tag,MPI_COMM_WORLD,         &
                       &         self%r_request(irank),ierr)
               !end if
               !PERFOFF
               !PERFON('isend')
               if (send_pe == n_procs-1) then
                  call MPI_Isend(arr_LMloc(llm,1+nr_per_rank*send_pe,1),          &
                       &         1,r_transfer_type_nr_end_cont(rank+1,self%count),&
                       &         send_pe,transfer_tag,MPI_COMM_WORLD,             &
                       &         self%s_request(irank),ierr)
               else
                  call MPI_Isend(arr_LMloc(llm,1+nr_per_rank*send_pe,1),   &
                       &         1,r_transfer_type_cont(rank+1,self%count),&
                       &         send_pe,transfer_tag,MPI_COMM_WORLD,      &
                       &         self%s_request(irank),ierr)
               end if
               !PERFOFF
            end if
         end do

         i=1
         do irank=1,n_procs-1
            self%final_wait_array(i)=self%s_request(irank)
            self%final_wait_array(i+1)=self%r_request(irank)
            i = i + 2
         end do
         !PRINT*,"Waiting for completion of nonblocking communication 1"
         !call mpi_waitall(2*n_procs,final_wait_array,array_of_statuses,ierr)
         !PRINT*,"Nonblocking communication 1 is done."
      else
         ! rank  ==  n_procs-1
         ! all receives are with the s_transfer_type_nr_end
         ! all sends are done with r_transfer_type_lm_end
         do irank=0,n_procs-1
            ! send_pe: send to this rank
            ! recv_pe: receive from this rank
            send_pe = modulo(rank+irank,n_procs)
            recv_pe = modulo(rank-irank+n_procs,n_procs)
            !PRINT*,"send to ",send_pe,",     recv from ",recv_pe
            if (rank == send_pe) then
               !PERFON('loc_copy')
               ! just copy
               do i=1,self%count
                  arr_Rloc(llm:ulm,nRstart:nRstop,i)= &
                        arr_LMloc(llm:ulm,nRstart:nRstop,i)
               end do
               !PERFOFF
            else
               !PERFON('irecv')
               call MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),nRstart,1),            &
                    &         1,s_transfer_type_nr_end_cont(recv_pe+1,self%count),&
                    &         recv_pe,transfer_tag,MPI_COMM_WORLD,                &
                    &         self%r_request(irank),ierr)
               !PERFOFF
               !PERFON('isend')
               call MPI_Isend(arr_LMloc(llm,1+nr_per_rank*send_pe,1),    &
                    &         1,r_transfer_type_cont(rank+1,self%count), &
                    &         send_pe,transfer_tag,MPI_COMM_WORLD,       &
                    &         self%s_request(irank),ierr)
               !PERFOFF
            end if
         end do
         i=1
         do irank=1,n_procs-1
            self%final_wait_array(i)=self%s_request(irank)
            self%final_wait_array(i+1)=self%r_request(irank)
            i = i + 2
         end do
         !PRINT*,"Waiting for completion of nonblocking communication 1"
         !call mpi_waitall(2*(n_procs-1),final_wait_array,array_of_statuses,ierr)
         !PRINT*,"Nonblocking communication 1 is done."
      end if
      !write(*,"(A,I3)") "lm2r_redist_start on n_procs=",n_procs
      !PERFOFF


#else
      do i=1,self%count
         arr_Rloc(llm:ulm,nRstart:nRstop,i)= arr_LMloc(llm:ulm,nRstart:nRstop,i)
      end do
#endif

   end subroutine lm2r_redist_start
!-------------------------------------------------------------------------------
   subroutine lm2r_redist_wait(self)

      type(lm2r_type) :: self
#ifdef WITH_MPI
      integer :: ierr
      integer :: array_of_statuses(MPI_STATUS_SIZE,2*n_procs)

      !PERFON('lm2r_wt')
      !write(*,"(A,I3)") "n_procs = ",n_procs
      !write(*,"(2(A,I3))") "Waiting for ",2*(n_procs-1)," requests,", &
      !   &             size(self%final_wait_array)
      call MPI_Waitall(2*(n_procs-1),self%final_wait_array,array_of_statuses,ierr)
      !PERFOFF
#endif

   end subroutine lm2r_redist_wait
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_start(self,arr_lo,arr_Rloc)

      type(lm2r_type) :: self
      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max,*)
      complex(cp), target, intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
  
  
      PERFON('lo2r_st')
      !call lm2r_redist(arr_lo,temp_lo)
      self%arr_Rloc(1:,nRstart:,1:) => arr_Rloc(1:lm_max,nRstart:nRstop,1:self%count)
      !self%arr_Rloc(1:,nRstart:) => arr_Rloc(1:,nRstart:)
      call lm2r_redist_start(self,arr_lo,self%temp_Rloc)
      PERFOFF

   end subroutine lo2r_redist_start
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_wait(self)

      type(lm2r_type) :: self
  
      ! Local variables
      integer :: nR,l,m,i
  
      !PERFON("lo2r_wt")
      call lm2r_redist_wait(self)
      ! now in self%temp_Rloc we do have the lo_ordered r-local part
      ! now reorder to the original ordering
      do i=1,self%count
         do nR=nRstart,nRstop
            do l=0,l_max
               do m=0,l,minc
                  self%arr_Rloc(st_map%lm2(l,m),nR,i) = &
                         self%temp_Rloc(lo_map%lm2(l,m),nR,i)
               end do
            end do
         end do
      end do
      !PERFOFF

   end subroutine lo2r_redist_wait
!-------------------------------------------------------------------------------
   subroutine r2lm_redist(arr_rloc,arr_LMloc)

      complex(cp), intent(in) :: arr_Rloc(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,n_r_max)
  
#ifdef WITH_MPI
      ! Local variables
      integer :: send_pe, recv_pe,i,irank
      integer :: transfer_tag=1111
      logical :: yetComplete(2*(n_procs-1))
      logical :: flag
      integer :: status(MPI_STATUS_SIZE)
      integer :: completeCounter
  
      !write(*,"(A)") "----------- start r2lm_redist -------------"
      !PERFON('r2lm_dst')
      if (rank < n_procs-1) then
         ! all the ranks from [0,n_procs-2]
         do irank=0,n_procs-1
            !if (rank == irank) then
            ! just copy
            !   arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
            !else
            ! send_pe: send to this rank
            ! recv_pe: receive from this rank
            send_pe = modulo(rank+irank,n_procs)
            recv_pe = modulo(rank-irank+n_procs,n_procs)
            if (rank == send_pe) then
               arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
            else
               call MPI_Isend(arr_Rloc(lmStartB(send_pe+1),nRstart),  &
                    &         1,s_transfer_type(send_pe+1),send_pe,   &
                    &         transfer_tag,MPI_COMM_WORLD,            &
                    &         s_request(irank),ierr)
               !write(*,"(2(A,I3))") "Sending s_transfer_type(",send_pe+1,") &
               !&                    to pe ",send_pe
               if (recv_pe == n_procs-1) then
                  call MPI_Irecv(arr_LMloc(llm,1+nr_per_rank*recv_pe),    &
                       &         1,r_transfer_type_nr_end(rank+1),recv_pe,&
                       &         transfer_tag,MPI_COMM_WORLD,             &
                       &         r_request(irank),ierr)
                  !write(*,"(2(A,I3))") "Receiving r_transfer_type_nr_end(",rank+1,&
                  !     &") from pe ",recv_pe
               else
                  call MPI_Irecv(arr_LMloc(llm,1+nr_per_rank*recv_pe),          &
                       &         1,r_transfer_type(rank+1),recv_pe,transfer_tag,&
                       &         MPI_COMM_WORLD,r_request(irank),ierr)
                  !write(*,"(2(A,I3))") "Receiving r_transfer_type(",rank+1,") &
                  !                      from pe ",recv_pe
               end if
            end if
         end do
  
         i=1
         do irank=1,n_procs-1
            final_wait_array(i)=s_request(irank)
            final_wait_array(i+1)=r_request(irank)
            i = i + 2
         end do
         !PRINT*,"Waiting for completion of nonblocking communication 1"
#if 0
         yetComplete=.false.
         completeCounter=0
         i=1
         do while (completeCounter < 2*(n_procs-1))
            if (.not.yetComplete(i)) call mpi_test(final_wait_array(i),flag,status,ierr)
            if (flag) then
               yetComplete(i)=.true.
               !write(*,"(A,I3,A,I10)") "status of request ",i," is ",status(MPI_ERROR)
               completeCounter = completeCounter + 1
            end if
            i = modulo(i,2*(n_procs-1))+1
         end do
#endif
         call mpi_waitall(2*(n_procs-1),final_wait_array,array_of_statuses,ierr)
         if (ierr /= MPI_SUCCESS) write(*,"(A)") "Error with nonblocking comm. 1"
         !PRINT*,"Nonblocking communication 1 is done."
      else
         ! rank  ==  n_procs-1
         ! all receives are with the r_transfer_type_lm_end
         ! all sends are done with s_transfer_type_nr_end
         do irank=0,n_procs-1
            ! send_pe: send to this rank
            ! recv_pe: receive from this rank
            send_pe = modulo(rank+irank,n_procs)
            recv_pe = modulo(rank-irank+n_procs,n_procs)
            !PRINT*,"send to ",send_pe,",     recv from ",recv_pe
            if (rank == send_pe) then
               ! just copy
               arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
            else
               call MPI_Irecv(arr_LMloc(llm,1+nr_per_rank*recv_pe),&
                    &         1,r_transfer_type(rank+1),recv_pe,   &
                    &         transfer_tag,MPI_COMM_WORLD,         &
                    &         r_request(irank),ierr)
               !write(*,"(2(A,I3))") "Receiving r_transfer_type(",rank+1,") &
               !                     from pe ",recv_pe
               call MPI_Isend(arr_Rloc(lmStartB(send_pe+1),nRstart),       &
                    &                  1,s_transfer_type_nr_end(send_pe+1),&
                    &                  send_pe,transfer_tag,MPI_COMM_WORLD,&
                    &                  s_request(irank),ierr)
               !write(*,"(2(A,I3))") "Sending s_transfer_type_nr_end(", &
               !                     send_pe+1,") to pe ",send_pe
  
            end if
         end do
         i=1
         do irank=1,n_procs-1
            final_wait_array(i)=s_request(irank)
            final_wait_array(i+1)=r_request(irank)
            i = i + 2
         end do
         !PRINT*,"Waiting for completion of nonblocking communication 2"
         call mpi_waitall(2*(n_procs-1),final_wait_array,array_of_statuses,ierr)
         if (ierr /= MPI_SUCCESS) write(*,"(A)") "Error with nonblocking comm. 2"
         !PRINT*,"Nonblocking communication 2 is done."
      end if
  
      !PERFOFF
      !write(*,"(A)") "----------- end   r2lm_redist -------------"
#else
      arr_LMLoc(llm:ulm,nRstart:nRstop)=arr_Rloc(llm:ulm,nRstart:nRstop)
#endif

   end subroutine r2lm_redist
!-------------------------------------------------------------------------------
   subroutine r2lo_redist(arr_Rloc,arr_lo)

      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max)
  
      ! Local variables
      integer :: nR,l,m
      ! temporary reordered array
      !complex(cp) :: temp_lo(1:lm_max,nRstart:nRstop)
  
      ! Just copy the array with permutation
      !PERFON('r2lo_dst')
      do nR=nRstart,nRstop
         do l=0,l_max
            do m=0,l,minc
               temp_r2lo(lo_map%lm2(l,m),nR) = arr_Rloc(st_map%lm2(l,m),nR)
            end do
         end do
      end do
  
      call r2lm_redist(temp_r2lo,arr_lo)
      !PERFOFF

   end subroutine r2lo_redist
!-------------------------------------------------------------------------------
   subroutine lm2lo_redist(arr_LMloc,arr_lo)

      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max)
  
      ! Local variables
      integer :: nR,l,m
  
      if (n_procs > 1) then
         PRINT*,"lm2lo not yet parallelized"
         stop
      end if
  
      do nR=1,n_r_max
         do l=0,l_max
            do m=0,l,minc
               arr_lo(lo_map%lm2(l,m),nR) = arr_LMloc(st_map%lm2(l,m),nR)
            end do
         end do
      end do
    
   end subroutine lm2lo_redist
!-------------------------------------------------------------------------------
   subroutine lo2lm_redist(arr_lo,arr_LMloc)

      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max)
  
      ! Local variables
      integer :: nR,l,m
  
      if (n_procs > 1) then
         PRINT*,"lo2lm not yet parallelized"
         stop
      end if
  
      do nR=1,n_r_max
         do l=0,l_max
            do m=0,l,minc
               arr_LMloc(st_map%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
            end do
         end do
      end do
      
   end subroutine lo2lm_redist
!-------------------------------------------------------------------------------
#ifdef WITH_MPI
#define STANDARD 1001
#define EXTENDED 1002
#define DATATYPE 1003

#define ALLGATHER STANDARD
#if (ALLGATHER==STANDARD)
!-------------------------------------------------------------------------------
   subroutine myAllGather(arr,dim1,dim2)

      use blocking
      use parallel_mod
  
      implicit none
  
      integer,     intent(in) :: dim1,dim2
      complex(cp), intent(inout) :: arr(dim1,dim2)
  
      integer :: lmStart_on_rank,lmStop_on_rank
      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank,nR
      !double precision :: local_sum, global_sum, recvd_sum
  
      !write(*,"(A,ES15.8)") "before: arr = ",sum(real(conjg(arr)*arr))
  
      lmStart_on_rank = lmStartB(1+rank*nLMBs_per_rank)
      lmStop_on_rank  = lmStopB(min((rank+1)*nLMBs_per_rank,nLMBs)) 
      sendcount  = lmStop_on_rank-lmStart_on_rank+1
      do irank=0,n_procs-1
         recvcounts(irank) = lmStopB ( min((irank+1)*nLMBs_per_rank,nLMBs) ) &
                           - lmStartB( 1+irank*nLMBs_per_rank ) + 1
      end do
      do irank=0,n_procs-1
         !displs(irank) = (nR-1)*dim1 + sum(recvcounts(0:irank-1))
         displs(irank) = SUM(recvcounts(0:irank-1))
      end do
      !write(*,"(4X,2I4,A,I6)") rank,nR," displs = ",displs(rank)
      !write(*,"(5(A,I4))") "LMBlocks ",1+rank*nLMBs_per_rank,"->", &
      !     &                (rank+1)*nLMBs_per_rank,&
      !     &", lm runs from ",lmStart_on_rank," to ",lmStop_on_rank,", &
      !     &  recvcounts = ",recvcounts(rank)
      do nR=1,dim2
         !local_sum = sum( real( conjg(arr(lmStart_on_rank:lmStop_on_rank,nR))*&
         !     &           arr(lmStart_on_rank:lmStop_on_rank,nR) ) )
         call MPI_AllGatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,     &
              &              arr(1,nR),recvcounts,displs,MPI_DOUBLE_COMPLEX,&
              &              MPI_COMM_WORLD,ierr)
         !recvd_sum = sum( real( &
         !     & conjg(arr( lmStartB(1+&
         !     &(modulo(rank+1,2)*nLMBs_per_rank)):&
         !     &lmStopB(1+(modulo(rank+1,2)+1)*nLMBs_per_rank-1),nR ))*&
         !     &       arr( lmStartB(1+(modulo(rank+1,2)*nLMBs_per_rank)):&
         !     &lmStopB(1+(modulo(rank+1,2)+1)*nLMBs_per_rank-1),nR ) ) )
         !global_sum = sum( real( conjg(arr(:,nR))*arr(:,nR) ) )
         !write(*,"(4X,A,I4,3(A,ES20.13))") "nR = ",nR,": l_sum = ",&
         !&      local_sum,", r_sum = ",recvd_sum,", g_sum = ", global_sum
      end do
      !write(*,"(A)") "---------------------------"

   end subroutine myAllGather
!-------------------------------------------------------------------------------
#elif (ALLGATHER==EXTENDED)
   subroutine myAllGather(arr,edim1,dim2)

      use blocking
      use parallel_mod

      implicit none

      integer,     intent(in) :: edim1,dim2
      complex(cp), intent(inout) :: arr(edim1, dim2)

      integer :: sendcount,recvcount
      integer :: irank,nR

      recvcount = edim1/n_procs
      do nR=1,dim2
         call MPI_AllGather(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
              &             arr(1,nR),recvcount,MPI_DOUBLE_COMPLEX,   &
              &             MPI_COMM_WORLD,ierr)
      end do

   end subroutine myAllGather
!-------------------------------------------------------------------------------
#elif (ALLGATHER==DATATYPE)
   ! now with a datatype to have one big transfer, not many small ones
   subroutine myAllGather(arr,dim1,dim2)

      use blocking
      use parallel_mod

      implicit none

      integer,     intent(in) :: dim1,dim2
      complex(cp), intent(inout) :: arr(dim1,dim2)

      integer :: lmStart_on_rank,lmStop_on_rank
      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank,nR
      integer :: sendtype, new_sendtype
      integer(kind=MPI_ADDRESS_KIND) :: lb,extent,extent_dcmplx

      PERFON('mk_dt')
      ! definition of the datatype (will later be pulled out of here)
      ! we assume dim1=lm_max and dim2=n_r_max
      call mpi_type_get_extent(MPI_DOUBLE_COMPLEX,lb,extent_dcmplx,ierr)
      call mpi_type_vector(dim2,1,dim1,MPI_DOUBLE_COMPLEX,sendtype,ierr)
      lb=0
      extent=extent_dcmplx
      call mpi_type_create_resized(sendtype,lb,extent,new_sendtype,ierr)
      call mpi_type_commit(new_sendtype,ierr)


      lmStart_on_rank = lmStartB(1+rank*nLMBs_per_rank)
      lmStop_on_rank  = lmStopB(1+(rank+1)*nLMBs_per_rank-1) 
      sendcount  = lmStop_on_rank-lmStart_on_rank+1
      do irank=0,n_procs-1
         recvcounts(irank) = lmStopB( (irank+1)*nLMBs_per_rank ) - &
                            lmStartB( 1+irank*nLMBs_per_rank ) + 1
      end do
      do irank=0,n_procs-1
         !displs(irank) = (nR-1)*dim1 + sum(recvcounts(0:irank-1))
         displs(irank) = sum(recvcounts(0:irank-1))
      end do
      PERFOFF
      PERFON('comm')
      call MPI_AllGatherV(MPI_IN_PLACE,sendcount,new_sendtype,&
           &              arr,recvcounts,displs,new_sendtype, &
           &              MPI_COMM_WORLD,ierr)
      PERFOFF

   end subroutine myAllGather
!------------------------------------------------------------------------------
#endif
#endif
!------------------------------------------------------------------------------
end module communications
