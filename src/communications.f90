#include "perflib_preproc.cpp"
module communications

#ifdef WITH_MPI
   use mpimod
#endif
   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use parallel_mod, only: rank, n_procs, ierr
   use truncation, only: l_max, lm_max, minc, n_r_max, n_r_ic_max, l_axi
   use blocking, only: st_map, lo_map, lm_balance, llm, ulm
   use radial_data, only: nRstart, nRstop, radial_balance
   use logic, only: l_mag, l_conv, l_heat, l_chemical_conv, l_finite_diff, &
       &            l_mag_kin, l_double_curl, l_save_out
   use useful, only: abortRun
   use output_data, only: n_log_file, log_file
   use iso_fortran_env, only: output_unit
   use mpi_ptop_mod, only: type_mpiptop
   use mpi_alltoall_mod, only: type_mpiatoav, type_mpiatoaw
   use charmanip, only: capitalize
   use num_param, only: mpi_transp
   use mpi_transp, only: type_mpitransp

   implicit none

   private

   interface get_global_sum
      module procedure get_global_sum_cmplx_2d, get_global_sum_cmplx_1d, &
      &                get_global_sum_real_2d
   end interface

   type, public :: gather_type
      integer, allocatable :: gather_mpi_type(:)
      integer :: dim2
   end type gather_type

   public :: gather_from_lo_to_rank0,scatter_from_rank0_to_lo, &
   &         gather_all_from_lo_to_rank0, gather_from_Rloc
   public :: get_global_sum, finalize_communications, initialize_communications

#ifdef WITH_MPI
   public :: myAllGather
#endif

   interface reduce_radial
      module procedure reduce_radial_1D
      module procedure reduce_radial_2D
   end interface

   public :: reduce_radial, reduce_scalar

   ! declaration of the types for the redistribution
   class(type_mpitransp), public, pointer :: lo2r_s, r2lo_s, lo2r_press
   class(type_mpitransp), public, pointer :: lo2r_flow, r2lo_flow
   class(type_mpitransp), public, pointer :: lo2r_field, r2lo_field
   class(type_mpitransp), public, pointer :: lo2r_xi, r2lo_xi

   type(gather_type), public :: gt_OC,gt_IC,gt_cheb

   complex(cp), allocatable :: temp_gather_lo(:)

contains

   subroutine initialize_communications

      integer(lip) :: local_bytes_used
      integer :: idx

      local_bytes_used=bytes_allocated

      call create_gather_type(gt_OC,n_r_max)
      call create_gather_type(gt_IC,n_r_ic_max)

      call capitalize(mpi_transp)
      if ( index(mpi_transp, 'AUTO') /= 0 ) then
         call find_faster_comm(idx)
      else
         if ( index(mpi_transp, 'P2P') /= 0 .or. index(mpi_transp, 'PTOP') /= 0 &
         &    .or. index(mpi_transp, 'POINTTOPOINT') /= 0  ) then
            idx = 1
         else if ( index(mpi_transp, 'ATOAV') /= 0 .or. index(mpi_transp, 'A2AV') /=0&
         &         .or. index(mpi_transp, 'ALLTOALLV') /= 0 .or. &
         &         index(mpi_transp, 'ALL2ALLV') /= 0 .or. &
         &         index(mpi_transp, 'ALL-TO-ALLV') /= 0 ) then
            idx = 2
         else if ( index(mpi_transp, 'ATOAW') /= 0 .or. index(mpi_transp, 'A2AW') /=0&
         &         .or. index(mpi_transp, 'ALLTOALLW') /= 0 .or. &
         &         index(mpi_transp, 'ALL2ALLW') /= 0 .or. &
         &         index(mpi_transp, 'ALL-TO-ALLW') /= 0 ) then
            idx = 3
         end if
      end if

      !call find_faster_block(idx)

      if ( idx == 1 ) then
         allocate( type_mpiptop :: lo2r_s )
         allocate( type_mpiptop :: r2lo_s )
         allocate( type_mpiptop :: lo2r_flow )
         allocate( type_mpiptop :: r2lo_flow )
         allocate( type_mpiptop :: lo2r_field )
         allocate( type_mpiptop :: r2lo_field )
         allocate( type_mpiptop :: lo2r_xi )
         allocate( type_mpiptop :: r2lo_xi )
         allocate( type_mpiptop :: lo2r_press )
      else if ( idx == 2 ) then
         allocate( type_mpiatoav :: lo2r_s )
         allocate( type_mpiatoav :: r2lo_s )
         allocate( type_mpiatoav :: lo2r_flow )
         allocate( type_mpiatoav :: r2lo_flow )
         allocate( type_mpiatoav :: lo2r_field )
         allocate( type_mpiatoav :: r2lo_field )
         allocate( type_mpiatoav :: lo2r_xi )
         allocate( type_mpiatoav :: r2lo_xi )
         allocate( type_mpiatoav :: lo2r_press )
      else if ( idx == 3 ) then
         allocate( type_mpiatoaw :: lo2r_s )
         allocate( type_mpiatoaw :: r2lo_s )
         allocate( type_mpiatoaw :: lo2r_flow )
         allocate( type_mpiatoaw :: r2lo_flow )
         allocate( type_mpiatoaw :: lo2r_field )
         allocate( type_mpiatoaw :: r2lo_field )
         allocate( type_mpiatoaw :: lo2r_xi )
         allocate( type_mpiatoaw :: r2lo_xi )
         allocate( type_mpiatoaw :: lo2r_press )
      end if

      if ( l_heat ) then
         call lo2r_s%create_comm(2)
         call r2lo_s%create_comm(2)
      end if
      if ( l_chemical_conv ) then
         call lo2r_xi%create_comm(2)
         call r2lo_xi%create_comm(2)
      end if
      if ( l_conv .or. l_mag_kin) then
         call lo2r_flow%create_comm(5)
         call lo2r_press%create_comm(2)
         if ( l_double_curl ) then
            call r2lo_flow%create_comm(4)
         else
            call r2lo_flow%create_comm(3)
         end if
      end if
      if ( l_mag ) then
         call lo2r_field%create_comm(5)
         call r2lo_field%create_comm(3)
      end if

      ! allocate a temporary array for the gather operations.
      if ( rank == 0 ) then
         allocate(temp_gather_lo(1:lm_max))
         bytes_allocated = bytes_allocated + lm_max*SIZEOF_DEF_COMPLEX
      else
         allocate(temp_gather_lo(1))
      end if

      local_bytes_used = bytes_allocated - local_bytes_used
      call memWrite('communications.f90', local_bytes_used)

   end subroutine initialize_communications
!-------------------------------------------------------------------------------
   subroutine finalize_communications

      call destroy_gather_type(gt_OC)
      call destroy_gather_type(gt_IC)

      if ( l_heat ) then
         call lo2r_s%destroy_comm()
         call r2lo_s%destroy_comm()
      end if
      if ( l_chemical_conv ) then
         call lo2r_xi%destroy_comm()
         call r2lo_xi%destroy_comm()
      end if
      if ( l_conv .or. l_mag_kin) then
         call lo2r_flow%destroy_comm()
         call r2lo_flow%destroy_comm()
         call lo2r_press%destroy_comm()
      end if
      if ( l_mag ) then
         call lo2r_field%destroy_comm()
         call r2lo_field%destroy_comm()
      end if

      deallocate( temp_gather_lo )

   end subroutine finalize_communications
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_2d(dwdt_local) result(global_sum)

      complex(cp), intent(in) :: dwdt_local(:,:)

#ifdef WITH_MPI
      integer :: ierr
      complex(cp) :: local_sum

      local_sum = sum( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_COMPLEX, &
           &             MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      global_sum= sum(dwdt_local)
#endif

   end function get_global_sum_cmplx_2d
!-------------------------------------------------------------------------------
   real(cp) function get_global_sum_real_2d(dwdt_local) result(global_sum)

      real(cp), intent(in) :: dwdt_local(:,:)

#ifdef WITH_MPI
      integer :: ierr
      real(cp) :: local_sum

      local_sum = sum( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_REAL,MPI_SUM, &
           &             MPI_COMM_WORLD,ierr)
#else
      global_sum= sum(dwdt_local)
#endif

   end function get_global_sum_real_2d
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_1d(arr_local) result(global_sum)
      !
      ! Kahan summation algorithm
      !
      ! .. code-block:: c
      !
      !    function KahanSum(input)
      !    var sum = 0.0
      !    var c = 0.0             //A running compensation for lost low-order bits.
      !    for i = 1 to input.length do
      !       y = input[i] - c    //So far, so good: c is zero.
      !       t = sum + y         //Alas, sum is big, y small,
      !                           //so low-order digits of y are lost.
      !       c = (t - sum) - y   //(t - sum) recovers the high-order part of y;
      !                           //subtracting y recovers -(low part of y)
      !       sum = t             //Algebraically, c should always be zero.
      !                           //Beware eagerly optimising compilers!
      !       //Next time around, the lost low part will be added to y in a fresh attempt.
      !    return sum
      !
      ! ..
      !

      complex(cp), intent(in) :: arr_local(:)

#ifdef WITH_MPI
      integer :: lb,ub,ierr,i
      complex(cp) :: local_sum,c,y,t

      lb = lbound(arr_local,1)
      ub = ubound(arr_local,1)
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

      !local_sum = sum( arr_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_COMPLEX, &
           &             MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      global_sum = sum( arr_local )
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
      integer :: gather_tag,status(MPI_STATUS_SIZE)

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
               call MPI_Recv(temp_lo(lm_balance(irank)%nStart,1),1,        &
                    &        self%gather_mpi_type(irank),irank,gather_tag, &
                    &        MPI_COMM_WORLD,status,ierr)
            end do

            ! copy the data on rank 0
            do nR=1,self%dim2
               temp_lo(llm:ulm,nR)=arr_lo(:,nR)
            end do
            !write(*,"(A,I3,A,2ES22.14)") "recving temp_lo(",1+irank*lm_per_rank,") &
            !     & = ",sum(temp_lo(1+irank*lm_per_rank:1+irank*lm_per_rank+        &
            !     & lm_on_last_rank,:))
         else
            ! Now send the data to rank 0
            !write(*,"(A,I5,A,I2)") "Sending ",(ulm-llm+1)*self%dim2," &
            !   &    dc from rank ",rank
            !write(*,"(A,2ES22.14)") "sending arr_lo = ", sum(arr_lo)
            call MPI_Send(arr_lo,self%dim2*(ulm-llm+1),MPI_DEF_COMPLEX, &
                 &        0,gather_tag,MPI_COMM_WORLD,ierr)
         end if
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end if

      if ( rank == 0 ) then
         ! reorder
         if ( .not. l_axi ) then
            do nR=1,self%dim2
               do l=0,l_max
                  do m=0,l,minc
                     arr_full(st_map%lm2(l,m),nR) = temp_lo(lo_map%lm2(l,m),nR)
                  end do
               end do
            end do
         else
            do nR=1,self%dim2
               do l=0,l_max
                  arr_full(st_map%lm2(l,0),nR) = temp_lo(lo_map%lm2(l,0),nR)
               end do
            end do
         end if
         deallocate(temp_lo)
      end if
#else
      if ( .not. l_axi ) then
         do nR=1,self%dim2
            do l=0,l_max
               do m=0,l,minc
                  arr_full(st_map%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
               end do
            end do
         end do
      else
         do nR=1,self%dim2
            do l=0,l_max
               arr_full(st_map%lm2(l,0),nR) = arr_lo(lo_map%lm2(l,0),nR)
            end do
         end do
      end if
#endif

   end subroutine gather_all_from_lo_to_rank0
!-------------------------------------------------------------------------------
   subroutine create_gather_type(self,dim2)
      !
      ! Define the datatypes for gather_all_from_lo_to_rank0
      ! the sending array has dimension (llm:ulm,1:dim2)
      ! receiving array has dimension (1:lm_max,1:dim2)
      !

      type(gather_type) :: self
      integer :: dim2

      integer :: proc

#ifdef WITH_MPI
      allocate(self%gather_mpi_type(0:n_procs-1))
      ! 1. Datatype for the data on one rank
      do proc=0,n_procs-1
         call MPI_type_vector(dim2,lm_balance(proc)%n_per_rank,    &
              &               lm_max,MPI_DEF_COMPLEX,              &
              &               self%gather_mpi_type(proc),ierr)
         call MPI_Type_commit(self%gather_mpi_type(proc),ierr)
      end do
#endif

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

      deallocate(self%gather_mpi_type)
#endif

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
         sendcounts(irank) = lm_balance(irank)%n_per_rank
         displs(irank) = lm_balance(irank)%nStart-1 !irank*lm_per_rank
      end do
      !sendcounts(n_procs-1) = lm_on_last_rank

      call MPI_GatherV(arr_lo,sendcounts(rank),MPI_DEF_COMPLEX,  &
           &           temp_gather_lo,sendcounts,displs,         &
           &           MPI_DEF_COMPLEX,0,MPI_COMM_WORLD,ierr)

      if ( rank == 0 ) then
         ! reorder
         if ( .not. l_axi ) then
            do l=0,l_max
               do m=0,l,minc
                  arr_full(st_map%lm2(l,m)) = temp_gather_lo(lo_map%lm2(l,m))
               end do
            end do
         else
            do l=0,l_max
               arr_full(st_map%lm2(l,0)) = temp_gather_lo(lo_map%lm2(l,0))
            end do
         end if
      end if
#else
      if ( .not. l_axi ) then
         do l=0,l_max
            do m=0,l,minc
               arr_full(st_map%lm2(l,m)) = arr_lo(lo_map%lm2(l,m))
            end do
         end do
      else
         do l=0,l_max
            arr_full(st_map%lm2(l,0)) = arr_lo(lo_map%lm2(l,0))
         end do
      end if
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
         sendcounts(irank) = lm_balance(irank)%n_per_rank
         displs(irank) = lm_balance(irank)%nStart-1
      end do

      if ( rank == 0 ) then
         ! reorder
         if ( .not. l_axi ) then
            do l=0,l_max
               do m=0,l,minc
                  temp_gather_lo(lo_map%lm2(l,m)) = arr_full(st_map%lm2(l,m))
               end do
            end do
         else
            do l=0,l_max
               temp_gather_lo(lo_map%lm2(l,0)) = arr_full(st_map%lm2(l,0))
            end do
         end if
      end if

      call MPI_ScatterV(temp_gather_lo,sendcounts,displs,MPI_DEF_COMPLEX,&
           &            arr_lo,sendcounts(rank),MPI_DEF_COMPLEX,0,       &
           &            MPI_COMM_WORLD,ierr)
#else
      if ( .not. l_axi ) then
         do l=0,l_max
            do m=0,l,minc
               arr_lo(lo_map%lm2(l,m)) = arr_full(st_map%lm2(l,m))
            end do
         end do
      else
         do l=0,l_max
            arr_lo(lo_map%lm2(l,0)) = arr_full(st_map%lm2(l,0))
         end do
      end if
#endif

   end subroutine scatter_from_rank0_to_lo
!-------------------------------------------------------------------------------
   subroutine lm2lo_redist(arr_LMloc,arr_lo)

      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max)

      ! Local variables
      integer :: nR,l,m

      if (n_procs > 1) then
         call abortRun('lm2lo not yet parallelized')
      end if

      if ( .not. l_axi ) then
         do nR=1,n_r_max
            do l=0,l_max
               do m=0,l,minc
                  arr_lo(lo_map%lm2(l,m),nR) = arr_LMloc(st_map%lm2(l,m),nR)
               end do
            end do
         end do
      else
         do nR=1,n_r_max
            do l=0,l_max
               arr_lo(lo_map%lm2(l,0),nR) = arr_LMloc(st_map%lm2(l,0),nR)
            end do
         end do
      end if

   end subroutine lm2lo_redist
!-------------------------------------------------------------------------------
   subroutine lo2lm_redist(arr_lo,arr_LMloc)

      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max)

      ! Local variables
      integer :: nR,l,m

      if (n_procs > 1) then
         call abortRun('lo2lm not yet parallelized')
      end if

      if ( .not. l_axi ) then
         do nR=1,n_r_max
            do l=0,l_max
               do m=0,l,minc
                  arr_LMloc(st_map%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
               end do
            end do
         end do
      else
         do nR=1,n_r_max
            do l=0,l_max
               arr_LMloc(st_map%lm2(l,0),nR) = arr_lo(lo_map%lm2(l,0),nR)
            end do
         end do
      end if

   end subroutine lo2lm_redist
!-------------------------------------------------------------------------------
   subroutine gather_from_Rloc(arr_Rloc, arr_glob, irank)
      !
      ! This subroutine gather a r-distributed array on rank=irank
      !

      !-- Input variable
      real(cp), intent(in) :: arr_Rloc(nRstart:nRstop)
      integer,  intent(in) :: irank

      !-- Output variable
      real(cp), intent(out) :: arr_glob(1:n_r_max)

#ifdef WITH_MPI
      !-- Local variables:
      integer :: p
      integer :: scount, rcounts(0:n_procs-1), rdisp(0:n_procs-1)

      scount = nRstop-nRstart+1
      do p=0,n_procs-1
         rcounts(p)=radial_balance(p)%n_per_rank
      end do
      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      call MPI_GatherV(arr_Rloc, scount, MPI_DEF_REAL, arr_glob, rcounts, &
           &           rdisp, MPI_DEF_REAL, irank, MPI_COMM_WORLD, ierr)
#else
      arr_glob(:)=arr_Rloc(:)
#endif

   end subroutine gather_from_Rloc
!-------------------------------------------------------------------------------
   subroutine reduce_radial_2D(arr_dist, arr_glob, irank)

      !-- Input variable
      integer,  intent(in) :: irank
      real(cp), intent(in) :: arr_dist(:,:)

      !-- Output variable
      real(cp), intent(out) :: arr_glob(:,:)

#ifdef WITH_MPI
      call MPI_Reduce(arr_dist, arr_glob, size(arr_dist), MPI_DEF_REAL, &
           &          MPI_SUM, irank, MPI_COMM_WORLD, ierr)
#else
      arr_glob(:,:)=arr_dist(:,:)
#endif

   end subroutine reduce_radial_2D
!-------------------------------------------------------------------------------
   subroutine reduce_radial_1D(arr_dist, arr_glob, irank)

      !-- input variable
      integer,  intent(in) :: irank
      real(cp), intent(in) :: arr_dist(:)

      !-- output variable
      real(cp), intent(inout) :: arr_glob(:)

#ifdef WITH_MPI
      call MPI_Reduce(arr_dist, arr_glob, size(arr_dist), MPI_DEF_REAL, &
           &          MPI_SUM, irank, MPI_COMM_WORLD, ierr)
#else
      arr_glob(:)=arr_dist(:)
#endif

   end subroutine reduce_radial_1D
!-------------------------------------------------------------------------------
   subroutine reduce_scalar(scal_dist, scal_glob, irank)

      !-- input variable
      integer,  intent(in) :: irank
      real(cp), intent(in) :: scal_dist

      !-- output variable
      real(cp), intent(inout) :: scal_glob

#ifdef WITH_MPI
      call MPI_Reduce(scal_dist, scal_glob, 1, MPI_DEF_REAL, &
           &          MPI_SUM, irank, MPI_COMM_WORLD, ierr)
#else
      scal_glob=scal_dist
#endif

   end subroutine reduce_scalar
!-------------------------------------------------------------------------------
   subroutine find_faster_block(idx_type)

      !-- Input variable
      integer, intent(in) :: idx_type

#ifdef WITH_MPI
      !-- Local variables
      complex(cp), allocatable :: arr_Rloc(:,:,:), arr_LMloc(:,:,:)
      class(type_mpitransp), pointer :: lo2r_test
      real(cp) :: rdm_real, rdm_imag, tStart, tStop
      real(cp) :: tBlock(6), tBlock_avg(6)
      integer, parameter :: n_transp=10
      integer :: block_size(6), ind(1)
      integer :: iblock, idx, n_f, n_r, lm, n, n_t, n_out
      character(len=80) :: message

      block_size = [1, 2, 3, 4, 5, 8]

      if ( idx_type == 1 ) then
         allocate( type_mpiptop :: lo2r_test )
      else if ( idx_type == 2 ) then
         allocate( type_mpiatoav :: lo2r_test )
      else if ( idx_type == 3 ) then
         allocate( type_mpiatoaw :: lo2r_test )
      end if
      do iblock=1,size(block_size)

         allocate( arr_RLoc(lm_max,nRstart:nRstop, block_size(iblock)) )
         allocate( arr_LMLoc(llm:ulm,n_r_max, block_size(iblock)) )

         !-- Try the all-to-allv strategy (10 back and forth transposes)
         call lo2r_test%create_comm(block_size(iblock))

         do n_f=1,block_size(iblock)
            do n_r=nRstart,nRstop
               do lm=1,lm_max
                  call random_number(rdm_real)
                  call random_number(rdm_imag)
                  arr_Rloc(lm,n_r,n_f)=cmplx(rdm_real,rdm_imag,kind=cp)
               end do
            end do
         end do

         tStart = MPI_Wtime()
         do n_t=1,n_transp
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call lo2r_test%transp_r2lm(arr_Rloc, arr_LMloc)
            call lo2r_test%transp_lm2r(arr_LMloc, arr_Rloc)
         end do
         tStop = MPI_Wtime()
         tBlock(iblock) = tStop-tStart


         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call lo2r_test%destroy_comm()

         deallocate( arr_Rloc, arr_LMloc )
      end do
      deallocate( lo2r_test )

      !-- Now determine the average over the ranks and send it to rank=0
      call MPI_Reduce(tBlock, tBlock_avg, size(tBlock), MPI_DEF_REAL, &
           &          MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         !-- Average over procs and number of transposes
         tBlock_avg(:) = tBlock_avg(:)/real(n_procs,cp)/real(n_transp,cp)

         tBlock_avg(:) = tBlock_avg(:)/real(block_size(:),cp)

         !-- Determine the fastest
         ind = minloc(tBlock_avg)
         idx = ind(1)

         do n=1,2
            if ( n==1 ) then
               n_out = n_log_file
               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown', &
                  &    position='append')
               end if
            else
               n_out = output_unit
            end if
            write(n_out,*)
            write(n_out,*) '! MPI blocking for the transposes :'
            write(message,'('' ! Container contains 1 field  ='', &
            &               ES10.3, '' s'')') tBlock_avg(1)
            write(n_out,'(A80)') message
            write(message,'('' ! Container contains 2 fields ='', &
            &               ES10.3, '' s'')') tBlock_avg(2)
            write(n_out,'(A80)') message
            write(message,'('' ! Container contains 3 fields ='', &
            &               ES10.3, '' s'')') tBlock_avg(3)
            write(n_out,'(A80)') message
            write(message,'('' ! Container contains 4 fields ='', &
            &               ES10.3, '' s'')') tBlock_avg(4)
            write(n_out,'(A80)') message
            write(message,'('' ! Container contains 5 fields ='', &
            &               ES10.3, '' s'')') tBlock_avg(5)
            write(n_out,'(A80)') message
            write(message,'('' ! Container contains 8 fields ='', &
            &               ES10.3, '' s'')') tBlock_avg(6)
            write(n_out,'(A80)') message
            write(n_out,*)
            if ( n==1 .and. l_save_out ) close(n_log_file)
         end do

      end if
#endif


   end subroutine find_faster_block
!-------------------------------------------------------------------------------
   subroutine find_faster_comm(idx)
      !
      ! This subroutine tests the two MPI transposition strategies and
      ! selects the fastest one.
      !

      !-- Output variable:
      integer, intent(out) :: idx

#ifdef WITH_MPI
      !-- Local variables
      class(type_mpitransp), pointer :: lo2r_test
      complex(cp) :: arr_Rloc(lm_max,nRstart:nRstop,5)
      complex(cp) :: arr_LMloc(llm:ulm,n_r_max,5)
      real(cp) :: tStart, tStop, tAlltoAllv, tPointtoPoint, tAlltoAllw
      real(cp) :: rdm_real, rdm_imag, timers(3)
      integer :: n_f, n_r, lm, n_t, n, n_out, ind(1)
      integer, parameter :: n_transp=10
      character(len=80) :: message

      !-- First fill an array with random numbers
      do n_f=1,5
         do n_r=nRstart,nRstop
            do lm=1,lm_max
               call random_number(rdm_real)
               call random_number(rdm_imag)
               arr_Rloc(lm,n_r,n_f)=cmplx(rdm_real,rdm_imag,kind=cp)
            end do
         end do
      end do

      !-- Try the all-to-allv strategy (10 back and forth transposes)
      allocate( type_mpiatoav :: lo2r_test )
      call lo2r_test%create_comm(5)
      tStart = MPI_Wtime()
      do n_t=1,n_transp
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call lo2r_test%transp_r2lm(arr_Rloc, arr_LMloc)
         call lo2r_test%transp_lm2r(arr_LMloc, arr_Rloc)
      end do
      tStop = MPI_Wtime()
      tAlltoAllv = tStop-tStart
      call lo2r_test%destroy_comm()
      deallocate( lo2r_test)

      !-- Try the point-to-point strategy (10 back and forth transposes)
      allocate( type_mpiptop :: lo2r_test )
      call lo2r_test%create_comm(5)
      tStart = MPI_Wtime()
      do n_t=1,n_transp
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call lo2r_test%transp_r2lm(arr_Rloc, arr_LMloc)
         call lo2r_test%transp_lm2r(arr_LMloc, arr_Rloc)
      end do
      tStop = MPI_Wtime()
      tPointtoPoint = tStop-tStart
      call lo2r_test%destroy_comm()
      deallocate ( lo2r_test )

      !-- Try the all-to-allw strategy (10 back and forth transposes)
      allocate( type_mpiatoaw :: lo2r_test )
      call lo2r_test%create_comm(5)
      tStart = MPI_Wtime()
      do n_t=1,n_transp
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call lo2r_test%transp_r2lm(arr_Rloc, arr_LMloc)
         call lo2r_test%transp_lm2r(arr_LMloc, arr_Rloc)
      end do
      tStop = MPI_Wtime()
      tAlltoAllw = tStop-tStart
      call lo2r_test%destroy_comm()
      deallocate( lo2r_test)

      !-- Now determine the average over the ranks and send it to rank=0
      call MPI_Reduce(tPointtoPoint, timers(1), 1, MPI_DEF_REAL, MPI_SUM, 0, &
           &          MPI_COMM_WORLD, ierr)
      call MPI_Reduce(tAlltoAllv, timers(2), 1, MPI_DEF_REAL, MPI_SUM, 0, &
           &          MPI_COMM_WORLD, ierr)
      call MPI_Reduce(tAlltoAllw, timers(3), 1, MPI_DEF_REAL, MPI_SUM, 0, &
           &          MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         !-- Average over procs and number of transposes
         timers(:) = timers(:)/real(n_procs,cp)/real(n_transp,cp)

         !-- Determine the fastest
         ind = minloc(timers)
         idx = ind(1)

         do n=2,2
            if ( n==1 ) then
               n_out = n_log_file
               if ( l_save_out ) open(newunit=n_log_file, file=log_file, &
               &                      status='unknown', position='append')
            else
               n_out = output_unit
            end if
            write(n_out,*)
            write(n_out,*) '! MPI transpose strategy :'
            write(message,'('' ! isend/irecv/waitall communicator='', &
            &               ES10.3, '' s'')') timers(1)
            write(n_out,'(A80)') message
            write(message,'('' ! alltoallv communicator          ='', &
            &               ES10.3, '' s'')') timers(2)
            write(n_out,'(A80)') message
            write(message,'('' ! alltoallw communicator          ='', &
            &               ES10.3, '' s'')') timers(3)
            write(n_out,'(A80)') message

            if ( idx == 1 ) then
               write(n_out,*) '! -> I choose isend/irecv/waitall'
            else if ( idx == 2 ) then
               write(n_out,*) '! -> I choose alltoallv'
            else if ( idx == 3 ) then
               write(n_out,*) '! -> I choose alltoallw'
            end if
            write(n_out,*)

            if ( n==1 .and. l_save_out ) close(n_log_file)
         end do

      end if

      call MPI_Bcast(idx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
      idx=1  ! In that case it does not matter
#endif

   end subroutine find_faster_comm
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

      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank,nR
      !double precision :: local_sum, global_sum, recvd_sum

      !write(*,"(A,ES15.8)") "before: arr = ",sum(real(conjg(arr)*arr))

      sendcount  = ulm-llm+1
      do irank=0,n_procs-1
         recvcounts(irank) = lm_balance(irank)%n_per_rank
      end do
      do irank=0,n_procs-1
         displs(irank) = sum(recvcounts(0:irank-1))
      end do

      do nR=1,dim2
         call MPI_AllGatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,     &
              &              arr(1,nR),recvcounts,displs,MPI_DEF_COMPLEX,&
              &              MPI_COMM_WORLD,ierr)
      end do

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
         call MPI_AllGather(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
              &             arr(1,nR),recvcount,MPI_DEF_COMPLEX,   &
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

      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank,nR
      integer :: sendtype, new_sendtype
      integer(kind=MPI_ADDRESS_KIND) :: lb,extent,extent_dcmplx

      PERFON('mk_dt')
      ! definition of the datatype (will later be pulled out of here)
      ! we assume dim1=lm_max and dim2=n_r_max
      call MPI_Type_get_extent(MPI_DEF_COMPLEX,lb,extent_dcmplx,ierr)
      call MPI_Type_vector(dim2,1,dim1,MPI_DEF_COMPLEX,sendtype,ierr)
      lb=0
      extent=extent_dcmplx
      call MPI_Type_create_resized(sendtype,lb,extent,new_sendtype,ierr)
      call MPI_Type_commit(new_sendtype,ierr)

      sendcount  = ulm-llm+1
      do irank=0,n_procs-1
         recvcounts(irank) = lm_balance(irank)%n_per_rank
      end do
      do irank=0,n_procs-1
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
