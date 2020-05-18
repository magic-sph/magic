module mod_mpiatoap
   ! 
   !-- First alltoall implementation, with zero padding.
   !>@TODO Test alltoall without padding, followed by bcast with "extra"
   !  points from/to ranks which have/need them
   !>@TODO use MPI_TYPES to create a strided array! We will still need the 
   !  reordering into the buffer, but we can force n_fields to remain the last
   !  dimension thus reducing considerably the memory reshuffling
   !>@TODO skip the reordering of the local data into the buffer; we can copy it
   ! into the the output array directly
   !  

   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp, only: type_mpitransp
   use fft
   use LMmapping

   implicit none

   private

   type, extends(type_mpitransp) :: type_mpiatoap
      integer :: ncount
      complex(cp), pointer :: buffer(:) ! actual buffer
      complex(cp), pointer :: send_buf(:,:,:) ! points to buffer; used to reorder
      complex(cp), pointer :: recv_buf(:,:,:,:)
      logical :: initialized = .false.

   contains
      procedure :: create_comm => create_mlo_atoap
      procedure :: destroy_comm => destroy_mlo_atoap
      procedure :: transp_lm2r => transp_lm2r_atoap_dummy
      procedure :: transp_r2lm => transp_r2lm_atoap_dummy
      procedure :: transp_lm2r_dist => transp_lm2r_atoap_dist
      procedure :: transp_r2lm_dist => transp_r2lm_atoap_dist
      procedure :: reorder_rloc2buffer
      procedure :: reorder_buffer2lmloc
      procedure :: reorder_buffer2rloc
      procedure :: reorder_lmloc2buffer
   end type type_mpiatoap
   
   logical, private :: mlo_atoap_initialized = .false.
   integer, private :: n_atoap_obj = 0
   
   integer, private, allocatable :: lmr2buf(:)
   integer, private, allocatable :: buf2mlo(:)
   integer, private :: send_count
   integer, private :: recv_count
   
   public :: type_mpiatoap


contains

   !-- Transposition from (n_lm_loc,nRStart) to (n_mlo_loc,n_r_max).
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   subroutine initialize_mlo_atoap
      !-- Initialize the MPI types for the transposition from ML to Radial
      !   
      !   This transposition assumes that all m's are kept within the comm_r
      !   Thus, comm_mlo is not used, just comm_r instead.
      !   
      !   Author: Rafael Lago (MPCDF) May 2020
      !    
      !-- Local variables
      integer :: i, j, l, m, mlo, icoord_mo, icoord_molo, icoord_lo, ierr
      
      !-- Help variables for computing lmr2buf
      integer :: send_count_array(0:n_ranks_r-1)
      integer :: dest(n_lm_loc) ! holds rank_R which will receive the corresponding lm point
      integer :: counter(0:n_ranks_r-1)
      
      if (mlo_atoap_initialized) return
      
      !-- Counts how many (l,m) tuples will be sent to each rank.
      !   The max will be taken afterwards
      send_count_array = 0
      do i=1,n_lm_loc
         l = map_dist_st%lm2l(i)
         m = map_dist_st%lm2m(i)
         
         icoord_lo   = map_mlo%ml2coord_lo(m,l)
         dest(i) = icoord_lo
         send_count_array(icoord_lo) = send_count_array(icoord_lo) + 1
      end do
      send_count = maxval(send_count_array)
     
! !       ! BUUUUG REPORTTT
! !       send_displacements = 0
! !       do icoord_lo=1,n_ranks_r-1
! !          send_displacements(icoord_lo) = sum(send_count_array(0:icoord_lo-1))
! ! !          print * , "SUM: ", icoord_lo, sum(send_count_array(0:icoord_lo-1))
! !       end do
! !       print *, "send_displacements: ", send_displacements
! !       ! END BUUUUG REPORTTT

      ! Builds the reordering into the send_buf
      counter = 1
      allocate(lmr2buf(n_lm_loc))
      allocate(buf2mlo(n_mlo_loc))
      do i=1,n_lm_loc
         icoord_lo = dest(i)
         lmr2buf(i) = icoord_lo*send_count + counter(icoord_lo)
         if (icoord_lo == coord_r) then
            l = map_dist_st%lm2l(i)
            m = map_dist_st%lm2m(i)
            mlo = map_mlo%ml2i(m,l)
            buf2mlo(counter(icoord_lo)) = mlo
         end if
         counter(icoord_lo) = counter(icoord_lo) + 1
      end do
      
      recv_count = maxval(dist_r(:,0))
      
      mlo_atoap_initialized = .true.
   end subroutine initialize_mlo_atoap
   
   !----------------------------------------------------------------------------
   subroutine reorder_rloc2buffer(this, arr_Rloc)
   !-- Reorder input variable into send buffer
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap) :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         this%send_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i)) = &
            arr_Rloc(i,nRstart:nRstop, 1:this%n_fields)
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2rloc(this, arr_Rloc)
   !-- Reorder receive buffer into output variable
   !   OBS: in the lm->r transposition, the buffer names are flipped.
   !     recv_buffer is actually send_buffer and vice versa.
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap) :: this
      complex(cp), intent(out)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         arr_Rloc(i,nRstart:nRstop, 1:this%n_fields) = &
         &  this%send_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i))
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2lmloc(this, arr_LMloc)
   !-- Reorder receive buffer onto output variable
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap) :: this
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      integer :: icoord_lo, lr, ur, lb, ub, i, nr
      
      !>@TODO optimize the order of this loop
      do i=1, n_mlo_loc
         do icoord_lo=0,n_ranks_r-1
            nr = dist_r(icoord_lo,0)
            lr = dist_r(icoord_lo,1)
            ur = dist_r(icoord_lo,2)
            arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields) = &
            &  this%recv_buf(1:nr, 1:this%n_fields, i, icoord_lo)
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_lmloc2buffer(this, arr_LMloc)
   !-- Reorder input variable into send buffer
   !   OBS: in the lm->r transposition, the buffer names are flipped.
   !     recv_buffer is actually send_buffer and vice versa.
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap) :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      integer :: icoord_lo, lr, ur, lb, ub, i, nr
      
      !>@TODO optimize the order of this loop
      do i=1, n_mlo_loc
         do icoord_lo=0,n_ranks_r-1
            nr = dist_r(icoord_lo,0)
            lr = dist_r(icoord_lo,1)
            ur = dist_r(icoord_lo,2)
            this%recv_buf(1:nr, 1:this%n_fields, i, icoord_lo) = &
            & arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields)
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine finalize_mlo_atoap
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
      integer :: i, ierr
      
      deallocate(lmr2buf)
      deallocate(buf2mlo)
      mlo_atoap_initialized = .false.
   end subroutine finalize_mlo_atoap
   

   !----------------------------------------------------------------------------
   subroutine create_mlo_atoap(this, n_fields)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpiatoap) :: this
      integer, intent(in) :: n_fields
      
      if (this%initialized) return 
      
      if (.not. mlo_atoap_initialized) call initialize_mlo_atoap
      this%n_fields = n_fields
      this%ncount = send_count*recv_count*n_fields
      
      allocate(this%buffer(recv_count*send_count*n_ranks_r*n_fields))
      this%send_buf(1:recv_count, 1:n_fields, 1:send_count*n_ranks_r) => this%buffer
      this%recv_buf(1:recv_count, 1:n_fields, 1:send_count, 0:n_ranks_r-1) => this%buffer
      this%buffer = cmplx(-1.0, -1.0)
      
      n_atoap_obj = n_atoap_obj + 1
      
   end subroutine create_mlo_atoap
   
   !----------------------------------------------------------------------------
   subroutine destroy_mlo_atoap(this)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   !   
      class(type_mpiatoap) :: this
      integer :: i, ierr
      
      if (.not. this%initialized) return
      
      nullify(this%send_buf)
      nullify(this%recv_buf)
      nullify(this%buffer)
      
      n_atoap_obj = n_atoap_obj - 1
      this%initialized = .false.
      
      if (n_atoap_obj==0) call finalize_mlo_atoap
      
   end subroutine destroy_mlo_atoap
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoap_dist(this, arr_LMloc, arr_Rloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap)     :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      !-- Local variables:
      integer :: ierr

      call reorder_lmloc2buffer(this, arr_LMloc)
      call mpi_alltoall(  MPI_IN_PLACE, 1, 1, this%buffer, this%ncount,  &
         &  MPI_DEF_COMPLEX, comm_r, ierr)
      call this%reorder_buffer2rloc(arr_Rloc)
      
   end subroutine transp_lm2r_atoap_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_atoap_dist(this, arr_Rloc, arr_LMloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap)     :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      !-- Local variables:
      integer :: ierr
      
      call this%reorder_rloc2buffer(arr_Rloc)
      call mpi_alltoall(  MPI_IN_PLACE, 1, 1, this%buffer, this%ncount,  &
         &  MPI_DEF_COMPLEX, comm_r, ierr)
      call this%reorder_buffer2lmloc(arr_LMloc)

   end subroutine transp_r2lm_atoap_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoap_dummy(this, arr_LMloc, arr_Rloc)
      class(type_mpiatoap) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      print*, "Dummy transp_lm2r_atoap_dummy, not yet implemented!"
   end subroutine transp_lm2r_atoap_dummy

   subroutine transp_r2lm_atoap_dummy(this, arr_Rloc, arr_LMloc)
      class(type_mpiatoap) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      print*, "Dummy transp_r2lm_atoap_dummy, not yet implemented!"
   end subroutine transp_r2lm_atoap_dummy
   
end module mod_mpiatoap
