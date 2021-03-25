#include "perflib_preproc.cpp"
module mod_mpiatoap_2d
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

   type, extends(type_mpitransp) :: type_mpiatoap_2d
      integer :: ncount
      complex(cp), pointer :: buffer(:) ! actual buffer
      complex(cp), pointer :: rloc_buf(:,:,:) ! points to buffer; used to reorder
      complex(cp), pointer :: lmloc_buf(:,:,:,:)
      logical :: initialized = .false.

   contains
      procedure :: create_comm => create_mlo_atoap_2d
      procedure :: destroy_comm => destroy_mlo_atoap_2d
      procedure :: transp_lm2r => transp_lm2r_atoap_2d_dummy
      procedure :: transp_r2lm => transp_r2lm_atoap_2d_dummy
      procedure :: transp_lm2r_dist => transp_lm2r_atoap_2d_dist
      procedure :: transp_r2lm_dist => transp_r2lm_atoap_2d_dist
      procedure :: reorder_rloc2buffer
      procedure :: reorder_buffer2lmloc
      procedure :: reorder_buffer2rloc
      procedure :: reorder_lmloc2buffer
   end type type_mpiatoap_2d
   
   logical, private :: mlo_atoap_2d_initialized = .false.
   integer, private :: n_atoap_2d_obj = 0
   
   integer, private, allocatable :: lmr2buf(:)
   integer, private, allocatable :: buf2mlo(:)
   integer, private :: rloc_count
   integer, private :: lmloc_count
   
   public :: type_mpiatoap_2d


contains

   !-- Transposition from (n_lm_loc,nRStart) to (n_mlo_loc,n_r_max).
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   subroutine initialize_mlo_atoap_2d
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
      integer :: rloc_count_array(0:n_ranks_r-1)
      integer :: dest(n_lm_loc) ! holds rank_R which will receive the corresponding lm point
      integer :: counter(0:n_ranks_r-1)
      
      if (mlo_atoap_2d_initialized) return
      
      !-- Counts how many (l,m) tuples will be sent to each rank.
      !   The max will be taken afterwards
      rloc_count_array = 0
      do i=1,n_lm_loc
         l = map_dist_st%lm2l(i)
         m = map_dist_st%lm2m(i)
         
         icoord_lo   = map_mlo%ml2coord_lo(m,l)
         dest(i) = icoord_lo
         rloc_count_array(icoord_lo) = rloc_count_array(icoord_lo) + 1
      end do
      rloc_count = maxval(rloc_count_array)
     
! !       ! BUUUUG REPORTTT
! !       rloc_displacements = 0
! !       do icoord_lo=1,n_ranks_r-1
! !          rloc_displacements(icoord_lo) = sum(rloc_count_array(0:icoord_lo-1))
! ! !          print * , "SUM: ", icoord_lo, sum(rloc_count_array(0:icoord_lo-1))
! !       end do
! !       print *, "rloc_displacements: ", rloc_displacements
! !       ! END BUUUUG REPORTTT

      ! Builds the reordering into the rloc_buf
      counter = 1
      allocate(lmr2buf(n_lm_loc))
      allocate(buf2mlo(n_mlo_loc))
      bytes_allocated=bytes_allocated+(n_lm_loc+n_mlo_loc)*SIZEOF_INTEGER
      do i=1,n_lm_loc
         icoord_lo = dest(i)
         lmr2buf(i) = icoord_lo*rloc_count + counter(icoord_lo)
         if (icoord_lo == coord_r) then
            l = map_dist_st%lm2l(i)
            m = map_dist_st%lm2m(i)
            mlo = map_mlo%ml2i(m,l)
            buf2mlo(counter(icoord_lo)) = mlo
         end if
         counter(icoord_lo) = counter(icoord_lo) + 1
      end do
      
      lmloc_count = maxval(dist_r(:,0))
      
      mlo_atoap_2d_initialized = .true.
   end subroutine initialize_mlo_atoap_2d
   
   !----------------------------------------------------------------------------
   subroutine reorder_rloc2buffer(this, arr_Rloc)
   !-- Reorder input variable into send buffer
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap_2d) :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         this%rloc_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i)) = &
            arr_Rloc(i,nRstart:nRstop, 1:this%n_fields)
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2rloc(this, arr_Rloc)
   !-- Reorder receive buffer into output variable
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap_2d) :: this
      complex(cp), intent(out)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         arr_Rloc(i,nRstart:nRstop, 1:this%n_fields) = &
         &  this%rloc_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i))
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2lmloc(this, arr_LMloc)
   !-- Reorder receive buffer onto output variable
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap_2d) :: this
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      integer :: icoord_lo, lr, ur, lb, ub, i, nr
      
      !>@TODO optimize the order of this loop
      do i=1, n_mlo_loc
         do icoord_lo=0,n_ranks_r-1
            nr = dist_r(icoord_lo,0)
            lr = dist_r(icoord_lo,1)
            ur = dist_r(icoord_lo,2)
            arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields) = &
            &  this%lmloc_buf(1:nr, 1:this%n_fields, i, icoord_lo)
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_lmloc2buffer(this, arr_LMloc)
   !-- Reorder input variable into send buffer
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap_2d) :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      integer :: icoord_lo, lr, ur, lb, ub, i, nr
      
      !>@TODO optimize the order of this loop
      do i=1, n_mlo_loc
         do icoord_lo=0,n_ranks_r-1
            nr = dist_r(icoord_lo,0)
            lr = dist_r(icoord_lo,1)
            ur = dist_r(icoord_lo,2)
            this%lmloc_buf(1:nr, 1:this%n_fields, i, icoord_lo) = &
            & arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields)
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine finalize_mlo_atoap_2d
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
      integer :: i, ierr
      
      deallocate(lmr2buf)
      deallocate(buf2mlo)
      mlo_atoap_2d_initialized = .false.
   end subroutine finalize_mlo_atoap_2d
   

   !----------------------------------------------------------------------------
   subroutine create_mlo_atoap_2d(this, n_fields)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpiatoap_2d) :: this
      integer, intent(in) :: n_fields
      
      if (this%initialized) return 
      
      if (.not. mlo_atoap_2d_initialized) call initialize_mlo_atoap_2d
      this%n_fields = n_fields
      this%ncount = rloc_count*lmloc_count*n_fields
      
      n_atoap_2d_obj = n_atoap_2d_obj + 1
      
   end subroutine create_mlo_atoap_2d
   
   !----------------------------------------------------------------------------
   subroutine destroy_mlo_atoap_2d(this)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   !   
      class(type_mpiatoap_2d) :: this
      integer :: i, ierr
      
      if (.not. this%initialized) return
      
      n_atoap_2d_obj = n_atoap_2d_obj - 1
      this%initialized = .false.
      
      if (n_atoap_2d_obj==0) call finalize_mlo_atoap_2d
      
   end subroutine destroy_mlo_atoap_2d
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoap_2d_dist(this, arr_LMloc, arr_Rloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap_2d)     :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      !-- Local variables:
      integer :: ierr
      complex(cp), allocatable, target :: buffer(:)
      
      PERFON('lm2rS')
      allocate(buffer(lmloc_count*rloc_count*n_ranks_r*this%n_fields))
      this%rloc_buf(1:lmloc_count, 1:this%n_fields, 1:rloc_count*n_ranks_r) => buffer
      this%lmloc_buf(1:lmloc_count, 1:this%n_fields, 1:rloc_count, 0:n_ranks_r-1) => buffer
      call this%reorder_lmloc2buffer(arr_LMloc)
      PERFOFF
      PERFON('lm2rW')
      call MPI_Alltoall( MPI_IN_PLACE, 1, 1, buffer, this%ncount,  &
           &             MPI_DEF_COMPLEX, comm_r, ierr)
      PERFOFF
      PERFON('lm2rS')
      call this%reorder_buffer2rloc(arr_Rloc)
      nullify(this%rloc_buf)
      nullify(this%lmloc_buf)
      deallocate(buffer)
      PERFOFF
   end subroutine transp_lm2r_atoap_2d_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_atoap_2d_dist(this, arr_Rloc, arr_LMloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoap_2d)     :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      !-- Local variables:
      integer :: ierr
      complex(cp), allocatable, target :: buffer(:)
      
      
      PERFON('r2lmS')
      allocate(buffer(lmloc_count*rloc_count*n_ranks_r*this%n_fields))
      this%rloc_buf(1:lmloc_count, 1:this%n_fields, 1:rloc_count*n_ranks_r) => buffer
      this%lmloc_buf(1:lmloc_count, 1:this%n_fields, 1:rloc_count, 0:n_ranks_r-1) => buffer
      call this%reorder_rloc2buffer(arr_Rloc)
      PERFOFF
      PERFON('r2lmW')
      call MPI_Alltoall( MPI_IN_PLACE, 1, 1, buffer, this%ncount,  &
           &             MPI_DEF_COMPLEX, comm_r, ierr)
      PERFOFF
      PERFON('r2lmS')
      call this%reorder_buffer2lmloc(arr_LMloc)
      nullify(this%rloc_buf)
      nullify(this%lmloc_buf)
      deallocate(buffer)
      PERFOFF

   end subroutine transp_r2lm_atoap_2d_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoap_2d_dummy(this, arr_LMloc, arr_Rloc)
      class(type_mpiatoap_2d) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      print*, "Dummy transp_lm2r_atoap_2d_dummy, not yet implemented!"
   end subroutine transp_lm2r_atoap_2d_dummy

   subroutine transp_r2lm_atoap_2d_dummy(this, arr_Rloc, arr_LMloc)
      class(type_mpiatoap_2d) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      print*, "Dummy transp_r2lm_atoap_2d_dummy, not yet implemented!"
   end subroutine transp_r2lm_atoap_2d_dummy
   
end module mod_mpiatoap_2d
