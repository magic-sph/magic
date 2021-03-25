#include "perflib_preproc.cpp"
module mod_mpiatoav_2d
   ! 
   !-- First alltoallv implementation
   !>@TODO skip the reordering of the local data into the buffer; we can copy it
   ! into the the output array directly
   !>@TODO use MPI_TYPES to create a strided array! We will still need the 
   !  reordering into the buffer, but we can force n_fields to remain the last
   !  dimension thus reducing considerably the memory reshuffling
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

   type, extends(type_mpitransp) :: type_mpiatoav_2d
      integer :: ncount
!       complex(cp), pointer :: rloc_buf(:,:,:)
!       complex(cp), pointer :: lmloc_buf(:)
      integer, allocatable :: rloc_count(:)
      integer, allocatable :: lmloc_count(:)
      integer, allocatable :: rloc_disp(:)
      integer, allocatable :: lmloc_disp(:)
      integer :: rloc_sum, lmloc_sum
      logical :: initialized = .false.

   contains
      procedure :: create_comm => create_mlo_atoav_2d
      procedure :: destroy_comm => destroy_mlo_atoav_2d
!       procedure :: transp_lm2r => transp_lm2r_atoav_2d_dummy
!       procedure :: transp_r2lm => transp_r2lm_atoav_2d_dummy
      procedure :: transp_lm2r => transp_lm2r_atoav_2d_dist
      procedure :: transp_r2lm => transp_r2lm_atoav_2d_dist
      procedure :: transp_lm2r_dist => transp_lm2r_atoav_2d_dist
      procedure :: transp_r2lm_dist => transp_r2lm_atoav_2d_dist
   end type type_mpiatoav_2d
   
   logical, private :: mlo_atoav_2d_initialized = .false.
   integer, private :: n_atoav_2d_obj = 0
   
   integer, private, allocatable :: lmr2buf(:)
   integer, private, allocatable :: buf2mlo(:)
   integer, private, allocatable :: send_lm_count(:)
   integer, private, allocatable :: recv_mlo_count(:)
   
   public :: type_mpiatoav_2d


contains

   !-- Transposition from (n_lm_loc,nRStart) to (n_mlo_loc,n_r_max).
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   subroutine initialize_mlo_atoav_2d
      !-- Initialize the MPI types for the transposition from ML to Radial
      !   
      !   This transposition assumes that all m's are kept within the comm_r
      !   Thus, comm_mlo is not used, just comm_r instead.
      !   
      !   Author: Rafael Lago (MPCDF) May 2020
      !    
      !-- Local variables
      integer :: i, j, l, m, mlo, icoord_mo, icoord_molo, icoord_lo, ierr
      integer :: disp
      
      !-- Help variables for computing lmr2buf
      integer :: dest(n_lm_loc) ! holds rank_R which will receive the corresponding lm point
      integer :: counter(0:n_ranks_r-1)
      
      if (mlo_atoav_2d_initialized) return
      
      allocate(send_lm_count(0:n_ranks_lo-1))
      allocate(recv_mlo_count(0:n_ranks_lo-1))
      bytes_allocated = bytes_allocated+2*n_ranks_lo*SIZEOF_INTEGER
      
      ! First we count only (l,m) points here
      send_lm_count = 0
      do i=1,n_lm_loc
         l = map_dist_st%lm2l(i)
         m = map_dist_st%lm2m(i)
         
         icoord_lo   = map_mlo%ml2coord_lo(m,l)
         dest(i) = icoord_lo
         send_lm_count(icoord_lo) = send_lm_count(icoord_lo) + 1
      end do
      
      ! Builds the reordering into the rloc_buf
      counter = 1
      allocate(lmr2buf(n_lm_loc))
      allocate(buf2mlo(n_mlo_loc))
      bytes_allocated = bytes_allocated+(n_lm_loc+n_mlo_loc)*SIZEOF_INTEGER
      do i=1,n_lm_loc
         icoord_lo = dest(i)
         if (icoord_lo>0) then
            lmr2buf(i) = SUM(send_lm_count(0:icoord_lo-1)) + counter(icoord_lo)
         else
            lmr2buf(i) = counter(icoord_lo)
         end if
         if (icoord_lo == coord_r) then
            l = map_dist_st%lm2l(i)
            m = map_dist_st%lm2m(i)
            mlo = map_mlo%ml2i(m,l)
            buf2mlo(counter(icoord_lo)) = mlo
         end if
         counter(icoord_lo) = counter(icoord_lo) + 1
      end do
      
      recv_mlo_count = dist_r(:,0)
      
      ! Some few sanity checks
      ! Once iimplement the bypassing of reodering local data, most of these
      ! checks will be invalid!
      if (minval(buf2mlo)/=1)  then
         write (*,'(A,I0,A)') ' * Wrong buf2mlo! min=', minval(buf2mlo), ', should be 1'
         call abortRun(' * initialize_mlo_atoav_2d: failed sanity check test')
      end if
      if (maxval(buf2mlo)/=n_mlo_loc)  then
         write (*,'(A,I0,A,I0)') ' * Wrong buf2mlo! max=', maxval(buf2mlo), ', should be ',n_mlo_loc
         call abortRun(' * initialize_mlo_atoav_2d: failed sanity check test')
      end if
      if (minval(lmr2buf)/=1)  then
         write (*,'(A,I0,A)') ' * Wrong lmr2buf! min=', minval(lmr2buf), ', should be 1'
         call abortRun(' * initialize_mlo_atoav_2d: failed sanity check test')
      end if
      if (maxval(lmr2buf)/=n_lm_loc)  then
         write (*,'(A,I0,A,I0)') ' * Wrong lmr2buf! max=', maxval(lmr2buf), ', should be ',n_lm_loc
         call abortRun(' * initialize_mlo_atoav_2d: failed sanity check test')
      end if
      if (sum(recv_mlo_count)/=n_r_max)  then
         write (*,'(A,I0,A,I0)') ' * Wrong recv_mlo_count! sum=', sum(recv_mlo_count), ', should be ',n_r_max
         call abortRun(' * initialize_mlo_atoav_2d: failed sanity check test')
      end if
      if (sum(send_lm_count)/=n_lm_loc)  then
         write (*,'(A,I0,A,I0)') ' * Wrong send_lm_count! sum=', sum(send_lm_count), ', should be ',n_lm_loc
         call abortRun(' * initialize_mlo_atoav_2d: failed sanity check test')
      end if
      
      mlo_atoav_2d_initialized = .true.
   end subroutine initialize_mlo_atoav_2d
   
   !----------------------------------------------------------------------------
   subroutine reorder_rloc2buffer(this, arr_Rloc, rloc_buf)
   !-- Reorder input variable into send buffer
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_2d)  :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: rloc_buf(:,:,:)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         rloc_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i)) = &
            arr_Rloc(i,nRstart:nRstop, 1:this%n_fields)
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2rloc(this, arr_Rloc, rloc_buf)
   !-- Reorder receive buffer into output variable
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_2d)  :: this
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(in)  :: rloc_buf(:,:,:)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         arr_Rloc(i,nRstart:nRstop, 1:this%n_fields) = &
         &  rloc_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i))
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2lmloc(this, arr_LMloc, lmloc_buf)
   !-- Reorder receive buffer onto output variable
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_2d)  :: this
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(in)  :: lmloc_buf(:)
      integer :: icoord_lo, lr, ur, lb, ub, i, n

      !>@TODO optimize the order of this loop
      do icoord_lo=0,n_ranks_r-1
         lr = dist_r(icoord_lo,1)
         ur = dist_r(icoord_lo,2)
         n  = dist_r(icoord_lo,0)
         lb = this%lmloc_disp(icoord_lo)+1
         ub = lb + n*this%n_fields-1
         do i=1, n_mlo_loc
            arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields) = reshape(lmloc_buf(lb:ub), (/n, this%n_fields/))
            lb = ub + 1
            ub = lb + n*this%n_fields-1
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_lmloc2buffer(this, arr_LMloc, lmloc_buf)
   !-- Reorder input variable into send buffer
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_2d)  :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: lmloc_buf(:)
      integer :: icoord_lo, lr, ur, lb, ub, i, n
      
      !>@TODO optimize the order of this loop
      do icoord_lo=0,n_ranks_r-1
         lr = dist_r(icoord_lo,1)
         ur = dist_r(icoord_lo,2)
         n  = dist_r(icoord_lo,0)
         lb = this%lmloc_disp(icoord_lo)+1
         ub = lb + n*this%n_fields-1
         do i=1, n_mlo_loc
            lmloc_buf(lb:ub) = reshape(arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields), (/n*this%n_fields/))
            lb = ub + 1
            ub = lb + n*this%n_fields-1
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine finalize_mlo_atoav_2d
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
      integer :: i, ierr
      
      deallocate(lmr2buf)
      deallocate(buf2mlo)
      mlo_atoav_2d_initialized = .false.
   end subroutine finalize_mlo_atoav_2d
   

   !----------------------------------------------------------------------------
   subroutine create_mlo_atoav_2d(this, n_fields)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpiatoav_2d) :: this
      integer, intent(in) :: n_fields
      integer :: i, buffer_size
      
      if (this%initialized) return 
      
      if (.not. mlo_atoav_2d_initialized) call initialize_mlo_atoav_2d
      this%n_fields = n_fields
      allocate(this%rloc_count(0:n_ranks_lo-1))
      allocate(this%lmloc_count(0:n_ranks_lo-1))
      allocate(this%rloc_disp(0:n_ranks_lo-1))
      allocate(this%lmloc_disp(0:n_ranks_lo-1))
      bytes_allocated = bytes_allocated+4*n_ranks_lo*SIZEOF_INTEGER
      
      this%rloc_count = send_lm_count*n_r_loc*n_fields
      this%lmloc_count = recv_mlo_count*n_mlo_loc*n_fields
      
      this%rloc_disp = 0
      this%lmloc_disp = 0
      do i=1,n_ranks_lo-1
         this%rloc_disp(i) = this%rloc_disp(i-1)+this%rloc_count(i-1)
         this%lmloc_disp(i) = this%lmloc_disp(i-1)+this%lmloc_count(i-1)
      end do
      this%rloc_sum = sum(this%rloc_count)
      this%lmloc_sum = sum(this%lmloc_count)
      buffer_size = max(this%rloc_sum, this%lmloc_sum)
      
      n_atoav_2d_obj = n_atoav_2d_obj + 1
      
   end subroutine create_mlo_atoav_2d
   
   !----------------------------------------------------------------------------
   subroutine destroy_mlo_atoav_2d(this)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   !   
      class(type_mpiatoav_2d) :: this
      integer :: i, ierr
      
      if (.not. this%initialized) return
      
      deallocate(this%rloc_count)
      deallocate(this%lmloc_count)
      deallocate(this%rloc_disp)
      deallocate(this%lmloc_disp)
      
      n_atoav_2d_obj = n_atoav_2d_obj - 1
      this%initialized = .false.
      
      if (n_atoav_2d_obj==0) call finalize_mlo_atoav_2d
      
   end subroutine destroy_mlo_atoav_2d
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoav_2d_dist(this, arr_LMloc, arr_Rloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_2d)     :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      !-- Local variables:
      integer :: ierr
      complex(cp), allocatable :: rloc_buf(:,:,:), lmloc_buf(:)
      

      PERFON('lm2rS')
      allocate(rloc_buf(n_r_loc, this%n_fields, n_lm_loc))
      allocate(lmloc_buf(n_r_max*this%n_fields*n_mlo_loc))
      call reorder_lmloc2buffer(this, arr_LMloc, lmloc_buf)
      PERFOFF
      PERFON('lm2rW')
      call mpi_alltoallv(lmloc_buf, this%lmloc_count, this%lmloc_disp, MPI_DEF_COMPLEX, &
                         rloc_buf, this%rloc_count, this%rloc_disp, MPI_DEF_COMPLEX, &
                         comm_r, ierr)
      PERFOFF
      PERFON('lm2rS')
      call reorder_buffer2rloc(this, arr_Rloc, rloc_buf)
      deallocate(rloc_buf)
      deallocate(lmloc_buf)
      PERFOFF
      
   end subroutine transp_lm2r_atoav_2d_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_atoav_2d_dist(this, arr_Rloc, arr_LMloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_2d)     :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      !-- Local variables:
      integer :: ierr
      complex(cp), allocatable :: rloc_buf(:,:,:), lmloc_buf(:)
      
      PERFON('r2lmS')
      allocate(rloc_buf(n_r_loc, this%n_fields, n_lm_loc))
      allocate(lmloc_buf(n_r_max*this%n_fields*n_mlo_loc))
      call reorder_rloc2buffer(this, arr_Rloc, rloc_buf)
      PERFOFF
      PERFON('r2lmW')
      call mpi_alltoallv(rloc_buf, this%rloc_count, this%rloc_disp, MPI_DEF_COMPLEX, &
                         lmloc_buf, this%lmloc_count, this%lmloc_disp, MPI_DEF_COMPLEX, &
                         comm_r, ierr)
      PERFOFF
      PERFON('r2lmS')
      call reorder_buffer2lmloc(this, arr_LMloc, lmloc_buf)
      deallocate(rloc_buf)
      deallocate(lmloc_buf)
      PERFOFF

   end subroutine transp_r2lm_atoav_2d_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoav_2d_dummy(this, arr_LMloc, arr_Rloc)
      class(type_mpiatoav_2d) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      print*, "Dummy transp_lm2r_atoav_2d_dummy, not yet implemented!"
   end subroutine transp_lm2r_atoav_2d_dummy

   subroutine transp_r2lm_atoav_2d_dummy(this, arr_Rloc, arr_LMloc)
      class(type_mpiatoav_2d) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      print*, "Dummy transp_r2lm_atoav_2d_dummy, not yet implemented!"
   end subroutine transp_r2lm_atoav_2d_dummy
   
end module mod_mpiatoav_2d
