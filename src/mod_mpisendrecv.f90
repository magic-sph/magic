module mod_mpisendrecv

   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp, only: type_mpitransp
   use fft
   use LMmapping


   type, extends(type_mpitransp) :: type_mpisendrecv
      integer, pointer :: rq(:)
      integer, pointer :: recvs(:)  ! pointer to rq
      integer, pointer :: sends(:)  ! pointer to rq
      logical :: initialized = .false.

   contains
      procedure :: create_comm => create_mlo_sendrecv
      procedure :: destroy_comm => destroy_mlo_sendrecv
      procedure :: transp_lm2r => transp_lm2r_sendrecv_dummy
      procedure :: transp_r2lm => transp_r2lm_sendrecv_dummy
      procedure :: transp_lm2r_dist => transp_lm2r_sendrecv_dist
      procedure :: transp_r2lm_dist => transp_r2lm_sendrecv_dist
      procedure :: lm2r_start => transp_lm2r_sendrecv_start
      procedure :: r2lm_start => transp_r2lm_sendrecv_start
      procedure :: lm2r_wait  => transp_lm2r_sendrecv_wait
      procedure :: r2lm_wait  => transp_r2lm_sendrecv_wait
   end type type_mpisendrecv
   
   !
   ! These are the same for all objects of the type type_mpisendrecv
   ! They were written for lm2r direction, but are also used for 
   ! r2lm direction, except that flipped (the position of 
   ! sendrecv_dests and sendrecv_sources and etc are inverted)
   !@>TODO maybe add a pointer within the object to these entities?
   integer, private, allocatable :: lm2r_s_type(:), lm2r_r_type(:)
   integer, private, allocatable :: lm2r_dests(:),  lm2r_sources(:)
   integer, private, allocatable :: lm2r_loc_dspl(:,:)
   logical, private :: mlo_sendrecv_initialized = .false.
   integer, private :: n_sendrecv_obj = 0
   
   public :: type_mpisendrecv
   
   contains
   
   subroutine initialize_mlo_sendrecv
      !-- Initialize the MPI types for the transposition from ML to Radial
      !   All type_mpisendrecv objects will use the MPI types defined here (hence 
      !   why we don't declare them once for each container)
      !   
      !   This is supposed to be a general-purpose transposition. It might be 
      !   possible to optimize it for further specialized data structures.
      !   The MPI datatype here is built as indexed->resized->vector->resized
      !   
      !   Author: Rafael Lago (MPCDF) January 2019
      !    
      !-- TODO
      !
      !-- Local variables
      integer :: i, j, l, m, icoord_m, icoord_mlo, icoord_r, in_r, lm, ierr
      integer :: col_type, ext_type, mtx_type
      integer :: send_displacements(n_mlo_array, 0:n_ranks_mlo-1)
      integer :: recv_displacements(n_lm_loc, 0:n_ranks_mlo-1)
      integer :: send_counter_i(0:n_ranks_mlo-1), recv_counter_i(0:n_ranks_mlo-1), inblocks
      integer :: blocklenghts(max(n_mlo_array, n_lm_loc))
      
      integer (kind=mpi_address_kind) :: lb, extend, bytesCMPLX
      integer :: nsends, nrecvs
      
      if (mlo_sendrecv_initialized) return
      
      call mpi_type_get_extent(MPI_DOUBLE_COMPLEX, lb, bytesCMPLX, ierr) 
      
      send_displacements = -1
      recv_displacements = -1
      send_counter_i     = 0
      recv_counter_i     = 0
      blocklenghts       = 1 ! This is actually fixed to 1!
      lb = 0
      
      !-- Loops over each (m,l) tuple to determine which rank in lmr will need it
      !   There is no good way of doing this; I could loop over the *local* tuples,
      !   but then I would get displacements in a funny order. The "easiest" 
      !   solution I see is to simply loop over *all* (m,l) tuples and skip those
      !   which do not belong here.
      do lm=1,lm_max
         l = map_glbl_st%lm2l(lm)
         m = map_glbl_st%lm2m(lm)
         i = map_mlo%ml2i(m,l)
         if (i<1) cycle
         
         icoord_m = m_tsid(m)
         
         do icoord_r=0, n_ranks_r-1
            icoord_mlo = mpi_map%lmr2mlo(icoord_m, icoord_r)
            inblocks = send_counter_i(icoord_mlo) + 1
            send_counter_i(icoord_mlo) = inblocks
            send_displacements(inblocks,icoord_mlo) = i-1
         end do
      end do
      
      !-- Counts how many ranks will receive from me, and then allocate the (global) 
      !   array of types. 
      !   We subtract 1, if this rank has to "send data to itself", because we will
      !   not create a datatype nor an MPI request for that. 
      !   I cannot fathom a realitic situation in which a rank would *not* have to
      !   send something to itself, but who knows, maybe I'm missing something
      nsends = count(send_counter_i>0)
      if (send_counter_i(coord_mlo)>0) nsends = nsends - 1
      allocate(lm2r_s_type(nsends))
      allocate(lm2r_dests(nsends))
      
      !-- Build the Send MPI types:
      !   First we create indexed type; it is a vector which looks more or less like
      !      send_buf({i1,i2,i3,...}, j, k)
      !   Then we resize it, so that it has the length of the first dimension. Next
      !   we build an vector type out of it obtaining
      !      send_buf({i1,i2,i3,...}, nRstart:nRstop, k)
      !   Finally, we resize that, so that it comprises the whole data for a single 
      !   variable. Since we send the whole containers at once, we will be sending
      !   and array of those last types
      j = 0
      do icoord_mlo=0,n_ranks_mlo-1
         
         inblocks = send_counter_i(icoord_mlo)
         if (inblocks>0 .and. icoord_mlo /= coord_mlo) then
            j=j+1
            lm2r_dests(j) = icoord_mlo
            icoord_r = mpi_map%mlo2lmr(icoord_mlo,2)
            in_r = dist_r(icoord_r,0)
            
            call MPI_Type_indexed(inblocks,blocklenghts(1:inblocks), send_displacements(1:inblocks,icoord_mlo), &
               MPI_DOUBLE_COMPLEX, col_type, ierr)
               
            extend = int(n_mlo_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(col_type, lb, extend, ext_type, ierr)
            
            call MPI_Type_vector(in_r, 1, 1, ext_type, mtx_type, ierr)
            
            extend = int(n_mlo_loc*n_r_max*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(mtx_type, lb, extend, lm2r_s_type(j), ierr)
            
            call MPI_Type_commit(lm2r_s_type(j),ierr)
            call MPI_Type_free(col_type,ierr)
            call MPI_Type_free(mtx_type,ierr)
            call MPI_Type_free(ext_type,ierr)
         end if
      end do
      
      !-- Loops over each local (l,m) tuple to figure out which rank will send it to me!
      do i=1,n_lm_loc
         m  = map_dist_st%lm2m(i)
         l  = map_dist_st%lm2l(i)
         icoord_mlo = mlo_tsid(m,l)
         
         inblocks = recv_counter_i(icoord_mlo) + 1
         recv_counter_i(icoord_mlo) = inblocks
         recv_displacements(inblocks,icoord_mlo) = i-1
      end do
      
      !-- Counts how many ranks will send to me (similar to nsends)
      nrecvs = count(recv_counter_i>0)
      if (send_counter_i(coord_mlo)>0) nrecvs = nrecvs - 1
      allocate(lm2r_r_type(nrecvs))
      allocate(lm2r_sources(nrecvs))
      
      !-- Build the Recv MPI types (similar to Send MPI Types)
      !
      j = 0
      do icoord_mlo=0,n_ranks_mlo-1

         inblocks = recv_counter_i(icoord_mlo)
         if (inblocks>0 .and. icoord_mlo /= coord_mlo) then
            j = j + 1
            lm2r_sources(j) = icoord_mlo
            call MPI_Type_indexed(inblocks,blocklenghts(1:inblocks),&
               recv_displacements(1:inblocks,icoord_mlo), MPI_DOUBLE_COMPLEX, col_type, ierr)
            
            extend = int(n_lm_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(col_type, lb, extend, ext_type, ierr)
            call MPI_Type_vector(n_r_loc, 1, 1, ext_type, mtx_type, ierr)
            
            extend = int(n_lm_loc*n_r_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(mtx_type, lb, extend, lm2r_r_type(j), ierr)
            
            call MPI_Type_commit(lm2r_r_type(j),ierr)
            call MPI_Type_free(col_type,ierr)
            call MPI_Type_free(mtx_type,ierr)
            call MPI_Type_free(ext_type,ierr)
         end if
      end do
      
      !-- Finally, we need to keep the info pertaining to the copying of 
      !   data which is already local
      inblocks = send_counter_i(coord_mlo)
      allocate(lm2r_loc_dspl(inblocks,2))
      lm2r_loc_dspl(:,1) = send_displacements(1:inblocks,coord_mlo) + 1
      lm2r_loc_dspl(:,2) = recv_displacements(1:inblocks,coord_mlo) + 1
      
      mlo_sendrecv_initialized = .true.
   end subroutine initialize_mlo_sendrecv
   
   !----------------------------------------------------------------------------
   subroutine finalize_mlo_sendrecv
      !
      !   Author: Rafael Lago (MPCDF) January 2018
      !   
      integer :: i, ierr
      
      do i=1,size(lm2r_s_type)
         call MPI_Type_free(lm2r_s_type(i),ierr)
      end do
      do i=1,size(lm2r_r_type)
         call MPI_Type_free(lm2r_r_type(i),ierr)
      end do
      
      deallocate(lm2r_s_type)
      deallocate(lm2r_r_type)
      deallocate(lm2r_dests)
      deallocate(lm2r_sources)
      deallocate(lm2r_loc_dspl)
      mlo_sendrecv_initialized = .false.
!       print *, "finalize_mlo_sendrecv" 
   end subroutine finalize_mlo_sendrecv
   

   !----------------------------------------------------------------------------
   subroutine create_mlo_sendrecv(this, n_fields)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpisendrecv) :: this
      integer, intent(in) :: n_fields
      integer :: nsends, nrecvs
      
      if (this%initialized) return
      if (.not. mlo_sendrecv_initialized) call initialize_mlo_sendrecv
      this%n_fields = n_fields
      
      nsends = size(lm2r_dests)
      nrecvs = size(lm2r_sources)
      allocate(this%rq(nsends+nrecvs))
      this%sends => this%rq(1:nsends)
      this%recvs => this%rq(nsends+1:nrecvs)
      
      this%rq = MPI_REQUEST_NULL
      
      n_sendrecv_obj = n_sendrecv_obj + 1
      this%initialized = .true.
      
   end subroutine create_mlo_sendrecv
   
   !----------------------------------------------------------------------------
   subroutine destroy_mlo_sendrecv(this)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
      
      class(type_mpisendrecv) :: this
      integer :: i, ierr
      
      if (.not. this%initialized) return
      do i=1,size(this%rq)
         if (this%rq(i)/=MPI_REQUEST_NULL) then
            call mpi_cancel(this%rq(i), ierr)
            call mpi_request_free(this%rq(i), ierr)
         end if
      end do
      
      nullify(this%sends)
      nullify(this%recvs)
      deallocate(this%rq)
      
      n_sendrecv_obj = n_sendrecv_obj - 1
      if (n_sendrecv_obj==0) call finalize_mlo_sendrecv
      
   end subroutine destroy_mlo_sendrecv
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_sendrecv_start(this, arr_LMloc, arr_Rloc)
   !-- General-purpose transposition. All the customization should happen 
   !   during the creation of the types!
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   transp_lm2r_sendrecv_start
      class(type_mpisendrecv)     :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      integer :: i, j, k, icoord_mlo, icoord_r, il_r, ierr
      
      !-- Starts the sends
      do j=1,size(lm2r_dests)
         icoord_mlo = lm2r_dests(j)
         icoord_r = mpi_map%mlo2lmr(icoord_mlo,2)
         il_r = dist_r(icoord_r,1)
         call mpi_isend(arr_LMloc(1,il_r,1), this%n_fields, lm2r_s_type(j), icoord_mlo, 1, comm_mlo, this%sends(j), ierr)
      end do
      
      !-- Starts the receives
      do j=1,size(lm2r_sources)
         icoord_mlo = lm2r_sources(j)
         call mpi_irecv(arr_Rloc, this%n_fields, lm2r_r_type(j), icoord_mlo, 1, comm_mlo, this%recvs(j), ierr)
      end do
      
      !-- Copies data which is already local
      do i=1,size(lm2r_loc_dspl,1)
         k = lm2r_loc_dspl(i,1)
         j = lm2r_loc_dspl(i,2)
         arr_Rloc(j,nRstart:nRstop,1:this%n_fields) = arr_LMloc(k,nRstart:nRstop,1:this%n_fields)
      end do
      
   end subroutine transp_lm2r_sendrecv_start
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_sendrecv_wait(this)
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
      class(type_mpisendrecv) :: this
      call mpi_waitall(size(this%rq),this%rq,MPI_STATUSES_IGNORE,ierr)
   end subroutine transp_lm2r_sendrecv_wait
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_sendrecv_dist(this, arr_LMloc, arr_Rloc)
   !   Author: Rafael Lago (MPCDF) May 2020
      class(type_mpisendrecv)     :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      call this%lm2r_start(arr_LMloc,arr_Rloc)
      call this%lm2r_wait()
   end subroutine transp_lm2r_sendrecv_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_sendrecv_start(this, arr_Rloc, arr_LMloc)
   !-- General-purpose transposition. All the customization should happen 
   !   during the creation of the types!
   !
   !   Author: Rafael Lago (MPCDF) October 2019
   !   
      class(type_mpisendrecv)     :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      integer :: i, j, k, icoord_mlo, icoord_r, il_r, ierr
      
      !-- Starts the receives
      do j=1,size(lm2r_sources)
         icoord_mlo = lm2r_sources(j)
         call mpi_isend(arr_Rloc, this%n_fields, lm2r_r_type(j), icoord_mlo, 1, comm_mlo, this%recvs(j), ierr)
      end do
      
      !-- Starts the sends
      do j=1,size(lm2r_dests)
         icoord_mlo = lm2r_dests(j)
         icoord_r = mpi_map%mlo2lmr(icoord_mlo,2)
         il_r = dist_r(icoord_r,1)
         call mpi_irecv(arr_LMloc(1,il_r,1), this%n_fields, lm2r_s_type(j), icoord_mlo, 1, comm_mlo, this%sends(j), ierr)
      end do
      
      !-- Copies data which is already local
      do i=1,size(lm2r_loc_dspl,1)
         k = lm2r_loc_dspl(i,1)
         j = lm2r_loc_dspl(i,2)
         arr_LMloc(k,nRstart:nRstop,1:this%n_fields) = arr_Rloc(j,nRstart:nRstop,1:this%n_fields)
      end do
   end subroutine transp_r2lm_sendrecv_start
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_sendrecv_wait(this)
   !
   !   Author: Rafael Lago (MPCDF) October 2019
   !   
      class(type_mpisendrecv) :: this
      call mpi_waitall(size(this%rq),this%rq,MPI_STATUSES_IGNORE,ierr)
   end subroutine transp_r2lm_sendrecv_wait
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_sendrecv_dist(this, arr_Rloc, arr_LMloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpisendrecv)     :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      call this%r2lm_start(arr_Rloc, arr_LMloc)
      call this%r2lm_wait()
   end subroutine transp_r2lm_sendrecv_dist
   
   subroutine transp_lm2r_sendrecv_dummy(this, arr_LMloc, arr_Rloc)
      class(type_mpisendrecv) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      print*, "Dummy transp_lm2r_sendrecv_dummy, not yet implemented!"
   end subroutine transp_lm2r_sendrecv_dummy

   subroutine transp_r2lm_sendrecv_dummy(this, arr_Rloc, arr_LMloc)
      class(type_mpisendrecv) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      print*, "Dummy transp_r2lm_sendrecv_dummy, not yet implemented!"
   end subroutine transp_r2lm_sendrecv_dummy
   
end module mod_mpisendrecv
