!----------------------------------------------------------------------------------
   !
module mpi_thetap_mod
   ! 
   ! This module contains the implementation of theta-parallel transpositions
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

   type, extends(type_mpitransp) :: type_mpisendrecv
      integer, pointer :: rq(:)
      integer, pointer :: recvs(:)  ! pointer to rq
      integer, pointer :: sends(:)  ! pointer to rq

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
   


  public :: transpose_m_theta, transpose_theta_m, transform_m2phi,            &
     & transform_phi2m, type_mpisendrecv, transform_new2old, transform_old2new, test_field

contains

   !-- Transposition from (n_lm_loc,nRStart) to (n_mlo_loc,n_r_max).
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   
   
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
      do i=1,size(lm2r_dests)
         call MPI_Type_free(lm2r_dests(i),ierr)
      end do
      
      deallocate(lm2r_s_type)
      deallocate(lm2r_r_type)
      deallocate(lm2r_dests)
      deallocate(lm2r_sources)
      deallocate(lm2r_loc_dspl)
      mlo_sendrecv_initialized = .false.
   end subroutine finalize_mlo_sendrecv
   

   !----------------------------------------------------------------------------
   subroutine create_mlo_sendrecv(this, n_fields)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpisendrecv) :: this
      integer, intent(in) :: n_fields
      integer :: nsends, nrecvs
           
      if (.not. mlo_sendrecv_initialized) call initialize_mlo_sendrecv
      this%n_fields = n_fields
      
      nsends = size(lm2r_dests)
      nrecvs = size(lm2r_sources)
      allocate(this%rq(nsends+nrecvs))
      this%sends => this%rq(1:nsends)
      this%recvs => this%rq(nsends+1:nrecvs)
      
      this%rq = MPI_REQUEST_NULL
      
      n_sendrecv_obj = n_sendrecv_obj + 1
      
   end subroutine create_mlo_sendrecv
   
   !----------------------------------------------------------------------------
   subroutine destroy_mlo_sendrecv(this)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpisendrecv) :: this
      integer :: i, ierr
      
      do i=1,size(this%rq)
         call mpi_cancel(this%rq(i), ierr)
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
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
      
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   !-- Transposition from (m_loc,θ_glb) to (θ_loc,m_glb).
   !   
   !   
   !   Author: Rafael Lago (MPCDF) August 2017
   !
   !-- TODO this with mpi_type to stride the data and check if performance 
   !--      improves
   !
   subroutine transpose_m_theta(f_m_theta, f_theta_m)
      complex(cp), intent(inout) :: f_m_theta(n_m_max, n_theta_loc)
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      
      complex(cp) :: sendbuf(n_m_max * n_theta_loc)
      complex(cp) :: recvbuf(n_m_loc, n_theta_max)
      
      integer :: sendcount(0:n_ranks_m-1)
      integer :: recvcount(0:n_ranks_m-1)
      integer :: senddispl(0:n_ranks_m-1)
      integer :: recvdispl(0:n_ranks_m-1)
      integer :: irank, j, itheta, m, pos
      
      pos = 1
      do irank=0,n_ranks_m-1
         !-- Copy each m which belongs to the irank-th coord_r into the send buffer
         !   column-wise. That will simplify a lot things later
         !
         !@>TODO check performance of this; implementing this with mpi_type
         !  striding the data could be faster
         senddispl(irank) = pos-1
         do itheta=1,n_theta_loc
            do j=1,dist_m(irank,0)
               m = dist_m(irank,j)/minc
               sendbuf(pos) = f_m_theta(m+1,itheta)
               pos = pos + 1
            end do
         end do
         
         sendcount(irank) = pos - senddispl(irank) - 1
         recvdispl(irank) = irank*n_m_loc*dist_theta(irank,0)
         recvcount(irank) =   n_m_loc*dist_theta(irank,0)
      end do
      
      call MPI_Alltoallv(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
                         comm_m, irank)
      f_theta_m = transpose(recvbuf)
      
   end subroutine transpose_m_theta
   
   !-- Transposition from (θ_loc,m_glb) to (m_loc,θ_glb)
   !   
   !   Author: Rafael Lago (MPCDF) August 2017
   !
   !-- TODO this with mpi_type to stride the data
   !
   subroutine transpose_theta_m(f_theta_m, f_m_theta)
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      complex(cp), intent(inout) :: f_m_theta(n_m_max, n_theta_loc)
      
      complex(cp) :: sendbuf(n_m_loc * n_theta_max)
      complex(cp) :: recvbuf(n_theta_loc,  n_m_max)
      
      integer :: sendcount(0:n_ranks_theta-1)
      integer :: recvcount(0:n_ranks_theta-1)
      integer :: senddispl(0:n_ranks_theta-1)
      integer :: recvdispl(0:n_ranks_theta-1)
      integer :: irank, j, pos, n_t, l_t, u_t
      integer :: m_arr(n_ranks_theta*n_m_array) 
      
      recvcount = 0
      pos = 1
      do irank=0,n_ranks_theta-1
         !-- Copy each theta chunk so that the send buffer is contiguous
         !-- TODO check performance of this; implementing this with mpi_type
         !   striding the data will probably be faster
         senddispl(irank) = pos-1
         n_t = dist_theta(irank,0)
         l_t = dist_theta(irank,1)
         u_t = dist_theta(irank,2)
         do j=1, n_m_loc
            sendbuf(pos:pos + n_t - 1) = f_theta_m(l_t:u_t,j)
            pos = pos + n_t
         end do
         
         sendcount(irank) = pos - senddispl(irank) - 1
         recvdispl(irank) = sum(recvcount)
         recvcount(irank) = dist_m(irank,0) * n_t
      end do
      
      call MPI_Alltoallv(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
                         comm_theta, irank)
      
      !-- Now we reorder the receiver buffer. If the m distribution looks like:
      !   coord_r 0: 0, 4,  8, 12, 16
      !   coord_r 1: 1, 5,  9, 13
      !   coord_r 2: 2, 6, 10, 14
      !   coord_r 3: 3, 7, 11, 15
      !   then the columns of recvbuf are ordered as 0,4,8,12,16,1,5,9,13(...)
      !   and so forth. m_arr will contain this ordering (+1):
      m_arr = reshape(transpose(dist_m(:,1:)), &
                      (/n_ranks_m*n_m_array/))/minc + 1
      j = 1
      do pos = 1, n_ranks_theta*n_m_array
         if (m_arr(pos) < 1) cycle
         f_m_theta(m_arr(pos),:) = recvbuf(:,j)
         j = j + 1
      end do
   end subroutine transpose_theta_m
   
   !-- Transforms from (θ,m) space into (φ,θ) space including transpositions 
   !   and FFT. 
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   !   TODO: some functions might requires this transformation for multiple
   !      fields at once (e.g. vector transform). In that case, a non-blocking 
   !      transpose_theta_m would make more sense. This would make this 
   !      function obsolete.
   !   TODO: there is a lot of room for immprovement here (e.g. in-place, 
   !     use just one intermediate buffer, vector transform, etc)
   !
   subroutine transform_m2phi(fL, f)
      
      !-- Input variables
      complex(cp), intent(inout) :: fL(n_theta_max,n_m_loc)
      
      !-- Output variables
      real(cp),    intent(out)   :: f(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: lF(n_m_max,n_theta_loc)
      complex(cp) :: Ff(n_phi_max/2+1,n_theta_loc)
   
      call transpose_theta_m(fL, lF)
      !-- TODO: The FFT must be performed for an array with the dimensions of 
      !   F_loc which may end up paded with zeroes.
      !   Is there any way to tell MKL to perform a "truncated" FFT?
      Ff = 0.0
      Ff(1:n_m_max,1:n_theta_loc) = lF
      
      call fft_phi_loc(f, Ff, -1)
   end subroutine transform_m2phi
   
   
   !-- Transforms from (φ,θ) space into (θ,m) space including transpositions 
   !   and FFT. 
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   !   TODO: some functions might requires this transformation for multiple
   !      fields at once (e.g. vector transform). In that case, a non-blocking 
   !      transpose_theta_m would make more sense. This would make this 
   !      function obsolete.
   !   TODO: there is a lot of room for immprovement here (e.g. in-place, 
   !     use just one intermediate buffer, vector transform, etc)
   !
   subroutine transform_phi2m(f, fL)
      
      !-- Input variables
      real(cp),    intent(inout) :: f(n_phi_max, n_theta_loc)
      
      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      
      !-- Local variables
      complex(cp) :: lF(n_m_max,n_theta_loc)
      complex(cp) :: Ff(n_phi_max/2+1,n_theta_loc)
   
      call fft_phi_loc(f, Ff, 1)
      lF(1:n_m_max,1:n_theta_loc) = Ff(1:n_m_max,1:n_theta_loc)
      call transpose_m_theta(lF, fL)
   end subroutine transform_phi2m
   
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!  From here on, are only temporary functions. They are not optimized
!  (actually they are pretty terrible). They 
!  should be deleted once the transition is over
!  
!  LM Loop transposes and Gathers and Etc
!
!-------------------------------------------------------------------------------
   subroutine transform_new2old(Fmlo_new, Fmlo_old, n_r)
      integer,     intent(in) :: n_r ! Needed since it can be either n_r_ic_max, or n_r_max
      complex(cp), intent(in) :: Fmlo_new(n_mlo_loc, n_r)
      complex(cp), intent(inout) :: Fmlo_old(llm:ulm,   n_r)
      
      complex(cp) :: recvbuff(n_r)
      integer :: irank, ierr, lm, l, m, lo
      
      do lm=1,lm_max
         m = map_glbl_st%lm2m(lm)
         l = map_glbl_st%lm2l(lm)
         irank = map_mlo%ml2coord(m,l)
         recvbuff = 0.0
         if (irank==coord_mlo) recvbuff = Fmlo_new(map_mlo%ml2i(m,l),:)
         call mpi_bcast(recvbuff, n_r, MPI_DOUBLE_COMPLEX, irank, comm_mlo, ierr)
         lo = lo_map%lm2(l,m)
         if (lo>=llm .and. lo<=ulm) Fmlo_old(lo,:) = recvbuff
      end do
   end subroutine transform_new2old
   
   
!-------------------------------------------------------------------------------
   subroutine transform_old2new(Fmlo_old, Fmlo_new, n_r)
!-------------------------------------------------------------------------------
      integer,     intent(in) :: n_r ! Needed since it can be either n_r_ic_max, or n_r_max
      complex(cp), intent(in) :: Fmlo_old(llm:ulm,   n_r)
      complex(cp), intent(inout) :: Fmlo_new(n_mlo_loc, n_r)
      
      complex(cp) :: recvbuff(n_r)
      integer :: old2coord(l_max, l_max)
      integer :: irank, ierr, lm, l, m, size_irank, i
      
      do irank=0,n_ranks_r-1
         ! lolololololol
         if (irank==coord_r) then
            call mpi_bcast(ulm-llm+1, 1, MPI_INTEGER, irank, comm_r, ierr)
            do lm=llm,ulm
               m = lo_map%lm2m(lm)
               l = lo_map%lm2l(lm)
               
               call mpi_bcast(m, 1, MPI_INTEGER, irank, comm_r, ierr)
               call mpi_bcast(l, 1, MPI_INTEGER, irank, comm_r, ierr)
               
               call mpi_bcast(Fmlo_old(lm,:), n_r, MPI_DOUBLE_COMPLEX, irank, comm_r, ierr)
               if (map_mlo%ml2coord(m,l)==coord_mlo) Fmlo_new(map_mlo%ml2i(m,l),:) = Fmlo_old(lm,:)
            end do
         else
            call mpi_bcast(size_irank, 1, MPI_INTEGER, irank, comm_r, ierr)
            do i=1, size_irank
               call mpi_bcast(m, 1, MPI_INTEGER, irank, comm_r, ierr)
               call mpi_bcast(l, 1, MPI_INTEGER, irank, comm_r, ierr)
               call mpi_bcast(recvbuff, n_r, MPI_DOUBLE_COMPLEX, irank, comm_r, ierr)
               if (map_mlo%ml2coord(m,l)==coord_mlo) Fmlo_new(map_mlo%ml2i(m,l),:) = recvbuff
            end do
         end if
         
      end do
   end subroutine transform_old2new
!--------------------------------------------------------------------------------
!@> Delete me after conversion!
   subroutine test_field(newfield, oldfield, name, n_r)
      integer,     intent(in) :: n_r ! Needed since it can be either n_r_ic_max, or n_r_max
      character(len=*), intent(in) :: name
      complex(cp), intent(in) :: newfield(n_mlo_loc,n_r)
      complex(cp), intent(in) :: oldfield(llm:ulm,n_r)
      complex(cp) :: test_old(llm:ulm,n_r)
      complex(cp) :: test_new(n_mlo_loc,n_r)
      real(cp)    :: test_norm, error_threshold
      
      error_threshold = 1e-16 !EPSILON(1.0_cp)
      
      call transform_new2old(newfield, test_old, n_r)
      test_norm = ABS(SUM(oldfield - test_old))
      IF (test_norm>error_threshold) print *, "||",name,"n2o|| : ", test_norm

      call transform_old2new(oldfield, test_new, n_r)
      test_norm = ABS(SUM(newfield - test_new))
      IF (test_norm>error_threshold) print *, "||",name,"o2n|| : ", test_norm
      
   end subroutine
   
end module mpi_thetap_mod
