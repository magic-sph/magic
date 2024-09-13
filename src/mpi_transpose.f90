#define KNL_BIG 0
module mpi_transp_mod
   !
   ! This is an abstract class that will be used to define MPI transposers
   ! The actual implementation is deferred to either point-to-point (MPI_Isend
   ! and MPI_IRecv) communications or all-to-all (MPI_AlltoAll)
   !

   use precision_mod
   use truncation, only: lm_max, n_r_max
   use radial_data, only: nRstart, nRstop
   use blocking, only: llm, ulm

   implicit none

   private

   type, abstract, public :: type_mpitransp
      integer :: n_fields
   contains
      procedure(create_if), deferred :: create_comm
      procedure(destroy_if), deferred :: destroy_comm
      procedure(transp_lm2r_if), deferred :: transp_lm2r
      procedure(transp_r2lm_if), deferred :: transp_r2lm
   end type type_mpitransp

   interface

      subroutine create_if(this, n_fields)
         import
         class(type_mpitransp) :: this
         integer, intent(in) :: n_fields
      end subroutine create_if

      subroutine destroy_if(this)
         import
         class(type_mpitransp) :: this
      end subroutine destroy_if

      subroutine transp_lm2r_if(this, arr_LMloc, arr_Rloc)
         import
         class(type_mpitransp) :: this
         complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
         complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      end subroutine transp_lm2r_if

      subroutine transp_r2lm_if(this, arr_Rloc, arr_LMloc)
         import
         class(type_mpitransp) :: this
         complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
         complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      end subroutine transp_r2lm_if

   end interface

end module mpi_transp_mod


!----------------------------------------------------------------------------------
module  mpi_alltoall_mod
   !
   ! This module contains the implementation of all-to-all global communicators
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use radial_data, only: radial_balance
   use truncation, only: lm_max, n_r_max, l_max, minc
   use radial_data, only: nRstart, nRstop
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp_mod, only: type_mpitransp
   use constants, only: zero

   implicit none

   private

   type, public, extends(type_mpitransp) :: type_mpiatoav
      integer, allocatable :: rcounts(:)
      integer, allocatable :: scounts(:)
      integer, allocatable :: rdisp(:)
      integer, allocatable :: sdisp(:)
      complex(cp), allocatable :: sbuff(:)
      complex(cp), allocatable :: rbuff(:)
   contains
      procedure :: create_comm => create_comm_alltoallv
      procedure :: destroy_comm => destroy_comm_alltoallv
      procedure :: transp_lm2r => transp_lm2r_alltoallv
      procedure :: transp_r2lm => transp_r2lm_alltoallv
   end type type_mpiatoav

   type, public, extends(type_mpitransp) :: type_mpiatoap
      integer :: rcounts
      integer :: scounts
      integer :: n_r_loc
      integer :: lm_loc
      complex(cp), allocatable :: buff(:)
   contains
      procedure :: create_comm => create_comm_alltoallp
      procedure :: destroy_comm => destroy_comm_alltoallp
      procedure :: transp_lm2r => transp_lm2r_alltoallp
      procedure :: transp_r2lm => transp_r2lm_alltoallp
   end type type_mpiatoap

   type, public, extends(type_mpitransp) :: type_mpiatoaw
      integer, allocatable :: rtype(:)
      integer, allocatable :: stype(:)
      integer, allocatable :: disp(:)
      integer, allocatable :: counts(:)
   contains
      procedure :: create_comm => create_comm_alltoallw
      procedure :: destroy_comm => destroy_comm_alltoallw
      procedure :: transp_lm2r => transp_lm2r_alltoallw
      procedure :: transp_r2lm => transp_r2lm_alltoallw
   end type type_mpiatoaw

contains

   subroutine create_comm_alltoallv(this, n_fields)

      class(type_mpiatoav) :: this
      integer, intent(in) :: n_fields

      !-- Local variables
      integer :: p, nlm_per_rank, my_lm_counts

      this%n_fields = n_fields

      allocate ( this%rcounts(0:n_procs-1), this%scounts(0:n_procs-1) )
      allocate ( this%rdisp(0:n_procs-1), this%sdisp(0:n_procs-1) )
      bytes_allocated = bytes_allocated +4*n_procs*SIZEOF_INTEGER

      do p=0,n_procs-1
         my_lm_counts = lm_balance(p)%n_per_rank
         nlm_per_rank = ulm-llm+1
         this%scounts(p)=nR_per_rank*my_lm_counts*this%n_fields
         this%rcounts(p)=radial_balance(p)%n_per_rank*nlm_per_rank*this%n_fields
      end do

      this%rdisp(0)=0
      this%sdisp(0)=0
      do p=1,n_procs-1
         this%sdisp(p)=this%sdisp(p-1)+this%scounts(p-1)
         this%rdisp(p)=this%rdisp(p-1)+this%rcounts(p-1)
      end do

      allocate( this%rbuff(sum(this%rcounts)), this%sbuff(sum(this%scounts)) )
      this%rbuff(:)=zero
      this%sbuff(:)=zero
      bytes_allocated=bytes_allocated+(sum(this%rcounts)+sum(this%scounts))* &
      &               SIZEOF_DEF_COMPLEX

   end subroutine create_comm_alltoallv
!----------------------------------------------------------------------------------
   subroutine create_comm_alltoallp(this, n_fields)

      class(type_mpiatoap) :: this
      integer, intent(in) :: n_fields

      this%n_r_loc = maxval(radial_balance(:)%n_per_rank)
      this%lm_loc = maxval(lm_balance(:)%n_per_rank)

      this%n_fields = n_fields
      this%scounts  = this%n_r_loc*this%lm_loc*this%n_fields
      this%rcounts  = this%n_r_loc*this%lm_loc*this%n_fields

      allocate( this%buff(n_procs*this%lm_loc*this%n_r_loc*this%n_fields) )
      this%buff(:)=zero
      bytes_allocated=bytes_allocated+n_procs*this%lm_loc*this%n_r_loc*this%n_fields*&
      &               SIZEOF_DEF_COMPLEX

   end subroutine create_comm_alltoallp
!----------------------------------------------------------------------------------
   subroutine create_comm_alltoallw(this, n_fields)

      class(type_mpiatoaw) :: this
      integer, intent(in) :: n_fields

      !-- Local variables
#ifdef WITH_MPI
      integer (kind=mpi_address_kind) :: lb, extend, bytesCMPLX
#endif
      integer :: arr_size(3), arr_loc_size(3), arr_start(3)
      integer :: p, my_lm_counts, nlm_per_rank
      integer :: k, col_type, ext_type, lm, l, m
      integer, allocatable :: indices(:), blocklengths(:)

      integer :: lm2l_balance(lm_max, 0:n_procs-1), tmp(lm_max)
      integer :: lm2m_balance(lm_max, 0:n_procs-1)

      tmp(:)=-1
      tmp(llm:ulm)=lo_map%lm2l(llm:ulm)
#ifdef WITH_MPI
      call MPI_Allgather(tmp, lm_max, MPI_INTEGER, lm2l_balance, lm_max, MPI_INTEGER, &
           &             MPI_COMM_WORLD, ierr)
#endif
      tmp(:)=-1
      tmp(llm:ulm)=lo_map%lm2m(llm:ulm)
#ifdef WITH_MPI
      call MPI_Allgather(tmp, lm_max, MPI_INTEGER, lm2m_balance, lm_max, MPI_INTEGER, &
           &             MPI_COMM_WORLD, ierr)
#endif

#ifdef WITH_MPI
      call MPI_Type_Get_Extent(MPI_DEF_COMPLEX, lb, bytesCMPLX, ierr)
#endif

      allocate(indices(maxval(lm_balance(:)%n_per_rank)))
      allocate(blocklengths(maxval(lm_balance(:)%n_per_rank)))
      blocklengths(:)=1
      this%n_fields = n_fields

      allocate ( this%counts(0:n_procs-1), this%disp(0:n_procs-1) )
      allocate ( this%rtype(0:n_procs-1), this%stype(0:n_procs-1) )
      bytes_allocated = bytes_allocated+4*n_procs*SIZEOF_INTEGER

      do p=0,n_procs-1
         my_lm_counts = lm_balance(p)%n_per_rank
         nlm_per_rank = ulm-llm+1

         this%counts(p)=1
         this%disp(p)  =0

         k = 0
         indices(:)=-1
         do lm=lm_balance(p)%nStart, lm_balance(p)%nStop
            l = lm2l_balance(lm,p)
            m = lm2m_balance(lm,p)
            if ( l >= 0 .and. m >= 0 ) then
               k = k+1
               indices(k) = st_map%lm2(l,m)-1
            end if
         end do

#ifdef WITH_MPI
         if ( k > 0 ) then
            call MPI_Type_Indexed(k, blocklengths(1:k), indices(1:k), MPI_DEF_COMPLEX, &
                 &                col_type, ierr)
            extend = int(lm_max*bytesCMPLX, kind=mpi_address_kind)
            call MPI_Type_create_resized(col_type, lb, extend, ext_type, ierr)
            call MPI_Type_contiguous((nRstop-nRstart+1)*this%n_fields, ext_type, &
                 &                   this%stype(p), ierr)

            call MPI_Type_commit(this%stype(p), ierr)

            call MPI_Type_free(col_type, ierr)
            call MPI_Type_free(ext_type, ierr)
         end if
#endif

         arr_size(1)=nlm_per_rank
         arr_size(2)=n_r_max
         arr_size(3)=this%n_fields
         arr_loc_size(1)=nlm_per_rank
         arr_loc_size(2)=radial_balance(p)%n_per_rank
         arr_loc_size(3)=this%n_fields
         arr_start(1)=0
         arr_start(2)=radial_balance(p)%nStart-1
         arr_start(3)=0
#ifdef WITH_MPI
         call MPI_Type_Create_Subarray(3, arr_size, arr_loc_size, arr_start, &
              &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,   &
              &                        this%rtype(p), ierr)
         call MPI_Type_Commit(this%rtype(p), ierr)
#endif
      end do

      deallocate(indices, blocklengths)

   end subroutine create_comm_alltoallw
!----------------------------------------------------------------------------------
   subroutine destroy_comm_alltoallv(this)

      class(type_mpiatoav) :: this

      deallocate( this%rbuff, this%sbuff )
      deallocate( this%sdisp, this%rdisp, this%scounts, this%rcounts )

   end subroutine destroy_comm_alltoallv
!----------------------------------------------------------------------------------
   subroutine destroy_comm_alltoallp(this)

      class(type_mpiatoap) :: this

      deallocate( this%buff )

   end subroutine destroy_comm_alltoallp
!----------------------------------------------------------------------------------
   subroutine destroy_comm_alltoallw(this)

      class(type_mpiatoaw) :: this

#ifdef WITH_MPI
      !-- Local variables
      integer :: p

      do p = 0, n_procs-1
         call MPI_Type_Free(this%rtype(p), ierr)
         call MPI_Type_Free(this%stype(p), ierr)
      end do
#endif
      deallocate( this%counts, this%disp, this%rtype, this%stype )

   end subroutine destroy_comm_alltoallw
!----------------------------------------------------------------------------------
   subroutine transp_lm2r_alltoallv(this, arr_LMloc, arr_Rloc)
      !
      ! This subroutine transposes a LM-distributed container of arrays into
      ! a r-distributed container of arrays
      !

      class(type_mpiatoav) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)

      !-- Local variables
      integer :: p, ii, n_r, lm, l, m, lm_st, n_f

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm)
      do p = 0, n_procs-1
         ii = this%rdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do lm=llm,ulm
                  this%rbuff(ii)=arr_LMloc(lm,n_r,n_f)
                  ii = ii+1
               end do
            end do
         end do
      end do
      !$omp end parallel do

#ifdef WITH_MPI
      call MPI_Alltoallv(this%rbuff, this%rcounts, this%rdisp, MPI_DEF_COMPLEX, &
           &             this%sbuff, this%scounts, this%sdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr)
#endif

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm,l,m,lm_st)
      do p = 0, n_procs-1
         ii = this%sdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
                  l = lo_map%lm2l(lm)
                  m = lo_map%lm2m(lm)
                  lm_st = st_map%lm2(l,m)
                  arr_Rloc(lm_st,n_r,n_f)=this%sbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine transp_lm2r_alltoallv
!----------------------------------------------------------------------------------
   subroutine transp_lm2r_alltoallp(this, arr_LMloc, arr_Rloc)
      !
      ! This subroutine transposes a LM-distributed container of arrays into
      ! a r-distributed container of arrays
      !

      class(type_mpiatoap) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)

      !-- Local variables
      integer :: p, ii, n_r, lm, l, m, lm_st, n_f

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm)
      !ii=1
      do p = 0, n_procs-1
         ii = p*this%rcounts+1
         do n_f=1,this%n_fields
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do lm=llm,ulm
                  this%buff(ii)=arr_LMloc(lm,n_r,n_f)
                  ii = ii+1
               end do
               !-- Padding happens here
               do lm = ulm-llm+1,this%lm_loc-1
                  ii = ii+1
               end do
            end do
            !-- Padding happens here
            do n_r=radial_balance(p)%n_per_rank,this%n_r_loc-1
               ii=ii+1
            end do
         end do
      end do
      !$omp end parallel do

#ifdef WITH_MPI
      call MPI_Alltoall(MPI_IN_PLACE, this%rcounts, MPI_DEF_COMPLEX, this%buff, &
           &            this%scounts, MPI_DEF_COMPLEX, MPI_COMM_WORLD, ierr)
#endif

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm,l,m,lm_st)
      !ii=1
      do p = 0, n_procs-1
         ii = p*this%scounts+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
                  l = lo_map%lm2l(lm)
                  m = lo_map%lm2m(lm)
                  lm_st = st_map%lm2(l,m)
                  arr_Rloc(lm_st,n_r,n_f)=this%buff(ii)
                  ii=ii+1
               end do
               do lm = lm_balance(p)%n_per_rank,this%lm_loc-1
                  ii = ii+1
               end do
            end do
            do n_r=nR_per_rank,this%n_r_loc-1
               ii=ii+1
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine transp_lm2r_alltoallp
!----------------------------------------------------------------------------------
   subroutine transp_lm2r_alltoallw(this, arr_LMloc, arr_Rloc)

      class(type_mpiatoaw) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)

#ifdef WITH_MPI
      call MPI_Alltoallw(arr_LMloc, this%counts, this%disp, this%rtype, &
           &             arr_Rloc, this%counts, this%disp, this%stype,  &
           &             MPI_COMM_WORLD, ierr)
#endif

   end subroutine transp_lm2r_alltoallw
!----------------------------------------------------------------------------------
   subroutine transp_r2lm_alltoallv(this, arr_Rloc, arr_LMloc)
      !
      ! This subroutine transposes a r-distributed container of arrays into
      ! a LM-distributed container of arrays
      !

      class(type_mpiatoav) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)

      !-- Local variables
#if (KNL_BIG==1)
      complex(cp) :: temp_Rloc(lm_max,nRstart:nRstop,this%n_fields)
      integer :: p, ii, n_r, lm, l, m, n_f

      !$omp parallel default(shared) private(p,ii,n_f,n_r,lm,l,m)
      !$omp do collapse(3)
      do n_f=1,this%n_fields
         do n_r=nRstart,nRstop
            do lm=1,lm_max
               l = lo_map%lm2l(lm)
               m = lo_map%lm2m(lm)
               temp_Rloc(lm,n_r,n_f)=arr_Rloc(st_map%lm2(l,m),n_r,n_f)
            end do
         end do
      end do
      !$omp end do

      !$omp do
      do p = 0, n_procs-1
         ii = this%sdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
                  this%sbuff(ii)=temp_Rloc(lm,n_r,n_f)
                  ii = ii +1
               end do
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

#elif (KNL_BIG==0)
      integer :: p, ii, n_r, lm, l, m, lm_st, n_f

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm,l,m,lm_st)
      do p = 0, n_procs-1
         ii = this%sdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
                  l = lo_map%lm2l(lm)
                  m = lo_map%lm2m(lm)
                  lm_st = st_map%lm2(l,m)
                  this%sbuff(ii)=arr_Rloc(lm_st,n_r,n_f)
                  ii = ii +1
               end do
            end do
         end do
      end do
      !$omp end parallel do
#endif

#ifdef WITH_MPI
      call MPI_Alltoallv(this%sbuff, this%scounts, this%sdisp, MPI_DEF_COMPLEX, &
           &             this%rbuff, this%rcounts, this%rdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr)
#endif

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm)
      do p = 0, n_procs-1
         ii = this%rdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do lm=llm,ulm
                  arr_LMloc(lm,n_r,n_f)=this%rbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine transp_r2lm_alltoallv
!----------------------------------------------------------------------------------
   subroutine transp_r2lm_alltoallp(this, arr_Rloc, arr_LMloc)
      !
      ! This subroutine transposes a r-distributed container of arrays into
      ! a LM-distributed container of arrays
      !

      class(type_mpiatoap) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)

      !-- Local variables
      integer :: p, ii, n_r, lm, l, m, lm_st, n_f

      !ii = 1
      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm,l,m,lm_st)
      do p = 0, n_procs-1
         ii = p*this%scounts+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
                  l = lo_map%lm2l(lm)
                  m = lo_map%lm2m(lm)
                  lm_st = st_map%lm2(l,m)

                  this%buff(ii)=arr_Rloc(lm_st,n_r,n_f)
                  ii = ii +1
               end do
               !-- padding happens here
               do lm = lm_balance(p)%n_per_rank, this%lm_loc-1
                  ii = ii+1
               end do
            end do
            !-- padding happens here
            do n_r=nR_per_rank,this%n_r_loc-1
               ii=ii+1
            end do
         end do
      end do
      !$omp end parallel do

#ifdef WITH_MPI
      call MPI_Alltoall(MPI_IN_PLACE, this%scounts, MPI_DEF_COMPLEX, this%buff, &
           &            this%rcounts, MPI_DEF_COMPLEX, MPI_COMM_WORLD, ierr)
#endif

      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_f,n_r,lm)
      do p = 0, n_procs-1
         ii = p*this%rcounts+1
         do n_f=1,this%n_fields
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do lm=llm,ulm
                  arr_LMloc(lm,n_r,n_f)=this%buff(ii)
                  ii=ii+1
               end do
               do lm = ulm-llm+1, this%lm_loc-1
                  ii = ii+1
               end do
            end do
            do n_r=radial_balance(p)%n_per_rank,this%n_r_loc-1
               ii=ii+1
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine transp_r2lm_alltoallp
!----------------------------------------------------------------------------------
   subroutine transp_r2lm_alltoallw(this, arr_Rloc, arr_LMloc)

      class(type_mpiatoaw) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)

#ifdef WITH_MPI
      call MPI_Alltoallw(arr_Rloc, this%counts, this%disp, this%stype,  &
           &             arr_LMloc, this%counts, this%disp, this%rtype, &
           &             MPI_COMM_WORLD, ierr)
#endif

   end subroutine transp_r2lm_alltoallw
!----------------------------------------------------------------------------------
end module mpi_alltoall_mod



!----------------------------------------------------------------------------------
module  mpi_ptop_mod
   !
   ! This module contains the implementation of MPI_Isend/MPI_Irecv global
   ! transpose
   !

   use precision_mod
   use mem_alloc
   use parallel_mod
   use truncation, only: l_max, minc
   use logic, only: l_finite_diff
   use truncation, only: lm_max, n_r_max
   use radial_data, only: nRstart, nRstop, radial_balance
   use blocking, only: lm_balance, st_map, lo_map, llm, ulm
   use mpi_transp_mod, only: type_mpitransp

   implicit none

   private

   type, public, extends(type_mpitransp) :: type_mpiptop
      integer, allocatable :: s_request(:)
      integer, allocatable :: r_request(:)
      integer, allocatable :: final_wait_array(:)
      complex(cp), allocatable :: temp_Rloc(:,:,:)
      integer, allocatable :: s_transfer_type_cont(:,:)
      integer, allocatable :: s_transfer_type_nr_end_cont(:,:)
      integer, allocatable :: r_transfer_type_cont(:,:)
      integer, allocatable :: r_transfer_type_nr_end_cont(:,:)
   contains
      procedure :: create_comm
      procedure :: destroy_comm
      procedure :: transp_lm2r
      procedure :: transp_r2lm
   end type type_mpiptop


contains

   subroutine initialize_comm(this, n_fields)

      type(type_mpiptop) :: this

      integer, intent(in) :: n_fields

#ifdef WITH_MPI
      integer :: proc, my_lm_per_rank, i, nR_main_ranks, nR_last_rank
      integer(kind=MPI_ADDRESS_KIND) :: zerolb, extent, sizeof_double_complex
      integer :: base_col_type,temptype
      integer :: blocklengths(n_fields),blocklengths_on_last(n_fields)
      integer :: displs(n_fields),displs_on_last(n_fields)


      if (.not. l_finite_diff ) then
         nR_main_ranks = (n_r_max-1)/n_procs
         nR_last_rank = nR_main_ranks+1
      else
         nR_main_ranks = n_r_max/n_procs
         nR_last_rank = nR_main_ranks+n_r_max-n_procs*nR_main_ranks
      end if

      allocate(this%s_transfer_type_cont(n_procs,n_fields))
      allocate(this%s_transfer_type_nr_end_cont(n_procs,n_fields))
      allocate(this%r_transfer_type_cont(n_procs,n_fields))
      allocate(this%r_transfer_type_nr_end_cont(n_procs,n_fields))
      bytes_allocated = bytes_allocated + 4*n_fields*n_procs*SIZEOF_INTEGER

      do proc=0,n_procs-1
         my_lm_per_rank=lm_balance(proc)%n_per_rank

         ! define the transfer types for the containers
         ! same schema as for the other types
         ! some temporary datatypes, not needed for communication
         ! but only for constructing the final datatypes
         call MPI_Type_get_extent(MPI_DEF_COMPLEX,zerolb,sizeof_double_complex,ierr)
         call MPI_Type_contiguous(my_lm_per_rank,MPI_DEF_COMPLEX,temptype,ierr)
         zerolb=0
         extent=lm_max*sizeof_double_complex
         call MPI_Type_create_resized(temptype,zerolb,extent,base_col_type,ierr)

         do i=1,n_fields
            blocklengths(i)         = nR_main_ranks
            displs(i)               = (i-1)*nR_main_ranks
            blocklengths_on_last(i) = nR_last_rank
            displs_on_last(i)       = (i-1)*nR_last_rank
         end do

         do i=1,n_fields
            call MPI_Type_vector(i,nR_main_ranks*my_lm_per_rank,        &
                 &               n_r_max*my_lm_per_rank,MPI_DEF_COMPLEX,&
                 &               this%r_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_commit(this%r_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_vector(i,nR_last_rank*my_lm_per_rank,            &
                 &               n_r_max*my_lm_per_rank,MPI_DEF_COMPLEX,   &
                 &               this%r_transfer_type_nr_end_cont(proc+1,i),ierr)
            call MPI_Type_commit(this%r_transfer_type_nr_end_cont(proc+1,i),ierr)

            call MPI_Type_indexed(i,blocklengths(1:i),displs(1:i),base_col_type, &
                 &                this%s_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_commit(this%s_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_indexed(i,blocklengths_on_last(1:i),             &
                 &                displs_on_last(1:i),base_col_type,       &
                 &                this%s_transfer_type_nr_end_cont(proc+1,i),ierr)
            call MPI_Type_commit(this%s_transfer_type_nr_end_cont(proc+1,i),ierr)
         end do

      end do
#endif

   end subroutine initialize_comm
!----------------------------------------------------------------------------------
   subroutine finalize_comm(this)

      type(type_mpiptop) :: this

#ifdef WITH_MPI
      deallocate( this%s_transfer_type_cont, this%s_transfer_type_nr_end_cont )
      deallocate( this%r_transfer_type_cont, this%r_transfer_type_nr_end_cont )
#endif

   end subroutine finalize_comm
!----------------------------------------------------------------------------------
   subroutine create_comm(this, n_fields)

      class(type_mpiptop) :: this
      integer, intent(in) :: n_fields

      this%n_fields = n_fields

#ifdef WITH_MPI
      allocate(this%s_request(n_procs-1))
      allocate(this%r_request(n_procs-1))
      allocate(this%final_wait_array(2*(n_procs-1)))
      bytes_allocated = bytes_allocated+4*(n_procs-1)*SIZEOF_INTEGER
#endif

      allocate(this%temp_Rloc(1:lm_max,nRstart:nRstop,1:this%n_fields))
      bytes_allocated = bytes_allocated+&
      &                 lm_max*(nRstop-nRstart+1)*this%n_fields*SIZEOF_DEF_COMPLEX

      call initialize_comm(this, n_fields)

   end subroutine create_comm
!----------------------------------------------------------------------------------
   subroutine destroy_comm(this)

      class(type_mpiptop) :: this

      call finalize_comm(this)

#ifdef WITH_MPI
      deallocate(this%final_wait_array, this%s_request, this%r_request)
#endif
      deallocate(this%temp_Rloc)

   end subroutine destroy_comm
!----------------------------------------------------------------------------------
   subroutine lm2r_redist_start(this,arr_LMloc,arr_Rloc)

      type(type_mpiptop) :: this
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
            if (rank == send_pe) then
               !PERFON('loc_copy')
               ! just copy
               do i=1,this%n_fields
                  arr_Rloc(llm:ulm,nRstart:nRstop,i)= &
                  &    arr_LMloc(llm:ulm,nRstart:nRstop,i)
               end do
               !PERFOFF
            else
               !PERFON('irecv')
               call MPI_Irecv(arr_Rloc(lm_balance(recv_pe)%nStart,nRstart,1),      &
                    &         1,this%s_transfer_type_cont(recv_pe+1,this%n_fields),&
                    &         recv_pe,transfer_tag,MPI_COMM_WORLD,                 &
                    &         this%r_request(irank),ierr)
               !PERFOFF
               !PERFON('isend')
               if (send_pe == n_procs-1) then
                  call MPI_Isend(arr_LMloc(llm,radial_balance(send_pe)%nStart,1),   &
                       &         1,this%r_transfer_type_nr_end_cont(rank+1,         &
                       &         this%n_fields),send_pe,transfer_tag,MPI_COMM_WORLD,&
                       &         this%s_request(irank),ierr)
               else
                  call MPI_Isend(arr_LMloc(llm,radial_balance(send_pe)%nStart,1),   &
                       &         1,this%r_transfer_type_cont(rank+1,this%n_fields), &
                       &         send_pe,transfer_tag,MPI_COMM_WORLD,               &
                       &         this%s_request(irank),ierr)
               end if
               !PERFOFF
            end if
         end do

         i=1
         do irank=1,n_procs-1
            this%final_wait_array(i)=this%s_request(irank)
            this%final_wait_array(i+1)=this%r_request(irank)
            i = i + 2
         end do
      else
         ! rank  ==  n_procs-1
         ! all receives are with the s_transfer_type_nr_end
         ! all sends are done with r_transfer_type_lm_end
         do irank=0,n_procs-1
            ! send_pe: send to this rank
            ! recv_pe: receive from this rank
            send_pe = modulo(rank+irank,n_procs)
            recv_pe = modulo(rank-irank+n_procs,n_procs)
            if (rank == send_pe) then
               !PERFON('loc_copy')
               ! just copy
               do i=1,this%n_fields
                  arr_Rloc(llm:ulm,nRstart:nRstop,i)= &
                  &     arr_LMloc(llm:ulm,nRstart:nRstop,i)
               end do
               !PERFOFF
            else
               !PERFON('irecv')
               call MPI_Irecv(arr_Rloc(lm_balance(recv_pe)%nStart,nRstart,1),1,  &
                    &         this%s_transfer_type_nr_end_cont(recv_pe+1,        &
                    &         this%n_fields),recv_pe,transfer_tag,MPI_COMM_WORLD,&
                    &         this%r_request(irank),ierr)
               !PERFOFF
               !PERFON('isend')
               call MPI_Isend(arr_LMloc(llm,radial_balance(send_pe)%nStart,1),   &
                    &         1,this%r_transfer_type_cont(rank+1,this%n_fields), &
                    &         send_pe,transfer_tag,MPI_COMM_WORLD,               &
                    &         this%s_request(irank),ierr)
               !PERFOFF
            end if
         end do
         i=1
         do irank=1,n_procs-1
            this%final_wait_array(i)=this%s_request(irank)
            this%final_wait_array(i+1)=this%r_request(irank)
            i = i + 2
         end do
      end if
      !PERFOFF

#else
      do i=1,this%n_fields
         arr_Rloc(llm:ulm,nRstart:nRstop,i)= arr_LMloc(llm:ulm,nRstart:nRstop,i)
      end do
#endif

   end subroutine lm2r_redist_start
!----------------------------------------------------------------------------------
   subroutine lm2r_redist_wait(this)

      type(type_mpiptop) :: this
#ifdef WITH_MPI
      integer :: ierr
      integer :: array_of_statuses(MPI_STATUS_SIZE,2*n_procs)

      call MPI_Waitall(2*(n_procs-1),this%final_wait_array,array_of_statuses,ierr)
#endif

   end subroutine lm2r_redist_wait
!----------------------------------------------------------------------------------
   subroutine lo2r_redist_start(this,arr_lo)

      type(type_mpiptop) :: this
      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max,*)

      call lm2r_redist_start(this,arr_lo,this%temp_Rloc)

   end subroutine lo2r_redist_start
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_wait(this,arr_Rloc)

      type(type_mpiptop) :: this

      !-- Output variable
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)

      ! Local variables
      integer :: nR, l, m, lm, i, start_lm, stop_lm

      !PERFON("lo2r_wt")
      call lm2r_redist_wait(this)

      !$omp parallel default(shared) private(start_lm,stop_lm,nR,i,lm,l,m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      ! now in this%temp_Rloc we do have the lo_ordered r-local part
      ! now reorder to the original ordering
      do i=1,this%n_fields
         do nR=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               m = st_map%lm2m(lm)
               arr_Rloc(lm,nR,i)=this%temp_Rloc(lo_map%lm2(l,m),nR,i)
            end do
         end do
      end do
      !$omp end parallel

      !PERFOFF

   end subroutine lo2r_redist_wait
!-------------------------------------------------------------------------------
   subroutine r2lo_redist_start(this,arr_Rloc,arr_lo)

      type(type_mpiptop) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max,*)

      ! Local variables
      integer :: nR,l,m,i,lm,start_lm,stop_lm

      !$omp parallel default(shared) private(start_lm,stop_lm,nR,i,lm,l,m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)

      ! Just copy the array with permutation
      !PERFON('r2lo_dst')
      do i=1,this%n_fields
         do nR=nRstart,nRstop
            do lm=start_lm,stop_lm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               this%temp_Rloc(lm,nR,i)=arr_Rloc(st_map%lm2(l,m),nR,i)
            end do
         end do
      end do
      !$omp end parallel

      call r2lm_redist_start(this,this%temp_Rloc,arr_lo)
      !PERFOFF

   end subroutine r2lo_redist_start
!-------------------------------------------------------------------------------
   subroutine r2lo_redist_wait(this)

      type(type_mpiptop) :: this

#ifdef WITH_MPI
      integer :: ierr
      integer :: array_of_statuses(MPI_STATUS_SIZE,2*n_procs)

      !PERFON('lm2r_wt')
      call MPI_Waitall(2*(n_procs-1),this%final_wait_array,array_of_statuses,ierr)
      !PERFOFF
#endif

   end subroutine r2lo_redist_wait
!-------------------------------------------------------------------------------
   subroutine r2lm_redist_start(this,arr_rloc,arr_LMloc)

      type(type_mpiptop) :: this
      complex(cp), intent(in) :: arr_Rloc(lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,n_r_max,*)

      integer :: i
#ifdef WITH_MPI
      ! Local variables
      integer :: send_pe, recv_pe, irank
      integer :: transfer_tag=1111

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
               do i=1,this%n_fields
                  arr_LMLoc(llm:ulm,nRstart:nRstop,i)= &
                  &         arr_Rloc(llm:ulm,nRstart:nRstop,i)
               end do
            else
               call MPI_Isend(arr_Rloc(lm_balance(send_pe)%nStart,nRstart,1),      &
                    &         1,this%s_transfer_type_cont(send_pe+1,this%n_fields),&
                    &         send_pe,transfer_tag,MPI_COMM_WORLD,                 &
                    &         this%s_request(irank),ierr)
               if (recv_pe == n_procs-1) then
                  call MPI_Irecv(arr_LMloc(llm,radial_balance(recv_pe)%nStart,1),   &
                       &         1,this%r_transfer_type_nr_end_cont(rank+1,         &
                       &         this%n_fields),recv_pe,transfer_tag,MPI_COMM_WORLD,&
                       &         this%r_request(irank),ierr)
               else
                  call MPI_Irecv(arr_LMloc(llm,radial_balance(recv_pe)%nStart,1),   &
                       &         1,this%r_transfer_type_cont(rank+1,this%n_fields), &
                       &         recv_pe,transfer_tag,MPI_COMM_WORLD,               &
                       &         this%r_request(irank),ierr)
               end if
            end if
         end do

         i=1
         do irank=1,n_procs-1
            this%final_wait_array(i)=this%s_request(irank)
            this%final_wait_array(i+1)=this%r_request(irank)
            i = i + 2
         end do
      else
         ! rank  ==  n_procs-1
         ! all receives are with the r_transfer_type_lm_end
         ! all sends are done with s_transfer_type_nr_end
         do irank=0,n_procs-1
            ! send_pe: send to this rank
            ! recv_pe: receive from this rank
            send_pe = modulo(rank+irank,n_procs)
            recv_pe = modulo(rank-irank+n_procs,n_procs)
            if (rank == send_pe) then
               ! just copy
               do i=1,this%n_fields
                  arr_LMLoc(llm:ulm,nRstart:nRstop,i)= &
                  &        arr_Rloc(llm:ulm,nRstart:nRstop,i)
               end do
            else
               call MPI_Irecv(arr_LMloc(llm,radial_balance(recv_pe)%nStart,1),   &
                    &         1,this%r_transfer_type_cont(rank+1,this%n_fields), &
                    &         recv_pe,transfer_tag,MPI_COMM_WORLD,               &
                    &         this%r_request(irank),ierr)
               call MPI_Isend(arr_Rloc(lm_balance(send_pe)%nStart,nRstart,1),    &
                    &         1,this%s_transfer_type_nr_end_cont(send_pe+1,      &
                    &         this%n_fields),send_pe,transfer_tag,MPI_COMM_WORLD,&
                    &         this%s_request(irank),ierr)

            end if
         end do
         i=1
         do irank=1,n_procs-1
            this%final_wait_array(i)=this%s_request(irank)
            this%final_wait_array(i+1)=this%r_request(irank)
            i = i + 2
         end do
      end if
      !PERFOFF
#else
      do i=1,this%n_fields
         arr_LMLoc(llm:ulm,nRstart:nRstop,i)=arr_Rloc(llm:ulm,nRstart:nRstop,i)
      end do
#endif

   end subroutine r2lm_redist_start
!-------------------------------------------------------------------------------
   subroutine transp_lm2r(this, arr_LMloc, arr_Rloc)
      !
      ! This subroutine transposes a LM-distributed container of arrays into
      ! a r-distributed container of arrays
      !

      class(type_mpiptop) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)

      call lo2r_redist_start(this, arr_LMloc)
      call lo2r_redist_wait(this, arr_Rloc)

   end subroutine transp_lm2r
!----------------------------------------------------------------------------------
   subroutine transp_r2lm(this, arr_Rloc, arr_LMloc)
      !
      ! This subroutine transposes a r-distributed container of arrays into
      ! a LM-distributed container of arrays
      !

      class(type_mpiptop) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)

      call r2lo_redist_start(this, arr_Rloc, arr_LMloc)
      call r2lo_redist_wait(this)

   end subroutine transp_r2lm
!----------------------------------------------------------------------------------
end module mpi_ptop_mod
