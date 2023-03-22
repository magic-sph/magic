module real_matrices

   use precision_mod
   use iso_c_binding

   implicit none

   private

   type, abstract, public :: type_realmat
      integer :: nrow ! Number of rows
      integer :: ncol ! Number of columns or number of bands
      logical :: l_pivot
      real(cp), pointer :: dat(:,:) ! Actual data
      integer,  pointer :: pivot(:)
      logical :: gpu_is_used
   contains
      procedure(initialize_if), deferred :: initialize
      procedure(finalize_if), deferred :: finalize
      procedure :: mat_add
      generic :: operator(+) => mat_add
      procedure(prepare_if), deferred :: prepare
      procedure(solve_real_multi_if), deferred :: solve_real_multi
      procedure(solve_real_single_if), deferred :: solve_real_single
      procedure(solve_complex_single_if), deferred :: solve_complex_single
      generic :: solve => solve_real_single, solve_complex_single, solve_real_multi
      procedure(set_data_if), deferred :: set_data
   end type type_realmat

   interface

#ifdef WITH_OMP_GPU
      subroutine initialize_if(this, nx, ny, l_pivot, use_gpu, nfull)
#else
      subroutine initialize_if(this, nx, ny, l_pivot, nfull)
#endif
         import
         class(type_realmat) :: this
         integer, intent(in) :: nx
         integer, intent(in) :: ny
         logical, intent(in) :: l_pivot
         integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
         logical, optional, intent(in) :: use_gpu
#endif
      end subroutine initialize_if

      subroutine finalize_if(this)
         import
         class(type_realmat) :: this
      end subroutine finalize_if

      subroutine prepare_if(this, info, handle, devInfo)
         import
         class(type_realmat) :: this
         integer, intent(out) :: info
         integer,        optional, intent(inout) :: devInfo(:)
         type(c_ptr),    optional, intent(inout) :: handle
      end subroutine prepare_if

      subroutine solve_real_multi_if(this, rhs, nRHS, handle, devInfo)
         import
         class(type_realmat) :: this
         integer,  intent(in) :: nRHS
         real(cp), intent(inout) :: rhs(:,:)
         integer,        optional, intent(inout) :: devInfo(:)
         type(c_ptr),    optional, intent(inout) :: handle
      end subroutine solve_real_multi_if

      subroutine solve_real_single_if(this, rhs, handle, devInfo)
         import
         class(type_realmat) :: this
         real(cp), intent(inout) :: rhs(:)
         type(c_ptr), optional,    intent(inout) :: handle
         integer,     optional,    intent(inout) :: devInfo(:)
      end subroutine solve_real_single_if

      subroutine solve_complex_single_if(this, rhs, tmpr, tmpi, handle, devInfo)
         import
         class(type_realmat) :: this
         complex(cp), intent(inout) :: rhs(:)
         real(cp), optional,    intent(inout) :: tmpr(:)
         real(cp), optional,    intent(inout) :: tmpi(:)
         type(c_ptr), optional, intent(inout) :: handle
         integer, optional,     intent(inout) :: devInfo(:)
      end subroutine solve_complex_single_if

      subroutine set_data_if(this, dat)
         import
         class(type_realmat) :: this
         real(cp), intent(in) :: dat(:,:)
      end subroutine set_data_if

   end interface

contains

   function mat_add(this, B)
      class(type_realmat), intent(in) :: this
      class(type_realmat), intent(in) :: B

      class(type_realmat), allocatable :: mat_add

!       if ( .not. allocated(mat_add%dat) ) then
!          call mat_add%initialize(this%nrow,this%ncol,this%l_pivot)
!       end if
!
!       mat_add%dat(:,:) = this%dat(:,:)+B%dat(:,:)

   end function mat_add
!------------------------------------------------------------------------------
end module real_matrices

module real_many_matrices

   use precision_mod
   use iso_c_binding

   implicit none

   private

   type, abstract, public :: type_mrealmat
      integer :: nrow ! Number of rows
      integer :: ncol ! Number of columns or number of bands
      integer :: nmat ! Number of matrices
      logical :: l_pivot
      real(cp), pointer :: dat(:,:,:) ! Actual data
      integer,  pointer :: pivot(:,:)
      logical :: gpu_is_used
   contains
      procedure(initialize_if), deferred :: initialize
      procedure(finalize_if), deferred :: finalize
      procedure(prepare_single_if), deferred :: prepare_single
      procedure(prepare_all_if), deferred :: prepare_all
      procedure(solve_real_multi_if), deferred :: solve_real_multi
      procedure(solve_real_single_if), deferred :: solve_real_single
      procedure(solve_complex_single_if), deferred :: solve_complex_single
      procedure(solve_all_mats_if), deferred :: solve_all_mats
      generic :: solve => solve_real_single, solve_complex_single, &
      &                   solve_real_multi, solve_all_mats
      generic :: prepare => prepare_single, prepare_all
      procedure(set_data_if), deferred :: set_data
   end type type_mrealmat

   interface

#ifdef WITH_OMP_GPU
      subroutine initialize_if(this, nx, ny, nmat, l_pivot, use_gpu, nfull)
#else
      subroutine initialize_if(this, nx, ny, nmat, l_pivot, nfull)
#endif
         import
         class(type_mrealmat) :: this
         integer, intent(in) :: nx
         integer, intent(in) :: ny
         integer, intent(in) :: nmat
         logical, intent(in) :: l_pivot
         integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
         logical, optional, intent(in) :: use_gpu
#endif
      end subroutine initialize_if

      subroutine finalize_if(this)
         import
         class(type_mrealmat) :: this
      end subroutine finalize_if

      subroutine prepare_single_if(this, idx, info, handle, devInfo)
         import
         class(type_mrealmat) :: this
         integer, intent(in) :: idx
         integer, intent(out) :: info
         integer,     optional, intent(inout) :: devInfo(:)
         type(c_ptr), optional, intent(inout) :: handle
      end subroutine prepare_single_if

      subroutine prepare_all_if(this, info)
         import
         class(type_mrealmat) :: this
         integer, intent(out) :: info
      end subroutine prepare_all_if

      subroutine solve_real_multi_if(this, rhs, nRHS, idx, handle, devInfo)
         import
         class(type_mrealmat) :: this
         integer,     intent(in) :: nRHS
         integer,     intent(in) :: idx
         real(cp),    intent(inout) :: rhs(:,:)
         integer,     optional, intent(inout) :: devInfo(:)
         type(c_ptr), optional, intent(inout) :: handle
      end subroutine solve_real_multi_if

      subroutine solve_real_single_if(this, rhs, handle, devInfo)
         import
         class(type_mrealmat) :: this
         real(cp),    intent(inout) :: rhs(:)
         type(c_ptr), optional, intent(inout) :: handle
         integer,     optional, intent(inout) :: devInfo(:)
      end subroutine solve_real_single_if

      subroutine solve_complex_single_if(this, rhs, tmpr, tmpi, handle, devInfo)
         import
         class(type_mrealmat) :: this
         complex(cp), intent(inout) :: rhs(:)
         real(cp), optional,    intent(inout) :: tmpr(:)
         real(cp), optional,    intent(inout) :: tmpi(:)
         type(c_ptr), optional, intent(inout) :: handle
         integer, optional,     intent(inout) :: devInfo(:)
      end subroutine solve_complex_single_if

      subroutine set_data_if(this, dat, idx)
         import
         class(type_mrealmat) :: this
         integer,  intent(in) :: idx
         real(cp), intent(in) :: dat(:,:)
      end subroutine set_data_if

      subroutine solve_all_mats_if(this, rhs, llm, ulm, lm2l, l2nLMB2)
         import
         class(type_mrealmat) :: this
         integer, intent(in) :: llm
         integer, intent(in) :: ulm
         integer, intent(in) :: lm2l(:)
         integer, intent(in) :: l2nLMB2(:)
         complex(cp), intent(inout) :: rhs(this%ncol,llm:ulm)
      end subroutine solve_all_mats_if

   end interface

end module real_many_matrices

module dense_matrices

   use precision_mod
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc
#endif
   use real_matrices, only: type_realmat
   use real_many_matrices, only: type_mrealmat
   use constants, only: one
   use algebra, only: solve_mat, prepare_mat
#ifdef WITH_OMP_GPU
   use algebra_hipfort, only: gpu_solve_mat, gpu_prepare_mat
#endif
   use iso_c_binding

   implicit none

   type, public, extends(type_realmat) :: type_densemat
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: prepare
      procedure :: solve_real_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
   end type type_densemat

   type, public, extends(type_mrealmat) :: type_mdensemat
#ifdef WITH_OMP_GPU
      integer, allocatable :: info_array(:,:)
#endif
   contains
      procedure :: initialize => initialize_
      procedure :: finalize => finalize_
      procedure :: prepare_single
      procedure :: prepare_all
      procedure :: solve_real_multi => solve_real_multi_
      procedure :: solve_real_single => solve_real_single_
      procedure :: solve_complex_single => solve_complex_single_
      procedure :: set_data => set_data_
      procedure :: solve_all_mats
   end type type_mdensemat

contains

#ifdef WITH_OMP_GPU
   subroutine initialize(this, nx, ny, l_pivot, use_gpu, nfull)
#else
   subroutine initialize(this, nx, ny, l_pivot, nfull)
#endif
      !
      ! Memory allocation
      !
      class(type_densemat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
      logical, optional, intent(in) :: use_gpu
#endif

      !--
      logical :: loc_use_gpu
      loc_use_gpu = .false.
      this%gpu_is_used=.false.
#ifdef WITH_OMP_GPU
      if ( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      end if
      if(loc_use_gpu) then
         this%gpu_is_used = .true.
      end if
#endif

      this%nrow = nx
      this%ncol = ny
      this%l_pivot = l_pivot
      allocate( this%dat(nx, ny) )
      this%dat(:,:) = 0.0_cp
      bytes_allocated = bytes_allocated+nx*ny*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
      if ( loc_use_gpu) then
         !$omp target enter data map(to : this%dat)
         gpu_bytes_allocated = gpu_bytes_allocated+nx*ny*SIZEOF_DEF_REAL
      end if
#endif

      if ( this%l_pivot ) then
         allocate( this%pivot(this%nrow) )
         this%pivot(:) = 0
         bytes_allocated = bytes_allocated+this%nrow*SIZEOF_INTEGER
#ifdef WITH_OMP_GPU
         if ( loc_use_gpu) then
            !$omp target enter data map(to : this%pivot)
            gpu_bytes_allocated = gpu_bytes_allocated+this%nrow*SIZEOF_INTEGER
         end if
#endif
      end if

   end subroutine initialize
!------------------------------------------------------------------------------
#ifdef WITH_OMP_GPU
   subroutine initialize_(this, nx, ny, nmat, l_pivot, use_gpu, nfull)
#else
   subroutine initialize_(this, nx, ny, nmat, l_pivot, nfull)
#endif
      !
      ! Memory allocation
      !
      class(type_mdensemat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      integer, intent(in) :: nmat
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
      logical, optional, intent(in) :: use_gpu
#endif

      !--
      logical :: loc_use_gpu
      loc_use_gpu = .false.
      this%gpu_is_used=.false.
#ifdef WITH_OMP_GPU
      if ( present(use_gpu) ) then
         loc_use_gpu = use_gpu
      end if
      if(loc_use_gpu) then
         this%gpu_is_used = .true.
      end if
#endif

      this%nrow = nx
      this%ncol = ny
      this%nmat = nmat
      this%l_pivot = l_pivot
      allocate( this%dat(nx, ny, nmat) )
      this%dat(:,:,:) = 0.0_cp
      bytes_allocated = bytes_allocated+nx*ny*nmat*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
      if ( loc_use_gpu) then
         gpu_bytes_allocated = gpu_bytes_allocated+nx*ny*nmat*SIZEOF_DEF_REAL
      end if
#endif

      if ( this%l_pivot ) then
         allocate( this%pivot(this%nrow, this%nmat) )
         this%pivot(:,:) = 0
         bytes_allocated = bytes_allocated+this%nrow*this%nmat*SIZEOF_INTEGER
#ifdef WITH_OMP_GPU
         if ( loc_use_gpu) then
            gpu_bytes_allocated = gpu_bytes_allocated+this%nrow*this%nmat*SIZEOF_INTEGER
         end if
#endif
      end if

#ifdef WITH_OMP_GPU
      if ( loc_use_gpu) then
         allocate( this%info_array(this%nrow, this%nmat) )
         this%info_array(:,:) = 0
         !$omp target enter data map(to : this)
      end if
#endif

   end subroutine initialize_
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(type_densemat) :: this

#ifdef WITH_OMP_GPU
      if ( this%gpu_is_used ) then
         !$omp target exit data map(delete : this%dat)
      end if
#endif
      deallocate( this%dat )

      if ( this%l_pivot ) then
#ifdef WITH_OMP_GPU
         if ( this%gpu_is_used ) then
            !$omp target exit data map(delete : this%pivot)
         end if
#endif
         deallocate (this%pivot)
      end if

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine finalize_(this)
      !
      ! Memory deallocation
      !
      class(type_mdensemat) :: this

#ifdef WITH_OMP_GPU
      if ( this%gpu_is_used ) then
         !$omp target exit data map(release : this)
      end if
#endif

      deallocate( this%dat )

      if ( this%l_pivot ) then
         deallocate (this%pivot)
      end if

#ifdef WITH_OMP_GPU
      if ( this%gpu_is_used ) then
         deallocate( this%info_array )
      end if
#endif

   end subroutine finalize_
!------------------------------------------------------------------------------
   subroutine prepare(this, info, handle, devInfo)

      class(type_densemat) :: this
      integer, intent(out) :: info
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         call gpu_prepare_mat(this%dat, this%nrow, this%nrow, this%pivot, info, &
              &               handle, devInfo)
#endif
      else
         call prepare_mat(this%dat, this%nrow, this%nrow, this%pivot, info)
      end if

   end subroutine prepare
!------------------------------------------------------------------------------
   subroutine prepare_single(this, idx, info, handle, devInfo)

      class(type_mdensemat) :: this
      integer, intent(in) :: idx
      integer, intent(out) :: info
      integer,     optional, intent(inout) :: devInfo(:)
      type(c_ptr), optional, intent(inout) :: handle

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         call gpu_prepare_mat(this%dat(:,:,idx), this%nrow, this%nrow, &
              &               this%pivot(:,idx), info, handle, devInfo)
#endif
      else
         call prepare_mat(this%dat(:,:,idx), this%nrow, this%nrow, &
              &           this%pivot(:,idx), info)
      end if

   end subroutine prepare_single
!------------------------------------------------------------------------------
   subroutine prepare_all(this, info)
      !
      ! LU decomposition for all the real matrix a(n,n,ell) via Gaussian elimination
      !

      class(type_mdensemat) :: this
      integer, intent(out) :: info ! Output diagnostic of success

      !-- Local variables:
      integer :: nm1,k,kp1,l,i,j,idx
      real(cp) :: help

      info=0
      nm1 =this%nrow-1

      !-- This external loop should be put on GPU
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(k,kp1,l,i,j,help)
#endif
      do idx=1,this%nmat
         do k=1,nm1

            kp1=k+1
            l  =k

            do i=kp1,this%nrow
               if ( abs(this%dat(i,k,idx)) > abs(this%dat(l,k,idx)) ) l=i
            end do
            this%pivot(k,idx)=l

            if ( abs(this%dat(l,k,idx)) > 10.0_cp*epsilon(0.0_cp) ) then
               if ( l /= k ) then
                  do i=1,this%nrow
                     help             =this%dat(k,i,idx)
                     this%dat(k,i,idx)=this%dat(l,i,idx)
                     this%dat(l,i,idx)=help
                  end do
               end if

               help=one/this%dat(k,k,idx)
               do i=kp1,this%nrow
                  this%dat(i,k,idx)=help*this%dat(i,k,idx)
               end do

               do j=kp1,this%nrow
                  do i=kp1,this%nrow
                     this%dat(i,j,idx)=this%dat(i,j,idx)-this%dat(k,j,idx)* &
                     &                 this%dat(i,k,idx)
                  end do
               end do
            else
#ifdef WITH_OMP_GPU
               this%info_array(k, idx)=k
#else
               info=k
#endif
            end if

         end do

         this%pivot(this%nrow,idx)=this%nrow
         if ( abs(this%dat(this%nrow,this%nrow,idx)) <= 10.0_cp*epsilon(0.0_cp) ) &
#ifdef WITH_OMP_GPU
          &   this%info_array(this%nrow,idx)=this%nrow
#else
          &   info=this%nrow
          if ( info > 0 ) return
#endif

         do i=1,this%nrow
            this%dat(i,i,idx)=one/this%dat(i,i,idx)
         end do

      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
      !$omp target teams distribute parallel do collapse(2) reduction(+: info)
      do idx=1,this%nmat
         do k=1,this%nrow
            if ( this%info_array(k, idx) > 0 ) then
               info = info + 1
            end if
         end do
      end do
      !$omp end target teams distribute parallel do
#endif

   end subroutine prepare_all
!------------------------------------------------------------------------------
   subroutine solve_real_single(this, rhs, handle, devInfo)

      class(type_densemat) :: this
      real(cp),                 intent(inout) :: rhs(:)
      type(c_ptr), optional,    intent(inout) :: handle
      integer, optional,        intent(inout) :: devInfo(:)

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         call gpu_solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs, &
              &             handle, devInfo)
#endif
      else
         call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs)
      end if

   end subroutine solve_real_single
!------------------------------------------------------------------------------
   subroutine solve_real_single_(this, rhs, handle, devInfo)

      class(type_mdensemat) :: this
      real(cp),                 intent(inout) :: rhs(:)
      type(c_ptr), optional,    intent(inout) :: handle
      integer, optional,        intent(inout) :: devInfo(:)

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         call gpu_solve_mat(this%dat(:,:,1), this%nrow, this%nrow, this%pivot(:,1), rhs, &
              &             handle, devInfo)
#endif
      else
         call solve_mat(this%dat(:,:,1), this%nrow, this%nrow, this%pivot(:,1), rhs)
      end if

   end subroutine solve_real_single_
!------------------------------------------------------------------------------
   subroutine solve_complex_single(this, rhs, tmpr, tmpi, handle, devInfo)

      class(type_densemat) :: this
      complex(cp), intent(inout) :: rhs(:)

      !-- Output variables
      real(cp), optional,    intent(inout) :: tmpr(:)
      real(cp), optional,    intent(inout) :: tmpi(:)
      type(c_ptr), optional, intent(inout) :: handle
      integer, optional,     intent(inout) :: devInfo(:)

      !--
#ifdef WITH_OMP_GPU
      integer :: n, i
      n =this%nrow
#endif

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         !-- Extract real and imag parts of input rhs matrix
         !$omp target teams distribute parallel do
         do i=1,n
            tmpr(i) = real(rhs(i))
            tmpi(i) = aimag(rhs(i))
         end do
         !$omp end target teams distribute parallel do
         call gpu_solve_mat(this%dat, this%nrow, this%nrow, this%pivot, tmpr, &
              &             tmpi, handle, devInfo)
         !$omp target teams distribute parallel do
         do i=1,n
            rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
         end do
         !$omp end target teams distribute parallel do
#endif
      else
         call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs)
      end if

   end subroutine solve_complex_single
!------------------------------------------------------------------------------
   subroutine solve_complex_single_(this, rhs, tmpr, tmpi, handle, devInfo)

      class(type_mdensemat) :: this
      complex(cp), intent(inout) :: rhs(:)

      !-- Output variables
      real(cp), optional,    intent(inout) :: tmpr(:)
      real(cp), optional,    intent(inout) :: tmpi(:)
      type(c_ptr), optional, intent(inout) :: handle
      integer, optional,     intent(inout) :: devInfo(:)

      !--
#ifdef WITH_OMP_GPU
      integer :: n, i
      n =this%nrow
#endif

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         !-- Extract real and imag parts of input rhs matrix
         !$omp target teams distribute parallel do
         do i=1,n
            tmpr(i) = real(rhs(i))
            tmpi(i) = aimag(rhs(i))
         end do
         !$omp end target teams distribute parallel do
         call gpu_solve_mat(this%dat(:,:,1), this%nrow, this%nrow, this%pivot(:,1), &
              &             tmpr, tmpi, handle, devInfo)
         !$omp target teams distribute parallel do
         do i=1,n
            rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
         end do
         !$omp end target teams distribute parallel do
#endif
      else
         call solve_mat(this%dat(:,:,1), this%nrow, this%nrow, this%pivot(:,1), rhs)
      end if

   end subroutine solve_complex_single_
!------------------------------------------------------------------------------
   subroutine solve_real_multi(this, rhs, nRHS, handle, devInfo)

      class(type_densemat) :: this
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         call gpu_solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs, &
              &             nRHS, handle, devInfo)
#endif
      else
         call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs, nRHS)
      end if

   end subroutine solve_real_multi
!------------------------------------------------------------------------------
   subroutine solve_real_multi_(this, rhs, nRHS, idx, handle, devInfo)

      class(type_mdensemat) :: this
      integer,  intent(in) :: idx
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         call gpu_solve_mat(this%dat(:,:,idx), this%nrow, this%nrow, &
              &             this%pivot(:,idx), rhs, nRHS, handle, devInfo)
#endif
      else
         call solve_mat(this%dat(:,:,idx), this%nrow, this%nrow, this%pivot(:,idx), &
              &         rhs, nRHS)
      end if

   end subroutine solve_real_multi_
!------------------------------------------------------------------------------
   subroutine set_data(this, dat)

      class(type_densemat) :: this
      real(cp), intent(in) :: dat(:,:)

#ifdef WITH_OMP_GPU
      integer :: i,j, row, col
      real(cp), pointer :: ptr_dat(:,:)
      ptr_dat => this%dat
      row = this%nrow
      col = this%ncol
#endif

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
         do j=1,col
            do i=1,row
               ptr_dat(i,j) = dat(i,j)
            end do
         end do
         !$omp end target teams distribute parallel do
#endif
      else
         this%dat(:,:) = dat(:,:)
      end if

   end subroutine set_data
!------------------------------------------------------------------------------
   subroutine set_data_(this, dat, idx)

      class(type_mdensemat) :: this
      integer,  intent(in) :: idx
      real(cp), intent(in) :: dat(:,:)

#ifdef WITH_OMP_GPU
      integer :: i,j, row, col
      real(cp), pointer :: ptr_dat(:,:)
      ptr_dat => this%dat(:,:,idx)
      row = this%nrow
      col = this%ncol
#endif

      if ( this%gpu_is_used ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
         do j=1,col
            do i=1,row
               ptr_dat(i,j) = dat(i,j)
            end do
         end do
         !$omp end target teams distribute parallel do
#endif
      else
         this%dat(:,:,idx) = dat(:,:)
      end if

   end subroutine set_data_
!------------------------------------------------------------------------------
   subroutine solve_all_mats(this, rhs, llm, ulm, lm2l, l2nLMB2)

      !use blocking, only: lo_map, lo_sub_map

      class(type_mdensemat) :: this

      !-- Input variables
      integer, intent(in) :: llm
      integer, intent(in) :: ulm
      integer, intent(in) :: lm2l(:)
      integer, intent(in) :: l2nLMB2(0:)

      !-- In/Out variables
      complex(cp), intent(inout) :: rhs(1:this%ncol,llm:ulm)

      !-- Local variables:
      integer :: nm1,nodd,i,m
      integer :: k,k1,nRHS,nLMB2,l,n
      complex(cp) :: help
      integer, pointer :: ptr_pivot(:,:)
      real(cp), pointer :: ptr_dat(:,:,:)

      ptr_dat   => this%dat
      ptr_pivot => this%pivot

      n = this%ncol

      nm1 =n-1
      nodd=mod(n,2)

      !-- Single loop over lm's
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do private(i,m,k,k1,nLMB2,l,help)
#endif
      do nRHS=llm,ulm

         l=lm2l(nRHS)
         nLMB2=l2nLMB2(l)

         !-- Permute vectors rhs
         do k=1,nm1
            m=ptr_pivot(k,nLMB2)
            help       =rhs(m,nRHS)
            rhs(m,nRHS) =rhs(k,nRHS)
            rhs(k,nRHS) =help
         end do

         !-- Solve  l * y = b
         do k=1,n-2,2
            k1=k+1
            rhs(k1,nRHS) =rhs(k1,nRHS)-rhs(k,nRHS)*ptr_dat(k1,k,nLMB2)
            !DIR$ CONCURRENT
            do i=k+2,n
               rhs(i,nRHS)=rhs(i,nRHS)-(rhs(k,nRHS)*ptr_dat(i,k,nLMB2) + &
               &                      rhs(k1,nRHS)*ptr_dat(i,k1,nLMB2))
            end do
         end do
         if ( nodd == 0 ) then
            rhs(n,nRHS) =rhs(n,nRHS)-rhs(nm1,nRHS)*ptr_dat(n,nm1,nLMB2)
         end if

         !-- Solve  u * x = y
         do k=n,3,-2
            k1=k-1
            rhs(k,nRHS)  =rhs(k,nRHS)*ptr_dat(k,k,nLMB2)
            rhs(k1,nRHS) =(rhs(k1,nRHS)-rhs(k,nRHS)*ptr_dat(k1,k,nLMB2)) * &
            &            ptr_dat(k1,k1,nLMB2)
            !DIR$ CONCURRENT
            do i=1,k-2
               rhs(i,nRHS)=rhs(i,nRHS)-rhs(k,nRHS)*ptr_dat(i,k,nLMB2) - &
               &          rhs(k1,nRHS)*ptr_dat(i,k1,nLMB2)
            end do
         end do
         if ( nodd == 0 ) then
            rhs(2,nRHS)=rhs(2,nRHS)*ptr_dat(2,2,nLMB2)
            rhs(1,nRHS)=(rhs(1,nRHS)-rhs(2,nRHS)*ptr_dat(1,2,nLMB2))*ptr_dat(1,1,nLMB2)
         else
            rhs(1,nRHS)=rhs(1,nRHS)*ptr_dat(1,1,nLMB2)
         end if

      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

   end subroutine solve_all_mats
!------------------------------------------------------------------------------
end module dense_matrices

module band_matrices

   use precision_mod
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc
#endif
   use real_matrices, only: type_realmat
   use real_many_matrices, only: type_mrealmat
   use algebra, only: solve_tridiag, prepare_tridiag, prepare_band, solve_band
   use iso_c_binding

   implicit none

   type, public, extends(type_realmat) :: type_bandmat
      real(cp), pointer :: du2(:)
      integer :: kl
      integer :: ku
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: prepare
      procedure :: solve_real_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
!      procedure :: mat_add
!      generic :: operator(+) => mat_add
   end type type_bandmat

   type, public, extends(type_mrealmat) :: type_mbandmat
      real(cp), pointer :: du2(:,:)
      integer :: kl
      integer :: ku
   contains
      procedure :: initialize => initialize_
      procedure :: finalize => finalize_
      procedure :: prepare_single
      procedure :: prepare_all
      procedure :: solve_real_multi => solve_real_multi_
      procedure :: solve_real_single => solve_real_single_
      procedure :: solve_complex_single => solve_complex_single_
      procedure :: set_data => set_data_
      procedure :: solve_all_mats
!      procedure :: mat_add
!      generic :: operator(+) => mat_add
   end type type_mbandmat

contains

#ifdef WITH_OMP_GPU
   subroutine initialize(this, nx, ny, l_pivot, use_gpu, nfull)
#else
   subroutine initialize(this, nx, ny, l_pivot, nfull)
#endif
      !
      ! Memory allocation
      !
      class(type_bandmat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
      logical, optional, intent(in) :: use_gpu
#endif

      !--
      logical :: loc_use_gpu
      loc_use_gpu = .false.
      this%gpu_is_used=.false.

      this%nrow = nx
      this%ncol = ny
      this%l_pivot = l_pivot

      this%kl = (nx-1)/2
      this%ku = this%kl

      if ( nx > 3 .and. this%l_pivot ) then
         allocate( this%dat(nx+(nx-1)/2, ny) )
         this%dat(:,:) = 0.0_cp
         bytes_allocated = bytes_allocated+(nx+(nx-1)/2)*ny*SIZEOF_DEF_REAL
      else
         allocate( this%dat(nx, ny) )
         this%dat(:,:) = 0.0_cp
         bytes_allocated = bytes_allocated+nx*ny*SIZEOF_DEF_REAL
      end if
      if ( this%l_pivot ) then
         allocate( this%pivot(this%ncol) )
         this%pivot(:) = 0
         bytes_allocated = bytes_allocated+this%ncol*SIZEOF_INTEGER
         if ( nx == 3 ) then ! Only require for tridiag arrays
            allocate( this%du2(this%ncol-2) ) ! Help array for tridiag
            this%du2(:) = 0
            bytes_allocated = bytes_allocated+(this%ncol-2)*SIZEOF_DEF_REAL
         end if
      end if

   end subroutine initialize
!------------------------------------------------------------------------------
#ifdef WITH_OMP_GPU
   subroutine initialize_(this, nx, ny, nmat, l_pivot, use_gpu, nfull)
#else
   subroutine initialize_(this, nx, ny, nmat, l_pivot, nfull)
#endif
      !
      ! Memory allocation
      !
      class(type_mbandmat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      integer, intent(in) :: nmat
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
      logical, optional, intent(in) :: use_gpu
#endif

      !--
      logical :: loc_use_gpu
      loc_use_gpu = .false.
      this%gpu_is_used=.false.

      this%nrow = nx
      this%ncol = ny
      this%nmat = nmat
      this%l_pivot = l_pivot

      this%kl = (nx-1)/2
      this%ku = this%kl

      if ( nx > 3 .and. this%l_pivot ) then
         allocate( this%dat(nx+(nx-1)/2, ny, nmat) )
         this%dat(:,:,:) = 0.0_cp
         bytes_allocated = bytes_allocated+(nx+(nx-1)/2)*ny*nmat*SIZEOF_DEF_REAL
      else
         allocate( this%dat(nx, ny, nmat) )
         this%dat(:,:,:) = 0.0_cp
         bytes_allocated = bytes_allocated+nx*ny*nmat*SIZEOF_DEF_REAL
      end if
      if ( this%l_pivot ) then
         allocate( this%pivot(this%ncol,this%nmat) )
         this%pivot(:,:) = 0
         bytes_allocated = bytes_allocated+this%ncol*this%nmat*SIZEOF_INTEGER
         if ( nx == 3 ) then ! Only require for tridiag arrays
            allocate( this%du2(this%ncol-2,this%nmat) ) ! Help array for tridiag
            this%du2(:,:) = 0
            bytes_allocated = bytes_allocated+(this%ncol-2)*this%nmat*SIZEOF_DEF_REAL
         end if
      end if

   end subroutine initialize_
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(type_bandmat) :: this

      deallocate( this%dat )

      if ( this%l_pivot ) then
         deallocate (this%pivot)
         if ( this%nrow == 3 ) then
            deallocate(this%du2)
         end if
      end if

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine finalize_(this)
      !
      ! Memory deallocation
      !
      class(type_mbandmat) :: this

      deallocate( this%dat )

      if ( this%l_pivot ) then
         deallocate (this%pivot)
         if ( this%nrow == 3 ) then
            deallocate(this%du2)
         end if
      end if

   end subroutine finalize_
!------------------------------------------------------------------------------
   subroutine prepare(this, info, handle, devInfo)

      class(type_bandmat) :: this
      integer, intent(out) :: info
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%nrow == 3 ) then
         call prepare_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:),  &
              &               this%dat(1,2:), this%du2, this%ncol,       &
              &               this%pivot, info)
      else
         call prepare_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, info)
      end if

   end subroutine prepare
!------------------------------------------------------------------------------
   subroutine prepare_single(this, idx, info, handle, devInfo)

      class(type_mbandmat) :: this
      integer, intent(in) :: idx
      integer, intent(out) :: info
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%nrow == 3 ) then
         call prepare_tridiag(this%dat(3,1:this%ncol-1,idx), this%dat(2,:,idx),  &
              &               this%dat(1,2:,idx), this%du2(:,idx), this%ncol,    &
              &               this%pivot(:,idx), info)
      else
         call prepare_band(this%dat(:,:,idx), this%ncol, this%kl, this%ku, &
              &            this%pivot(:,idx), info)
      end if

   end subroutine prepare_single
!------------------------------------------------------------------------------
   subroutine prepare_all(this, info)

      class(type_mbandmat) :: this
      integer, intent(out) :: info

      print*, 'Not implemented at this stage'

   end subroutine prepare_all
!------------------------------------------------------------------------------
   subroutine solve_real_single(this, rhs, handle, devInfo)

      class(type_bandmat) :: this
      real(cp), intent(inout) :: rhs(:)
      type(c_ptr), optional,    intent(inout) :: handle
      integer, optional,        intent(inout) :: devInfo(:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:), &
              &             this%dat(1,2:), this%du2, this%ncol,      &
              &             this%pivot, rhs)
      else
         call solve_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, rhs)
      end if

   end subroutine solve_real_single
!------------------------------------------------------------------------------
   subroutine solve_real_single_(this, rhs, handle, devInfo)

      class(type_mbandmat) :: this
      real(cp), intent(inout) :: rhs(:)
      type(c_ptr), optional,    intent(inout) :: handle
      integer, optional,        intent(inout) :: devInfo(:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1,1), this%dat(2,:,1), &
              &             this%dat(1,2:,1), this%du2(:,1), this%ncol,   &
              &             this%pivot(:,1), rhs)
      else
         call solve_band(this%dat(:,:,1), this%ncol, this%kl, this%ku, &
              &          this%pivot(:,1), rhs)
      end if

   end subroutine solve_real_single_
!------------------------------------------------------------------------------
   subroutine solve_complex_single(this, rhs, tmpr, tmpi, handle, devInfo)

      class(type_bandmat) :: this
      complex(cp), intent(inout) :: rhs(:)

      !-- Output variables
      real(cp), optional,    intent(inout) :: tmpr(:)
      real(cp), optional,    intent(inout) :: tmpi(:)
      type(c_ptr), optional, intent(inout) :: handle
      integer, optional,     intent(inout) :: devInfo(:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:), &
              &             this%dat(1,2:), this%du2, this%ncol,      &
              &             this%pivot, rhs)
      else
         call solve_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, rhs)
      end if

   end subroutine solve_complex_single
!------------------------------------------------------------------------------
   subroutine solve_complex_single_(this, rhs, tmpr, tmpi, handle, devInfo)

      class(type_mbandmat) :: this
      complex(cp), intent(inout) :: rhs(:)

      !-- Output variables
      real(cp), optional,    intent(inout) :: tmpr(:)
      real(cp), optional,    intent(inout) :: tmpi(:)
      type(c_ptr), optional, intent(inout) :: handle
      integer, optional,     intent(inout) :: devInfo(:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1,1), this%dat(2,:,1), &
              &             this%dat(1,2:,1), this%du2(:,1), this%ncol,   &
              &             this%pivot(:,1), rhs)
      else
         call solve_band(this%dat(:,:,1), this%ncol, this%kl, this%ku, &
              &          this%pivot(:,1), rhs)
      end if

   end subroutine solve_complex_single_
!------------------------------------------------------------------------------
   subroutine solve_real_multi(this, rhs, nRHS, handle, devInfo)

      class(type_bandmat) :: this
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:),   &
              &             this%dat(1,2:), this%du2, this%ncol,        &
              &             this%pivot, rhs, nRHS)
      else
         call solve_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, &
              &          rhs, nRHS)
      end if

   end subroutine solve_real_multi
!------------------------------------------------------------------------------
   subroutine solve_real_multi_(this, rhs, nRHS, idx, handle, devInfo)

      class(type_mbandmat) :: this
      integer,  intent(in) :: nRHS
      integer,  intent(in) :: idx
      real(cp), intent(inout) :: rhs(:,:)
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1,idx), this%dat(2,:,idx),   &
              &             this%dat(1,2:,idx), this%du2(:,idx), this%ncol,     &
              &             this%pivot(:,idx), rhs, nRHS)
      else
         call solve_band(this%dat(:,:,idx), this%ncol, this%kl, this%ku, &
              &          this%pivot(:,idx), rhs, nRHS)
      end if

   end subroutine solve_real_multi_
!------------------------------------------------------------------------------
   subroutine set_data(this, dat)

      class(type_bandmat) :: this

      !-- Input array
      real(cp), intent(in) :: dat(:,:)

      !-- Local variables
      integer :: i,j

      if ( this%nrow == 3 ) then
         do j=1,this%ncol
            do i=max(1,j-this%ku),min(this%ncol,j+this%kl)
               this%dat(this%ku+1+i-j,j)=dat(i,j)
            end do
         end do
      else
         do j=1,this%ncol
            do i=max(1,j-this%ku),min(this%ncol,j+this%kl)
               this%dat(this%kl+this%ku+1+i-j,j)=dat(i,j)
            end do
         end do
      end if

   end subroutine set_data
!------------------------------------------------------------------------------
   subroutine set_data_(this, dat, idx)

      class(type_mbandmat) :: this
      
      !-- Input variables
      integer,  intent(in) :: idx
      real(cp), intent(in) :: dat(:,:)

      !-- Local variables
      integer :: i,j

      if ( this%nrow == 3 ) then
         do j=1,this%ncol
            do i=max(1,j-this%ku),min(this%ncol,j+this%kl)
               this%dat(this%ku+1+i-j,j,idx)=dat(i,j)
            end do
         end do
      else
         do j=1,this%ncol
            do i=max(1,j-this%ku),min(this%ncol,j+this%kl)
               this%dat(this%kl+this%ku+1+i-j,j,idx)=dat(i,j)
            end do
         end do
      end if

   end subroutine set_data_
!------------------------------------------------------------------------------
   subroutine solve_all_mats(this, rhs, llm, ulm, lm2l, l2nLMB2)

      class(type_mbandmat) :: this

      !-- Input variables
      !-- Input variables
      integer, intent(in) :: llm
      integer, intent(in) :: ulm
      integer, intent(in) :: lm2l(:)
      integer, intent(in) :: l2nLMB2(0:)

      !-- In/Out variables
      complex(cp), intent(inout) :: rhs(this%ncol, llm:ulm)

      print*, 'Implementation to be completed...'

   end subroutine solve_all_mats
!------------------------------------------------------------------------------
!   function mat_add(this, B)
!      class(type_bandmat), intent(in) :: this
!      class(type_bandmat), intent(in) :: B
!
!      type(type_bandmat), allocatable :: mat_add
!
!      if ( .not. allocated(mat_add%dat) ) then
!         call mat_add%initialize(shape(this%dat))
!      end if
!
!      mat_add%dat(:,:) = this%dat(:,:)+this%B(:,:)
!
!   end function mat_add
!------------------------------------------------------------------------------
end module band_matrices

module bordered_matrices

   use precision_mod
   use mem_alloc
   use real_matrices, only: type_realmat
   use algebra, only: solve_bordered, prepare_bordered
   use iso_c_binding

   implicit none

   type, public, extends(type_realmat) :: type_bordmat
      real(cp), allocatable :: A1(:,:)
      real(cp), allocatable :: A2(:,:)
      real(cp), allocatable :: A3(:)
      real(cp), allocatable :: A4(:,:)
      integer, allocatable :: pivA1(:)
      integer, allocatable :: pivA4(:)
      integer :: kl
      integer :: ku
      integer :: nfull
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: prepare
      procedure :: solve_real_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
!      procedure :: mat_add
!      generic :: operator(+) => mat_add
   end type type_bordmat

contains

#ifdef WITH_OMP_GPU
   subroutine initialize(this, nx, ny, l_pivot, use_gpu, nfull)
#else
   subroutine initialize(this, nx, ny, l_pivot, nfull)
#endif
      !
      ! Memory allocation
      !
      class(type_bordmat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull
#ifdef WITH_OMP_GPU
      logical, optional, intent(in) :: use_gpu
#endif

      !--
      logical :: loc_use_gpu
      loc_use_gpu = .false.
      this%gpu_is_used=.false.

      this%nrow = nx
      this%ncol = ny
      this%l_pivot = l_pivot

      this%kl = (nx-1)/2
      this%ku = this%kl
      this%nfull = nfull

      allocate( this%A1(nx+(nx-1)/2, ny) )
      allocate( this%A2(ny,nfull) )
      allocate( this%A3(ny) )
      allocate( this%A4(nfull,nfull) )
      this%A1(:,:)=0.0_cp
      this%A2(:,:)=0.0_cp
      this%A3(:)  =0.0_cp
      this%A4(:,:)=0.0_cp
      bytes_allocated = bytes_allocated+(nfull*nfull+nfull*ny+ny+ &
      &                 (nx+(nx-1)/2)*ny)*SIZEOF_DEF_REAL

      if ( this%l_pivot ) then
         allocate( this%pivA1(ny) )
         allocate( this%pivA4(nfull) )
         bytes_allocated = bytes_allocated+(ny+nfull)*SIZEOF_INTEGER
         this%pivA1(:) = 0
         this%pivA4(:) = 0
      end if

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(type_bordmat) :: this

      deallocate( this%A1, this%A2, this%A3, this%A4 )
      if ( this%l_pivot  ) deallocate( this%pivA1, this%pivA4 )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine prepare(this, info, handle, devInfo)

      class(type_bordmat) :: this
      integer, intent(out) :: info
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      call prepare_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &                this%kl,this%ku,this%pivA1,this%pivA4,info)

   end subroutine prepare
!------------------------------------------------------------------------------
   subroutine solve_real_single(this, rhs, handle, devInfo)

      class(type_bordmat) :: this
      real(cp), intent(inout) :: rhs(:)
      type(c_ptr), optional,    intent(inout) :: handle
      integer, optional,        intent(inout) :: devInfo(:)

      !-- Local variable :
      integer :: lenRhs

      lenRhs = this%nfull+this%ncol
      call solve_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &              this%kl,this%ku,this%pivA1,this%pivA4,rhs,lenRhs)

   end subroutine solve_real_single
!------------------------------------------------------------------------------
   subroutine solve_complex_single(this, rhs, tmpr, tmpi, handle, devInfo)

      class(type_bordmat) :: this
      complex(cp), intent(inout) :: rhs(:)

      !-- Output variables
      real(cp), optional,    intent(inout) :: tmpr(:)
      real(cp), optional,    intent(inout) :: tmpi(:)
      type(c_ptr), optional, intent(inout) :: handle
      integer, optional,     intent(inout) :: devInfo(:)

      !-- Local variable :
      integer :: lenRhs

      lenRhs = this%nfull+this%ncol
      call solve_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &              this%kl,this%ku,this%pivA1,this%pivA4,rhs,lenRhs)

   end subroutine solve_complex_single
!------------------------------------------------------------------------------
   subroutine solve_real_multi(this, rhs, nRHS, handle, devInfo)

      class(type_bordmat) :: this
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)
      integer,        optional, intent(inout) :: devInfo(:)
      type(c_ptr),    optional, intent(inout) :: handle

      call solve_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &              this%kl,this%ku,this%pivA1,this%pivA4,rhs,nRHS)

   end subroutine solve_real_multi
!------------------------------------------------------------------------------
   subroutine set_data(this, dat)

      class(type_bordmat) :: this
      real(cp), intent(in) :: dat(:,:)

      !-- Local variables
      integer :: i, j

      !-- A1 = band matrix
      do j=1,this%ncol
         do i=max(1,j-this%ku),min(this%ncol,j+this%kl)
            this%A1(this%kl+this%ku+1+i-j,j)=dat(i,j)
         end do
      end do

      do j=1,this%nfull
         do i=1,this%ncol
            this%A2(i,j)=dat(i,this%ncol+j)
         end do
      end do

      do j=1,this%ncol
         this%A3(j)=dat(this%ncol+1,j)
      end do

      do j=1,this%nfull
         do i=1,this%nfull
            this%A4(i,j)=dat(this%ncol+i,this%ncol+j)
         end do
      end do

   end subroutine set_data
!------------------------------------------------------------------------------
end module bordered_matrices
