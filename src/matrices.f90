module real_matrices

   use precision_mod

   implicit none

   private

   type, abstract, public :: type_realmat
      integer :: nrow ! Number of rows
      integer :: ncol ! Number of columns or number of bands
      logical :: l_pivot
      real(cp), allocatable :: dat(:,:) ! Actual data
      integer, allocatable :: pivot(:)
   contains
      procedure(initialize_if), deferred :: initialize
      procedure(finalize_if), deferred :: finalize
      procedure :: mat_add
      generic :: operator(+) => mat_add
      procedure(prepare_if), deferred :: prepare
      procedure(solve_complex_multi_if), deferred :: solve_complex_multi
      procedure(solve_real_single_if), deferred :: solve_real_single
      procedure(solve_complex_single_if), deferred :: solve_complex_single
      generic :: solve => solve_complex_multi, solve_real_single, &
      &                   solve_complex_single
      procedure(set_data_if), deferred :: set_data
   end type type_realmat

   interface

      subroutine initialize_if(this, nx, ny, l_pivot)
         import
         class(type_realmat) :: this
         integer, intent(in) :: nx
         integer, intent(in) :: ny
         logical, intent(in) :: l_pivot
      end subroutine initialize_if

      subroutine finalize_if(this)
         import
         class(type_realmat) :: this
      end subroutine finalize_if

      subroutine prepare_if(this, info)
         import
         class(type_realmat) :: this
         integer, intent(out) :: info
      end subroutine prepare_if

      subroutine solve_complex_multi_if(this, rhs, nRHS)
         import
         class(type_realmat) :: this
         integer,     intent(in) :: nRHS
         complex(cp), intent(inout) :: rhs(:,:)
      end subroutine solve_complex_multi_if

      subroutine solve_real_single_if(this, rhs)
         import
         class(type_realmat) :: this
         real(cp), intent(inout) :: rhs(:)
      end subroutine solve_real_single_if

      subroutine solve_complex_single_if(this, rhs)
         import
         class(type_realmat) :: this
         complex(cp), intent(inout) :: rhs(:)
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

      if ( .not. allocated(mat_add%dat) ) then
         call mat_add%initialize(this%nrow,this%ncol,this%l_pivot)
      end if

      mat_add%dat(:,:) = this%dat(:,:)+B%dat(:,:)

   end function mat_add
!------------------------------------------------------------------------------
end module real_matrices

module dense_matrices

   use precision_mod
   use mem_alloc
   use real_matrices, only: type_realmat
   use algebra, only: solve_mat, prepare_mat

   implicit none

   type, public, extends(type_realmat) :: type_densemat
   contains 
      procedure :: initialize
      procedure :: finalize
      procedure :: prepare
      procedure :: solve_complex_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
   end type type_densemat

contains 

   subroutine initialize(this, nx, ny, l_pivot)
      !
      ! Memory allocation
      !
      class(type_densemat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot

      this%nrow = nx
      this%ncol = ny
      this%l_pivot = l_pivot
      allocate( this%dat(nx, ny) )
      bytes_allocated = bytes_allocated+nx*ny*SIZEOF_DEF_REAL
      if ( this%l_pivot ) then
         allocate( this%pivot(this%nrow) )
         bytes_allocated = bytes_allocated+this%nrow*SIZEOF_INTEGER
      end if
   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(type_densemat) :: this

      deallocate( this%dat )
      if ( this%l_pivot ) deallocate (this%pivot)

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine prepare(this, info)

      class(type_densemat) :: this
      integer, intent(out) :: info

      call prepare_mat(this%dat, this%nrow, this%nrow, this%pivot, info)

   end subroutine prepare
!------------------------------------------------------------------------------
   subroutine solve_real_single(this, rhs)

      class(type_densemat) :: this
      real(cp), intent(inout) :: rhs(:)

      call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs)

   end subroutine solve_real_single
!------------------------------------------------------------------------------
   subroutine solve_complex_single(this, rhs)

      class(type_densemat) :: this
      complex(cp), intent(inout) :: rhs(:)

      call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs)

   end subroutine solve_complex_single
!------------------------------------------------------------------------------
   subroutine solve_complex_multi(this, rhs, nRHS)

      class(type_densemat) :: this
      integer,     intent(in) :: nRHS
      complex(cp), intent(inout) :: rhs(:,:)

      call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs, nRHS)

   end subroutine solve_complex_multi
!------------------------------------------------------------------------------
   subroutine set_data(this, dat)

      class(type_densemat) :: this
      real(cp), intent(in) :: dat(:,:)

      this%dat(:,:) = dat(:,:)

   end subroutine set_data
!------------------------------------------------------------------------------
end module dense_matrices

module band_matrices

   use precision_mod
   use mem_alloc
   use real_matrices, only: type_realmat
   use algebra, only: solve_tridiag, prepare_tridiag, prepare_band, solve_band

   implicit none

   type, public, extends(type_realmat) :: type_bandmat
      real(cp), allocatable :: du2(:)
      integer :: kl
      integer :: ku
   contains 
      procedure :: initialize
      procedure :: finalize
      procedure :: prepare
      procedure :: solve_complex_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
!      procedure :: mat_add
!      generic :: operator(+) => mat_add
   end type type_bandmat

contains 

   subroutine initialize(this, nx, ny, l_pivot)
      !
      ! Memory allocation
      !
      class(type_bandmat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot

      this%nrow = nx
      this%ncol = ny
      this%l_pivot = l_pivot

      this%kl = (nx-1)/2
      this%ku = this%kl

      if ( nx > 3 .and. this%l_pivot ) then
         allocate( this%dat(nx+(nx-1)/2, ny) )
         bytes_allocated = bytes_allocated+(nx+(nx-1)/2)*ny*SIZEOF_DEF_REAL
      else
         allocate( this%dat(nx, ny) )
         bytes_allocated = bytes_allocated+nx*ny*SIZEOF_DEF_REAL
      end if
      if ( this%l_pivot ) then
         allocate( this%pivot(this%ncol) )
         bytes_allocated = bytes_allocated+this%ncol*SIZEOF_INTEGER
         if ( nx == 3 ) then ! Only require for tridiag arrays
            allocate( this%du2(this%ncol-2) ) ! Help array for tridiag
            bytes_allocated = bytes_allocated+(this%ncol-2)*SIZEOF_DEF_REAL
         end if
      end if
   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(type_bandmat) :: this

      deallocate( this%dat )
      if ( this%l_pivot ) then
         deallocate (this%pivot)
         if ( this%nrow == 3 ) deallocate(this%du2)
      end if

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine prepare(this, info)

      class(type_bandmat) :: this
      integer, intent(out) :: info

      if ( this%nrow == 3 ) then
         call prepare_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:),  &
              &               this%dat(1,2:), this%du2, this%ncol,       &
              &               this%pivot, info)
      else
         call prepare_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, info)
      end if

   end subroutine prepare
!------------------------------------------------------------------------------
   subroutine solve_real_single(this, rhs)

      class(type_bandmat) :: this
      real(cp), intent(inout) :: rhs(:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:), &
              &             this%dat(1,2:), this%du2, this%ncol,      &
              &             this%pivot, rhs)
      else
         call solve_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, rhs)
      end if

   end subroutine solve_real_single
!------------------------------------------------------------------------------
   subroutine solve_complex_single(this, rhs)

      class(type_bandmat) :: this
      complex(cp), intent(inout) :: rhs(:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:), &
              &             this%dat(1,2:), this%du2, this%ncol,      &
              &             this%pivot, rhs)
      else
         call solve_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, rhs)
      end if

   end subroutine solve_complex_single
!------------------------------------------------------------------------------
   subroutine solve_complex_multi(this, rhs, nRHS)

      class(type_bandmat) :: this
      integer,     intent(in) :: nRHS
      complex(cp), intent(inout) :: rhs(:,:)

      if ( this%nrow == 3 ) then
         call solve_tridiag(this%dat(3,1:this%ncol-1), this%dat(2,:),   &
              &             this%dat(1,2:), this%du2, this%ncol,        &
              &             this%pivot, rhs, nRHS)
      else 
         call solve_band(this%dat, this%ncol, this%kl, this%ku, this%pivot, &
              &          rhs, nRHS)
      end if

   end subroutine solve_complex_multi
!------------------------------------------------------------------------------
   subroutine set_data(this, dat)

      class(type_bandmat) :: this
      real(cp), intent(in) :: dat(:,:)

      !-- Local variables 
      integer :: i, j

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
