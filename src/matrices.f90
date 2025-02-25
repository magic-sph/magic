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
      procedure(solve_real_multi_if), deferred :: solve_real_multi
      procedure(solve_real_single_if), deferred :: solve_real_single
      procedure(solve_complex_single_if), deferred :: solve_complex_single
      generic :: solve => solve_real_single, solve_complex_single, solve_real_multi
      procedure(set_data_if), deferred :: set_data
   end type type_realmat

   interface

      subroutine initialize_if(this, nx, ny, l_pivot, nfull)
         import
         class(type_realmat) :: this
         integer, intent(in) :: nx
         integer, intent(in) :: ny
         logical, intent(in) :: l_pivot
         integer, optional, intent(in) :: nfull
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

      subroutine solve_real_multi_if(this, rhs, nRHS)
         import
         class(type_realmat) :: this
         integer,  intent(in) :: nRHS
         real(cp), intent(inout) :: rhs(:,:)
      end subroutine solve_real_multi_if

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
      procedure :: solve_real_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
   end type type_densemat

contains

   subroutine initialize(this, nx, ny, l_pivot, nfull)
      !
      ! Memory allocation
      !
      class(type_densemat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull

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
   subroutine solve_real_multi(this, rhs, nRHS)

      class(type_densemat) :: this
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)

      call solve_mat(this%dat, this%nrow, this%nrow, this%pivot, rhs, nRHS)

   end subroutine solve_real_multi
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
      procedure :: solve_real_multi
      procedure :: solve_real_single
      procedure :: solve_complex_single
      procedure :: set_data
!      procedure :: mat_add
!      generic :: operator(+) => mat_add
   end type type_bandmat

contains

   subroutine initialize(this, nx, ny, l_pivot, nfull)
      !
      ! Memory allocation
      !
      class(type_bandmat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull

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
   subroutine solve_real_multi(this, rhs, nRHS)

      class(type_bandmat) :: this
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)

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

module bordered_matrices

   use precision_mod
   use mem_alloc
   use real_matrices, only: type_realmat
   use algebra, only: solve_bordered, prepare_bordered

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

   subroutine initialize(this, nx, ny, l_pivot, nfull)
      !
      ! Memory allocation
      !
      class(type_bordmat) :: this
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      logical, intent(in) :: l_pivot
      integer, optional, intent(in) :: nfull

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
   subroutine prepare(this, info)

      class(type_bordmat) :: this
      integer, intent(out) :: info

      call prepare_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &                this%kl,this%ku,this%pivA1,this%pivA4,info)

   end subroutine prepare
!------------------------------------------------------------------------------
   subroutine solve_real_single(this, rhs)

      class(type_bordmat) :: this
      real(cp), intent(inout) :: rhs(:)

      !-- Local variable :
      integer :: lenRhs

      lenRhs = this%nfull+this%ncol
      call solve_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &              this%kl,this%ku,this%pivA1,this%pivA4,rhs,lenRhs)

   end subroutine solve_real_single
!------------------------------------------------------------------------------
   subroutine solve_complex_single(this, rhs)

      class(type_bordmat) :: this
      complex(cp), intent(inout) :: rhs(:)

      !-- Local variable :
      integer :: lenRhs

      lenRhs = this%nfull+this%ncol
      call solve_bordered(this%A1,this%A2,this%A3,this%A4,this%ncol,this%nfull, &
           &              this%kl,this%ku,this%pivA1,this%pivA4,rhs,lenRhs)

   end subroutine solve_complex_single
!------------------------------------------------------------------------------
   subroutine solve_real_multi(this, rhs, nRHS)

      class(type_bordmat) :: this
      integer,  intent(in) :: nRHS
      real(cp), intent(inout) :: rhs(:,:)

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
