module time_array

   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_tarray
      complex(cp), allocatable :: impl(:,:,:)
      complex(cp), pointer :: expl(:,:,:)
      complex(cp), allocatable :: old(:,:,:)
      logical :: l_exp
   contains
      procedure :: initialize
      procedure :: finalize
   end type type_tarray

   type, public :: type_tscalar
      real(cp), allocatable :: impl(:)
      real(cp), allocatable :: old(:)
      real(cp), allocatable :: expl(:)
   contains
      procedure :: initialize => initialize_scalar
      procedure :: finalize => finalize_scalar
   end type type_tscalar

contains

   subroutine initialize(this, nMstart, nMstop, n_r_max, nold, norder_exp, &
              &          norder_imp_lin, l_allocate_exp)

      class(type_tarray) :: this

      !-- Input variables
      integer,           intent(in) :: nMstart
      integer,           intent(in) :: nMstop
      integer,           intent(in) :: n_r_max
      integer,           intent(in) :: nold
      integer,           intent(in) :: norder_exp
      integer,           intent(in) :: norder_imp_lin
      logical, optional, intent(in) :: l_allocate_exp

      !-- Local variable
      logical :: l_allocate

      if ( present(l_allocate_exp) ) then
         l_allocate = l_allocate_exp
      else
         l_allocate = .false.
      end if

      this%l_exp = l_allocate

      allocate( this%impl(nMstart:nMstop,n_r_max,norder_imp_lin-1) )
      allocate( this%old(nMstart:nMstop,n_r_max,nold) )

      bytes_allocated = bytes_allocated + (nMstop-nMstart+1)*n_r_max*(&
      &                 nold+norder_imp_lin-1)*SIZEOF_DEF_COMPLEX

      this%old(:,:,:) =zero
      this%impl(:,:,:)=zero

      if ( l_allocate ) then
         allocate( this%expl(nMstart:nMstop,n_r_max,norder_exp) )
         bytes_allocated = bytes_allocated + (nMstop-nMstart+1)*n_r_max*&
         &                 norder_exp*SIZEOF_DEF_COMPLEX
         this%expl(:,:,:)=zero
      end if

   end subroutine initialize
!----------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tarray) :: this

      deallocate( this%old, this%impl )
      if ( this%l_exp ) deallocate( this%expl )

   end subroutine finalize
!----------------------------------------------------------------------------------
   subroutine initialize_scalar(this, nold, norder_exp, norder_imp_lin )

      class(type_tscalar) :: this

      !-- Input variables
      integer, intent(in) :: nold
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      allocate( this%expl(norder_exp) )
      allocate( this%old(nold) )
      allocate( this%impl(norder_imp_lin-1) )

      bytes_allocated = bytes_allocated + (norder_exp+nold+norder_imp_lin-1)*&
      &                 SIZEOF_DEF_REAL

      this%expl(:)=0.0_cp
      this%old(:) =0.0_cp
      this%impl(:)=0.0_cp

   end subroutine initialize_scalar
!----------------------------------------------------------------------------------
   subroutine finalize_scalar(this)

      class(type_tscalar) :: this

      deallocate( this%old, this%expl )

   end subroutine finalize_scalar
!----------------------------------------------------------------------------------
end module time_array
