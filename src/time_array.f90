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

   subroutine initialize(this, nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
              &          norder_imp_lin)

      class(type_tarray) :: this

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      allocate( this%impl(nMstart:nMstop,n_r_max,norder_imp_lin-1) )
      !allocate( this%expl(nMstart:nMstop,n_r_max,norder_exp) )
      allocate( this%old(nMstart:nMstop,n_r_max,norder_imp-1) )

      bytes_allocated = bytes_allocated + (nMstop-nMstart+1)*n_r_max*(&
      &                 norder_imp+norder_imp_lin-2)*SIZEOF_DEF_COMPLEX

      this%old(:,:,:) =zero
      !this%expl(:,:,:)=zero
      this%impl(:,:,:)=zero

   end subroutine initialize
!----------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tarray) :: this

      deallocate( this%old, this%impl )

   end subroutine finalize
!----------------------------------------------------------------------------------
   subroutine initialize_scalar(this, norder_imp, norder_exp, norder_imp_lin )

      class(type_tscalar) :: this

      !-- Input variables
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      allocate( this%expl(norder_exp) )
      allocate( this%old(norder_imp-1) )
      allocate( this%impl(norder_imp_lin-1) )

      bytes_allocated = bytes_allocated + (norder_exp+norder_imp+norder_imp_lin-2)*&
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
