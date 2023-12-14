module dtB_arrays_mod

   use truncation, only: lm_max
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use precision_mod
   use constants, only: zero

   implicit none

   private

#ifdef WITH_OMP_GPU
   !$omp declare target (dtB_arrays_t)
#endif

   type, public :: dtB_arrays_t
      !----- Local dtB output stuff:
      complex(cp), allocatable :: BtVrLM(:), BpVrLM(:), BrVtLM(:)
      complex(cp), allocatable :: BtVpLM(:), BpVtLM(:), BrVpLM(:)
      complex(cp), allocatable :: BpVtBtVpCotLM(:), BpVtBtVpSn2LM(:)
      complex(cp), allocatable :: BrVZLM(:), BtVZLM(:), BtVZsn2LM(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
   end type dtB_arrays_t

contains

   subroutine initialize(this)

      class(dtB_arrays_t) :: this

      allocate( this%BtVrLM(lm_max), this%BpVrLM(lm_max) )
      allocate( this%BrVtLM(lm_max), this%BrVpLM(lm_max) )
      allocate( this%BtVpLM(lm_max), this%BpVtLM(lm_max) )
      allocate( this%BpVtBtVpCotLM(lm_max), this%BpVtBtVpSn2LM(lm_max) )
      allocate( this%BrVZLM(lm_max), this%BtVZLM(lm_max) )
      allocate( this%BtVZsn2LM(lm_max) )
      bytes_allocated = bytes_allocated+ 11*lm_max*SIZEOF_DEF_COMPLEX

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: this)
      gpu_bytes_allocated = gpu_bytes_allocated+ 11*lm_max_dtB*SIZEOF_DEF_COMPLEX
#endif

      !--
      call this%set_zero()

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(dtB_arrays_t) :: this

#ifdef WITH_OMP_GPU
      !$omp target exit data map(release: this)
#endif

      deallocate( this%BtVrLM, this%BpVrLM, this%BrVtLM )
      deallocate( this%BrVpLM, this%BtVpLM, this%BpVtLM )
      deallocate( this%BpVtBtVpCotLM, this%BpVtBtVpSn2LM )
      deallocate( this%BrVZLM, this%BtVZLM, this%BtVZsn2LM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(dtB_arrays_t) :: this

      this%BtVrLM(:) = zero
      this%BpVrLM(:) = zero
      this%BrVtLM(:) = zero
      this%BrVpLM(:) = zero
      this%BtVpLM(:) = zero
      this%BpVtLM(:) = zero
      this%BrVZLM(:) = zero
      this%BtVZLM(:) = zero
      this%BpVtBtVpCotLM(:) = zero
      this%BpVtBtVpSn2LM(:) = zero
      this%BtVZsn2LM(:) = zero

#ifdef WITH_OMP_GPU
      !$omp target update to(this)
#endif

   end subroutine set_zero
!----------------------------------------------------------------------------
end module dtB_arrays_mod
