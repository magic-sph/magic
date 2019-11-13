module num_param
   !
   !  Module containing numerical and control parameters
   !

   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max
   use timing, only: timer_type
   use precision_mod

   implicit none

   private

   !-- Time step control:
   integer, public :: n_time_steps     ! Total number of time steps requested in the name list
   real(cp), public :: alpha           ! Weight for implicit time step
   real(cp), public :: dtMin           ! Minimum allowed time step
   real(cp), public :: dtMax           ! Maximum allowed time step
   real(cp), public :: timeStart       ! Numerical time where run should start
   real(cp), public :: tEND            ! Numerical time where run should end

   !-- Z-angular momentum at start of integration:
   real(cp), public :: AMstart

   !-- Courant criteria:
   real(cp), public :: courfac         ! Value to scale velocity in courant criteria
   real(cp), public :: alffac          ! Value to scale Alfen-velocity in courant criteria
   real(cp), public :: intfac          ! Value to re-scale dtMax during simulation
   real(cp), public, allocatable :: delxr2(:) ! Auxiliary arrays containing effective Courant grid intervals
   real(cp), public, allocatable :: delxh2(:) ! Auxiliary arrays containing effective Courant grid intervals

   !-- Hyperdiffusivity:
   integer, public :: ldif             ! Degree where hyperdiffusion starts to act
   integer, public :: ldifexp          ! Exponent for hyperdiffusion function
   real(cp), public :: difeta          ! Amplitude of magnetic hyperdiffusion
   real(cp), public :: difnu           ! Amplitude of viscous hyperdiffusion
   real(cp), public :: difkap          ! Amplitude of thermal hyperdiffusion
   real(cp), public :: difchem         ! Amplitude of chemical hyperdiffusion

   !-- Scalings:
   real(cp), public :: tScale          ! Time scale
   real(cp), public :: lScale          ! Length scale
   real(cp), public :: vScale          ! Velocity scale
   real(cp), public :: pScale
   real(cp), public :: eScale          ! Energy scale
   real(cp), public :: enscale         ! Energies scale
   integer, public :: n_tScale         ! Control time scale
   integer, public :: n_lScale         ! Control length scale

   character(len=72), public :: anelastic_flavour ! version of the anelastic approximation
   character(len=72), public :: polo_flow_eq ! form of the poloidal flow equation: Pressure or Double Curl

   character(len=72), public :: mpi_transp ! Form of the MPI transpose (point to point or alltoall)

   !-- Stop signal:
   integer, public :: istop            ! Variable used in FFT soubroutine

   !-- Controlling run time:
   real(cp), public :: run_time_limit

   !-- For the nonlinear mapping)
   real(cp), public :: alph1  ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alph2  ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
   character(len=72), public :: map_function ! Mapping family: either tangent or arcsin

   type(timer_type), public :: dct_counter ! Time counter for discrete cosine transforms
   type(timer_type), public :: solve_counter ! Time counter for linear solves

   type(timer_type), public :: lm2phy_counter
   type(timer_type), public :: phy2lm_counter
   type(timer_type), public :: nl_counter
   type(timer_type), public :: td_counter
   character(len=72), public :: time_scheme ! Time scheme

   public :: initialize_num_param      ! Subroutine that allocates auxiliary arrays delxr2 and delxh2
   public :: finalize_num_param        ! Subroutine that deallocates arrays

contains

   subroutine initialize_num_param

      allocate( delxr2(n_r_max),delxh2(n_r_max) )
      bytes_allocated = bytes_allocated+2*n_r_max*SIZEOF_DEF_REAL
      call solve_counter%initialize()
      call dct_counter%initialize()
      call lm2phy_counter%initialize()
      call phy2lm_counter%initialize()
      call nl_counter%initialize()
      call td_counter%initialize()

   end subroutine initialize_num_param
!-------------------------------------------------------------------------------
   subroutine finalize_num_param

      deallocate( delxr2, delxh2 )

   end subroutine finalize_num_param
!-------------------------------------------------------------------------------
end module num_param
