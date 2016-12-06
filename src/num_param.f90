module num_param
   !
   !  Module containing numerical and control parameters
   !

   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max
   use precision_mod

   implicit none

   private
 
   !-- Time step control:
   integer, public :: n_time_steps     ! Total number of time steps requested in the name list
   real(cp), public :: alpha           ! Weight for implicit time step
   real(cp), public :: dtstart         ! Initial time step if start solution is initialized
   real(cp), public :: dtMin           ! Minimum allowed time step
   real(cp), public :: dtMax           ! Maximum allowed time step
   real(cp), public :: timeStart       ! Numerical time where run should start
   real(cp), public :: tEND            ! Numerical time where run should end
 
   !-- Z-angular momentum at start of integration:
   real(cp), public :: AMstart
 
   !-- Courant criteria:
   integer, public :: n_cour_step      ! Step for controlling  Courant criteria
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
   character(len=72), public :: thermo_variable ! thermodynamic variable: S or T

 
   !-- Stop signal:
   integer, public :: istop            ! Variable used in FFT soubroutine
 
   !-- Controlling run time:
   integer, public :: runTimeLimit(4)  ! Maximum running time
   integer, public :: runTime(4)       ! Running time
   integer, public :: runTimeStart(4)  ! Wall clock time of start of the run

   public :: initialize_num_param      ! Subroutine that allocates auxiliary arrays delxr2 and delxh2
   public :: finalize_num_param        ! Subroutine that deallocates arrays

contains

   subroutine initialize_num_param

      allocate( delxr2(n_r_max),delxh2(n_r_max) )
      bytes_allocated = bytes_allocated+2*n_r_max*SIZEOF_DEF_REAL

   end subroutine initialize_num_param
!-------------------------------------------------------------------------------
   subroutine finalize_num_param

      deallocate( delxr2, delxh2 )

   end subroutine finalize_num_param
!-------------------------------------------------------------------------------
end module num_param
