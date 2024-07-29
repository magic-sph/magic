module output_data
   !
   !  This module contains the parameters for output control
   !

   use constants, only: one
   use mem_alloc, only: bytes_allocated
   use precision_mod

   implicit none

   private
 
   !----- Identification of run:
   character(len=64), public :: runid
 
   !----- Output time control:
   real(cp), public :: t_graph_start,t_graph_stop,dt_graph
   real(cp), public :: t_rst_start,t_rst_stop,dt_rst
   real(cp), public :: t_log_start,t_log_stop,dt_log
   real(cp), public :: t_spec_start,t_spec_stop,dt_spec
   real(cp), public :: t_cmb_start,t_cmb_stop,dt_cmb
   real(cp), public :: t_r_field_start,t_r_field_stop,dt_r_field
   real(cp), public :: t_TO_start,t_TO_stop,dt_TO
   real(cp), public :: t_movie_start,t_movie_stop,dt_movie
   real(cp), public :: t_TOmovie_start,t_TOmovie_stop,dt_TOmovie
   real(cp), public :: t_pot_start,t_pot_stop,dt_pot
   real(cp), public :: t_probe_start,t_probe_stop,dt_probe
   integer, public :: n_graph_step,n_graphs
   integer, public :: n_rst_step,n_rsts,n_stores
   integer, public :: n_log_step,n_logs
   integer, public :: n_spec_step,n_specs
   integer, public :: n_cmb_step,n_cmbs
   integer, public :: n_r_field_step,n_r_fields
   integer, public :: n_movie_step,n_movie_frames
   integer, public :: n_TO_step,n_TOs
   integer, public :: n_TOmovie_step,n_TOmovie_frames
   integer, public :: n_pot_step,n_pots
   integer, public :: n_probe_step,n_probe_out
   integer, public, parameter :: n_time_hits=30 ! Maximum number of specific times for I/O in input namelist
   real(cp), allocatable, public ::  t_graph(:)
   real(cp), allocatable, public ::  t_rst(:)
   real(cp), allocatable, public ::  t_log(:)
   real(cp), allocatable, public ::  t_spec(:)
   real(cp), allocatable, public ::  t_cmb(:)
   real(cp), allocatable, public ::  t_r_field(:)
   real(cp), allocatable, public ::  t_movie(:)
   real(cp), allocatable, public ::  t_TO(:)
   real(cp), allocatable, public ::  t_TOmovie(:)
   real(cp), allocatable, public ::  t_pot(:)
   real(cp), allocatable, public ::  t_probe(:)
 
   !----- Output radii and degrees for coeff files:
   integer, public :: n_coeff_r_max
   integer, public, allocatable :: n_coeff_r(:)
   integer, public :: n_r_array(100)
   integer, public :: l_max_cmb
   integer, public :: l_max_comp ! Maximum spherical harmonic degree to estimate Earth-likeness
   integer, public :: l_geo ! max degree for geomagnetic field seen on Earth
   integer, public :: l_max_r
   integer, public :: n_r_step
   integer, public :: m_max_modes
 
   !----- Output files:
   integer, public :: n_log_file
 
   character(len=55), public :: tag
   character(len=72), public :: log_file
   character(len=72), public :: lp_file
   !----- Z-integrated output:
   real(cp), public :: sDens ! Density in s when using z-integration
   real(cp), public :: zDens ! Density in z when using z-integration

   !----- RMS cut radius and dealiasing:
   real(cp), public :: rCut, rDea

   public :: initialize_output_data, finalize_output_data

contains

   subroutine initialize_output_data()
      !
      ! This subroutine allocates the arrays used in the input namelist
      ! to store the times for diagnostics
      !

      allocate(t_graph(n_time_hits), t_rst(n_time_hits), t_log(n_time_hits))
      allocate(t_spec(n_time_hits), t_cmb(n_time_hits), t_r_field(n_time_hits))
      allocate(t_movie(n_time_hits), t_pot(n_time_hits), t_TO(n_time_hits))
      allocate(t_TOmovie(n_time_hits), t_probe(n_time_hits))
      bytes_allocated=bytes_allocated+11*n_time_hits*SIZEOF_DEF_REAL

      !-- Fill with negative values
      t_graph(:)  =-one; t_rst(:)    =-one; t_log(:)    =-one
      t_spec(:)   =-one; t_cmb(:)    =-one; t_r_field(:)=-one
      t_movie(:)  =-one; t_pot(:)    =-one; t_TO(:)     =-one
      t_TOmovie(:)=-one; t_probe(:)  =-one

   end subroutine initialize_output_data
!-----------------------------------------------------------------------------------
   subroutine finalize_output_data()
      !
      ! This subroutine deallocates the arrays used in the input namelist
      ! to store the times for diagnostics
      !

      deallocate(t_graph, t_rst, t_log, t_spec, t_cmb, t_r_field)
      deallocate(t_movie, t_pot, t_TO, t_TOmovie, t_probe)

   end subroutine finalize_output_data
!-----------------------------------------------------------------------------------
end module output_data
