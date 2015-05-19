!$Id$
!***************************************************************
!  Common block containing parameters for output control
!****************************************************************
!
!     !------------ This is release 1 level 1  --------------!
!     !------------ Created on 1/17/02  by JW. --------------!
!

MODULE output_data
  USE truncation
  USE movie_data,ONLY: n_movies_max
  IMPLICIT NONE

  !----- Identification of run:
  character(len=64) :: runid

  !----- Information for graphic output grid:
  INTEGER :: ngform
  LOGICAL :: l_graph_time

  !----- Output time control:
  REAL(kind=8) :: t_graph_start,t_graph_stop,dt_graph
  REAL(kind=8) :: t_rst_start,t_rst_stop,dt_rst
  REAL(kind=8) :: t_log_start,t_log_stop,dt_log
  REAL(kind=8) :: t_p_start,t_p_stop,dt_p
  REAL(kind=8) :: t_spec_start,t_spec_stop,dt_spec
  REAL(kind=8) :: t_diag_start,t_diag_stop,dt_diag
  REAL(kind=8) :: t_cmb_start,t_cmb_stop,dt_cmb
  REAL(kind=8) :: t_r_field_start,t_r_field_stop,dt_r_field
  REAL(kind=8) :: t_TO_start,t_TO_stop,dt_TO
  REAL(kind=8) :: t_TOZ_start,t_TOZ_stop,dt_TOZ
  REAL(kind=8) :: t_movie_start,t_movie_stop,dt_movie
  REAL(kind=8) :: t_TOmovie_start,t_TOmovie_stop,dt_TOmovie
  REAL(kind=8) :: t_Bpot_start,t_Bpot_stop,dt_Bpot
  REAL(kind=8) :: t_Vpot_start,t_Vpot_stop,dt_Vpot
  REAL(kind=8) :: t_Tpot_start,t_Tpot_stop,dt_Tpot
  REAL(kind=8) :: t_pot_start,t_pot_stop,dt_pot
  INTEGER :: n_graph_step,n_graphs,n_t_graph
  INTEGER :: n_rst_step,n_rsts,n_t_rst,n_stores
  INTEGER :: n_log_step,n_logs,n_t_log
  INTEGER :: n_p_step,n_ps,n_t_p
  INTEGER :: n_spec_step,n_specs,n_t_spec
  INTEGER :: n_diag_step,n_diags,n_t_diag
  INTEGER :: n_cmb_step,n_cmbs,n_t_cmb
  INTEGER :: n_r_field_step,n_r_fields,n_t_r_field
  INTEGER :: n_movie_step,n_movie_frames,n_t_movie
  INTEGER :: n_TO_step,n_TOs,n_t_TO
  INTEGER :: n_TOZ_step,n_TOZs,n_t_TOZ
  INTEGER :: n_TOmovie_step,n_TOmovie_frames,n_t_TOmovie
  INTEGER :: n_Bpot_step,n_Bpots,n_t_Bpot      
  INTEGER :: n_Vpot_step,n_Vpots,n_t_Vpot      
  INTEGER :: n_Tpot_step,n_Tpots,n_t_Tpot      
  INTEGER :: n_pot_step,n_pots,n_t_pot      
  INTEGER,PARAMETER :: n_time_hits=5000
  REAL(kind=8) ::  t_graph(n_time_hits)
  REAL(kind=8) ::  t_rst(n_time_hits)
  REAL(kind=8) ::  t_log(n_time_hits)
  REAL(kind=8) ::  t_p(n_time_hits)
  REAL(kind=8) ::  t_spec(n_time_hits)
  REAL(kind=8) ::  t_diag(n_time_hits)
  REAL(kind=8) ::  t_cmb(n_time_hits)
  REAL(kind=8) ::  t_r_field(n_time_hits)
  REAL(kind=8) ::  t_movie(n_time_hits)
  REAL(kind=8) ::  t_TO(n_time_hits)
  REAL(kind=8) ::  t_TOZ(n_time_hits)
  REAL(kind=8) ::  t_TOmovie(n_time_hits)
  REAL(kind=8) ::  t_Bpot(n_time_hits)
  REAL(kind=8) ::  t_Vpot(n_time_hits)
  REAL(kind=8) ::  t_Tpot(n_time_hits)
  REAL(kind=8) ::  t_pot(n_time_hits)

  !----- Output radii and degrees for coeff files:
  INTEGER :: n_coeff_r_max
  INTEGER,ALLOCATABLE :: n_coeff_r(:)
  INTEGER :: n_r_array(100)
  INTEGER :: l_max_cmb
  INTEGER :: l_max_r
  INTEGER :: n_r_step

  !----- Output files:
  INTEGER :: n_log_file,nLF
#ifdef WITH_MPI
  integer :: graph_mpi_fh
  integer :: rst_mpi_fh
#endif
  INTEGER :: n_graph_file
  integer :: n_lp_file,n_rst_file
  INTEGER :: n_e_mag_oc_file,n_e_mag_ic_file,n_e_kin_file
  INTEGER :: n_u_square_file,n_par_file,n_angular_file
  INTEGER :: n_dtvrms_file,n_dtvasrms_file
  INTEGER :: n_dtbrms_file,n_dtdrms_file
  INTEGER :: n_mag_spec_file,n_kin_spec_file,n_u2_spec_file
  INTEGER :: n_rot_file
  INTEGER :: n_perpPar_file
  INTEGER :: n_dipole_file
  INTEGER :: n_cmb_file,n_dt_cmb_file
  INTEGER :: n_cmbMov_file
  INTEGER :: n_misc_file
  INTEGER :: n_SRIC_file
  INTEGER :: n_SRMA_file
  INTEGER,ALLOCATABLE :: n_v_r_file(:)
  INTEGER,ALLOCATABLE :: n_t_r_file(:)
  INTEGER,ALLOCATABLE:: n_b_r_file(:)
  INTEGER :: n_power_file
  INTEGER :: n_signal_file

  CHARACTER(len=55) :: tag
#ifdef WITH_MPI
  CHARACTER(len=55) :: tag_wo_rank
#endif
  CHARACTER(len=72) :: log_file
  CHARACTER(len=72) :: graph_file
  CHARACTER(len=72) :: lp_file
  CHARACTER(len=72) :: rst_file
  CHARACTER(len=72) :: e_mag_oc_file,e_mag_ic_file
  CHARACTER(len=72) :: e_kin_file
  CHARACTER(len=72) :: u_square_file
  CHARACTER(len=72) :: perpPar_file
  CHARACTER(len=72) :: par_file
  CHARACTER(len=72) :: angular_file
  CHARACTER(len=72) :: dtvrms_file,dtvasrms_file
  CHARACTER(len=72) :: dtbrms_file,dtdrms_file
  CHARACTER(len=72) :: rot_file
  CHARACTER(len=72) :: dipole_file
  CHARACTER(len=72) :: cmb_file,dt_cmb_file
  CHARACTER(len=72) :: cmbMov_file
  CHARACTER(len=72) :: misc_file
  CHARACTER(len=72) :: movie_file(n_movies_max)
  CHARACTER(len=72) :: SRIC_file
  CHARACTER(len=72) :: SRMA_file
  CHARACTER(len=72),ALLOCATABLE :: v_r_file(:)
  CHARACTER(len=72),ALLOCATABLE :: t_r_file(:)
  CHARACTER(len=72),ALLOCATABLE :: b_r_file(:)
  CHARACTER(len=72) :: power_file
  
  !----- Z-integrated output:
  REAL(kind=8) :: zDens,sDens
  !COMMON/zInts/zDens,sDens

  !----- RMS cut radius and dealiasing:
  REAL(kind=8) :: rCut,rDea
  !COMMON/RMScut/rCut,rDea

END MODULE output_data
