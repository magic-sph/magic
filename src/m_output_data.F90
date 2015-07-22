!$Id$
module output_data
   !***************************************************************
   !  Parameters for output control
   !****************************************************************

   use movie_data, only: n_movies_max

   implicit none
 
   !----- Identification of run:
   character(len=64) :: runid
 
   !----- Information for graphic output grid:
   integer :: ngform
   logical :: l_graph_time
 
   !----- Output time control:
   real(kind=8) :: t_graph_start,t_graph_stop,dt_graph
   real(kind=8) :: t_rst_start,t_rst_stop,dt_rst
   real(kind=8) :: t_log_start,t_log_stop,dt_log
   real(kind=8) :: t_p_start,t_p_stop,dt_p
   real(kind=8) :: t_spec_start,t_spec_stop,dt_spec
   real(kind=8) :: t_diag_start,t_diag_stop,dt_diag
   real(kind=8) :: t_cmb_start,t_cmb_stop,dt_cmb
   real(kind=8) :: t_r_field_start,t_r_field_stop,dt_r_field
   real(kind=8) :: t_TO_start,t_TO_stop,dt_TO
   real(kind=8) :: t_TOZ_start,t_TOZ_stop,dt_TOZ
   real(kind=8) :: t_movie_start,t_movie_stop,dt_movie
   real(kind=8) :: t_TOmovie_start,t_TOmovie_stop,dt_TOmovie
   real(kind=8) :: t_Bpot_start,t_Bpot_stop,dt_Bpot
   real(kind=8) :: t_Vpot_start,t_Vpot_stop,dt_Vpot
   real(kind=8) :: t_Tpot_start,t_Tpot_stop,dt_Tpot
   real(kind=8) :: t_pot_start,t_pot_stop,dt_pot
   integer :: n_graph_step,n_graphs,n_t_graph
   integer :: n_rst_step,n_rsts,n_t_rst,n_stores
   integer :: n_log_step,n_logs,n_t_log
   integer :: n_p_step,n_ps,n_t_p
   integer :: n_spec_step,n_specs,n_t_spec
   integer :: n_diag_step,n_diags,n_t_diag
   integer :: n_cmb_step,n_cmbs,n_t_cmb
   integer :: n_r_field_step,n_r_fields,n_t_r_field
   integer :: n_movie_step,n_movie_frames,n_t_movie
   integer :: n_TO_step,n_TOs,n_t_TO
   integer :: n_TOZ_step,n_TOZs,n_t_TOZ
   integer :: n_TOmovie_step,n_TOmovie_frames,n_t_TOmovie
   integer :: n_Bpot_step,n_Bpots,n_t_Bpot      
   integer :: n_Vpot_step,n_Vpots,n_t_Vpot      
   integer :: n_Tpot_step,n_Tpots,n_t_Tpot      
   integer :: n_pot_step,n_pots,n_t_pot      
   integer, parameter :: n_time_hits=5000
   real(kind=8) ::  t_graph(n_time_hits)
   real(kind=8) ::  t_rst(n_time_hits)
   real(kind=8) ::  t_log(n_time_hits)
   real(kind=8) ::  t_p(n_time_hits)
   real(kind=8) ::  t_spec(n_time_hits)
   real(kind=8) ::  t_diag(n_time_hits)
   real(kind=8) ::  t_cmb(n_time_hits)
   real(kind=8) ::  t_r_field(n_time_hits)
   real(kind=8) ::  t_movie(n_time_hits)
   real(kind=8) ::  t_TO(n_time_hits)
   real(kind=8) ::  t_TOZ(n_time_hits)
   real(kind=8) ::  t_TOmovie(n_time_hits)
   real(kind=8) ::  t_Bpot(n_time_hits)
   real(kind=8) ::  t_Vpot(n_time_hits)
   real(kind=8) ::  t_Tpot(n_time_hits)
   real(kind=8) ::  t_pot(n_time_hits)
 
   !----- Output radii and degrees for coeff files:
   integer :: n_coeff_r_max
   integer, allocatable :: n_coeff_r(:)
   integer :: n_r_array(100)
   integer :: l_max_cmb
   integer :: l_max_r
   integer :: n_r_step
 
   !----- Output files:
   integer :: n_log_file,nLF
#ifdef WITH_MPI
   integer :: graph_mpi_fh
   integer :: rst_mpi_fh
#endif
   integer :: n_graph_file
   integer :: n_lp_file,n_rst_file
   integer :: n_e_mag_oc_file,n_e_mag_ic_file,n_e_kin_file
   integer :: n_u_square_file,n_par_file,n_angular_file
   integer :: n_dtvrms_file,n_dtvasrms_file
   integer :: n_dtbrms_file,n_dtdrms_file
   integer :: n_mag_spec_file,n_kin_spec_file,n_u2_spec_file
   integer :: n_rot_file
   integer :: n_perpPar_file
   integer :: n_dipole_file
   integer :: n_cmb_file,n_dt_cmb_file
   integer :: n_cmbMov_file
   integer :: n_misc_file
   integer :: n_SRIC_file
   integer :: n_SRMA_file
   integer,allocatable :: n_v_r_file(:)
   integer,allocatable :: n_t_r_file(:)
   integer,allocatable:: n_b_r_file(:)
   integer :: n_power_file
   integer :: n_signal_file
 
   character(len=55) :: tag
#ifdef WITH_MPI
   character(len=55) :: tag_wo_rank
#endif
   character(len=72) :: log_file
   character(len=72) :: graph_file
   character(len=72) :: lp_file
   character(len=72) :: rst_file
   character(len=72) :: e_mag_oc_file,e_mag_ic_file
   character(len=72) :: e_kin_file
   character(len=72) :: u_square_file
   character(len=72) :: perpPar_file
   character(len=72) :: par_file
   character(len=72) :: angular_file
   character(len=72) :: dtvrms_file,dtvasrms_file
   character(len=72) :: dtbrms_file,dtdrms_file
   character(len=72) :: rot_file
   character(len=72) :: dipole_file
   character(len=72) :: cmb_file,dt_cmb_file
   character(len=72) :: cmbMov_file
   character(len=72) :: misc_file
   character(len=72) :: movie_file(n_movies_max)
   character(len=72) :: SRIC_file
   character(len=72) :: SRMA_file
   character(len=72), allocatable :: v_r_file(:)
   character(len=72), allocatable :: t_r_file(:)
   character(len=72), allocatable :: b_r_file(:)
   character(len=72) :: power_file
   
   !----- Z-integrated output:
   real(kind=8) :: zDens,sDens
 
   !----- RMS cut radius and dealiasing:
   real(kind=8) :: rCut,rDea

end module output_data
