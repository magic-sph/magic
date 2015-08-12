!$Id$
module output_data
   !-------------------------------------------------------
   !  Parameters for output control
   !-------------------------------------------------------

   use precision_mod, only: cp
   use logic, only: l_mag, l_anel, l_perpPar, l_r_field, l_r_fieldT, &
                    l_RMS, l_RMStest, l_save_out, l_cmb_field,       &
                    l_rot_ic, l_rot_ma, l_power, l_SRIC, l_SRMA,     &
                    l_dt_cmb_field, l_AM, l_movie
   use parallel_mod, only: rank
   use charmanip, only: length_to_blank

   implicit none

   private
 
   !----- Identification of run:
   character(len=64), public :: runid
 
   !----- Information for graphic output grid:
   integer, public :: ngform
   logical, public :: l_graph_time
 
   !----- Output time control:
   real(cp), public :: t_graph_start,t_graph_stop,dt_graph
   real(cp), public :: t_rst_start,t_rst_stop,dt_rst
   real(cp), public :: t_log_start,t_log_stop,dt_log
   real(cp), public :: t_p_start,t_p_stop,dt_p
   real(cp), public :: t_spec_start,t_spec_stop,dt_spec
   real(cp), public :: t_diag_start,t_diag_stop,dt_diag
   real(cp), public :: t_cmb_start,t_cmb_stop,dt_cmb
   real(cp), public :: t_r_field_start,t_r_field_stop,dt_r_field
   real(cp), public :: t_TO_start,t_TO_stop,dt_TO
   real(cp), public :: t_TOZ_start,t_TOZ_stop,dt_TOZ
   real(cp), public :: t_movie_start,t_movie_stop,dt_movie
   real(cp), public :: t_TOmovie_start,t_TOmovie_stop,dt_TOmovie
   real(cp), public :: t_Bpot_start,t_Bpot_stop,dt_Bpot
   real(cp), public :: t_Vpot_start,t_Vpot_stop,dt_Vpot
   real(cp), public :: t_Tpot_start,t_Tpot_stop,dt_Tpot
   real(cp), public :: t_pot_start,t_pot_stop,dt_pot
   integer, public :: n_graph_step,n_graphs,n_t_graph
   integer, public :: n_rst_step,n_rsts,n_t_rst,n_stores
   integer, public :: n_log_step,n_logs,n_t_log
   integer, public :: n_p_step,n_ps,n_t_p
   integer, public :: n_spec_step,n_specs,n_t_spec
   integer, public :: n_diag_step,n_diags,n_t_diag
   integer, public :: n_cmb_step,n_cmbs,n_t_cmb
   integer, public :: n_r_field_step,n_r_fields,n_t_r_field
   integer, public :: n_movie_step,n_movie_frames,n_t_movie
   integer, public :: n_TO_step,n_TOs,n_t_TO
   integer, public :: n_TOZ_step,n_TOZs,n_t_TOZ
   integer, public :: n_TOmovie_step,n_TOmovie_frames,n_t_TOmovie
   integer, public :: n_Bpot_step,n_Bpots,n_t_Bpot      
   integer, public :: n_Vpot_step,n_Vpots,n_t_Vpot      
   integer, public :: n_Tpot_step,n_Tpots,n_t_Tpot      
   integer, public :: n_pot_step,n_pots,n_t_pot      
   integer, public, parameter :: n_time_hits=5000
   real(cp), public ::  t_graph(n_time_hits)
   real(cp), public ::  t_rst(n_time_hits)
   real(cp), public ::  t_log(n_time_hits)
   real(cp), public ::  t_p(n_time_hits)
   real(cp), public ::  t_spec(n_time_hits)
   real(cp), public ::  t_diag(n_time_hits)
   real(cp), public ::  t_cmb(n_time_hits)
   real(cp), public ::  t_r_field(n_time_hits)
   real(cp), public ::  t_movie(n_time_hits)
   real(cp), public ::  t_TO(n_time_hits)
   real(cp), public ::  t_TOZ(n_time_hits)
   real(cp), public ::  t_TOmovie(n_time_hits)
   real(cp), public ::  t_Bpot(n_time_hits)
   real(cp), public ::  t_Vpot(n_time_hits)
   real(cp), public ::  t_Tpot(n_time_hits)
   real(cp), public ::  t_pot(n_time_hits)
 
   !----- Output radii and degrees for coeff files:
   integer, public :: n_coeff_r_max
   integer, public, allocatable :: n_coeff_r(:)
   integer, public :: n_r_array(100)
   integer, public :: l_max_cmb
   integer, public :: l_max_r
   integer, public :: n_r_step
 
   !----- Output files:
   integer, public :: n_log_file,nLF
#ifdef WITH_MPI
   integer, public :: graph_mpi_fh
   integer, public :: rst_mpi_fh
#endif
   integer, public :: n_graph_file
   integer, public :: n_lp_file,n_rst_file
   integer, public :: n_e_mag_oc_file,n_e_mag_ic_file,n_e_kin_file
   integer, public :: n_u_square_file,n_par_file,n_angular_file
   integer, public :: n_dtvrms_file,n_dtvasrms_file
   integer, public :: n_dtbrms_file,n_dtdrms_file
   integer, public :: n_mag_spec_file,n_kin_spec_file,n_u2_spec_file
   integer, public :: n_rot_file
   integer, public :: n_perpPar_file
   integer, public :: n_dipole_file
   integer, public :: n_cmb_file,n_dt_cmb_file
   integer, public :: n_cmbMov_file
   integer, public :: n_misc_file
   integer, public :: n_SRIC_file
   integer, public :: n_SRMA_file
   integer, public, allocatable :: n_v_r_file(:)
   integer, public, allocatable :: n_t_r_file(:)
   integer, public, allocatable:: n_b_r_file(:)
   integer, public :: n_power_file
   integer, public :: n_signal_file
 
   character(len=55), public :: tag
#ifdef WITH_MPI
   character(len=55), public :: tag_wo_rank
#endif
   character(len=72), public :: log_file
   character(len=72), public :: graph_file
   character(len=72), public :: lp_file
   character(len=72), public :: rst_file
   character(len=72), public :: e_mag_oc_file,e_mag_ic_file
   character(len=72), public :: e_kin_file
   character(len=72), public :: u_square_file
   character(len=72), public :: perpPar_file
   character(len=72), public :: par_file
   character(len=72), public :: angular_file
   character(len=72), public :: dtvrms_file,dtvasrms_file
   character(len=72), public :: dtbrms_file,dtdrms_file
   character(len=72), public :: rot_file
   character(len=72), public :: dipole_file
   character(len=72), public :: cmb_file,dt_cmb_file
   character(len=72), public :: cmbMov_file
   character(len=72), public :: misc_file
   character(len=72), public :: SRIC_file
   character(len=72), public :: SRMA_file
   character(len=72), public, allocatable :: v_r_file(:)
   character(len=72), public, allocatable :: t_r_file(:)
   character(len=72), public, allocatable :: b_r_file(:)
   character(len=72), public :: power_file
   !----- Z-integrated output:
   real(cp), public :: zDens,sDens
 
   !----- RMS cut radius and dealiasing:
   real(cp), public :: rCut,rDea

   public :: openFiles, closeFiles

contains

   subroutine openFiles
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  | Defines names and unit for output files and opens then.           |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      character(len=72) :: string
      integer ::  length
      integer :: n

      !-- Define output unit numbers:
      n_log_file         =7
      n_rst_file         =9
      n_graph_file       =10
      n_e_kin_file       =11
      n_kin_spec_file    =12
      n_e_mag_ic_file    =13
      n_e_mag_oc_file    =14
      n_mag_spec_file    =15
      n_u2_spec_file     =121
      n_lp_file          =16
      n_rot_file         =17
      n_dipole_file      =18
      n_signal_file      =19
      n_cmb_file         =20
      n_misc_file        =21
      n_cmbMov_file      =24
      n_SRIC_file        =25
      n_SRMA_file        =26
      n_dt_cmb_file      =27
      n_power_file       =28
      n_u_square_file    =29  
      n_perpPar_file     =290
      n_par_file         =300
      n_angular_file     =301
      n_dtvrms_file      =302
      n_dtvasrms_file    =303
      n_dtbrms_file      =304
      n_dtdrms_file      =305
      do n=1,n_coeff_r_max
         n_v_r_file(n)=40+n
      end do
      if ( l_r_fieldT ) then
         do n=1,n_coeff_r_max
            n_t_r_file(n)=60+n
         end do
      end if
      if ( l_mag ) then
         do n=1,n_coeff_r_max
            n_b_r_file(n)=50+n
         end do
      end if
      nLF=n_log_file

      !-- Define output file names:
      e_kin_file='e_kin.'//tag
      log_file='log.'//tag
      par_file='par.'//tag
      if ( l_mag ) then
         e_mag_ic_file='e_mag_ic.'//tag
         e_mag_oc_file='e_mag_oc.'//tag
         dipole_file='dipole.'//tag
         if ( l_RMS .or. l_RMStest) then
            dtbrms_file='dtBrms.'//tag
            dtdrms_file='dtDrms.'//tag
         end if
      end if
      if ( l_AM ) then
         angular_file='AM.'//tag
      end if
      if ( l_anel ) then
         u_square_file='u_square.'//tag
      end if
      if ( l_perpPar ) then
         perpPar_file='perpPar.'//tag
      end if
      if ( l_RMS .or. l_RMStest) then
         dtvrms_file='dtVrms.'//tag
         dtvasrms_file='dtVAsRms.'//tag
      end if
      if ( l_rot_ic .or. l_rot_ma ) then
         rot_file='rot.'//tag
      end if
      if ( l_cmb_field ) then
         cmb_file   ='B_coeff_cmb.'//tag
         cmbMov_file='B_coeff_cmbMov.'//tag
      end if
      if ( l_dt_cmb_field ) then
         dt_cmb_file   ='B_coeff_dt_cmb.'//tag
      end if
      if ( l_r_field ) then
         do n=1,n_coeff_r_max
            write(string,'(''V_coeff_r'',i1,''.'')') n
            length=length_to_blank(string)
            v_r_file(n)=string(1:length)//tag
            if ( l_mag ) then
               write(string,'(''B_coeff_r'',i1,''.'')') n
               length=length_to_blank(string)
               B_r_file(n)=string(1:length)//tag
            end if
         end do
      end if
      if ( l_r_fieldT ) then
         do n=1,n_coeff_r_max
            write(string,'(''T_coeff_r'',i1,''.'')') n
            length=length_to_blank(string)
            t_r_file(n)=string(1:length)//tag
         end do
      end if
      misc_file ='misc.'//tag
      SRIC_file ='SRIC.'//tag
      SRMA_file ='SRMA.'//tag
      power_file='power.'//tag

      !-- Open various output files that will be used throughout the run:
      if ( .not. l_save_out ) then

         if ( rank == 0 ) then
            open(n_log_file, file=log_file, status='unknown')
            open(n_e_kin_file, file=e_kin_file, status='new')
            open(n_par_file, file=par_file, status='new')
            if ( l_RMS .or. l_RMStest) then
               open(n_dtvrms_file, file=dtvrms_file, status='new')
               open(n_dtvasrms_file, file=dtvasrms_file, status='new')
            end if
            if ( l_anel ) then
               open(n_u_square_file, file=u_square_file, status='new')
            end if
            if ( l_perpPar ) then
               open(n_perpPar_file, file=perpPar_file, status='new')
            end if
            if ( l_AM ) then
               open(n_angular_file, file=angular_file, status='new')
            end if
            if ( l_r_field ) then
               do n=1,n_coeff_r_max
                  open(n_v_r_file(n), file=v_r_file(n), &
                       status='new', form='unformatted')
               end do
            end if
            if ( l_r_fieldT ) then
               do n=1,n_coeff_r_max
                  open(n_t_r_file(n), file=t_r_file(n), &
                       status='new', form='unformatted')
               end do
            end if
            if ( l_mag ) then
               open(n_e_mag_oc_file, file=e_mag_oc_file, status='new')
               open(n_e_mag_ic_file, file=e_mag_ic_file, status='new')
               open(n_dipole_file, file=dipole_file, status='new')
               if ( l_RMS .or. l_RMStest) then
                  open(n_dtbrms_file, file=dtbrms_file, status='new')
                  open(n_dtdrms_file, file=dtdrms_file, status='new')
               end if
               if ( l_cmb_field ) then
                  open(n_cmb_file, file=cmb_file, &
                       status='new', form='unformatted')
                  if ( l_movie ) then
                     open(n_cmbMov_file, file=cmbMov_file, &
                          status='new', form='unformatted')
                  end if
               end if
               if ( l_dt_cmb_field )                    &
                  open(n_dt_cmb_file, file=dt_cmb_file, &
                       status='new', form='unformatted')
               if ( l_r_field ) then
                  do n=1,n_coeff_r_max
                     open(n_b_r_file(n), file=b_r_file(n), &
                          status='new', form='unformatted')
                  end do
               end if
            end if
            if ( .not. l_SRIC .and. .not. l_SRMA ) then
               if ( l_rot_ic .or. l_rot_ma ) &
                    open(n_rot_file, file=rot_file, status="NEW")
            end if
            open(n_misc_file, file=misc_file, status="NEW")
            if ( l_power ) open(n_power_file, file=power_file, status="NEW")
         end if
      end if

   end subroutine openFiles
!-------------------------------------------------------------------
   subroutine closeFiles
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  | Defines names and unit for output files and opens then.           |
      !  | MPI: called only by the processor responsible for output !        |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      integer :: n
    
      !-- Close various output files:
      if ( .not. l_save_out ) then
    
         close(n_log_file)
    
         if ( rank == 0 ) then
            close(n_e_kin_file)
            close(n_par_file)
            if ( l_AM ) then
               close(n_angular_file)
            end if
            if ( l_anel ) then
               close(n_u_square_file)
            end if
            if ( l_perpPar ) then
               close(n_perpPar_file)
            end if
            if ( l_RMS .or. l_RMStest ) then
               close(n_dtvrms_file)
               close(n_dtvasrms_file)
            end if
            if ( l_r_field ) then
               do n=1,n_coeff_r_max
                  close(n_v_r_file(n))
               end do
            end if
            if ( l_r_fieldT ) then
               do n=1,n_coeff_r_max
                  close(n_t_r_file(n))
               end do
            end if
            if ( l_mag ) then
               close(n_e_mag_oc_file)
               close(n_e_mag_ic_file)
               close(n_dipole_file)
               if ( l_RMS .or. l_RMStest) then
                  close(n_dtbrms_file)
                  close(n_dtdrms_file)
               end if
               if ( l_cmb_field ) close(n_cmb_file)
               if ( l_cmb_field .and. l_movie ) close(n_cmbMov_file)
               if ( l_r_field ) then
                  do n=1,n_coeff_r_max
                     close(n_b_r_file(n))
                  end do
               end if
               if ( l_dt_cmb_field ) close(n_dt_cmb_file)
            end if
            if ( l_rot_ic .or. l_rot_ma .and.      &
                 .not. l_SRIC .and. .not. l_SRMA ) &
                 close(n_rot_file)
            close(n_misc_file)
            close(n_power_file)
            !end if
         end if
      end if

   end subroutine closeFiles
!-------------------------------------------------------------------
end module output_data
