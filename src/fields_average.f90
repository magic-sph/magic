module fields_average_mod
   !
   ! This module is used when one wants to store time-averaged quantities
   !

   use truncation
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use communications, only: lo2r_one
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: chebt_ic, chebt_ic_even, r, dr_fac_ic, &
       &                       rscheme_oc, l_R
   use grid_blocking, only: n_phys_space
   use blocking,only: lm2, llm, ulm, llmMag, ulmMag
   use logic, only: l_mag, l_conv, l_save_out, l_heat, l_cond_ic, l_onset, &
       &            l_chemical_conv, l_phase_field, l_batched_shts
   use kinetic_energy, only: get_e_kin
   use magnetic_energy, only: get_e_mag
   use output_data, only: tag, n_log_file, log_file, n_graphs, l_max_cmb
   use parallel_mod, only: rank
   use sht, only: torpol_to_spat, scal_to_spat, sht_l
   use constants, only: zero, vol_oc, vol_ic, one
   use communications, only: get_global_sum, gather_from_lo_to_rank0,&
       &                     gather_all_from_lo_to_rank0,gt_OC,gt_IC
   use out_coeff, only: write_Bcmb
#ifdef WITH_MPI
   use out_coeff, only: write_Pot_mpi
#else
   use out_coeff, only: write_Pot
#endif
   use spectra, only: spectrum
   use graphOut_mod, only: graphOut_IC, open_graph_file, close_graph_file
#ifdef WITH_MPI
   use graphOut_mod, only: graphOut_mpi, graphOut_mpi_header
#else
   use graphOut_mod, only: graphOut, graphOut_header
#endif
   use radial_der_even, only: get_drNS_even, get_ddrNS_even
   use radial_der, only: get_dr
   use fieldsLast, only: dwdt, dpdt, dzdt, dsdt, dxidt, dbdt, djdt, dbdt_ic, &
       &                 djdt_ic, domega_ma_dt, domega_ic_dt, dphidt
   use storeCheckPoints
   use time_schemes, only: type_tscheme

   implicit none

   private

   complex(cp), allocatable :: w_ave_LMloc(:,:), w_ave_Rloc(:,:)
   complex(cp), allocatable :: z_ave_LMloc(:,:), z_ave_Rloc(:,:)
   complex(cp), allocatable :: s_ave_LMloc(:,:), s_ave_Rloc(:,:)
   complex(cp), allocatable :: xi_ave_LMloc(:,:), xi_ave_Rloc(:,:)
   complex(cp), allocatable :: phi_ave_LMloc(:,:), phi_ave_Rloc(:,:)
   complex(cp), allocatable :: p_ave_LMloc(:,:), p_ave_Rloc(:,:)
   complex(cp), allocatable :: b_ave_LMloc(:,:), b_ave_Rloc(:,:)
   complex(cp), allocatable :: aj_ave_LMloc(:,:), aj_ave_Rloc(:,:)
   complex(cp), allocatable :: b_ic_ave(:,:)
   complex(cp), allocatable :: aj_ic_ave(:,:)
   ! on rank 0 we also allocate the following fields
   complex(cp), allocatable :: bICB(:)

   public :: initialize_fields_average_mod, fields_average, &
   &         finalize_fields_average_mod

contains

   subroutine initialize_fields_average_mod

      allocate( w_ave_LMloc(llm:ulm,n_r_max), z_ave_LMloc(llm:ulm,n_r_max) )
      allocate( s_ave_LMloc(llm:ulm,n_r_max), p_ave_LMloc(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+4*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      allocate( w_ave_Rloc(lm_max,nRstart:nRstop), z_ave_Rloc(lm_max,nRstart:nRstop) )
      allocate( s_ave_Rloc(lm_max,nRstart:nRstop), p_ave_Rloc(lm_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      if ( l_mag ) then
         allocate( b_ave_LMloc(llm:ulm,n_r_max), aj_ave_LMloc(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      allocate( b_ave_Rloc(lm_max,nRstart:nRstop), aj_ave_Rloc(lm_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         allocate( b_ic_ave(llm:ulm,n_r_ic_max), aj_ic_ave(llm:ulm,n_r_ic_max) )
         bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_ic_max*SIZEOF_DEF_COMPLEX
      else
         allocate( b_ave_LMloc(1,1), aj_ave_LMloc(1,1) )
         allocate( b_ave_Rloc(1,1), aj_ave_Rloc(1,1) )
         allocate( b_ic_ave(1,1), aj_ic_ave(1,1) )
      end if

      if ( l_chemical_conv ) then
         allocate( xi_ave_LMloc(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         allocate( xi_ave_Rloc(lm_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_ave_LMloc(1,1) )
         allocate( xi_ave_Rloc(1,1) )
      end if

      if ( l_phase_field ) then
         allocate( phi_ave_LMloc(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         allocate( phi_ave_Rloc(lm_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( phi_ave_LMloc(1,1) )
         allocate( phi_ave_Rloc(1,1) )
      end if

      !-- Set initial values to zero
      if ( l_conv ) then
         w_ave_LMloc(:,:)=zero
         z_ave_LMloc(:,:)=zero
         p_ave_LMloc(:,:)=zero
      end if
      if ( l_heat ) s_ave_LMloc(:,:)=zero
      if ( l_chemical_conv ) xi_ave_LMloc(:,:)=zero
      if ( l_phase_field ) phi_ave_LMloc(:,:)=zero
      if ( l_mag ) then
         b_ave_LMloc(:,:) =zero
         aj_ave_LMloc(:,:)=zero
         if ( l_cond_ic ) then
            b_ic_ave(:,:) =zero
            aj_ic_ave(:,:)=zero
         end if
      end if

      if ( rank == 0 ) then
         allocate( bICB(1:lm_max) )
         bytes_allocated = bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX
      else
         allocate( bICB(1) )
      end if

   end subroutine initialize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine finalize_fields_average_mod

      deallocate( w_ave_LMloc, z_ave_LMloc, s_ave_LMloc, p_ave_LMloc )
      deallocate( w_ave_Rloc, z_ave_Rloc, s_ave_Rloc, p_ave_Rloc )
      deallocate( b_ave_LMloc, aj_ave_LMloc, b_ic_ave )
      deallocate( b_ave_Rloc, aj_ave_Rloc, xi_ave_Rloc )
      deallocate( aj_ic_ave, bICB, xi_ave_LMloc, phi_ave_LMloc, phi_ave_Rloc )

   end subroutine finalize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine fields_average(simtime,tscheme,nAve,l_stop_time,        &
      &                      time_passed,time_norm,omega_ic,omega_ma, &
      &                      w,z,p,s,xi,phi,b,aj,b_ic,aj_ic)
      !
      ! This subroutine averages fields b and v over time.
      !

      !-- Input of variables:
      integer,             intent(in) :: nAve         ! number for averaged time steps
      logical,             intent(in) :: l_stop_time  ! true if this is the last time step
      real(cp),            intent(in) :: time_passed  ! time passed since last log
      real(cp),            intent(in) :: time_norm    ! time passed since start of time loop
      real(cp),            intent(in) :: omega_ic,omega_ma
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: simtime
      complex(cp),         intent(in) :: w(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: z(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: p(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: s(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: phi(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),         intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local stuff:
      ! fields for the gathering
      complex(cp) :: b_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: db_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: aj_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)

      !----- Time averaged fields:
      complex(cp) :: dw_ave_LMloc(llm:ulm,n_r_max)
      complex(cp) :: dw_ave_Rloc(lm_max,nRstart:nRstop)
      complex(cp) :: ds_ave_LMloc(llm:ulm,n_r_max)
      complex(cp) :: dxi_ave_LMloc(llm:ulm,n_r_max)
      complex(cp) :: db_ave_LMloc(llm:ulm,n_r_maxMag)
      complex(cp) :: db_ave_Rloc(lm_maxMag,nRstart:nRstop)
      complex(cp) :: db_ic_ave(llm:ulm,n_r_ic_max)
      complex(cp) :: ddb_ic_ave(llm:ulm,n_r_ic_max)
      complex(cp) :: dj_ic_ave(llm:ulm,n_r_ic_max)

      !----- Work array:
      complex(cp) :: workA_LMloc(llm:ulm,n_r_max)

      !----- Fields in grid space:
      real(cp) :: Br(n_phys_space),Bt(n_phys_space)
      real(cp) :: Bp(n_phys_space),Vr(n_phys_space)
      real(cp) :: Vt(n_phys_space),Vp(n_phys_space)
      real(cp) :: Sr(n_phys_space),Prer(n_phys_space)
      real(cp) :: Xir(n_phys_space),Phir(n_phys_space)

      !----- Energies of time average field:
      real(cp) :: e_kin_p_ave,e_kin_t_ave
      real(cp) :: e_kin_p_as_ave,e_kin_t_as_ave
      real(cp) :: e_mag_p_ave,e_mag_t_ave
      real(cp) :: e_mag_p_as_ave,e_mag_t_as_ave
      real(cp) :: e_mag_p_ic_ave,e_mag_t_ic_ave
      real(cp) :: e_mag_p_as_ic_ave,e_mag_t_as_ic_ave
      real(cp) :: e_mag_os_ave,e_mag_as_os_ave
      real(cp) :: Dip,DipCMB,e_cmb,elsAnel

      integer :: nR, n_graph_handle
      integer :: n_e_sets

      character(len=80) :: outFile
      integer :: nOut,n_cmb_sets,nPotSets

      real(cp) :: time
      real(cp) :: dt_norm

      !-- Add new time step:
      if ( l_conv ) then
         w_ave_LMloc(:,:)=w_ave_LMloc(:,:) + time_passed*w(:,:)
         z_ave_LMloc(:,:)=z_ave_LMloc(:,:) + time_passed*z(:,:)
         p_ave_LMloc(:,:)=p_ave_LMloc(:,:) + time_passed*p(:,:)
      end if
      if ( l_phase_field ) phi_ave_LMloc(:,:)=phi_ave_LMloc(:,:) + &
      &                                       time_passed*phi(:,:)
      if ( l_heat ) s_ave_LMloc(:,:)=s_ave_LMloc(:,:) + time_passed*s(:,:)
      if ( l_chemical_conv ) xi_ave_LMloc(:,:)=xi_ave_LMloc(:,:) + &
      &                                        time_passed*xi(:,:)
      if ( l_mag ) then
         b_ave_LMloc(:,:) =b_ave_LMloc(:,:)  + time_passed*b(:,:)
         aj_ave_LMloc(:,:)=aj_ave_LMloc(:,:) + time_passed*aj(:,:)
         if ( l_cond_ic ) then
            b_ic_ave(:,:) =b_ic_ave(:,:) + time_passed*b_ic(:,:)
            aj_ic_ave(:,:)=aj_ic_ave(:,:)+ time_passed*aj_ic(:,:)
         end if
      end if

      !--- Output, intermediate output every 10th averaging to save result
      !    will be overwritten.
      if ( l_stop_time .or. mod(nAve,10) == 0 ) then

         time   =-one  ! This signifies averaging in output files!
         dt_norm=one/time_norm

         if ( l_conv ) then
            w_ave_LMloc(:,:)=dt_norm*w_ave_LMloc(:,:)
            z_ave_LMloc(:,:)=dt_norm*z_ave_LMloc(:,:)
            p_ave_LMloc(:,:)=dt_norm*p_ave_LMloc(:,:)
         end if
         if ( l_phase_field ) phi_ave_LMloc(:,:)=dt_norm*phi_ave_LMloc(:,:)
         if ( l_heat ) s_ave_LMloc(:,:)=dt_norm*s_ave_LMloc(:,:)
         if ( l_chemical_conv ) xi_ave_LMloc(:,:)=dt_norm*xi_ave_LMloc(:,:)
         if ( l_mag ) then
            b_ave_LMloc(:,:) =dt_norm*b_ave_LMloc(:,:)
            aj_ave_LMloc(:,:)=dt_norm*aj_ave_LMloc(:,:)
            if ( l_cond_ic ) then
               b_ic_ave(:,:) =dt_norm*b_ic_ave(:,:)
               aj_ic_ave(:,:)=dt_norm*aj_ic_ave(:,:)
            end if
         end if

         !----- Get the radial derivatives:
         call get_dr(w_ave_LMloc,dw_ave_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max, &
              &      rscheme_oc,nocopy=.true.)
         if ( l_mag ) then
            call get_dr(b_ave_LMloc,db_ave_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max, &
                 &      rscheme_oc,nocopy=.true.)
         end if
         if ( l_heat ) then
            call get_dr(s_ave_LMloc,ds_ave_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max, &
                 &      rscheme_oc,nocopy=.true.)
         end if
         if ( l_chemical_conv ) then
            call get_dr(xi_ave_LMloc,dxi_ave_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max, &
                 &      rscheme_oc,nocopy=.true.)
         end if
         if ( l_cond_ic ) then
            call get_ddrNS_even(b_ic_ave,db_ic_ave,ddb_ic_ave,ulm-llm+1,1,     &
                 &              ulm-llm+1,n_r_ic_max,n_cheb_ic_max,dr_fac_ic,  &
                 &              workA_LMloc,chebt_ic, chebt_ic_even)
            call get_drNS_even(aj_ic_ave,dj_ic_ave,ulm-llm+1,1,ulm-llm+1,      &
                 &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA_LMloc, &
                 &             chebt_ic,chebt_ic_even)
         end if

         !----- Get averaged spectra:
         !      Note: average spectra will be in file no 0
         call spectrum(0,time,.false.,nAve,l_stop_time,time_passed,        &
              &        time_norm,s_ave_LMloc,ds_ave_LMloc,xi_ave_LMloc,    &
              &        dxi_ave_LMloc,phi_ave_LMloc,w_ave_LMloc,            &
              &        dw_ave_LMloc,z_ave_LMloc,b_ave_LMloc,db_ave_LMloc,  &
              &        aj_ave_LMloc,b_ic_ave,db_ic_ave,aj_ic_ave)

         if ( rank==0 .and. l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if

         !----- Write averaged energies into log-file at end of run:
         if ( l_stop_time ) then
            !----- Calculate energies of averaged field:
            n_e_sets=1
            call get_e_kin(time,.false.,.true.,n_e_sets,         &
                 &         w_ave_LMloc,dw_ave_LMloc,z_ave_LMloc, &
                 &         e_kin_p_ave,e_kin_t_ave,              &
                 &         e_kin_p_as_ave,e_kin_t_as_ave)

            call get_e_mag(time,.false.,.true.,n_e_sets,                  &
                 &         b_ave_LMloc,db_ave_LMloc,aj_ave_LMloc,         &
                 &         b_ic_ave,db_ic_ave,aj_ic_ave,                  &
                 &         e_mag_p_ave,e_mag_t_ave,                       &
                 &         e_mag_p_as_ave,e_mag_t_as_ave,                 &
                 &         e_mag_p_ic_ave,e_mag_t_ic_ave,                 &
                 &         e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,           &
                 &         e_mag_os_ave,e_mag_as_os_ave,e_cmb,Dip,DipCMB, &
                 &         elsAnel)

            if ( rank == 0 ) then
               !----- Output of energies of averaged field:
               write(n_log_file,'(/,A)')                           &
               &    ' ! ENERGIES OF TIME AVERAGED FIELD'
               write(n_log_file,                                   &
               &    '('' !  (total,poloidal,toroidal,total density)'')')
               write(n_log_file,'(1P,'' !  Kinetic energies:'',4ES16.6)') &
               &    (e_kin_p_ave+e_kin_t_ave),e_kin_p_ave,e_kin_t_ave,    &
               &    (e_kin_p_ave+e_kin_t_ave)/vol_oc
               write(n_log_file,'(1P,'' !  OC Mag  energies:'',4ES16.6)') &
               &    (e_mag_p_ave+e_mag_t_ave),e_mag_p_ave,e_mag_t_ave,    &
               &    (e_mag_p_ave+e_mag_t_ave)/vol_oc
               write(n_log_file,'(1P,'' !  IC Mag  energies:'',4ES16.6)') &
               &    (e_mag_p_ic_ave+e_mag_t_ic_ave),e_mag_p_ic_ave,       &
               &     e_mag_t_ic_ave,(e_mag_p_ic_ave+e_mag_t_ic_ave)/vol_ic
               write(n_log_file,'(1P,'' !  OS Mag  energies:'',ES16.6)')  &
               &     e_mag_os_ave
               write(n_log_file,'(/,'' !  AXISYMMETRIC PARTS:'')')
               write(n_log_file,                                          &
               &     '('' !  (total,poloidal,toroidal,total density)'')')
               write(n_log_file,'(1P,'' !  Kinetic AS energies:'',4ES16.6)') &
               &    (e_kin_p_as_ave+e_kin_t_as_ave),e_kin_p_as_ave,          &
               &     e_kin_t_as_ave,(e_kin_p_as_ave+e_kin_t_as_ave)/vol_oc
               write(n_log_file,'(1P,'' !  OC Mag  AS energies:'',4ES16.6)') &
               &    (e_mag_p_as_ave+e_mag_t_as_ave),e_mag_p_as_ave,          &
               &     e_mag_t_as_ave,(e_mag_p_as_ave+e_mag_t_as_ave)/vol_oc
               write(n_log_file,'(1P,'' !  IC Mag  AS energies:'',4ES16.6)') &
               &    (e_mag_p_as_ic_ave+e_mag_t_as_ic_ave),e_mag_p_as_ic_ave, &
               &     e_mag_t_as_ic_ave,(e_mag_p_as_ic_ave+e_mag_t_as_ic_ave) &
               &     /vol_ic
               write(n_log_file,'(1P,'' !  OC Mag  AS energies:'',ES16.6)')  &
               &     e_mag_os_ave
               write(n_log_file,'(1P,'' !  Relative ax. dip. E:'',ES16.6)')  &
               &     Dip
            end if
         end if ! End of run ?

         !----- Construct name of graphic file and open it:
         ! For the graphic file of the average fields, we gather them
         ! on rank 0 and use the old serial output routine.
         if ( .not. l_onset ) then
            call open_graph_file(0, time, .true., n_graph_handle)
            !----- Write header into graphic file:
#ifdef WITH_MPI
            call graphOut_mpi_header(time, n_graph_handle)
#else
            call graphOut_header(time, n_graph_handle)
#endif

            !-- This will be needed for the inner core
            if ( l_mag ) then
               call gather_from_lo_to_rank0(b_ave_LMloc(llm,n_r_icb),bICB)
            end if

            !-- MPI transposes from LMloc to Rloc
            call lo2r_one%transp_lm2r(z_ave_LMloc, z_ave_Rloc)
            call lo2r_one%transp_lm2r(w_ave_LMloc, w_ave_Rloc)
            call lo2r_one%transp_lm2r(dw_ave_LMloc, dw_ave_Rloc)
            call lo2r_one%transp_lm2r(s_ave_LMloc, s_ave_Rloc)
            call lo2r_one%transp_lm2r(p_ave_LMloc, p_ave_Rloc)

            if ( l_mag ) then
               call lo2r_one%transp_lm2r(aj_ave_LMloc, aj_ave_Rloc)
               call lo2r_one%transp_lm2r(aj_ave_LMloc, b_ave_Rloc)
               call lo2r_one%transp_lm2r(db_ave_LMloc, db_ave_Rloc)
            end if

            if ( l_chemical_conv ) call lo2r_one%transp_lm2r(xi_ave_LMloc, xi_ave_Rloc)
            if ( l_phase_field ) call lo2r_one%transp_lm2r(phi_ave_LMloc, phi_ave_Rloc)

            !----- Outer core:
            if ( .not. l_batched_shts ) then
               do nR=nRstart,nRstop
                  if ( l_mag ) then
                     call torpol_to_spat(b_ave_Rloc(:,nR), db_ave_Rloc(:,nR), &
                          &              aj_ave_Rloc(:,nR), Br, Bt, Bp, l_R(nR))
                  end if
                  call torpol_to_spat(w_ave_Rloc(:,nR), dw_ave_Rloc(:,nR), &
                       &              z_ave_Rloc(:,nR), Vr, Vt, Vp, l_R(nR))
                  call scal_to_spat(sht_l, p_ave_Rloc(:,nR), Prer, l_R(nR))
                  if ( l_heat ) then
                     call scal_to_spat(sht_l, s_ave_Rloc(:,nR), Sr, l_R(nR))
                  end if
                  if ( l_chemical_conv ) then
                     call scal_to_spat(sht_l, xi_ave_Rloc(:,nR), Xir, l_R(nR))
                  end if
                  if ( l_phase_field ) then
                     call scal_to_spat(sht_l, phi_ave_Rloc(:,nR), Phir, l_R(nR))
                  end if
#ifdef WITH_MPI
                  call graphOut_mpi(nR, Vr, Vt, Vp, Br, Bt, Bp, Sr, Prer, Xir, Phir, &
                       &            n_graph_handle)
#else
                  call graphOut(nR, Vr, Vt, Vp, Br, Bt, Bp, Sr, Prer, Xir, Phir, &
                       &        n_graph_handle)
#endif
               end do
            else
               if ( l_mag ) then
                  call torpol_to_spat(b_ave_Rloc, db_ave_Rloc, &
                       &              aj_ave_Rloc, Br, Bt, Bp, l_R(1))
               end if
               call torpol_to_spat(w_ave_Rloc, dw_ave_Rloc, &
                    &              z_ave_Rloc, Vr, Vt, Vp, l_R(1))
               call scal_to_spat(sht_l, p_ave_Rloc, Prer, l_R(1))
               if ( l_heat ) then
                  call scal_to_spat(sht_l, s_ave_Rloc, Sr, l_R(1))
               end if
               if ( l_chemical_conv ) then
                  call scal_to_spat(sht_l, xi_ave_Rloc, Xir, l_R(1))
               end if
               if ( l_phase_field ) then
                  call scal_to_spat(sht_l, phi_ave_Rloc, Phir, l_R(1))
               end if
               do nR=nRstart,nRstop
#ifdef WITH_MPI
                  call graphOut_mpi(nR, Vr, Vt, Vp, Br, Bt, Bp, Sr, Prer, Xir, Phir, &
                       &            n_graph_handle)
#else
                  call graphOut(nR, Vr, Vt, Vp, Br, Bt, Bp, Sr, Prer, Xir, Phir, &
                       &        n_graph_handle)
#endif
               end do

            end if

            !----- Inner core: Transform is included in graphOut_IC!
            if ( l_mag .and. n_r_ic_max > 0 ) then
               if ( l_cond_ic ) then
                  call gather_all_from_lo_to_rank0(gt_IC,b_ic_ave,b_ic_ave_global)
                  call gather_all_from_lo_to_rank0(gt_IC,db_ic_ave,db_ic_ave_global)
                  call gather_all_from_lo_to_rank0(gt_IC,aj_ic_ave,aj_ic_ave_global)
               else
                  if ( rank == 0 ) then
                     db_ic_ave_global(:,1)=zero
                     aj_ic_ave_global(:,1)=zero
                  end if
               end if

               call graphOut_IC(b_ic_ave_global,db_ic_ave_global,   &
                    &           aj_ic_ave_global,bICB,n_graph_handle)
            end if

            call close_graph_file(n_graph_handle)

            !----- Write info about graph-file into STDOUT and log-file:
            if ( l_stop_time .and. rank == 0 ) then
               write(n_log_file,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')
            end if
         end if

         !--- Store time averaged poloidal magnetic coeffs at cmb
         if ( l_mag) then
            outFile='B_coeff_cmb_ave.'//tag
            nOut   =93
            n_cmb_sets=-1
            !call write_Bcmb(time,b(1,n_r_cmb),lm_max,l_max,           &
            !     &           l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
            call write_Bcmb(time,b_ave_LMloc(:,n_r_cmb),l_max_cmb, &
                 &          n_cmb_sets,outFile,nOut)
         end if

         !--- Store potentials of averaged field:
         !    dw_ave_LMloc and db_ave_LMloc used as work arrays here.
         nPotSets=0
#ifdef WITH_MPI
         call write_Pot_mpi(time,w_ave_Rloc,z_ave_Rloc,b_ic_ave,aj_ic_ave,nPotSets,  &
              &             'V_lmr_ave.',omega_ma,omega_ic)
         if ( l_mag) then
            call write_Pot_mpi(time,b_ave_Rloc,aj_ave_Rloc,b_ic_ave,aj_ic_ave,nPotSets,&
                 &             'B_lmr_ave.',omega_ma,omega_ic)
         end if
         if ( l_heat ) then
            call write_Pot_mpi(time,s_ave_Rloc,z_ave_Rloc,b_ic_ave,aj_ic_ave,nPotSets, &
                 &             'T_lmr_ave.',omega_ma,omega_ic)
         end if
         if ( l_chemical_conv ) then
            call write_Pot_mpi(time,xi_ave_Rloc,z_ave_Rloc,b_ic_ave,aj_ic_ave,nPotSets,&
                 &             'Xi_lmr_ave.',omega_ma,omega_ic)
         end if
#else
         call write_Pot(time,w_ave_LMloc,z_ave_LMloc,b_ic_ave,aj_ic_ave,nPotSets,      &
              &        'V_lmr_ave.',omega_ma,omega_ic)
         if ( l_mag) then
            call write_Pot(time,b_ave_LMloc,aj_ave_LMloc,b_ic_ave,aj_ic_ave,nPotSets,  &
                 &        'B_lmr_ave.',omega_ma,omega_ic)
         end if
         if ( l_heat ) then
            call write_Pot(time,s_ave_LMloc,z_ave_LMloc,b_ic_ave,aj_ic_ave,nPotSets,   &
                 &        'T_lmr_ave.',omega_ma,omega_ic)
         end if
         if ( l_chemical_conv ) then
            call write_Pot(time,xi_ave_LMloc,z_ave_LMloc,b_ic_ave,aj_ic_ave,nPotSets,  &
                 &        'Xi_lmr_ave.',omega_ma,omega_ic)
         end if
#endif

         if ( rank==0 .and. l_save_out ) close(n_log_file)

         !--- Store checkpoint file
#ifdef WITH_MPI
         call store_mpi(simtime,tscheme,-1,l_stop_time,.false.,.true.,        &
              &         w_ave_Rloc,z_ave_Rloc,p_ave_Rloc,s_ave_Rloc,          &
              &         xi_ave_Rloc,phi_ave_Rloc,b_ave_Rloc,aj_ave_Rloc,      &
              &         b_ic_ave,aj_ic_ave,dwdt,dzdt,dpdt,dsdt,dxidt,dphidt,  &
              &         dbdt,djdt,dbdt_ic,djdt_ic,domega_ma_dt,domega_ic_dt)
#else
         call store(simtime,tscheme,-1,l_stop_time,.false.,.true.,        &
              &     w_ave_LMloc,z_ave_LMloc,p_ave_LMloc,s_ave_LMloc,      &
              &     xi_ave_LMloc,phi_ave_LMloc,b_ave_LMloc,aj_ave_LMloc,  &
              &     b_ic_ave,aj_ic_ave,dwdt,dzdt,dpdt,dsdt,dxidt,dphidt,  &
              &     dbdt,djdt,dbdt_ic,djdt_ic,domega_ma_dt,domega_ic_dt)
#endif

         ! now correct the stored average fields by the factor which has been
         ! applied before
         if ( l_conv ) then
            w_ave_LMloc(:,:)=w_ave_LMloc(:,:)*time_norm
            z_ave_LMloc(:,:)=z_ave_LMloc(:,:)*time_norm
            p_ave_LMloc(:,:)=p_ave_LMloc(:,:)*time_norm
         end if
         if ( l_chemical_conv ) phi_ave_LMloc(:,:)=phi_ave_LMloc(:,:)*time_norm
         if ( l_heat ) s_ave_LMloc(:,:)=s_ave_LMloc(:,:)*time_norm
         if ( l_chemical_conv ) xi_ave_LMloc(:,:)=xi_ave_LMloc(:,:)*time_norm
         if ( l_mag ) then
            b_ave_LMloc(:,:) =b_ave_LMloc(:,:)*time_norm
            aj_ave_LMloc(:,:)=aj_ave_LMloc(:,:)*time_norm
            if ( l_cond_ic ) then
               b_ic_ave(:,:) =b_ic_ave(:,:)*time_norm
               aj_ic_ave(:,:)=aj_ic_ave(:,:)*time_norm
            end if
         end if

      end if ! last time step ?

   end subroutine fields_average
!------------------------------------------------------------------------------
end module fields_average_mod
