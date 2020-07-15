module fields_average_mod
   !
   ! This module is used when one wants to store time-averaged quantities
   !

   use truncation
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use fields, only: work_LMdist
   use radial_functions, only: chebt_ic, chebt_ic_even, r, dr_fac_ic, &
       &                       rscheme_oc, l_R
   use logic, only: l_mag, l_conv, l_save_out, l_heat, l_cond_ic, &
       &            l_chemical_conv
   use kinetic_energy, only: get_e_kin
   use magnetic_energy, only: get_e_mag
   use output_data, only: tag, n_log_file, log_file, n_graphs, l_max_cmb
   use sht, only: torpol_to_spat_loc, scal_to_spat_loc
   use constants, only: zero, vol_oc, vol_ic, one
   use communications, only: gather_from_mlo_to_master
   use out_coeff, only: write_Bcmb, write_Pot
   use spectra, only: spectrum, spectrum_temp
   use graphOut_mod, only: graphOut, graphOut_IC, n_graph_file, graphOut_header
   use radial_der_even, only: get_drNS_even, get_ddrNS_even
   use radial_der, only: get_dr
   use fieldsLast, only: dwdt_dist, dpdt_dist, dzdt_dist, dsdt_dist, dxidt_dist, &
       &                 dbdt_dist, djdt_dist, dbdt_ic_dist, djdt_ic_dist,       &
       &                 domega_ma_dt, domega_ic_dt, lorentz_torque_ic_dt,       &
       &                 lorentz_torque_ma_dt
   use storeCheckPoints, only: store
   use time_schemes, only: type_tscheme

   implicit none

   private

   complex(cp), allocatable :: w_ave(:,:)
   complex(cp), allocatable :: z_ave(:,:)
   complex(cp), allocatable :: s_ave(:,:)
   complex(cp), allocatable :: xi_ave(:,:)
   complex(cp), allocatable :: p_ave(:,:)
   complex(cp), allocatable :: b_ave(:,:)
   complex(cp), allocatable :: aj_ave(:,:)
   complex(cp), allocatable :: b_ic_ave(:,:)
   complex(cp), allocatable :: aj_ic_ave(:,:)
   ! on coord_r 0 we also allocate the following fields
   complex(cp), allocatable :: b_ave_global(:), bICB(:)
   complex(cp), allocatable :: db_ave_global(:), aj_ave_global(:)
   complex(cp), allocatable :: w_ave_global(:), dw_ave_global(:)
   complex(cp), allocatable :: z_ave_global(:), s_ave_global(:)
   complex(cp), allocatable :: p_ave_global(:), xi_ave_global(:)

   public :: initialize_fields_average_mod, fields_average, &
   &         finalize_fields_average_mod


contains

   subroutine initialize_fields_average_mod

      allocate( w_ave(n_mlo_loc,n_r_max) )
      allocate( z_ave(n_mlo_loc,n_r_max) )
      allocate( s_ave(n_mlo_loc,n_r_max) )
      allocate( p_ave(n_mlo_loc,n_r_max) )
      allocate( b_ave(n_mloMag_loc,n_r_maxMag) )
      allocate( aj_ave(n_mloMag_loc,n_r_maxMag) )
      bytes_allocated = bytes_allocated+4*n_mlo_loc*n_r_max*SIZEOF_DEF_COMPLEX+&
      &                 2*n_mloMag_loc*n_r_maxMag*SIZEOF_DEF_COMPLEX
      allocate( b_ic_ave(n_mloMag_loc,n_r_ic_max) )
      allocate( aj_ic_ave(n_mloMag_loc,n_r_ic_max) )
      bytes_allocated = bytes_allocated+2*n_mloMag_loc*n_r_ic_max*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( xi_ave(n_mlo_loc,n_r_max) )
         bytes_allocated = bytes_allocated+n_mlo_loc*n_r_max*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_ave(1,1) )
      end if

      if ( l_master_rank ) then
         allocate( bICB(1:lm_max) )
         allocate( b_ave_global(1:lm_max) )
         allocate( db_ave_global(1:lm_max) )
         allocate( aj_ave_global(1:lm_max) )
         allocate( w_ave_global(1:lm_max) )
         allocate( dw_ave_global(1:lm_max) )
         allocate( z_ave_global(1:lm_max) )
         allocate( s_ave_global(1:lm_max) )
         allocate( p_ave_global(1:lm_max) )
         bytes_allocated = bytes_allocated+9*lm_max*SIZEOF_DEF_COMPLEX
         if ( l_chemical_conv ) then
            allocate( xi_ave_global(1:lm_max) )
            bytes_allocated = bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX
         end if
      else
         allocate( bICB(1) )
         allocate( b_ave_global(1) )
         allocate( db_ave_global(1) )
         allocate( aj_ave_global(1) )
         allocate( w_ave_global(1) )
         allocate( dw_ave_global(1) )
         allocate( z_ave_global(1) )
         allocate( s_ave_global(1) )
         allocate( p_ave_global(1) )
         if ( l_chemical_conv ) then
            allocate( xi_ave_global(1) )
         end if
      end if

   end subroutine initialize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine finalize_fields_average_mod

      deallocate( w_ave, z_ave, s_ave, p_ave, b_ave, aj_ave, b_ic_ave )
      deallocate( aj_ic_ave, db_ave_global, aj_ave_global, w_ave_global )
      deallocate( dw_ave_global, z_ave_global, s_ave_global, p_ave_global )
      deallocate( b_ave_global, bICB )

      if ( l_chemical_conv ) deallocate( xi_ave, xi_ave_global )

   end subroutine finalize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine fields_average(simtime,tscheme,nAve,l_stop_time,        &
      &                      time_passed,time_norm,omega_ic,omega_ma, &
      &                      w,z,p,s,xi,b,aj,b_ic,aj_ic)
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
      complex(cp),         intent(in) :: w(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: z(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: p(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: s(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: xi(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: b(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(in) :: aj(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(in) :: b_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp),         intent(in) :: aj_ic(n_mloMag_loc,n_r_ic_maxMag)

      !-- Local stuff:
      ! fields for the gathering
      complex(cp) :: b_ic_ave_global(lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: db_ic_ave_global(lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: ddb_ic_ave_global(lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: aj_ic_ave_global(lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: dj_ic_ave_global(lm_maxMag,n_r_ic_maxMag)

      !----- Time averaged fields:
      complex(cp) :: dw_ave(n_mlo_loc,n_r_max)
      complex(cp) :: ds_ave(n_mlo_loc,n_r_max)
      complex(cp) :: dxi_ave(n_mlo_loc,n_r_max)
      complex(cp) :: db_ave(n_mloMag_loc,n_r_maxMag)
      complex(cp) :: db_ic_ave(n_mloMag_loc,n_r_ic_max)
      complex(cp) :: ddb_ic_ave(n_mloMag_loc,n_r_ic_max)
      complex(cp) :: dj_ic_ave(n_mloMag_loc,n_r_ic_max)

      !----- Fields in grid space:
      real(cp) :: Br(nlat_padded,n_phi_max),Bt(nlat_padded,n_phi_max)
      real(cp) :: Bp(nlat_padded,n_phi_max),Vr(nlat_padded,n_phi_max)
      real(cp) :: Vt(nlat_padded,n_phi_max),Vp(nlat_padded,n_phi_max) 
      real(cp) :: Sr(nlat_padded,n_phi_max),Prer(nlat_padded,n_phi_max)
      real(cp) :: Xir(nlat_padded,n_phi_max)

      !----- Energies of time average field:
      real(cp) :: e_kin_p_ave,e_kin_t_ave,e_kin_p_as_ave,e_kin_t_as_ave
      real(cp) :: e_mag_p_ave,e_mag_t_ave,e_mag_p_as_ave,e_mag_t_as_ave
      real(cp) :: e_mag_p_ic_ave,e_mag_t_ic_ave,e_mag_p_as_ic_ave,e_mag_t_as_ic_ave
      real(cp) :: e_mag_os_ave,e_mag_as_os_ave,Dip,DipCMB,e_cmb,elsAnel
      real(cp) :: time,dt_norm

      integer :: nR,n_e_sets,n_spec,nOut,n_cmb_sets,nPotSets

      character(len=72) :: graph_file
      character(len=80) :: outFile

      !-- Initialise average for first time step:

      if ( nAve == 1 ) then

         !zero=zero
         if ( n_graphs > 0 ) then
            if ( l_conv ) then
               w_ave(:,:)=zero
               z_ave(:,:)=zero
               p_ave(:,:)=zero
            end if
            if ( l_heat ) then
               s_ave(:,:)=zero
            end if
            if ( l_chemical_conv ) then
               xi_ave(:,:)=zero
            end if
            if ( l_mag ) then
               b_ave(:,:) =zero
               aj_ave(:,:)=zero
               if ( l_cond_ic ) then
                  b_ic_ave(:,:) =zero
                  aj_ic_ave(:,:)=zero
               end if
            end if
         end if

      end if  ! First step

      !-- Add new time step:

      if ( l_conv ) then
         w_ave(:,:)=w_ave(:,:) + time_passed*w(:,:)
         z_ave(:,:)=z_ave(:,:) + time_passed*z(:,:)
         p_ave(:,:)=p_ave(:,:) + time_passed*p(:,:)
      end if
      if ( l_heat ) then
         s_ave(:,:)=s_ave(:,:) + time_passed*s(:,:)
      end if
      if ( l_chemical_conv ) then
         xi_ave(:,:)=xi_ave(:,:) + time_passed*xi(:,:)
      end if
      if ( l_mag ) then
         b_ave(:,:) =b_ave(:,:)  + time_passed*b(:,:)
         aj_ave(:,:)=aj_ave(:,:) + time_passed*aj(:,:)
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
            w_ave(:,:)=dt_norm*w_ave(:,:)
            z_ave(:,:)=dt_norm*z_ave(:,:)
            p_ave(:,:)=dt_norm*p_ave(:,:)
         end if
         if ( l_heat ) then
            s_ave(:,:)=dt_norm*s_ave(:,:)
         end if
         if ( l_chemical_conv ) then
            xi_ave(:,:)=dt_norm*xi_ave(:,:)
         end if
         if ( l_mag ) then
            b_ave(:,:) =dt_norm*b_ave(:,:)
            aj_ave(:,:)=dt_norm*aj_ave(:,:)
         end if
         if ( l_cond_ic ) then
            b_ic_ave(:,:) =dt_norm*b_ic_ave(:,:)
            aj_ic_ave(:,:)=dt_norm*aj_ic_ave(:,:)
         end if

         !----- Get the radial derivatives:
         call get_dr(w_ave,dw_ave,n_mlo_loc,1,n_mlo_loc,n_r_max,rscheme_oc, &
              &      nocopy=.true.)
         if ( l_mag ) then
            call get_dr(b_ave,db_ave,n_mlo_loc,1,n_mlo_loc,n_r_max,rscheme_oc, &
                 &      nocopy=.true.)
         end if
         if ( l_heat ) then
            call get_dr(s_ave,ds_ave,n_mlo_loc,1,n_mlo_loc,n_r_max,rscheme_oc, &
                 &      nocopy=.true.)
         end if
         if ( l_chemical_conv ) then
            call get_dr(xi_ave,dxi_ave,n_mlo_loc,1,n_mlo_loc,n_r_max,rscheme_oc, &
                 &      nocopy=.true.)
         end if
         if ( l_cond_ic ) then
            call get_ddrNS_even(b_ic_ave,db_ic_ave,ddb_ic_ave,n_mlo_loc,1,     &
                 &              n_mlo_loc,n_r_ic_max,n_cheb_ic_max,dr_fac_ic,  &
                 &              work_LMdist,chebt_ic, chebt_ic_even)
            call get_drNS_even(aj_ic_ave,dj_ic_ave,n_mlo_loc,1,n_mlo_loc,      &
                 &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,work_LMdist, &
                 &             chebt_ic,chebt_ic_even)
         end if

         !----- Get averaged spectra:
         !      Note: average spectra will be in file no 0
         n_spec=0
         call spectrum(n_spec,time,.false.,nAve,l_stop_time,time_passed, &
              &        time_norm,w_ave,dw_ave,z_ave,b_ave,db_ave,        &
              &        aj_ave,b_ic_ave,db_ic_ave,aj_ic_ave)

         if ( l_heat ) then
            call spectrum_temp(n_spec,time,.false.,0,l_stop_time,     &
                 &             0.0_cp,0.0_cp,s_ave,ds_ave)
         end if
         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if

         !----- Write averaged energies into log-file at end of run:
         if ( l_stop_time ) then
            !----- Calculate energies of averaged field:
            n_e_sets=1
            call get_e_kin(time,.false.,.true.,n_e_sets,w_ave,dw_ave,z_ave,   &
                 &         e_kin_p_ave,e_kin_t_ave,e_kin_p_as_ave,e_kin_t_as_ave)

            call get_e_mag(time,.false.,.true.,n_e_sets,b_ave,db_ave,aj_ave,  &
                 &         b_ic_ave,db_ic_ave,aj_ic_ave,e_mag_p_ave,          &
                 &         e_mag_t_ave,e_mag_p_as_ave,e_mag_t_as_ave,         &
                 &         e_mag_p_ic_ave,e_mag_t_ic_ave,e_mag_p_as_ic_ave,   &
                 &         e_mag_t_as_ic_ave,e_mag_os_ave,e_mag_as_os_ave,    &
                 &         e_cmb,Dip,DipCMB,elsAnel)

            if ( l_master_rank ) then
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
         ! on coord_r 0 and use the old serial output routine.

         if ( l_master_rank ) then
            graph_file='G_ave.'//tag
            open(newunit=n_graph_file, file=graph_file, status='unknown', &
            &    form='unformatted', access='stream')

            !----- Write header into graphic file:
            call graphOut_header(time)
         end if

         !-- This will be needed for the inner core
         if ( l_mag ) then
            call gather_from_mlo_to_master(b_ave(:,n_r_icb),bICB)
         end if

         !----- Outer core:
         do nR=1,n_r_max
            if ( l_mag ) then
               call gather_from_mlo_to_master(b_ave(:,nR),b_ave_global)
               call gather_from_mlo_to_master(db_ave(:,nR),db_ave_global)
               call gather_from_mlo_to_master(aj_ave(:,nR),aj_ave_global)
            end if
            call gather_from_mlo_to_master(w_ave(:,nR),w_ave_global)
            call gather_from_mlo_to_master(dw_ave(:,nR),dw_ave_global)
            call gather_from_mlo_to_master(z_ave(:,nR),z_ave_global)
            call gather_from_mlo_to_master(p_ave(:,nR),p_ave_global)
            if ( l_heat ) then
               call gather_from_mlo_to_master(s_ave(:,nR),s_ave_global)
            end if
            if ( l_chemical_conv ) then
               call gather_from_mlo_to_master(xi_ave(:,nR),xi_ave_global)
            end if

            if ( l_master_rank ) then
               if ( l_mag ) then
                  call torpol_to_spat_loc(b_ave_global, db_ave_global, &
                       &              aj_ave_global, Br, Bt, Bp, l_R(nR))
               end if
               call torpol_to_spat_loc(w_ave_global, dw_ave_global, &
                    &              z_ave_global, Vr, Vt, Vp, l_R(nR))
               call scal_to_spat_loc(p_ave_global, Prer, l_R(nR))
               call scal_to_spat_loc(s_ave_global, Sr, l_R(nR))
               if ( l_chemical_conv ) then
                  call scal_to_spat_loc(xi_ave_global, Xir, l_R(nR))
               end if
               call graphOut(nR, Vr, Vt, Vp, Br, Bt, Bp, Sr, Prer, Xir)
            end if
         end do

         !----- Inner core: Transform is included in graphOut_IC!
         if ( l_mag .and. n_r_ic_max > 0 ) then
            do nR=1,n_r_ic_max 
               call gather_from_mlo_to_master(b_ic_ave(:,nR),b_ic_ave_global(:,nR))
               call gather_from_mlo_to_master(db_ic_ave(:,nR),db_ic_ave_global(:,nR))
               call gather_from_mlo_to_master(ddb_ic_ave(:,nR),ddb_ic_ave_global(:,nR))
               call gather_from_mlo_to_master(aj_ic_ave(:,nR),aj_ic_ave_global(:,nR))
               call gather_from_mlo_to_master(dj_ic_ave(:,nR),dj_ic_ave_global(:,nR))
            end do

            if ( l_master_rank ) then
               call graphOut_IC(b_ic_ave_global,db_ic_ave_global,   &
                    &           aj_ic_ave_global,bICB,l_avg=.true.)
            end if
         end if

         if ( l_master_rank ) close(n_graph_file)  ! close graphic output file !

         !----- Write info about graph-file into STDOUT and log-file:
         if ( l_stop_time .and. l_master_rank )  &
            &  write(n_log_file,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')

         !--- Store time averaged poloidal magnetic coeffs at cmb
         if ( l_mag) then
            outFile='B_coeff_cmb_ave.'//tag
            nOut   =93
            n_cmb_sets=-1
            call write_Bcmb(time,b_ave(:,n_r_cmb),l_max_cmb,n_cmb_sets,outFile,nOut)
         end if

         !--- Store potentials of averaged field:
         !    dw_ave and db_ave used as work arrays here.
         nPotSets=-1
         call write_Pot(time,w_ave,z_ave,b_ic_ave,aj_ic_ave,nPotSets,      &
              &        'V_lmr_ave.',omega_ma,omega_ic)
         if ( l_mag) then
            call write_Pot(time,b_ave,aj_ave,b_ic_ave,aj_ic_ave,nPotSets,  &
                 &        'B_lmr_ave.',omega_ma,omega_ic)
         end if
         if ( l_heat ) then
            call write_Pot(time,s_ave,z_ave,b_ic_ave,aj_ic_ave,nPotSets,   &
                 &        'T_lmr_ave.',omega_ma,omega_ic)
         end if
         if ( l_chemical_conv ) then
            call write_Pot(time,xi_ave,z_ave,b_ic_ave,aj_ic_ave,nPotSets,  &
                 &        'Xi_lmr_ave.',omega_ma,omega_ic)
         end if

         if ( l_save_out ) close(n_log_file)

         !--- Store checkpoint file
         call store(simtime,tscheme,-1,l_stop_time,.false.,.true.,            &
              &     w_ave,z_ave,p_ave,s_ave,xi_ave,b_ave,aj_ave,b_ic_ave,     &
              &     aj_ic_ave,dwdt_dist,dzdt_dist,dpdt_dist,dsdt_dist,        &
              &     dxidt_dist,dbdt_dist,djdt_dist,dbdt_ic_dist,djdt_ic_dist, &
              &      domega_ma_dt,domega_ic_dt,lorentz_torque_ma_dt,          &
              &     lorentz_torque_ic_dt)

         ! now correct the stored average fields by the factor which has been
         ! applied before
         if ( l_conv ) then
            w_ave(:,:)=w_ave(:,:)*time_norm
            z_ave(:,:)=z_ave(:,:)*time_norm
            p_ave(:,:)=p_ave(:,:)*time_norm
         end if
         if ( l_heat ) s_ave(:,:)=s_ave(:,:)*time_norm
         if ( l_chemical_conv ) xi_ave(:,:)=xi_ave(:,:)*time_norm
         if ( l_mag ) then
            b_ave(:,:) =b_ave(:,:)*time_norm
            aj_ave(:,:)=aj_ave(:,:)*time_norm
         end if
         if ( l_cond_ic ) then
            b_ic_ave(:,:) =b_ic_ave(:,:)*time_norm
            aj_ic_ave(:,:)=aj_ic_ave(:,:)*time_norm
         end if

      end if ! last time step ?

   end subroutine fields_average
!------------------------------------------------------------------------------
end module fields_average_mod
