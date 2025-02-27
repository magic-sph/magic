module start_fields
   !
   ! This module is used to set-up the initial starting fields. They can be obtained
   ! by reading a starting checkpoint file or by setting some starting conditions.
   !

   use truncation
   use precision_mod
   use parallel_mod
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use communications, only: lo2r_one
   use radial_functions, only: rscheme_oc, r, or1, alpha0, dLtemp0,      &
       &                       dLalpha0, beta, orho1, temp0, rho0, r_cmb,&
       &                       otemp1, ogrun, dentropy0, dxicond, r_icb
   use physical_parameters, only: interior_model, epsS, impS, n_r_LCR,   &
       &                          ktopv, kbotv, LFfac, imagcon, ThExpNb, &
       &                          ViscHeatFac, impXi
   use num_param, only: dtMax, alpha
   use special, only: lGrenoble, ampForce
   use output_data, only: log_file, n_log_file
   use blocking, only: lo_map, st_map, llm, ulm, ulmMag, llmMag
   use logic, only: l_conv, l_mag, l_cond_ic, l_heat, l_SRMA, l_SRIC,    &
       &            l_mag_kin, l_mag_LF, l_temperature_diff, l_onset,    &
       &            l_chemical_conv, l_anelastic_liquid, l_save_out,     &
       &            l_parallel_solve, l_mag_par_solve, l_phase_field,    &
       &            l_single_matrix, l_non_adia, l_tidal
   use init_fields, only: l_start_file, init_s1, init_b1, tops, pt_cond,  &
       &                  initV, initS, initB, initXi, ps_cond,           &
       &                  start_file, init_xi1, topxi, xi_cond, omega_ic1,&
       &                  omega_ma1, initPhi, init_phi, initF, initTidal  &
       &                 ,topv, force_z_vol, botv, pertur_w
   use fields ! The entire module is required
   use fieldsLast ! The entire module is required
   use timing, only: timer_type
   use constants, only: zero, c_lorentz_ma, c_lorentz_ic, osq4pi, one, two, sq4pi&
        ,c_z10_omega_ma,c_z10_omega_ic
   use useful, only: cc2real, logWrite
   use radial_der, only: get_dr, exch_ghosts, bulk_to_ghost
   use readCheckPoints, only: readStartFields_old, readStartFields
   use time_schemes, only: type_tscheme
#ifdef WITH_MPI
   use readCheckPoints, only: readStartFields_mpi
#endif
   use updateWPS_mod, only: get_single_rhs_imp
   use updateWP_mod, only: get_pol_rhs_imp, get_pol_rhs_imp_ghost, w_ghost, &
       &                   fill_ghosts_W, p0_ghost
   use updateS_mod, only: get_entropy_rhs_imp, get_entropy_rhs_imp_ghost, s_ghost, &
       &                  fill_ghosts_S
   use updateXI_mod, only: get_comp_rhs_imp, get_comp_rhs_imp_ghost, xi_ghost, &
       &                   fill_ghosts_Xi
   use updatePhi_mod, only: get_phase_rhs_imp, get_phase_rhs_imp_ghost, phi_ghost, &
       &                   fill_ghosts_Phi
   use updateZ_mod, only: get_tor_rhs_imp, get_tor_rhs_imp_ghost, z_ghost, &
       &                  fill_ghosts_Z
   use updateB_mod, only: get_mag_rhs_imp, get_mag_ic_rhs_imp, b_ghost, aj_ghost, &
       &                  get_mag_rhs_imp_ghost, fill_ghosts_B


   implicit none

   private

   real(cp), public :: topcond ! Conducting heat flux at the outer boundary
   real(cp), public :: botcond ! Conducting heat flux at the inner boundary
   real(cp), public :: deltacond ! Temperature or entropy difference between boundaries
   real(cp), public :: topxicond ! Conducting mass flux at the outer boundary
   real(cp), public :: botxicond ! Conducting mass flux at the inner boundary
   real(cp), public :: deltaxicond ! Composition difference between boundaries

   public :: getStartFields

contains

   subroutine getStartFields(time,tscheme,n_time_step)
      !
      !  Purpose of this subroutine is to initialize the fields and
      !  other auxiliary parameters.
      !

      !---- Output variables:
      real(cp),            intent(out) :: time ! Time of the restart
      integer,             intent(out) :: n_time_step ! Number of past iterations
      class(type_tscheme), intent(inout) :: tscheme

      !-- Local variables:
      integer :: l, m, lm, n_r, filehandle, lm00, l1m0, ierr
      character(len=76) :: message
      real(cp) :: sEA,sES,sAA,xiEA,xiES,xiAA,topval,botval
      real(cp) :: s0(n_r_max),p0(n_r_max),ds0(n_r_max),dp0(n_r_max)
      logical :: lMat
      type(timer_type) :: t_reader
      complex(cp) :: local_topv(0:l_max,0:m_max)
      complex(cp) :: local_botv(0:l_max,0:m_max)

      call t_reader%initialize()

      !---- Computations for the Nusselt number if we are anelastic
      !     Can be done before setting the fields
      if ( l_heat ) then

         if ( rank == 0 ) open(newunit=filehandle, file='pscond.dat')

         if ( l_anelastic_liquid ) then ! temperature

            call pt_cond(s0,p0)

            if ( rank == 0 ) then
               do n_r=1,n_r_max
                  write(filehandle,'(5ES20.12)') r(n_r), osq4pi*otemp1(n_r)*&
                  &            (s0(n_r)-ViscHeatFac*ThExpNb*alpha0(n_r)*    &
                  &            temp0(n_r)*orho1(n_r)*p0(n_r)),              &
                  &            osq4pi*p0(n_r), osq4pi*s0(n_r),              &
                  &            osq4pi*alpha0(n_r)*(-rho0(n_r)*s0(n_r)+      &
                  &            ViscHeatFac*ThExpNb*(alpha0(n_r)*temp0(n_r)  &
                  &            +ogrun(n_r))*p0(n_r))
               end do
            end if

            call get_dr(s0,ds0,n_r_max,rscheme_oc)

            if ( l_temperature_diff ) then ! temperature diffusion
               topcond=-osq4pi*ds0(1)
               botcond=-osq4pi*ds0(n_r_max)
               deltacond=osq4pi*(s0(n_r_max)-s0(1))
            else ! entropy diffusion
               call get_dr(p0,dp0,n_r_max,rscheme_oc)

               topcond = -osq4pi*(otemp1(1)*( -dLtemp0(1)*s0(1)+ds0(1))- &
               &        ViscHeatFac*ThExpNb*alpha0(1)*orho1(1)*(         &
               &        (dLalpha0(1)-beta(1))*p0(1) + dp0(1)) )
               botcond = -osq4pi*(otemp1(n_r_max)*( -dLtemp0(n_r_max)*    &
               &                   s0(n_r_max) + ds0(n_r_max))-           &
               &      ViscHeatFac*ThExpNb*alpha0(n_r_max)*orho1(n_r_max)*(&
               &             (dLalpha0(n_r_max)-beta(n_r_max))*           &
               &                   p0(n_r_max) + dp0(n_r_max)) )
               deltacond=osq4pi*(otemp1(n_r_max)*s0(n_r_max)-otemp1(1)*s0(1)- &
               &               ViscHeatFac*ThExpNb*( alpha0(n_r_max)*         &
               &               orho1(n_r_max)*p0(n_r_max)-                    &
               &               alpha0(1)*orho1(1)*p0(1)) )
            end if

         else ! entropy is the thermodynamic variable

            call ps_cond(s0,p0)

            if ( rank == 0 ) then
               do n_r=1,n_r_max
                  write(filehandle,'(5ES20.12)') r(n_r), s0(n_r)*osq4pi, &
                  &            p0(n_r)*osq4pi, osq4pi*temp0(n_r)*(       &
                  &            s0(n_r)+alpha0(n_r)*orho1(n_r)*p0(n_r)*   &
                  &            ThExpNb*ViscHeatFac), osq4pi*alpha0(n_r)* &
                  &            ThExpNb*(-rho0(n_r)*temp0(n_r)*s0(n_r)+   &
                  &            ViscHeatFac*ogrun(n_r)*p0(n_r))
               end do
            end if

            call get_dr(s0,ds0,n_r_max,rscheme_oc)

            if ( l_temperature_diff ) then
               call get_dr(p0,dp0,n_r_max,rscheme_oc)

               topcond = -osq4pi*temp0(1)*( dLtemp0(1)*s0(1)+ds0(1)+   &
               &        ViscHeatFac*ThExpNb*alpha0(1)*orho1(1)*(       &
               &        (dLalpha0(1)+dLtemp0(1)-beta(1))*p0(1) +       &
               &                            dp0(1)) )
               botcond = -osq4pi*temp0(n_r_max)*( dLtemp0(n_r_max)*       &
               &                   s0(n_r_max) + ds0(n_r_max)+            &
               &      ViscHeatFac*ThExpNb*alpha0(n_r_max)*orho1(n_r_max)*(&
               &    (dLtemp0(n_r_max)+dLalpha0(n_r_max)-beta(n_r_max))*   &
               &                   p0(n_r_max) + dp0(n_r_max)) )
               deltacond=osq4pi*(temp0(n_r_max)*s0(n_r_max)-temp0(1)*s0(1)+ &
               &               ViscHeatFac*ThExpNb*( alpha0(n_r_max)*       &
               &               temp0(n_r_max)*orho1(n_r_max)*p0(n_r_max)-   &
               &               alpha0(1)*temp0(1)*orho1(1)*p0(1)) )
            else ! entropy diffusion
               topcond=-osq4pi*ds0(1)
               botcond=-osq4pi*ds0(n_r_max)
               deltacond=osq4pi*(s0(n_r_max)-s0(1))
            end if

         end if

         if ( l_onset .and. ( .not. l_non_adia ) ) dentropy0(:) = ds0(:) * osq4pi

         if ( rank == 0 ) close(filehandle)

      else
         topcond  =0.0_cp
         botcond  =0.0_cp
         deltacond=0.0_cp
      end if

      if ( l_chemical_conv ) then
         call xi_cond(s0)
         call get_dr(s0,ds0,n_r_max,rscheme_oc)
         topxicond=-osq4pi*ds0(1)
         botxicond=-osq4pi*ds0(n_r_max)
         deltaxicond=osq4pi*(s0(n_r_max)-s0(1))
         if ( l_onset ) dxicond(:)=ds0(:) * osq4pi
      else
         topxicond  =0.0_cp
         botxicond  =0.0_cp
         deltaxicond=0.0_cp
      end if

      !-- Start with setting fields to zero:
      !   Touching the fields with the appropriate processor
      !   for the LM-distribute parallel region (LMLoop) makes
      !   sure that they are located close the individual
      !   processors in memory:

      if ( l_start_file ) then

         call t_reader%start_count()
         if ( index(start_file, 'rst_') /= 0 ) then
            call readStartFields_old( w_LMloc,dwdt,z_LMloc,dzdt,p_LMloc,dpdt,      &
                 &                    s_LMloc,dsdt,xi_LMloc,dxidt,phi_LMloc,       &
                 &                    dphidt,b_LMloc,dbdt,aj_LMloc,djdt,           &
                 &                    b_ic_LMloc,dbdt_ic,aj_ic_LMloc,djdt_ic,      &
                 &                    omega_ic,omega_ma,domega_ic_dt,domega_ma_dt, &
                 &                    time,tscheme,n_time_step )
         else
#ifdef WITH_MPI
            call readStartFields_mpi( w_LMloc,dwdt,z_LMloc,dzdt,p_LMloc,dpdt,   &
                 &                    s_LMloc,dsdt,xi_LMloc,dxidt,phi_LMloc,    &
                 &                    dphidt,b_LMloc,dbdt,aj_LMloc,djdt,        &
                 &                    b_ic_LMloc,dbdt_ic,aj_ic_LMloc,djdt_ic,   &
                 &                    omega_ic,omega_ma,domega_ic_dt,           &
                 &                    domega_ma_dt,time,tscheme,n_time_step )
#else
            call readStartFields( w_LMloc,dwdt,z_LMloc,dzdt,p_LMloc,dpdt,s_LMloc,&
                 &                dsdt,xi_LMloc,dxidt,phi_LMloc,dphidt,b_LMloc,  &
                 &                dbdt,aj_LMloc,djdt,b_ic_LMloc,dbdt_ic,         &
                 &                aj_ic_LMloc,djdt_ic,omega_ic,omega_ma,         &
                 &                domega_ic_dt,domega_ma_dt,time,tscheme,n_time_step )
#endif
         end if

         call t_reader%stop_count()
         if ( rank == 0 .and. l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if
         call t_reader%finalize('! Time taken to read the checkpoint file:', &
              &                 n_log_file)
         if ( rank == 0 .and. l_save_out ) close(n_log_file)

         if ( tscheme%dt(1) > 0.0_cp ) then
            if ( rank==0 ) write(message,'(''! Using old time step:'',ES16.6)') tscheme%dt(1)
         else
            tscheme%dt(1)=dtMax
            if ( rank==0 ) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
         end if

         if ( .not. l_heat ) s_LMloc(:,:)=zero

         if (ktopv == 3 .and. .not. pertur_w==-1) then
            l1m0=lo_map%lm2(1,0)
            local_topv(:,:)= zero
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if(lm == l1m0) then
                  local_topv(l,m)=local_topv(l,m)+ z_LMloc(lm,n_r_cmb)*c_z10_omega_ma
               else
                  local_topv(l,m)=local_topv(l,m)+ z_LMloc(lm,n_r_cmb)!/rscheme_oc%rnorm
               end if
            end do
#ifdef WITH_MPI
             call MPI_Allreduce(local_topv,topv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
         end if
         if (kbotv == 3) then
            l1m0=lo_map%lm2(1,0)
            local_botv(:,:)= zero
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if(lm == l1m0) then
                  local_botv(l,m)=local_botv(l,m)+ z_LMloc(lm,n_r_icb)*c_z10_omega_ic
               else
                  local_botv(l,m)=local_botv(l,m)+ z_LMloc(lm,n_r_icb)!/rscheme_oc%rnorm
               end if
            end do
#ifdef WITH_MPI
             call MPI_Allreduce(local_botv,botv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
         end if

         lm00 = lo_map%lm2(0,0)
         !-- In case temperature/entropy was imposed on both boundaries, simply
         !-- translate it to (0,1) instead of legacy contrast
         if ( l_heat .and. llm <= lm00 .and. ulm >= lm00 ) then
            topval=-r_icb**2/(r_icb**2+r_cmb**2)*sq4pi
            botval= r_cmb**2/(r_icb**2+r_cmb**2)*sq4pi
            if ( abs(s_LMloc(lm00,n_r_cmb)-topval) <= 10000.0_cp*epsilon(one) .and. &
            &    abs(s_LMloc(lm00,n_r_icb)-botval) <= 10000.0_cp*epsilon(one) ) then
               s_LMloc(lm00,:)=s_LMloc(lm00,:)-cmplx(topval,0.0_cp,cp)
            end if
         end if

         !-- In case chemical composition was imposed on both boundaries, simply
         !-- translate it to (0,1) instead of legacy contrast
         if ( l_chemical_conv .and. llm <= lm00 .and. ulm >= lm00 ) then
            topval=-r_icb**2/(r_icb**2+r_cmb**2)*sq4pi
            botval= r_cmb**2/(r_icb**2+r_cmb**2)*sq4pi
            if ( abs(xi_LMloc(lm00,n_r_cmb)-topval) <= 10000.0_cp*epsilon(one) .and. &
            &    abs(xi_LMloc(lm00,n_r_icb)-botval) <= 10000.0_cp*epsilon(one) ) then
               xi_LMloc(lm00,:)=xi_LMloc(lm00,:)-cmplx(topval,0.0_cp,cp)
            end if
         end if

      else ! If there's no restart file

         ! Initialize with zero
         if ( l_conv .or. l_mag_kin ) then
            w_LMloc(:,:)=zero
            z_LMloc(:,:)=zero
            p_LMloc(:,:)=zero
         end if
         if ( l_heat ) s_LMloc(:,:)=zero
         if ( l_chemical_conv ) xi_LMloc(:,:)=zero
         if ( l_phase_field ) phi_LMloc(:,:)=zero
         if ( l_mag ) then
            b_LMloc(:,:) =zero
            aj_LMloc(:,:)=zero
         end if
         if ( l_cond_ic ) then
            b_ic_LMloc(:,:) =zero
            aj_ic_LMloc(:,:)=zero
         end if

         time         =0.0_cp
         tscheme%dt(:)=dtMax !dt   =0.01*dtMax!0.001*dtMax
         n_time_step  =0
         if (rank == 0) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
      end if
      call logWrite(message)

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

      !-- Initialize the weights of the time scheme
      call tscheme%set_weights(lMat)

      !-- Initialize/add fields
      !----- Initialize/add magnetic field:
      if ( ( imagcon /= 0 .or. init_b1 /= 0 .or. lGrenoble ) &
      &      .and. ( l_mag .or. l_mag_LF ) ) then
         call initB(b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc)
      end if

      !----- Initialize/add velocity, set IC and ma rotation:
      if ( l_conv .or. l_mag_kin .or. l_SRIC .or. l_SRMA ) then
         call initV(w_LMloc,z_LMloc,omega_ic,omega_ma)
      end if

      !----- Initialize/add entropy:
      if ( ( init_s1 /= 0 .or. impS /= 0 ) .and. l_heat ) call initS(s_LMloc,p_LMloc)

      !----- Initialize/add chemical convection:
      if ( ( init_xi1 /= 0 .or. impXi /= 0 ) .and. l_chemical_conv ) call initXi(xi_LMloc)

      !----- Initialize/add phase field:
      if ( init_phi /= 0 .and. l_phase_field ) call initPhi(s_LMloc, phi_LMloc)

      !---- For now fiels initialized in R-distributed arrays: now transpose them if needed
      if ( l_parallel_solve ) then
         call lo2r_one%transp_lm2r(w_LMloc, w_Rloc)
         call lo2r_one%transp_lm2r(z_LMloc, z_Rloc)
         if ( l_chemical_conv ) call lo2r_one%transp_lm2r(xi_LMloc, xi_Rloc)
         if ( l_phase_field ) call lo2r_one%transp_lm2r(phi_LMloc, phi_Rloc)
         if ( l_heat ) call lo2r_one%transp_lm2r(s_LMloc, s_Rloc)
         call lo2r_one%transp_lm2r(p_LMloc, p_Rloc)
         if ( l_mag .and. l_mag_par_solve ) then
            call lo2r_one%transp_lm2r(b_LMloc, b_Rloc)
            call lo2r_one%transp_lm2r(aj_LMloc, aj_Rloc)
         end if
      end if


      !----- Assemble initial implicit terms
      if ( l_chemical_conv ) then
         if ( l_parallel_solve ) then
            call bulk_to_ghost(xi_Rloc, xi_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
            call exch_ghosts(xi_ghost, lm_max, nRstart, nRstop, 1)
            call fill_ghosts_Xi(xi_ghost)
            call get_comp_rhs_imp_ghost(xi_ghost, dxidt, 1, .true.)
         else
            call get_comp_rhs_imp(xi_LMloc, dxi_LMloc, dxidt, 1, .true.)
         end if
      end if

      if ( l_phase_field ) then
         if ( l_parallel_solve ) then
            call bulk_to_ghost(phi_Rloc, phi_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
            call exch_ghosts(phi_ghost, lm_max, nRstart, nRstop, 1)
            call fill_ghosts_Phi(phi_ghost)
            call get_phase_rhs_imp_ghost(phi_ghost, dphidt, 1, .true.)
         else
            call get_phase_rhs_imp(phi_LMloc, dphidt, 1, .true.)
         end if
      end if

      if ( l_single_matrix ) then
         call get_single_rhs_imp(s_LMloc, ds_LMloc, w_LMloc, dw_LMloc,     &
              &                  ddw_LMloc, p_LMloc, dp_LMloc, dsdt, dwdt, &
              &                  dpdt, tscheme, 1, .true., .false., time=time)
      else
         if ( l_heat ) then
            if ( l_parallel_solve ) then
               call bulk_to_ghost(s_Rloc, s_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
               call exch_ghosts(s_ghost, lm_max, nRstart, nRstop, 1)
               call fill_ghosts_S(s_ghost)
               call get_entropy_rhs_imp_ghost(s_ghost, ds_Rloc, dsdt, phi_Rloc, &
                    &                         1, .true.)
            else
               call get_entropy_rhs_imp(s_LMloc, ds_LMloc, dsdt, phi_LMloc, 1, .true.)
            end if
         end if
         if ( l_parallel_solve ) then
            call bulk_to_ghost(w_Rloc, w_ghost, 2, nRstart, nRstop, lm_max, 1, lm_max)
            call bulk_to_ghost(p_Rloc(1,:), p0_ghost, 1, nRstart, nRstop, 1, 1, 1)
            call exch_ghosts(w_ghost, lm_max, nRstart, nRstop, 2)
            call fill_ghosts_W(w_ghost, p0_ghost, .true.)
            call get_pol_rhs_imp_ghost(w_ghost, dw_Rloc, ddw_Rloc, p_Rloc, dp_Rloc,  &
                 &                     dwdt, tscheme, 1, .true., .false., .false.,   &
                 &                     dwdt%expl(:,:,1)) ! Work array
         else
            call get_pol_rhs_imp(s_LMloc, xi_LMloc, w_LMloc, dw_LMloc, ddw_LMloc,  &
                 &               p_LMloc, dp_LMloc, dwdt, dpdt, tscheme, 1, .true.,&
                 &               .false., .false., work_LMloc)
         end if
      end if
      if ( l_parallel_solve ) then
         call bulk_to_ghost(z_Rloc, z_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
         call exch_ghosts(z_ghost, lm_max, nRstart, nRstop, 1)
         call fill_ghosts_Z(z_ghost)
         call get_tor_rhs_imp_ghost(time, z_ghost, dz_Rloc, dzdt, domega_ma_dt,  &
              &                     domega_ic_dt, omega_ic, omega_ma, omega_ic1, &
              &                     omega_ma1, tscheme, 1, .true., .false.)
      else
         call get_tor_rhs_imp(time, z_LMloc, dz_LMloc, dzdt, domega_ma_dt, &
              &               domega_ic_dt, omega_ic, omega_ma, omega_ic1, &
              &               omega_ma1, tscheme, 1, .true., .false.)
      end if
      call get_tor_rhs_imp(time, z_LMloc, dz_LMloc, dzdt, domega_ma_dt, &
           &               domega_ic_dt, omega_ic, omega_ma, omega_ic1, &
           &               omega_ma1, tscheme, 1, .true., .false.)

      if ( l_mag .or. l_mag_kin  ) then
         if ( l_mag_par_solve ) then
            call bulk_to_ghost(b_Rloc, b_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
            call bulk_to_ghost(aj_Rloc, aj_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
            call exch_ghosts(aj_ghost, lm_max, nRstart, nRstop, 1)
            call exch_ghosts(b_ghost, lm_max, nRstart, nRstop, 1)
            call fill_ghosts_B(b_ghost, aj_ghost)
            call get_mag_rhs_imp_ghost(b_ghost, db_Rloc, ddb_Rloc, aj_ghost,     &
                 &                     dj_Rloc, ddj_Rloc, dbdt, djdt, tscheme, 1,&
                 &                     .true., .false.)
         else
            call get_mag_rhs_imp(b_LMloc, db_LMloc, ddb_LMloc, aj_LMloc,     &
                 &               dj_LMloc, ddj_LMloc, dbdt, djdt, tscheme, 1,&
                 &               .true., .false.)
         end if
      end if
      if ( l_cond_ic ) then
         call get_mag_ic_rhs_imp(b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc,    &
              &                  aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc,   &
              &                  dbdt_ic, djdt_ic, 1, .true.)
      end if

      !--- Get symmetry properties of tops excluding l=m=0:
      sES=0.0_cp
      sEA=0.0_cp
      sAA=0.0_cp
      do m=m_min,m_max,minc
         do l=m,l_max
            if ( l > 0 ) then
               if ( mod(l+m,2) == 0 ) then
                  sES=sES+cc2real(tops(l,m),m)
               else
                  sEA=sEA+cc2real(tops(l,m),m)
               end if
               if ( m /= 0 ) sAA=sAA+cc2real(tops(l,m),m)
            end if
         end do
      end do
      if ( sEA+sES == 0 ) then
         write(message,'(''! Only l=m=0 comp. in tops:'')')
         call logWrite(message)
      else
         sEA=sqrt(sEA/(sEA+sES))
         sAA=sqrt(sAA/(sEA+sES))
         write(message,'(''! Rel. RMS equ. asym. tops:'',ES16.6)') sEA
         call logWrite(message)
         write(message,'(''! Rel. RMS axi. asym. tops:'',ES16.6)') sAA
         call logWrite(message)
      end if

      !--- Get symmetry properties of topxi excluding l=m=0:
      if ( l_chemical_conv ) then
         xiES=0.0_cp
         xiEA=0.0_cp
         xiAA=0.0_cp
         do m=m_min,m_max,minc
            do l=m,l_max
               if ( l > 0 ) then
                  if ( mod(l+m,2) == 0 ) then
                     xiES=xiES+cc2real(topxi(l,m),m)
                  else
                     xiEA=xiEA+cc2real(topxi(l,m),m)
                  end if
                  if ( m /= 0 ) xiAA=xiAA+cc2real(topxi(l,m),m)
               end if
            end do
         end do
         if ( xiEA+xiES == 0 ) then
            write(message,'(''! Only l=m=0 comp. in topxi:'')')
            call logWrite(message)
         else
            xiEA=sqrt(xiEA/(xiEA+xiES))
            xiAA=sqrt(xiAA/(xiEA+xiES))
            write(message,'(''! Rel. RMS equ. asym. topxi:'',ES16.6)') xiEA
            call logWrite(message)
            write(message,'(''! Rel. RMS axi. asym. topxi:'',ES16.6)') xiAA
            call logWrite(message)
         end if
      end if

      if ( ampForce /= 0.0_cp ) then
         call initF(bodyForce_LMloc)
         if ( l_parallel_solve ) then
            call lo2r_one%transp_lm2r(bodyForce_LMloc, bodyForce_Rloc)
         end if
      end if

      if (l_tidal) then
         wtidal_LMloc(:,:)=zero
         dwtidal_LMloc(:,:)=zero
         ddwtidal_LMloc(:,:)=zero
         wtidal_Rloc(:,:)=zero
         dwtidal_Rloc(:,:)=zero
         ddwtidal_Rloc(:,:)=zero
         call initTidal(wtidal_LMloc,dwtidal_LMloc,ddwtidal_LMloc &
              & ,wtidal_Rloc,dwtidal_Rloc,ddwtidal_Rloc)
      end if

   end subroutine getStartFields
!------------------------------------------------------------------------------
end module start_fields
