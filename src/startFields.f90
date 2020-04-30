#include "perflib_preproc.cpp"
module start_fields

#ifdef WITH_MPI
   use mpimod
#endif
   use truncation
   use precision_mod
   use radial_functions, only: rscheme_oc, r, or1, alpha0, dLtemp0,      &
       &                       dLalpha0, beta, orho1, temp0, rho0,       &
       &                       otemp1, ogrun
   use physical_parameters, only: interior_model, epsS, impS, n_r_LCR,   &
       &                          ktopv, kbotv, LFfac, imagcon, ThExpNb, &
       &                          ViscHeatFac, impXi
   use mpi_thetap_mod, only: test_field, transform_new2old !@> TODO: remove!!
   use num_param, only: dtMax, alpha
   use special, only: lGrenoble
   use output_data, only: log_file, n_log_file
   use logic, only: l_conv, l_mag, l_cond_ic, l_heat, l_SRMA, l_SRIC,    &
       &            l_mag_kin, l_mag_LF, l_rot_ic, l_z10Mat, l_LCR,      &
       &            l_rot_ma, l_temperature_diff, l_single_matrix,       &
       &            l_chemical_conv, l_anelastic_liquid, l_save_out
   use init_fields, only: l_start_file, init_s1, init_b1, tops, pt_cond,  &
       &                  initV, initS, initB, initXi, ps_cond,           &
       &                  start_file, init_xi1, topxi, xi_cond, omega_ic1,&
       &                  omega_ma1
   use fields ! The entire module is required
   use fieldsLast ! The entire module is required
   use timing, only: timer_type
   use constants, only: zero, c_lorentz_ma, c_lorentz_ic, osq4pi, one, two
   use useful, only: cc2real, logWrite
   use radial_der, only: get_dr
   use readCheckPoints, only: readStartFields_old, readStartFields
   use time_schemes, only: type_tscheme
#ifdef WITH_MPI
   use readCheckPoints, only: readStartFields_mpi
#endif
   use updateWPS_mod, only: get_single_rhs_imp
   use updateWP_mod, only: get_pol_rhs_imp
   use updateS_mod, only: get_entropy_rhs_imp
   use updateXI_mod, only: get_comp_rhs_imp
   use updateZ_mod, only: get_rot_rates, get_tor_rhs_imp
   use updateB_mod, only: get_mag_rhs_imp, get_mag_ic_rhs_imp


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
      integer :: l, m, n_r
      character(len=76) :: message

      real(cp) :: sEA,sES,sAA
      real(cp) :: xiEA,xiES,xiAA

      real(cp) :: s0(n_r_max),p0(n_r_max),ds0(n_r_max),dp0(n_r_max)

      logical :: lMat
      type(timer_type) :: t_reader
      integer :: ierr, filehandle

      call t_reader%initialize()

      !---- Computations for the Nusselt number if we are anelastic
      !     Can be done before setting the fields
      if ( l_heat ) then

         if ( l_master_rank ) open(newunit=filehandle, file='pscond.dat')

         if ( l_anelastic_liquid ) then ! temperature

            call pt_cond(s0,p0)

            if ( l_master_rank ) then
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

            if ( l_master_rank ) then
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

         if ( l_master_rank ) close(filehandle)

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
            call readStartFields_old( w_LMdist,dwdt_dist,z_LMdist,dzdt_dist,p_LMdist,dpdt_dist,   &
                 &                    s_LMdist,dsdt_dist,xi_LMdist,dxidt_dist,b_LMdist,      &
                 &                    dbdt_dist,aj_LMdist,djdt_dist,b_ic_LMdist,dbdt_ic_dist,    &
                 &                    aj_ic_LMdist,djdt_ic_dist,omega_ic,omega_ma,    &
                 &                    domega_ic_dt,domega_ma_dt,                &
                 &                    lorentz_torque_ic_dt,lorentz_torque_ma_dt,&
                 &                    time,tscheme,n_time_step )
         else
#ifdef WITH_MPI
            call readStartFields_mpi( w_LMdist,dwdt_dist,z_LMdist,dzdt_dist,p_LMdist,dpdt_dist,   &
                 &                    s_LMdist,dsdt_dist,xi_LMdist,dxidt_dist,b_LMdist,dbdt_dist, &
                 &                    aj_LMdist,djdt_dist,b_ic_LMdist,dbdt_ic_dist,         &
                 &                    aj_ic_LMdist,djdt_ic_dist,omega_ic,omega_ma,    &
                 &                    domega_ic_dt,domega_ma_dt,                &
                 &                    lorentz_torque_ic_dt,lorentz_torque_ma_dt,&
                 &                    time,tscheme,n_time_step )
#else
            call readStartFields( w_LMdist,dwdt_dist,z_LMdist,dzdt_dist,p_LMdist,dpdt_dist,s_LMdist,&
                 &                dsdt_dist,xi_LMdist,dxidt_dist,b_LMdist,dbdt_dist,aj_LMdist,djdt_dist,&
                 &                b_ic_LMdist,dbdt_ic_dist,aj_ic_LMdist,djdt_ic_dist,        &
                 &                omega_ic,omega_ma,domega_ic_dt,domega_ma_dt,   &
                 &                lorentz_torque_ic_dt,lorentz_torque_ma_dt,     &
                 &                time,tscheme,n_time_step )
#endif

         end if
         call t_reader%stop_count()
         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if
         call t_reader%finalize('! Time taken to read the checkpoint file:', &
              &                 n_log_file)
         if ( l_save_out ) close(n_log_file)

         if ( tscheme%dt(1) > 0.0_cp ) then
            if ( l_master_rank ) write(message,'(''! Using old time step:'',ES16.6)') tscheme%dt(1)
         else
            tscheme%dt(1)=dtMax
            if ( l_master_rank ) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
         end if

         if ( .not. l_heat ) s_LMdist(:,:)=zero

      else ! If there's no restart file

         ! Initialize with zero
         if ( l_conv .or. l_mag_kin ) then
            w_LMdist(:,:)=zero
            z_LMdist(:,:)=zero
            p_LMdist(:,:)=zero
         end if
         if ( l_heat ) s_LMdist(:,:)=zero
         if ( l_chemical_conv ) xi_LMdist(:,:)=zero
         if ( l_mag ) then
            b_LMdist(:,:) =zero
            aj_LMdist(:,:)=zero
         end if
         if ( l_cond_ic ) then
            b_ic_LMdist(:,:) =zero
            aj_ic_LMdist(:,:)=zero
         end if

         time         =0.0_cp
         tscheme%dt(:)=dtMax
         n_time_step  =0
         if (l_master_rank) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
      end if
      call logWrite(message)

#ifdef WITH_MPI
      call MPI_Barrier(comm_r, ierr)
#endif

      !-- Initialize the weights of the time scheme
      call tscheme%set_weights(lMat)


      !-- Initialize/add fields
      !----- Initialize/add magnetic field:
      if ( ( imagcon /= 0 .or. init_b1 /= 0 .or. lGrenoble ) &
      &      .and. ( l_mag .or. l_mag_LF ) ) then
         call initB(b_LMdist,aj_LMdist,b_ic_LMdist,aj_ic_LMdist)
      end if

      !----- Initialize/add velocity, set IC and ma rotation:
      if ( l_conv .or. l_mag_kin .or. l_SRIC .or. l_SRMA ) then
         call initV(w_LMdist,z_LMdist,omega_ic,omega_ma)
      end if

      !----- Initialize/add entropy:
      if ( ( init_s1 /= 0 .or. impS /= 0 ) .and. l_heat ) then
         call initS(s_LMdist,p_LMdist)
      end if

      !----- Initialize/add chemical convection:
      if ( ( init_xi1 /= 0 .or. impXi /= 0 ) .and. l_chemical_conv ) then
         call initXi(xi_LMdist)
      end if

      !----- Assemble initial implicit terms
      if ( l_chemical_conv ) then
         call get_comp_rhs_imp(xi_LMdist, dxi_LMdist, dxidt_dist, 1, .true.)
      end if

      if ( l_single_matrix ) then
         call get_single_rhs_imp(s_LMdist, ds_LMdist, w_LMdist, dw_LMdist,   &
              &                  ddw_LMdist, p_LMdist, dp_LMdist, dsdt_dist, &
              &                  dwdt_dist, dpdt_dist, tscheme, 1, .true., .false.)
      else
         if ( l_heat ) then
            call get_entropy_rhs_imp(s_LMdist, ds_LMdist, dsdt_dist, 1, .true.)
         end if
         call get_pol_rhs_imp(s_LMdist, xi_LMdist, w_LMdist, dw_LMdist, ddw_LMdist,  &
              &               p_LMdist, dp_LMdist, dwdt_dist, dpdt_dist, tscheme, 1, .true.,&
              &               .false., .false., work_LMdist)
      end if
      call get_rot_rates(omega_ma, lorentz_torque_ma_dt%old(1))
      call get_rot_rates(omega_ic, lorentz_torque_ic_dt%old(1))
      call get_tor_rhs_imp(time, z_LMdist, dz_LMdist, dzdt_dist, domega_ma_dt, &
           &               domega_ic_dt, omega_ic, omega_ma, omega_ic1,        &
           &               omega_ma1, tscheme, 1, .true., .false.)

      if ( l_mag .or. l_mag_kin  ) then
         call get_mag_rhs_imp(b_LMdist, db_LMdist, ddb_LMdist, aj_LMdist,   &
              &               dj_LMdist, ddj_LMdist, dbdt_dist, djdt_dist,  &
              &               tscheme, 1, .true., .false.)
      end if
      if ( l_cond_ic ) then
         call get_mag_ic_rhs_imp(b_ic_LMdist, db_ic_LMdist, ddb_ic_LMdist,    &
              &                  aj_ic_LMdist, dj_ic_LMdist, ddj_ic_LMdist,   &
              &                  dbdt_ic_dist, djdt_ic_dist, 1, .true.)
      end if

      !~~~~~~~~~~~~~~~~~~~~~~~ Conversion Loc > Dist ~~~~~~~~~~~~~~~~~~~~~~
      call transform_new2old(w_LMdist, w_LMloc, n_r_max)
      call transform_new2old(dw_LMdist, dw_LMloc, n_r_max)
      call transform_new2old(ddw_LMdist, ddw_LMloc, n_r_max)
      call transform_new2old(p_LMdist, p_LMloc, n_r_max)
      call transform_new2old(dp_LMdist, dp_LMloc, n_r_max)
      call dwdt_dist%gather_all(dwdt)
      call dpdt_dist%gather_all(dpdt)
      call transform_new2old(z_LMdist, z_LMloc, n_r_max)
      call transform_new2old(dz_LMdist, dz_LMloc, n_r_max)
      call dzdt_dist%gather_all(dzdt)

      if ( l_heat ) then
         call transform_new2old(s_LMdist, s_LMloc, n_r_max)
         call transform_new2old(ds_LMdist, ds_LMloc, n_r_max)
         call dsdt_dist%gather_all(dsdt)
      end if
      if ( l_chemical_conv ) then
         call transform_new2old(xi_LMdist, xi_LMloc, n_r_max)
         call transform_new2old(dxi_LMdist, dxi_LMloc, n_r_max)
         call dxidt_dist%gather_all(dxidt)
      end if
      if ( l_mag ) then
         call transform_new2old(b_LMdist, b_LMloc, n_r_max)
         call transform_new2old(db_LMdist, db_LMloc, n_r_max)
         call transform_new2old(ddb_LMdist, ddb_LMloc, n_r_max)
         call transform_new2old(aj_LMdist, aj_LMloc, n_r_max)
         call transform_new2old(dj_LMdist, dj_LMloc, n_r_max)
         call transform_new2old(ddj_LMdist, ddj_LMloc, n_r_max)
         call dbdt_dist%gather_all(dbdt)
         call djdt_dist%gather_all(djdt)
         if ( l_cond_ic ) then
            call transform_new2old(b_ic_LMdist, b_ic_LMloc, n_r_ic_max)
            call transform_new2old(db_ic_LMdist, db_ic_LMloc, n_r_ic_max)
            call transform_new2old(ddb_ic_LMdist, ddb_ic_LMloc, n_r_ic_max)
            call transform_new2old(aj_ic_LMdist, aj_ic_LMloc, n_r_ic_max)
            call transform_new2old(dj_ic_LMdist, dj_ic_LMloc, n_r_ic_max)
            call transform_new2old(ddj_ic_LMdist, ddj_ic_LMloc, n_r_ic_max)
            call dbdt_ic_dist%gather_all(dbdt_ic)
            call djdt_ic_dist%gather_all(djdt_ic)
         end if
      end if 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      !--- Get symmetry properties of tops excluding l=m=0:
      sES=0.0_cp
      sEA=0.0_cp
      sAA=0.0_cp
      if ( .not. l_axi ) then
         do m=0,l_max,minc
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
      else
         do l=0,l_max
            if ( l > 0 ) then
               if ( mod(l,2) == 0 ) then
                  sES=sES+cc2real(tops(l,0),0)
               else
                  sEA=sEA+cc2real(tops(l,0),0)
               end if
            end if
         end do
      end if
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
         if ( .not. l_axi ) then
            do m=0,l_max,minc
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
         else
            do l=0,l_max
               if ( l > 0 ) then
                  if ( mod(l,2) == 0 ) then
                     xiES=xiES+cc2real(topxi(l,0),0)
                  else
                     xiEA=xiEA+cc2real(topxi(l,0),0)
                  end if
               end if
            end do
         end if
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

   end subroutine getStartFields
!------------------------------------------------------------------------------
end module start_fields
