#include "perflib_preproc.cpp"
module start_fields

#ifdef WITH_MPI
   use mpimod
#endif
   use truncation
   use precision_mod
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: dr_fac_ic, chebt_ic, rscheme_oc,          &
       &                       chebt_ic_even, r, or1, alpha0, dLtemp0,   &
       &                       dLalpha0, beta, orho1, temp0, rho0,       &
       &                       otemp1, ogrun
   use physical_parameters, only: interior_model, epsS, impS, n_r_LCR,   &
       &                          ktopv, kbotv, LFfac, imagcon, ThExpNb, &
       &                          ViscHeatFac, impXi
   use num_param, only: dtMax, alpha
   use special, only: lGrenoble
   use blocking, only: lo_map
   use logic, only: l_conv, l_mag, l_cond_ic, l_heat, l_SRMA, l_SRIC,    &
       &            l_mag_kin, l_mag_LF, l_rot_ic, l_z10Mat, l_LCR,      &
       &            l_rot_ma, l_temperature_diff, l_single_matrix,       &
       &            l_chemical_conv, l_TP_form, l_anelastic_liquid
   use init_fields, only: l_start_file, init_s1, init_b1, tops, pt_cond, &
       &                  initV, initS, initB, initXi, ps_cond,          &
       &                  start_file, init_xi1, topxi, xi_cond
   use fields ! The entire module is required
   use fieldsLast ! The entire module is required
   use constants, only: zero, c_lorentz_ma, c_lorentz_ic, osq4pi, &
       &            one, two
   use useful, only: cc2real, logWrite
   use LMLoop_data, only: llm, ulm, ulmMag, llmMag
   use parallel_mod, only: rank, n_procs
   use radial_der, only: get_dr, get_ddr
   use radial_der_even, only: get_ddr_even
   use readCheckPoints, only: readStartFields_old, readStartFields

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

   subroutine getStartFields(time,dt,dtNew,n_time_step)
      !
      !  Purpose of this subroutine is to initialize the fields and
      !  other auxiliary parameters.
      !

      !---- Output variables:
      real(cp), intent(out) :: time,dt,dtNew
      integer,  intent(out) :: n_time_step

      !-- Local variables:
      integer :: nR,l1m0,l,m
      integer :: lm, n_r
      real(cp) :: coex
      real(cp) :: d_omega_ma_dt,d_omega_ic_dt
      character(len=76) :: message

      real(cp) :: sEA,sES,sAA
      real(cp) :: xiEA,xiES,xiAA

      real(cp) :: s0(n_r_max),p0(n_r_max),ds0(n_r_max),dp0(n_r_max)

      complex(cp), allocatable :: workA_LMloc(:,:),workB_LMloc(:,:)

      integer :: ierr, filehandle

      !PERFON('getFlds')
      !print*,"Starting getStartFields"
      !write(*,"(2(A,L1))") "l_conv=",l_conv,", l_heat=",l_heat
      !---- Computations for the Nusselt number if we are anelastic
      !     Can be done before setting the fields
      if ( l_heat ) then

         if ( rank == 0 ) open(newunit=filehandle, file='pscond.dat')

         if ( l_TP_form .or. l_anelastic_liquid ) then ! temperature

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

         if ( index(start_file, 'rst_') /= 0 ) then
            call readStartFields_old( w_LMloc,dwdtLast_LMloc,z_LMloc,dzdtLast_lo, &
                 &                p_LMloc,dpdtLast_LMloc,s_LMloc,dsdtLast_LMloc,  &
                 &                xi_LMloc,dxidtLast_LMloc,b_LMloc,dbdtLast_LMloc,&
                 &                aj_LMloc,djdtLast_LMloc,b_ic_LMloc,             &
                 &                dbdt_icLast_LMloc,aj_ic_LMloc,                  &
                 &                djdt_icLast_LMloc,omega_ic,omega_ma,            &
                 &                lorentz_torque_icLast,lorentz_torque_maLast,    &
                 &                time,dt,dtNew,n_time_step )
         else
            call readStartFields( w_LMloc,dwdtLast_LMloc,z_LMloc,dzdtLast_lo,     &
                 &                p_LMloc,dpdtLast_LMloc,s_LMloc,dsdtLast_LMloc,  &
                 &                xi_LMloc,dxidtLast_LMloc,b_LMloc,dbdtLast_LMloc,&
                 &                aj_LMloc,djdtLast_LMloc,b_ic_LMloc,             &
                 &                dbdt_icLast_LMloc,aj_ic_LMloc,                  &
                 &                djdt_icLast_LMloc,omega_ic,omega_ma,            &
                 &                lorentz_torque_icLast,lorentz_torque_maLast,    &
                 &                time,dt,dtNew,n_time_step )
         end if

         if ( dt > 0.0_cp ) then
            if ( rank==0 ) write(message,'(''! Using old time step:'',ES16.6)') dt
         else
            dt=dtMax
            if ( rank==0 ) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
         end if

         if ( .not. l_heat ) then
            s_LMloc(:,:)       =zero
            dsdtLast_LMloc(:,:)=zero
         end if

            !PERFOFF
      else ! If there's no restart file

         ! Initialize with zero
         if ( l_conv .or. l_mag_kin ) then
            w_LMloc(:,:)       =zero
            dwdtLast_LMloc(:,:)=zero
            z_LMloc(:,:)       =zero
            dzdtLast_lo(:,:)   =zero
            p_LMloc(:,:)       =zero
            dpdtLast_LMloc(:,:)=zero
         end if
         if ( l_heat ) then
            s_LMloc(:,:)       =zero
            dsdtLast_LMloc(:,:)=zero
         end if
         if ( l_chemical_conv ) then
            xi_LMloc(:,:)       =zero
            dxidtLast_LMloc(:,:)=zero
         end if
         if ( l_mag ) then
            b_LMloc(:,:)       =zero
            dbdtLast_LMloc(:,:)=zero
            aj_LMloc(:,:)      =zero
            djdtLast_LMloc(:,:)=zero
         end if
         if ( l_cond_ic ) then
            b_ic_LMloc(:,:)       =zero
            dbdt_icLast_LMloc(:,:)=zero
            aj_ic_LMloc(:,:)      =zero
            djdt_icLast_LMloc(:,:)=zero
         end if

         time =0.0_cp
         dt   =dtMax
         dtNew=dtMax
         n_time_step=0
         if (rank == 0) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
      end if
      call logWrite(message)

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

      allocate( workA_LMloc(llm:ulm,n_r_max) )
      allocate( workB_LMloc(llm:ulm,n_r_max) )

      !-- Initialize/add fields
      !----- Initialize/add magnetic field:
      if ( ( imagcon /= 0 .or. init_b1 /= 0 .or. lGrenoble ) &
           & .and. ( l_mag .or. l_mag_LF ) ) then
         call initB(b_LMloc,aj_LMloc,b_ic_LMloc,aj_ic_LMloc,      &
              &     lorentz_torque_icLast, lorentz_torque_maLast)
      end if

      !----- Initialize/add velocity, set IC and ma rotation:
      if ( l_conv .or. l_mag_kin .or. l_SRIC .or. l_SRMA ) then
         call initV(w_LMloc,z_LMloc,omega_ic,omega_ma)
      end if

      !----- Initialize/add entropy:
      if ( ( init_s1 /= 0 .or. impS /= 0 ) .and. l_heat ) then
         call initS(s_LMloc,p_LMloc)
      end if

      !----- Initialize/add chemical convection:
      if ( ( init_xi1 /= 0 .or. impXi /= 0 ) .and. l_chemical_conv ) then
         call initXi(xi_LMloc)
      end if

      !  Computing derivatives
      if ( l_conv .or. l_mag_kin ) then
         call get_ddr( w_LMloc,dw_LMloc,ddw_LMloc,ulm-llm+1,1, &
              &        ulm-llm+1,n_r_max,rscheme_oc )
         call get_dr( z_LMloc,dz_LMloc,ulm-llm+1, 1,ulm-llm+1, &
              &       n_r_max,rscheme_oc )
      end if

      if ( l_mag .or. l_mag_kin  ) then
         call get_ddr( b_LMloc,db_LMloc,ddb_LMloc,ulmMag-llmMag+1,  &
              &        1,ulmMag-llmMag+1,n_r_max,rscheme_oc )
         call get_ddr( aj_LMloc,dj_LMloc,ddj_LMloc,ulmMag-llmMag+1, &
              &        1,ulmMag-llmMag+1,n_r_max,rscheme_oc )
      end if
      if ( l_cond_ic ) then
         call get_ddr_even(b_ic_LMloc,db_ic_LMLoc,ddb_ic_LMloc,ulmMag-llmMag+1, &
              &            1,ulmMag-llmMag+1,n_r_ic_max,n_cheb_ic_max,          &
              &            dr_fac_ic,workA_LMloc,workB_LMloc,chebt_ic,          &
              &            chebt_ic_even)
         call get_ddr_even(aj_ic_LMloc,dj_ic_LMloc,ddj_ic_LMloc,ulmMag-llmMag+1,&
              &            1,ulmMag-llmMag+1,n_r_ic_max,n_cheb_ic_max,          &
              &            dr_fac_ic,workA_LMloc,workB_LMloc,chebt_ic,          &
              &            chebt_ic_even)
      end if

      if ( l_LCR ) then
         do nR=n_r_cmb,n_r_icb-1
            if ( nR<=n_r_LCR ) then
               do lm=llm,ulm
                  l=lo_map%lm2l(lm)
                  m=lo_map%lm2m(lm)

                  b_LMloc(lm,nR)=(r(n_r_LCR)/r(nR))**real(l,cp)* &
                  &               b_LMloc(lm,n_r_LCR)
                  db_LMloc(lm,nR)=-real(l,cp)*(r(n_r_LCR))**real(l,cp)/ &
                  &               (r(nR))**real(l+1,cp)*b_LMloc(lm,n_r_LCR)
                  ddb_LMloc(lm,nR)=real(l,cp)*real(l+1,cp)*    &
                  &                (r(n_r_LCR))**(real(l,cp))/ &
                  &                (r(nR))**real(l+2,cp)*b_LMloc(lm,n_r_LCR)
                  aj_LMloc(lm,nR)=zero
                  dj_LMloc(lm,nR)=zero
                  ddj_LMloc(lm,nR)=zero
               end do
            end if
         end do
      end if


      if ( l_heat ) then
         !-- Get radial derivatives of entropy:
         call get_dr( s_LMloc,ds_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc )
         if ( l_single_matrix ) then
            call get_dr( p_LMloc,dp_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc )
         end if
      end if

      if ( l_chemical_conv ) then
         !-- Get radial derivatives of chemical composition:
         call get_dr( xi_LMloc,dxi_LMloc,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc )
      end if

      deallocate(workA_LMloc)
      deallocate(workB_LMloc)

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

      !----- Get changes in mantle and ic rotation rate:
      if ( .not. l_mag_LF ) then
         lorentz_torque_icLast=0.0_cp
         lorentz_torque_maLast=0.0_cp
      end if
      if ( l_z10mat ) then
         l1m0=lo_map%lm2(1,0)
         coex=-two*(alpha-one)
         if ( ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ).and.&
              & (l1m0 >= llm .and.l1m0 <= ulm) ) then
            d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_maLast
            d_omega_ma_dtLast=d_omega_ma_dt -                              &
            &                 coex * ( two*or1(1)*real( z_LMloc(l1m0,1)) - &
            &                                     real(dz_LMloc(l1m0,1)) )
         end if
         if ( ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ).and.&
              & (l1m0 >= llm .and. l1m0 <= ulm) ) then
            d_omega_ic_dt=LFfac*c_lorentz_ic*lorentz_torque_icLast
            d_omega_ic_dtLast= d_omega_ic_dt +coex * (                         &
            &                  two*or1(n_r_max)*real( z_LMloc(l1m0,n_r_max)) - &
            &                                   real(dz_LMloc(l1m0,n_r_max)) )
         end if
      else
         d_omega_ma_dtLast=0.0_cp
         d_omega_ic_dtLast=0.0_cp
      end if

   end subroutine getStartFields
!------------------------------------------------------------------------------
end module start_fields
