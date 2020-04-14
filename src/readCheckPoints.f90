module readCheckPoints
   !
   ! This module contains the functions that can help reading and
   ! mapping of the restart files
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use parallel_mod
   use communications, only: scatter_from_rank0_to_lo
   use fields, only: dw_LMloc, ddw_LMloc, ds_LMloc, dp_LMloc, dz_LMloc,   &
       &             dxi_LMloc, db_LMloc, ddb_LMloc, dj_LMloc, ddj_LMloc, &
       &             db_ic_LMloc, ddb_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc
   use truncation, only: n_r_max,lm_max,n_r_maxMag,lm_maxMag,n_r_ic_max, &
       &                 n_r_ic_maxMag,nalias,n_phi_tot,l_max,m_max,     &
       &                 minc,lMagMem,fd_stretch,fd_ratio,n_r_icb,       &
       &                 n_r_cmb, load, getBlocks
   use logic, only: l_rot_ma,l_rot_ic,l_SRIC,l_SRMA,l_cond_ic,l_heat,l_mag, &
       &            l_mag_LF, l_chemical_conv, l_AB1, l_bridge_step,        &
       &            l_double_curl, l_z10Mat, l_single_matrix
   use blocking, only: lo_map, lm2l, lm2m, lm_balance, llm, ulm, llmMag, &
       &               ulmMag, st_map
   use init_fields, only: start_file,inform,tOmega_ic1,tOmega_ic2,             &
       &                  tOmega_ma1,tOmega_ma2,omega_ic1,omegaOsz_ic1,        &
       &                  omega_ic2,omegaOsz_ic2,omega_ma1,omegaOsz_ma1,       &
       &                  omega_ma2,omegaOsz_ma2,tShift_ic1,tShift_ic2,        &
       &                  tShift_ma1,tShift_ma2,tipdipole, scale_b, scale_v,   &
       &                  scale_s,scale_xi
   use radial_functions, only: rscheme_oc, chebt_ic, cheb_norm_ic, r
   use num_param, only: alph1, alph2, alpha
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio, &
       &                          kbotv, ktopv, sc, raxi, LFfac
   use constants, only: c_z10_omega_ic, c_z10_omega_ma, pi, zero, two
   use chebyshev, only: type_cheb_odd
   use radial_scheme, only: type_rscheme
   use finite_differences, only: type_fd
   use cosine_transform_odd, only: costf_odd_t
   use useful, only: polynomial_interpolation, abortRun
   use constants, only: one, c_lorentz_ma, c_lorentz_ic
   use updateWP_mod, only: get_pol_rhs_imp
   use updateZ_mod, only: get_tor_rhs_imp
   use updateS_mod, only: get_entropy_rhs_imp
   use updateXI_mod, only: get_comp_rhs_imp
   use updateB_mod, only: get_mag_rhs_imp, get_mag_ic_rhs_imp
   use updateWPS_mod, only: get_single_rhs_imp
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray, type_tscalar

   implicit none

   private

   logical :: lreadS, lreadXi, lreadR
   logical :: l_axi_old

   integer :: n_start_file
   integer(lip) :: bytes_allocated=0
   class(type_rscheme), pointer :: rscheme_oc_old
   real(cp) :: ratio1_old, ratio2_old, ratio1, ratio2

   public :: readStartFields_old, readStartFields
#ifdef WITH_MPI
   public :: readStartFields_mpi
#endif

contains

   subroutine readStartFields_old(w,dwdt,z,dzdt,p,dpdt,s,dsdt,xi,dxidt,b,     &
              &                   dbdt,aj,djdt,b_ic,dbdt_ic,aj_ic,djdt_ic,    &
              &                   omega_ic,omega_ma,domega_ic_dt,domega_ma_dt,&
              &                   lorentz_torque_ic_dt,lorentz_torque_ma_dt,  &
              &                   time,tscheme,n_time_step)
      !
      ! This subroutine is used to read the old restart files produced
      ! by MagIC. This is now deprecated with the change of the file format.
      ! This is still needed to read old files.
      !

      !-- Output:
      real(cp),            intent(out) :: time
      class(type_tscheme), intent(inout) :: tscheme
      integer,             intent(out) :: n_time_step
      real(cp),            intent(out) :: omega_ic,omega_ma
      complex(cp),         intent(out) :: w(llm:ulm,n_r_max),z(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: s(llm:ulm,n_r_max),p(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),         intent(out) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      type(type_tarray),   intent(inout) :: dwdt, dzdt, dpdt, dsdt, dxidt
      type(type_tarray),   intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ma_dt, domega_ic_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ma_dt,lorentz_torque_ic_dt

      !-- Local:
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old,n_r_ic_max_old,lm,nR,n_o
      real(cp) :: pr_old,ra_old,pm_old,raxi_old,sc_old
      real(cp) :: ek_old,radratio_old,sigma_ratio_old
      logical :: l_mag_old,startfile_does_exist
      integer :: informOld,ioerr,n_r_maxL,n_r_ic_maxL,lm_max_old
      integer, allocatable :: lm2lmo(:)

      real(cp) :: omega_ic1Old,omegaOsz_ic1Old,omega_ic2Old,omegaOsz_ic2Old
      real(cp) :: omega_ma1Old,omegaOsz_ma1Old,omega_ma2Old,omegaOsz_ma2Old

      character(len=72) :: rscheme_version_old
      real(cp) :: r_icb_old, r_cmb_old, dom_ic, dom_ma, coex
      integer :: n_in, n_in_2, l1m0

      complex(cp), allocatable :: wo(:,:),zo(:,:),po(:,:),so(:,:),xio(:,:)
      complex(cp), allocatable :: workA(:,:),workB(:,:),workC(:,:)
      complex(cp), allocatable :: workD(:,:),workE(:,:)
      real(cp), allocatable :: r_old(:), dt_array_old(:)

      if ( rscheme_oc%version == 'cheb') then
         ratio1 = alph1
         ratio2 = alph2
      else
         ratio1 = fd_stretch
         ratio2 = fd_ratio
      end if

      allocate( dt_array_old(max(2,tscheme%nexp)) )
      dt_array_old(:)=0.0_cp

      if ( coord_r ==0 ) then
         inquire(file=start_file, exist=startfile_does_exist)

         if ( startfile_does_exist ) then
            open(newunit=n_start_file, file=start_file, status='old', &
            &    form='unformatted')
         else
            call abortRun('! The restart file does not exist !')
         end if

         sigma_ratio_old=0.0_cp  ! assume non conducting inner core !
         if ( inform == -1 ) then ! This is default !
            read(n_start_file)                                         &
            &    time,dt_array_old(2),ra_old,pr_old,pm_old,ek_old,     &
            &    radratio_old,informOld,n_r_max_old,n_theta_max_old,   &
            &    n_phi_tot_old,minc_old,nalias_old,n_r_ic_max_old,     &
            &    sigma_ratio_old
            n_time_step=0
         else if ( inform == 0 ) then
            read(n_start_file)                                         &
            &    time,dt_array_old(2),ra_old,pr_old,pm_old,ek_old,     &
            &    radratio_old,n_time_step,n_r_max_old,n_theta_max_old, &
            &    n_phi_tot_old,minc_old,nalias_old
         else if ( inform == 1 ) then
            read(n_start_file)                                         &
            &    time,dt_array_old(2),ra_old,pr_old,pm_old,ek_old,     &
            &    radratio_old,n_time_step,n_r_max_old,n_theta_max_old, &
            &    n_phi_tot_old,minc_old
            nalias_old=nalias
         else if ( inform >= 2 ) then
            read(n_start_file)                                         &
            &    time,dt_array_old(2),ra_old,pr_old,pm_old,ek_old,     &
            &    radratio_old,n_time_step,n_r_max_old,n_theta_max_old, &
            &    n_phi_tot_old,minc_old,nalias_old,n_r_ic_max_old,     &
            &    sigma_ratio_old
         end if
         if ( inform == -1 ) inform=informOld

         dt_array_old(3:tscheme%nexp)=dt_array_old(2)

         !---- Compare parameters:
         if ( l_master_rank ) then
            if ( ra /= ra_old )                                                         &
            &    write(output_unit,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') &
            &    ra_old,ra
            if ( ek /= ek_old )                                                      &
            &    write(output_unit,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') &
            &    ek_old,ek
            if ( pr /= pr_old )                                                        &
            &    write(output_unit,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') &
            &    pr_old,pr
            if ( prmag /= pm_old )                                                    &
            &    write(output_unit,'(/,'' ! New mag Pr.number (old/new):'',2ES16.6)') &
            &    pm_old,prmag
            if ( radratio /= radratio_old )                                              &
            &    write(output_unit,'(/,'' ! New mag aspect ratio (old/new):'',2ES16.6)') &
            &    radratio_old,radratio
            if ( sigma_ratio /= sigma_ratio_old )                                       &
            &    write(output_unit,'(/,'' ! New mag cond. ratio (old/new):'',2ES16.6)') &
            &    sigma_ratio_old,sigma_ratio
         end if

         if ( n_phi_tot_old == 1 ) then ! Axisymmetric restart file
            l_max_old=nalias_old*n_theta_max_old/30
            l_axi_old=.true.
         else
            l_max_old=nalias_old*n_phi_tot_old/60
            l_axi_old=.false.
         end if
         l_mag_old=.false.
         if ( pm_old /= 0.0_cp ) l_mag_old= .true.

         if ( l_master_rank ) then
            if ( n_phi_tot_old /= n_phi_tot) &
            &    write(output_unit,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
            if ( nalias_old /= nalias) &
            &    write(output_unit,*) '! New nalias (old,new)   :',nalias_old,nalias
            if ( l_max_old /= l_max ) &
            &    write(output_unit,*) '! New l_max (old,new)    :',l_max_old,l_max
         end if

         if ( inform==6 .or. inform==7 .or. inform==9 .or. inform==11 .or. &
         &    inform==13 .or. inform==21 .or. inform==23) then
            lreadS=.false.
         else
            lreadS=.true.
         end if

         if ( inform==13 .or. inform==14 .or. inform==23 .or. inform==24 ) then
            lreadXi=.true.
         else
            lreadXi=.false.
         end if

         !-- Radius is now stored since it can also handle finite differences
         if ( inform > 20 ) then
            lreadR=.true.
         else
            lreadR=.false.
         end if

         allocate( r_old(n_r_max_old) )

         if ( lreadR ) then
            read(n_start_file) rscheme_version_old, n_in, n_in_2, &
            &                  ratio1_old, ratio2_old
            if ( rscheme_version_old == 'cheb' ) then
               allocate ( type_cheb_odd :: rscheme_oc_old )
            else
               allocate ( type_fd :: rscheme_oc_old )
            end if
         else
            rscheme_version_old='cheb'
            n_in  =n_r_max_old-2 ! Just a guess here
            n_in_2=0 ! Regular grid
            ratio1_old=ratio1
            ratio2_old=0.0_cp
            allocate ( type_cheb_odd :: rscheme_oc_old )
         end if


         r_icb_old=radratio_old/(one-radratio_old)
         r_cmb_old=one/(one-radratio_old)

         call rscheme_oc_old%initialize(n_r_max_old, n_in, n_in_2)

         call rscheme_oc_old%get_grid(n_r_max_old, r_icb_old, r_cmb_old, &
              &                       ratio1_old, ratio2_old, r_old)

         if ( l_master_rank ) then
            if ( rscheme_oc%version /= rscheme_oc_old%version )                       &
            &    write(output_unit,'(/,'' ! New radial scheme (old/new):'',A4,A1,A4)')&
            &    rscheme_oc_old%version,'/', rscheme_oc%version
         end if

         allocate( lm2lmo(lm_max) )

         call getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old,m_max, &
              &         minc,minc_old,lm_max,lm_max_old,lm2lmo)

         ! allocation of local arrays.
         ! if this becomes a performance bottleneck, one can make a module
         ! and allocate the array only once in the initialization
         allocate( wo(lm_max_old,n_r_max_old),zo(lm_max_old,n_r_max_old) )
         allocate( po(lm_max_old,n_r_max_old),so(lm_max_old,n_r_max_old) )
         bytes_allocated = bytes_allocated + 4*lm_max_old*n_r_max_old* &
         &                 SIZEOF_DEF_COMPLEX
         ! end of allocation

         !PERFON('mD_rd')
         if ( lreadXi ) then
            allocate(xio(lm_max_old,n_r_max_old))
            bytes_allocated = bytes_allocated + lm_max_old*n_r_max_old* &
            &                 SIZEOF_DEF_COMPLEX
            if ( lreadS ) then
               read(n_start_file) wo, zo, po, so, xio
            else
               read(n_start_file) wo, zo, po, xio
            end if
         else
            allocate(xio(1,1))
            if ( lreadS ) then
               read(n_start_file) wo, zo, po, so
            else
               read(n_start_file) wo, zo, po
            end if
         end if
         !PERFOFF

      end if ! l_master_rank

#ifdef WITH_MPI
      call MPI_Bcast(l_mag_old,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(lreadS,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(lreadXi,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(minc_old,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(inform,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(sigma_ratio_old,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(time,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(n_time_step,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(dt_array_old,size(dt_array_old),MPI_DEF_REAL,0, &
           &         comm_r,ierr)
#endif

      coex = two*(one-alpha)

      if ( coord_r == 0 ) then
         allocate( workA(lm_max,n_r_max), workB(lm_max,n_r_max) )
         allocate( workC(lm_max,n_r_max) )
         allocate( workD(lm_max,n_r_max) )
         allocate( workE(lm_max,n_r_max) )
         bytes_allocated = bytes_allocated + 5*lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      else
         allocate( workA(1,n_r_max), workB(1,n_r_max), workC(1,n_r_max) )
         allocate( workD(1,n_r_max), workE(1,n_r_max) )
      end if

      workA(:,:)=zero
      workB(:,:)=zero
      workC(:,:)=zero
      workD(:,:)=zero
      workE(:,:)=zero

      if ( coord_r == 0 ) then
         n_r_maxL = max(n_r_max,n_r_max_old)

         call mapDataHydro( wo,zo,po,so,xio,r_old,lm2lmo,n_r_max_old,  &
              &             lm_max_old,n_r_maxL,.false.,.false.,       &
              &             .false.,.false.,.false.,workA,workB,workC, &
              &             workD,workE )
      end if

      !-- Scatter everything
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(workA(:,nR),w(llm:ulm,nR))
         call scatter_from_rank0_to_lo(workB(:,nR),z(llm:ulm,nR))
         call scatter_from_rank0_to_lo(workC(:,nR),p(llm:ulm,nR))
         if ( l_heat ) then
            call scatter_from_rank0_to_lo(workD(:,nR),s(llm:ulm,nR))
         end if
         if ( l_chemical_conv ) then
            call scatter_from_rank0_to_lo(workE(:,nR),xi(llm:ulm,nR))
         end if
      end do

      if ( l_heat .and. .not. lreadS ) then ! No entropy before
         s(:,:)=zero
      end if
      if ( l_chemical_conv .and. .not. lreadXi ) then ! No composition before
         xi(:,:)=zero
      end if

      workA(:,:)=zero
      workB(:,:)=zero
      workC(:,:)=zero
      workD(:,:)=zero
      workE(:,:)=zero

      !-- Read the d?dt fields
      if ( coord_r == 0 ) then
         if ( lreadXi ) then
            if ( lreadS ) then
               read(n_start_file) so,wo,zo,po,xio
            else
               read(n_start_file) wo,zo,po,xio
            end if
         else
            if ( lreadS ) then
               read(n_start_file) so,wo,zo,po
            else
               read(n_start_file) wo,zo,po
            end if
         end if
         !PERFOFF

         call mapDataHydro( wo,zo,po,so,xio,r_old,lm2lmo,n_r_max_old,  &
              &             lm_max_old,n_r_maxL,.true.,.true.,.true.,  &
              &             .true.,.true.,workA,workB,workC,workD,workE )
      end if

      !-- Scatter everything
      if ( tscheme%family == 'MULTISTEP' .and. tscheme%nexp >=2 ) then
         do nR=1,n_r_max
            call scatter_from_rank0_to_lo(workA(:,nR),dwdt%expl(llm:ulm,nR,2))
            call scatter_from_rank0_to_lo(workB(:,nR),dzdt%expl(llm:ulm,nR,2))
            call scatter_from_rank0_to_lo(workC(:,nR),dpdt%expl(llm:ulm,nR,2))
            if ( l_heat ) then
               call scatter_from_rank0_to_lo(workD(:,nR),dsdt%expl(llm:ulm,nR,2))
            end if

            if ( l_chemical_conv ) then
               call scatter_from_rank0_to_lo(workE(:,nR),dxidt%expl(llm:ulm,nR,2))
            end if
         end do

         if ( l_single_matrix ) then
            call get_single_rhs_imp(s, ds_LMloc, w, dw_LMloc, ddw_LMloc, p,     &
                 &                  dp_LMloc, dsdt, dwdt, dpdt, tscheme, 1,     &
                 &                  .true., .false.)
         else
            call get_pol_rhs_imp(s, xi, w, dw_LMloc, ddw_LMloc, p, dp_LMloc, &
                 &               dwdt, dpdt, tscheme, 1, .true., .false.,    &
                 &               .false., z)
            !-- z is a work array in the above expression
            if ( l_heat ) call get_entropy_rhs_imp(s, ds_LMloc, dsdt, 1, .true.)
         end if
         dwdt%expl(:,:,2)=dwdt%expl(:,:,2)+coex*dwdt%impl(:,:,1)
         if ( .not. l_double_curl ) dpdt%expl(:,:,2)=dpdt%expl(:,:,2)+coex*dpdt%impl(:,:,1)
         if ( l_heat ) dsdt%expl(:,:,2)=dsdt%expl(:,:,2)+coex*dsdt%impl(:,:,1)

         call get_tor_rhs_imp(time, z, dz_LMloc, dzdt, domega_ma_dt, domega_ic_dt, &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1, tscheme, 1,&
              &               .true., .false.)
         dzdt%expl(:,:,2)=dzdt%expl(:,:,2)+coex*dzdt%impl(:,:,1)

         if ( l_chemical_conv ) then
            call get_comp_rhs_imp(xi, dxi_LMloc, dxidt, 1, .true.)
            dxidt%expl(:,:,2)=dxidt%expl(:,:,2)+coex*dxidt%impl(:,:,1)
         end if
      end if

      deallocate(workA, workB, workC, workD, workE)

      if ( coord_r == 0 ) then
         allocate( workA(lm_maxMag,n_r_maxMag), workB(lm_maxMag,n_r_maxMag) )
         allocate( workC(lm_maxMag,n_r_maxMag) )
         allocate( workD(lm_maxMag, n_r_maxMag) )
         bytes_allocated = bytes_allocated - 5*lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + 4*lm_maxMag*n_r_maxMag*SIZEOF_DEF_COMPLEX
      else
         allocate( workA(1,n_r_maxMag), workB(1,n_r_maxMag), workC(1,n_r_maxMag) )
         allocate( workD(1,n_r_maxMag) )
      end if

      workA(:,:)=zero
      workB(:,:)=zero
      workC(:,:)=zero
      workD(:,:)=zero

      if ( coord_r == 0 ) then
         if ( lreadXi ) then
            read(n_start_file) raxi_old, sc_old
            if ( raxi /= raxi_old ) &
               write(output_unit,'(/,'' ! New composition-based Rayleigh number (old/new):'',2ES16.6)') raxi_old,raxi
            if ( sc /= sc_old ) &
               write(output_unit,'(/,'' ! New Schmidt number (old/new):'',2ES16.6)') sc_old,sc
         end if

         if ( l_mag_old ) then
            read(n_start_file) so,wo,zo,po

            if ( l_mag ) then
               call mapDataMag( wo,zo,po,so,r_old,n_r_max,n_r_max_old,   &
                     &           lm_max_old,n_r_maxL,lm2lmo,n_r_maxMag,   &
                     &           .false.,workA,workB,workC,workD )
            end if
         else
            write(output_unit,*) '! No magnetic data in input file!'
         end if
      end if

      !-- Scatter everything
      if ( l_mag_old .and. l_mag ) then
         do nR=1,n_r_maxMag
            call scatter_from_rank0_to_lo(workA(:,nR),aj(llm:ulm,nR))
            call scatter_from_rank0_to_lo(workD(:,nR),b(llm:ulm,nR))
         end do
      end if
      if ( tscheme%family == 'MULTISTEP' .and. tscheme%nexp >=2 ) then
         if ( l_mag_old .and. l_mag ) then
            do nR=1,n_r_maxMag
               call scatter_from_rank0_to_lo(workB(:,nR),dbdt%expl(llm:ulm,nR,2))
               call scatter_from_rank0_to_lo(workC(:,nR),djdt%expl(llm:ulm,nR,2))
            end do

            call get_mag_rhs_imp(b, db_LMloc, ddb_LMloc, aj, dj_LMloc, ddj_LMloc, &
                 &               dbdt, djdt, tscheme, 1, .true., .false.)
            dbdt%expl(:,:,2)=dbdt%expl(:,:,2)+coex*dbdt%impl(:,:,1)
            djdt%expl(:,:,2)=djdt%expl(:,:,2)+coex*djdt%impl(:,:,1)
         end if
      end if

      deallocate( workA, workB, workC, workD )

      !-- Inner core part
      if ( l_mag_old ) then
         if ( coord_r == 0 ) then
            allocate( workA(lm_max,n_r_ic_max), workB(lm_max,n_r_ic_max) )
            allocate( workC(lm_max,n_r_ic_max), workD(lm_max,n_r_ic_max) )
            bytes_allocated = bytes_allocated - 4*lm_maxMag*n_r_maxMag* &
            &                 SIZEOF_DEF_COMPLEX
            bytes_allocated = bytes_allocated + 4*lm_max*n_r_ic_max* &
            &                 SIZEOF_DEF_COMPLEX
         else
            allocate( workA(1,n_r_ic_max), workB(1,n_r_ic_max) )
            allocate( workC(1,n_r_ic_max), workD(1,n_r_ic_max) )
         end if

         workA(:,:)=zero
         workB(:,:)=zero
         workC(:,:)=zero
         workD(:,:)=zero
      end if

      if ( coord_r == 0 ) then
         ! deallocation of the local arrays
         deallocate( wo,zo,po,so )
         bytes_allocated = bytes_allocated - 4*(lm_max_old*n_r_max_old)*&
         &                 SIZEOF_DEF_COMPLEX

         if ( lreadXi ) then
            bytes_allocated = bytes_allocated - lm_max_old*n_r_max_old* &
            &                 SIZEOF_DEF_COMPLEX
         end if
      end if


      !-- Inner core fields:
      if ( l_mag .or. l_mag_LF ) then
         if ( l_mag_old ) then

            if ( inform >= 2 .and. sigma_ratio_old /= 0.0_cp ) then
               if ( coord_r == 0 ) then
                  n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)
                  allocate( wo(lm_max_old,n_r_ic_max_old) )
                  allocate( zo(lm_max_old,n_r_ic_max_old) )
                  allocate( po(lm_max_old,n_r_ic_max_old) )
                  allocate( so(lm_max_old,n_r_ic_max_old) )

                  read(n_start_file) so,wo,zo,po
                  if ( l_mag ) then
                     call mapDataMag( wo,zo,po,so,r_old,n_r_ic_max,          &
                          &           n_r_ic_max_old,lm_max_old,n_r_ic_maxL, &
                          &           lm2lmo,n_r_ic_maxMag,.true.,workA,     &
                          &           workB,workC,workD )
                  end if

                  deallocate( wo,zo,po,so )
               end if

               if ( l_cond_ic ) then
                  do nR=1,n_r_ic_max
                     call scatter_from_rank0_to_lo(workA(:,nR),aj_ic(llm:ulm,nR))
                     call scatter_from_rank0_to_lo(workD(:,nR),b_ic(llm:ulm,nR))
                  end do
               end if

               if ( l_cond_ic .and. tscheme%family == 'MULTISTEP' .and. &
               &    tscheme%nexp >=2 ) then
                  do nR=1,n_r_ic_max
                     call scatter_from_rank0_to_lo(workB(:,nR), &
                          &                        dbdt_ic%expl(llm:ulm,nR,2))
                     call scatter_from_rank0_to_lo(workC(:,nR), &
                          &                        djdt_ic%expl(llm:ulm,nR,2))
                  end do

                  call get_mag_ic_rhs_imp(b_ic, db_ic_LMloc, ddb_ic_LMloc, aj_ic,  & 
                       &                  dj_ic_LMloc, ddj_ic_LMloc, dbdt_ic,      &
                       &                  djdt_ic, 1, .true.)
                  dbdt_ic%expl(:,:,2)=dbdt_ic%expl(:,:,2)+coex*dbdt_ic%impl(:,:,1)
                  djdt_ic%expl(:,:,2)=djdt_ic%expl(:,:,2)+coex*djdt_ic%impl(:,:,1)
               end if

            else if ( l_cond_ic ) then
               !----- No inner core fields provided by start_file, we thus assume that
               !      simple the inner core field decays like r**(l+1) from
               !      the ICB to r=0:
               if ( l_master_rank ) write(output_unit,'(/,'' ! USING POTENTIAL IC fields!'')')

               do lm=llm,ulm
                  do nR=1,n_r_ic_max
                     b_ic(lm,nR)   =b(lm,n_r_CMB)
                     aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                  end do
               end do
            end if
         end if
      end if

      if ( l_mag_old ) deallocate( workA, workB, workC, workD )

      omega_ic1Old     =0.0_cp
      omegaOsz_ic1Old  =0.0_cp
      tOmega_ic1       =0.0_cp
      omega_ic2Old     =0.0_cp
      omegaOsz_ic2Old  =0.0_cp
      tOmega_ic2       =0.0_cp
      omega_ma1Old     =0.0_cp
      omegaOsz_ma1Old  =0.0_cp
      tOmega_ma1       =0.0_cp
      omega_ma2Old     =0.0_cp
      omegaOsz_ma2Old  =0.0_cp
      tOmega_ma2       =0.0_cp
      dt_array_old(1)  =dt_array_old(2)

      if ( coord_r == 0 ) then
         deallocate( r_old, lm2lmo )
         call rscheme_oc_old%finalize() ! deallocate old radial scheme

         !-- Lorentz-torques:
         !   NOTE: If lMagMem=.false. the memory required to read
         !         magnetic field is not available. The code therefore
         !         cannot read lorentz torques and rotations that
         !         are stored after the magnetic fields.
         !         In this case I set the lorentz torques to zero and
         !         calculate the rotation from the speed at the
         !         boundaries in the case of no slip conditions.
         if ( inform == 3 .and. l_mag_old .and. lMagMem == 1 ) then
            read(n_start_file,iostat=ioerr) dom_ic, dom_ma
            if( ioerr/=0 ) then
               write(output_unit,*) '! Could not read last line in input file!'
               write(output_unit,*) '! Data missing or wrong format!'
               write(output_unit,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else if ( inform >= 4 .and. inform <= 6 .and. lMagMem == 1 )then
            read(n_start_file,iostat=ioerr) dom_ic,dom_ma,omega_ic,omega_ma
            omega_ic1Old=omega_ic
            omega_ic2Old=0.0_cp
            omega_ma1Old=omega_ma
            omega_ma2Old=0.0_cp
            if( ioerr/=0 .and. l_master_rank ) then
               write(output_unit,*) '! Could not read last line in input file!'
               write(output_unit,*) '! Data missing or wrong format!'
               write(output_unit,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else if ( inform == 7 .or. inform == 8 ) then
            read(n_start_file,iostat=ioerr) dom_ic, dom_ma,    &
            &    omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1,      &
            &    omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2,      &
            &    omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1,      &
            &    omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2
            if( ioerr/=0 ) then
               write(output_unit,*) '! Could not read last line in input file!'
               write(output_unit,*) '! Data missing or wrong format!'
               write(output_unit,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else if ( inform > 8 ) then
            read(n_start_file,iostat=ioerr) dom_ic, dom_ma,    &
            &    omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1,      &
            &    omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2,      &
            &    omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1,      &
            &    omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2,      &
            &    dt_array_old(1)
            if( ioerr/=0 ) then
               write(output_unit,*) '! Could not read last line in input file!'
               write(output_unit,*) '! Data missing or wrong format!'
               write(output_unit,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else
            !-- These could possibly be calcualted from the B-field
            dom_ic=0.0_cp
            dom_ma=0.0_cp
         end if
         if ( inform < 11 ) then
            dom_ic=pm_old*dom_ic
            dom_ma=pm_old*dom_ma
         end if

         if ( l_master_rank ) then
            if ( l_SRIC ) then
               if ( omega_ic1Old /= omega_ic1 )                                 &
               &    write(output_unit,*) '! New IC rotation rate 1 (old/new):', &
               &    omega_ic1Old,omega_ic1
               if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                                &
               &    write(output_unit,*) '! New IC rotation osz. rate 1 (old/new):', &
               &    omegaOsz_ic1Old,omegaOsz_ic1
               if ( omega_ic2Old /= omega_ic2 )                                 &
               &    write(output_unit,*) '! New IC rotation rate 2 (old/new):', &
               &    omega_ic2Old,omega_ic2
               if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                                &
               &    write(output_unit,*) '! New IC rotation osz. rate 2 (old/new):', &
               &    omegaOsz_ic2Old,omegaOsz_ic2
            end if
            if ( l_SRMA ) then
               if ( omega_ma1Old /= omega_ma1 )                                 &
               &    write(output_unit,*) '! New MA rotation rate 1 (old/new):', &
               &    omega_ma1Old,omega_ma1
               if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                                &
               &    write(output_unit,*) '! New MA rotation osz. rate 1 (old/new):', &
               &    omegaOsz_ma1Old,omegaOsz_ma1
               if ( omega_ma2Old /= omega_ma2 )                                 &
               &    write(output_unit,*) '! New MA rotation rate 2 (old/new):', &
               &    omega_ma2Old,omega_ma2
               if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                                &
               &    write(output_unit,*) '! New MA rotation osz. rate 2 (old/new):', &
               &    omegaOsz_ma2Old,omegaOsz_ma2
            end if
         end if ! l_master_rank
      end if ! coord_r == 0

#ifdef WITH_MPI
      call MPI_Bcast(omega_ic1Old,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(omega_ma1Old,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ic1,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ic2,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ma1,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ma2,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(dom_ic,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(dom_ma,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(dt_array_old,size(dt_array_old),MPI_DEF_REAL,0, &
           &         comm_r,ierr)
#endif

      !-- Fill the time step array
      do n_o=1,size(tscheme%dt)
         !-- If the new scheme has higher order one fill the missing dt values
         !-- with the oldest
         if ( n_o > size(dt_array_old) ) then
            tscheme%dt(n_o)=dt_array_old(size(dt_array_old))
         else
            tscheme%dt(n_o)=dt_array_old(n_o)
         end if
      end do

      !-- If old and new schemes differ in precision, one has to use a bridging step
      if ( tscheme%family == 'MULTISTEP' ) then
         if ( tscheme%nold > 1 .or. tscheme%nimp > 1 ) then
            l_bridge_step = .true.
         else
            l_bridge_step = .false.
         end if
      else if ( tscheme%family == 'DIRK' ) then
         l_bridge_step = .false.
      end if

      if (coord_r == 0) close(n_start_file)

      deallocate ( dt_array_old )

      !-- Put the torques correctly in the time step array
      if ( tscheme%family == 'MULTISTEP' .and. tscheme%nexp >=2 ) then
         lorentz_torque_ic_dt%expl(2) = dom_ic
         lorentz_torque_ma_dt%expl(2) = dom_ma
      end if

      !-- Finish computation to restart
      call finish_start_fields(time, minc_old, l_mag_old, omega_ic1Old, &
           &                   omega_ma1Old, z, s, xi, b, omega_ic, omega_ma)

      !----- Get changes in mantle and ic rotation rate:
      if ( tscheme%family == 'MULTISTEP' .and. tscheme%nexp >=2 ) then
         if ( .not. l_mag_LF ) then
            lorentz_torque_ic_dt%expl(2)=0.0_cp
            lorentz_torque_ma_dt%expl(2)=0.0_cp
         end if
         if ( l_z10mat ) then
            l1m0=lo_map%lm2(1,0)
            if ( ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ).and.&
            &     (l1m0 >= llm .and.l1m0 <= ulm) ) then
               domega_ma_dt%expl(2)=LFfac*c_lorentz_ma*dom_ma
            end if
            if ( ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ).and.&
                 & (l1m0 >= llm .and. l1m0 <= ulm) ) then
               domega_ic_dt%expl(2)=LFfac*c_lorentz_ic*dom_ic
            end if
         else
            domega_ma_dt%expl(2)=0.0_cp
            domega_ic_dt%expl(2)=0.0_cp
         end if
      end if

   end subroutine readStartFields_old
!------------------------------------------------------------------------------
   subroutine readStartFields(w,dwdt,z,dzdt,p,dpdt,s,dsdt,xi,dxidt,b,dbdt, &
              &               aj,djdt,b_ic,dbdt_ic,aj_ic,djdt_ic,omega_ic, &
              &               omega_ma,domega_ic_dt,domega_ma_dt,          &
              &               lorentz_torque_ic_dt,lorentz_torque_ma_dt,   &
              &               time,tscheme,n_time_step)
      !
      ! This subroutine is used to read the restart files produced
      ! by MagIC.
      !

      !-- Output:
      !-- Output:
      real(cp),            intent(out) :: time
      class(type_tscheme), intent(inout) :: tscheme
      integer,             intent(out) :: n_time_step
      real(cp),            intent(out) :: omega_ic,omega_ma
      complex(cp),         intent(out) :: w(llm:ulm,n_r_max),z(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: s(llm:ulm,n_r_max),p(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),         intent(out) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      type(type_tarray),   intent(inout) :: dwdt, dzdt, dpdt, dsdt, dxidt
      type(type_tarray),   intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ma_dt, domega_ic_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ma_dt,lorentz_torque_ic_dt

      !-- Local:
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old,n_r_ic_max_old, io_status,lm,nR
      real(cp) :: pr_old,ra_old,pm_old,raxi_old,sc_old
      real(cp) :: ek_old,radratio_old,sigma_ratio_old,coex
      logical :: l_mag_old, l_heat_old, l_cond_ic_old, l_chemical_conv_old
      logical :: startfile_does_exist
      integer :: nimp_old, nexp_old, nold_old
      logical :: l_press_store_old
      integer :: n_r_maxL,n_r_ic_maxL,lm_max_old, n_o
      integer, allocatable :: lm2lmo(:)
      real(cp) :: dom_ic, dom_ma
      real(cp) :: omega_ic1Old,omegaOsz_ic1Old,omega_ic2Old,omegaOsz_ic2Old
      real(cp) :: omega_ma1Old,omegaOsz_ma1Old,omega_ma2Old,omegaOsz_ma2Old

      character(len=10) :: tscheme_family_old
      character(len=72) :: rscheme_version_old
      real(cp) :: r_icb_old, r_cmb_old
      integer :: n_in, n_in_2, version, l1m0

      complex(cp), allocatable :: workOld(:,:), work(:,:)
      real(cp), allocatable :: r_old(:), dt_array_old(:)

      if ( rscheme_oc%version == 'cheb') then
         ratio1 = alph1
         ratio2 = alph2
      else
         ratio1 = fd_stretch
         ratio2 = fd_ratio
      end if

      if ( coord_r ==0 ) then
         inquire(file=start_file, exist=startfile_does_exist)

         if ( startfile_does_exist ) then
            !-- First try without stream
            open(newunit=n_start_file, file=start_file, status='old', &
            &    form='unformatted')

            read(n_start_file, iostat=io_status) version

            if ( io_status /= 0 .or. abs(version) > 100 ) then
               write(output_unit,*) '! The checkpoint file does not have record markers'
               write(output_unit,*) '! I try to read it with a stream access...'
               close(n_start_file)

               !-- Second attempt without stream
               open(newunit=n_start_file, file=start_file, status='old', &
               &    form='unformatted', access='stream')

               read(n_start_file, iostat=io_status) version
               if ( io_status /= 0 ) then
                  call abortRun('! The restart file has a wrong formatting !')
               end if
            end if
         else
            call abortRun('! The restart file does not exist !')
         end if

         if ( index(start_file, 'checkpoint_ave') /= 0 .and. l_master_rank) then
            write(output_unit,*) '! This is a time-averaged checkpoint'
            write(output_unit,*) '! d#dt arrays will be set to zero'
            l_AB1=.true.
         end if

         if ( version == 1 ) then ! This was CN/AB2 in the initial version
            allocate( dt_array_old(max(2,tscheme%nexp)) )
            dt_array_old(:)=0.0_cp
            nimp_old = 1
            nold_old = 1
            nexp_old = 2
            tscheme_family_old = 'MULTISTEP'
            read(n_start_file) time, dt_array_old(2), n_time_step
            dt_array_old(nexp_old+1:)=dt_array_old(nexp_old)
         else
            read(n_start_file) time
            read(n_start_file) tscheme_family_old
            read(n_start_file) nexp_old
            read(n_start_file) nimp_old
            read(n_start_file) nold_old

            if ( tscheme_family_old == 'MULTISTEP' ) then
               allocate( dt_array_old(max(nexp_old,tscheme%nexp) ) )
               dt_array_old(:)=0.0_cp
               read(n_start_file) dt_array_old(1:nexp_old)
               dt_array_old(nexp_old+1:)=dt_array_old(nexp_old)
            else if ( tscheme_family_old == 'DIRK' ) then
               allocate( dt_array_old(max(1,size(tscheme%dt))) )
               dt_array_old(:)=0.0_cp
               read(n_start_file) dt_array_old(1)
               dt_array_old(2:size(tscheme%dt))=dt_array_old(1)
            end if

            read(n_start_file) n_time_step
         end if


         read(n_start_file) ra_old,pr_old,raxi_old,sc_old,pm_old, &
         &                  ek_old,radratio_old,sigma_ratio_old
         read(n_start_file) n_r_max_old,n_theta_max_old,n_phi_tot_old,&
         &                  minc_old,nalias_old,n_r_ic_max_old

         if ( n_phi_tot_old == 1 ) then ! Axisymmetric restart file
            l_max_old=nalias_old*n_theta_max_old/30
            l_axi_old=.true.
         else
            l_max_old=nalias_old*n_phi_tot_old/60
            l_axi_old=.false.
         end if

         !---- Compare parameters:
         call print_info(ra_old,ek_old,pr_old,sc_old,raxi_old,pm_old, &
              &          radratio_old,sigma_ratio_old,n_phi_tot_old,  &
              &          nalias_old, l_max_old)

         allocate( r_old(n_r_max_old) )

         read(n_start_file) rscheme_version_old, n_in, n_in_2, &
         &                  ratio1_old, ratio2_old
         if ( rscheme_version_old == 'cheb' ) then
            allocate ( type_cheb_odd :: rscheme_oc_old )
         else
            allocate ( type_fd :: rscheme_oc_old )
         end if

         r_icb_old=radratio_old/(one-radratio_old)
         r_cmb_old=one/(one-radratio_old)

         call rscheme_oc_old%initialize(n_r_max_old, n_in, n_in_2)

         if ( version == 1 ) then
            call rscheme_oc_old%get_grid(n_r_max_old, r_icb_old, r_cmb_old, &
                 &                       ratio1_old, ratio2_old, r_old)
         else
            read(n_start_file) r_old(:)
         end if

         if ( rscheme_oc%version /= rscheme_oc_old%version )             &
         &    write(output_unit,'(/,'' ! New radial scheme (old/new):'',A4,A1,A4)')&
         &    rscheme_oc_old%version,'/', rscheme_oc%version

         !-- Determine the old mapping
         allocate( lm2lmo(lm_max) )
         call getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old,m_max,minc, &
              &         minc_old,lm_max,lm_max_old,lm2lmo)
         n_r_maxL = max(n_r_max,n_r_max_old)

         !-- Read Lorentz torques and rotation rates:
         if ( version == 1 ) then

            read(n_start_file) dom_ic, dom_ma, omega_ic1Old,            &
            &                  omegaOsz_ic1Old,tOmega_ic1,              &
            &                  omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
            &                  omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
            &                  omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2, &
            &                  dt_array_old(1)
            if ( tscheme%nexp >=2 .and. tscheme%family=='MULTISTEP' ) then
               lorentz_torque_ic_dt%expl(2)=dom_ic
               lorentz_torque_ma_dt%expl(2)=dom_ma
            end if

         else ! New version

            call read_map_one_scalar(n_start_file, tscheme, nexp_old, nimp_old, &
                 &                   nold_old, tscheme_family_old, domega_ic_dt)
            call read_map_one_scalar(n_start_file, tscheme, nexp_old, nimp_old, &
                 &                   nold_old, tscheme_family_old, domega_ma_dt)
            call read_map_one_scalar(n_start_file, tscheme, nexp_old,       &
                 &                   nimp_old, nold_old, tscheme_family_old,&
                 &                   lorentz_torque_ic_dt)
            call read_map_one_scalar(n_start_file, tscheme, nexp_old,       &
                 &                   nimp_old, nold_old, tscheme_family_old,&
                 &                   lorentz_torque_ma_dt)

            read(n_start_file) omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
            &                  omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
            &                  omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
            &                  omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2

         end if

         if ( l_SRIC ) then
            if ( omega_ic1Old /= omega_ic1 )                       &
            &    write(output_unit,*) '! New IC rotation rate 1 (old/new):', &
            &    omega_ic1Old,omega_ic1
            if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                      &
            &    write(output_unit,*) '! New IC rotation osz. rate 1 (old/new):', &
            &    omegaOsz_ic1Old,omegaOsz_ic1
            if ( omega_ic2Old /= omega_ic2 )                       &
            &    write(output_unit,*) '! New IC rotation rate 2 (old/new):', &
            &    omega_ic2Old,omega_ic2
            if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                      &
            &    write(output_unit,*) '! New IC rotation osz. rate 2 (old/new):', &
            &    omegaOsz_ic2Old,omegaOsz_ic2
         end if
         if ( l_SRMA ) then
            if ( omega_ma1Old /= omega_ma1 )                       &
            &    write(output_unit,*) '! New MA rotation rate 1 (old/new):', &
            &    omega_ma1Old,omega_ma1
            if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                      &
            &    write(output_unit,*) '! New MA rotation osz. rate 1 (old/new):', &
            &    omegaOsz_ma1Old,omegaOsz_ma1
            if ( omega_ma2Old /= omega_ma2 )                       &
            &    write(output_unit,*) '! New MA rotation rate 2 (old/new):', &
            &    omega_ma2Old,omega_ma2
            if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                      &
            &    write(output_unit,*) '! New MA rotation osz. rate 2 (old/new):', &
            &    omegaOsz_ma2Old,omegaOsz_ma2
         end if

         if ( version == 1 ) then
            read(n_start_file) l_heat_old, l_chemical_conv_old, l_mag_old, &
            &                  l_cond_ic_old
            l_press_store_old = .true.
         else
            read(n_start_file) l_heat_old, l_chemical_conv_old, l_mag_old, &
            &                  l_press_store_old, l_cond_ic_old
         end if

      end if ! coord_r == 0

#ifdef WITH_MPI
      call MPI_Bcast(version,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(nexp_old,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(nold_old,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(nimp_old,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(tscheme_family_old,len(tscheme_family_old),MPI_CHARACTER,0, &
           &         comm_r,ierr)
      if ( tscheme_family_old == 'MULTISTEP' ) then
         if ( coord_r /= 0 ) allocate( dt_array_old(max(nexp_old,tscheme%nexp)) )
      else if ( tscheme_family_old == 'DIRK' ) then
         if ( coord_r /= 0 ) allocate( dt_array_old(max(1,size(tscheme%dt))) )
      end if
      call MPI_Bcast(dt_array_old,size(dt_array_old),MPI_DEF_REAL,0, &
           &         comm_r,ierr)
      call MPI_Bcast(l_mag_old,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(l_heat_old,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(l_press_store_old,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(l_chemical_conv_old,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(l_cond_ic_old,1,MPI_LOGICAL,0,comm_r,ierr)
      call MPI_Bcast(minc_old,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(time,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(omega_ic1Old,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(omega_ma1Old,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ic1,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ic2,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ma1,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(tOmega_ma2,1,MPI_DEF_REAL,0,comm_r,ierr)
      call MPI_Bcast(n_time_step,1,MPI_INTEGER,0,comm_r,ierr)
      call MPI_Bcast(domega_ic_dt%expl, tscheme%nexp, MPI_DEF_REAL, 0, &
           &         comm_r, ierr)
      call MPI_Bcast(domega_ic_dt%impl, tscheme%nimp, MPI_DEF_REAL, 0, &
           &         comm_r, ierr)
      call MPI_Bcast(domega_ic_dt%old, tscheme%nold, MPI_DEF_REAL, 0, &
           &         comm_r, ierr)
      call MPI_Bcast(domega_ma_dt%expl, tscheme%nexp, MPI_DEF_REAL, 0, &
           &         comm_r, ierr)
      call MPI_Bcast(domega_ma_dt%impl, tscheme%nimp, MPI_DEF_REAL, 0, &
           &         comm_r, ierr)
      call MPI_Bcast(domega_ma_dt%old, tscheme%nold, MPI_DEF_REAL, 0, &
           &         comm_r, ierr)
      call MPI_Bcast(lorentz_torque_ic_dt%expl, tscheme%nexp, MPI_DEF_REAL, &
           &         0, comm_r, ierr)
      call MPI_Bcast(lorentz_torque_ic_dt%impl, tscheme%nimp, &
           &         MPI_DEF_REAL, 0, comm_r, ierr)
      call MPI_Bcast(lorentz_torque_ic_dt%old, tscheme%nold, MPI_DEF_REAL, &
           &         0, comm_r, ierr)
      call MPI_Bcast(lorentz_torque_ma_dt%expl, tscheme%nexp, MPI_DEF_REAL, &
           &         0, comm_r, ierr)
      call MPI_Bcast(lorentz_torque_ma_dt%impl, tscheme%nimp, &
           &         MPI_DEF_REAL, 0, comm_r, ierr)
      call MPI_Bcast(lorentz_torque_ma_dt%old, tscheme%nold, MPI_DEF_REAL, &
           &         0, comm_r, ierr)
#endif

      !-- Fill the time step array
      do n_o=1,size(tscheme%dt)
         !-- If the new scheme has higher order one fill the missing dt values
         !-- with the oldest
         if ( n_o > size(dt_array_old) ) then
            tscheme%dt(n_o)=dt_array_old(size(dt_array_old))
         else
            tscheme%dt(n_o)=dt_array_old(n_o)
         end if
      end do

      !-- If old and new schemes differ in precision, one has to use a bridging step
      if ( tscheme%family == 'MULTISTEP' ) then
         if ( tscheme_family_old == 'DIRK' ) then
            l_bridge_step = .true.
         else
            if ( tscheme%nold > nold_old .or. tscheme%nimp > nimp_old ) then
               l_bridge_step = .true.
            else
               l_bridge_step = .false.
            end if
         end if
      else if ( tscheme%family == 'DIRK' ) then
         l_bridge_step = .false.
      end if

      !-- Allocate arrays
      if ( coord_r==0 ) then
         allocate( work(lm_max,n_r_max), workOld(lm_max_old,n_r_max_old) )
      else
         allocate( work(1,n_r_max))
         allocate(workOld(1,1))
         allocate(r_old(1))
         allocate(lm2lmo(1) )
      end if

      !-- Read the poloidal flow
      call read_map_one_field( n_start_file, tscheme, workOld, work, scale_v,    &
           &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,    &
           &                   nexp_old, nimp_old, nold_old, tscheme_family_old, &
           &                   w, dwdt, .true.)

      !-- Read the toroidal flow
      call read_map_one_field( n_start_file, tscheme, workOld, work, scale_v,    &
           &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,    &
           &                   nexp_old, nimp_old, nold_old, tscheme_family_old, &
           &                   z, dzdt, .true.)

      !-- Read the pressure
      if ( l_press_store_old ) then
         call read_map_one_field( n_start_file, tscheme, workOld, work, scale_v,   &
              &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,   &
              &                   nexp_old, nimp_old, nold_old, tscheme_family_old,&
              &                   p, dpdt, .not. l_double_curl)
      end if

      !-- Read the entropy
      if ( l_heat_old ) then
         call read_map_one_field( n_start_file, tscheme, workOld, work, scale_s,    &
              &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,    &
              &                   nexp_old, nimp_old, nold_old,  tscheme_family_old,&
              &                   s, dsdt, l_heat)
      end if

      !-- Read the chemical composition
      if ( l_chemical_conv_old ) then
         call read_map_one_field( n_start_file, tscheme, workOld, work, scale_xi,  &
              &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,   &
              &                   nexp_old, nimp_old, nold_old, tscheme_family_old,&
              &                   xi, dxidt, l_chemical_conv)
      end if
      if ( l_heat .and. .not. l_heat_old ) s(:,:)=zero
      if ( l_chemical_conv .and. .not. l_chemical_conv_old ) xi(:,:)=zero
      if ( .not. l_double_curl .and. .not. l_press_store_old ) p(:,:)=zero

      if ( (l_mag .or. l_mag_LF) .and.  l_mag_old ) then

         !-- Read the poloidal magnetic field
         call read_map_one_field( n_start_file, tscheme, workOld, work, scale_b,    &
              &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,    &
              &                   nexp_old, nimp_old, nold_old,  tscheme_family_old,&
              &                   b, dbdt, .true. )

         !-- Read the toroidal magnetic field
         call read_map_one_field( n_start_file, tscheme, workOld, work, scale_b,    &
              &                   r_old, lm2lmo, n_r_max_old, n_r_maxL, n_r_max,    &
              &                   nexp_old, nimp_old, nold_old,  tscheme_family_old,&
              &                   aj, djdt, .true. )

         if ( l_cond_ic ) then

            if ( l_cond_ic_old ) then
               deallocate( work, workOld )

               if ( coord_r==0 ) then
                  n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)
                  allocate( work(lm_max,n_r_ic_max) )
                  allocate( workOld(lm_max_old,n_r_ic_max_old) )
               else
                  allocate( work(1,n_r_ic_max), workOld(1,1) )
               end if

               !-- Read the inner core poloidal magnetic field
               if ( coord_r==0 ) then
                  work(:,:)=zero
                  read(n_start_file) workOld
                  call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                       &            n_r_ic_max_old,n_r_ic_maxL,       &
                       &            n_r_ic_max,.false.,.true.,work )
                  !-- Cancel the spherically-symmetric part
                  work(1,:)=zero
               end if
               do nR=1,n_r_ic_max
                  call scatter_from_rank0_to_lo(work(:,nR),b_ic(llm:ulm,nR))
               end do

               !-- Read dbdt_ic
               if ( tscheme_family_old == 'MULTISTEP' ) then
                  do n_o=2,nexp_old
                     if ( coord_r==0 ) then
                        work(:,:)=zero
                        read(n_start_file) workOld
                        call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                             &            n_r_ic_max_old,n_r_ic_maxL,      &
                             &            n_r_ic_max,.true.,.true.,work )
                        !-- Cancel the spherically-symmetric part
                        work(1,:)=zero
                     end if
                     if ( n_o <= tscheme%nexp .and. tscheme%family=='MULTISTEP' ) then
                        do nR=1,n_r_ic_max
                           call scatter_from_rank0_to_lo(work(:,nR),  &
                                &                dbdt_ic%expl(llm:ulm,nR,n_o))
                        end do
                     end if
                  end do
                  do n_o=2,nimp_old
                     if ( coord_r==0 ) then
                        work(:,:)=zero
                        read(n_start_file) workOld
                        call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                             &            n_r_ic_max_old,n_r_ic_maxL,      &
                             &            n_r_ic_max,.true.,.true.,work )
                        !-- Cancel the spherically-symmetric part
                        work(1,:)=zero
                     end if
                     if ( n_o <= tscheme%nimp .and. tscheme%family=='MULTISTEP' ) then
                        do nR=1,n_r_ic_max
                           call scatter_from_rank0_to_lo(work(:,nR),  &
                                &                dbdt_ic%impl(llm:ulm,nR,n_o))
                        end do
                     end if
                  end do
                  do n_o=2,nold_old
                     if ( coord_r==0 ) then
                        work(:,:)=zero
                        read(n_start_file) workOld
                        call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                             &            n_r_ic_max_old,n_r_ic_maxL,      &
                             &            n_r_ic_max,.true.,.true.,work )
                        !-- Cancel the spherically-symmetric part
                        work(1,:)=zero
                     end if
                     if ( n_o <= tscheme%nold .and. &
                     &   tscheme%family=='MULTISTEP' ) then
                        do nR=1,n_r_ic_max
                           call scatter_from_rank0_to_lo(work(:,nR),  &
                                &                dbdt_ic%old(llm:ulm,nR,n_o))
                        end do
                     end if
                  end do
               end if

               !-- Read the inner core toroidal magnetic field
               if ( coord_r==0 ) then
                  work(:,:)=zero
                  read(n_start_file) workOld
                  call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                       &            n_r_ic_max_old,n_r_ic_maxL,       &
                       &            n_r_ic_max,.false.,.true.,work )
                  !-- Cancel the spherically-symmetric part
                  work(1,:)=zero
               end if
               do nR=1,n_r_ic_max
                  call scatter_from_rank0_to_lo(work(:,nR),aj_ic(llm:ulm,nR))
               end do

               !-- Read djdt_ic
               if ( tscheme_family_old == 'MULTISTEP' ) then
                  do n_o=2,nexp_old
                     if ( coord_r==0 ) then
                        work(:,:)=zero
                        read(n_start_file) workOld
                        call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                             &            n_r_ic_max_old,n_r_ic_maxL,      &
                             &            n_r_ic_max,.true.,.true.,work )
                        !-- Cancel the spherically-symmetric part
                        work(1,:)=zero
                     end if
                     if ( n_o <= tscheme%nexp  .and. tscheme%family=='MULTISTEP' ) then
                        do nR=1,n_r_ic_max
                           call scatter_from_rank0_to_lo(work(:,nR),  &
                                &                djdt_ic%expl(llm:ulm,nR,n_o))
                        end do
                     end if
                  end do
                  do n_o=2,nimp_old
                     if ( coord_r==0 ) then
                        work(:,:)=zero
                        read(n_start_file) workOld
                        call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                             &            n_r_ic_max_old,n_r_ic_maxL,      &
                             &            n_r_ic_max,.true.,.true.,work )
                        !-- Cancel the spherically-symmetric part
                        work(1,:)=zero
                     end if
                     if ( n_o <= tscheme%nimp .and. tscheme%family=='MULTISTEP' ) then
                        do nR=1,n_r_ic_max
                           call scatter_from_rank0_to_lo(work(:,nR),  &
                                &                djdt_ic%impl(llm:ulm,nR,n_o))
                        end do
                     end if
                  end do
                  do n_o=2,nold_old
                     if ( coord_r==0 ) then
                        work(:,:)=zero
                        read(n_start_file) workOld
                        call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                             &            n_r_ic_max_old,n_r_ic_maxL,      &
                             &            n_r_ic_max,.true.,.true.,work )
                        !-- Cancel the spherically-symmetric part
                        work(1,:)=zero
                     end if
                     if ( n_o <= tscheme%nold .and. &
                     &   tscheme%family=='MULTISTEP' ) then
                        do nR=1,n_r_ic_max
                           call scatter_from_rank0_to_lo(work(:,nR),  &
                                &                djdt_ic%old(llm:ulm,nR,n_o))
                        end do
                     end if
                  end do
               end if

            else
               !-- No inner core fields provided by start_file, we thus assume that
               !   simple the inner core field decays like r**(l+1) from
               !   the ICB to r=0:
               if ( l_master_rank ) write(output_unit,'(/,'' ! USING POTENTIAL IC fields!'')')

               do lm=llm,ulm
                  do nR=1,n_r_ic_max
                     b_ic(lm,nR) =b(lm,n_r_CMB)
                     aj_ic(lm,nR)=aj(lm,n_r_CMB)
                     !dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                     !djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
                  end do
               end do

            end if

         end if

      end if

      !-- Free memory
      deallocate( work, workOld, dt_array_old, r_old, lm2lmo )

      !-- Close file
      if ( coord_r==0 ) then
         call rscheme_oc_old%finalize() ! deallocate old radial scheme
         close(n_start_file)
      end if ! coord_r==0

      !-- Finish computation to restart
      call finish_start_fields(time, minc_old, l_mag_old, omega_ic1Old, &
           &                   omega_ma1Old, z, s, xi, b, omega_ic, omega_ma)

      !-- Correct explicit arrays if old version with CNAB2 was stored
      !-- Back then the d?dtLast arrays did not carry the same meaning
      if ( tscheme%family == 'MULTISTEP' .and. tscheme%nexp >= 2 .and. &
      &    version == 1 ) then
         coex = two*(one-alpha)

         if ( l_single_matrix ) then
            call get_single_rhs_imp(s, ds_LMloc, w, dw_LMloc, ddw_LMloc, p,     &
                 &                  dp_LMloc, dsdt, dwdt, dpdt, tscheme, 1,     &
                 &                  .true., .false.)
         else
            call get_pol_rhs_imp(s, xi, w, dw_LMloc, ddw_LMloc, p, dp_LMloc, &
                 &               dwdt, dpdt, tscheme, 1, .true., .false.,    &
                 &               .false., z)
            !-- z is a work array in the above expression
            if ( l_heat ) call get_entropy_rhs_imp(s, ds_LMloc, dsdt, 1, .true.)
         end if

         dwdt%expl(:,:,2)=dwdt%expl(:,:,2)+coex*dwdt%impl(:,:,1)
         if ( .not. l_double_curl ) then
            dpdt%expl(:,:,2)=dpdt%expl(:,:,2)+coex*dpdt%impl(:,:,1)
         end if
         if ( l_heat) dsdt%expl(:,:,2)=dsdt%expl(:,:,2)+coex*dsdt%impl(:,:,1)

         call get_tor_rhs_imp(time, z, dz_LMloc, dzdt, domega_ma_dt, domega_ic_dt, &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1, tscheme, 1,&
              &               .true., .false.)
         dzdt%expl(:,:,2)=dzdt%expl(:,:,2)+coex*dzdt%impl(:,:,1)


         if ( l_chemical_conv ) then
            call get_comp_rhs_imp(xi, dxi_LMloc, dxidt, 1, .true.)
            dxidt%expl(:,:,2)=dxidt%expl(:,:,2)+coex*dxidt%impl(:,:,1)
         end if

         if ( l_mag ) then
            call get_mag_rhs_imp(b, db_LMloc, ddb_LMloc, aj, dj_LMloc, ddj_LMloc, &
                 &               dbdt, djdt, tscheme, 1, .true., .false.)
            dbdt%expl(:,:,2)=dbdt%expl(:,:,2)+coex*dbdt%impl(:,:,1)
            djdt%expl(:,:,2)=djdt%expl(:,:,2)+coex*djdt%impl(:,:,1)
         end if

         if ( l_cond_ic ) then
            call get_mag_ic_rhs_imp(b_ic, db_ic_LMloc, ddb_ic_LMloc, aj_ic,  &
                 &                  dj_ic_LMloc, ddj_ic_LMloc, dbdt_ic,      &
                 &                  djdt_ic, 1, .true.)
            dbdt_ic%expl(:,:,2)=dbdt_ic%expl(:,:,2)+coex*dbdt_ic%impl(:,:,1)
            djdt_ic%expl(:,:,2)=djdt_ic%expl(:,:,2)+coex*djdt_ic%impl(:,:,1)
         end if

         if ( .not. l_mag_LF ) then
            lorentz_torque_ic_dt%expl(2)=0.0_cp
            lorentz_torque_ma_dt%expl(2)=0.0_cp
         end if
         if ( l_z10mat ) then
            l1m0=lo_map%lm2(1,0)
            if ( ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ) .and. &
            &     (l1m0 >= llm .and.l1m0 <= ulm) ) then
               domega_ma_dt%expl(2)=LFfac*c_lorentz_ma*lorentz_torque_ma_dt%expl(2)
            end if
            if ( ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ) .and. &
            &     (l1m0 >= llm .and. l1m0 <= ulm) ) then
               domega_ic_dt%expl(2)=LFfac*c_lorentz_ic*lorentz_torque_ic_dt%expl(2)
            end if
         else
            domega_ma_dt%expl(2)=0.0_cp
            domega_ic_dt%expl(2)=0.0_cp
         end if

      end if

   end subroutine readStartFields
!------------------------------------------------------------------------------
   subroutine read_map_one_scalar(fh, tscheme, nexp_old, nimp_old, nold_old,&
              &                   tscheme_family_old, dscal_dt)

      !-- Input variables
      integer,             intent(in) :: fh, nold_old
      integer,             intent(in) :: nexp_old, nimp_old
      character(len=*),    intent(in) :: tscheme_family_old
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variable
      type(type_tscalar), intent(inout) :: dscal_dt

      !-- Local variables
      integer :: n_o
      real(cp) :: scal

      if ( tscheme_family_old == 'MULTISTEP' ) then

         do n_o=2,nexp_old
            read(fh) scal
            if ( n_o <= tscheme%nexp .and.  &
            &    tscheme%family=='MULTISTEP') dscal_dt%expl(n_o)=scal
         end do
         do n_o=2,nimp_old
            read(fh) scal
            if ( n_o <= tscheme%nimp .and. &
            &    tscheme%family=='MULTISTEP' ) dscal_dt%impl(n_o)=scal
         end do
         do n_o=2,nold_old
            read(fh) scal
            if ( n_o <= tscheme%nold .and. &
            &    tscheme%family=='MULTISTEP') dscal_dt%old(n_o)=scal
         end do

      end if

   end subroutine read_map_one_scalar
!------------------------------------------------------------------------------
   subroutine read_map_one_field( fh, tscheme, wOld, work, scale_w, r_old, lm2lmo,&
              &                   n_r_max_old,  n_r_maxL, dim1, nexp_old,   &
              &                   nimp_old, nold_old, tscheme_family_old, w, dwdt,&
              &                   l_map)

      !--- Input variables
      logical,             intent(in) :: l_map
      integer,             intent(in) :: fh, nold_old
      integer,             intent(in) :: nexp_old, nimp_old
      character(len=*),    intent(in) :: tscheme_family_old
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_r_max_old,dim1
      integer,             intent(in) :: n_r_maxL
      integer,             intent(in) :: lm2lmo(lm_max)
      real(cp),            intent(in) :: r_old(:)
      real(cp),            intent(in) :: scale_w

      !--- Output variables
      complex(cp),       intent(inout) :: wOld(:,:)
      complex(cp),       intent(inout) :: work(:,:)
      complex(cp),       intent(out) :: w(llm:ulm,dim1)
      type(type_tarray), intent(inout) :: dwdt

      !-- Local variable
      integer :: n_o, nR

      if ( l_master_rank ) then
         work(:,:)=zero
         read(fh) wOld
         call mapOneField( wOld,scale_w,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,dim1,.false.,.false.,work )
      end if
      if ( l_map ) then
         do nR=1,dim1
            call scatter_from_rank0_to_lo(work(:,nR),w(llm:ulm,nR))
         end do
      end if

      if ( tscheme_family_old == 'MULTISTEP' ) then
         do n_o=2,nexp_old
            if ( l_master_rank ) then
               work(:,:)=zero
               read(fh) wOld
               call mapOneField( wOld,scale_w,r_old,lm2lmo,n_r_max_old, &
                    &            n_r_maxL,dim1,.true.,.false.,work )
            end if
            if ( n_o <= tscheme%nexp .and. l_map .and. tscheme%family == 'MULTISTEP') then
               do nR=1,n_r_max
                  call scatter_from_rank0_to_lo(work(:,nR),dwdt%expl(llm:ulm,nR,n_o))
               end do
            end if
         end do
         do n_o=2,nimp_old
            if ( l_master_rank ) then
               work(:,:)=zero
               read(fh) wOld
               call mapOneField( wOld,scale_w,r_old,lm2lmo,n_r_max_old, &
                    &            n_r_maxL,dim1,.true.,.false.,work )
            end if
            if ( n_o <= tscheme%nimp .and. l_map .and. tscheme%family=='MULTISTEP') then
               do nR=1,n_r_max
                  call scatter_from_rank0_to_lo(work(:,nR),dwdt%impl(llm:ulm,nR,n_o))
               end do
            end if
         end do
         do n_o=2,nold_old
            if ( l_master_rank ) then
               work(:,:)=zero
               read(fh) wOld
               call mapOneField( wOld,scale_w,r_old,lm2lmo,n_r_max_old, &
                    &            n_r_maxL,dim1,.true.,.false.,work )
            end if
            if ( n_o <= tscheme%nold .and. l_map .and. &
            &    tscheme%family == 'MULTISTEP' ) then
               do nR=1,n_r_max
                  call scatter_from_rank0_to_lo(work(:,nR),dwdt%old(llm:ulm,nR,n_o))
               end do
            end if
         end do
      end if

   end subroutine read_map_one_field
!------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine readStartFields_mpi(w,dwdt,z,dzdt,p,dpdt,s,dsdt,xi,dxidt,b,   &
              &                   dbdt,aj,djdt,b_ic,dbdt_ic,aj_ic,djdt_ic,  &
              &                   omega_ic,omega_ma,domega_ic_dt,           &
              &                   domega_ma_dt,lorentz_torque_ic_dt,        &
              &                   lorentz_torque_ma_dt,time,tscheme,        &
              &                   n_time_step)   
      !
      ! This subroutine is used to read the restart files produced
      ! by MagIC using MPI-IO
      !

      !-- Output:
      real(cp),            intent(out) :: time
      class(type_tscheme), intent(inout) :: tscheme
      integer,             intent(out) :: n_time_step
      real(cp),            intent(out) :: omega_ic,omega_ma
      complex(cp),         intent(out) :: w(llm:ulm,n_r_max),z(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: s(llm:ulm,n_r_max),p(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),         intent(out) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      type(type_tarray),   intent(inout) :: dwdt, dzdt, dpdt, dsdt, dxidt
      type(type_tarray),   intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ma_dt, domega_ic_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ma_dt,lorentz_torque_ic_dt

      !-- Local:
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old,lm,nR,n_r_ic_max_old
      real(cp) :: pr_old,ra_old,pm_old,raxi_old,sc_old,coex
      real(cp) :: ek_old,radratio_old,sigma_ratio_old
      logical :: l_mag_old, l_heat_old, l_cond_ic_old, l_chemical_conv_old
      logical :: startfile_does_exist
      integer :: n_r_maxL,n_r_ic_maxL,lm_max_old
      integer, allocatable :: lm2lmo(:)

      real(cp) :: omega_ic1Old,omegaOsz_ic1Old,omega_ic2Old,omegaOsz_ic2Old
      real(cp) :: omega_ma1Old,omegaOsz_ma1Old,omega_ma2Old,omegaOsz_ma2Old

      character(len=72) :: rscheme_version_old
      character(len=10) :: tscheme_family_old
      real(cp) :: r_icb_old, r_cmb_old
      integer :: n_in, n_in_2, version, info, fh, nRStart_old, nRStop_old, n_o
      integer :: nR_per_rank_old, datatype, l1m0
      integer :: istat(MPI_STATUS_SIZE)
      integer :: nimp_old, nexp_old, nold_old
      logical :: l_press_store_old
      integer(lip) :: disp, offset, size_old

      complex(cp), allocatable :: workOld(:,:)
      complex(cp), allocatable :: work(:,:)
      real(cp), allocatable :: r_old(:), dt_array_old(:)
      type(load), allocatable :: radial_balance_old(:)

      if ( rscheme_oc%version == 'cheb') then
         ratio1 = alph1
         ratio2 = alph2
      else
         ratio1 = fd_stretch
         ratio2 = fd_ratio
      end if

      !-- Default setup of MPI-IO
      call mpiio_setup(info)

      inquire(file=start_file, exist=startfile_does_exist)

      if ( startfile_does_exist ) then
         !-- Open file
         call MPI_File_Open(comm_r, start_file, MPI_MODE_RDONLY, &
              &             info, fh, ierr)
      else
         call abortRun('! The restart file does not exist !')
      end if

      disp = 0
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE, "native", &
           &                 info, ierr)

      !-- Read the header
      call MPI_File_Read(fh, version, 1, MPI_INTEGER, istat, ierr)
      !-- Little trick if wrong-endianness is detected 
      !-- version gets crazy large, so flip back to default reader then
      if ( abs(version) > 100 ) then
         call MPI_File_close(fh, ierr)
         if ( l_master_rank ) then
            write(output_unit,*) '! I cannot read it with MPI-IO'
            write(output_unit,*) '! I try to fall back on serial reader...'
         end if
         call readStartFields(w,dwdt,z,dzdt,p,dpdt,s,dsdt,xi,dxidt,b,dbdt, &
              &               aj,djdt,b_ic,dbdt_ic,aj_ic,djdt_ic,omega_ic, &
              &               omega_ma,domega_ic_dt, domega_ma_dt,         &
              &               lorentz_torque_ic_dt, lorentz_torque_ma_dt,  &
              &               time,tscheme,n_time_step)
         return
      end if
      call MPI_File_Read(fh, time, 1, MPI_DEF_REAL, istat, ierr)
      if ( version == 1 ) then ! This was CN/AB2 in the initial version
         allocate( dt_array_old(max(2,tscheme%nexp)) )
         dt_array_old(:)=0.0_cp
         nimp_old = 1
         nold_old = 1
         nexp_old = 2
         tscheme_family_old = 'MULTISTEP'
         call MPI_File_Read(fh, dt_array_old(2), 1, MPI_DEF_REAL, istat, ierr)
         dt_array_old(nexp_old+1:)=dt_array_old(nexp_old)
      else
         call MPI_File_Read(fh, tscheme_family_old, len(tscheme%family), &
              &              MPI_CHARACTER, istat, ierr)
         call MPI_File_Read(fh, nexp_old, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Read(fh, nimp_old, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Read(fh, nold_old, 1, MPI_INTEGER, istat, ierr)
         if ( tscheme_family_old == 'MULTISTEP' ) then
            allocate( dt_array_old(max(nexp_old,tscheme%nexp) ) )
            dt_array_old(:)=0.0_cp
            call MPI_File_Read(fh, dt_array_old(1:nexp_old), nexp_old, &
                 &             MPI_DEF_REAL, istat, ierr)
            dt_array_old(nexp_old+1:)=dt_array_old(nexp_old)
         else if ( tscheme_family_old == 'DIRK' ) then
            allocate( dt_array_old(max(1,size(tscheme%dt))) )
            dt_array_old(:)=0.0_cp
            call MPI_File_Read(fh, dt_array_old(1), 1, MPI_DEF_REAL, istat, ierr)
            dt_array_old(2:size(tscheme%dt))=dt_array_old(1)
         end if
      end if
      call MPI_File_Read(fh, n_time_step, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, ra_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, pr_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, raxi_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, sc_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, pm_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, ek_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, radratio_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, sigma_ratio_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, n_r_max_old, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, n_theta_max_old, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, n_phi_tot_old, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, minc_old, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, nalias_old, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, n_r_ic_max_old, 1, MPI_INTEGER, istat, ierr)

      if ( n_phi_tot_old == 1 ) then ! Axisymmetric restart file
         l_max_old=nalias_old*n_theta_max_old/30
         l_axi_old=.true.
      else
         l_max_old=nalias_old*n_phi_tot_old/60
         l_axi_old=.false.
      end if

      !---- Compare parameters:
      if ( l_master_rank ) then
         call print_info(ra_old,ek_old,pr_old,sc_old,raxi_old,pm_old, &
              &          radratio_old,sigma_ratio_old,n_phi_tot_old,  &
              &          nalias_old, l_max_old)
      end if

      allocate( r_old(n_r_max_old) )
      !-- Read scheme version (FD or CHEB)
      call MPI_File_Read(fh, rscheme_version_old, len(rscheme_version_old), &
           &             MPI_CHARACTER, istat, ierr)
      call MPI_File_Read(fh, n_in, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, n_in_2, 1, MPI_INTEGER, istat, ierr)
      call MPI_File_Read(fh, ratio1_old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, ratio2_old, 1, MPI_DEF_REAL, istat, ierr)

      if ( rscheme_version_old == 'cheb' ) then
         allocate ( type_cheb_odd :: rscheme_oc_old )
      else
         allocate ( type_fd :: rscheme_oc_old )
      end if

      r_icb_old=radratio_old/(one-radratio_old)
      r_cmb_old=one/(one-radratio_old)

      call rscheme_oc_old%initialize(n_r_max_old, n_in, n_in_2)

      if ( version == 1 ) then
         call rscheme_oc_old%get_grid(n_r_max_old, r_icb_old, r_cmb_old, &
              &                       ratio1_old, ratio2_old, r_old)
      else
         call MPI_File_Read(fh, r_old(:), n_r_max_old, MPI_DEF_REAL, istat, ierr)
      end if

      if ( l_master_rank .and. rscheme_oc%version /= rscheme_oc_old%version )      &
      &    write(output_unit,'(/,'' ! New radial scheme (old/new):'',A4,A1,A4)')   &
      &    rscheme_oc_old%version,'/', rscheme_oc%version

      !-- Determine the old mapping
      allocate( lm2lmo(lm_max) )
      call getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old,m_max,minc, &
           &         minc_old,lm_max,lm_max_old,lm2lmo)
      n_r_maxL = max(n_r_max,n_r_max_old)

      !-- Read Lorentz-torques and rotation rates:
      if ( version > 1 ) then
         call read_map_one_scalar_mpi(fh, tscheme, nexp_old, nimp_old, nold_old,  &
              &                       tscheme_family_old, domega_ic_dt)

         call read_map_one_scalar_mpi(fh, tscheme, nexp_old, nimp_old, nold_old,  &
              &                       tscheme_family_old, domega_ma_dt)
      end if

      call read_map_one_scalar_mpi(fh, tscheme, nexp_old, nimp_old, nold_old,  &
           &                       tscheme_family_old, lorentz_torque_ic_dt)

      call read_map_one_scalar_mpi(fh, tscheme, nexp_old, nimp_old, nold_old,  &
           &                       tscheme_family_old, lorentz_torque_ma_dt)

      call MPI_File_Read(fh, omega_ic1Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omegaOsz_ic2Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, tOmega_ic1, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omega_ic2Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omegaOsz_ic2Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, tOmega_ic2, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omega_ma1Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omegaOsz_ma1Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, tOmega_ma1, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omega_ma2Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, omegaOsz_ma2Old, 1, MPI_DEF_REAL, istat, ierr)
      call MPI_File_Read(fh, tOmega_ma2, 1, MPI_DEF_REAL, istat, ierr)
      if ( version == 1 ) then
         call MPI_File_Read(fh, dt_array_old(1), 1, MPI_DEF_REAL, istat, ierr)
      end if

      !-- Fill the time step array
      do n_o=1,size(tscheme%dt)
         !-- If the new scheme has higher order one fill the missing dt values
         !-- with the oldest
         if ( n_o > size(dt_array_old) ) then
            tscheme%dt(n_o)=dt_array_old(nexp_old)
         else
            tscheme%dt(n_o)=dt_array_old(n_o)
         end if
      end do

      !-- If old and new schemes differ in precision, one has to use a bridging step
      if ( tscheme%family == 'MULTISTEP' ) then
         if ( tscheme_family_old == 'DIRK' ) then
            l_bridge_step = .true.
         else
            if ( tscheme%nold > nold_old .or. tscheme%nimp > nimp_old ) then
               l_bridge_step = .true.
            else
               l_bridge_step = .false.
            end if
         end if
      else if ( tscheme%family == 'DIRK' ) then
         l_bridge_step = .false.
      end if

      if ( l_SRIC ) then
         if ( omega_ic1Old /= omega_ic1 )                       &
         &    write(output_unit,*) '! New IC rotation rate 1 (old/new):', &
         &    omega_ic1Old,omega_ic1
         if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                      &
         &    write(output_unit,*) '! New IC rotation osz. rate 1 (old/new):', &
         &    omegaOsz_ic1Old,omegaOsz_ic1
         if ( omega_ic2Old /= omega_ic2 )                       &
         &    write(output_unit,*) '! New IC rotation rate 2 (old/new):', &
         &    omega_ic2Old,omega_ic2
         if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                      &
         &    write(output_unit,*) '! New IC rotation osz. rate 2 (old/new):', &
         &    omegaOsz_ic2Old,omegaOsz_ic2
      end if
      if ( l_SRMA ) then
         if ( omega_ma1Old /= omega_ma1 )                       &
         &    write(output_unit,*) '! New MA rotation rate 1 (old/new):', &
         &    omega_ma1Old,omega_ma1
         if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                      &
         &    write(output_unit,*) '! New MA rotation osz. rate 1 (old/new):', &
         &    omegaOsz_ma1Old,omegaOsz_ma1
         if ( omega_ma2Old /= omega_ma2 )                       &
         &    write(output_unit,*) '! New MA rotation rate 2 (old/new):', &
         &    omega_ma2Old,omega_ma2
         if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                      &
         &    write(output_unit,*) '! New MA rotation osz. rate 2 (old/new):', &
         &    omegaOsz_ma2Old,omegaOsz_ma2
      end if

      !-- Read logical to know how many fields are stored
      call MPI_File_Read(fh, l_heat_old, 1, MPI_LOGICAL, istat, ierr)
      call MPI_File_Read(fh, l_chemical_conv_old, 1, MPI_LOGICAL, istat, ierr)
      call MPI_File_Read(fh, l_mag_old, 1, MPI_LOGICAL, istat, ierr)
      if ( version > 1 ) then
         call MPI_File_Read(fh, l_press_store_old, 1, MPI_LOGICAL, istat, ierr)
      else
         l_press_store_old = .true.
      end if
      call MPI_File_Read(fh, l_cond_ic_old, 1, MPI_LOGICAL, istat, ierr)

      !-- Measure offset 
      call MPI_File_get_position(fh, offset, ierr)
      call MPI_File_get_byte_offset(fh, offset, disp, ierr)

      !-- Allocate and determine old radial balance
      allocate( radial_balance_old(0:n_ranks_r-1) )
      call getBlocks(radial_balance_old, n_r_max_old, n_ranks_r)
      nRstart_old = radial_balance_old(coord_r)%nStart
      nRstop_old = radial_balance_old(coord_r)%nStop
      nR_per_rank_old = radial_balance_old(coord_r)%n_per_rank

      allocate( workOld(lm_max_old, nRstart_old:nRstop_old) )

      !-- Define a MPI type to handle the reading
      call MPI_Type_Vector(1,lm_max_old*nR_per_rank_old, lm_max_old*n_r_max_old,&
           &               MPI_DEF_COMPLEX, datatype, ierr)
      call MPI_Type_Commit(datatype, ierr)

      !-- Set-up the first displacement
      size_old = int(nRstart_old-1,kind=lip)*int(lm_max_old,kind=lip)* &
      &          int(SIZEOF_DEF_COMPLEX,kind=lip)
      disp = disp+size_old
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      !-- Poloidal potential: w
      call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
           &                      lm_max_old, n_r_max_old, nRstart_old,   &
           &                      nRstop_old, radial_balance_old, lm2lmo, &
           &                      r_old, n_r_maxL, n_r_max, scale_v,      &
           &                      nexp_old, nimp_old, nold_old,           &
           &                      tscheme_family_old, w, dwdt, disp, .true. )

      !-- Toroidal potential: z
      call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
           &                      lm_max_old, n_r_max_old, nRstart_old,   &
           &                      nRstop_old, radial_balance_old, lm2lmo, &
           &                      r_old, n_r_maxL, n_r_max, scale_v,      &
           &                      nexp_old, nimp_old, nold_old,           &
           &                      tscheme_family_old, z, dzdt, disp, .true. )

      !-- Pressure: p
      if ( l_press_store_old ) then
         call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
              &                      lm_max_old, n_r_max_old, nRstart_old,   &
              &                      nRstop_old, radial_balance_old, lm2lmo, &
              &                      r_old, n_r_maxL, n_r_max, scale_v,      &
              &                      nexp_old, nimp_old, nold_old,           &
              &                      tscheme_family_old, p, dpdt, disp,      &
              &                      .not. l_double_curl )
      end if

      !-- Entropy: s
      if ( l_heat_old ) then
         call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
              &                      lm_max_old, n_r_max_old, nRstart_old,   &
              &                      nRstop_old, radial_balance_old, lm2lmo, &
              &                      r_old, n_r_maxL, n_r_max, scale_s,      &
              &                      nexp_old, nimp_old, nold_old,           &
              &                      tscheme_family_old, s, dsdt, disp, l_heat )
      end if

      !-- Chemical composition: xi
      if ( l_chemical_conv_old ) then
         call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
              &                      lm_max_old, n_r_max_old, nRstart_old,   &
              &                      nRstop_old, radial_balance_old, lm2lmo, &
              &                      r_old, n_r_maxL, n_r_max, scale_xi,     &
              &                      nexp_old, nimp_old, nold_old,           &
              &                      tscheme_family_old, xi, dxidt, disp,    &
              &                      l_chemical_conv )
      end if

      if ( .not. l_double_curl .and. .not. l_press_store_old ) p(:,:)=zero
      if ( l_heat .and. (.not. l_heat_old) ) s(:,:)=zero
      if ( l_chemical_conv .and. .not. l_chemical_conv_old ) xi(:,:)=zero

      if ( (l_mag .or. l_mag_LF) .and. l_mag_old ) then
         !-- Read poloidal potential: b
         call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
              &                      lm_max_old, n_r_max_old, nRstart_old,   &
              &                      nRstop_old, radial_balance_old, lm2lmo, &
              &                      r_old, n_r_maxL, n_r_max, scale_b,      &
              &                      nexp_old, nimp_old, nold_old,           &
              &                      tscheme_family_old, b, dbdt, disp, .true. )

         !-- Read toroidal potential: aj
         call read_map_one_field_mpi(fh, info, datatype, tscheme, workOld,   &
              &                      lm_max_old, n_r_max_old, nRstart_old,   &
              &                      nRstop_old, radial_balance_old, lm2lmo, &
              &                      r_old, n_r_maxL, n_r_max, scale_b,      &
              &                      nexp_old, nimp_old, nold_old,           &
              &                      tscheme_family_old, aj, djdt, disp, .true. )
      end if

      deallocate(workOld)

      !-- Inner core: for the Inner core right now only coord_r 0 can read
      !-- and broadcast
      if ( (l_mag .or. l_mag_LF) .and. l_cond_ic ) then

         if ( l_cond_ic_old ) then
            n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)

            if ( coord_r==0 ) then
               allocate( workOld(lm_max_old, n_r_ic_max_old) )
               allocate( work(lm_max, n_r_ic_max) )
            else
               allocate( work(1,n_r_ic_max), workOld(1,1) )
            end if

            !-- Read the inner core poloidal magnetic field
            if ( coord_r==0 ) then
               work(:,:)=zero
               call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                    &             MPI_DEF_COMPLEX, istat, ierr)
               call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                    &            n_r_ic_max_old,n_r_ic_maxL,       &
                    &            n_r_ic_max,.false.,.true.,work )
               !-- Cancel the spherically-symmetric part
               work(1,:)=zero
            end if
            do nR=1,n_r_ic_max
               call scatter_from_rank0_to_lo(work(:,nR),b_ic(llm:ulm,nR))
            end do

            !-- Read dbdt_ic
            if ( tscheme_family_old == 'MULTISTEP' ) then
               do n_o=2,nexp_old
                  if ( coord_r==0 ) then
                     work(:,:)=zero
                     call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                          &             MPI_DEF_COMPLEX, istat, ierr)
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  if ( n_o <= tscheme%nexp .and. tscheme%family=='MULTISTEP' ) then
                     do nR=1,n_r_ic_max
                        call scatter_from_rank0_to_lo(work(:,nR), &
                             &                        dbdt_ic%expl(llm:ulm,nR,n_o))
                     end do
                  end if
               end do

               do n_o=2,nimp_old
                  if ( coord_r==0 ) then
                     work(:,:)=zero
                     call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                          &             MPI_DEF_COMPLEX, istat, ierr)
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  if ( n_o <= tscheme%nimp .and. tscheme%family=='MULTISTEP' ) then
                     do nR=1,n_r_ic_max
                        call scatter_from_rank0_to_lo(work(:,nR), &
                             &                        dbdt_ic%impl(llm:ulm,nR,n_o))
                     end do
                  end if
               end do

               do n_o=2,nold_old
                  if ( coord_r==0 ) then
                     work(:,:)=zero
                     call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                          &             MPI_DEF_COMPLEX, istat, ierr)
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  if ( n_o <= tscheme%nold .and. &
                  &    tscheme%family=='MULTISTEP' ) then
                     do nR=1,n_r_ic_max
                        call scatter_from_rank0_to_lo(work(:,nR), &
                             &                        dbdt_ic%old(llm:ulm,nR,n_o))
                     end do
                  end if
               end do
            end if

            !-- Read the inner core toroidal magnetic field
            if ( coord_r==0 ) then
               work(:,:)=zero
               call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                    &             MPI_DEF_COMPLEX, istat, ierr)
               call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                    &            n_r_ic_max_old,n_r_ic_maxL,       &
                    &            n_r_ic_max,.false.,.true.,work )
               !-- Cancel the spherically-symmetric part
               work(1,:)=zero
            end if
            do nR=1,n_r_ic_max
               call scatter_from_rank0_to_lo(work(:,nR),aj_ic(llm:ulm,nR))
            end do

            !-- Read djdt_ic
            if ( tscheme_family_old == 'MULTISTEP' ) then
               do n_o=2,nexp_old
                  if ( coord_r==0 ) then
                     work(:,:)=zero
                     call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                          &             MPI_DEF_COMPLEX, istat, ierr)
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  if ( n_o <= tscheme%nexp .and. tscheme%family=='MULTISTEP' ) then
                     do nR=1,n_r_ic_max
                        call scatter_from_rank0_to_lo(work(:,nR), &
                             &                        djdt_ic%expl(llm:ulm,nR,n_o))
                     end do
                  end if
               end do

               do n_o=2,nimp_old
                  if ( coord_r==0 ) then
                     work(:,:)=zero
                     call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                          &             MPI_DEF_COMPLEX, istat, ierr)
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  if ( n_o <= tscheme%nimp .and. tscheme%family=='MULTISTEP' ) then
                     do nR=1,n_r_ic_max
                        call scatter_from_rank0_to_lo(work(:,nR), &
                             &                        djdt_ic%impl(llm:ulm,nR,n_o))
                     end do
                  end if
               end do

               do n_o=2,nold_old
                  if ( coord_r==0 ) then
                     work(:,:)=zero
                     call MPI_File_Read(fh, workOld, lm_max_old*n_r_ic_max_old, &
                          &             MPI_DEF_COMPLEX, istat, ierr)
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  if ( n_o <= tscheme%nold .and.  &
                  &    tscheme%family=='MULTISTEP' ) then
                     do nR=1,n_r_ic_max
                        call scatter_from_rank0_to_lo(work(:,nR), &
                             &                        djdt_ic%old(llm:ulm,nR,n_o))
                     end do
                  end if
               end do

            end if ! only if multistep is the old one

            deallocate( workOld, work )

         else
            !-- No inner core fields provided by start_file, we thus assume that
            !   simple the inner core field decays like r**(l+1) from
            !   the ICB to r=0:
            if ( l_master_rank ) write(output_unit,'(/,'' ! USING POTENTIAL IC fields!'')')

            do lm=llm,ulm
               do nR=1,n_r_ic_max
                  b_ic(lm,nR)   =b(lm,n_r_CMB)
                  aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                  !dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                  !djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
               end do
            end do
         end if

      end if

      call MPI_Type_Free(datatype, ierr)
      call MPI_File_close(fh, ierr)

      !-- Deallocate work_arrays
      deallocate(radial_balance_old, r_old, dt_array_old)
      call rscheme_oc_old%finalize()

      !-- Correct explicit arrays if old version with CNAB2 was stored
      !-- Back then the d?dtLast arrays did not carry the same meaning
      if ( tscheme%family == 'MULTISTEP' .and. tscheme%nexp >= 2 .and. &
      &    version == 1 ) then
         coex = two*(one-alpha)
         if ( l_single_matrix ) then
            call get_single_rhs_imp(s, ds_LMloc, w, dw_LMloc, ddw_LMloc, p,     &
                 &                  dp_LMloc, dsdt, dwdt, dpdt, tscheme, 1,     &
                 &                  .true., .false.)
         else
            call get_pol_rhs_imp(s, xi, w, dw_LMloc, ddw_LMloc, p, dp_LMloc, &
                 &               dwdt, dpdt, tscheme, 1, .true., .false.,    &
                 &               .false., z)
            !-- z is a work array in the above expression
            if ( l_heat ) call get_entropy_rhs_imp(s, ds_LMloc, dsdt, 1, .true.)
         end if
         dwdt%expl(:,:,2)=dwdt%expl(:,:,2)+coex*dwdt%impl(:,:,1)
         if ( .not. l_double_curl ) dpdt%expl(:,:,2)=dpdt%expl(:,:,2)+coex*dpdt%impl(:,:,1)
         if ( l_heat ) dsdt%expl(:,:,2)=dsdt%expl(:,:,2)+coex*dsdt%impl(:,:,1)

         call get_tor_rhs_imp(time, z, dz_LMloc, dzdt, domega_ma_dt, domega_ic_dt, &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1, tscheme, 1,&
              &               .true., .false.)
         dzdt%expl(:,:,2)=dzdt%expl(:,:,2)+coex*dzdt%impl(:,:,1)

         if ( l_chemical_conv ) then
            call get_comp_rhs_imp(xi, dxi_LMloc, dxidt, 1, .true.)
            dxidt%expl(:,:,2)=dxidt%expl(:,:,2)+coex*dxidt%impl(:,:,1)
         end if

         if ( l_mag ) then
            call get_mag_rhs_imp(b, db_LMloc, ddb_LMloc, aj, dj_LMloc, ddj_LMloc, &
                 &               dbdt, djdt, tscheme, 1, .true., .false.)
            dbdt%expl(:,:,2)=dbdt%expl(:,:,2)+coex*dbdt%impl(:,:,1)
            djdt%expl(:,:,2)=djdt%expl(:,:,2)+coex*djdt%impl(:,:,1)
         end if

         if ( l_cond_ic ) then
            call get_mag_ic_rhs_imp(b_ic, db_ic_LMloc, ddb_ic_LMloc, aj_ic,  &
                 &                  dj_ic_LMloc, ddj_ic_LMloc, dbdt_ic,      &
                 &                  djdt_ic, 1, .true.)
            dbdt_ic%expl(:,:,2)=dbdt_ic%expl(:,:,2)+coex*dbdt_ic%impl(:,:,1)
            djdt_ic%expl(:,:,2)=djdt_ic%expl(:,:,2)+coex*djdt_ic%impl(:,:,1)
         end if

         if ( .not. l_mag_LF ) then
            lorentz_torque_ic_dt%expl(2)=0.0_cp
            lorentz_torque_ma_dt%expl(2)=0.0_cp
         end if
         if ( l_z10mat ) then
            l1m0=lo_map%lm2(1,0)
            if ( ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ) .and. &
            &     (l1m0 >= llm .and.l1m0 <= ulm) ) then
               domega_ma_dt%expl(2)=LFfac*c_lorentz_ma*lorentz_torque_ma_dt%expl(2)
            end if
            if ( ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ) .and. &
            &      (l1m0 >= llm .and. l1m0 <= ulm) ) then
               domega_ic_dt%expl(2)=LFfac*c_lorentz_ic*lorentz_torque_ic_dt%expl(2)
            end if
         else
            domega_ma_dt%expl(2)=0.0_cp
            domega_ic_dt%expl(2)=0.0_cp
         end if

      end if

      !-- Finish computation to restart
      call finish_start_fields(time, minc_old, l_mag_old, omega_ic1Old, &
           &                   omega_ma1Old, z, s, xi, b, omega_ic, omega_ma)

   end subroutine readStartFields_mpi
!------------------------------------------------------------------------------
   subroutine read_map_one_scalar_mpi(fh, tscheme, nexp_old, nimp_old, nold_old, &
              &                       tscheme_family_old, dscal_dt)

      !-- Input variables
      integer,             intent(in) :: fh, nold_old
      integer,             intent(in) :: nexp_old, nimp_old
      character(len=*),    intent(in) :: tscheme_family_old
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variable
      type(type_tscalar), intent(inout) :: dscal_dt

      !-- Local variables
      integer :: n_o
      integer :: istat(MPI_STATUS_SIZE)
      real(cp) :: scal_exp(nexp_old-1), scal_old(nold_old-1), scal_imp(nimp_old-1)

      if ( tscheme_family_old == 'MULTISTEP' ) then

         if ( nexp_old >= 2 ) then
            call MPI_File_Read(fh, scal_exp, nexp_old-1, MPI_DEF_REAL, istat, ierr)
            do n_o=2,nexp_old
               if ( n_o <= tscheme%nexp .and.  &
               &    tscheme%family=='MULTISTEP') dscal_dt%expl(n_o)=scal_exp(n_o-1)
            end do
         end if
         if ( nimp_old >= 2 ) then
            call MPI_File_Read(fh, scal_imp, nimp_old-1, MPI_DEF_REAL, istat, ierr)
            do n_o=2,nimp_old
               if ( n_o <= tscheme%nimp .and. &
               &    tscheme%family=='MULTISTEP' ) dscal_dt%impl(n_o)=scal_imp(n_o-1)
            end do
         end if
         if(  nold_old  >= 2 ) then
            call MPI_File_Read(fh, scal_old, nold_old-1, MPI_DEF_REAL, istat, ierr)
            do n_o=2,nold_old
               if ( n_o <= tscheme%nold .and. &
               &    tscheme%family=='MULTISTEP') dscal_dt%old(n_o)=scal_old(n_o-1)
            end do
         end if

      end if

   end subroutine read_map_one_scalar_mpi
!------------------------------------------------------------------------------
   subroutine read_map_one_field_mpi(fh, info, datatype, tscheme, wOld,        &
              &                      lm_max_old, n_r_max_old, nRstart_old,     &
              &                      nRstop_old, radial_balance_old, lm2lmo,   &
              &                      r_old, n_r_maxL, dim1, scale_w, nexp_old, &
              &                      nimp_old, nold_old, tscheme_family_old,   &
              &                      w, dwdt, disp, l_map )

      !--- Input variables
      logical,             intent(in) :: l_map
      integer,             intent(in) :: nexp_old, nimp_old
      integer,             intent(in) :: nold_old
      character(len=*),    intent(in) :: tscheme_family_old
      integer,             intent(in) :: fh, info, datatype
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_r_max_old, lm_max_old, nRstart_old
      integer,             intent(in) :: nRstop_old, dim1, n_r_maxL
      real(cp),            intent(in) :: r_old(:)
      integer,             intent(in) :: lm2lmo(lm_max)
      type(load),          intent(in) :: radial_balance_old(0:n_ranks_r-1)
      complex(cp),         intent(in) :: wOld(lm_max_old,nRstart_old:nRstop_old)
      real(cp),            intent(in) :: scale_w

      !--- Output variables
      integer(lip),      intent(inout) :: disp
      complex(cp),       intent(out) :: w(llm:ulm,dim1)
      type(type_tarray), intent(inout) :: dwdt

      !-- Local variables:
      integer(lip) :: size_old
      integer :: istat(MPI_STATUS_SIZE)
      integer :: n_o, nR_per_rank_old

      nR_per_rank_old = nRstop_old-nRstart_old+1

      !-- Poloidal potential: w
      call MPI_File_Read_All(fh, wOld, lm_max_old*nR_per_rank_old, &
           &                 MPI_DEF_COMPLEX, istat, ierr)
      size_old = int(n_r_max_old,kind=lip)*int(lm_max_old,kind=lip)* &
      &          int(SIZEOF_DEF_COMPLEX,kind=lip)
      disp = disp+size_old
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      if ( l_map ) then
         call mapOneField_mpi( wOld, lm_max_old, n_r_max_old, nRstart_old, &
              &                nRstop_old, radial_balance_old, lm2lmo,     &
              &                r_old, n_r_maxL, n_r_max, .false., .false., &
              &                scale_w, w )
      end if

      !-- dwdt
      if ( tscheme_family_old == 'MULTISTEP' ) then
         do n_o=2,nexp_old
            call MPI_File_Read_All(fh, wOld, lm_max_old*nR_per_rank_old, &
                 &                 MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+size_old
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
                 &                 info, ierr)
            if ( n_o <= tscheme%nexp .and. l_map .and. tscheme%family=='MULTISTEP' ) then
               call mapOneField_mpi( wOld, lm_max_old, n_r_max_old, nRstart_old, &
                    &                nRstop_old, radial_balance_old, lm2lmo,     &
                    &                r_old, n_r_maxL, n_r_max, .true., .false.,  &
                    &                scale_w, dwdt%expl(:,:,n_o) )
            end if
         end do
         do n_o=2,nimp_old
            call MPI_File_Read_All(fh, wOld, lm_max_old*nR_per_rank_old, &
                 &                 MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+size_old
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
                 &                 info, ierr)
            if ( n_o <= tscheme%nimp .and. l_map .and. tscheme%family=='MULTISTEP' ) then
               call mapOneField_mpi( wOld, lm_max_old, n_r_max_old, nRstart_old, &
                    &                nRstop_old, radial_balance_old, lm2lmo,     &
                    &                r_old, n_r_maxL, n_r_max, .true., .false.,  &
                    &                scale_v, dwdt%impl(:,:,n_o) )
            end if
         end do
         do n_o=2,nold_old
            call MPI_File_Read_All(fh, wOld, lm_max_old*nR_per_rank_old, &
                 &                 MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+size_old
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
                 &                 info, ierr)
            if ( n_o <= tscheme%nold .and. l_map .and. & 
            &    tscheme%family=='MULTISTEP' ) then
               call mapOneField_mpi( wOld, lm_max_old, n_r_max_old, nRstart_old, &
                    &                nRstop_old, radial_balance_old, lm2lmo,     &
                    &                r_old, n_r_maxL, n_r_max, .true., .false.,  &
                    &                scale_v, dwdt%old(:,:,n_o) )
            end if
         end do
      end if

   end subroutine read_map_one_field_mpi
#endif
!------------------------------------------------------------------------------
   subroutine getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old,m_max,minc, &
              &         minc_old,lm_max,lm_max_old,lm2lmo)

      !--- Input variables
      integer, intent(in) :: n_r_max,l_max,m_max,minc
      integer, intent(in) :: n_r_max_old,l_max_old,minc_old
      integer, intent(in) :: lm_max

      !--- Output variables
      integer,intent(out) :: lm2lmo(lm_max)
      integer,intent(out) :: lm_max_old

      !--- Local variables
      integer :: m_max_old,l,m,lm,lmo,lo,mo

      if ( .not. l_axi_old ) then
         m_max_old=(l_max_old/minc_old)*minc_old
      else
         m_max_old=0
      end if

      if ( l_max==l_max_old .and. minc==minc_old .and. n_r_max==n_r_max_old &
      &    .and. m_max==m_max_old ) then

         !----- Direct reading of fields, grid not changed:
         if ( l_master_rank ) write(output_unit,'(/,'' ! Reading fields directly.'')')
         lm_max_old=lm_max
      else

         !----- Mapping onto new grid !
         if ( l_master_rank ) write(output_unit,'(/,'' ! Mapping onto new grid.'')')

         if ( mod(minc_old,minc) /= 0 .and. l_master_rank)                &
         &     write(output_unit,'('' ! Warning: Incompatible old/new minc= '',2i3)')

         lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
         &          m_max_old*(m_max_old-minc_old)/(2*minc_old) +     &
         &          l_max_old-m_max_old+1

         !-- Write info to STdoUT:
         if ( l_master_rank ) then
            write(output_unit,'('' ! Old/New  l_max= '',2I4,''  m_max= '',2I4, &
            &            ''  minc= '',2I3,''  lm_max= '',2I8/)')               &
            &                l_max_old,l_max,m_max_old,m_max,                  &
            &                minc_old,minc,lm_max_old,lm_max
         end if
         if ( n_r_max_old /= n_r_max .and. l_master_rank )                   &
         &   write(output_unit,'('' ! Old/New n_r_max='',2i4)') n_r_max_old,n_r_max

      end if

      if ( .not. l_axi_old ) then
         do lm=1,lm_max
            l=lm2l(lm)
            m=lm2m(lm)
            lm2lmo(lm)=-1 ! -1 means that there is no data in the startfile
            lmo=0
            do mo=0,l_max_old,minc_old
               do lo=mo,l_max_old
                  lmo=lmo+1
                  if ( lo==l .and. mo==m ) then
                     lm2lmo(lm)=lmo ! data found in startfile
                     cycle
                  end if
               end do
            end do
         end do
      else
         do lm=1,lm_max
            l=lm2l(lm)
            m=lm2m(lm)
            lm2lmo(lm)=-1 ! -1 means that there is no data in the startfile
            lmo=0
            do lo=0,l_max_old
               lmo=lmo+1
               if ( lo==l .and. m==0 ) then
                  lm2lmo(lm)=lmo ! data found in startfile
                  cycle
               end if
            end do
         end do
      end if

   end subroutine getLm2lmO
!------------------------------------------------------------------------------
   subroutine mapOneField_mpi( wOld, lm_max_old, n_r_max_old, nRstart_old, &
              &                nRstop_old, radial_balance_old, lm2lmo,     &
              &                r_old, n_r_maxL, dim1, lBc1, l_IC, scale_w, w )

      !--- Input variables
      integer,     intent(in) :: n_r_max_old
      integer,     intent(in) :: lm_max_old
      integer,     intent(in) :: nRstart_old
      integer,     intent(in) :: nRstop_old
      integer,     intent(in) :: dim1
      integer,     intent(in) :: n_r_maxL
      logical,     intent(in) :: lBc1,l_IC
      real(cp),    intent(in) :: r_old(:)
      integer,     intent(in) :: lm2lmo(lm_max)
      type(load),  intent(in) :: radial_balance_old(0:n_ranks_r-1)
      complex(cp), intent(in) :: wOld(lm_max_old,nRstart_old:nRstop_old)
      real(cp),    intent(in) :: scale_w

      !--- Output variables
      complex(cp), intent(out) :: w(llm:ulm,dim1)

      !--- Local variables
      complex(cp) :: wtmp_Rloc(lm_max,nRstart_old:nRstop_old)
      complex(cp) :: wtmp_LMloc(llm:ulm,n_r_max_old)
      complex(cp), allocatable :: sbuff(:), rbuff(:)
      integer :: rcounts(0:n_ranks_r-1), scounts(0:n_ranks_r-1)
      integer :: rdisp(0:n_ranks_r-1), sdisp(0:n_ranks_r-1)
      integer :: n_r, lm, lmo, p, nlm_per_rank, my_lm_counts, ii
      integer :: max_send, max_recv, l, m, lm_st
      complex(cp) :: woR(n_r_maxL)
      !integer :: lm,lmo,lmStart,lmStop,n_proc


      !-- First step: remap the lm keeping the old radii distributed
      !$omp parallel do default(shared) private(n_r,lm,lmo)
      do n_r=nRstart_old, nRstop_old
         do lm=1,lm_max
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               wtmp_Rloc(lm,n_r)=scale_w*wOld(lmo,n_r)
            else
               wtmp_Rloc(lm,n_r)=zero
            end if
         end do
      end do
      !$omp end parallel do

      !-- Second step: perform a MPI transpose from
      !--  (lm_max,nRstart_old:nRstop_old) to (llm:ulm,n_r_max_old)
      do p=0,n_ranks_r-1
         my_lm_counts = lm_balance(p)%n_per_rank
         nlm_per_rank = ulm-llm+1
         scounts(p)=(nRstop_old-nRstart_old+1)*my_lm_counts
         rcounts(p)=radial_balance_old(p)%n_per_rank*nlm_per_rank
      end do

      rdisp(0)=0
      sdisp(0)=0
      do p=1,n_ranks_r-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      max_send = sum(scounts)
      max_recv = sum(rcounts)

      allocate( sbuff(1:max_send), rbuff(1:max_recv) )

      !$omp barrier
      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_r,lm,l,m,lm_st)
      do p = 0, n_ranks_r-1
         ii = sdisp(p)+1
         do n_r=nRstart_old,nRstop_old
            do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
               l = lo_map%lm2l(lm)
               m = lo_map%lm2m(lm)
               lm_st = st_map%lm2(l,m)
               sbuff(ii)=wtmp_Rloc(lm_st,n_r)
               ii = ii +1
            end do
         end do
      end do
      !$omp end parallel do

#ifdef WITH_MPI
      call MPI_Alltoallv(sbuff, scounts, sdisp, MPI_DEF_COMPLEX, &
           &             rbuff, rcounts, rdisp, MPI_DEF_COMPLEX, &
           &             comm_r, ierr)
#endif


      !$omp barrier
      !$omp parallel do default(shared) &
      !$omp private(p,ii,n_r,lm)
      do p = 0, n_ranks_r-1
         ii = rdisp(p)+1
         do n_r=radial_balance_old(p)%nStart,radial_balance_old(p)%nStop
            do lm=llm,ulm
               wtmp_LMloc(lm,n_r)=rbuff(ii)
               ii=ii+1
            end do
         end do
      end do
      !$omp end parallel do

      deallocate( sbuff, rbuff )

      !-- Third step: map the radius now that all the old radii are
      !-- in processor while the lm are already distributed
      !$omp parallel do default(shared) private(lm,woR)
      do lm=llm,ulm
         if ( dim1 /= n_r_max_old .or. ratio1 /= ratio1_old .or.        &
         &    ratio2 /= ratio2_old .or.                                 &
         &    rscheme_oc%order_boundary /= rscheme_oc_old%order_boundary&
         &    .or. rscheme_oc%version /= rscheme_oc_old%version ) then

            woR(1:n_r_max_old)=wtmp_LMloc(lm,:)
            call mapDataR(woR,r_old,dim1,n_r_max_old,n_r_maxL,lBc1,l_IC)
            w(lm,:)=woR(1:dim1)
         else
            w(lm,:)=wtmp_LMloc(lm,:)
         end if
      end do
      !$omp end parallel do

   end subroutine mapOneField_mpi
!------------------------------------------------------------------------------
   subroutine mapOneField( wo,scale_w,r_old,lm2lmo,n_r_max_old,n_r_maxL,dim1,&
              &            lBc1,l_IC,w )

      !--- Input variables
      integer,     intent(in) :: n_r_max_old,dim1
      integer,     intent(in) :: n_r_maxL
      logical,     intent(in) :: lBc1,l_IC
      integer,     intent(in) :: lm2lmo(lm_max)
      complex(cp), intent(in) :: wo(:,:)
      real(cp),    intent(in) :: r_old(:)
      real(cp),    intent(in) :: scale_w

      !--- Output variables
      complex(cp), intent(out) :: w(lm_max,dim1)

      !--- Local variables
      integer :: lm,lmo,lmStart,lmStop,n_proc
      complex(cp) :: woR(n_r_maxL)

      !$omp parallel do default(shared) private(n_proc,lmStart,lmStop,lm,lmo,woR)
      do n_proc=0,n_ranks_r-1 ! Blocking of loop over all (l,m)
         lmStart=lm_balance(n_proc)%nStart
         lmStop =lm_balance(n_proc)%nStop

         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( dim1 /= n_r_max_old .or. ratio1 /= ratio1_old .or.        &
               &    ratio2 /= ratio2_old .or.                                 &
               &    rscheme_oc%order_boundary /= rscheme_oc_old%order_boundary&
               &    .or. rscheme_oc%version /= rscheme_oc_old%version ) then

                  woR(1:n_r_max_old)=wo(lmo,:)
                  call mapDataR(woR,r_old,dim1,n_r_max_old,n_r_maxL,lBc1,l_IC)
                  w(lm,:)=scale_w*woR(1:dim1)
               else
                  w(lm,:)=scale_w*wo(lmo,:)
               end if
            else
               w(lm,:)=zero
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine mapOneField
!------------------------------------------------------------------------------
   subroutine mapDataHydro( wo,zo,po,so,xio,r_old,lm2lmo,n_r_max_old,   &
              &             lm_max_old,n_r_maxL,lBc1,lBc2,lBc3,lBc4,    &
              &             lBc5,w,z,p,s,xi )

      !--- Input variables
      integer,     intent(in) :: n_r_max_old,lm_max_old
      integer,     intent(in) :: n_r_maxL
      logical,     intent(in) :: lBc1,lBc2,lBc3,lBc4,lBc5
      integer,     intent(in) :: lm2lmo(lm_max)
      real(cp),    intent(in) :: r_old(:)
      complex(cp), intent(in) :: wo(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: zo(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: po(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: so(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: xio(lm_max_old,n_r_max_old)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_max,n_r_max),z(lm_max,n_r_max)
      complex(cp), intent(out) :: p(lm_max,n_r_max),s(lm_max,n_r_max)
      complex(cp), intent(out) :: xi(lm_max,n_r_max)

      !--- Local variables
      integer :: lm,lmo,lmStart,lmStop,n_proc
      complex(cp), allocatable :: woR(:),zoR(:)
      complex(cp), allocatable :: poR(:),soR(:)
      complex(cp), allocatable :: xioR(:)

      allocate( woR(n_r_maxL),zoR(n_r_maxL),poR(n_r_maxL) )
      bytes_allocated = bytes_allocated + 3*n_r_maxL*SIZEOF_DEF_COMPLEX
      if ( lreadS .and. l_heat ) then
         allocate( soR(n_r_maxL) )
         bytes_allocated = bytes_allocated + n_r_maxL*SIZEOF_DEF_COMPLEX
      endif
      if ( lreadXi .and. l_chemical_conv ) then
         allocate( xioR(n_r_maxL) )
         bytes_allocated = bytes_allocated + n_r_maxL*SIZEOF_DEF_COMPLEX
      end if

      !PERFON('mD_map')
      do n_proc=0,n_ranks_r-1 ! Blocking of loop over all (l,m)
         lmStart=lm_balance(n_proc)%nStart
         lmStop =lm_balance(n_proc)%nStop

         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_r_max /= n_r_max_old .or. ratio1 /= ratio1_old .or.     &
               &    ratio2 /= ratio2_old .or.                                 &
               &    rscheme_oc%order_boundary /= rscheme_oc_old%order_boundary&
               &    .or. rscheme_oc%version /= rscheme_oc_old%version ) then

                  woR(1:n_r_max_old)=wo(lmo,:)
                  zoR(1:n_r_max_old)=zo(lmo,:)
                  poR(1:n_r_max_old)=po(lmo,:)
                  if ( lreadS .and. l_heat ) soR(1:n_r_max_old)=so(lmo,:)
                  if ( lreadXi .and. l_chemical_conv ) &
                  &            xioR(1:n_r_max_old)=xio(lmo,:)
                  call mapDataR(woR,r_old,n_r_max,n_r_max_old,n_r_maxL,lBc1, &
                       &        .false.)
                  call mapDataR(zoR,r_old,n_r_max,n_r_max_old,n_r_maxL,lBc2, &
                       &        .false.)
                  call mapDataR(poR,r_old,n_r_max,n_r_max_old,n_r_maxL,lBc3, &
                       &        .false.)
                  if ( lreadS .and. l_heat ) then
                     call mapDataR(soR,r_old,n_r_max,n_r_max_old,n_r_maxL,lBc4,&
                          &        .false.)
                  end if
                  if ( lreadXi .and. l_chemical_conv ) then
                     call mapDataR(xioR,r_old,n_r_max,n_r_max_old,n_r_maxL,lBc5,&
                          &        .false.)
                  end if
                  if ( lm > 1 ) then
                     w(lm,:)=scale_v*woR(:)
                     z(lm,:)=scale_v*zoR(:)
                  else
                     w(1,:)=zero
                     z(1,:)=zero
                  end if
                  p(lm,:)=scale_v*poR(:)
                  if ( lreadS .and. l_heat ) s(lm,:)=scale_s*soR(:)
                  if ( lreadXi .and. l_chemical_conv ) &
                  &     xi(lm,:)=scale_xi*xioR(:)
               else
                  if ( lm > 1 ) then
                     w(lm,:)=scale_v*wo(lmo,:)
                     z(lm,:)=scale_v*zo(lmo,:)
                  else
                     w(1,:)=zero
                     z(1,:)=zero
                  end if
                  p(lm,:)=scale_v*po(lmo,:)
                  if ( lreadS .and. l_heat ) s(lm,:)=scale_s*so(lmo,:)
                  if ( lreadXi .and. l_chemical_conv ) &
                  &    xi(lm,:)=scale_xi*xio(lmo,:)
               end if
            else
               w(lm,:)=zero
               z(lm,:)=zero
               p(lm,:)=zero
               if ( lreadS .and. l_heat ) s(lm,:)=zero
               if ( lreadXi .and. l_chemical_conv ) xi(lm,:)=zero
            end if
         end do
      end do
      !PERFOFF
      deallocate(woR,zoR,poR)
      bytes_allocated = bytes_allocated - 3*n_r_maxL*SIZEOF_DEF_COMPLEX
      if ( lreadS .and. l_heat ) then
         deallocate(soR)
         bytes_allocated = bytes_allocated - n_r_maxL*SIZEOF_DEF_COMPLEX
      end if
      if ( lreadXi .and. l_chemical_conv ) then
         deallocate(xioR)
         bytes_allocated = bytes_allocated - n_r_maxL*SIZEOF_DEF_COMPLEX
      end if

   end subroutine mapDataHydro
!------------------------------------------------------------------------------
   subroutine mapDataMag( wo,zo,po,so,r_old,n_rad_tot,n_r_max_old, &
              &           lm_max_old,n_r_maxL,lm2lmo,dim1,l_IC,w,z,p,s )

      !--- Input variables
      integer,     intent(in) :: n_rad_tot,n_r_max_old,lm_max_old
      integer,     intent(in) :: n_r_maxL,dim1
      integer,     intent(in) :: lm2lmo(lm_max)
      logical,     intent(in) :: l_IC
      real(cp),    intent(in) :: r_old(:)
      complex(cp), intent(in) :: wo(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: zo(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: po(lm_max_old,n_r_max_old)
      complex(cp), intent(in) :: so(lm_max_old,n_r_max_old)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_maxMag,dim1),z(lm_maxMag,dim1)
      complex(cp), intent(out) :: p(lm_maxMag,dim1),s(lm_maxMag,dim1)

      !--- Local variables
      integer :: lm,lmo,lmStart,lmStop,n_proc
      complex(cp), allocatable :: woR(:),zoR(:),poR(:),soR(:)

      allocate( woR(n_r_maxL),zoR(n_r_maxL) )
      allocate( poR(n_r_maxL),soR(n_r_maxL) )
      bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DEF_COMPLEX

      !PERFON('mD_map')
      do n_proc=0,n_ranks_r-1 ! Blocking of loop over all (l,m)
         lmStart=lm_balance(n_proc)%nStart
         lmStop =lm_balance(n_proc)%nStop
         lmStart=max(2,lmStart)
         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_rad_tot /= n_r_max_old .or. ratio1 /= ratio1_old .or.    &
               &    ratio2 /= ratio2_old .or.                                  &
               &    rscheme_oc%order_boundary /= rscheme_oc_old%order_boundary &
               &    .or. rscheme_oc%version /= rscheme_oc_old%version ) then
                  woR(1:n_r_max_old)=wo(lmo,:)
                  zoR(1:n_r_max_old)=zo(lmo,:)
                  poR(1:n_r_max_old)=po(lmo,:)
                  soR(1:n_r_max_old)=so(lmo,:)
                  call mapDataR(woR,r_old,dim1,n_r_max_old,n_r_maxL,.false.,l_IC)
                  call mapDataR(zoR,r_old,dim1,n_r_max_old,n_r_maxL,.true.,l_IC)
                  call mapDataR(poR,r_old,dim1,n_r_max_old,n_r_maxL,.true.,l_IC)
                  call mapDataR(soR,r_old,dim1,n_r_max_old,n_r_maxL,.false.,l_IC)
                  w(lm,:)=scale_b*woR(:)
                  z(lm,:)=scale_b*zoR(:)
                  p(lm,:)=scale_b*poR(:)
                  s(lm,:)=scale_b*soR(:)
               else
                  w(lm,:)=scale_b*wo(lmo,:)
                  z(lm,:)=scale_b*zo(lmo,:)
                  p(lm,:)=scale_b*po(lmo,:)
                  s(lm,:)=scale_b*so(lmo,:)
               end if
            else
               w(lm,:)=zero
               z(lm,:)=zero
               p(lm,:)=zero
               s(lm,:)=zero
            end if
         end do
      end do
      !PERFOFF
      deallocate(woR,zoR,poR,soR)
      bytes_allocated = bytes_allocated - 4*n_r_maxL*SIZEOF_DEF_COMPLEX

   end subroutine mapDataMag
!------------------------------------------------------------------------------
   subroutine mapDataR(dataR,r_old,n_rad_tot,n_r_max_old,n_r_maxL,lBc,l_IC)
      !
      !
      !  Copy (interpolate) data (read from disc file) from old grid structure
      !  to new grid. Linear interploation is used in r if the radial grid
      !  structure differs
      !
      !  called in mapdata
      !
      !

      !--- Input variables
      integer,  intent(in) :: n_r_max_old
      integer,  intent(in) :: n_r_maxL,n_rad_tot
      real(cp), intent(in) :: r_old(:)
      logical,  intent(in) :: lBc,l_IC

      !--- Output variables
      complex(cp), intent(inout) :: dataR(:)  ! old data

      !-- Local variables
      integer :: nR, nR_old, n_r_index_start
      real(cp) :: xold(4)
      complex(cp) :: yold(4)
      complex(cp), allocatable :: work(:)
      real(cp) :: cheb_norm_old,scale
      type(costf_odd_t) :: chebt_oc_old


      !-- If **both** the old and the new schemes are Chebyshev, we can
      !-- use costf to get the new data
      !-- Both have also to use the same mapping
      !-- (order_boundary is a proxy of l_map
      if ( rscheme_oc%version == 'cheb' .and. rscheme_oc_old%version == 'cheb' &
      &   .and. rscheme_oc%order_boundary == rscheme_oc_old%order_boundary     &
      &   .and. ratio1 == ratio1_old .and. ratio2 == ratio2_old ) then

         !-- Guess the boundary values, since they have not been stored:
         if ( .not. l_IC .and. lBc ) then
            dataR(1)=two*dataR(2)-dataR(3)
            dataR(n_r_max_old)=two*dataR(n_r_max_old-1)-dataR(n_r_max_old-2)
         end if

         !----- Transform old data to cheb space:
         if ( l_IC ) then
            allocate( work(n_r_maxL) )
            call chebt_oc_old%initialize(n_r_max_old, 2*n_r_maxL+2, 2*n_r_maxL+5)
            call chebt_oc_old%costf1(dataR, work)
            call chebt_oc_old%finalize()
         else
            call rscheme_oc_old%costf1(dataR)
         end if

         !----- Fill up cheb polynomial with zeros:
         if ( n_rad_tot>n_r_max_old ) then
            if ( l_IC) then
               n_r_index_start=n_r_max_old
            else
               n_r_index_start=n_r_max_old+1
            end if
            do nR=n_r_index_start,n_rad_tot
               dataR(nR)=zero
            end do
         end if

         !----- Now transform to new radial grid points:
         if ( l_IC ) then
            call chebt_ic%costf1(dataR,work)
            !----- Rescale :
            cheb_norm_old=sqrt(two/real(n_r_max_old-1,kind=cp))
            scale=cheb_norm_old/cheb_norm_ic

            deallocate( work )
         else
            call rscheme_oc%costf1(dataR)
            !----- Rescale :
            cheb_norm_old=sqrt(two/real(n_r_max_old-1,kind=cp))
            scale=cheb_norm_old/rscheme_oc%rnorm
         end if
         do nR=1,n_rad_tot
            dataR(nR)=scale*dataR(nR)
         end do

      !-- If either the old grid or the new grid is FD, we use a
      !-- polynomial interpolation
      else

         !----- Now transform the inner core:
         if ( l_IC ) then

            allocate( work(n_r_maxL) )

            !-- This is needed for the inner core
            call chebt_oc_old%initialize(n_r_max_old, 2*n_r_maxL+2, 2*n_r_maxL+5)

            !----- Transform old data to cheb space:
            call chebt_oc_old%costf1(dataR, work)

            call chebt_oc_old%finalize()

            !----- Fill up cheb polynomial with zeros:
            if ( n_rad_tot>n_r_max_old ) then
               n_r_index_start=n_r_max_old
               do nR=n_r_index_start,n_rad_tot
                  dataR(nR)=zero
               end do
            end if


            call chebt_ic%costf1(dataR,work)
            !----- Rescale :
            cheb_norm_old=sqrt(two/real(n_r_max_old-1,kind=cp))
            scale=cheb_norm_old/cheb_norm_ic

            deallocate( work )

            do nR=1,n_rad_tot
               dataR(nR)=scale*dataR(nR)
            end do

         else

            allocate( work(n_r_max) )

            !-- Interpolate data and store into a work array
            do nR=1,n_r_max

               nR_old=minloc(abs(r_old-r(nR)),1)
               if ( nR_old < 3 ) nR_old=3
               if ( nR_old == n_r_max_old ) nR_old=n_r_max_old-1

               xold(1)=r_old(nR_old-2)
               xold(2)=r_old(nR_old-1)
               xold(3)=r_old(nR_old)
               xold(4)=r_old(nR_old+1)

               yold(1)=dataR(nR_old-2)
               yold(2)=dataR(nR_old-1)
               yold(3)=dataR(nR_old)
               yold(4)=dataR(nR_old+1)

               call polynomial_interpolation(xold, yold, r(nR), work(nR))

            end do

            !-- Copy interpolated data
            do nR=1,n_r_max
               dataR(nR)=work(nR)
            end do

            deallocate( work )

         end if

      end if

   end subroutine mapDataR
!---------------------------------------------------------------------
   subroutine finish_start_fields(time, minc_old, l_mag_old, omega_ic1Old, &
              &                   omega_ma1Old, z, s, xi, b, omega_ic, omega_ma)

      !-- Input variables
      real(cp), intent(in) :: time
      integer,  intent(in) :: minc_old
      logical,  intent(in) :: l_mag_old
      real(cp), intent(in) :: omega_ic1Old
      real(cp), intent(in) :: omega_ma1Old

      !-- In/Out variables
      complex(cp), intent(inout) :: z(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: b(llmMag:ulmMag,n_r_maxMag)

      !-- Output variables
      real(cp), intent(out) :: omega_ic
      real(cp), intent(out) :: omega_ma

      !-- Local variables
      real(cp) :: fr
      integer :: lm, nR, l1m0


      !-- If mapping towards reduced symmetry, add thermal perturbation in
      !   mode (l,m)=(minc,minc) if parameter tipdipole  /=  0
      if ( l_heat .and. minc<minc_old .and. tipdipole>0.0_cp ) then
         lm=l_max+2
         if ( llm<=lm .and. ulm>=lm ) then
            do nR=1,n_r_max+1
               fr=sin(pi*(r(nR)-r(n_r_max)))
               s(lm,nR)=tipdipole*fr
            end do
         end if
      end if

      if ( l_chemical_conv .and. minc<minc_old .and. tipdipole>0.0_cp ) then
         lm=l_max+2
         if ( llm<=lm .and. ulm>=lm ) then
            do nR=1,n_r_max+1
               fr=sin(pi*(r(nR)-r(n_r_max)))
               xi(lm,nR)=tipdipole*fr
            end do
         end if
      end if

      !-- If starting from data file with longitudinal symmetry, add
      !   weak non-axisymmetric dipole component if tipdipole  /=  0
      if ( ( l_mag .or. l_mag_LF )                                &
      &            .and. minc==1 .and. minc_old/=1 .and.          &
      &            tipdipole>0.0_cp .and. l_mag_old ) then
         lm=l_max+2
         if ( llm<=lm .and. ulm>=lm ) then
            do nR=1,n_r_max+1
               b(lm,nR)=tipdipole
            end do
         end if
      end if

      !----- Set IC and mantle rotation rates:
      !      Following cases are covered:
      !       1) Prescribed inner-core rotation omega_ic_pre
      !       2) Rotation has been read above ( inform >= 4)
      !       3) Rotation calculated from flow field z(l=1,m=0)
      !       4) No rotation
      !       5) Flow driven by prescribed inner core rotation
      !       l_SRIC=.true. (spherical Couette case)
      l1m0=lo_map%lm2(1,0)
      if ( l_rot_ic ) then
         if ( l_SRIC .or. omega_ic1 /= 0.0_cp ) then
            if ( tShift_ic1 == 0.0_cp ) tShift_ic1=tOmega_ic1-time
            if ( tShift_ic2 == 0.0_cp ) tShift_ic2=tOmega_ic2-time
            tOmega_ic1=time+tShift_ic1
            tOmega_ic2=time+tShift_ic2
            omega_ic=omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
            &        omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
            if ( l_master_rank ) then
               write(output_unit,*)
               write(output_unit,*) '! I use prescribed inner core rotation rate:'
               write(output_unit,*) '! omega_ic=',omega_ic
            end if
            if ( kbotv == 2 ) then
               if ( llm<=l1m0 .and. ulm>=l1m0 ) then
                  z(l1m0,n_r_icb)=cmplx(omega_ic/c_z10_omega_ic,0.0_cp,kind=cp)
               end if
            end if
         else
            omega_ic=omega_ic1Old
         end if
      else
         omega_ic=0.0_cp
      end if

      !----- Mantle rotation, same as for inner core (see above)
      !      exept the l_SRIC case.
      if ( l_rot_ma ) then
         if ( l_SRMA .or. omega_ma1 /= 0.0_cp ) then
            if ( tShift_ma1 == 0.0_cp ) tShift_ma1=tOmega_ma1-time
            if ( tShift_ma2 == 0.0_cp ) tShift_ma2=tOmega_ma2-time
            tOmega_ma1=time+tShift_ma1
            tOmega_ma2=time+tShift_ma2
            omega_ma=omega_ma1*cos(omegaOsz_ma1*tOmega_ma1) + &
            &        omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
            if ( l_master_rank ) then
               write(output_unit,*)
               write(output_unit,*) '! I use prescribed mantle rotation rate:'
               write(output_unit,*) '! omega_ma =',omega_ma
               write(output_unit,*) '! omega_ma1=',omega_ma1
            end if
            if ( ktopv == 2 ) then
               if ( llm<=l1m0 .and. ulm>=l1m0 ) then
                  z(l1m0,n_r_cmb)=cmplx(omega_ma/c_z10_omega_ma,0.0_cp,kind=cp)
               end if
            end if
         else
            omega_ma=omega_ma1Old
         end if
      else
         omega_ma=0.0_cp
      end if

   end subroutine finish_start_fields
!------------------------------------------------------------------------------
   subroutine print_info(ra_old,ek_old,pr_old,sc_old,raxi_old,pm_old, &
              &          radratio_old,sigma_ratio_old,n_phi_tot_old,  &
              &          nalias_old, l_max_old)

      !-- Input variables
      real(cp), intent(in) :: ra_old, ek_old, pr_old, sc_old, raxi_old
      real(cp), intent(in) :: pm_old, radratio_old, sigma_ratio_old
      integer,  intent(in) :: nalias_old, l_max_old, n_phi_tot_old

      !---- Compare parameters:
      if ( ra /= ra_old ) &
      &    write(output_unit,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') ra_old,ra
      if ( ek /= ek_old ) &
      &    write(output_unit,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') ek_old,ek
      if ( pr /= pr_old ) &
      &    write(output_unit,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') pr_old,pr
      if ( prmag /= pm_old )                                          &
      &    write(output_unit,'(/,'' ! New mag Pr number (old/new):'',2ES16.6)') &
      &    pm_old,prmag
      if ( raxi /= raxi_old ) &
      &    write(output_unit,'(/,'' ! New composition-based Rayleigh number (old/new):'',2ES16.6)') raxi_old,raxi
      if ( sc /= sc_old ) &
      &    write(output_unit,'(/,'' ! New Schmidt number (old/new):'',2ES16.6)') sc_old,sc
      if ( radratio /= radratio_old )                                    &
      &    write(output_unit,'(/,'' ! New radius ratio (old/new):'',2ES16.6)') &
      &    radratio_old,radratio
      if ( sigma_ratio /= sigma_ratio_old )                             &
      &    write(output_unit,'(/,'' ! New mag cond. ratio (old/new):'',2ES16.6)') &
      &    sigma_ratio_old,sigma_ratio

      if ( n_phi_tot_old /= n_phi_tot) &
      &    write(output_unit,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
      if ( nalias_old /= nalias) &
      &    write(output_unit,*) '! New nalias (old,new)   :',nalias_old,nalias
      if ( l_max_old /= l_max ) &
      &    write(output_unit,*) '! New l_max (old,new)    :',l_max_old,l_max

   end subroutine print_info
!------------------------------------------------------------------------------
end module readCheckPoints
