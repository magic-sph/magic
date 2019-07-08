module readCheckPoints
   !
   ! This module contains the functions that can help reading and
   ! mapping of the restart files
   !

   use precision_mod
   use parallel_mod
   use communications, only: scatter_from_rank0_to_lo
   use truncation, only: n_r_max,lm_max,n_r_maxMag,lm_maxMag,n_r_ic_max, &
       &                 n_r_ic_maxMag,nalias,n_phi_tot,l_max,m_max,     &
       &                 minc,lMagMem,fd_stretch,fd_ratio
   use logic, only: l_rot_ma,l_rot_ic,l_SRIC,l_SRMA,l_cond_ic,l_heat,l_mag, &
       &            l_mag_LF, l_chemical_conv, l_AB1
   use blocking, only: lo_map, lm2l, lm2m, lm_balance, llm, ulm, llmMag, ulmMag
   use init_fields, only: start_file,inform,tOmega_ic1,tOmega_ic2,             &
       &                  tOmega_ma1,tOmega_ma2,omega_ic1,omegaOsz_ic1,        &
       &                  omega_ic2,omegaOsz_ic2,omega_ma1,omegaOsz_ma1,       &
       &                  omega_ma2,omegaOsz_ma2,tShift_ic1,tShift_ic2,        &
       &                  tShift_ma1,tShift_ma2,tipdipole, scale_b, scale_v,   &
       &                  scale_s,scale_xi
   use radial_functions, only: rscheme_oc, chebt_ic, cheb_norm_ic, r
   use num_param, only: alph1, alph2
   use radial_data, only: n_r_icb, n_r_cmb
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio, &
       &                          kbotv, ktopv, sc, raxi
   use constants, only: c_z10_omega_ic, c_z10_omega_ma, pi, zero, two
   use chebyshev, only: type_cheb_odd
   use radial_scheme, only: type_rscheme
   use finite_differences, only: type_fd
   use cosine_transform_odd, only: costf_odd_t
   use useful, only: polynomial_interpolation, abortRun
   use constants, only: one


   implicit none

   private

   logical :: lreadS, lreadXi, lreadR
   logical :: l_axi_old

   integer :: n_start_file
   integer(lip) :: bytes_allocated=0
   class(type_rscheme), pointer :: rscheme_oc_old
   real(cp) :: ratio1_old, ratio2_old, ratio1, ratio2

   public :: readStartFields_old, readStartFields

contains

   subroutine readStartFields_old(w,dwdt,z,dzdt,p,dpdt,s,dsdt,   &
              &                   xi,dxidt,b,dbdt,aj,djdt,b_ic,  &
              &                   dbdt_ic,aj_ic,djdt_ic,omega_ic,&
              &                   omega_ma,lorentz_torque_ic,    &
              &                   lorentz_torque_ma,time,dt_old, &
              &                   dt_new,n_time_step)
      !
      ! This subroutine is used to read the old restart files produced
      ! by MagIC. This is now deprecated with the change of the file format.
      ! This is still needed to read old files.
      !

      !-- Output:
      real(cp),    intent(out) :: time,dt_old,dt_new
      integer,     intent(out) :: n_time_step
      real(cp),    intent(out) :: omega_ic,omega_ma
      real(cp),    intent(out) :: lorentz_torque_ic,lorentz_torque_ma
      complex(cp), intent(out) :: w(llm:ulm,n_r_max),z(llm:ulm,n_r_max)
      complex(cp), intent(out) :: s(llm:ulm,n_r_max),p(llm:ulm,n_r_max)
      complex(cp), intent(out) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dwdt(llm:ulm,n_r_max),dzdt(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dsdt(llm:ulm,n_r_max),dpdt(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dxidt(llm:ulm,n_r_max)
      complex(cp), intent(out) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: djdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: dbdt_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: djdt_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local:
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old
      integer :: n_r_ic_max_old
      real(cp) :: pr_old,ra_old,pm_old
      real(cp) :: raxi_old,sc_old
      real(cp) :: ek_old,radratio_old
      real(cp) :: sigma_ratio_old
      integer :: lm,nR,l1m0
      logical :: l_mag_old
      logical :: startfile_does_exist
      integer :: informOld,ioerr
      integer :: n_r_maxL,n_r_ic_maxL,n_data_oldP,lm_max_old
      integer, allocatable :: lm2lmo(:)

      real(cp) :: fr
      real(cp) :: omega_ic1Old,omegaOsz_ic1Old
      real(cp) :: omega_ic2Old,omegaOsz_ic2Old
      real(cp) :: omega_ma1Old,omegaOsz_ma1Old
      real(cp) :: omega_ma2Old,omegaOsz_ma2Old

      character(len=72) :: rscheme_version_old
      real(cp) :: r_icb_old, r_cmb_old
      integer :: n_in, n_in_2

      complex(cp), allocatable :: wo(:),zo(:),po(:),so(:),xio(:)
      complex(cp), allocatable :: workA(:,:),workB(:,:),workC(:,:)
      complex(cp), allocatable :: workD(:,:),workE(:,:)
      real(cp), allocatable :: r_old(:)

      if ( rscheme_oc%version == 'cheb') then
         ratio1 = alph1
         ratio2 = alph2
      else
         ratio1 = fd_stretch
         ratio2 = fd_ratio
      end if

      if ( rank == 0 ) then
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
            &    time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
            &    informOld,n_r_max_old,n_theta_max_old,n_phi_tot_old,  &
            &    minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
            n_time_step=0
         else if ( inform == 0 ) then
            read(n_start_file)                                         &
            &    time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
            &    n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
            &    minc_old,nalias_old
         else if ( inform == 1 ) then
            read(n_start_file)                                         &
            &    time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
            &    n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
            &    minc_old
            nalias_old=nalias
         else if ( inform >= 2 ) then
            read(n_start_file)                                         &
            &    time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
            &    n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
            &    minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
         end if
         if ( inform == -1 ) inform=informOld

         !---- Compare parameters:
         if ( ra /= ra_old ) &
         &    write(*,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') ra_old,ra
         if ( ek /= ek_old ) &
         &    write(*,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') ek_old,ek
         if ( pr /= pr_old ) &
         &    write(*,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') pr_old,pr
         if ( prmag /= pm_old )                                          &
         &    write(*,'(/,'' ! New mag Pr.number (old/new):'',2ES16.6)') &
         &    pm_old,prmag
         if ( radratio /= radratio_old )                                    &
         &    write(*,'(/,'' ! New mag aspect ratio (old/new):'',2ES16.6)') &
         &    radratio_old,radratio
         if ( sigma_ratio /= sigma_ratio_old )                             &
         &    write(*,'(/,'' ! New mag cond. ratio (old/new):'',2ES16.6)') &
         &    sigma_ratio_old,sigma_ratio

         if ( n_phi_tot_old == 1 ) then ! Axisymmetric restart file
            l_max_old=nalias_old*n_theta_max_old/30
            l_axi_old=.true.
         else
            l_max_old=nalias_old*n_phi_tot_old/60
            l_axi_old=.false.
         end if
         l_mag_old=.false.
         if ( pm_old /= 0.0_cp ) l_mag_old= .true.

         if ( n_phi_tot_old /= n_phi_tot) &
         &    write(*,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
         if ( nalias_old /= nalias) &
         &    write(*,*) '! New nalias (old,new)   :',nalias_old,nalias
         if ( l_max_old /= l_max ) &
         &    write(*,*) '! New l_max (old,new)    :',l_max_old,l_max


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
            ratio1_old=0.0_cp
            ratio2_old=0.0_cp
            allocate ( type_cheb_odd :: rscheme_oc_old )
         end if


         r_icb_old=radratio_old/(one-radratio_old)
         r_cmb_old=one/(one-radratio_old)

         call rscheme_oc_old%initialize(n_r_max_old, n_in, n_in_2)

         call rscheme_oc_old%get_grid(n_r_max_old, r_icb_old, r_cmb_old, &
              &                       ratio1_old, ratio2_old, r_old)

         if ( rscheme_oc%version /= rscheme_oc_old%version )         &
         &    write(*,'(/,'' ! New radial scheme (old/new):'',2A4)') &
         &    rscheme_oc_old%version, rscheme_oc%version

         allocate( lm2lmo(lm_max) )

         call getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
              &         m_max,minc,minc_old,inform,lm_max,   &
              &         lm_max_old,n_data_oldP,lm2lmo)

         ! allocation of local arrays.
         ! if this becomes a performance bottleneck, one can make a module
         ! and allocate the array only once in the initialization
         allocate( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP),so(n_data_oldP) )
         bytes_allocated = bytes_allocated + 4*n_data_oldP*SIZEOF_DEF_COMPLEX
         ! end of allocation

         !PERFON('mD_rd')
         if ( lreadXi ) then
            allocate(xio(n_data_oldP))
            bytes_allocated = bytes_allocated + n_data_oldP*SIZEOF_DEF_COMPLEX
            if ( lreadS ) then
               read(n_start_file) wo, zo, po, so, xio
            else
               read(n_start_file) wo, zo, po, xio
            end if
         else
            allocate(xio(1))
            if ( lreadS ) then
               read(n_start_file) wo, zo, po, so
            else
               read(n_start_file) wo, zo, po
            end if
         end if
         !PERFOFF

      end if ! rank == 0

#ifdef WITH_MPI
      call MPI_Bcast(l_mag_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lreadS,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lreadXi,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(minc_old,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(inform,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(sigma_ratio_old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(time,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt_old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(n_time_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      if ( rank == 0 ) then
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

      if ( rank == 0 ) then
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
      if ( rank == 0 ) then
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
      do nR=1,n_r_max
         !write(*,"(8X,A,I4)") "nR = ",nR
         call scatter_from_rank0_to_lo(workA(:,nR),dwdt(llm:ulm,nR))
         call scatter_from_rank0_to_lo(workB(:,nR),dzdt(llm:ulm,nR))
         call scatter_from_rank0_to_lo(workC(:,nR),dpdt(llm:ulm,nR))
         if ( l_heat ) then
            call scatter_from_rank0_to_lo(workD(:,nR),dsdt(llm:ulm,nR))
         end if

         if ( l_chemical_conv ) then
            call scatter_from_rank0_to_lo(workE(:,nR),dxidt(llm:ulm,nR))
         end if
      end do

      if ( l_heat .and. .not. lreadS ) then ! No entropy before
         dsdt(:,:)=zero
      end if
      if ( l_chemical_conv .and. .not. lreadXi ) then ! No composition before
         dxidt(:,:)=zero
      end if

      deallocate(workA, workB, workC, workD, workE)

      if ( rank == 0 ) then
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

      if ( rank == 0 ) then
         if ( lreadXi ) then
            read(n_start_file) raxi_old, sc_old
            if ( raxi /= raxi_old ) &
              write(*,'(/,'' ! New composition-based Rayleigh number (old/new):'',2ES16.6)') raxi_old,raxi
            if ( sc /= sc_old ) &
              write(*,'(/,'' ! New Schmidt number (old/new):'',2ES16.6)') sc_old,sc
         end if

         if ( l_mag_old ) then
            read(n_start_file) so,wo,zo,po

            if ( l_mag ) then
               call mapDataMag( wo,zo,po,so,r_old,n_data_oldP,n_r_max,n_r_max_old, &
                    &           lm_max_old,n_r_maxL,lm2lmo,n_r_maxMag,             &
                    &           .false.,workA,workB,workC,workD )
            end if
         else
            write(*,*) '! No magnetic data in input file!'
         end if
      end if


      !-- Scatter everything
      if ( l_mag_old .and. l_mag ) then
         do nR=1,n_r_maxMag
            call scatter_from_rank0_to_lo(workA(:,nR),aj(llm:ulm,nR))
            call scatter_from_rank0_to_lo(workB(:,nR),dbdt(llm:ulm,nR))
            call scatter_from_rank0_to_lo(workC(:,nR),djdt(llm:ulm,nR))
            call scatter_from_rank0_to_lo(workD(:,nR),b(llm:ulm,nR))
         end do
      end if

      deallocate( workA, workB, workC, workD )


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


      !-- Inner core part
      !
      !
      if ( l_mag_old ) then
         if ( rank == 0 ) then
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

      if ( rank == 0 ) then
         ! deallocation of the local arrays
         deallocate( lm2lmo )
         deallocate( wo,zo,po,so )
         bytes_allocated = bytes_allocated - 4*n_data_oldP*SIZEOF_DEF_COMPLEX

         if ( lreadXi ) then
            bytes_allocated = bytes_allocated - n_data_oldP*SIZEOF_DEF_COMPLEX
         end if
      end if

      !-- Inner core fields:
      if ( l_mag .or. l_mag_LF ) then
         if ( l_mag_old ) then

            if ( inform >= 2 .and. sigma_ratio_old /= 0.0_cp ) then
               if ( rank == 0 ) then
                  allocate( lm2lmo(lm_max) )
                  call getLm2lmO(n_r_ic_max,n_r_ic_max_old,l_max,l_max_old, &
                       &         m_max,minc,minc_old,inform,lm_max,         &
                       &         lm_max_old,n_data_oldP,lm2lmo)

                  n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)
                  allocate( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP), &
                            so(n_data_oldP) )

                  read(n_start_file) so,wo,zo,po
                  if ( l_mag ) then
                     call mapDataMag( wo,zo,po,so,r_old,n_data_oldP,n_r_ic_max, &
                          &           n_r_ic_max_old,lm_max_old,n_r_ic_maxL,    &
                          &           lm2lmo,n_r_ic_maxMag,.true.,workA,        &
                          &           workB,workC,workD )
                  end if

                  deallocate( lm2lmo )
                  deallocate( wo,zo,po,so )
               end if

               do nR=1,n_r_ic_max
                  call scatter_from_rank0_to_lo(workA(:,nR),aj_ic(llm:ulm,nR))
                  call scatter_from_rank0_to_lo(workB(:,nR),dbdt_ic(llm:ulm,nR))
                  call scatter_from_rank0_to_lo(workC(:,nR),djdt_ic(llm:ulm,nR))
                  call scatter_from_rank0_to_lo(workD(:,nR),b_ic(llm:ulm,nR))
               end do

            else if ( l_cond_ic ) then
               !----- No inner core fields provided by start_file, we thus assume that
               !      simple the inner core field decays like r**(l+1) from
               !      the ICB to r=0:
               if ( rank == 0 ) write(*,'(/,'' ! USING POTENTIAL IC fields!'')')

               do lm=llm,ulm
                  do nR=1,n_r_ic_max
                     b_ic(lm,nR)   =b(lm,n_r_CMB)
                     aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                     dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                     djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
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
      dt_new           =dt_old

      if ( rank == 0 ) then
         deallocate( r_old )
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
            read(n_start_file,iostat=ioerr) lorentz_torque_ic, lorentz_torque_ma
            if( ioerr/=0 ) then
               write(*,*) '! Could not read last line in input file!'
               write(*,*) '! Data missing or wrong format!'
               write(*,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else if ( inform >= 4 .and. inform <= 6 .and. lMagMem == 1 )then
            read(n_start_file,iostat=ioerr) lorentz_torque_ic, &
            &    lorentz_torque_ma,omega_ic,omega_ma
            if( ioerr/=0 ) then
               write(*,*) '! Could not read last line in input file!'
               write(*,*) '! Data missing or wrong format!'
               write(*,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else if ( inform == 7 .or. inform == 8 ) then
            read(n_start_file,iostat=ioerr) lorentz_torque_ic, &
            &    lorentz_torque_ma,                            &
            &    omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1,      &
            &    omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2,      &
            &    omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1,      &
            &    omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2
            if( ioerr/=0 ) then
               write(*,*) '! Could not read last line in input file!'
               write(*,*) '! Data missing or wrong format!'
               write(*,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else if ( inform > 8 ) then
            read(n_start_file,iostat=ioerr) lorentz_torque_ic, &
            &    lorentz_torque_ma,                            &
            &    omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1,      &
            &    omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2,      &
            &    omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1,      &
            &    omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2,      &
            &    dt_new
            if( ioerr/=0 ) then
               write(*,*) '! Could not read last line in input file!'
               write(*,*) '! Data missing or wrong format!'
               write(*,*) '! Change inform accordingly!'
               call abortRun('! Stop run in readStartFields')
            end if
         else
            !-- These could possibly be calcualted from the B-field
            lorentz_torque_ic=0.0_cp
            lorentz_torque_ma=0.0_cp
         end if
         if ( inform < 11 ) then
            lorentz_torque_ic=pm_old*lorentz_torque_ic
            lorentz_torque_ma=pm_old*lorentz_torque_ma
         end if

         if ( l_SRIC ) then
            if ( omega_ic1Old /= omega_ic1 )                       &
            &    write(*,*) '! New IC rotation rate 1 (old/new):', &
            &    omega_ic1Old,omega_ic1
            if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                      &
            &    write(*,*) '! New IC rotation osz. rate 1 (old/new):', &
            &    omegaOsz_ic1Old,omegaOsz_ic1
            if ( omega_ic2Old /= omega_ic2 )                       &
            &    write(*,*) '! New IC rotation rate 2 (old/new):', &
            &    omega_ic2Old,omega_ic2
            if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                      &
            &    write(*,*) '! New IC rotation osz. rate 2 (old/new):', &
            &    omegaOsz_ic2Old,omegaOsz_ic2
         end if
         if ( l_SRMA ) then
            if ( omega_ma1Old /= omega_ma1 )                       &
            &    write(*,*) '! New MA rotation rate 1 (old/new):', &
            &    omega_ma1Old,omega_ma1
            if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                      &
            &    write(*,*) '! New MA rotation osz. rate 1 (old/new):', &
            &    omegaOsz_ma1Old,omegaOsz_ma1
            if ( omega_ma2Old /= omega_ma2 )                       &
            &    write(*,*) '! New MA rotation rate 2 (old/new):', &
            &    omega_ma2Old,omega_ma2
            if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                      &
            &    write(*,*) '! New MA rotation osz. rate 2 (old/new):', &
            &    omegaOsz_ma2Old,omegaOsz_ma2
         end if
      end if ! rank == 0

#ifdef WITH_MPI
      call MPI_Bcast(omega_ic1Old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ma1Old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ic1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ic2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ma1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ma2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt_new,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
#endif

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
            if ( rank == 0 ) then
               write(*,*)
               write(*,*) '! I use prescribed inner core rotation rate:'
               write(*,*) '! omega_ic=',omega_ic
            end if
            if ( kbotv == 2 ) then
               if ( llm<=l1m0 .and. ulm>=l1m0 ) then
                  z(l1m0,n_r_icb)=cmplx(omega_ic/c_z10_omega_ic,0.0_cp,kind=cp)
               end if
            end if
         else if ( inform >= 7 ) then
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
                     omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
            if ( rank == 0 ) then
               write(*,*)
               write(*,*) '! I use prescribed mantle rotation rate:'
               write(*,*) '! omega_ma =',omega_ma
               write(*,*) '! omega_ma1=',omega_ma1
            end if
            if ( ktopv == 2 ) then
               if ( llm<=l1m0 .and. ulm>=l1m0 ) then
                  z(l1m0,n_r_cmb)=cmplx(omega_ma/c_z10_omega_ma,0.0_cp,kind=cp)
               end if
            end if
         else if ( inform >= 7 ) then
            omega_ma=omega_ma1Old
         end if
      else
         omega_ma=0.0_cp
      end if

      if (rank == 0) close(n_start_file)

      ! The broadcast is still needed for old file format
#ifdef WITH_MPI
      call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lorentz_torque_ic,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lorentz_torque_ma,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
#endif

   end subroutine readStartFields_old
!------------------------------------------------------------------------------
   subroutine readStartFields(w,dwdt,z,dzdt,p,dpdt,s,dsdt,   &
              &               xi,dxidt,b,dbdt,aj,djdt,b_ic,  &
              &               dbdt_ic,aj_ic,djdt_ic,omega_ic,&
              &               omega_ma,lorentz_torque_ic,    &
              &               lorentz_torque_ma,time,dt_old, &
              &               dt_new,n_time_step)
      !
      ! This subroutine is used to read the restart files produced
      ! by MagIC.
      !

      !-- Output:
      real(cp),    intent(out) :: time,dt_old,dt_new
      integer,     intent(out) :: n_time_step
      real(cp),    intent(out) :: omega_ic,omega_ma
      real(cp),    intent(out) :: lorentz_torque_ic,lorentz_torque_ma
      complex(cp), intent(out) :: w(llm:ulm,n_r_max),z(llm:ulm,n_r_max)
      complex(cp), intent(out) :: s(llm:ulm,n_r_max),p(llm:ulm,n_r_max)
      complex(cp), intent(out) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dwdt(llm:ulm,n_r_max),dzdt(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dsdt(llm:ulm,n_r_max),dpdt(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dxidt(llm:ulm,n_r_max)
      complex(cp), intent(out) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: djdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: dbdt_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: djdt_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local:
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old
      integer :: n_r_ic_max_old, io_status
      real(cp) :: pr_old,ra_old,pm_old
      real(cp) :: raxi_old,sc_old
      real(cp) :: ek_old,radratio_old
      real(cp) :: sigma_ratio_old
      integer :: lm,nR,l1m0
      logical :: l_mag_old, l_heat_old
      logical :: l_cond_ic_old, l_chemical_conv_old
      logical :: startfile_does_exist
      integer :: n_r_maxL,n_r_ic_maxL,n_data_oldP,lm_max_old
      integer, allocatable :: lm2lmo(:)

      real(cp) :: fr
      real(cp) :: omega_ic1Old,omegaOsz_ic1Old
      real(cp) :: omega_ic2Old,omegaOsz_ic2Old
      real(cp) :: omega_ma1Old,omegaOsz_ma1Old
      real(cp) :: omega_ma2Old,omegaOsz_ma2Old

      character(len=72) :: rscheme_version_old
      real(cp) :: r_icb_old, r_cmb_old
      integer :: n_in, n_in_2, version

      complex(cp), allocatable :: workOld(:,:)
      complex(cp), allocatable :: work(:,:)
      real(cp), allocatable :: r_old(:)

      if ( rscheme_oc%version == 'cheb') then
         ratio1 = alph1
         ratio2 = alph2
      else
         ratio1 = fd_stretch
         ratio2 = fd_ratio
      end if

      if ( rank == 0 ) then
         inquire(file=start_file, exist=startfile_does_exist)

         if ( startfile_does_exist ) then

            !-- First try without stream
            open(newunit=n_start_file, file=start_file, status='old', &
            &    form='unformatted')

            read(n_start_file, iostat=io_status) version

            if ( io_status /= 0 ) then
               write(*,*) '! The checkpoint file does not have record markers'
               write(*,*) '! I try to read it with a stream access...'
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

         if ( index(start_file, 'checkpoint_ave') /= 0 ) then
            write(*,*) '! This is a time-averaged checkpoint'
            write(*,*) '! d#dt arrays will be set to zero'
            l_AB1=.true.
         end if

         read(n_start_file) time, dt_old, n_time_step
         read(n_start_file) ra_old,pr_old,raxi_old,sc_old,pm_old, &
         &                  ek_old,radratio_old,sigma_ratio_old
         read(n_start_file) n_r_max_old,n_theta_max_old,n_phi_tot_old,&
         &                  minc_old,nalias_old,n_r_ic_max_old

         !---- Compare parameters:
         if ( ra /= ra_old ) &
              write(*,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') ra_old,ra
         if ( ek /= ek_old ) &
              write(*,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') ek_old,ek
         if ( pr /= pr_old ) &
              write(*,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') pr_old,pr
         if ( prmag /= pm_old )                                          &
              write(*,'(/,'' ! New mag Pr number (old/new):'',2ES16.6)') &
              pm_old,prmag
            if ( raxi /= raxi_old ) &
              write(*,'(/,'' ! New composition-based Rayleigh number (old/new):'',2ES16.6)') raxi_old,raxi
            if ( sc /= sc_old ) &
              write(*,'(/,'' ! New Schmidt number (old/new):'',2ES16.6)') sc_old,sc
         if ( radratio /= radratio_old )                                    &
              write(*,'(/,'' ! New mag aspect ratio (old/new):'',2ES16.6)') &
              radratio_old,radratio
         if ( sigma_ratio /= sigma_ratio_old )                             &
              write(*,'(/,'' ! New mag cond. ratio (old/new):'',2ES16.6)') &
              sigma_ratio_old,sigma_ratio

         if ( n_phi_tot_old == 1 ) then ! Axisymmetric restart file
            l_max_old=nalias_old*n_theta_max_old/30
            l_axi_old=.true.
         else
            l_max_old=nalias_old*n_phi_tot_old/60
            l_axi_old=.false.
         end if

         if ( n_phi_tot_old /= n_phi_tot) &
              write(*,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
         if ( nalias_old /= nalias) &
              write(*,*) '! New nalias (old,new)   :',nalias_old,nalias
         if ( l_max_old /= l_max ) &
              write(*,*) '! New l_max (old,new)    :',l_max_old,l_max


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

         call rscheme_oc_old%get_grid(n_r_max_old, r_icb_old, r_cmb_old, &
              &                       ratio1_old, ratio2_old, r_old)

         if ( rscheme_oc%version /= rscheme_oc_old%version )         &
         &    write(*,'(/,'' ! New radial scheme (old/new):'',2A4)') &
         &    rscheme_oc_old%version, rscheme_oc%version

         !-- Determine the old mapping
         allocate( lm2lmo(lm_max) )
         call getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
              &         m_max,minc,minc_old,4,lm_max,        &
              &         lm_max_old,n_data_oldP,lm2lmo)
         n_r_maxL = max(n_r_max,n_r_max_old)

         !-- Read Lorentz torques and rotation rates:
         read(n_start_file) lorentz_torque_ic, lorentz_torque_ma,    &
         &                  omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
         &                  omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
         &                  omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
         &                  omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2, &
         &                  dt_new

         if ( l_SRIC ) then
            if ( omega_ic1Old /= omega_ic1 )                       &
            &    write(*,*) '! New IC rotation rate 1 (old/new):', &
            &    omega_ic1Old,omega_ic1
            if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                      &
            &    write(*,*) '! New IC rotation osz. rate 1 (old/new):', &
            &    omegaOsz_ic1Old,omegaOsz_ic1
            if ( omega_ic2Old /= omega_ic2 )                       &
            &    write(*,*) '! New IC rotation rate 2 (old/new):', &
            &    omega_ic2Old,omega_ic2
            if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                      &
            &    write(*,*) '! New IC rotation osz. rate 2 (old/new):', &
            &    omegaOsz_ic2Old,omegaOsz_ic2
         end if
         if ( l_SRMA ) then
            if ( omega_ma1Old /= omega_ma1 )                       &
            &    write(*,*) '! New MA rotation rate 1 (old/new):', &
            &    omega_ma1Old,omega_ma1
            if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                      &
            &    write(*,*) '! New MA rotation osz. rate 1 (old/new):', &
            &    omegaOsz_ma1Old,omegaOsz_ma1
            if ( omega_ma2Old /= omega_ma2 )                       &
            &    write(*,*) '! New MA rotation rate 2 (old/new):', &
            &    omega_ma2Old,omega_ma2
            if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                      &
            &    write(*,*) '! New MA rotation osz. rate 2 (old/new):', &
            &    omegaOsz_ma2Old,omegaOsz_ma2
         end if

         read(n_start_file) l_heat_old, l_chemical_conv_old, l_mag_old, &
         &                  l_cond_ic_old

      end if ! rank == 0

#ifdef WITH_MPI
      call MPI_Bcast(l_mag_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_heat_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_chemical_conv_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_cond_ic_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(minc_old,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(time,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ic1Old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ma1Old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ic1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ic2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ma1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ma2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt_old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(n_time_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt_new,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lorentz_torque_ic,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lorentz_torque_ma,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
#endif

      !-- Allocate arrays
      if ( rank == 0 ) then
         allocate( work(lm_max,n_r_max), workOld(lm_max_old,n_r_max_old) )
      else
         allocate( work(1,n_r_max), workOld(1,1) )
      end if

      !-- Read the poloidal flow
      if ( rank == 0 ) then
         work(:,:)=zero
         read(n_start_file) workOld
         call mapOneField( workOld,scale_v,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,n_r_max,.false.,.false.,work )
         !-- Cancel the spherically symmetric part for poloidal flow
         work(1,:)=zero
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(work(:,nR),w(llm:ulm,nR))
      end do

      !-- Read dwdt
      if ( rank == 0 ) then
         work(:,:)=zero
         read(n_start_file) workOld
         call mapOneField( workOld,scale_v,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,n_r_max,.true.,.false.,work )
         !-- Cancel the spherically symmetric part for poloidal flow
         work(1,:)=zero
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(work(:,nR),dwdt(llm:ulm,nR))
      end do

      !-- Read the toroidal flow
      if ( rank == 0 ) then
         work(:,:)=zero
         read(n_start_file) workOld
         call mapOneField( workOld,scale_v,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,n_r_max,.false.,.false.,work )
         !-- Cancel the spherically symmetric part for toroidal flow
         work(1,:)=zero
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(work(:,nR),z(llm:ulm,nR))
      end do

      !-- Read dzdt
      if ( rank == 0 ) then
         work(:,:)=zero
         read(n_start_file) workOld
         call mapOneField( workOld,scale_v,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,n_r_max,.true.,.false.,work )
         !-- Cancel the spherically symmetric part for toroidal flow
         work(1,:)=zero
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(work(:,nR),dzdt(llm:ulm,nR))
      end do

      !-- Read the pressure
      if ( rank == 0 ) then
         work(:,:)=zero
         read(n_start_file) workOld
         call mapOneField( workOld,scale_v,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,n_r_max,.false.,.false.,work )
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(work(:,nR),p(llm:ulm,nR))
      end do

      !-- Read dpdt
      if ( rank == 0 ) then
         work(:,:)=zero
         read(n_start_file) workOld
         call mapOneField( workOld,scale_v,r_old,lm2lmo,n_r_max_old, &
              &            n_r_maxL,n_r_max,.true.,.false.,work )
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(work(:,nR),dpdt(llm:ulm,nR))
      end do

      if ( l_heat_old ) then
         !-- Read the entropy
         if ( rank == 0 ) then
            work(:,:)=zero
            read(n_start_file) workOld
            call mapOneField( workOld,scale_s,r_old,lm2lmo,n_r_max_old, &
                 &            n_r_maxL,n_r_max,.false.,.false.,work )
         end if
         if ( l_heat ) then
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),s(llm:ulm,nR))
            end do
         end if

         !-- Read dsdt
         if ( rank == 0 ) then
            work(:,:)=zero
            read(n_start_file) workOld
            call mapOneField( workOld,scale_s,r_old,lm2lmo,n_r_max_old, &
                 &            n_r_maxL,n_r_max,.true.,.false.,work )
         end if
         if ( l_heat ) then
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),dsdt(llm:ulm,nR))
            end do
         end if
      end if

      if ( l_chemical_conv_old ) then
         !-- Read the chemical composition
         if ( rank == 0 ) then
            work(:,:)=zero
            read(n_start_file) workOld
            call mapOneField( workOld,scale_xi,r_old,lm2lmo,n_r_max_old, &
                 &            n_r_maxL,n_r_max,.false.,.false.,work )
         end if
         if ( l_chemical_conv ) then
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),xi(llm:ulm,nR))
            end do
         end if

         !-- Read dsdt
         if ( rank == 0 ) then
            work(:,:)=zero
            read(n_start_file) workOld
            call mapOneField( workOld,scale_xi,r_old,lm2lmo,n_r_max_old, &
                 &            n_r_maxL,n_r_max,.true.,.false.,work )
         end if
         if ( l_chemical_conv ) then
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),dxidt(llm:ulm,nR))
            end do
         end if
      end if

      if ( l_heat .and. .not. l_heat_old ) then ! No entropy before
         s(:,:)   =zero
         dsdt(:,:)=zero
      end if
      if ( l_chemical_conv .and. .not. l_chemical_conv_old ) then ! No composition before
         xi(:,:)   =zero
         dxidt(:,:)=zero
      end if

      if ( l_mag .or. l_mag_LF ) then

         if ( l_mag_old ) then
            !-- Read the poloidal magnetic field
            if ( rank == 0 ) then
               work(:,:)=zero
               read(n_start_file) workOld
               call mapOneField(workOld,scale_b,r_old,lm2lmo,n_r_max_old, &
                    &           n_r_maxL,n_r_max,.false.,.false.,work)
               !-- Cancel the spherically-symmetric part
               work(1,:)=zero
            end if
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),b(llm:ulm,nR))
            end do

            !-- Read dbdt
            if ( rank == 0 ) then
               work(:,:)=zero
               read(n_start_file) workOld
               call mapOneField(workOld,scale_b,r_old,lm2lmo,n_r_max_old,   &
                    &           n_r_maxL,n_r_max,.true.,.false.,work)
               !-- Cancel the spherically-symmetric part
               work(1,:)=zero
            end if
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),dbdt(llm:ulm,nR))
            end do

            !-- Read the toroidal magnetic field
            if ( rank == 0 ) then
               work(:,:)=zero
               read(n_start_file) workOld
               call mapOneField(workOld,scale_b,r_old,lm2lmo,n_r_max_old,   &
                    &           n_r_maxL,n_r_max,.false.,.false.,work)
               !-- Cancel the spherically-symmetric part
               work(1,:)=zero
            end if
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),aj(llm:ulm,nR))
            end do

            !-- Read djdt
            if ( rank == 0 ) then
               work(:,:)=zero
               read(n_start_file) workOld
               call mapOneField(workOld,scale_b,r_old,lm2lmo,n_r_max_old,   &
                    &           n_r_maxL,n_r_max,.true.,.false.,work)
               !-- Cancel the spherically-symmetric part
               work(1,:)=zero
            end if
            do nR=1,n_r_max
               call scatter_from_rank0_to_lo(work(:,nR),djdt(llm:ulm,nR))
            end do

            if ( l_cond_ic ) then

               if ( l_cond_ic_old ) then
                  deallocate( work, workOld )

                  if ( rank == 0 ) then
                     deallocate( lm2lmo )
                     allocate( lm2lmo(lm_max) )
                     call getLm2lmO(n_r_ic_max,n_r_ic_max_old,l_max,l_max_old, &
                          &         m_max,minc,minc_old,4,lm_max,              &
                          &         lm_max_old,n_data_oldP,lm2lmo)

                     n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)
                     allocate( work(lm_max,n_r_ic_max) )
                     allocate( workOld(lm_max_old,n_r_ic_max_old) )
                  else
                     allocate( work(1,n_r_ic_max), workOld(1,1) )
                  end if

                  !-- Read the inner core poloidal magnetic field
                  if ( rank == 0 ) then
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
                  if ( rank == 0 ) then
                     work(:,:)=zero
                     read(n_start_file) workOld
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,    &
                          &            n_r_ic_max_old,n_r_ic_maxL,      &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  do nR=1,n_r_ic_max
                     call scatter_from_rank0_to_lo(work(:,nR),dbdt_ic(llm:ulm,nR))
                  end do

                  !-- Read the inner core toroidal magnetic field
                  if ( rank == 0 ) then
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
                  if ( rank == 0 ) then
                     work(:,:)=zero
                     read(n_start_file) workOld
                     call mapOneField( workOld,scale_b,r_old,lm2lmo,     &
                          &            n_r_ic_max_old,n_r_ic_maxL,       &
                          &            n_r_ic_max,.true.,.true.,work )
                     !-- Cancel the spherically-symmetric part
                     work(1,:)=zero
                  end if
                  do nR=1,n_r_ic_max
                     call scatter_from_rank0_to_lo(work(:,nR),djdt_ic(llm:ulm,nR))
                  end do

                  if ( rank == 0 ) deallocate( lm2lmo )

               else
                  !-- No inner core fields provided by start_file, we thus assume that
                  !   simple the inner core field decays like r**(l+1) from
                  !   the ICB to r=0:
                  if ( rank == 0 ) write(*,'(/,'' ! USING POTENTIAL IC fields!'')')

                  do lm=llm,ulm
                     do nR=1,n_r_ic_max
                        b_ic(lm,nR)   =b(lm,n_r_CMB)
                        aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                        dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                        djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
                     end do
                  end do

               end if

            end if

         else
            if ( rank == 0 ) write(*,*) '! No magnetic data in input file!'
         end if

      end if


      !-- Free memory
      deallocate( work, workOld )

      !-- Close file
      if ( rank == 0 ) then
         deallocate( r_old )
         call rscheme_oc_old%finalize() ! deallocate old radial scheme
         close(n_start_file)
      end if ! rank == 0


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
            if ( rank == 0 ) then
               write(*,*)
               write(*,*) '! I use prescribed inner core rotation rate:'
               write(*,*) '! omega_ic=',omega_ic
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
            if ( rank == 0 ) then
               write(*,*)
               write(*,*) '! I use prescribed mantle rotation rate:'
               write(*,*) '! omega_ma =',omega_ma
               write(*,*) '! omega_ma1=',omega_ma1
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

   end subroutine readStartFields
!------------------------------------------------------------------------------
   subroutine getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
              &         m_max,minc,minc_old,inform,lm_max,   &
              &         lm_max_old,n_data_oldP,lm2lmo)

      !--- Input variables
      integer, intent(in) :: n_r_max,l_max,m_max,minc
      integer, intent(in) :: n_r_max_old,l_max_old,minc_old
      integer, intent(in) :: inform,lm_max

      !--- Output variables
      integer,intent(out) :: lm2lmo(lm_max)
      integer,intent(out) :: n_data_oldP
      integer,intent(out) :: lm_max_old

      !--- Local variables
      integer :: n_data,n_data_old
      integer :: m_max_old
      integer :: l,m,lm,lmo,lo,mo

      !-- Outer core fields:
      n_data  = lm_max*n_r_max

      if ( .not. l_axi_old ) then
         m_max_old=(l_max_old/minc_old)*minc_old
      else
         m_max_old=0
      end if

      if ( l_max==l_max_old .and. minc==minc_old .and. n_r_max==n_r_max_old &
      &    .and. m_max==m_max_old ) then

         !----- Direct reading of fields, grid not changed:
         write(*,'(/,'' ! Reading fields directly.'')')

         n_data_old=n_data
         if ( inform>2 ) then
            n_data_oldP=n_data
         else
            !----- In the past an 'extra' radial grid point has been
            !      stored which was not really necessary
            n_data_oldP=lm_max*(n_r_max+1)
         end if

         lm_max_old=lm_max
      else

         !----- Mapping onto new grid !
         write(*,'(/,'' ! Mapping onto new grid.'')')

         if ( mod(minc_old,minc) /= 0 )                                &
         &     write(*,'('' ! Warning: Incompatible old/new minc= '',2i3)')

         lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
         &          m_max_old*(m_max_old-minc_old)/(2*minc_old) +     &
         &          l_max_old-m_max_old+1

         n_data_old=lm_max_old*n_r_max_old
         if ( inform>2 ) then
            n_data_oldP=n_data_old
         else
            n_data_oldP=lm_max_old*(n_r_max_old+1)
         end if

         !-- Write info to STdoUT:
         write(*,'('' ! Old/New  l_max= '',2I4,''  m_max= '',2I4,     &
         &            ''  minc= '',2I3,''  lm_max= '',2I8/)')         &
         &                l_max_old,l_max,m_max_old,m_max,            &
         &                minc_old,minc,lm_max_old,lm_max
         if ( n_r_max_old /= n_r_max )                                &
         &   write(*,'('' ! Old/New n_r_max='',2i4)') n_r_max_old,n_r_max

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
      do n_proc=0,n_procs-1 ! Blocking of loop over all (l,m)
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
      complex(cp), intent(in) :: wo(:),zo(:)
      complex(cp), intent(in) :: po(:),so(:)
      real(cp),    intent(in) :: r_old(:)
      complex(cp), intent(in) :: xio(:)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_max,n_r_max),z(lm_max,n_r_max)
      complex(cp), intent(out) :: p(lm_max,n_r_max),s(lm_max,n_r_max)
      complex(cp), intent(out) :: xi(lm_max,n_r_max)

      !--- Local variables
      integer :: lm,lmo,n,nR,lmStart,lmStop,n_proc
      complex(cp),allocatable :: woR(:),zoR(:)
      complex(cp),allocatable :: poR(:),soR(:)
      complex(cp),allocatable :: xioR(:)

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
      do n_proc=0,n_procs-1 ! Blocking of loop over all (l,m)
         lmStart=lm_balance(n_proc)%nStart
         lmStop =lm_balance(n_proc)%nStop

         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_r_max /= n_r_max_old .or. ratio1 /= ratio1_old .or.     &
               &    ratio2 /= ratio2_old .or.                                 &
               &    rscheme_oc%order_boundary /= rscheme_oc_old%order_boundary&
               &    .or. rscheme_oc%version /= rscheme_oc_old%version ) then

                  do nR=1,n_r_max_old  ! copy on help arrays
                     n=lmo+(nR-1)*lm_max_old
                     woR(nR)=wo(n)
                     zoR(nR)=zo(n)
                     poR(nR)=po(n)
                     if ( lreadS .and. l_heat ) soR(nR)=so(n)
                     if ( lreadXi .and. l_chemical_conv ) xioR(nR)=xio(n)
                  end do
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
                  do nR=1,n_r_max
                     if ( lm > 1 ) then
                        w(lm,nR)=scale_v*woR(nR)
                        z(lm,nR)=scale_v*zoR(nR)
                     else
                        w(1,nR)=zero
                        z(1,nR)=zero
                     end if
                     p(lm,nR)=scale_v*poR(nR)
                     if ( lreadS .and. l_heat ) s(lm,nR)=scale_s*soR(nR)
                     if ( lreadXi .and. l_chemical_conv ) &
                     &     xi(lm,nR)=scale_xi*xioR(nR)
                  end do
               else
                  do nR=1,n_r_max
                     n=lmo+(nR-1)*lm_max_old
                     if ( lm > 1 ) then
                        w(lm,nR)=scale_v*wo(n)
                        z(lm,nR)=scale_v*zo(n)
                     else
                        w(1,nR)=zero
                        z(1,nR)=zero
                     end if
                     p(lm,nR)=scale_v*po(n)
                     if ( lreadS .and. l_heat ) s(lm,nR)=scale_s*so(n)
                     if ( lreadXi .and. l_chemical_conv ) &
                     &    xi(lm,nR)=scale_xi*xio(n)
                  end do
               end if
            else
               do nR=1,n_r_max
                  w(lm,nR)=zero
                  z(lm,nR)=zero
                  p(lm,nR)=zero
                  if ( lreadS .and. l_heat ) s(lm,nR)=zero
                  if ( lreadXi .and. l_chemical_conv ) xi(lm,nR)=zero

               end do
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
   subroutine mapDataMag( wo,zo,po,so,r_old,n_data_oldP,n_rad_tot,n_r_max_old, &
              &           lm_max_old,n_r_maxL,lm2lmo,dim1,l_IC,w,z,p,s )

      !--- Input variables
      integer,     intent(in) :: n_rad_tot,n_r_max_old,lm_max_old
      integer,     intent(in) :: n_r_maxL,n_data_oldP,dim1
      integer,     intent(in) :: lm2lmo(lm_max)
      logical,     intent(in) :: l_IC
      complex(cp), intent(in) :: wo(n_data_oldP),zo(n_data_oldP)
      real(cp),    intent(in) :: r_old(:)
      complex(cp), intent(in) :: po(n_data_oldP),so(n_data_oldP)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_maxMag,dim1),z(lm_maxMag,dim1)
      complex(cp), intent(out) :: p(lm_maxMag,dim1),s(lm_maxMag,dim1)

      !--- Local variables
      integer :: lm,lmo,n,nR,lmStart,lmStop,n_proc
      complex(cp), allocatable :: woR(:),zoR(:),poR(:),soR(:)

      allocate( woR(n_r_maxL),zoR(n_r_maxL) )
      allocate( poR(n_r_maxL),soR(n_r_maxL) )
      bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DEF_COMPLEX

      !PERFON('mD_map')
      do n_proc=0,n_procs-1 ! Blocking of loop over all (l,m)
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
                  do nR=1,n_r_max_old  ! copy on help arrays
                     n=lmo+(nR-1)*lm_max_old
                     woR(nR)=wo(n)
                     zoR(nR)=zo(n)
                     poR(nR)=po(n)
                     soR(nR)=so(n)
                  end do
                  call mapDataR(woR,r_old,dim1,n_r_max_old,n_r_maxL,.false.,l_IC)
                  call mapDataR(zoR,r_old,dim1,n_r_max_old,n_r_maxL,.true.,l_IC)
                  call mapDataR(poR,r_old,dim1,n_r_max_old,n_r_maxL,.true.,l_IC)
                  call mapDataR(soR,r_old,dim1,n_r_max_old,n_r_maxL,.false.,l_IC)
                  do nR=1,n_rad_tot
                     w(lm,nR)=scale_b*woR(nR)
                     z(lm,nR)=scale_b*zoR(nR)
                     p(lm,nR)=scale_b*poR(nR)
                     s(lm,nR)=scale_b*soR(nR)
                  end do
               else
                  do nR=1,n_rad_tot
                     n=lmo+(nR-1)*lm_max_old
                     w(lm,nR)=scale_b*wo(n)
                     z(lm,nR)=scale_b*zo(n)
                     p(lm,nR)=scale_b*po(n)
                     s(lm,nR)=scale_b*so(n)
                  end do
               end if
            else
               do nR=1,n_rad_tot
                  w(lm,nR)=zero
                  z(lm,nR)=zero
                  p(lm,nR)=zero
                  s(lm,nR)=zero
               end do
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
end module readCheckPoints
