module readCheckPoints
   !
   ! This module contains the functions that can help reading and 
   ! mapping of the restart files
   !

   use precision_mod
   use truncation, only: n_r_max,lm_max,n_r_maxMag,lm_maxMag,n_r_ic_max, &
       &                 n_r_ic_maxMag,nalias,n_phi_tot,l_max,m_max,     &
       &                 minc,lMagMem
   use logic, only: l_rot_ma,l_rot_ic,l_SRIC,l_SRMA,l_cond_ic,l_heat,l_mag, &
       &            l_mag_LF, l_chemical_conv
   use blocking, only: lm2,lmStartB,lmStopB,nLMBs,lm2l,lm2m
   use init_fields, only: start_file,n_start_file,inform,tOmega_ic1,tOmega_ic2,&
       &                  tOmega_ma1,tOmega_ma2,omega_ic1,omegaOsz_ic1,        &
       &                  omega_ic2,omegaOsz_ic2,omega_ma1,omegaOsz_ma1,       &
       &                  omega_ma2,omegaOsz_ma2,tShift_ic1,tShift_ic2,        &
       &                  tShift_ma1,tShift_ma2,tipdipole, scale_b, scale_v,   &
       &                  scale_s,scale_xi
   use radial_functions, only: chebt_oc, cheb_norm, chebt_ic, cheb_norm_ic, r
   use radial_data, only: n_r_icb, n_r_cmb
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio, &
       &                          kbotv, ktopv, sc, raxi
   use constants, only: c_z10_omega_ic, c_z10_omega_ma, pi, zero, two
   use cosine_transform_odd


   implicit none

   private

   logical :: lreadS,lreadXi
   logical :: l_axi_old

   integer(lip) :: bytes_allocated=0

#ifdef WITH_HDF5
   public :: readStartFields, readHdf5_serial
#else
   public :: readStartFields
#endif

contains

   subroutine readStartFields(w,dwdt,z,dzdt,p,dpdt,s,dsdt,   &
              &               xi,dxidt,b,dbdt,aj,djdt,b_ic,  &
              &               dbdt_ic,aj_ic,djdt_ic,omega_ic,&
              &               omega_ma,lorentz_torque_ic,    &
              &               lorentz_torque_ma,time,dt_old, &
              &               dt_new,n_time_step)
      !
      !   read initial condition from restart file 
      !

      !-- Output:
      real(cp),    intent(out) :: time,dt_old,dt_new
      integer,     intent(out) :: n_time_step
      real(cp),    intent(out) :: omega_ic,omega_ma
      real(cp),    intent(out) :: lorentz_torque_ic,lorentz_torque_ma
      complex(cp), intent(out) :: w(lm_max,n_r_max),z(lm_max,n_r_max)
      complex(cp), intent(out) :: s(lm_max,n_r_max),p(lm_max,n_r_max)
      complex(cp), intent(out) :: xi(lm_max,n_r_max)
      complex(cp), intent(out) :: dwdt(lm_max,n_r_max),dzdt(lm_max,n_r_max)
      complex(cp), intent(out) :: dsdt(lm_max,n_r_max),dpdt(lm_max,n_r_max)
      complex(cp), intent(out) :: dxidt(lm_max,n_r_max)
      complex(cp), intent(out) :: b(lm_maxMag,n_r_maxMag),aj(lm_maxMag,n_r_maxMag)
      complex(cp), intent(out) :: dbdt(lm_maxMag,n_r_maxMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,n_r_maxMag)
      complex(cp), intent(out) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: dbdt_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: djdt_ic(lm_maxMag,n_r_ic_maxMag)

      !-- Local:
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old
      integer :: n_r_ic_max_old
      real(cp) :: pr_old,ra_old,pm_old
      real(cp) :: raxi_old,sc_old
      real(cp) :: ek_old,radratio_old
      real(cp) :: sigma_ratio_old
      integer :: nLMB,lm,lmStart,lmStop,nR,l1m0
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

      complex(cp), allocatable :: wo(:),zo(:),po(:),so(:),xio(:)

      inquire(file=start_file, exist=startfile_does_exist)
    
      if ( startfile_does_exist ) then
         open(newunit=n_start_file, file=start_file, status='old', &
         &    form='unformatted')
      else
         write(*,*)
         write(*,*) '! The restart file does not exist !'
         stop
      end if
    
      sigma_ratio_old=0.0_cp  ! assume non conducting inner core !
      if ( inform == -1 ) then ! This is default !
         read(n_start_file)                                         &
              time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
              informOld,n_r_max_old,n_theta_max_old,n_phi_tot_old,  &
              minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
         n_time_step=0
      else if ( inform == 0 ) then
         read(n_start_file)                                         &
              time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
              n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
              minc_old,nalias_old
      else if ( inform == 1 ) then
         read(n_start_file)                                         &
              time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
              n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
              minc_old
         nalias_old=nalias
      else if ( inform >= 2 ) then
         read(n_start_file)                                         &
              time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
              n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
              minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
      end if
      if ( inform == -1 ) inform=informOld
    
      !---- Compare parameters:
      if ( ra /= ra_old ) &
           write(*,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') ra_old,ra
      if ( ek /= ek_old ) &
           write(*,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') ek_old,ek
      if ( pr /= pr_old ) &
           write(*,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') pr_old,pr
      if ( prmag /= pm_old )                                          &
           write(*,'(/,'' ! New mag Pr.number (old/new):'',2ES16.6)') &
           pm_old,prmag
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
      l_mag_old=.false.
      if ( pm_old /= 0.0_cp ) l_mag_old= .true. 
    
      if ( n_phi_tot_old /= n_phi_tot) &
           write(*,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
      if ( nalias_old /= nalias) &
           write(*,*) '! New nalias (old,new)   :',nalias_old,nalias
      if ( l_max_old /= l_max ) &
           write(*,*) '! New l_max (old,new)    :',l_max_old,l_max

    
      if ( inform==6 .or. inform==7 .or. inform==9 .or. inform==11 .or. &
           inform==13) then
         lreadS=.false.
      else
         lreadS=.true.
      end if

      if ( inform==13 .or. inform==14 ) then
         lreadXi=.true.
      else
         lreadXi=.false.
      end if
    
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
         if ( lreadS ) then
            read(n_start_file) wo, zo, po, so
         else
            read(n_start_file) wo, zo, po
         end if
      end if
      !PERFOFF
    
      n_r_maxL = max(n_r_max,n_r_max_old)
    
      call mapDataHydro( wo,zo,po,so,xio,n_data_oldP,lm2lmo,  &
           &            n_r_max_old,lm_max_old,n_r_maxL,      &
           &            .false.,.false.,.false.,.false.,      &
           &            .false.,w,z,p,s,xi )

      if ( l_heat .and. .not. lreadS ) then ! No entropy before
         s(:,:)=zero
      end if
      if ( l_chemical_conv .and. .not. lreadXi ) then ! No composition before
         xi(:,:)=zero
      end if
    
      !PERFON('mD_rd_dt')
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
    
      call mapDataHydro( wo,zo,po,so,xio,n_data_oldP,lm2lmo,     &
           &             n_r_max_old,lm_max_old,n_r_maxL,.true., &
           &             .true.,.true.,.true.,.true.,dwdt,dzdt,  &
           &             dpdt,dsdt,dxidt )

      if ( l_heat .and. .not. lreadS ) then ! No entropy before
         dsdt(:,:)=zero
      end if
      if ( l_chemical_conv .and. .not. lreadXi ) then ! No composition before
         dxidt(:,:)=zero
      end if

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
            call mapDataMag( wo,zo,po,so,n_data_oldP,n_r_max,n_r_max_old, &
                 &           lm_max_old,n_r_maxL,lm2lmo,n_r_maxMag,    &
                 &           .false.,aj,dbdt,djdt,b )
         end if
      else
         write(*,*) '! No magnetic data in input file!'
      end if
    
    
      !-- If mapping towards reduced symmetry, add thermal perturbation in
      !   mode (l,m)=(minc,minc) if parameter tipdipole  /=  0
      if ( l_heat .and. minc<minc_old .and. tipdipole>0.0_cp ) then
         do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
            lmStart=lmStartB(nLMB)
            lmStop =lmStopB(nLMB)
            lm=l_max+2
            if ( lmStart<=lm .and. lmStop>=lm ) then
               do nR=1,n_r_max+1
                  fr=sin(pi*(r(nR)-r(n_r_max)))
                  s(lm,nR)=tipdipole*fr
               end do
            end if
         end do
      end if

      if ( l_chemical_conv .and. minc<minc_old .and. tipdipole>0.0_cp ) then
         do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
            lmStart=lmStartB(nLMB)
            lmStop =lmStopB(nLMB)
            lm=l_max+2
            if ( lmStart<=lm .and. lmStop>=lm ) then
               do nR=1,n_r_max+1
                  fr=sin(pi*(r(nR)-r(n_r_max)))
                  xi(lm,nR)=tipdipole*fr
               end do
            end if
         end do
      end if
    
      !-- If starting from data file with longitudinal symmetry, add
      !   weak non-axisymmetric dipole component if tipdipole  /=  0
      if ( ( l_mag .or. l_mag_LF )                                &
           &       .and. minc==1 .and. minc_old/=1 .and.          &
           &       tipdipole>0.0_cp .and. l_mag_old ) then
         do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
            lmStart=lmStartB(nLMB)
            lmStop =lmStopB(nLMB)
            lm=l_max+2
            if ( lmStart<=lm .and. lmStop>=lm ) then
               do nR=1,n_r_max+1
                  b(lm,nR)=tipdipole
               end do
            end if
         end do
      end if
    
      ! deallocation of the local arrays
      deallocate( lm2lmo )
      deallocate( wo,zo,po,so )
      bytes_allocated = bytes_allocated - 4*n_data_oldP*SIZEOF_DEF_COMPLEX

      if ( lreadXi ) then
         bytes_allocated = bytes_allocated - n_data_oldP*SIZEOF_DEF_COMPLEX
      end if
    
      !call mapData(n_r_max_old,l_max_old,minc_old,l_mag_old, &
      !     w,dwdt,z,dzdt,p,dpdt,s,dsdt,b,dbdt,aj,djdt)
    
      !-- Inner core fields:
      if ( l_mag_old ) then
         if ( inform >= 2 .and. sigma_ratio_old /= 0.0_cp ) then
            allocate( lm2lmo(lm_max) )
            call getLm2lmO(n_r_ic_max,n_r_ic_max_old,l_max,l_max_old, &
                 &         m_max,minc,minc_old,inform,lm_max,         &
                 &         lm_max_old,n_data_oldP,lm2lmo)
    
            n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)
            allocate( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP), &
                      so(n_data_oldP) )
    
            read(n_start_file) so,wo,zo,po
            if ( l_mag ) then
               call mapDataMag( wo,zo,po,so,n_data_oldP,n_r_ic_max,    &
                    &           n_r_ic_max_old,lm_max_old,n_r_ic_maxL, &
                    &           lm2lmo,n_r_ic_maxMag,.true.,aj_ic,     &
                    &           dbdt_ic,djdt_ic,b_ic )
            end if
    
            deallocate( lm2lmo )
            deallocate( wo,zo,po,so )
         else if ( l_cond_ic ) then
            !----- No inner core fields provided by start_file, we thus assume that
            !      simple the inner core field decays like r**(l+1) from
            !      the ICB to r=0:
            write(*,'(/,'' ! USING POTENTIAL IC fields!'')')
    
            do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
               lmStart=lmStartB(nLMB)
               lmStop =lmStopB(nLMB)
               do lm=lmStart,lmStop
                  do nR=1,n_r_ic_max
                     b_ic(lm,nR)   =b(lm,n_r_CMB)
                     aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                     dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                     djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
                  end do
               end do
            end do
         end if
      end if
    
      !-- Lorentz-torques:
      !   NOTE: If lMagMem=.false. the memory required to read
      !         magnetic field is not available. The code therefore
      !         cannot read lorentz torques and rotations that
      !         are stored after the magnetic fields.
      !         In this case I set the lorentz torques to zero and
      !         calculate the rotation from the speed at the
      !         boundaries in the case of no slip conditions.
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
      if ( inform == 3 .and. l_mag_old .and. lMagMem == 1 ) then
         read(n_start_file,iostat=ioerr) lorentz_torque_ic, lorentz_torque_ma
         if( ioerr/=0 ) then
            write(*,*) '! Could not read last line in input file!'
            write(*,*) '! Data missing or wrong format!'
            write(*,*) '! Change inform accordingly!'
            stop
         end if
      else if ( inform >= 4 .and. inform <= 6 .and. lMagMem == 1 )then
         read(n_start_file,iostat=ioerr) lorentz_torque_ic, &
              lorentz_torque_ma,omega_ic,omega_ma
         if( ioerr/=0 ) then
            write(*,*) '! Could not read last line in input file!'
            write(*,*) '! Data missing or wrong format!'
            write(*,*) '! Change inform accordingly!'
            stop
         end if
      else if ( inform == 7 .or. inform == 8 ) then
         read(n_start_file,iostat=ioerr) lorentz_torque_ic, &
              lorentz_torque_ma, &
              omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
              omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
              omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
              omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2
         if( ioerr/=0 ) then
            write(*,*) '! Could not read last line in input file!'
            write(*,*) '! Data missing or wrong format!'
            write(*,*) '! Change inform accordingly!'
            stop
         end if
      else if ( inform > 8 ) then
         read(n_start_file,iostat=ioerr) lorentz_torque_ic, &
              lorentz_torque_ma,                       &
              omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
              omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
              omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
              omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2, &
              dt_new
         if( ioerr/=0 ) then
            write(*,*) '! Could not read last line in input file!'
            write(*,*) '! Data missing or wrong format!'
            write(*,*) '! Change inform accordingly!'
            stop
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
         if ( omega_ic1Old /= omega_ic1 )                     &
              write(*,*) '! New IC rotation rate 1 (old/new):', &
              omega_ic1Old,omega_ic1
         if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                    &
              write(*,*) '! New IC rotation osz. rate 1 (old/new):', &
              omegaOsz_ic1Old,omegaOsz_ic1
         if ( omega_ic2Old /= omega_ic2 )                     &
              write(*,*) '! New IC rotation rate 2 (old/new):', &
              omega_ic2Old,omega_ic2
         if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                    &
              write(*,*) '! New IC rotation osz. rate 2 (old/new):', &
              omegaOsz_ic2Old,omegaOsz_ic2
      end if
      if ( l_SRMA ) then
         if ( omega_ma1Old /= omega_ma1 )                     &
              write(*,*) '! New MA rotation rate 1 (old/new):', &
              omega_ma1Old,omega_ma1
         if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                    &
              write(*,*) '! New MA rotation osz. rate 1 (old/new):', &
              omegaOsz_ma1Old,omegaOsz_ma1
         if ( omega_ma2Old /= omega_ma2 )                     &
              write(*,*) '! New MA rotation rate 2 (old/new):', &
              omega_ma2Old,omega_ma2
         if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                    &
              write(*,*) '! New MA rotation osz. rate 2 (old/new):', &
              omegaOsz_ma2Old,omegaOsz_ma2
      end if
    
    
      !----- Set IC and mantle rotation rates:
      !      Following cases are covered:
      !       1) Prescribed inner-core rotation omega_ic_pre
      !       2) Rotation has been read above ( inform >= 4)
      !       3) Rotation calculated from flow field z(l=1,m=0)
      !       4) No rotation
      !       5) Flow driven by prescribed inner core rotation
      !       l_SRIC=.true. (spherical Couette case)
      l1m0=lm2(1,0)
      if ( l_rot_ic ) then
         if ( l_SRIC .or. omega_ic1 /= 0.0_cp ) then
            if ( tShift_ic1 == 0.0_cp ) tShift_ic1=tOmega_ic1-time
            if ( tShift_ic2 == 0.0_cp ) tShift_ic2=tOmega_ic2-time
            tOmega_ic1=time+tShift_ic1
            tOmega_ic2=time+tShift_ic2
            omega_ic=omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
                     omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
            write(*,*)
            write(*,*) '! I use prescribed inner core rotation rate:'
            write(*,*) '! omega_ic=',omega_ic
            if ( kbotv == 2 ) &
                 z(l1m0,n_r_icb)=cmplx(omega_ic/c_z10_omega_ic,0.0_cp,kind=cp)
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
            write(*,*)
            write(*,*) '! I use prescribed mantle rotation rate:'
            write(*,*) '! omega_ma =',omega_ma
            write(*,*) '! omega_ma1=',omega_ma1
            if ( ktopv == 2 ) &
                 z(l1m0,n_r_cmb)=cmplx(omega_ma/c_z10_omega_ma,0.0_cp,kind=cp)
         else if ( inform >= 7 ) then
            omega_ma=omega_ma1Old
         end if
      else
         omega_ma=0.0_cp
      end if
    
      close(n_start_file)
 
   end subroutine readStartFields
!------------------------------------------------------------------------------
#ifdef WITH_HDF5
   subroutine readHdf5_serial(w,dwdt,z,dzdt,p,dpdt,s,dsdt,xi,dxidt, &
              &               b,dbdt,aj,djdt,b_ic,dbdt_ic,aj_ic,    &
              &               djdt_ic,omega_ic,omega_ma,            &
              &               lorentz_torque_ic,lorentz_torque_ma,  &
              &               time,dt_old,dt_new)

      use hdf5
      use hdf5Helpers, only: readHdf5_attribute

      !--- Output variables
      real(cp),    intent(out) :: time,dt_old,dt_new
      real(cp),    intent(out) :: omega_ic,omega_ma
      real(cp),    intent(out) :: lorentz_torque_ic,lorentz_torque_ma
      complex(cp), intent(out) :: w(lm_max,n_r_max),z(lm_max,n_r_max)
      complex(cp), intent(out) :: s(lm_max,n_r_max),p(lm_max,n_r_max)
      complex(cp), intent(out) :: xi(lm_max,n_r_max)
      complex(cp), intent(out) :: dwdt(lm_max,n_r_max),dzdt(lm_max,n_r_max)
      complex(cp), intent(out) :: dsdt(lm_max,n_r_max),dpdt(lm_max,n_r_max)
      complex(cp), intent(out) :: dxidt(lm_max,n_r_max)
      complex(cp), intent(out) :: b(lm_maxMag,n_r_maxMag),aj(lm_maxMag,n_r_maxMag)
      complex(cp), intent(out) :: dbdt(lm_maxMag,n_r_maxMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,n_r_maxMag)
      complex(cp), intent(out) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: dbdt_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: djdt_ic(lm_maxMag,n_r_ic_maxMag)

      !--- Local variables
      integer :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
      integer :: l_max_old,n_r_max_old
      integer :: n_r_ic_max_old,n_data_oldP,n_r_maxL,n_r_ic_maxL
      integer :: l1m0,nLMB,lm,lmStart,lmStop,nR,lm_max_old
      logical :: l_mag_old

      real(cp) :: fr
      real(cp) :: pr_old,ra_old,pm_old
      real(cp) :: ek_old,radratio_old
      real(cp) :: sigma_ratio_old
      real(cp) :: omega_ic1Old,omegaOsz_ic1Old
      real(cp) :: omega_ic2Old,omegaOsz_ic2Old
      real(cp) :: omega_ma1Old,omegaOsz_ma1Old
      real(cp) :: omega_ma2Old,omegaOsz_ma2Old

      complex(cp), allocatable, target :: so(:),wo(:),zo(:),po(:),xio(:)

      type(C_PTR) :: f_ptr

      integer, allocatable :: lm2lmo(:)

      !--- HDF5 file
      integer(HID_T) :: file_id

      !--- HDF5 Attributes
      integer(HID_T) :: attr_id
      logical :: attr_exists, link_exists

      !--- HDF5 Groups
      character(len=12) :: grpname,grpname1     ! Group names
      integer(HID_T) :: grp_id,grp1_id

      !--- HDF5 Datasets
      integer(HID_T) :: dset_id

      !--- HDF5 Type
      integer(HSIZE_T) :: re_size,im_size,complex_t_size,offset
      integer(HID_T) :: type_id

      integer     ::   i,error

      inform = 12

      ! Initialize FORTRAN interface.
      call h5open_f(error)

      ! Create a new file using default properties.
      call h5fopen_f(start_file, H5F_ACC_RDONLY_F, file_id, error)

      ! Open group for control parameters and read attributes
      grpname = '/Params'
      call h5gopen_f(file_id, grpname, grp_id, error)
      call readHdf5_attribute(grp_id,'Ek',ek_old)
      call readHdf5_attribute(grp_id,'Ra',ra_old)
      call readHdf5_attribute(grp_id,'Pr',pr_old)
      call readHdf5_attribute(grp_id,'Prmag',pm_old)
      call readHdf5_attribute(grp_id,'Radratio',radratio_old)
      call readHdf5_attribute(grp_id,'Time',time)
      call readHdf5_attribute(grp_id,'dt',dt_old)
      call readHdf5_attribute(grp_id,'sigma_ratio',sigma_ratio_old)

      call h5gclose_f(grp_id, error)

      !---- Compare parameters:
      if ( ra /= ra_old ) &
         write(*,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') ra_old,ra
      if ( ek /= ek_old ) &
         write(*,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') ek_old,ek
      if ( pr /= pr_old ) &
         write(*,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') pr_old,pr
      if ( prmag /= pm_old )                                        &
         write(*,'(/,'' ! New mag Pr.number (old/new):'',2ES16.6)') &
         pm_old,prmag
      if ( radratio /= radratio_old )                                  &
         write(*,'(/,'' ! New mag aspect ratio (old/new):'',2ES16.6)') &
         radratio_old,radratio
      if ( sigma_ratio /= sigma_ratio_old )                           &
         write(*,'(/,'' ! New mag cond. ratio (old/new):'',2ES16.6)') &
         sigma_ratio_old,sigma_ratio

      ! Create a HF compound type  to store Fortran complex
      call h5tget_size_f(H5T_NATIVE_DOUBLE,re_size,error)
      call h5tget_size_f(H5T_NATIVE_DOUBLE,im_size,error)
      complex_t_size = re_size+im_size
      call h5tcreate_f(H5T_COMPOUND_F,complex_t_size,type_id,error)
      offset = 0
      call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, error)
      offset = offset + re_size
      call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, error)

      grpname = '/Fields'
      grpname1 = '/dtFields'

      call h5gopen_f(file_id, grpname, grp_id, error)
      call h5gopen_f(file_id, grpname1, grp1_id, error)

      call readHdf5_attribute(grp_id,'n_phi_tot',n_phi_tot_old)
      call readHdf5_attribute(grp_id,'minc',minc_old)
      call readHdf5_attribute(grp_id,'n_r_ic_max',n_r_ic_max_old)
      call readHdf5_attribute(grp_id,'n_r_max',n_r_max_old)
      call readHdf5_attribute(grp_id,'n_theta_max',n_theta_max_old)
      call readHdf5_attribute(grp_id,'nalias',nalias_old)

      l_max_old=nalias_old*n_phi_tot_old/60
      l_mag_old=.false.
      if ( pm_old /= 0.0_cp ) l_mag_old= .true.

      if ( n_phi_tot_old /= n_phi_tot) &
         write(*,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
      if ( nalias_old /= nalias) &
         write(*,*) '! New nalias (old,new)   :',nalias_old,nalias
      if ( l_max_old /= l_max ) &
         write(*,*) '! New l_max (old,new)    :',l_max_old,l_max

      allocate( lm2lmo(lm_max) )

      call getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
           &         m_max,minc,minc_old,inform,lm_max,   &
           &         lm_max_old,n_data_oldP,lm2lmo)
      n_data_oldP=lm_max*(n_r_max+1)

      allocate( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP),so(n_data_oldP) )
      bytes_allocated = bytes_allocated + 4*n_data_oldP*SIZEOF_DEF_COMPLEX

      call h5oexists_by_name_f(file_id, '/Fields/entropy', link_exists, error)
      if ( link_exists ) then
         call h5dopen_f(grp_id, 'entropy', dset_id, error)
         f_ptr=C_LOC(so)
         call h5dread_f(dset_id, type_id, f_ptr, error)
         call h5dclose_f(dset_id, error)
      end if
      call h5dopen_f(grp_id, 'w_pol', dset_id, error)
      f_ptr=C_LOC(wo)
      call h5dread_f(dset_id, type_id, f_ptr, error)
      call h5dclose_f(dset_id, error)

      call h5dopen_f(grp_id, 'z_tor', dset_id, error)
      f_ptr=C_LOC(zo)
      call h5dread_f(dset_id, type_id, f_ptr, error)
      call h5dclose_f(dset_id, error)

      call h5dopen_f(grp_id, 'pressure', dset_id, error)
      f_ptr=C_LOC(po)
      call h5dread_f(dset_id, type_id, f_ptr, error)
      call h5dclose_f(dset_id, error)

      n_r_maxL = max(n_r_max,n_r_max_old)
      call mapDataHydro( wo,zo,po,so,xio,n_data_oldP,lm2lmo,  &
           &             n_r_max_old,lm_max_old,n_r_maxL,     &
           &             .false.,.false.,.false.,.false.,     &
           &             .false.,w,z,p,s,xi )

      call h5oexists_by_name_f(file_id, '/dtFields/dsdtLast', link_exists, error)
      if ( link_exists ) then
         call h5dopen_f(grp1_id, 'dsdtLast', dset_id, error)
         f_ptr=C_LOC(so)
         call h5dread_f(dset_id, type_id, f_ptr, error)
         call h5dclose_f(dset_id, error)
      end if
      call h5dopen_f(grp1_id, 'dwdtLast', dset_id, error)
      f_ptr=C_LOC(wo)
      call h5dread_f(dset_id, type_id, f_ptr, error)
      call h5dclose_f(dset_id, error)

      call h5dopen_f(grp1_id, 'dzdtLast', dset_id, error)
      f_ptr=C_LOC(zo)
      call h5dread_f(dset_id, type_id, f_ptr, error)
      call h5dclose_f(dset_id, error)

      call h5dopen_f(grp1_id, 'dpdtLast', dset_id, error)
      f_ptr=C_LOC(po)
      call h5dread_f(dset_id, type_id, f_ptr, error)
      call h5dclose_f(dset_id, error)

      call mapDataHydro( wo,zo,po,so,xio,n_data_oldP,lm2lmo,     &
           &             n_r_max_old,lm_max_old,n_r_maxL,.true., &
           &            .true.,.true.,.true.,.true.,dwdt,dzdt,   &
           &            dpdt,dsdt,dxidt )

      if ( l_mag_old ) then
         call h5dopen_f(grp_id, 'b_pol', dset_id, error)
         f_ptr=C_LOC(so)
         call h5dread_f(dset_id, type_id, f_ptr, error)
         call h5dclose_f(dset_id, error)

         call h5dopen_f(grp_id, 'aj_tor', dset_id, error)
         f_ptr=C_LOC(wo)
         call h5dread_f(dset_id, type_id, f_ptr, error)
         call h5dclose_f(dset_id, error)

         call h5dopen_f(grp1_id, 'dbdtLast', dset_id, error)
         f_ptr=C_LOC(zo)
         call h5dread_f(dset_id, type_id, f_ptr, error)
         call h5dclose_f(dset_id, error)

         call h5dopen_f(grp1_id, 'djdtLast', dset_id, error)
         f_ptr=C_LOC(po)
         call h5dread_f(dset_id, type_id, f_ptr, error)
         call h5dclose_f(dset_id, error)

        call mapDataMag( wo,zo,po,so,n_data_oldP,n_r_max,n_r_max_old, &
             &              lm_max_old,n_r_maxL,lm2lmo,n_r_maxMag,    &
             &              .false.,aj,dbdt,djdt,b )
      else
        write(*,*) '! No magnetic data in input file!'
      end if

      !-- If mapping towards reduced symmetry, add thermal perturbation in
      !   mode (l,m)=(minc,minc) if parameter tipdipole /= 0
      if ( l_heat .and.  minc < minc_old .and. tipdipole > 0.0_cp ) then
         do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
            lmStart=lmStartB(nLMB)
            lmStop =lmStopB(nLMB)
            lm=l_max+2
            if ( lmStart <= lm .and. lmStop >= lm ) then
               do nR=1,n_r_max+1
                  fr=sin(pi*(r(nR)-r(n_r_max)))
                  s(lm,nR)=tipdipole*fr
               end do
            end if
         end do
      end if

      !-- If starting from data file with longitudinal symmetry, add
      !   weak non-axisymmetric dipole component if tipdipole /= 0
      if ( ( l_mag .or. l_mag_LF )                                    &
          &       .and. minc==1 .and. minc_old/=1 .and.               &
          &       tipdipole>0.0_cp .and. l_mag_old ) then
         do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
            lmStart=lmStartB(nLMB)
            lmStop =lmStopB(nLMB)
            lm=l_max+2
            if ( lmStart<=lm .and. lmStop>=lm ) then
               do nR=1,n_r_max+1
                  b(lm,nR)=tipdipole
               end do
            end if
         end do
      end if

      ! deallocation of the local arrays
      deallocate( lm2lmo )
      deallocate( wo,zo,po,so )
      !bytes_allocated = bytes_allocated - 4*n_data_oldP*SIZEOF_DEF_COMPLEX

      if ( l_mag_old ) then
         if ( sigma_ratio_old /= 0.0_cp ) then
            allocate( lm2lmo(lm_max) )
            call getLm2lmO(n_r_ic_max,n_r_ic_max_old,l_max,l_max_old, &
                 &         m_max,minc,minc_old,inform,lm_max,         &
                 &         lm_max_old,n_data_oldP,lm2lmo)

            n_r_ic_maxL = max(n_r_ic_max,n_r_ic_max_old)
            allocate( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP), &
                      so(n_data_oldP) )

            call h5dopen_f(grp_id, 'b_ic_pol', dset_id, error)
            f_ptr=C_LOC(so)
            call h5dread_f(dset_id, type_id, f_ptr, error)
            call h5dclose_f(dset_id, error)

            call h5dopen_f(grp_id, 'aj_ic_tor', dset_id, error)
            f_ptr=C_LOC(wo)
            call h5dread_f(dset_id, type_id, f_ptr, error)
            call h5dclose_f(dset_id, error)

            call h5dopen_f(grp1_id, 'dbdt_icLast', dset_id, error)
            f_ptr=C_LOC(zo)
            call h5dread_f(dset_id, type_id, f_ptr, error)
            call h5dclose_f(dset_id, error)

            call h5dopen_f(grp1_id, 'djdt_icLast', dset_id, error)
            f_ptr=C_LOC(po)
            call h5dread_f(dset_id, type_id, f_ptr, error)
            call h5dclose_f(dset_id, error)

            call mapDataMag( wo,zo,po,so,n_data_oldP,n_r_ic_max,n_r_ic_max_old, &
                             lm_max_old,n_r_ic_maxL,lm2lmo,n_r_ic_maxMag,       &
                             .true.,aj_ic,dbdt_ic,djdt_ic,b_ic )

            deallocate( lm2lmo )
            deallocate( wo,zo,po,so )
         else if ( l_cond_ic ) then
            !----- No inner core fields provided by start_file, we thus assume that
            !      simple the inner core field decays like r**(l+1) from
            !      the ICB to r=0:
            write(*,'(/,'' ! USING POTENTIAL IC fields!'')')

            do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
               lmStart=lmStartB(nLMB)
               lmStop =lmStopB(nLMB)
               do lm=lmStart,lmStop
                  do nR=1,n_r_ic_max
                     b_ic(lm,nR)   =b(lm,n_r_CMB)
                     aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                     dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                     djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
                  end do
               end do
            end do
         end if
      end if

      call h5gclose_f(grp1_id, error)
      call h5gclose_f(grp_id, error)

      !-- Lorentz-torques:
      !   NOTE: If lMagMem=.false. the memory required to read
      !         magnetic field is not available. The code therefore
      !         cannot read lorentz torques and rotations that
      !         are stored after the magnetic fields.
      !         In this case I set the lorentz torques to zero and
      !         calculate the rotation from the speed at the
      !         boundaries in the case of no slip conditions.
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

      ! Open group for control parameters and read attributes
      grpname = '/Torque'
      call h5gopen_f(file_id, grpname, grp_id, error)

      call readHdf5_attribute(grp_id,'lorentz_torque_ic',lorentz_torque_ic)
      call readHdf5_attribute(grp_id,'lorentz_torque_ma',lorentz_torque_ma)
      call readHdf5_attribute(grp_id,'omega_ic1',omega_ic1Old)
      call readHdf5_attribute(grp_id,'omegaOsz_ic1',omegaOsz_ic1Old)
      call readHdf5_attribute(grp_id,'tOmega_ic1',tOmega_ic1)
      call readHdf5_attribute(grp_id,'omega_ic2',omega_ic2Old)
      call readHdf5_attribute(grp_id,'omegaOsz_ic2',omegaOsz_ic2Old)
      call readHdf5_attribute(grp_id,'tOmega_ic2',tOmega_ic2)
      call readHdf5_attribute(grp_id,'omega_ma1',omega_ma1Old)
      call readHdf5_attribute(grp_id,'omegaOsz_ma1',omegaOsz_ma1Old)
      call readHdf5_attribute(grp_id,'tOmega_ma1',tOmega_ma1)
      call readHdf5_attribute(grp_id,'omega_ma2',omega_ma2Old)
      call readHdf5_attribute(grp_id,'omegaOsz_ma2',omegaOsz_ma2Old)
      call readHdf5_attribute(grp_id,'tOmega_ma2',tOmega_ma2)
      call readHdf5_attribute(grp_id,'dtNew',dt_new)

      call h5gclose_f(grp_id, error)

      if ( l_SRIC ) then
         if ( omega_ic1Old /= omega_ic1 )                     &
            write(*,*) '! New IC rotation rate 1 (old/new):', &
            omega_ic1Old,omega_ic1
         if ( omegaOsz_ic1Old /= omegaOsz_ic1 )                    &
            write(*,*) '! New IC rotation osz. rate 1 (old/new):', &
            omegaOsz_ic1Old,omegaOsz_ic1
         if ( omega_ic2Old /= omega_ic2 )                     &
            write(*,*) '! New IC rotation rate 2 (old/new):', &
            omega_ic2Old,omega_ic2
         if ( omegaOsz_ic2Old /= omegaOsz_ic2 )                    &
            write(*,*) '! New IC rotation osz. rate 2 (old/new):', &
            omegaOsz_ic2Old,omegaOsz_ic2
      end if
      if ( l_SRMA ) then
         if ( omega_ma1Old /= omega_ma1 )                     &
            write(*,*) '! New MA rotation rate 1 (old/new):', &
            omega_ma1Old,omega_ma1
         if ( omegaOsz_ma1Old /= omegaOsz_ma1 )                    &
            write(*,*) '! New MA rotation osz. rate 1 (old/new):', &
            omegaOsz_ma1Old,omegaOsz_ma1
         if ( omega_ma2Old /= omega_ma2 )                     &
            write(*,*) '! New MA rotation rate 2 (old/new):', &
            omega_ma2Old,omega_ma2
         if ( omegaOsz_ma2Old /= omegaOsz_ma2 )                    &
            write(*,*) '! New MA rotation osz. rate 2 (old/new):', &
            omegaOsz_ma2Old,omegaOsz_ma2
       end if


      ! Close the file.
      call h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      call h5close_f(error)

      !----- Set IC and mantle rotation rates:
      l1m0=lm2(1,0)
      if ( l_rot_ic ) then
         if ( l_SRIC .or. omega_ic1 /= 0.0_cp ) then
            if ( tShift_ic1 == 0.0_cp ) tShift_ic1=tOmega_ic1-time
            if ( tShift_ic2 == 0.0_cp ) tShift_ic2=tOmega_ic2-time
            tOmega_ic1=time+tShift_ic1
            tOmega_ic2=time+tShift_ic2
            omega_ic=omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
                     omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
            write(*,*)
            write(*,*) '! I use prescribed inner core rotation rate:'
            write(*,*) '! omega_ic=',omega_ic
            if ( kbotv == 2 ) &
               z(l1m0,n_r_icb)=cmplx(omega_ic/c_z10_omega_ic,0.0_cp,kind=cp)
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
            write(*,*)
            write(*,*) '! I use prescribed mantle rotation rate:'
            write(*,*) '! omega_ma =',omega_ma
            write(*,*) '! omega_ma1=',omega_ma1
            if ( ktopv == 2 ) &
               z(l1m0,n_r_cmb)=cmplx(omega_ma/c_z10_omega_ma,0.0_cp,kind=cp)
         else if ( inform >= 7 ) then
            omega_ma=omega_ma1Old
         end if
      else
         omega_ma=0.0_cp
      end if

   end subroutine readHdf5_serial
#endif
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

      if ( l_max==l_max_old .and. minc==minc_old .and. n_r_max==n_r_max_old .and. m_max==m_max_old ) then

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
              &     write(6,'('' ! Warning: Incompatible old/new minc= '',2i3)')

         lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
              &     m_max_old*(m_max_old-minc_old)/(2*minc_old) +     &
              &     l_max_old-m_max_old+1

         n_data_old=lm_max_old*n_r_max_old
         if ( inform>2 ) then
            n_data_oldP=n_data_old
         else
            n_data_oldP=lm_max_old*(n_r_max_old+1)
         end if

         !-- Write info to STdoUT:
         write(*,'('' ! Old/New  l_max= '',2I4,''  m_max= '',2I4,     &
              &       ''  minc= '',2I3,''  lm_max= '',2I5/)')         &
              &           l_max_old,l_max,m_max_old,m_max,            &
              &           minc_old,minc,lm_max_old,lm_max
         if ( n_r_max_old /= n_r_max )                                &
              &        write(*,'('' ! Old/New n_r_max='',2i4)')       &
             &              n_r_max_old,n_r_max

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
   subroutine mapDataHydro( wo,zo,po,so,xio,n_data_oldP,lm2lmo, &
              &             n_r_max_old,lm_max_old,n_r_maxL,    &
              &             lBc1,lBc2,lBc3,lBc4,lBc5,w,z,p,s,xi )

      !--- Input variables
      integer,     intent(in) :: n_r_max_old,lm_max_old
      integer,     intent(in) :: n_r_maxL,n_data_oldP
      logical,     intent(in) :: lBc1,lBc2,lBc3,lBc4,lBc5
      integer,     intent(in) :: lm2lmo(lm_max)
      complex(cp), intent(in) :: wo(:),zo(:)
      complex(cp), intent(in) :: po(:),so(:)
      complex(cp), intent(in) :: xio(:)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_max,n_r_max),z(lm_max,n_r_max)
      complex(cp), intent(out) :: p(lm_max,n_r_max),s(lm_max,n_r_max)
      complex(cp), intent(out) :: xi(lm_max,n_r_max)

      !--- Local variables
      integer :: lm,lmo,n,nR,lmStart,lmStop,nLMB
      complex(cp),allocatable :: woR(:),zoR(:)
      complex(cp),allocatable :: poR(:),soR(:)
      complex(cp),allocatable :: xioR(:)

      !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
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
      write(*,"(A,I12)") "maximal allocated bytes in mapData are ",bytes_allocated

      !PERFON('mD_map')
      do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)

         !PRINT*,nLMB,lmStart,lmStop
         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_r_max /= n_r_max_old ) then
                  do nR=1,n_r_max_old  ! copy on help arrays
                     n=lmo+(nR-1)*lm_max_old
                     woR(nR)=wo(n)
                     zoR(nR)=zo(n)
                     poR(nR)=po(n)
                     if ( lreadS .and. l_heat ) soR(nR)=so(n)
                     if ( lreadXi .and. l_chemical_conv ) xioR(nR)=xio(n)
                  end do
                  call mapDataR(woR,n_r_max,n_r_max_old,n_r_maxL,lBc1,.false.)
                  call mapDataR(zoR,n_r_max,n_r_max_old,n_r_maxL,lBc2,.false.)
                  call mapDataR(poR,n_r_max,n_r_max_old,n_r_maxL,lBc3,.false.)
                  if ( lreadS .and. l_heat ) call mapDataR(soR,n_r_max, &
                                                  n_r_max_old,n_r_maxL, &
                                                  lBc4,.false.)
                  if ( lreadXi .and. l_chemical_conv ) call mapDataR(xioR, &
                                              n_r_max,n_r_max_old,n_r_maxL,&
                                              lBc5,.false.)
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
                     if ( lreadXi .and. l_chemical_conv ) xi(lm,nR)=scale_xi*xioR(nR)
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
                     if ( lreadXi .and. l_chemical_conv ) xi(lm,nR)=scale_xi*xio(n)
                  end do
               end if
            end if
         end do
      end do
      !PERFOFF
      !PRINT*,omp_get_thread_num(),": After nLMB loop"
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
   subroutine mapDataMag( wo,zo,po,so,n_data_oldP,n_rad_tot,n_r_max_old, &
                          lm_max_old,n_r_maxL,lm2lmo,dim1,l_IC,          & 
                          w,z,p,s )

      !--- Input variables
      integer,     intent(in) :: n_rad_tot,n_r_max_old,lm_max_old
      integer,     intent(in) :: n_r_maxL,n_data_oldP,dim1
      integer,     intent(in) :: lm2lmo(lm_max)
      logical,     intent(in) :: l_IC
      complex(cp), intent(in) :: wo(n_data_oldP),zo(n_data_oldP)
      complex(cp), intent(in) :: po(n_data_oldP),so(n_data_oldP)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_maxMag,dim1),z(lm_maxMag,dim1)
      complex(cp), intent(out) :: p(lm_maxMag,dim1),s(lm_maxMag,dim1)

      !--- Local variables
      integer :: lm,lmo,n,nR,lmStart,lmStop,nLMB
      complex(cp), allocatable :: woR(:),zoR(:),poR(:),soR(:)

      !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
      allocate( woR(n_r_maxL),zoR(n_r_maxL) )
      allocate( poR(n_r_maxL),soR(n_r_maxL) )
      bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DEF_COMPLEX
      write(*,"(A,I12)") "maximal allocated bytes in mapData are ",bytes_allocated

      !PERFON('mD_map')
      do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)
         lmStart=max(2,lmStart)
         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_rad_tot /= n_r_max_old ) then
                  do nR=1,n_r_max_old  ! copy on help arrays
                     n=lmo+(nR-1)*lm_max_old
                     woR(nR)=wo(n)
                     zoR(nR)=zo(n)
                     poR(nR)=po(n)
                     soR(nR)=so(n)
                  end do
                  call mapDataR(woR,dim1,n_r_max_old,n_r_maxL,.false.,l_IC)
                  call mapDataR(zoR,dim1,n_r_max_old,n_r_maxL,.true.,l_IC)
                  call mapDataR(poR,dim1,n_r_max_old,n_r_maxL,.true.,l_IC)
                  call mapDataR(soR,dim1,n_r_max_old,n_r_maxL,.false.,l_IC)
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
            end if
         end do
      end do
      !PERFOFF
      !PRINT*,omp_get_thread_num(),": After nLMB loop"
      deallocate(woR,zoR,poR,soR)
      bytes_allocated = bytes_allocated - 4*n_r_maxL*SIZEOF_DEF_COMPLEX

   end subroutine mapDataMag
!------------------------------------------------------------------------------
   subroutine mapDataR(dataR,n_rad_tot,n_r_max_old,n_r_maxL,lBc,l_IC)
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
      integer, intent(in) :: n_r_max_old
      integer, intent(in) :: n_r_maxL,n_rad_tot
      logical, intent(in) :: lBc,l_IC

      !--- Output variables
      complex(cp), intent(out) :: dataR(:)  ! old data 

      !-- Local variables
      integer :: nR, n_r_index_start
      type(costf_odd_t) :: chebt_oc_old
      complex(cp), allocatable :: work(:)
      real(cp) :: cheb_norm_old,scale

      allocate( work(n_r_maxL) )

      !----- Initialize transform to cheb space:
      call chebt_oc_old%initialize(n_r_max_old, 2*n_r_maxL+2,2*n_r_maxL+5)

      !-- Guess the boundary values, since they have not been stored:
      if ( .not. l_IC .and. lBc ) then
         dataR(1)=two*dataR(2)-dataR(3)
         dataR(n_r_max_old)=two*dataR(n_r_max_old-1)-dataR(n_r_max_old-2)
      end if

      !----- Transform old data to cheb space:
      call chebt_oc_old%costf1(dataR,work)

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
      else
         call chebt_oc%costf1(dataR,work)
         !----- Rescale :
         cheb_norm_old=sqrt(two/real(n_r_max_old-1,kind=cp))
         scale=cheb_norm_old/cheb_norm
      end if
      do nR=1,n_rad_tot
         dataR(nR)=scale*dataR(nR)
      end do

      call chebt_oc_old%finalize()
      deallocate( work )

   end subroutine mapDataR
!---------------------------------------------------------------------
end module readCheckPoints
