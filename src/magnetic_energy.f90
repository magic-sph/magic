module magnetic_energy

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use LMmapping, only: map_mlo, map_glbl_st
   use truncation, only: n_r_maxMag, n_r_ic_maxMag, n_r_max, n_r_ic_max, &
       &                 lm_max, minc, n_r_cmb, n_mloMag_loc, n_mlo_loc
   use radial_functions, only: r_icb, r_cmb, r_ic, dr_fac_ic, chebt_ic, &
       &                       sigma, orho1, r, or2, rscheme_oc
   use physical_parameters, only: LFfac, kbotb, ktopb
   use num_param, only: eScale, tScale
   use logic, only: l_cond_ic, l_mag, l_mag_LF, l_save_out, l_earth_likeness, &
       &            l_full_sphere
   use movie_data, only: movieDipColat, movieDipLon, movieDipStrength, &
       &                 movieDipStrengthGeo
   use output_data, only: tag, l_max_comp
   use constants, only: pi, zero, one, two, half, four, osq4pi
   use special, only: n_imp, rrMP
   use integration, only: rInt_R,rIntIC
   use useful, only: cc2real,cc22real
   use plms_theta, only: plm_theta
   use communications, only: gather_from_mlo_to_master, reduce_to_master, &
       &                     reduce_scalar, send_lm_pair_to_master

   implicit none

   private

   integer :: n_theta_max_comp
   integer :: n_phi_max_comp
   integer :: lm_max_comp
   real(cp), allocatable :: Plm_comp(:,:)

   real(cp), allocatable :: e_dipA(:)    ! Time-averaged dipole (l=1) energy
   real(cp), allocatable :: e_pA(:)      ! Time-averaged poloidal energy
   real(cp), allocatable :: e_p_asA(:)   ! Time-averaged axisymmetric poloidal energy
   real(cp), allocatable :: e_tA(:)      ! Time-averaged toroidal energy
   real(cp), allocatable :: e_t_asA(:)   ! Time-averaged axisymmetric toroidal energy
   complex(cp), allocatable :: bCMB(:)

   integer :: n_dipole_file, n_e_mag_ic_file, n_e_mag_oc_file
   integer :: n_compliance_file
   character(len=72) :: dipole_file, e_mag_ic_file, e_mag_oc_file
   character(len=72) :: earth_compliance_file

   public :: initialize_magnetic_energy, get_e_mag, finalize_magnetic_energy

contains

   subroutine initialize_magnetic_energy
      !
      ! Open diagnostic files and allocate memory
      !

      !-- Local variables
      integer :: lm, nTheta, l, m
      real(cp) :: theta
      real(cp), allocatable :: plma(:), dtheta_plma(:)

      allocate( e_dipA(n_r_max) )
      allocate( e_pA(n_r_max),e_p_asA(n_r_max) )
      allocate( e_tA(n_r_max),e_t_asA(n_r_max) )
      bytes_allocated = bytes_allocated+5*n_r_max*SIZEOF_DEF_REAL

      e_mag_ic_file='e_mag_ic.'//tag
      e_mag_oc_file='e_mag_oc.'//tag
      dipole_file  ='dipole.'//tag
      if ( l_master_rank .and. .not. l_save_out ) then
         open(newunit=n_e_mag_oc_file, file=e_mag_oc_file, status='new')
         if ( .not. l_full_sphere ) then
            open(newunit=n_e_mag_ic_file, file=e_mag_ic_file, status='new')
         end if
         open(newunit=n_dipole_file, file=dipole_file, status='new')
      end if

      if ( l_earth_likeness ) then
         if ( coord_r == 0 ) then
            earth_compliance_file='earth_like.'//tag
            if ( l_master_rank .and. .not. l_save_out ) then
               open(newunit=n_compliance_file, file=earth_compliance_file, &
               &    status='new')
            end if
            allocate( bCMB(lm_max) )
            bytes_allocated = bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX

            !-- Define a regularly spaced theta-grid
            !-- to estimate Br skewness
            n_theta_max_comp = 5*l_max_comp
            n_phi_max_comp = 2*n_theta_max_comp
            lm_max_comp=l_max_comp*(l_max_comp+1)/minc-&
            &           l_max_comp*(l_max_comp-minc)/(2*minc)+1
            !-- Plm_comp are values of associated Legendre functions for equidistant
            !-- grid points in theta  (theta = pi*(ith-0.5)/nth,  ith= 1  ... nth)
            allocate( Plm_comp(lm_max_comp,1:n_theta_max_comp) )

            !-- Temporary arrays
            allocate( plma(lm_max_comp), dtheta_plma(lm_max_comp) )
            do nTheta=1,n_theta_max_comp
               theta = pi*(nTheta-half)/real(n_theta_max_comp,cp)
               !-- Schmidt normalisation --!
               call plm_theta(theta, l_max_comp, l_max_comp, minc, &
                    &         plma, dtheta_plma, lm_max_comp, 1)
               lm = 1
               do m=0,l_max_comp,minc
                  do l=m,l_max_comp
                     Plm_comp(lm,nTheta)=plma(lm)
                     lm = lm+1
                  end do
               end do
            end do
            !
            deallocate( plma, dtheta_plma )
         else
            allocate( bCMB(1) )
         end if
      end if

   end subroutine initialize_magnetic_energy
!----------------------------------------------------------------------------
   subroutine finalize_magnetic_energy
      !
      ! Close file and deallocates global arrays
      !

      if ( l_earth_likeness ) then
         if ( coord_r == 0 ) deallocate( Plm_comp )
         deallocate( bCMB )
         if ( l_master_rank .and. .not. l_save_out ) close(n_compliance_file)
      end if
      deallocate( e_dipA, e_pA, e_p_asA, e_tA, e_t_asA)

      if ( l_master_rank .and. .not. l_save_out ) then
         close(n_e_mag_oc_file)
         if ( .not. l_full_sphere ) close(n_e_mag_ic_file)
         close(n_dipole_file)
      end if

   end subroutine finalize_magnetic_energy
!----------------------------------------------------------------------------
   subroutine get_e_mag(time,l_write,l_stop_time,n_e_sets,        &
              &         b,db,aj,b_ic,db_ic,aj_ic,                 &
              &         e_p,e_t,e_p_as,e_t_as,                    &
              &         e_p_ic,e_t_ic,e_p_as_ic,e_t_as_ic,        &
              &         e_p_os,e_p_as_os,e_cmb,Dip,DipCMB,        &
              &         elsAnel)
      !
      !  calculates magnetic energy  = 1/2 Integral(B^2 dV)
      !  integration in theta,phi by summation over harmonic coeffs.
      !  integration in r by Chebyshev integrals or Simpson rule
      !  depending whether FD or Cheb is used.
      !

      !-- Input of variables:
      integer,     intent(in) :: n_e_sets  ! Switch for time-average and to determine first time step
      real(cp),    intent(in) :: time      ! Current time
      logical,     intent(in) :: l_write   ! Switch to write output
      logical,     intent(in) :: l_stop_time ! Indicates when last time step of the run is reached for radial output
      complex(cp), intent(in) :: b(n_mloMag_loc,n_r_maxMag)        ! Array containing magnetic field poloidal potential
      complex(cp), intent(in) :: db(n_mloMag_loc,n_r_maxMag)       ! Array containing radial derivative of b
      complex(cp), intent(in) :: aj(n_mloMag_loc,n_r_maxMag)       ! Array containing magnetic field toroidal potential
      complex(cp), intent(in) :: b_ic(n_mloMag_loc,n_r_ic_maxMag)  ! Array containing IC magnetic field poloidal potential
      complex(cp), intent(in) :: db_ic(n_mloMag_loc,n_r_ic_maxMag) ! Array containing radial derivative of IC b
      complex(cp), intent(in) :: aj_ic(n_mloMag_loc,n_r_ic_maxMag) ! Array containing IC magnetic field toroidal potential

      !-- Output variables:
      real(cp), intent(out) :: e_p          ! Volume averaged poloidal magnetic energy
      real(cp), intent(out) :: e_t          ! Volume averaged toroidal magnetic energy
      real(cp), intent(out) :: e_p_as       ! Volume averaged axisymmetric poloidal magnetic energy
      real(cp), intent(out) :: e_t_as       ! Volume averaged axisymmetric toroidal magnetic energy
      real(cp), intent(out) :: e_p_ic       ! IC poloidal magnetic energy
      real(cp), intent(out) :: e_t_ic       ! IC toroidal magnetic energy
      real(cp), intent(out) :: e_p_as_ic    ! IC axisymmetric poloidal magnetic energy
      real(cp), intent(out) :: e_t_as_ic    ! IC axisymmetric toroidal magnetic energy
      real(cp), intent(out) :: e_p_os       ! Outside poloidal magnetic energy
      real(cp), intent(out) :: e_p_as_os    ! Outside axisymmetric poloidal magnetic energy
      real(cp), intent(out) :: e_cmb        ! Magnetic energy at the CMB
      real(cp), intent(out) :: Dip          ! Relative magnetic energy of axial dipole
      real(cp), intent(out) :: DipCMB       ! Relative magnetic energy of axial dipole at the CMB
      real(cp), intent(out) :: elsAnel      ! Radially averaged Elsasser number

      !-- local:
      integer :: nR,lm,l,m,l_geo

      real(cp) :: e_p_r(n_r_max), e_p_r_global(n_r_max)
      real(cp) :: e_t_r(n_r_max), e_t_r_global(n_r_max)
      real(cp) :: els_r(n_r_max), els_r_global(n_r_max)
      real(cp) :: e_p_as_r(n_r_max), e_p_as_r_global(n_r_max)
      real(cp) :: e_t_as_r(n_r_max), e_t_as_r_global(n_r_max)
      real(cp) :: e_p_es_r(n_r_max), e_p_es_r_global(n_r_max)
      real(cp) :: e_t_es_r(n_r_max), e_t_es_r_global(n_r_max)
      real(cp) :: e_p_eas_r(n_r_max), e_p_eas_r_global(n_r_max)
      real(cp) :: e_t_eas_r(n_r_max), e_t_eas_r_global(n_r_max)
      real(cp) :: e_dipole_r(n_r_max), e_dipole_r_global(n_r_max)
      real(cp) :: e_dipole_ax_r(n_r_max), e_dipole_ax_r_global(n_r_max)

      real(cp) :: e_p_ic_r(n_r_ic_max), e_p_ic_r_global(n_r_ic_max)
      real(cp) :: e_t_ic_r(n_r_ic_max), e_t_ic_r_global(n_r_ic_max)
      real(cp) :: e_p_as_ic_r(n_r_ic_max), e_p_as_ic_r_global(n_r_ic_max)
      real(cp) :: e_t_as_ic_r(n_r_ic_max), e_t_as_ic_r_global(n_r_ic_max)

      real(cp) :: e_geo,e_es_geo,e_as_geo,e_eas_geo
      real(cp) :: e_geo_global,e_es_geo_global,e_as_geo_global,e_eas_geo_global
      real(cp) :: e_p_ic_global, e_p_as_ic_global, e_p_os_global, e_p_as_os_global
      real(cp) :: e_p_e,e_p_as_e,e_p_e_global,e_p_as_e_global

      real(cp) :: r_ratio,fac,e_p_temp,e_t_temp,dLh
      real(cp) :: e_dipole, e_dipole_ax, e_dipole_ax_cmb
      real(cp) :: e_dipole_e, e_dipole_e_global
      real(cp) :: e_p_e_ratio,O_r_icb_E_2,rad
      real(cp) :: e_p_es,e_t_es,e_es_cmb,e_as_cmb
      real(cp) :: e_p_eas,e_t_eas,e_eas_cmb

      real(cp) :: ad, nad, sym, asym, zon, nzon, br_skew
      real(cp) :: fluxConcentration
      real(cp) :: ad_global, nad_global, sym_global, asym_global
      real(cp) :: zon_global, nzon_global

      real(cp) :: e_dip_cmb,eTot,eDR,theta_dip,phi_dip

      complex(cp) :: r_dr_b,b10,b11

      !-- time averaging of e(r):
      character(len=80) :: filename
      real(cp) :: dt, osurf
      real(cp), save :: timeLast,timeTot
      integer :: fileHandle

      l_geo=11   ! max degree for geomagnetic field seen on Earth

      e_p      =0.0_cp
      e_t      =0.0_cp
      e_p_as   =0.0_cp
      e_t_as   =0.0_cp
      e_p_ic   =0.0_cp
      e_t_ic   =0.0_cp
      e_p_as_ic=0.0_cp
      e_t_as_ic=0.0_cp
      e_p_os   =0.0_cp
      e_p_as_os=0.0_cp
      e_geo    =0.0_cp
      e_es_geo =0.0_cp
      e_as_geo =0.0_cp
      e_eas_geo=0.0_cp
      Dip      =0.0_cp
      DipCMB   =0.0_cp
      ad       =0.0_cp
      nad      =0.0_cp
      sym      =0.0_cp
      asym     =0.0_cp
      zon      =0.0_cp
      nzon     =0.0_cp

      if ( .not.( l_mag .or. l_mag_LF ) ) return

      do nR=1,n_r_max

         e_p_r(nR)     =0.0_cp
         e_t_r(nR)     =0.0_cp
         e_p_as_r(nR)  =0.0_cp
         e_t_as_r(nR)  =0.0_cp
         e_p_es_r(nR)  =0.0_cp
         e_t_es_r(nR)  =0.0_cp
         e_p_eas_r(nR) =0.0_cp
         e_t_eas_r(nR) =0.0_cp
         e_dipole_r(nR)=0.0_cp
         e_dipole_ax_r(nR)=0.0_cp

         do lm=1,n_mlo_loc
            l=map_mlo%i2l(lm)
            m=map_mlo%i2m(lm)
            if ( l == 0 ) cycle
            dLh = real(l*(l+1),cp)

            e_p_temp= dLh*( dLh*or2(nR)*cc2real( b(lm,nR),m) &
            &                         + cc2real(db(lm,nR),m) )
            e_t_temp= dLh * cc2real(aj(lm,nR),m)

            if ( m == 0 ) then  ! axisymmetric part
               e_p_as_r(nR)=e_p_as_r(nR) + e_p_temp
               e_t_as_r(nR)=e_t_as_r(nR) + e_t_temp
               if ( mod(l,2) == 1 ) then
                  e_p_eas_r(nR)=e_p_eas_r(nR)+e_p_temp
               else
                  e_t_eas_r(nR)=e_t_eas_r(nR)+e_t_temp
               end if
            else
               e_p_r(nR)=e_p_r(nR) + e_p_temp
               e_t_r(nR)=e_t_r(nR) + e_t_temp
            end if
            if ( mod(l+m,2) == 1 ) then
               e_p_es_r(nR)=e_p_es_r(nR) + e_p_temp
            else
               e_t_es_r(nR)=e_t_es_r(nR) + e_t_temp
            end if
            if ( l <= l_geo .and. nR == n_r_cmb ) then
               e_geo    =e_geo    +e_p_temp
               if ( mod(l+m,2) == 1 ) e_es_geo = e_es_geo + e_p_temp
               if ( m == 0 )          e_as_geo = e_as_geo + e_p_temp
               if ( mod(l+m,2) == 1 .and. m == 0 )                    &
               &                      e_eas_geo=e_eas_geo+e_p_temp
            end if

            if ( l == 1 .and. m == 0 ) then
               e_dipole_ax_r(nR)=e_p_temp
               if ( nR == n_r_cmb ) ad = e_dipole_ax_r(nR)
            end if
            if ( l == 1 ) e_dipole_r(nR)=e_dipole_r(nR)+e_p_temp

            if ( l_earth_likeness ) then
               if ( l <= l_max_comp .and. nR == n_r_cmb ) then
                  if ( l /= 1 .or. m /= 0 ) then
                     nad = nad + e_p_temp
                  end if
               end if

               if ( nR ==  n_r_cmb .and. l >= 2 .and. l <= l_max_comp ) then
                  if ( mod(l+m,2) == 1 ) then
                     sym = sym + e_p_temp
                  else
                     asym = asym+e_p_temp
                  end if
                  if ( m == 0 ) then
                     zon = zon + e_p_temp
                  else
                     nzon = nzon + e_p_temp
                  end if
               end if
            end if

         end do    ! do loop over lms in block

         e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
         e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)

         ! In anelastic models it is also interesting to have Lambda
         els_r(nR)=(e_p_r(nR)+e_t_r(nR))*orho1(nR)*sigma(nR)

      end do    ! radial grid points

      call reduce_to_master(e_p_r, e_p_r_global, 0)
      call reduce_to_master(e_t_r, e_t_r_global, 0)
      call reduce_to_master(e_p_as_r, e_p_as_r_global, 0)
      call reduce_to_master(e_t_as_r, e_t_as_r_global, 0)
      call reduce_to_master(e_p_es_r, e_p_es_r_global, 0)
      call reduce_to_master(e_t_es_r, e_t_es_r_global, 0)
      call reduce_to_master(e_p_eas_r, e_p_eas_r_global, 0)
      call reduce_to_master(e_t_eas_r, e_t_eas_r_global, 0)
      call reduce_to_master(e_dipole_ax_r, e_dipole_ax_r_global, 0)
      call reduce_to_master(e_dipole_r, e_dipole_r_global, 0)
      call reduce_to_master(els_r, els_r_global, 0)

      ! reduce some scalars
      call reduce_scalar(e_geo, e_geo_global, 0)
      call reduce_scalar(e_es_geo, e_es_geo_global, 0)
      call reduce_scalar(e_as_geo, e_as_geo_global, 0)
      call reduce_scalar(e_eas_geo, e_eas_geo_global, 0)
      if ( l_earth_likeness ) then
         call reduce_scalar(ad, ad_global, 0)
         call reduce_scalar(nad, nad_global, 0)
         call reduce_scalar(sym, sym_global, 0)
         call reduce_scalar(asym, asym_global, 0)
         call reduce_scalar(zon, zon_global, 0)
         call reduce_scalar(nzon, nzon_global, 0)
      end if

      if ( l_earth_likeness ) call gather_from_mlo_to_master(b(:,n_r_cmb), bCMB)

      if ( l_master_rank ) then
         !-- Get Values at CMB:
         e_cmb          =e_p_r_global(n_r_cmb)+e_t_r_global(n_r_cmb)
         e_dip_cmb      =e_dipole_r_global(n_r_cmb)
         e_dipole_ax_cmb=e_dipole_ax_r_global(n_r_cmb)
         e_es_cmb       =e_p_es_r_global(n_r_cmb)
         e_as_cmb       =e_p_as_r_global(n_r_cmb)
         e_eas_cmb      =e_p_eas_r_global(n_r_cmb)

         ! NOTE: n_e_sets=0 prevents averaging
         if ( n_e_sets == 1 ) then
            timeTot=one
            do nR=1,n_r_max
               e_dipA(nR) =e_dipole_r_global(nR)
               e_pA(nR)   =e_p_r_global(nR)
               e_p_asA(nR)=e_p_as_r_global(nR)
               e_tA(nR)   =e_t_r_global(nR)
               e_t_asA(nR)=e_t_as_r_global(nR)
            end do
         else if ( n_e_sets == 2 ) then
            dt=time-timeLast
            timeTot=two*dt
            do nR=1,n_r_max
               e_dipA(nR) =dt*(e_dipA(nR) +e_dipole_r_global(nR))
               e_pA(nR)   =dt*(e_pA(nR)   +e_p_r_global(nR)     )
               e_p_asA(nR)=dt*(e_p_asA(nR)+e_p_as_r_global(nR)  )
               e_tA(nR)   =dt*(e_tA(nR)   +e_t_r_global(nR)     )
               e_t_asA(nR)=dt*(e_t_asA(nR)+e_t_as_r_global(nR)  )
            end do
         else
            dt=time-timeLast
            timeTot=timeTot+dt
            do nR=1,n_r_max
               e_dipA(nR) =e_dipA(nR) +dt*e_dipole_r_global(nR)
               e_pA(nR)   =e_pA(nR)   +dt*e_p_r_global(nR)
               e_p_asA(nR)=e_p_asA(nR)+dt*e_p_as_r_global(nR)
               e_tA(nR)   =e_tA(nR)   +dt*e_t_r_global(nR)
               e_t_asA(nR)=e_t_asA(nR)+dt*e_t_as_r_global(nR)
            end do
         end if
         if ( ktopb == 1) then
            e_tA(1)   =0.0_cp
            e_t_asA(1)=0.0_cp
         endif
         if ( kbotb == 1 ) then
            e_tA(n_r_max)   =0.0_cp
            e_t_asA(n_r_max)=0.0_cp
         endif
         if ( l_stop_time ) then
            fac=half*LFfac*eScale
            filename='eMagR.'//tag
            open(newunit=fileHandle, file=filename, status='unknown')
            do nR=1,n_r_max
               eTot=e_pA(nR)+e_tA(nR)
               if ( e_dipA(nR)  <  1.e-6_cp*eTot ) then
                  eDR=0.0_cp
               else
                  eDR=e_dipA(nR)/eTot
               end if
               osurf=0.25_cp/pi*or2(nR)
               write(fileHandle,'(ES20.10,9ES15.7)') r(nR),         &
               &                    fac*e_pA(nR)/timetot,           &
               &                    fac*e_p_asA(nR)/timetot,        &
               &                    fac*e_tA(nR)/timetot,           &
               &                    fac*e_t_asA(nR)/timetot,        &
               &                    fac*e_pA(nR)/timetot*osurf,     &
               &                    fac*e_p_asA(nR)/timetot*osurf,  &
               &                    fac*e_tA(nR)/timetot*osurf,     &
               &                    fac*e_t_asA(nR)/timetot*osurf,  &
               &                    eDR
            end do
            close(fileHandle)
         end if
         timeLast=time

         !-- Radial integrals:
         e_p        =rInt_R(e_p_r_global,r,rscheme_oc)
         e_t        =rInt_R(e_t_r_global,r,rscheme_oc)
         e_p_as     =rInt_R(e_p_as_r_global,r,rscheme_oc)
         e_t_as     =rInt_R(e_t_as_r_global,r,rscheme_oc)
         e_p_es     =rInt_R(e_p_es_r_global,r,rscheme_oc)
         e_t_es     =rInt_R(e_t_es_r_global,r,rscheme_oc)
         e_p_eas    =rInt_R(e_p_eas_r_global,r,rscheme_oc)
         e_t_eas    =rInt_R(e_t_eas_r_global,r,rscheme_oc)
         e_dipole   =rInt_R(e_dipole_r_global,r,rscheme_oc)
         e_dipole_ax=rInt_R(e_dipole_ax_r_global,r,rscheme_oc)
         elsAnel    =rInt_R(els_r_global,r,rscheme_oc)

         fac=half*LFfac*eScale
         e_p        =fac*e_p
         e_t        =fac*e_t
         e_p_as     =fac*e_p_as
         e_t_as     =fac*e_t_as
         e_p_es     =fac*e_p_es
         e_t_es     =fac*e_t_es
         e_p_eas    =fac*e_p_eas
         e_t_eas    =fac*e_t_eas
         e_dipole   =fac*e_dipole
         e_dipole_ax=fac*e_dipole_ax
         ! Elsasser must not change with scalings
         !elsAnel    =eScale*elsAnel

         e_cmb          =fac*e_cmb
         e_dip_cmb      =fac*e_dip_cmb
         e_dipole_ax_cmb=fac*e_dipole_ax_cmb
         e_es_cmb       =fac*e_es_cmb
         e_as_cmb       =fac*e_as_cmb
         e_eas_cmb      =fac*e_eas_cmb
         e_geo          =fac*e_geo_global
         e_es_geo       =fac*e_es_geo_global
         e_as_geo       =fac*e_as_geo_global
         e_eas_geo      =fac*e_eas_geo_global
      end if

      !-- Inner core:

      if ( l_cond_ic ) then

         O_r_icb_E_2=one/(r(n_r_max)*r(n_r_max))

         do nR=1,n_r_ic_max

            r_ratio=r_ic(nR)/r_ic(1)

            e_p_ic_r(nR)   =0.0_cp
            e_t_ic_r(nR)   =0.0_cp
            e_p_as_ic_r(nR)=0.0_cp
            e_t_as_ic_r(nR)=0.0_cp

            do lm=1,n_mlo_loc
               l=map_mlo%i2l(lm)
               m=map_mlo%i2m(lm)
               if ( l == 0 ) cycle
               dLh = real(l*(l+1),cp)
               r_dr_b=r_ic(nR)*db_ic(lm,nR)

               e_p_temp=     dLh*O_r_icb_E_2*r_ratio**(2*l) * (               &
               &         real((l+1)*(2*l+1),cp)*cc2real(b_ic(lm,nR),m)     +  &
               &         real(2*(l+1),cp)*cc22real(b_ic(lm,nR),r_dr_b,m)   +  &
               &                                 cc2real(r_dr_b,m)            )
               e_t_temp=  dLh*r_ratio**(2*l+2)*cc2real(aj_ic(lm,nR),m)

               if ( m == 0 ) then  ! axisymmetric part
                  e_p_as_ic_r(nR)=e_p_as_ic_r(nR) + e_p_temp
                  e_t_as_ic_r(nR)=e_t_as_ic_r(nR) + e_t_temp
               else
                  e_p_ic_r(nR)   =e_p_ic_r(nR) + e_p_temp
                  e_t_ic_r(nR)   =e_t_ic_r(nR) + e_t_temp
               end if

            end do    ! do loop over lms in block

            e_p_ic_r(nR)=e_p_ic_r(nR)+e_p_as_ic_r(nR)
            e_t_ic_r(nR)=e_t_ic_r(nR)+e_t_as_ic_r(nR)

         end do    ! radial grid points

         ! reduce over the ranks
         call reduce_to_master(e_p_ic_r, e_p_ic_r_global, 0)
         call reduce_to_master(e_t_ic_r, e_t_ic_r_global, 0)
         call reduce_to_master(e_p_as_ic_r, e_p_as_ic_r_global, 0)
         call reduce_to_master(e_t_as_ic_r, e_t_as_ic_r_global, 0)

         if ( l_master_rank ) then
            e_p_ic   =rIntIC(e_p_ic_r_global,n_r_ic_max,dr_fac_ic,chebt_ic)
            e_t_ic   =rIntIC(e_t_ic_r_global,n_r_ic_max,dr_fac_ic,chebt_ic)
            e_p_as_ic=rIntIC(e_p_as_ic_r_global,n_r_ic_max,dr_fac_ic,chebt_ic)
            e_t_as_ic=rIntIC(e_t_as_ic_r_global,n_r_ic_max,dr_fac_ic,chebt_ic)
            fac=half*LFfac*eScale
            e_p_ic   =fac*e_p_ic
            e_t_ic   =fac*e_t_ic
            e_p_as_ic=fac*e_p_as_ic
            e_t_as_ic=fac*e_t_as_ic
         end if

      else if ( (.not. l_cond_ic) .and. (.not. l_full_sphere) ) then

         !do lm=2,lm_max
         do lm=1,n_mlo_loc
            l=map_mlo%i2l(lm)
            m=map_mlo%i2m(lm)
            if ( l == 0 ) cycle
            fac=real(l*(l+1)*(l+1),cp)
            e_p_temp=fac*cc2real(b(lm,n_r_max),m)
            e_p_ic=e_p_ic + e_p_temp
            if ( m == 0 ) e_p_as_ic=e_p_as_ic+e_p_temp
         end do    ! do loop over lms in block

         call reduce_scalar(e_p_ic, e_p_ic_global, 0)
         call reduce_scalar(e_p_as_ic, e_p_as_ic_global, 0)

         if ( l_master_rank ) then
            fac      =half*LFfac/r_icb*eScale
            e_p_ic   =fac*e_p_ic_global
            e_t_ic   =0.0_cp
            e_p_as_ic=fac*e_p_as_ic_global
            e_t_as_ic=0.0_cp
         end if

      end if  ! conducting inner core ?

      !-- Outside energy:
      nR=n_r_cmb
      e_p_os   =0.0_cp
      e_p_as_os=0.0_cp
      !do lm=2,lm_max
      do lm=1,n_mlo_loc
         l=map_mlo%i2l(lm)
         m=map_mlo%i2m(lm)
         if ( l == 0 ) cycle
         fac=real( l*l*(l+1),cp)
         e_p_temp=fac*cc2real(b(lm,nR),m)
         e_p_os  =e_p_os + e_p_temp
         if ( m == 0 ) e_p_as_os=e_p_as_os + e_p_temp
      end do

      call reduce_scalar(e_p_os, e_p_os_global, 0)
      call reduce_scalar(e_p_as_os, e_p_as_os_global, 0)

      if ( l_master_rank ) then
         fac      =half*LFfac/r_cmb*eScale
         e_p_os   =fac*e_p_os_global
         e_p_as_os=fac*e_p_as_os_global
      end if

      !-- External potential field energy in Uli case (n_imp=1)
      e_p_e     =0.0_cp
      e_p_as_e  =0.0_cp
      e_dipole_e=0.0_cp
      if ( n_imp == 1 ) then
         do lm=1,n_mlo_loc
            l=map_mlo%i2l(lm)
            m=map_mlo%i2m(lm)
            if ( l == 0 ) cycle
            fac=real(l*(l+1)**2*(2*l+1),cp)*one/(rrMP**(2*l+1)-one)
            e_p_temp=fac*cc2real(b(lm,nR),m)
            e_p_e   =e_p_e  + e_p_temp
            if ( m == 0 ) e_p_as_e =e_p_as_e  + e_p_temp
            if ( l == 1 ) e_dipole_e=e_dipole_e+e_p_temp
         end do

         call reduce_scalar(e_p_e, e_p_e_global, 0)
         call reduce_scalar(e_p_as_e, e_p_as_e_global, 0)
         call reduce_scalar(e_dipole_e, e_dipole_e_global, 0)

         if ( l_master_rank ) then
            fac       =half*LFfac/r_cmb**2*eScale
            e_p_e     =fac*e_p_e_global
            e_p_as_e  =fac*e_p_as_e_global
            e_dipole_e=fac*e_dipole_e_global
         end if
      end if


      !-- Output of OC and outside energies:
      if ( l_master_rank ) then
         if ( l_write ) then
            if ( l_save_out ) then
               open(newunit=n_e_mag_oc_file, file=e_mag_oc_file, &
               &    status='unknown', position='append')
            end if
            write(n_e_mag_oc_file,'(1P,ES20.12,12ES16.8)')   &
            &                             time*tScale,       &! 1
            &                             e_p,e_t,           &! 2,3
            &                             e_p_as,e_t_as,     &! 4,5
            &                             e_p_os,e_p_as_os,  &! 6,7
            &                             e_p_es,e_t_es,     &! 8,9
            &                             e_p_eas,e_t_eas,   &! 10,11
            &                             e_p_e,e_p_as_e      ! 12,13
            if ( l_save_out ) close(n_e_mag_oc_file)
         end if

         !-- Output of IC energies:
         if ( l_write .and. .not. l_full_sphere ) then
            if ( l_save_out ) then
               open(newunit=n_e_mag_ic_file, file=e_mag_ic_file, &
               &    status='unknown', position='append')
            end if
            write(n_e_mag_ic_file,'(1P,ES20.12,4ES16.8)')             &
            &                            time*tScale,                 &
            &                            e_p_ic,e_t_ic,               &
            &                            e_p_as_ic,e_t_as_ic
            if ( l_save_out ) close(n_e_mag_ic_file)
         end if
      end if
      b10=zero
      b11=zero

      call send_lm_pair_to_master(b(:,n_r_cmb),1,0,b10)
      call send_lm_pair_to_master(b(:,n_r_cmb),1,1,b11)

      if ( l_master_rank) then
         !-- Calculate pole position:
         rad =180.0_cp/pi

         theta_dip= rad*atan2(sqrt(two)*abs(b11),real(b10))
         if ( theta_dip < 0.0_cp ) theta_dip=180.0_cp+theta_dip
         if ( abs(b11) < 1.e-20_cp ) then
            phi_dip=0.0_cp
         else
            phi_dip=-rad*atan2(aimag(b11),real(b11))
         end if
         Dip      =e_dipole_ax/(e_p+e_t)
         DipCMB   =e_dipole_ax_cmb/e_cmb

         !-- Output of pole position:
         if ( l_write ) then
            if ( l_save_out ) then
               open(newunit=n_dipole_file, file=dipole_file,  &
               &    status='unknown', position='append')
            end if
            if ( e_p_e == 0 ) then
               e_p_e_ratio=0.0_cp
            else
               e_p_e_ratio=e_dipole_e/e_p_e
            end if

            write(n_dipole_file,'(1P,ES20.12,19ES14.6)')   &
            &            time*tScale,                      &! 1
            &            theta_dip,phi_dip,                &! 2,3
            &            Dip,                              &! 4
            &            DipCMB,                           &! 5
            &            e_dipole_ax_cmb/e_geo,            &! 6
            &            e_dipole/(e_p+e_t),               &! 7
            &            e_dip_cmb/e_cmb,                  &! 8
            &            e_dip_cmb/e_geo,                  &! 9
            &            e_dip_cmb,e_dipole_ax_cmb,        &! 10,11
            &            e_dipole,e_dipole_ax,             &! 12,13
            &            e_cmb,e_geo,                      &! 14,15
            &            e_p_e_ratio,                      &! 16
            &            (e_cmb-e_es_cmb)/e_cmb,           &! 17
            &            (e_cmb-e_as_cmb)/e_cmb,           &! 18
            &            (e_geo-e_es_geo)/e_geo,           &! 19
            &            (e_geo-e_as_geo)/e_geo             ! 20
            if ( l_save_out ) close(n_dipole_file)

            if ( l_earth_likeness ) then
               if ( l_save_out ) then
                  open(newunit=n_compliance_file, file=earth_compliance_file,  &
                  &    status='unknown', position='append')
               end if

               call get_Br_skew(bCMB,br_skew)
               fluxConcentration = br_skew-one

               if ( nad_global /= 0.0_cp ) then
                  ad=ad_global/nad_global
               else
                  ad=0.0_cp
               end if
               if ( asym_global /= 0.0_cp ) then
                  sym=sym_global/asym_global
               else
                  sym=0.0_cp
               end if
               if ( nzon_global /= 0.0_cp ) then
                  zon=zon_global/nzon_global
               else
                  zon=0.0_cp
               end if
               write(n_compliance_file,'(1P,ES20.12,4ES16.8)')    &
               &            time*tScale, ad,                      &
               &            sym, zon, fluxConcentration
               if ( l_save_out ) close(n_compliance_file)
            end if

         end if
         ! Store values needed for movie output:
         movieDipColat      =theta_dip
         movieDipLon        =phi_dip
         movieDipStrength   =e_dip_cmb/e_cmb
         movieDipStrengthGeo=e_dip_cmb/e_geo
      end if

   end subroutine get_e_mag
!----------------------------------------------------------------------------
   subroutine get_Br_skew(bCMB,Br_skew)
      !
      ! This subroutine calculates the skewness of the radial magnetic field
      ! at the outer boundary.
      ! bskew :=  <B^4> / <B^2>^2   where <..> is average over the surface
      !

      !-- Input variable
      complex(cp), intent(in) :: bCMB(lm_max)

      !-- Output variable
      real(cp), intent(out) :: Br_skew

      !-- Local variables
      real(cp) :: br(n_theta_max_comp, n_phi_max_comp)
      real(cp) :: glm, hlm, phi, fac
      real(cp) :: sinTh, sinTh_m, Br2tmp, Br2s, Br4s, Br2, Br4
      integer :: nPhi, m, l, nTheta, lm, lm1

      br(:,:) = 0.0_cp
      !-- First compute Br on the equidistant physical grid
      do nPhi=1,n_phi_max_comp
         phi = two*pi*(nPhi-1)/real(n_phi_max_comp,kind=cp)
         lm1 = 1
         do m=0,l_max_comp,minc
            do l=m,l_max_comp
               lm=map_glbl_st%lm2(l,m)
               fac = (-one)**(m)*l*sqrt(two*l+one)*osq4pi
               if ( m > 0 ) fac = fac*sqrt(two)
               glm =  fac*real(bCMB(lm))
               hlm = -fac*aimag(bCMB(lm))

               do nTheta=1,n_theta_max_comp
                  br(nTheta,nPhi)=br(nTheta,nPhi)+Plm_comp(lm1,nTheta)*( &
                  &               glm*cos(m*phi)+hlm*sin(m*phi))*(l+1)/r(n_r_cmb)**2
               end do
               lm1 = lm1+1
            end do
         end do
      end do

      !-- Then calculate skewness
      Br2 = 0.0_cp
      Br4 = 0.0_cp
      sinTh_m = 0.0_cp
      do nTheta=1,n_theta_max_comp
         sinTh=sin(pi*(nTheta-half)/real(n_theta_max_comp,cp))
         sinTh_m=sinTh_m+sinTh

         Br2s = 0.0_cp
         Br4s = 0.0_cp
         do nPhi=1,n_phi_max_comp
            Br2tmp = br(nTheta,nPhi)*br(nTheta,nPhi)
            Br2s = Br2s+Br2tmp
            Br4s = Br4s+Br2tmp*Br2tmp
         end do
         Br2=Br2+Br2s*sinTh
         Br4=Br4+Br4s*sinTh
      end do
      Br2=Br2/(sinTh_m*real(n_phi_max_comp,kind=cp))
      Br4=Br4/(sinTh_m*real(n_phi_max_comp,kind=cp))

      Br_skew= Br4/(Br2*Br2)

   end subroutine get_Br_skew
!----------------------------------------------------------------------------
end module magnetic_energy
