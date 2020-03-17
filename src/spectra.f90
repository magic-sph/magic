module spectra

   use parallel_mod
   use precision_mod
   use communications, only: reduce_radial
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_ic_maxMag, n_r_maxMag, &
       &                 n_r_ic_max, l_max, minc, n_r_cmb, n_r_icb
   use radial_functions, only: orho1, orho2, r_ic, chebt_ic, r, r_cmb,  &
       &                       rscheme_oc, or2, r_icb, dr_fac_ic
   use physical_parameters, only: LFfac
   use num_param, only: eScale, tScale
   use blocking, only: lo_map, st_map, llm, ulm, llmMag, ulmMag
   use horizontal_data, only: dLh
   use logic, only: l_mag, l_anel, l_cond_ic, l_heat, l_save_out, &
       &            l_energy_modes, l_2D_spectra
   use output_data, only: tag, log_file, n_log_file, m_max_modes
   use useful, only: cc2real, cc22real, abortRun
   use integration, only: rInt_R, rIntIC
   use constants, only: pi, vol_oc, half, one, four
   use mean_sd, only: mean_sd_type, mean_sd_2D_type

   implicit none
 
   private

   type(mean_sd_type) :: e_mag_p_l_ave, e_mag_p_m_ave
   type(mean_sd_type) :: e_mag_t_l_ave, e_mag_t_m_ave
   type(mean_sd_type) :: e_mag_cmb_l_ave, e_mag_cmb_m_ave

   type(mean_sd_type) :: u2_p_l_ave, u2_p_m_ave
   type(mean_sd_type) :: u2_t_l_ave, u2_t_m_ave

   type(mean_sd_type) :: e_kin_p_l_ave, e_kin_p_m_ave
   type(mean_sd_type) :: e_kin_t_l_ave, e_kin_t_m_ave

   type(mean_sd_2D_type) :: e_mag_p_r_l_ave, e_mag_p_r_m_ave
   type(mean_sd_2D_type) :: e_mag_t_r_l_ave, e_mag_t_r_m_ave

   type(mean_sd_2D_type) :: e_kin_p_r_l_ave, e_kin_p_r_m_ave
   type(mean_sd_2D_type) :: e_kin_t_r_l_ave, e_kin_t_r_m_ave

   type(mean_sd_type) :: T_l_ave, T_ICB_l_ave, dT_ICB_l_ave
   type(mean_sd_type) :: T_m_ave, T_ICB_m_ave, dT_ICB_m_ave

   integer :: n_am_kpol_file, n_am_ktor_file
   integer :: n_am_mpol_file, n_am_mtor_file
   character(len=72) :: am_kpol_file, am_ktor_file
   character(len=72) :: am_mpol_file, am_mtor_file

 
   public :: initialize_spectra, spectrum, spectrum_temp, &
   &         get_amplitude, finalize_spectra

contains

   subroutine initialize_spectra

      if ( l_mag ) then
         call e_mag_p_l_ave%initialize(1,l_max)
         call e_mag_p_m_ave%initialize(0,l_max)
         call e_mag_t_l_ave%initialize(1,l_max)
         call e_mag_t_m_ave%initialize(0,l_max)
         call e_mag_cmb_l_ave%initialize(1,l_max)
         call e_mag_cmb_m_ave%initialize(0,l_max)
      end if

      if ( l_anel ) then
         call u2_p_l_ave%initialize(1,l_max)
         call u2_p_m_ave%initialize(0,l_max)
         call u2_t_l_ave%initialize(1,l_max)
         call u2_t_m_ave%initialize(0,l_max)
      end if

      call e_kin_p_l_ave%initialize(1,l_max)
      call e_kin_p_m_ave%initialize(0,l_max)
      call e_kin_t_l_ave%initialize(1,l_max)
      call e_kin_t_m_ave%initialize(0,l_max)

      if ( l_2D_spectra ) then

         if ( l_mag ) then
            call e_mag_p_r_l_ave%initialize(1,n_r_max,1,l_max,.false.)
            call e_mag_p_r_m_ave%initialize(1,n_r_max,1,l_max+1,.false.)
            call e_mag_t_r_l_ave%initialize(1,n_r_max,1,l_max,.false.)
            call e_mag_t_r_m_ave%initialize(1,n_r_max,1,l_max+1,.false.)
         end if

         call e_kin_p_r_l_ave%initialize(1,n_r_max,1,l_max,.false.)
         call e_kin_p_r_m_ave%initialize(1,n_r_max,1,l_max+1,.false.)
         call e_kin_t_r_l_ave%initialize(1,n_r_max,1,l_max,.false.)
         call e_kin_t_r_m_ave%initialize(1,n_r_max,1,l_max+1,.false.)

      end if

      if ( l_heat ) then
         call T_l_ave%initialize(1,l_max+1)
         call T_ICB_l_ave%initialize(1,l_max+1)
         call dT_ICB_l_ave%initialize(1,l_max+1)
         call T_m_ave%initialize(1,l_max+1)
         call T_ICB_m_ave%initialize(1,l_max+1)
         call dT_ICB_m_ave%initialize(1,l_max+1)
      end if

      am_kpol_file='am_kin_pol.'//tag
      am_ktor_file='am_kin_tor.'//tag
      am_mpol_file='am_mag_pol.'//tag
      am_mtor_file='am_mag_tor.'//tag

      if ( l_master_rank .and. (.not. l_save_out) ) then
         if ( l_mag .and. l_energy_modes ) then
            open(newunit=n_am_kpol_file,file=am_kpol_file,status='new', &
            &    form='unformatted')
            open(newunit=n_am_ktor_file,file=am_ktor_file,status='new', &
            &    form='unformatted')
            open(newunit=n_am_mpol_file,file=am_mpol_file,status='new', &
            &    form='unformatted')
            open(newunit=n_am_mtor_file,file=am_mtor_file,status='new', &
            &    form='unformatted')
         end if
      end if

   end subroutine initialize_spectra
!----------------------------------------------------------------------------
   subroutine finalize_spectra

      if ( l_mag ) then
         call e_mag_p_l_ave%finalize()
         call e_mag_p_m_ave%finalize()
         call e_mag_t_l_ave%finalize()
         call e_mag_t_m_ave%finalize()
         call e_mag_cmb_l_ave%finalize()
         call e_mag_cmb_m_ave%finalize()
      end if

      if ( l_anel ) then
         call u2_p_l_ave%finalize()
         call u2_p_m_ave%finalize()
         call u2_t_l_ave%finalize()
         call u2_t_m_ave%finalize()
      end if

      call e_kin_p_l_ave%finalize()
      call e_kin_p_m_ave%finalize()
      call e_kin_t_l_ave%finalize()
      call e_kin_t_m_ave%finalize()

      if ( l_2D_spectra ) then
         if ( l_mag ) then
            call e_mag_p_r_l_ave%finalize() 
            call e_mag_p_r_m_ave%finalize()
            call e_mag_t_r_l_ave%finalize()
            call e_mag_t_r_m_ave%finalize()
         end if
         call e_kin_p_r_l_ave%finalize()
         call e_kin_p_r_m_ave%finalize()
         call e_kin_t_r_l_ave%finalize()
         call e_kin_t_r_m_ave%finalize()
      end if

      if ( l_heat ) then
         call T_l_ave%finalize()
         call T_ICB_l_ave%finalize()
         call dT_ICB_l_ave%finalize()
         call T_m_ave%finalize()
         call T_ICB_m_ave%finalize()
         call dT_ICB_m_ave%finalize()
      end if

      if ( l_master_rank .and. (.not. l_save_out) ) then
         if ( l_mag .and. l_energy_modes ) then
            close(n_am_kpol_file)
            close(n_am_ktor_file)
            close(n_am_mpol_file)
            close(n_am_mtor_file)
         end if
      end if

   end subroutine finalize_spectra
!----------------------------------------------------------------------------
   subroutine spectrum(n_spec,time,l_avg,n_time_ave,l_stop_time,time_passed, &
              &        time_norm,w,dw,z,b,db,aj,b_ic,db_ic,aj_ic)
      !
      !  calculates magnetic energy  = 1/2 Integral(B^2 dV)
      !  integration in theta,phi by summation over harmonic coeffs.
      !  integration in r by Chebycheff integrals
      !
      !  Output:
      !  enbp: Total poloidal        enbt: Total toroidal
      !  apome: Axisym. poloidal     atome: Axisym. toroidal
      !
    
      !-- Input of variables:
      integer,     intent(in) :: n_spec     ! number of spectrum/call, file
      integer,     intent(in) :: n_time_ave
      real(cp),    intent(in) :: time
      real(cp),    intent(in) :: time_passed
      real(cp),    intent(in) :: time_norm
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      logical,     intent(in) :: l_avg, l_stop_time

      !-- Output:
      real(cp) :: b_rms
      real(cp) :: e_mag_p_l(l_max),e_mag_t_l(l_max)
      real(cp) :: e_kin_p_l(l_max),e_kin_t_l(l_max)
      real(cp) :: e_mag_p_ic_l(l_max),e_mag_t_ic_l(l_max)
      real(cp) :: u2_p_l(l_max),u2_t_l(l_max)
    
      real(cp) :: e_mag_p_m(l_max+1),e_mag_t_m(l_max+1)
      real(cp) :: e_kin_p_m(l_max+1),e_kin_t_m(l_max+1)
      real(cp) :: e_mag_p_ic_m(l_max+1),e_mag_t_ic_m(l_max+1)
      real(cp) :: u2_p_m(l_max+1),u2_t_m(l_max+1)
    
      real(cp) :: e_mag_cmb_l(l_max)
      real(cp) :: e_mag_cmb_m(l_max+1)
      real(cp) :: e_kin_nearSurf_l(l_max)
      real(cp) :: e_kin_nearSurf_m(l_max+1)
    
      real(cp) :: eCMB(l_max),eCMB_global(l_max)
    
      !-- local:
      character(len=14) :: string
      character(len=72) :: mag_spec_file,kin_spec_file,u2_spec_file
      integer :: n_mag_spec_file,n_kin_spec_file,n_u2_spec_file
      integer :: n_r,lm,ml,l,mc,m,n_const
    
      real(cp) :: r_ratio,O_r_icb_E_2
      real(cp) :: e_mag_p_temp,e_mag_t_temp
      real(cp) :: e_kin_p_temp,e_kin_t_temp
      real(cp) :: u2_p_temp,u2_t_temp
      real(cp) :: O_surface
      real(cp) :: fac_mag,fac_kin
      real(cp) :: nearSurfR

      real(cp) :: e_mag_p_r_l(n_r_max,l_max),e_mag_p_r_l_global(n_r_max,l_max)
      real(cp) :: e_mag_t_r_l(n_r_max,l_max),e_mag_t_r_l_global(n_r_max,l_max)
      real(cp) :: e_kin_p_r_l(n_r_max,l_max),e_kin_p_r_l_global(n_r_max,l_max)
      real(cp) :: e_kin_t_r_l(n_r_max,l_max),e_kin_t_r_l_global(n_r_max,l_max)
      real(cp) :: u2_p_r_l(n_r_max,l_max),u2_p_r_l_global(n_r_max,l_max)
      real(cp) :: u2_t_r_l(n_r_max,l_max),u2_t_r_l_global(n_r_max,l_max)
      real(cp) :: e_mag_p_r_m(n_r_max,l_max+1),e_mag_p_r_m_global(n_r_max,l_max+1)
      real(cp) :: e_mag_t_r_m(n_r_max,l_max+1),e_mag_t_r_m_global(n_r_max,l_max+1)
      real(cp) :: e_kin_p_r_m(n_r_max,l_max+1),e_kin_p_r_m_global(n_r_max,l_max+1)
      real(cp) :: e_kin_t_r_m(n_r_max,l_max+1),e_kin_t_r_m_global(n_r_max,l_max+1)
      real(cp) :: u2_p_r_m(n_r_max,l_max+1),u2_p_r_m_global(n_r_max,l_max+1)
      real(cp) :: u2_t_r_m(n_r_max,l_max+1),u2_t_r_m_global(n_r_max,l_max+1)
    
      real(cp) :: e_mag_p_ic_r_l(n_r_ic_max,l_max)
      real(cp) :: e_mag_p_ic_r_l_global(n_r_ic_max,l_max)
      real(cp) :: e_mag_t_ic_r_l(n_r_ic_max,l_max)
      real(cp) :: e_mag_t_ic_r_l_global(n_r_ic_max,l_max)
      real(cp) :: e_mag_p_ic_r_m(n_r_ic_max,l_max+1)
      real(cp) :: e_mag_p_ic_r_m_global(n_r_ic_max,l_max+1)
      real(cp) :: e_mag_t_ic_r_m(n_r_ic_max,l_max+1)
      real(cp) :: e_mag_t_ic_r_m_global(n_r_ic_max,l_max+1)
    
      complex(cp) :: r_dr_b
    
    
      eCMB(:)=0.0_cp
    
      do n_r=1,n_r_max
    
         do l=1,l_max
            if ( l_mag ) then
               e_mag_p_r_l(n_r,l)=0.0_cp
               e_mag_t_r_l(n_r,l)=0.0_cp
            end if
            if ( l_anel ) then
               u2_p_r_l(n_r,l)=0.0_cp
               u2_t_r_l(n_r,l)=0.0_cp
            end if
            e_kin_p_r_l(n_r,l)=0.0_cp
            e_kin_t_r_l(n_r,l)=0.0_cp
         end do
         do mc=1,l_max+1
            if ( l_mag ) then
               e_mag_p_r_m(n_r,mc)=0.0_cp
               e_mag_t_r_m(n_r,mc)=0.0_cp
            end if
            if ( l_anel ) then
               u2_p_r_m(n_r,mc)=0.0_cp
               u2_t_r_m(n_r,mc)=0.0_cp
            end if
            e_kin_p_r_m(n_r,mc)=0.0_cp
            e_kin_t_r_m(n_r,mc)=0.0_cp
         end do
    
         !do lm=2,lm_max
         do lm=max(llm,2),ulm
    
            l  =lo_map%lm2l(lm)
            m  =lo_map%lm2m(lm)
            mc=m+1
    
            if ( l_mag ) then
               e_mag_p_temp= dLh(st_map%lm2(l,m)) * ( dLh(st_map%lm2(l,m))*     &
               &             or2(n_r)*cc2real(b(lm,n_r),m) + cc2real(db(lm,n_r),m) )
               e_mag_t_temp=dLh(st_map%lm2(l,m))*cc2real(aj(lm,n_r),m)
            end if
            if ( l_anel ) then
               u2_p_temp=  orho2(n_r)*dLh(st_map%lm2(l,m)) *  (                   &
               &             dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
               &             cc2real(dw(lm,n_r),m) )
               u2_t_temp=orho2(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)
            end if
            e_kin_p_temp= orho1(n_r)*dLh(st_map%lm2(l,m)) *  (                   &
            &               dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
            &               cc2real(dw(lm,n_r),m) )
            e_kin_t_temp=orho1(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)
    
            !----- l-spectra:
            if ( l_mag ) then
               e_mag_p_r_l(n_r,l) = e_mag_p_r_l(n_r,l) + e_mag_p_temp
               e_mag_t_r_l(n_r,l) = e_mag_t_r_l(n_r,l) + e_mag_t_temp
               if ( m == 0 .and. n_r == n_r_cmb ) eCMB(l)=e_mag_p_temp
            end if
            if ( l_anel ) then
               u2_p_r_l(n_r,l) = u2_p_r_l(n_r,l) + u2_p_temp
               u2_t_r_l(n_r,l) = u2_t_r_l(n_r,l) + u2_t_temp
            end if
            e_kin_p_r_l(n_r,l) = e_kin_p_r_l(n_r,l) + e_kin_p_temp
            e_kin_t_r_l(n_r,l) = e_kin_t_r_l(n_r,l) + e_kin_t_temp
    
            !----- m-spectra:
            if ( l_mag ) then
               e_mag_p_r_m(n_r,mc) = e_mag_p_r_m(n_r,mc) + e_mag_p_temp
               e_mag_t_r_m(n_r,mc) = e_mag_t_r_m(n_r,mc) + e_mag_t_temp
            end if
            if ( l_anel ) then
               u2_p_r_m(n_r,mc) = u2_p_r_m(n_r,mc) + u2_p_temp
               u2_t_r_m(n_r,mc) = u2_t_r_m(n_r,mc) + u2_t_temp
            end if
            e_kin_p_r_m(n_r,mc)=e_kin_p_r_m(n_r,mc) + e_kin_p_temp
            e_kin_t_r_m(n_r,mc)=e_kin_t_r_m(n_r,mc) + e_kin_t_temp
    
         end do    ! do loop over lms in block
    
      end do    ! radial grid points
    
      ! ----------- We need a reduction here ----------------
      ! first the l-spectra
      if ( l_mag ) then
         call reduce_radial(e_mag_p_r_l, e_mag_p_r_l_global, 0)
         call reduce_radial(e_mag_t_r_l, e_mag_t_r_l_global, 0)
         call reduce_radial(e_mag_p_r_m, e_mag_p_r_m_global, 0)
         call reduce_radial(e_mag_t_r_m, e_mag_t_r_m_global, 0)
         call reduce_radial(eCMB, eCMB_global, 0)
      end if
      if ( l_anel ) then
         call reduce_radial(u2_p_r_l, u2_p_r_l_global, 0)
         call reduce_radial(u2_t_r_l, u2_t_r_l_global, 0)
         call reduce_radial(u2_p_r_m, u2_p_r_m_global, 0)
         call reduce_radial(u2_t_r_m, u2_t_r_m_global, 0)
      end if
      call reduce_radial(e_kin_p_r_l, e_kin_p_r_l_global, 0)
      call reduce_radial(e_kin_t_r_l, e_kin_t_r_l_global, 0)
      call reduce_radial(e_kin_p_r_m, e_kin_p_r_m_global, 0)
      call reduce_radial(e_kin_t_r_m, e_kin_t_r_m_global, 0)

      ! now switch to coord_r 0 for the postprocess
    
      if ( coord_r == 0 ) then
         ! Getting appropriate radius index for e_kin_nearSurf spectra
         nearSurfR = r_cmb-0.01_cp
         do n_r=2,n_r_max
            if ( r(n_r-1) > nearSurfR .and. r(n_r)  <= nearSurfR ) then
               if ( r(n_r-1)-nearSurfR < nearSurfR-r(n_r) ) then
                  n_const=n_r-1
               else
                  n_const=n_r
               end if
            end if
         end do
 
         !-- Save CMB energy spectra:
         O_surface=one/(four*pi*r(1)*r(1))
    
         if ( l_mag ) then
            b_rms=0.0_cp
            do l=1,l_max
               e_mag_cmb_l(l)=e_mag_p_r_l_global(1,l)
               b_rms=b_rms + e_mag_cmb_l(l)
            end do
            b_rms=sqrt(b_rms*O_surface)
            do mc=1,l_max+1
               e_mag_cmb_m(mc)=e_mag_p_r_m_global(1,mc)
            end do
         end if
    
         !-- Save nearSurf kin energy spectra:
         do l=1,l_max
            e_kin_nearSurf_l(l)=e_kin_p_r_l_global(n_const,l)
         end do
         do mc=1,l_max+1
            e_kin_nearSurf_m(mc)=e_kin_p_r_m_global(n_const,mc)
         end do
    
         !-- Radial Integrals:
         fac_mag=half*LFfac*eScale
         fac_kin=half*eScale
         do l=1,l_max
            if ( l_mag ) then
               e_mag_p_l(l)=fac_mag*rInt_R(e_mag_p_r_l_global(:,l),r,rscheme_oc)
               e_mag_t_l(l)=fac_mag*rInt_R(e_mag_t_r_l_global(:,l),r,rscheme_oc)
               e_mag_cmb_l(l)=fac_mag*e_mag_cmb_l(l)
            end if
            if ( l_anel ) then
               u2_p_l(l)  =fac_kin*rInt_R(u2_p_r_l_global(:,l),r,rscheme_oc)
               u2_t_l(l)  =fac_kin*rInt_R(u2_t_r_l_global(:,l),r,rscheme_oc)
            end if
            e_kin_p_l(l)  =fac_kin*rInt_R(e_kin_p_r_l_global(:,l),r,rscheme_oc)
            e_kin_t_l(l)  =fac_kin*rInt_R(e_kin_t_r_l_global(:,l),r,rscheme_oc)
            e_kin_nearSurf_l(l)=fac_kin*e_kin_nearSurf_l(l)
         end do
         do m=1,l_max+1 ! Note: counter m is actual order+1
            if ( l_mag )  then
               e_mag_p_m(m)=fac_mag*rInt_R(e_mag_p_r_m_global(:,m),r,rscheme_oc)
               e_mag_t_m(m)=fac_mag*rInt_R(e_mag_t_r_m_global(:,m),r,rscheme_oc)
               e_mag_cmb_m(m)=fac_mag*e_mag_cmb_m(m)
            end if
            if ( l_anel ) then
               u2_p_m(m)   =fac_kin*rInt_R(u2_p_r_m_global(:,m),r,rscheme_oc)
               u2_t_m(m)   =fac_kin*rInt_R(u2_t_r_m_global(:,m),r,rscheme_oc)
            end if
            e_kin_p_m(m)   =fac_kin*rInt_R(e_kin_p_r_m_global(:,m),r,rscheme_oc)
            e_kin_t_m(m)   =fac_kin*rInt_R(e_kin_t_r_m_global(:,m),r,rscheme_oc)
            e_kin_nearSurf_m(m)=fac_kin*e_kin_nearSurf_m(m)
         end do
      end if

      !-- Averaging:
      if ( coord_r == 0 .and. l_avg ) then
         if ( l_mag ) then
            call e_mag_p_l_ave%compute(e_mag_p_l, n_time_ave, time_passed, time_norm)
            call e_mag_t_l_ave%compute(e_mag_t_l, n_time_ave, time_passed, time_norm)
            call e_mag_p_m_ave%compute(e_mag_p_m, n_time_ave, time_passed, time_norm)
            call e_mag_t_m_ave%compute(e_mag_t_m, n_time_ave, time_passed, time_norm)
            call e_mag_cmb_l_ave%compute(e_mag_cmb_l, n_time_ave, time_passed, &
                 &                       time_norm)
            call e_mag_cmb_m_ave%compute(e_mag_cmb_m, n_time_ave, time_passed, &
                 &                       time_norm)

            !-- Averaging of 2D spectra
            if ( l_2D_spectra ) then
               call e_mag_p_r_l_ave%compute(e_mag_p_r_l_global, n_time_ave, &
                    &                       time_passed, time_norm)
               call e_mag_t_r_l_ave%compute(e_mag_t_r_l_global, n_time_ave, &
                    &                       time_passed, time_norm)
               call e_mag_p_r_m_ave%compute(e_mag_p_r_m_global, n_time_ave, &
                    &                       time_passed, time_norm)
               call e_mag_t_r_m_ave%compute(e_mag_t_r_m_global, n_time_ave, &
                    &                       time_passed, time_norm)
            end if
         end if

         if ( l_anel ) then
            call u2_p_l_ave%compute(u2_p_l, n_time_ave, time_passed, time_norm)
            call u2_t_l_ave%compute(u2_t_l, n_time_ave, time_passed, time_norm)
            call u2_p_m_ave%compute(u2_p_m, n_time_ave, time_passed, time_norm)
            call u2_t_m_ave%compute(u2_t_m, n_time_ave, time_passed, time_norm)
         end if

         call e_kin_p_l_ave%compute(e_kin_p_l, n_time_ave, time_passed, time_norm)
         call e_kin_t_l_ave%compute(e_kin_t_l, n_time_ave, time_passed, time_norm)
         call e_kin_p_m_ave%compute(e_kin_p_m, n_time_ave, time_passed, time_norm)
         call e_kin_t_m_ave%compute(e_kin_t_m, n_time_ave, time_passed, time_norm)

         !-- Averaging of 2D spectra
         if ( l_2D_spectra ) then
            call e_kin_p_r_l_ave%compute(e_kin_p_r_l_global, n_time_ave, &
                 &                       time_passed, time_norm)
            call e_kin_t_r_l_ave%compute(e_kin_t_r_l_global, n_time_ave, &
                 &                       time_passed, time_norm)
            call e_kin_p_r_m_ave%compute(e_kin_p_r_m_global, n_time_ave, &
                 &                       time_passed, time_norm)
            call e_kin_t_r_m_ave%compute(e_kin_t_r_m_global, n_time_ave, &
                 &                       time_passed, time_norm)
         end if
      end if
    
      !-- inner core:
    
      if ( l_cond_ic ) then
    
         O_r_icb_E_2=one/(r_ic(1)*r_ic(1))
         do n_r=1,n_r_ic_max
            r_ratio=r_ic(n_r)/r_ic(1)
            do mc=1,l_max+1
               e_mag_p_ic_r_m(n_r,mc)=0.0_cp
               e_mag_t_ic_r_m(n_r,mc)=0.0_cp
            end do
            do l=1,l_max
               e_mag_p_ic_r_l(n_r,l)=0.0_cp
               e_mag_t_ic_r_l(n_r,l)=0.0_cp
            end do
            !do lm=2,lm_max
            do lm=max(llm,2),ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)
               mc=m+1
               r_dr_b=r_ic(n_r)*db_ic(lm,n_r)
    
               e_mag_p_temp=dLh(st_map%lm2(l,m))*O_r_icb_E_2*r_ratio**(2*l) * ( &
               &            real((2*l+1)*(l+1),cp)*cc2real(b_ic(lm,n_r),m)   +  &
               &            real(2*(l+1),cp)*cc22real(b_ic(lm,n_r),r_dr_b,m) +  &
               &            cc2real(r_dr_b,m) )
               e_mag_t_temp=dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) * &
               &            cc2real(aj_ic(lm,n_r),m)
    
               e_mag_p_ic_r_l(n_r,l)=e_mag_p_ic_r_l(n_r,l) + e_mag_p_temp
               e_mag_t_ic_r_l(n_r,l)=e_mag_t_ic_r_l(n_r,l) + e_mag_t_temp
               e_mag_p_ic_r_m(n_r,mc)=e_mag_p_ic_r_m(n_r,mc) + e_mag_p_temp
               e_mag_t_ic_r_m(n_r,mc)=e_mag_t_ic_r_m(n_r,mc) + e_mag_t_temp
            end do  ! loop over lm's
         end do ! loop over radial levels
    
         call reduce_radial(e_mag_p_ic_r_l, e_mag_p_ic_r_l_global, 0)
         call reduce_radial(e_mag_t_ic_r_l, e_mag_t_ic_r_l_global, 0)
         call reduce_radial(e_mag_p_ic_r_m, e_mag_p_ic_r_m_global, 0)
         call reduce_radial(e_mag_t_ic_r_m, e_mag_t_ic_r_m_global, 0)
    
         if ( coord_r == 0 ) then
            !----- Radial Integrals:
            fac_mag=LFfac*half*eScale
            do l=1,l_max
               e_mag_p_ic_l(l)=fac_mag*rIntIC(e_mag_p_ic_r_l_global(1,l),    &
               &                              n_r_ic_max,dr_fac_ic,chebt_ic)
               e_mag_t_ic_l(l)=fac_mag*rIntIC(e_mag_t_ic_r_l_global(1,l),    &
               &                              n_r_ic_max,dr_fac_ic,chebt_ic)
            end do
            do m=1,l_max+1
               e_mag_p_ic_m(m)=fac_mag*rIntIC(e_mag_p_ic_r_m_global(1,m),    &
               &                              n_r_ic_max,dr_fac_ic,chebt_ic)
               e_mag_t_ic_m(m)=fac_mag*rIntIC(e_mag_t_ic_r_m_global(1,m),    &
               &                              n_r_ic_max,dr_fac_ic,chebt_ic)
            end do
         end if
      else
         do l=1,l_max
            e_mag_p_ic_l(l)=0.0_cp
            e_mag_t_ic_l(l)=0.0_cp
         end do
         do mc=1,l_max+1
            e_mag_p_ic_m(mc)=0.0_cp
            e_mag_t_ic_m(mc)=0.0_cp
         end do
      end if  ! conducting inner core ?

      if ( l_master_rank ) then
         !-- Output into files:
         if ( l_mag ) then
            write(string, *) n_spec
            mag_spec_file='mag_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=n_mag_spec_file, file=mag_spec_file, status='unknown')
            if ( n_spec == 0 ) then
               write(n_mag_spec_file,'(1x, &
               &           ''Magnetic energy spectra of time averaged field:'')')
            else
               write(n_mag_spec_file,'(1x, &
               &           ''Magnetic energy spectra at time:'', &
               &           ES20.12)') time*tScale
            end if
            write(n_mag_spec_file,'(1p,i4,11ES16.8)')            &
            &     0,0.0_cp,e_mag_p_m(1)   ,0.0_cp,e_mag_t_m(1),  &
            &     0.0_cp,e_mag_p_ic_m(1),0.0_cp,e_mag_t_ic_m(1), &
            &     0.0_cp,e_mag_cmb_m(1),0.0_cp
            do ml=1,l_max
               write(n_mag_spec_file,'(1p,i4,11ES16.8)')  &
               &     ml,e_mag_p_l(ml),   e_mag_p_m(ml+1), &
               &     e_mag_t_l(ml),   e_mag_t_m(ml+1),    &
               &     e_mag_p_ic_l(ml),e_mag_p_ic_m(ml+1), &
               &     e_mag_t_ic_l(ml),e_mag_t_ic_m(ml+1), &
               &     e_mag_cmb_l(ml), e_mag_cmb_m(ml+1),  &
               &     eCMB_global(ml)
            end do
            close(n_mag_spec_file)

            if ( l_2D_spectra ) then    
               mag_spec_file='2D_mag_spec_'//trim(adjustl(string))//'.'//tag
               open(newunit=n_mag_spec_file, file=mag_spec_file, &
               &    status='unknown', form='unformatted')
               write(n_mag_spec_file) time*tScale,n_r_max,l_max,minc
               write(n_mag_spec_file) r
               write(n_mag_spec_file) fac_mag*e_mag_p_r_l_global
               write(n_mag_spec_file) fac_mag*e_mag_p_r_m_global
               write(n_mag_spec_file) fac_mag*e_mag_t_r_l_global
               write(n_mag_spec_file) fac_mag*e_mag_t_r_m_global
               close(n_mag_spec_file)
            end if
            
            if ( l_avg .and. l_stop_time ) then
               !-- Output: at end of run
               mag_spec_file='mag_spec_ave.'//tag
               open(newunit=n_mag_spec_file, file=mag_spec_file, &
               &    status='unknown')
               call e_mag_p_l_ave%finalize_SD(time_norm)
               call e_mag_t_l_ave%finalize_SD(time_norm)
               call e_mag_p_m_ave%finalize_SD(time_norm)
               call e_mag_t_m_ave%finalize_SD(time_norm)
               call e_mag_cmb_l_ave%finalize_SD(time_norm)
               call e_mag_cmb_m_ave%finalize_SD(time_norm)
               write(n_mag_spec_file,'(2X,1P,I4,12ES16.8)') 0,                  &
               &     0.0_cp,e_mag_p_m_ave%mean(0),0.0_cp,e_mag_t_m_ave%mean(0), &
               &     0.0_cp,e_mag_cmb_m_ave%mean(0),0.0_cp,e_mag_p_m_ave%SD(0), &
               &     0.0_cp,e_mag_t_m_ave%SD(0),0.0_cp,e_mag_cmb_m_ave%SD(0)
               do l=1,l_max
                  write(n_mag_spec_file,'(2X,1P,I4,12ES16.8)') l,         &
                  &     e_mag_p_l_ave%mean(l), e_mag_p_m_ave%mean(l),     &
                  &     e_mag_t_l_ave%mean(l), e_mag_t_m_ave%mean(l),     &
                  &     e_mag_cmb_l_ave%mean(l), e_mag_cmb_m_ave%mean(l), &
                  &     e_mag_p_l_ave%SD(l), e_mag_p_m_ave%SD(l),         &
                  &     e_mag_t_l_ave%SD(l), e_mag_t_m_ave%SD(l),         &
                  &     e_mag_cmb_l_ave%SD(l), e_mag_cmb_m_ave%SD(l)
               end do
               close(n_mag_spec_file)

               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown',&
                  &    position='append')
               end if
               write(n_log_file,"(/,A,A)")                             &
               &     ' ! TIME AVERAGED EMAG SPECTRA STORED IN FILE: ', &
               &     mag_spec_file 
               write(n_log_file,"(A,I5)")                         &
               &     ' !              No. of averaged spectra: ', &
               &     n_time_ave
               if ( l_save_out ) close(n_log_file)
              
               if ( l_2D_spectra ) then
                  mag_spec_file='2D_mag_spec_ave.'//tag
                  open(newunit=n_mag_spec_file, file=mag_spec_file, &
                  &    status='unknown', form='unformatted')
                  write(n_mag_spec_file) n_r_max,l_max,minc
                  write(n_mag_spec_file) r
                  write(n_mag_spec_file) fac_mag*e_mag_p_r_l_ave%mean
                  write(n_mag_spec_file) fac_mag*e_mag_p_r_m_ave%mean
                  write(n_mag_spec_file) fac_mag*e_mag_t_r_l_ave%mean
                  write(n_mag_spec_file) fac_mag*e_mag_t_r_m_ave%mean
                  close(n_mag_spec_file)
               end if
            end if
         end if
    
         if ( l_anel ) then
            write(string, *) n_spec
            u2_spec_file='u2_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=n_u2_spec_file, file=u2_spec_file, status='unknown')
            if ( n_spec == 0 ) then
               write(n_u2_spec_file,'(1x, &
               &          ''Velocity square spectra of time averaged field:'')')
            else
               write(n_u2_spec_file,'(1x,                       &
               &          ''Velocity square spectra at time:'', &
               &          ES20.12)') time*tScale
            end if
            write(n_u2_spec_file,'(1p,i4,4ES16.8)') &
            &        0,0.0_cp,u2_p_m(1),0.0_cp,u2_t_m(1)
            do ml=1,l_max
               write(n_u2_spec_file,'(1p,i4,4ES16.8)') &
               &     ml,u2_p_l(ml),u2_p_m(ml+1),u2_t_l(ml),u2_t_m(ml+1)
            end do
            close(n_u2_spec_file)
    
            if ( l_avg .and. l_stop_time ) then
               !-- Output: at end of run
               u2_spec_file='u2_spec_ave.'//tag
               open(newunit=n_u2_spec_file, file=u2_spec_file, &
               &    status='unknown')
               call u2_p_l_ave%finalize_SD(time_norm)
               call u2_t_l_ave%finalize_SD(time_norm)
               call u2_p_m_ave%finalize_SD(time_norm)
               call u2_t_m_ave%finalize_SD(time_norm)
               write(n_u2_spec_file,'(2X,1P,I4,8ES16.8)') 0,              &
               &     0.0_cp,u2_p_m_ave%mean(0),0.0_cp,u2_t_m_ave%mean(0), &
               &     0.0_cp,u2_p_m_ave%SD(0),0.0_cp,u2_t_m_ave%SD(0)
               do l=1,l_max
                  write(n_u2_spec_file,'(2X,1P,I4,8ES16.8)') l, &
                  &     u2_p_l_ave%mean(l), u2_p_m_ave%mean(l), &
                  &     u2_t_l_ave%mean(l), u2_t_m_ave%mean(l), &
                  &     u2_p_l_ave%SD(l), u2_p_m_ave%SD(l),     &
                  &     u2_t_l_ave%SD(l), u2_t_m_ave%SD(l)
               end do
               close(n_u2_spec_file)

               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown',&
                  &    position='append')
               end if
               write(n_log_file,"(/,A,A)")                           &
               &     ' ! TIME AVERAGED U2 SPECTRA STORED IN FILE: ', &
               &     u2_spec_file
               write(n_log_file,"(A,I5)")                         &
               &     ' !              No. of averaged spectra: ', &
               &     n_time_ave
               if ( l_save_out ) close(n_log_file)
            end if
         end if
    
         write(string, *) n_spec
         kin_spec_file='kin_spec_'//trim(adjustl(string))//'.'//tag
         open(newunit=n_kin_spec_file, file=kin_spec_file, status='unknown')
         if ( n_spec == 0 ) then
            write(n_kin_spec_file,'(1x, &
            &           ''Kinetic energy spectra of time averaged field:'')')
         else
            write(n_kin_spec_file,'(1x,                      &
            &           ''Kinetic energy spectra at time:'', &
            &           ES20.12)') time*tScale
         end if
         write(n_kin_spec_file,'(1p,i4,6ES16.8)')            &
         &     0,0.0_cp,e_kin_p_m(1),0.0_cp,e_kin_t_m(1),    &
         &     0.0_cp, e_kin_nearSurf_m(1)
         do ml=1,l_max
            write(n_kin_spec_file,'(1p,i4,6ES16.8)')    &
            &     ml,e_kin_p_l(ml),e_kin_p_m(ml+1),     &
            &     e_kin_t_l(ml),e_kin_t_m(ml+1),        &
            &     e_kin_nearSurf_l(ml), e_kin_nearSurf_m(ml+1)
         end do
         close(n_kin_spec_file)
    
         if ( l_2D_spectra ) then
            kin_spec_file='2D_kin_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=n_kin_spec_file, file=kin_spec_file, status='unknown', &
            &    form='unformatted')
            write(n_kin_spec_file) time*tScale,n_r_max,l_max,minc
            write(n_kin_spec_file) r
            write(n_kin_spec_file) fac_kin*e_kin_p_r_l_global
            write(n_kin_spec_file) fac_kin*e_kin_p_r_m_global
            write(n_kin_spec_file) fac_kin*e_kin_t_r_l_global
            write(n_kin_spec_file) fac_kin*e_kin_t_r_m_global
            close(n_kin_spec_file)
         end if

         if ( l_avg .and. l_stop_time ) then
            !-- Output: at end of run
            kin_spec_file='kin_spec_ave.'//tag
            open(newunit=n_kin_spec_file, file=kin_spec_file, status='unknown')
            call e_kin_p_l_ave%finalize_SD(time_norm)
            call e_kin_t_l_ave%finalize_SD(time_norm)
            call e_kin_p_m_ave%finalize_SD(time_norm)
            call e_kin_t_m_ave%finalize_SD(time_norm)
            write(n_kin_spec_file,'(2X,1P,I4,8ES16.8)') 0, &
            &     0.0_cp,e_kin_p_m_ave%mean(0),0.0_cp,e_kin_t_m_ave%mean(0), &
            &     0.0_cp,e_kin_p_m_ave%SD(0),0.0_cp,e_kin_t_m_ave%SD(0)
            do l=1,l_max
               write(n_kin_spec_file,'(2X,1P,I4,8ES16.8)') l,      &
               &     e_kin_p_l_ave%mean(l), e_kin_p_m_ave%mean(l), &
               &     e_kin_t_l_ave%mean(l), e_kin_t_m_ave%mean(l), &
               &     e_kin_p_l_ave%SD(l), e_kin_p_m_ave%SD(l),     &
               &     e_kin_t_l_ave%SD(l), e_kin_t_m_ave%SD(l)
            end do
            close(n_kin_spec_file)

            if ( l_save_out ) then
               open(newunit=n_log_file, file=log_file, status='unknown', &
               &    position='append')
            end if
            write(n_log_file,"(/,A,A)")                             &
            &     ' ! TIME AVERAGED EKIN SPECTRA STORED IN FILE: ', &
            &     kin_spec_file
            write(n_log_file,"(A,I5)")                         &
            &     ' !              No. of averaged spectra: ', &
            &     n_time_ave
            if ( l_save_out ) close(n_log_file)

            if ( l_2D_spectra ) then
               kin_spec_file='2D_kin_spec_ave.'//tag
               open(newunit=n_kin_spec_file, file=kin_spec_file, &
               &    status='unknown', form='unformatted')
               write(n_kin_spec_file) n_r_max,l_max,minc
               write(n_kin_spec_file) r
               write(n_kin_spec_file) fac_kin*e_kin_p_r_l_ave%mean
               write(n_kin_spec_file) fac_kin*e_kin_p_r_m_ave%mean
               write(n_kin_spec_file) fac_kin*e_kin_t_r_l_ave%mean
               write(n_kin_spec_file) fac_kin*e_kin_t_r_m_ave%mean
               close(n_kin_spec_file)
            end if
         end if
      end if
    
   end subroutine spectrum
!----------------------------------------------------------------------------
   subroutine spectrum_temp(n_spec,time,l_avg,n_time_ave,l_stop_time,   &
              &             time_passed,time_norm,s,ds)

      !-- Direct input:
      real(cp),     intent(in) :: time
      integer,     intent(in) :: n_time_ave
      integer,     intent(in) :: n_spec
      logical,     intent(in) :: l_stop_time
      real(cp),    intent(in) :: time_passed
      real(cp),    intent(in) :: time_norm
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)
      logical,     intent(in) :: l_avg

      !-- Local:
      character(len=14) :: string
      character(len=72) :: spec_file
      integer :: n_temp_spec_file
      integer :: n_r,lm,l,m,lc,mc
      real(cp) :: T_temp
      real(cp) :: dT_temp
      real(cp) :: surf_ICB
      real(cp) :: fac,facICB

      real(cp) :: T_r_l(n_r_max,l_max+1),T_r_l_global(n_r_max,l_max+1)
      real(cp) :: T_r_m(n_r_max,l_max+1),T_r_m_global(n_r_max,l_max+1)
      real(cp) :: T_l(l_max+1), T_m(l_max+1)
      real(cp) :: T_ICB_l(l_max+1), T_ICB_l_global(l_max+1)
      real(cp) :: dT_ICB_l(l_max+1), dT_ICB_l_global(l_max+1)
      real(cp) :: T_ICB_m(l_max+1), T_ICB_m_global(l_max+1)
      real(cp) :: dT_ICB_m(l_max+1), dT_ICB_m_global(l_max+1)

      integer :: nOut

      T_l(:)     =0.0_cp
      T_ICB_l(:) =0.0_cp
      dT_ICB_l(:)=0.0_cp
      T_m(:)     =0.0_cp
      T_ICB_m(:) =0.0_cp
      dT_ICB_m(:)=0.0_cp

      do n_r=1,n_r_max
         do l=1,l_max+1
            T_r_l(n_r,l)=0.0_cp
            T_ICB_l(l)  =0.0_cp
            dT_ICB_l(l) =0.0_cp
            T_r_m(n_r,l)=0.0_cp
            T_ICB_m(l)  =0.0_cp
            dT_ICB_m(l) =0.0_cp
         end do
         do lm=llm,ulm
            l =lo_map%lm2l(lm)
            m =lo_map%lm2m(lm)
            lc=l+1
            mc=m+1

            T_temp =sqrt(cc2real(s(lm,n_r),m))/or2(n_r)
            dT_temp=sqrt(cc2real(ds(lm,n_r),m))/or2(n_r)

            !----- l-spectra:
            T_r_l(n_r,lc)=T_r_l(n_r,lc) + T_temp
            !----- m-spectra:
            T_r_m(n_r,mc)=T_r_m(n_r,mc) + T_temp

            !----- ICB spectra:
            if ( n_r == n_r_icb ) then
               T_ICB_l(lc) =T_ICB_l(lc) +T_temp
               T_ICB_m(mc) =T_ICB_m(mc) +T_temp
               dT_ICB_l(lc)=dT_ICB_l(lc)+dT_temp
               dT_ICB_m(mc)=dT_ICB_m(mc)+dT_temp
            end if

         end do    ! do loop over lms in block 
      end do    ! radial grid points 

      ! Reduction over all ranks
      call reduce_radial(T_r_l, T_r_l_global, 0)
      call reduce_radial(T_r_m, T_r_m_global, 0)
      call reduce_radial(T_ICB_l, T_ICB_l_global, 0)
      call reduce_radial(T_ICB_m, T_ICB_m_global, 0)
      call reduce_radial(dT_ICB_l, dT_ICB_l_global, 0)
      call reduce_radial(dT_ICB_m, dT_ICB_m_global, 0)

      if ( coord_r == 0 .and. l_heat ) then
         !-- Radial Integrals:
         surf_ICB=four*pi*r_icb*r_icb
         fac      =one/vol_oc
         facICB   =one/surf_ICB
         do l=1,l_max+1
            T_l(l)=fac*rInt_R(T_r_l_global(:,l),r,rscheme_oc)
            T_ICB_l(l)=facICB*T_ICB_l_global(l)
            dT_ICB_l(l)=facICB*dT_ICB_l_global(l)
         end do
         do m=1,l_max+1 ! Note: counter m is actual order+1
            T_m(m)=fac*rInt_R(T_r_m_global(:,m),r,rscheme_oc)
            T_ICB_m(m)=facICB*T_ICB_m_global(m)
            dT_ICB_m(m)=facICB*dT_ICB_m_global(m)
         end do

         if ( l_avg ) then
            !-- Averaging:
            call T_l_ave%compute(T_l, n_time_ave, time_passed, time_norm)
            call T_ICB_l_ave%compute(T_ICB_l, n_time_ave, time_passed, time_norm)
            call dT_ICB_l_ave%compute(dT_ICB_l, n_time_ave, time_passed, time_norm)
            call T_m_ave%compute(T_m, n_time_ave, time_passed, time_norm)
            call T_ICB_m_ave%compute(T_ICB_m, n_time_ave, time_passed, time_norm)
            call dT_ICB_m_ave%compute(dT_ICB_m, n_time_ave, time_passed, time_norm)

            !-- Output:
            if ( l_stop_time ) then
               
               call T_l_ave%finalize_SD(time_norm)
               call T_ICB_l_ave%finalize_SD(time_norm)
               call dT_ICB_l_ave%finalize_SD(time_norm)
               call T_m_ave%finalize_SD(time_norm)
               call T_ICB_m_ave%finalize_SD(time_norm)
               call dT_ICB_m_ave%finalize_SD(time_norm)

               !------ Output:
               spec_file='T_spec_ave.'//tag
               open(newunit=nOut,file=spec_file,status='unknown')
               do l=1,l_max+1
                  write(nOut,'(2X,1P,I4,12ES16.8)') l-1, T_l_ave%mean(l), &
                  &                T_m_ave%mean(l),  T_ICB_l_ave%mean(l), &
                  &            T_ICB_m_ave%mean(l), dT_ICB_l_ave%mean(l), &
                  &           dT_ICB_m_ave%mean(l),        T_l_ave%SD(l), &
                  &                  T_m_ave%SD(l),    T_ICB_l_ave%SD(l), &
                  &              T_ICB_m_ave%SD(l),   dT_ICB_l_ave%SD(l), &
                  &             dT_ICB_m_ave%SD(l)
               end do
               close(nOut)

               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown', &
                  &    position='append')
               end if
               write(n_log_file,"(/,A,A)")  &
               &    ' ! TIME AVERAGED T/C SPECTRA STORED IN FILE: ', spec_file
               write(n_log_file,"(A,I5)")  &
               &    ' !              No. of averaged spectra: ', n_time_ave
               if ( l_save_out ) close(n_log_file)

            end if

         else ! Just one spectrum

            !-- Output into files:
            write(string, *) n_spec
            spec_file='T_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=n_temp_spec_file, file=spec_file, status='unknown')
            write(n_temp_spec_file,'(1x,''TC spectra at time:'', ES20.12)')  &
            &     time*tScale
            do l=0,l_max
               write(n_temp_spec_file,'(1P,I4,6ES12.4)') l, T_l(l+1), T_m(l+1), &
               &                                    T_ICB_l(l+1), T_ICB_m(l+1), &
               &                                  dT_ICB_l(l+1), dT_ICB_m(l+1)
            end do
            close(n_temp_spec_file)

         end if

      end if

   end subroutine spectrum_temp
!----------------------------------------------------------------------------
   subroutine get_amplitude(time,w,dw,z,b,db,aj)

      !-- Input of variables:
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)

      !-- Output: 
      real(cp) :: e_mag_p_m(0:l_max),e_mag_t_m(0:l_max)
      real(cp) :: e_kin_p_m(0:l_max),e_kin_t_m(0:l_max)

      !-- Local variables:
      integer :: n_r,lm,l,m

      real(cp) :: e_mag_p_temp,e_mag_t_temp
      real(cp) :: e_kin_p_temp,e_kin_t_temp
      real(cp) :: fac_mag,fac_kin

      real(cp) :: e_mag_p_r_m(n_r_max,0:l_max),e_mag_p_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_mag_t_r_m(n_r_max,0:l_max),e_mag_t_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_kin_p_r_m(n_r_max,0:l_max),e_kin_p_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_kin_t_r_m(n_r_max,0:l_max),e_kin_t_r_m_global(n_r_max,0:l_max)

      do n_r=1,n_r_max

         do m=0,l_max
            if ( l_mag ) then
               e_mag_p_r_m(n_r,m)=0.0_cp
               e_mag_t_r_m(n_r,m)=0.0_cp
            end if
            e_kin_p_r_m(n_r,m)=0.0_cp
            e_kin_t_r_m(n_r,m)=0.0_cp
         end do

         do lm=max(llm,2),ulm

            l  =lo_map%lm2l(lm)
            m  =lo_map%lm2m(lm)

            if ( l_mag ) then
               e_mag_p_temp=dLh(st_map%lm2(l,m)) * ( dLh(st_map%lm2(l,m))*    &
               &            or2(n_r)*cc2real(b(lm,n_r),m) + cc2real(db(lm,n_r),m) )
               e_mag_t_temp=dLh(st_map%lm2(l,m))*cc2real(aj(lm,n_r),m)     
            end if

            e_kin_p_temp=orho1(n_r)*dLh(st_map%lm2(l,m)) *  (                 &
            &            dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
            &            cc2real(dw(lm,n_r),m) )
            e_kin_t_temp=orho1(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)

            !----- m-spectra:
            if ( l_mag ) then
               e_mag_p_r_m(n_r,m)=e_mag_p_r_m(n_r,m)+e_mag_p_temp
               e_mag_t_r_m(n_r,m)=e_mag_t_r_m(n_r,m)+e_mag_t_temp
            end if
            e_kin_p_r_m(n_r,m)=e_kin_p_r_m(n_r,m)+e_kin_p_temp
            e_kin_t_r_m(n_r,m)=e_kin_t_r_m(n_r,m)+e_kin_t_temp

         end do    ! do loop over lms in block 

      end do    ! radial grid points 

      if ( l_mag ) then
         call reduce_radial(e_mag_p_r_m, e_mag_p_r_m_global, 0)
         call reduce_radial(e_mag_t_r_m, e_mag_t_r_m_global, 0)
      end if
      call reduce_radial(e_kin_p_r_m, e_kin_p_r_m_global, 0)
      call reduce_radial(e_kin_t_r_m, e_kin_t_r_m_global, 0)

      if ( coord_r == 0 ) then

         !-- Radial Integrals:
         fac_mag=0.5*LFfac*eScale
         fac_kin=0.5*eScale
         do m=0,l_max ! Note: counter m is actual order+1
            if ( l_mag ) then
               e_mag_p_m(m)= fac_mag*rInt_R(e_mag_p_r_m_global(:,m),r,rscheme_oc)
               e_mag_t_m(m)= fac_mag*rInt_R(e_mag_t_r_m_global(:,m),r,rscheme_oc)
            end if
            e_kin_p_m(m)   =fac_kin*rInt_R(e_kin_p_r_m_global(:,m),r,rscheme_oc) 
            e_kin_t_m(m)   =fac_kin*rInt_R(e_kin_t_r_m_global(:,m),r,rscheme_oc)
         end do

         !-- Output
         if ( l_save_out ) then
            open(newunit=n_am_kpol_file,file=am_kpol_file,status='unknown', &
            &    form='unformatted',position='append')
         end if

         write(n_am_kpol_file) time,(e_kin_p_m(m),m=0,m_max_modes)

         if ( l_save_out ) then
            close(n_am_kpol_file)
         end if

         if ( l_save_out ) then
            open(newunit=n_am_ktor_file,file=am_ktor_file,status='unknown', &
            &    form='unformatted',position='append')
         end if

         write(n_am_ktor_file) time,(e_kin_t_m(m),m=0,m_max_modes)

         if ( l_save_out ) then
            close(n_am_ktor_file)
         end if

         if ( l_mag ) then
            if ( l_save_out ) then
               open(newunit=n_am_mpol_file,file=am_mpol_file,status='unknown', &
               &    form='unformatted',position='append')
            end if

            write(n_am_mpol_file) time,(e_mag_p_m(m),m=0,m_max_modes)

            if ( l_save_out ) then
               close(n_am_mpol_file)
            end if

            if ( l_save_out ) then
               open(newunit=n_am_mtor_file,file=am_mtor_file,status='unknown', &
               &    form='unformatted',position='append')
            end if

            write(n_am_mtor_file) time,(e_mag_t_m(m),m=0,m_max_modes)

            if ( l_save_out ) then
               close(n_am_mtor_file)
            end if
         end if

      end if ! l_master_rank
    
   end subroutine get_amplitude 
!------------------------------------------------------------------------------
end module spectra
