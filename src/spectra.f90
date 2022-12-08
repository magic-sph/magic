module spectra
   !
   ! This module handles the computation and the writing of spectra. It handles
   ! both 2-D spectra in (r,l) and (r,m) spaces and usual spectra integrated
   ! over all radii in (l) or (m) spaces.
   !

   use parallel_mod
   use precision_mod
   use communications, only: reduce_radial
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_ic_maxMag, n_r_maxMag, n_r_ic_max, l_max, minc
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: orho1, orho2, r_ic, chebt_ic, r, r_cmb,  &
       &                       rscheme_oc, or2, r_icb, dr_fac_ic
   use physical_parameters, only: LFfac
   use num_param, only: eScale, tScale
   use blocking, only: lo_map, llm, ulm, llmMag, ulmMag
   use logic, only: l_mag, l_anel, l_cond_ic, l_heat, l_save_out, l_chemical_conv, &
       &            l_energy_modes, l_2D_spectra, l_full_sphere
   use output_data, only: tag, m_max_modes
   use useful, only: cc2real, cc22real, logWrite, round_off
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

   type(mean_sd_type) :: Xi_l_ave, Xi_ICB_l_ave, dXi_ICB_l_ave
   type(mean_sd_type) :: Xi_m_ave, Xi_ICB_m_ave, dXi_ICB_m_ave

   integer :: n_am_kpol_file, n_am_ktor_file
   integer :: n_am_mpol_file, n_am_mtor_file
   character(len=72) :: am_kpol_file, am_ktor_file
   character(len=72) :: am_mpol_file, am_mtor_file

   public :: initialize_spectra, spectrum, get_amplitude, finalize_spectra

contains

   subroutine initialize_spectra()
      !
      ! This subroutine allocates the arrays employed to generate spectra.
      !

      if ( l_mag ) then
         call e_mag_p_l_ave%initialize(0,l_max)
         call e_mag_p_m_ave%initialize(0,l_max)
         call e_mag_t_l_ave%initialize(0,l_max)
         call e_mag_t_m_ave%initialize(0,l_max)
         call e_mag_cmb_l_ave%initialize(0,l_max)
         call e_mag_cmb_m_ave%initialize(0,l_max)
      end if

      if ( l_anel ) then
         call u2_p_l_ave%initialize(0,l_max)
         call u2_p_m_ave%initialize(0,l_max)
         call u2_t_l_ave%initialize(0,l_max)
         call u2_t_m_ave%initialize(0,l_max)
      end if

      call e_kin_p_l_ave%initialize(0,l_max)
      call e_kin_p_m_ave%initialize(0,l_max)
      call e_kin_t_l_ave%initialize(0,l_max)
      call e_kin_t_m_ave%initialize(0,l_max)

      if ( l_2D_spectra ) then

         if ( l_mag ) then
            call e_mag_p_r_l_ave%initialize(1,n_r_max,0,l_max,.false.)
            call e_mag_p_r_m_ave%initialize(1,n_r_max,0,l_max,.false.)
            call e_mag_t_r_l_ave%initialize(1,n_r_max,0,l_max,.false.)
            call e_mag_t_r_m_ave%initialize(1,n_r_max,0,l_max,.false.)
         end if

         call e_kin_p_r_l_ave%initialize(1,n_r_max,0,l_max,.false.)
         call e_kin_p_r_m_ave%initialize(1,n_r_max,0,l_max,.false.)
         call e_kin_t_r_l_ave%initialize(1,n_r_max,0,l_max,.false.)
         call e_kin_t_r_m_ave%initialize(1,n_r_max,0,l_max,.false.)

      end if

      if ( l_heat ) then
         call T_l_ave%initialize(0,l_max)
         call T_ICB_l_ave%initialize(0,l_max)
         call dT_ICB_l_ave%initialize(0,l_max)
         call T_m_ave%initialize(0,l_max)
         call T_ICB_m_ave%initialize(0,l_max)
         call dT_ICB_m_ave%initialize(0,l_max)
      end if

      if ( l_chemical_conv ) then
         call Xi_l_ave%initialize(0,l_max)
         call Xi_ICB_l_ave%initialize(0,l_max)
         call dXi_ICB_l_ave%initialize(0,l_max)
         call Xi_m_ave%initialize(0,l_max)
         call Xi_ICB_m_ave%initialize(0,l_max)
         call dXi_ICB_m_ave%initialize(0,l_max)
      end if

      am_kpol_file='am_kin_pol.'//tag
      am_ktor_file='am_kin_tor.'//tag
      am_mpol_file='am_mag_pol.'//tag
      am_mtor_file='am_mag_tor.'//tag

      if ( rank == 0 .and. (.not. l_save_out) ) then
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
   subroutine finalize_spectra()
      !
      ! This subroutine terminates the memory allocation associated with
      ! spectra.
      !

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

      if ( l_chemical_conv ) then
         call Xi_l_ave%finalize()
         call Xi_ICB_l_ave%finalize()
         call dXi_ICB_l_ave%finalize()
         call Xi_m_ave%finalize()
         call Xi_ICB_m_ave%finalize()
         call dXi_ICB_m_ave%finalize()
      end if

      if ( l_heat ) then
         call T_l_ave%finalize()
         call T_ICB_l_ave%finalize()
         call dT_ICB_l_ave%finalize()
         call T_m_ave%finalize()
         call T_ICB_m_ave%finalize()
         call dT_ICB_m_ave%finalize()
      end if

      if ( rank == 0 .and. (.not. l_save_out) ) then
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
              &        time_norm,s,ds,xi,dxi,w,dw,z,b,db,aj,b_ic,db_ic,aj_ic)
      !
      ! This routine handles the computation and the writing of kinetic energy,
      ! magnetic energy and temperture spectra, depending on the field of interest.
      !

      !-- Input of variables:
      integer,     intent(in) :: n_spec     ! number of spectrum/call, file
      integer,     intent(in) :: n_time_ave
      real(cp),    intent(in) :: time
      real(cp),    intent(in) :: time_passed
      real(cp),    intent(in) :: time_norm
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dxi(llm:ulm,n_r_max)
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

      !-- Local arrays
      real(cp) :: e_mag_p_l(0:l_max),e_mag_t_l(0:l_max)
      real(cp) :: e_kin_p_l(0:l_max),e_kin_t_l(0:l_max)
      real(cp) :: e_mag_p_ic_l(0:l_max),e_mag_t_ic_l(0:l_max)
      real(cp) :: u2_p_l(0:l_max),u2_t_l(0:l_max)
      real(cp) :: e_mag_p_m(0:l_max),e_mag_t_m(0:l_max)
      real(cp) :: e_kin_p_m(0:l_max),e_kin_t_m(0:l_max)
      real(cp) :: e_mag_p_ic_m(0:l_max),e_mag_t_ic_m(0:l_max)
      real(cp) :: u2_p_m(0:l_max),u2_t_m(0:l_max)
      real(cp) :: e_mag_cmb_l(0:l_max), e_mag_cmb_m(0:l_max)
      real(cp) :: e_kin_nearSurf_l(0:l_max), e_kin_nearSurf_m(0:l_max)
      real(cp) :: eCMB(0:l_max),eCMB_global(0:l_max)
      real(cp) :: T_l(0:l_max), T_m(0:l_max), T_ICB_l(0:l_max), T_ICB_m(0:l_max)
      real(cp) :: dT_ICB_l(0:l_max), dT_ICB_m(0:l_max)
      real(cp) :: Xi_l(0:l_max), Xi_m(0:l_max), Xi_ICB_l(0:l_max), Xi_ICB_m(0:l_max)
      real(cp) :: dXi_ICB_l(0:l_max), dXi_ICB_m(0:l_max)

      !-- Local variables
      character(len=14) :: string
      character(len=72) :: file_name
      character(len=255) :: message
      integer :: file_handle
      integer :: n_r,lm,ml,l,m,n_const
      real(cp), parameter :: cut=10.0_cp
      real(cp) :: e_mag_p_tmp
      real(cp) :: fac_mag, fac_kin, nearSurfR

      !-- Local 2D spectra
      real(cp) :: e_mag_p_r_l(n_r_max,0:l_max), e_mag_t_r_l(n_r_max,0:l_max)
      real(cp) :: e_kin_p_r_l(n_r_max,0:l_max), e_kin_t_r_l(n_r_max,0:l_max)
      real(cp) :: u2_p_r_l(n_r_max,0:l_max), u2_t_r_l(n_r_max,0:l_max)
      real(cp) :: e_mag_p_r_m(n_r_max,0:l_max), e_mag_t_r_m(n_r_max,0:l_max)
      real(cp) :: e_kin_p_r_m(n_r_max,0:l_max), e_kin_t_r_m(n_r_max,0:l_max)
      real(cp) :: u2_p_r_m(n_r_max,0:l_max), u2_t_r_m(n_r_max,0:l_max)

      fac_mag=half*LFfac*eScale
      fac_kin=half*eScale
      !-- Compute spectra
      call spectrum_vec(w,dw,z,e_kin_p_r_l,e_kin_t_r_l,e_kin_p_r_m,e_kin_t_r_m, &
           &            e_kin_p_l,e_kin_t_l,e_kin_p_m,e_kin_t_m,fac_kin,orho1)

      if ( l_anel ) call spectrum_vec(w,dw,z,u2_p_r_l,u2_t_r_l,u2_p_r_m,u2_t_r_m, &
                         &            u2_p_l,u2_t_l,u2_p_m,u2_t_m,fac_kin,orho2)

      if ( l_heat ) call spectrum_scal(s, ds, T_l, T_m, T_ICB_l, T_ICB_m, &
                         &             dT_ICB_l, dT_ICB_m)

      if ( l_chemical_conv ) call spectrum_scal(xi, dxi, Xi_l, Xi_m, Xi_ICB_l, &
                                  &             Xi_ICB_m, dXi_ICB_l, dXi_ICB_m)

      !-- Compute the kinetic energy spectra near the external boundary
      if ( rank == 0 ) then
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

         !-- Save nearSurf kin energy spectra:
         e_kin_nearSurf_l(:)=fac_kin*e_kin_p_r_l(n_const,:)
         e_kin_nearSurf_m(:)=fac_kin*e_kin_p_r_m(n_const,:)
      end if

      if ( l_mag ) then
         call spectrum_vec(b,db,aj,e_mag_p_r_l,e_mag_t_r_l,e_mag_p_r_m,e_mag_t_r_m,&
              &            e_mag_p_l,e_mag_t_l,e_mag_p_m,e_mag_t_m,fac_mag)
         !-- inner core:
         if ( l_cond_ic ) then
            call spectrum_vec_IC(b_ic,db_ic,aj_ic,e_mag_p_ic_l,e_mag_t_ic_l, &
                 &               e_mag_p_ic_m,e_mag_t_ic_m,fac_mag)
         else
            e_mag_p_ic_l(:)=0.0_cp; e_mag_t_ic_l(:)=0.0_cp
            e_mag_p_ic_m(:)=0.0_cp; e_mag_t_ic_m(:)=0.0_cp
         end if

         !-- Extra outputs specific to magnetic spectra 
         eCMB(:)=0.0_cp
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            if ( l==0 ) cycle
            m=lo_map%lm2m(lm)
            e_mag_p_tmp= real(l*(l+1),cp) * ( real(l*(l+1),cp)*     &
            &            or2(n_r_cmb)*cc2real(b(lm,n_r_cmb),m) +    &
            &            cc2real(db(lm,n_r_cmb),m) )
            if ( m == 0 ) eCMB(l)=e_mag_p_tmp
         end do
         call reduce_radial(eCMB, eCMB_global, 0)

         if ( rank == 0 ) then
            !-- Save CMB energy spectra:
            e_mag_cmb_l(:)=fac_mag*e_mag_p_r_l(n_r_cmb,:)
            e_mag_cmb_m(:)=fac_mag*e_mag_p_r_m(n_r_cmb,:)
         end if
      end if ! l_mag ?


      !-- Averaging:
      if ( l_avg .and. (rank==0) ) then

         call e_kin_p_l_ave%compute(e_kin_p_l, n_time_ave, time_passed, time_norm)
         call e_kin_t_l_ave%compute(e_kin_t_l, n_time_ave, time_passed, time_norm)
         call e_kin_p_m_ave%compute(e_kin_p_m, n_time_ave, time_passed, time_norm)
         call e_kin_t_m_ave%compute(e_kin_t_m, n_time_ave, time_passed, time_norm)

         !-- Averaging of 2D spectra
         if ( l_2D_spectra ) then
            call e_kin_p_r_l_ave%compute(e_kin_p_r_l, n_time_ave, time_passed, time_norm)
            call e_kin_t_r_l_ave%compute(e_kin_t_r_l, n_time_ave, time_passed, time_norm)
            call e_kin_p_r_m_ave%compute(e_kin_p_r_m, n_time_ave, time_passed, time_norm)
            call e_kin_t_r_m_ave%compute(e_kin_t_r_m, n_time_ave, time_passed, time_norm)
         end if

         if ( l_anel ) then
            call u2_p_l_ave%compute(u2_p_l, n_time_ave, time_passed, time_norm)
            call u2_t_l_ave%compute(u2_t_l, n_time_ave, time_passed, time_norm)
            call u2_p_m_ave%compute(u2_p_m, n_time_ave, time_passed, time_norm)
            call u2_t_m_ave%compute(u2_t_m, n_time_ave, time_passed, time_norm)
         end if

         if ( l_heat ) then
            call T_l_ave%compute(T_l, n_time_ave, time_passed, time_norm)
            call T_ICB_l_ave%compute(T_ICB_l, n_time_ave, time_passed, time_norm)
            call dT_ICB_l_ave%compute(dT_ICB_l, n_time_ave, time_passed, time_norm)
            call T_m_ave%compute(T_m, n_time_ave, time_passed, time_norm)
            call T_ICB_m_ave%compute(T_ICB_m, n_time_ave, time_passed, time_norm)
            call dT_ICB_m_ave%compute(dT_ICB_m, n_time_ave, time_passed, time_norm)
         end if

         if ( l_chemical_conv ) then
            call Xi_l_ave%compute(Xi_l, n_time_ave, time_passed, time_norm)
            call Xi_ICB_l_ave%compute(Xi_ICB_l, n_time_ave, time_passed, time_norm)
            call dXi_ICB_l_ave%compute(dXi_ICB_l, n_time_ave, time_passed, time_norm)
            call Xi_m_ave%compute(Xi_m, n_time_ave, time_passed, time_norm)
            call Xi_ICB_m_ave%compute(Xi_ICB_m, n_time_ave, time_passed, time_norm)
            call dXi_ICB_m_ave%compute(dXi_ICB_m, n_time_ave, time_passed, time_norm)
         end if

         if ( l_mag ) then
            call e_mag_p_l_ave%compute(e_mag_p_l, n_time_ave, time_passed, time_norm)
            call e_mag_t_l_ave%compute(e_mag_t_l, n_time_ave, time_passed, time_norm)
            call e_mag_p_m_ave%compute(e_mag_p_m, n_time_ave, time_passed, time_norm)
            call e_mag_t_m_ave%compute(e_mag_t_m, n_time_ave, time_passed, time_norm)
            call e_mag_cmb_l_ave%compute(e_mag_cmb_l, n_time_ave, time_passed, time_norm)
            call e_mag_cmb_m_ave%compute(e_mag_cmb_m, n_time_ave, time_passed, time_norm)

            !-- Averaging of 2D spectra
            if ( l_2D_spectra ) then
               call e_mag_p_r_l_ave%compute(e_mag_p_r_l, n_time_ave, time_passed, time_norm)
               call e_mag_t_r_l_ave%compute(e_mag_t_r_l, n_time_ave, time_passed, time_norm)
               call e_mag_p_r_m_ave%compute(e_mag_p_r_m, n_time_ave, time_passed, time_norm)
               call e_mag_t_r_m_ave%compute(e_mag_t_r_m, n_time_ave, time_passed, time_norm)
            end if
         end if
      end if

      !-- Writing of snapshot spectra
      if ( (n_spec >= 0) .and. (rank==0) ) then
         write(string, *) n_spec

         file_name='kin_spec_'//trim(adjustl(string))//'.'//tag
         open(newunit=file_handle, file=file_name, status='unknown')
         if ( n_spec == 0 ) then
            write(file_handle,'(1x,''Kinetic energy spectra of time averaged field:'')')
         else
            write(file_handle,'(1x,''Kinetic energy spectra at time:'',ES20.12)')  &
            &     time*tScale
         end if
         do ml=0,l_max
            write(file_handle,'(1p,i4,6ES16.8)')                                &
            &     ml,round_off(e_kin_p_l(ml),maxval(e_kin_p_l),cut),            &
            &     round_off(e_kin_p_m(ml),maxval(e_kin_p_m),cut),               &
            &     round_off(e_kin_t_l(ml),maxval(e_kin_t_l),cut),               &
            &     round_off(e_kin_t_m(ml),maxval(e_kin_t_m),cut),               & 
            &     round_off(e_kin_nearSurf_l(ml),maxval(e_kin_nearSurf_l),cut), &
            &     round_off(e_kin_nearSurf_m(ml),maxval(e_kin_nearSurf_m),cut)
         end do
         close(file_handle)

         if ( l_2D_spectra ) then
            file_name='2D_kin_spec_'//trim(adjustl(string))//'.'//tag
            call write_2D_spectra(file_name,e_kin_p_r_l,e_kin_t_r_l, e_kin_p_r_m,   &
                 &                e_kin_t_r_m,fac_kin,time)
         end if

         if ( l_anel ) then
            file_name='u2_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            if ( n_spec == 0 ) then
               write(file_handle,'(1x,''Velocity square spectra of time averaged field:'')')
            else
               write(file_handle,'(1x,''Velocity square spectra at time:'', &
               &          ES20.12)') time*tScale
            end if
            do ml=0,l_max
               write(file_handle,'(1p,i4,4ES16.8)')               &
               &     ml,round_off(u2_p_l(ml),maxval(u2_p_l),cut), &
               &     round_off(u2_p_m(ml),maxval(u2_p_m),cut),    &
               &     round_off(u2_t_l(ml),maxval(u2_t_l),cut),    &
               &     round_off(u2_t_m(ml),maxval(u2_t_m),cut)
            end do
            close(file_handle)
         end if

         if ( l_heat ) then
            file_name='T_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            if ( n_spec == 0 ) then
               write(file_handle,'(1x,''Temperature spectra of time averaged field:'')')
            else
               write(file_handle,'(1x,''Temperature spectra at time:'', ES20.12)')  &
               &     time*tScale
            end if
            do l=0,l_max
               write(file_handle,'(1P,I4,6ES12.4)') l,            &
               &     round_off(T_l(l),maxval(T_l),cut),           &
               &     round_off(T_m(l),maxval(T_m),cut),           &
               &     round_off(T_ICB_l(l),maxval(T_ICB_l),cut),   &
               &     round_off(T_ICB_m(l),maxval(T_ICB_m),cut),   &
               &     round_off(dT_ICB_l(l),maxval(dT_ICB_l),cut), &
               &     round_off(dT_ICB_m(l),maxval(dT_ICB_m),cut)
            end do
            close(file_handle)
         end if

         if ( l_chemical_conv ) then
            file_name='Xi_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            if ( n_spec == 0 ) then
               write(file_handle,'(1x,''Chemical comp. spectra of time averaged field:'')')
            else
               write(file_handle,'(1x,''Chemical comp. spectra at time:'', ES20.12)')  &
               &     time*tScale
            end if
            do l=0,l_max
               write(file_handle,'(1P,I4,6ES12.4)') l,            &
               &     round_off(Xi_l(l),maxval(Xi_l),cut),           &
               &     round_off(Xi_m(l),maxval(Xi_m),cut),           &
               &     round_off(Xi_ICB_l(l),maxval(Xi_ICB_l),cut),   &
               &     round_off(Xi_ICB_m(l),maxval(Xi_ICB_m),cut),   &
               &     round_off(dXi_ICB_l(l),maxval(dXi_ICB_l),cut), &
               &     round_off(dXi_ICB_m(l),maxval(dXi_ICB_m),cut)
            end do
            close(file_handle)
         end if

         if ( l_mag ) then
            file_name='mag_spec_'//trim(adjustl(string))//'.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            if ( n_spec == 0 ) then
               write(file_handle,'(1x,''Magnetic energy spectra of time averaged field:'')')
            else
               write(file_handle,'(1x,''Magnetic energy spectra at time:'', &
               &           ES20.12)') time*tScale
            end if
            do ml=0,l_max
               write(file_handle,'(1p,i4,11ES16.8)')                         &
               &     ml, round_off(e_mag_p_l(ml),maxval(e_mag_p_l),cut),     &
               &     round_off(e_mag_p_m(ml),maxval(e_mag_p_m),cut),         &
               &     round_off(e_mag_t_l(ml),maxval(e_mag_t_l),cut),         &
               &     round_off(e_mag_t_m(ml),maxval(e_mag_t_m),cut),         &
               &     round_off(e_mag_p_ic_l(ml),maxval(e_mag_p_ic_l),cut),   &
               &     round_off(e_mag_p_ic_m(ml),maxval(e_mag_p_ic_m),cut),   &
               &     round_off(e_mag_t_ic_l(ml),maxval(e_mag_t_ic_l),cut),   &
               &     round_off(e_mag_t_ic_m(ml),maxval(e_mag_t_ic_m),cut),   &
               &     round_off(e_mag_cmb_l(ml),maxval(e_mag_cmb_l),cut),     &
               &     round_off(e_mag_cmb_m(ml),maxval(e_mag_cmb_m),cut),     &
               &     round_off(eCMB_global(ml),maxval(eCMB_global),cut)
            end do
            close(file_handle)

            if ( l_2D_spectra ) then
               file_name='2D_mag_spec_'//trim(adjustl(string))//'.'//tag
               call write_2D_spectra(file_name,e_mag_p_r_l,e_mag_t_r_l, e_mag_p_r_m, &
                    &                e_mag_t_r_m, fac_mag, time)
            end if
         end if

      end if

      !-- Writing of time-averaged spectra at the end of the run
      if ( l_avg .and. l_stop_time .and. (rank==0) ) then
         write(message,"(A,I5)") ' !              No. of averaged spectra: ',n_time_ave
         call logWrite(message)

         !-- Terminates computation of time-averaged quantities
         call e_kin_p_l_ave%finalize_SD(time_norm)
         call e_kin_t_l_ave%finalize_SD(time_norm)
         call e_kin_p_m_ave%finalize_SD(time_norm)
         call e_kin_t_m_ave%finalize_SD(time_norm)
         file_name='kin_spec_ave.'//tag
         open(newunit=file_handle, file=file_name, status='unknown')
         do l=0,l_max
            write(file_handle,'(2X,1P,I4,8ES16.8)') l,                             &
            &     round_off(e_kin_p_l_ave%mean(l),maxval(e_kin_p_l_ave%mean),cut), &
            &     round_off(e_kin_p_m_ave%mean(l),maxval(e_kin_p_m_ave%mean),cut), &
            &     round_off(e_kin_t_l_ave%mean(l),maxval(e_kin_t_l_ave%mean),cut), &
            &     round_off(e_kin_t_m_ave%mean(l),maxval(e_kin_t_m_ave%mean),cut), &
            &     round_off(e_kin_p_l_ave%SD(l),maxval(e_kin_p_l_ave%SD),cut),     &
            &     round_off(e_kin_p_m_ave%SD(l),maxval(e_kin_p_m_ave%SD),cut),     &
            &     round_off(e_kin_t_l_ave%SD(l),maxval(e_kin_t_l_ave%SD),cut),     &
            &     round_off(e_kin_t_m_ave%SD(l),maxval(e_kin_t_m_ave%SD),cut)
         end do
         close(file_handle)

         if ( l_2D_spectra ) then
            file_name='2D_kin_spec_ave.'//tag
            call write_2D_spectra(file_name,e_kin_p_r_l_ave%mean,e_kin_t_r_l_ave%mean, &
                 &                e_kin_p_r_m_ave%mean,e_kin_t_r_m_ave%mean,fac_kin)
         end if

         if ( l_anel ) then
            !-- Output: at end of run
            call u2_p_l_ave%finalize_SD(time_norm)
            call u2_t_l_ave%finalize_SD(time_norm)
            call u2_p_m_ave%finalize_SD(time_norm)
            call u2_t_m_ave%finalize_SD(time_norm)
            file_name='u2_spec_ave.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            do l=0,l_max
               write(file_handle,'(2X,1P,I4,8ES16.8)') l,                       &
               &     round_off(u2_p_l_ave%mean(l),maxval(u2_p_l_ave%mean),cut), &
               &     round_off(u2_p_m_ave%mean(l),maxval(u2_p_m_ave%mean),cut), &
               &     round_off(u2_t_l_ave%mean(l),maxval(u2_t_l_ave%mean),cut), &
               &     round_off(u2_t_m_ave%mean(l),maxval(u2_t_m_ave%mean),cut), &
               &     round_off(u2_p_l_ave%SD(l),maxval(u2_p_l_ave%SD),cut),     &
               &     round_off(u2_p_m_ave%SD(l),maxval(u2_p_m_ave%SD),cut),     &
               &     round_off(u2_t_l_ave%SD(l),maxval(u2_t_l_ave%SD),cut),     &
               &     round_off(u2_t_m_ave%SD(l),maxval(u2_t_m_ave%SD),cut)
            end do
            close(file_handle)
         end if

         if ( l_heat ) then
            call T_l_ave%finalize_SD(time_norm)
            call T_ICB_l_ave%finalize_SD(time_norm)
            call dT_ICB_l_ave%finalize_SD(time_norm)
            call T_m_ave%finalize_SD(time_norm)
            call T_ICB_m_ave%finalize_SD(time_norm)
            call dT_ICB_m_ave%finalize_SD(time_norm)

            file_name='T_spec_ave.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            do l=0,l_max
               write(file_handle,'(2X,1P,I4,12ES16.8)') l,                          &
               &     round_off(T_l_ave%mean(l),maxval(T_l_ave%mean),cut),           &
               &     round_off(T_m_ave%mean(l),maxval(T_m_ave%mean),cut),           &
               &     round_off(T_ICB_l_ave%mean(l),maxval(T_ICB_l_ave%mean),cut),   &
               &     round_off(T_ICB_m_ave%mean(l),maxval(T_ICB_m_ave%mean),cut),   &
               &     round_off(dT_ICB_l_ave%mean(l),maxval(dT_ICB_l_ave%mean),cut), &
               &     round_off(dT_ICB_m_ave%mean(l),maxval(dT_ICB_m_ave%mean),cut), &
               &     round_off(T_l_ave%SD(l),maxval(T_l_ave%SD),cut),               &
               &     round_off(T_m_ave%SD(l),maxval(T_m_ave%SD),cut),               &
               &     round_off(T_ICB_l_ave%SD(l),maxval(T_ICB_l_ave%SD),cut),       &
               &     round_off(T_ICB_m_ave%SD(l),maxval(T_ICB_m_ave%SD),cut),       &
               &     round_off(dT_ICB_l_ave%SD(l),maxval(dT_ICB_l_ave%SD),cut),     &
               &     round_off(dT_ICB_m_ave%SD(l),maxval(dT_ICB_m_ave%SD),cut)
            end do
            close(file_handle)
         end if ! l_heat ?

         if ( l_chemical_conv ) then
            call Xi_l_ave%finalize_SD(time_norm)
            call Xi_ICB_l_ave%finalize_SD(time_norm)
            call dXi_ICB_l_ave%finalize_SD(time_norm)
            call Xi_m_ave%finalize_SD(time_norm)
            call Xi_ICB_m_ave%finalize_SD(time_norm)
            call dXi_ICB_m_ave%finalize_SD(time_norm)

            file_name='Xi_spec_ave.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            do l=0,l_max
               write(file_handle,'(2X,1P,I4,12ES16.8)') l,                            &
               &     round_off(Xi_l_ave%mean(l),maxval(Xi_l_ave%mean),cut),           &
               &     round_off(Xi_m_ave%mean(l),maxval(Xi_m_ave%mean),cut),           &
               &     round_off(Xi_ICB_l_ave%mean(l),maxval(Xi_ICB_l_ave%mean),cut),   &
               &     round_off(Xi_ICB_m_ave%mean(l),maxval(Xi_ICB_m_ave%mean),cut),   &
               &     round_off(dXi_ICB_l_ave%mean(l),maxval(dXi_ICB_l_ave%mean),cut), &
               &     round_off(dXi_ICB_m_ave%mean(l),maxval(dXi_ICB_m_ave%mean),cut), &
               &     round_off(Xi_l_ave%SD(l),maxval(Xi_l_ave%SD),cut),               &
               &     round_off(Xi_m_ave%SD(l),maxval(Xi_m_ave%SD),cut),               &
               &     round_off(Xi_ICB_l_ave%SD(l),maxval(Xi_ICB_l_ave%SD),cut),       &
               &     round_off(Xi_ICB_m_ave%SD(l),maxval(Xi_ICB_m_ave%SD),cut),       &
               &     round_off(dXi_ICB_l_ave%SD(l),maxval(dXi_ICB_l_ave%SD),cut),     &
               &     round_off(dXi_ICB_m_ave%SD(l),maxval(dXi_ICB_m_ave%SD),cut)
            end do
            close(file_handle)
         end if ! l_chemical_conv ?

         if ( l_mag ) then
            !-- Output: at end of run
            call e_mag_p_l_ave%finalize_SD(time_norm)
            call e_mag_t_l_ave%finalize_SD(time_norm)
            call e_mag_p_m_ave%finalize_SD(time_norm)
            call e_mag_t_m_ave%finalize_SD(time_norm)
            call e_mag_cmb_l_ave%finalize_SD(time_norm)
            call e_mag_cmb_m_ave%finalize_SD(time_norm)
            file_name='mag_spec_ave.'//tag
            open(newunit=file_handle, file=file_name, status='unknown')
            do l=0,l_max
               write(file_handle,'(2X,1P,I4,12ES16.8)') l,                               &
               &     round_off(e_mag_p_l_ave%mean(l),maxval(e_mag_p_l_ave%mean),cut),    &
               &     round_off(e_mag_p_m_ave%mean(l),maxval(e_mag_p_m_ave%mean),cut),    &
               &     round_off(e_mag_t_l_ave%mean(l),maxval(e_mag_t_l_ave%mean),cut),    &
               &     round_off(e_mag_t_m_ave%mean(l),maxval(e_mag_t_m_ave%mean),cut),    &
               &     round_off(e_mag_cmb_l_ave%mean(l),maxval(e_mag_cmb_l_ave%mean),cut),&
               &     round_off(e_mag_cmb_m_ave%mean(l),maxval(e_mag_cmb_m_ave%mean),cut),&
               &     round_off(e_mag_p_l_ave%SD(l),maxval(e_mag_p_l_ave%SD),cut),        &
               &     round_off(e_mag_p_m_ave%SD(l),maxval(e_mag_p_m_ave%SD),cut),        &
               &     round_off(e_mag_t_l_ave%SD(l),maxval(e_mag_t_l_ave%SD),cut),        &
               &     round_off(e_mag_t_m_ave%SD(l),maxval(e_mag_t_m_ave%SD),cut),        &
               &     round_off(e_mag_cmb_l_ave%SD(l),maxval(e_mag_cmb_l_ave%SD),cut),    &
               &     round_off(e_mag_cmb_m_ave%SD(l),maxval(e_mag_cmb_m_ave%SD),cut)
            end do
            close(file_handle)

            if ( l_2D_spectra ) then
               file_name='2D_mag_spec_ave.'//tag
               call write_2D_spectra(file_name,e_mag_p_r_l_ave%mean,e_mag_t_r_l_ave%mean,&
                    &                e_mag_p_r_m_ave%mean,e_mag_t_r_m_ave%mean, fac_mag)
            end if
         end if

      end if

   end subroutine spectrum
!----------------------------------------------------------------------------
   subroutine write_2D_spectra(file_name,e_p_r_l,e_t_r_l,e_p_r_m,e_t_r_m,fac,time)

      !-- Input variables
      character(len=*), intent(in) :: file_name ! name of the output file
      real(cp), intent(in) :: e_p_r_l(n_r_max,0:l_max) ! Poloidal spectrum (r,l) space
      real(cp), intent(in) :: e_t_r_l(n_r_max,0:l_max) ! Toroidal spectrum (r,l) space
      real(cp), intent(in) :: e_p_r_m(n_r_max,0:l_max) ! Poloidal spectrum (r,m) space
      real(cp), intent(in) :: e_t_r_m(n_r_max,0:l_max) ! Toroidal spectrum (r,m) space
      real(cp), intent(in) :: fac ! normalisation factor
      real(cp), optional, intent(in) :: time ! Time for a snapshot

      !-- Local variables
      integer :: file_handle

      open(newunit=file_handle, file=file_name, status='unknown', form='unformatted')
      if ( present(time) ) then
         write(file_handle) time*tScale,n_r_max,l_max,minc
      else
         write(file_handle) n_r_max, l_max, minc
      end if
      write(file_handle) r
      write(file_handle) fac*e_p_r_l(1:n_r_max,1:l_max)
      write(file_handle) fac*e_p_r_m
      write(file_handle) fac*e_t_r_l(1:n_r_max,1:l_max)
      write(file_handle) fac*e_t_r_m
      close(file_handle)

   end subroutine write_2D_spectra
!----------------------------------------------------------------------------
   subroutine spectrum_scal(scal,dscal,T_l,T_m,T_ICB_l,T_ICB_m,dT_ICB_l,dT_ICB_m)
      !
      ! This routine is used to compute the spectra of one scalar field such
      ! as temperature or chemical composition.
      !

      !-- Input arrays:
      complex(cp), intent(in) :: scal(llm:ulm,n_r_max) ! The scalar field in l,m space
      complex(cp), intent(in) :: dscal(llm:ulm,n_r_max) ! The radial derivative of the scalar field

      !-- Outputs arrays:
      real(cp), intent(out) :: T_l(0:l_max) ! Spectrum as a function of degree l
      real(cp), intent(out) :: T_m(0:l_max) ! Spectrum as a funtion of order m
      real(cp), intent(out) :: T_ICB_l(0:l_max) ! Spectrum at ICB as a function of l
      real(cp), intent(out) :: T_ICB_m(0:l_max) ! Spectrum at ICB as a function of m
      real(cp), intent(out) :: dT_ICB_l(0:l_max) ! Spectrum of radial der. at ICB as a function of l
      real(cp), intent(out) :: dT_ICB_m(0:l_max) ! Spectrum of radial der. at ICB as a function of m

      !-- Local:
      integer :: n_r, lm, l, m
      real(cp) :: T_temp, dT_temp, surf_ICB, fac, facICB
      real(cp) :: T_r_l(n_r_max,0:l_max),T_r_l_global(n_r_max,0:l_max)
      real(cp) :: T_r_m(n_r_max,0:l_max),T_r_m_global(n_r_max,0:l_max)
      real(cp) ::  T_ICB_l_global(0:l_max), dT_ICB_l_global(0:l_max) 
      real(cp) :: T_ICB_m_global(0:l_max), dT_ICB_m_global(0:l_max)

      T_l(:)     =0.0_cp
      T_ICB_l(:) =0.0_cp
      dT_ICB_l(:)=0.0_cp
      T_m(:)     =0.0_cp
      T_ICB_m(:) =0.0_cp
      dT_ICB_m(:)=0.0_cp

      do n_r=1,n_r_max
         T_r_l(n_r,:)=0.0_cp
         T_ICB_l(:)  =0.0_cp
         dT_ICB_l(:) =0.0_cp
         T_r_m(n_r,:)=0.0_cp
         T_ICB_m(:)  =0.0_cp
         dT_ICB_m(:) =0.0_cp
         do lm=llm,ulm
            l =lo_map%lm2l(lm)
            m =lo_map%lm2m(lm)

            T_temp =cc2real(scal(lm,n_r),m)*r(n_r)*r(n_r)
            dT_temp=cc2real(dscal(lm,n_r),m)*r(n_r)*r(n_r)

            !----- l-spectra:
            T_r_l(n_r,l)=T_r_l(n_r,l) + T_temp
            !----- m-spectra:
            T_r_m(n_r,m)=T_r_m(n_r,m) + T_temp

            !----- ICB spectra:
            if ( n_r == n_r_icb ) then
               T_ICB_l(l) =T_ICB_l(l) +T_temp
               T_ICB_m(m) =T_ICB_m(m) +T_temp
               dT_ICB_l(l)=dT_ICB_l(l)+dT_temp
               dT_ICB_m(m)=dT_ICB_m(m)+dT_temp
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

      if ( rank == 0 ) then ! Only rank==0 handles the radial integrals
         !-- Radial Integrals:
         surf_ICB=four*pi*r_icb*r_icb
         fac      =one/vol_oc
         if ( l_full_sphere ) then
            facICB=0.0_cp
         else
            facICB=one/surf_ICB
         end if
         do l=0,l_max
            T_l(l)=fac*rInt_R(T_r_l_global(:,l),r,rscheme_oc)
            T_ICB_l(l)=facICB*T_ICB_l_global(l)
            dT_ICB_l(l)=facICB*dT_ICB_l_global(l)
         end do
         do m=0,l_max
            T_m(m)=fac*rInt_R(T_r_m_global(:,m),r,rscheme_oc)
            T_ICB_m(m)=facICB*T_ICB_m_global(m)
            dT_ICB_m(m)=facICB*dT_ICB_m_global(m)
         end do
      end if

   end subroutine spectrum_scal
!----------------------------------------------------------------------------
   subroutine spectrum_vec(w,dw,z,e_p_r_l,e_t_r_l,e_p_r_m,e_t_r_m, &
              &            e_p_l,e_t_l,e_p_m,e_t_m,fac_scal,fac)
      !
      ! This routine handles the computation of spectra of a solenoidal vector
      ! field.
      !

      complex(cp), intent(in) :: w(llm:ulm,n_r_max) ! Poloidal potential
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)! Radial derivative of poloidal potential
      complex(cp), intent(in) :: z(llm:ulm,n_r_max) ! Toroidal potential
      real(cp), optional, intent(in) :: fac(n_r_max) ! Factor which depends on radius (like density)
      real(cp),    intent(in) :: fac_scal ! Constant factor (like 1/2)

      real(cp), intent(out) :: e_p_r_l(n_r_max,0:l_max) ! Poloidal spectrum (r,l) space
      real(cp), intent(out) :: e_t_r_l(n_r_max,0:l_max) ! Toroidal spectrum (r,l) space
      real(cp), intent(out) :: e_p_r_m(n_r_max,0:l_max) ! Poloidal spectrum (r,m) space
      real(cp), intent(out) :: e_t_r_m(n_r_max,0:l_max) ! Toroidal spectrum (r,m) space
      real(cp), intent(out) :: e_p_l(0:l_max) ! Poloidal spectrum as a function of l
      real(cp), intent(out) :: e_t_l(0:l_max) ! Toroidal spectrum as a function of l
      real(cp), intent(out) :: e_p_m(0:l_max) ! Poloidal spectrum as a function of m
      real(cp), intent(out) :: e_t_m(0:l_max) ! Toroidal spectrum as a function of m

      !-- Local variables
      integer :: n_r, lm, l, m
      real(cp) :: e_p_tmp, e_t_tmp, facR(n_r_max)
      real(cp) :: e_p_r_l_loc(n_r_max,0:l_max), e_t_r_l_loc(n_r_max,0:l_max)
      real(cp) :: e_p_r_m_loc(n_r_max,0:l_max), e_t_r_m_loc(n_r_max,0:l_max)

      !eCMB(:)=0.0_cip
      if ( present(fac) ) then
         facR(:)=fac(:)
      else
         facR(:)=one
      end if

      e_p_r_l_loc(:,:)=0.0_cp; e_t_r_l_loc(:,:)=0.0_cp
      e_p_r_m_loc(:,:)=0.0_cp; e_t_r_m_loc(:,:)=0.0_cp
      e_p_l(:)=0.0_cp; e_t_l(:)=0.0_cp; e_p_m(:)=0.0_cp; e_t_m(:)=0.0_cp

      do n_r=1,n_r_max
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            if ( l == 0 ) cycle
            m=lo_map%lm2m(lm)

            e_p_tmp=facR(n_r)* real(l*(l+1),cp) * ( real(l*(l+1),cp) *     &
            &       or2(n_r)*cc2real(w(lm,n_r),m) + cc2real(dw(lm,n_r),m) )
            e_t_tmp=facR(n_r)*real(l*(l+1),cp)*cc2real(z(lm,n_r),m)

            !----- l-spectra:
            e_p_r_l_loc(n_r,l) = e_p_r_l_loc(n_r,l) + e_p_tmp
            e_t_r_l_loc(n_r,l) = e_t_r_l_loc(n_r,l) + e_t_tmp
            !if ( m == 0 .and. n_r == n_r_cmb ) eCMB(l)=e_mag_p_temp

            !----- m-spectra:
            e_p_r_m_loc(n_r,m) = e_p_r_m_loc(n_r,m) + e_p_tmp
            e_t_r_m_loc(n_r,m) = e_t_r_m_loc(n_r,m) + e_t_tmp
         end do    ! do loop over lms in block
      end do    ! radial grid points

      ! ----------- We need a reduction here ----------------
      call reduce_radial(e_p_r_l_loc, e_p_r_l, 0)
      call reduce_radial(e_t_r_l_loc, e_t_r_l, 0)
      call reduce_radial(e_p_r_m_loc, e_p_r_m, 0)
      call reduce_radial(e_t_r_m_loc, e_t_r_m, 0)
      !call reduce_radial(eCMB, eCMB_global, 0)

      if ( rank == 0 ) then
         !-- Radial integration only computed by rank==0
         do l=1,l_max
            e_p_l(l)=fac_scal*rInt_R(e_p_r_l(:,l),r,rscheme_oc)
            e_t_l(l)=fac_scal*rInt_R(e_t_r_l(:,l),r,rscheme_oc)
         end do
         do m=0,l_max
            e_p_m(m)=fac_scal*rInt_R(e_p_r_m(:,m),r,rscheme_oc)
            e_t_m(m)=fac_scal*rInt_R(e_t_r_m(:,m),r,rscheme_oc)
         end do
      end if

   end subroutine spectrum_vec
!----------------------------------------------------------------------------
   subroutine spectrum_vec_IC(b_ic,db_ic,aj_ic,e_p_l,e_t_l,e_p_m,e_t_m,fac)
      !
      ! This routine handles the computation of spectra of a solenoidal vector
      ! field in the Inner Core.
      !

      complex(cp), intent(in) :: b_ic(llm:ulm,n_r_ic_max) ! Poloidal potential
      complex(cp), intent(in) :: db_ic(llm:ulm,n_r_ic_max)! Radial derivative of poloidal potential
      complex(cp), intent(in) :: aj_ic(llm:ulm,n_r_ic_max) ! Toroidal potential
      real(cp),    intent(in) :: fac ! Constant normalisation factor

      real(cp), intent(out) :: e_p_l(0:l_max) ! Poloidal spectrum as a function of l
      real(cp), intent(out) :: e_t_l(0:l_max) ! Toroidal spectrum as a function of l
      real(cp), intent(out) :: e_p_m(0:l_max) ! Poloidal spectrum as a function of m
      real(cp), intent(out) :: e_t_m(0:l_max) ! Toroidal spectrum as a function of m

      !-- Local variables
      integer :: n_r, lm, l, m
      real(cp) :: e_p_tmp, e_t_tmp, O_r_icb_E_2, r_ratio
      complex(cp) :: r_dr_b
      real(cp) :: e_p_r_l_global(n_r_max,0:l_max), e_t_r_l_global(n_r_max,0:l_max)
      real(cp) :: e_p_r_m_global(n_r_max,0:l_max), e_t_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_p_r_l(n_r_max,0:l_max), e_t_r_l(n_r_max,0:l_max)
      real(cp) :: e_p_r_m(n_r_max,0:l_max), e_t_r_m(n_r_max,0:l_max)

      O_r_icb_E_2=one/(r_ic(1)*r_ic(1))
      e_p_r_l(:,:)=0.0_cp; e_t_r_l(:,:)=0.0_cp
      e_p_r_m(:,:)=0.0_cp; e_t_r_m(:,:)=0.0_cp
      do n_r=1,n_r_ic_max
         r_ratio=r_ic(n_r)/r_ic(1)
         do lm=llm,ulm
            l =lo_map%lm2l(lm)
            if ( l == 0 ) cycle
            m =lo_map%lm2m(lm)
            r_dr_b=r_ic(n_r)*db_ic(lm,n_r)

            e_p_tmp=real(l*(l+1),cp)*O_r_icb_E_2*r_ratio**(2*l) * (  &
            &       real((2*l+1)*(l+1),cp)*cc2real(b_ic(lm,n_r),m)   +  &
            &       real(2*(l+1),cp)*cc22real(b_ic(lm,n_r),r_dr_b,m) +  &
            &       cc2real(r_dr_b,m) )
            e_t_tmp=real(l*(l+1),cp)*r_ratio**(2*l+2)*cc2real(aj_ic(lm,n_r),m)

            e_p_r_l(n_r,l)=e_p_r_l(n_r,l) + e_p_tmp
            e_t_r_l(n_r,l)=e_t_r_l(n_r,l) + e_t_tmp
            e_p_r_m(n_r,m)=e_p_r_m(n_r,m) + e_p_tmp
            e_t_r_m(n_r,m)=e_t_r_m(n_r,m) + e_t_tmp
         end do  ! loop over lm's
      end do ! loop over radial levels

      call reduce_radial(e_p_r_l, e_p_r_l_global, 0)
      call reduce_radial(e_t_r_l, e_t_r_l_global, 0)
      call reduce_radial(e_p_r_m, e_p_r_m_global, 0)
      call reduce_radial(e_t_r_m, e_t_r_m_global, 0)

      if ( rank == 0 ) then
         !----- Radial Integrals:
         do l=1,l_max
            e_p_l(l)=fac*rIntIC(e_p_r_l_global(1,l),n_r_ic_max,dr_fac_ic,chebt_ic)
            e_t_l(l)=fac*rIntIC(e_t_r_l_global(1,l),n_r_ic_max,dr_fac_ic,chebt_ic)
         end do
         do m=0,l_max
            e_p_m(m)=fac*rIntIC(e_p_r_m_global(1,m),n_r_ic_max,dr_fac_ic,chebt_ic)
            e_t_m(m)=fac*rIntIC(e_t_r_m_global(1,m),n_r_ic_max,dr_fac_ic,chebt_ic)
         end do
      end if

   end subroutine spectrum_vec_IC
!----------------------------------------------------------------------------
   subroutine get_amplitude(time,w,dw,z,b,db,aj)
      !
      ! This routine is used to generate times series of magnetic and kinetic
      ! energy spectra as a function of the spherical harmonic order m.
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)

      !-- Output:
      real(cp) :: e_mag_p_m(0:l_max), e_mag_t_m(0:l_max)
      real(cp) :: e_kin_p_m(0:l_max), e_kin_t_m(0:l_max)

      !-- Local variables:
      integer :: n_r, lm, l, m
      real(cp) :: e_mag_p_tmp, e_mag_t_tmp, e_kin_p_tmp, e_kin_t_tmp
      real(cp) :: fac_mag, fac_kin
      real(cp) :: e_mag_p_r_m(n_r_max,0:l_max),e_mag_p_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_mag_t_r_m(n_r_max,0:l_max),e_mag_t_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_kin_p_r_m(n_r_max,0:l_max),e_kin_p_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_kin_t_r_m(n_r_max,0:l_max),e_kin_t_r_m_global(n_r_max,0:l_max)

      e_kin_p_r_m(:,:)=0.0_cp
      e_kin_t_r_m(:,:)=0.0_cp
      if ( l_mag ) then
         e_mag_p_r_m(:,:)=0.0_cp
         e_mag_t_r_m(:,:)=0.0_cp
      end if
      do n_r=1,n_r_max
         do lm=llm,ulm
            l  =lo_map%lm2l(lm)
            if ( l == 0 ) cycle
            m  =lo_map%lm2m(lm)

            if ( l_mag ) then
               e_mag_p_tmp=real(l*(l+1),cp) * ( real(l*(l+1),cp)*    &
               &           or2(n_r)*cc2real(b(lm,n_r),m) + cc2real(db(lm,n_r),m) )
               e_mag_t_tmp=real(l*(l+1),cp)*cc2real(aj(lm,n_r),m)
            end if

            e_kin_p_tmp=orho1(n_r)*real(l*(l+1),cp) *  (                 &
            &           real(l*(l+1),cp)*or2(n_r)*cc2real(w(lm,n_r),m) + &
            &           cc2real(dw(lm,n_r),m) )
            e_kin_t_tmp=orho1(n_r)*real(l*(l+1),cp)*cc2real(z(lm,n_r),m)

            !----- m-spectra:
            if ( l_mag ) then
               e_mag_p_r_m(n_r,m)=e_mag_p_r_m(n_r,m)+e_mag_p_tmp
               e_mag_t_r_m(n_r,m)=e_mag_t_r_m(n_r,m)+e_mag_t_tmp
            end if
            e_kin_p_r_m(n_r,m)=e_kin_p_r_m(n_r,m)+e_kin_p_tmp
            e_kin_t_r_m(n_r,m)=e_kin_t_r_m(n_r,m)+e_kin_t_tmp

         end do    ! do loop over lms in block
      end do    ! radial grid points

      if ( l_mag ) then
         call reduce_radial(e_mag_p_r_m, e_mag_p_r_m_global, 0)
         call reduce_radial(e_mag_t_r_m, e_mag_t_r_m_global, 0)
      end if
      call reduce_radial(e_kin_p_r_m, e_kin_p_r_m_global, 0)
      call reduce_radial(e_kin_t_r_m, e_kin_t_r_m_global, 0)

      if ( rank == 0 ) then

         !-- Radial Integrals:
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

      end if ! rank == 0

   end subroutine get_amplitude
!------------------------------------------------------------------------------
end module spectra
