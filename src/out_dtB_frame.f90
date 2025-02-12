module out_dtB_frame
   !
   ! This module handles the final assembly of the arrays produced to store
   ! the movies that diagnose the different contributions of the induction
   ! equations. It produces arrays that are stored in the corresponding
   ! movie files.
   !

   use precision_mod
   use parallel_mod
   use truncation, only: n_r_max, lm_max, n_r_ic_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag, l_max, m_max, minc, n_theta_max,     &
       &                 n_r_tot, nlat_padded, n_phi_max, n_cheb_ic_max
   use communications, only: lo2r_one, gather_from_lo_to_rank0
   use fields, only: work_LMloc, b_Rloc, db_Rloc, aj_Rloc, dj_Rloc, b_LMloc, &
       &             b_ic_LMloc, aj_ic_LMloc, db_ic_LMloc, dj_ic_LMloc
   use radial_data, only: nRstart, nRstop
   use radial_functions, only: r, or1, chebt_ic, r_ic, rscheme_oc, r_icb, &
       &                       dr_fac_ic, chebt_ic_even, l_R
   use blocking, only: lm2l, lm2m, lm2, lo_map, llmMag, ulmMag, llm, ulm
   use horizontal_data, only: cosTheta, n_theta_cal2ord, sinTheta, &
       &                      dLh, O_sin_theta
   use dtB_mod, only: PadvLMIC_LMloc, PdifLMIC_LMloc, TadvLMIC_LMloc,          &
       &              TdifLMIC_LMloc, PstrLM_LMloc, PstrLM_Rloc, PdifLM_LMloc, &
       &              PadvLM_Rloc, PadvLM_LMloc, TstrLM_Rloc, TstrLM_LMloc,    &
       &              TadvLM_Rloc, TadvLM_LMloc, TdifLM_LMloc, TomeLM_Rloc,    &
       &              TomeLM_LMloc
   use movie_data, only: n_movie_fields, n_movie_fields_ic, n_movie_file,   &
       &                 n_movie_const, n_movie_surface, movie_const,       &
       &                 n_movie_field_type, frames, n_movies, n_movie_field_start
   use logic, only: l_cond_ic
   use constants, only: zero, one
   use radial_der_even, only: get_drNS_even
   use radial_der, only: get_dr
   use sht, only: torpol_to_spat_single, sphtor_to_spat, scal_to_spat, &
       &          toraxi_to_spat, sht_l_single

   implicit none

   private

   public :: calc_dtB_frame, calc_dtB_frame_IC

contains

   subroutine calc_dtB_frame()
      !
      !  Controls output of specific movie frames related to magnetic field
      !  production and diffusion.
      !

      !-- Local variables:
      integer :: n_out, n_movie
      integer :: n_fields, n_field, n_const, n_r, n_store_last
      integer :: n_theta, n_theta_cal, n_phi, n_o, n_or
      integer :: lm, l, m
      integer :: n_surface, n_field_type

      real(cp) :: dtB(n_theta_max)
      real(cp) :: dtBr(nlat_padded,n_phi_max),dtBt(nlat_padded,n_phi_max)
      real(cp) :: dtBp(nlat_padded,n_phi_max)
      real(cp) :: const,rMov

      complex(cp) :: tmp(lm_max,nRstart:nRstop), work_Rloc(lm_max,nRstart:nRstop)

      do n_movie=1,n_movies

         n_fields =n_movie_fields(n_movie)
         n_out    =n_movie_file(n_movie)
         n_const  =n_movie_const(n_movie)
         n_surface=n_movie_surface(n_movie)
         const    =movie_const(n_movie)

         do n_field=1,n_fields

            n_field_type=n_movie_field_type(n_field,n_movie)
            n_store_last=n_movie_field_start(n_field,n_movie)-1

            if ( n_surface == 4 ) then ! Axisymmetric movies

               if ( n_field_type == 20 ) then
                  call lo2r_one%transp_lm2r(PstrLM_LMloc,PstrLM_Rloc)
               else if ( n_field_type == 21 ) then
                  call lo2r_one%transp_lm2r(PadvLM_LMloc,PadvLM_Rloc)
               else if ( n_field_type == 22 ) then
                  call lo2r_one%transp_lm2r(PdifLM_LMloc,work_Rloc)
               else if ( n_field_type == 23 ) then
                  call lo2r_one%transp_lm2r(TstrLM_LMloc,TstrLM_Rloc)
               else if ( n_field_type == 25 ) then
                  call lo2r_one%transp_lm2r(TadvLM_LMloc,TadvLM_Rloc)
               else if ( n_field_type == 26 ) then
                  call lo2r_one%transp_lm2r(TdifLM_LMloc,work_Rloc)
               else
                  cycle ! Leave the n_field loop
               end if

               do n_r=nRstart,nRstop

                  if ( n_field_type == 20 ) then
                     call get_dtB(dtB,PstrLM_Rloc,1,lm_max,nRstart,nRstop,n_r,.false.)
                  else if ( n_field_type == 21 ) then
                     call get_dtB(dtB,PadvLM_Rloc,1,lm_max,nRstart,nRstop,n_r,.false.)
                  else if ( n_field_type == 22 ) then
                     call get_dtB(dtB,work_Rloc,1,lm_max,nRstart,nRstop,n_r,.false.)
                  else if ( n_field_type == 23 ) then
                     call get_dtB(dtB,TstrLM_Rloc,1,lm_max,nRstart,nRstop,n_r,.false.)
                  else if ( n_field_type == 25 ) then
                     call get_dtB(dtB,TadvLM_Rloc,1,lm_max,nRstart,nRstop,n_r,.false.)
                  else if ( n_field_type == 26 ) then
                     call get_dtB(dtB,work_Rloc,1,lm_max,nRstart,nRstop,n_r,.false.)
                  end if

                  n_or=n_store_last+(n_r-1)*n_theta_max
                  !--- Store in frame field:
                  do n_theta_cal=1,n_theta_max
                     n_theta=n_theta_cal2ord(n_theta_cal)
                     frames(n_or+n_theta)=dtB(n_theta_cal)
                  end do

               end do

            else if ( n_surface == 0 .or. n_surface == 1 ) then ! 3D or rcut movies

               if ( n_field_type == 24 ) then
                  ! This reduces omega effect to field production of axisymm. toroidal field:
                  do n_r=1,n_r_max
                     do lm=llm,ulm
                        m=lo_map%lm2m(lm)
                        if ( m == 0 ) then
                           work_LMloc(lm,n_r)=TomeLM_LMloc(lm,n_r)
                        else
                           work_LMloc(lm,n_r)=zero
                        end if
                     end do
                  end do
               else if ( n_field_type == 81 ) then
                  ! This reduces poloidal field to the dipole contribution:
                  do n_r=1,n_r_max
                     do lm=llm,ulm
                        l=lo_map%lm2l(lm)
                        if ( l == 1 ) then
                           work_LMloc(lm,n_r)=b_LMloc(lm,n_r)
                        else
                           work_LMloc(lm,n_r)=zero
                        end if
                     end do
                  end do
               end if

               !------ Outer core contribution:
               !------ Calculate needed radial derivatives:
               if ( n_field_type == 35 ) then
                  call get_dr(PstrLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               else if ( n_field_type == 36 ) then
                  call get_dr(PadvLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               else if ( n_field_type == 37 ) then
                  call get_dr(PdifLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               else if ( n_field_type == 38 ) then
                  call get_dr(TstrLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               else if ( n_field_type == 39 ) then
                  call get_dr(TomeLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               else if ( n_field_type == 40 ) then
                  call get_dr(TadvLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               else if ( n_field_type == 41 ) then
                  call get_dr(TdifLM_LMloc,work_LMloc,ulm-llm+1,1,ulm-llm+1, &
                       &      n_r_max,rscheme_oc,nocopy=.true.)
               end if

               !-- Br:
               if ( n_field_type == 27 ) then
                  call lo2r_one%transp_lm2r(PstrLM_LMloc,PstrLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 28 ) then
                  call lo2r_one%transp_lm2r(PadvLM_LMloc,PadvLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 29 ) then
                  call lo2r_one%transp_lm2r(PdifLM_LMloc,tmp)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 81 ) then
                  call lo2r_one%transp_lm2r(PdifLM_LMloc,tmp)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               !-- Jr:
               else if ( n_field_type == 30 ) then
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 31 ) then
                  call lo2r_one%transp_lm2r(TstrLM_LMloc,TstrLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 32 ) then
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 33 ) then
                  call lo2r_one%transp_lm2r(TadvLM_LMloc,TadvLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 34 ) then
                  call lo2r_one%transp_lm2r(TdifLM_LMloc,tmp)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               !-- Bz poloidal
               else if ( n_field_type == 35 ) then
                  call lo2r_one%transp_lm2r(PstrLM_LMloc,PstrLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 36 ) then
                  call lo2r_one%transp_lm2r(PadvLM_LMloc,PadvLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 37 ) then
                  call lo2r_one%transp_lm2r(PdifLM_LMloc,tmp)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               !-- Jz poloidal
               else if ( n_field_type == 38 ) then
                  call lo2r_one%transp_lm2r(TstrLM_LMloc,TstrLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 39 ) then
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 40 ) then
                  call lo2r_one%transp_lm2r(TadvLM_LMloc,TadvLM_Rloc)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 41 ) then
                  call lo2r_one%transp_lm2r(TdifLM_LMloc,tmp)
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 24 ) then
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 43 ) then
                  call lo2r_one%transp_lm2r(TstrLM_LMloc,TstrLM_Rloc)
               else if ( n_field_type == 44 ) then
                  call lo2r_one%transp_lm2r(work_LMloc,work_Rloc)
               else if ( n_field_type == 45 ) then
                  call lo2r_one%transp_lm2r(TadvLM_LMloc,TadvLM_Rloc)
               else if ( n_field_type == 46 ) then
                  call lo2r_one%transp_lm2r(TdifLM_LMloc,tmp)
               else if ( n_field_type == 49 ) then
                  call lo2r_one%transp_lm2r(TomeLM_LMloc,TomeLM_Rloc)
               end if

               do n_r=nRstart,nRstop
                  if ( n_surface == 0 ) then
                     rMov=r(n_r)
                     n_or=n_store_last+(n_r-1)*n_theta_max*n_phi_max
                  else if ( n_surface == 1 ) then
                     rMov=const
                     n_or=n_store_last
                     if ( n_r /= n_const ) cycle  ! not desired radius
                  end if

                  !-- Br:
                  if ( n_field_type == 27 ) then
                     call get_Bpol(PstrLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov, &
                          &        .false.)
                  else if ( n_field_type == 28 ) then
                     call get_Bpol(PadvLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 29 ) then
                     call get_Bpol(tmp(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov,.false.)
                  else if ( n_field_type == 81 ) then
                     !---------- get radial field diffusion and radial dipole field:
                     call get_Bpol(tmp(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                     call get_Bpol(work_Rloc(:,n_r),work_Rloc(:,n_r),dtBt,dtBp,dtBp,  &
                          &        rMov,.false.)
                  !-- Jr:
                  else if ( n_field_type == 30 ) then
                     call get_Bpol(aj_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp,  &
                          &        rMov,.false.)
                  else if ( n_field_type == 31 ) then
                     call get_Bpol(TstrLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 32 ) then
                     call get_Bpol(work_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 33 ) then
                     call get_Bpol(TadvLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 34 ) then
                     call get_Bpol(tmp(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov,.false.)
                  !-- Bz poloidal
                  else if ( n_field_type == 13 ) then
                     call get_Bpol(b_Rloc(:,n_r),db_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov, &
                          &        .false.)
                  else if ( n_field_type == 35 ) then
                     call get_Bpol(PstrLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 36 ) then
                     call get_Bpol(PadvLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 37 ) then
                     call get_Bpol(tmp(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov,.false.)

                  !-- Jz poloidal
                  else if ( n_field_type == 14 ) then
                     call get_Bpol(aj_Rloc(:,n_r),dj_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov, &
                          &        .false.)
                  else if ( n_field_type == 38 ) then
                     call get_Bpol(TstrLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 39 ) then
                     call get_Bpol(work_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 40 ) then
                     call get_Bpol(TadvLM_Rloc(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.false.)
                  else if ( n_field_type == 41 ) then
                     call get_Bpol(tmp(:,n_r),work_Rloc(:,n_r),dtBr,dtBt,dtBp,rMov,.false.)
                  !--- Bp toroidal:
                  else if ( n_field_type == 24 ) then
                     call get_Btor(work_Rloc(:,n_r),dtBt,dtBr,rMov,.false.)
                  else if ( n_field_type == 42 ) then
                     call get_Btor(aj_Rloc(:,n_r),dtBt,dtBr,rMov,.false.)
                  else if ( n_field_type == 43 ) then
                     call get_Btor(TstrLM_Rloc(:,n_r),dtBt,dtBr,rMov,.false.)
                  else if ( n_field_type == 44 ) then
                     call get_Btor(work_Rloc(:,n_r),dtBt,dtBr,rMov,.false.)
                  else if ( n_field_type == 45 ) then
                     call get_Btor(TadvLM_Rloc(:,n_r),dtBt,dtBr,rMov,.false.)
                  else if ( n_field_type == 46 ) then
                     call get_Btor(tmp(:,n_r),dtBt,dtBr,rMov,.false.)
                  else if ( n_field_type == 49 ) then
                     call get_Btor(TomeLM_Rloc(:,n_r),dtBt,dtBr,rMov,.false.)
                  !--- Bt toroidal
                  else if ( n_field_type == 50 ) then
                     call get_Btor(aj_Rloc(:,n_r),dtBr,dtBt,rMov,.false.)
                  !--- Toroidal Potential
                  else if ( n_field_type == 51 ) then
                     call lm2pt(aj_Rloc(:,n_r),dtBr,rMov,.false.,.false.)
                  !--- Fieldlines for theta=const.
                  else if ( n_field_type == 52 ) then
                     call get_Btor(b_Rloc(:,n_r),dtBr,dtBt,rMov,.false.)
                  else
                     cycle ! Do not store anything in frames(*)
                  end if

                  !--- Now store the stuff to frames(*):
                  do n_phi=1,n_phi_max
                     do n_theta_cal=1,n_theta_max
                        n_theta=n_theta_cal2ord(n_theta_cal)
                        n_o=n_or+(n_theta-1)*n_phi_max
                        if ( n_field_type == 13 .or. n_field_type == 35 .or. &
                        &    n_field_type == 36 .or. n_field_type == 37 .or. &
                        &    n_field_type == 14 .or. n_field_type == 38 .or. &
                        &    n_field_type == 39 .or. n_field_type == 40 .or. &
                        &    n_field_type == 41 ) then
                           frames(n_phi+n_o) =                                  &
                           &    cosTheta(n_theta_cal)*dtBr(n_theta_cal,n_phi) - &
                           &    sinTheta(n_theta_cal)*dtBt(n_theta_cal,n_phi)
                        else if ( n_field_type == 81 ) then
                           if ( dtBr(n_theta_cal,n_phi) * &
                                dtBt(n_theta_cal,n_phi) < 0.0_cp ) then
                              frames(n_phi+n_o)=dtBr(n_theta_cal,n_phi)
                           else
                              frames(n_phi+n_o)=0.0_cp
                           end if
                        else
                           frames(n_phi+n_o)=dtBr(n_theta_cal,n_phi)
                        end if ! Br (Bp) or Bz ?
                     end do ! phi Loop
                  end do    ! theta loop

               end do          ! r loop

            end if ! which surface ?

         end do  ! LOOP over fields in movie file

      end do ! Loop over movies

   end subroutine calc_dtB_frame
!----------------------------------------------------------------------------
   subroutine calc_dtB_frame_IC()
      !
      !  Controls output of specific movie frames related to magnetic field
      !  production and diffusion in the inner core.
      !

      !-- Local variables:
      integer :: n_out, n_fields_oc, n_fields_ic, n_movie
      integer :: n_fields, n_field, n_const, n_r, n_store_last
      integer :: n_theta, n_theta_cal, n_phi, n_o, n_or
      integer :: n_surface, n_field_type

      real(cp) :: dtB(n_theta_max)
      real(cp) :: dtBr(nlat_padded,n_phi_max),dtBt(nlat_padded,n_phi_max)
      real(cp) :: dtBp(nlat_padded,n_phi_max)
      real(cp) :: const,rMov

      complex(cp) :: workA(llm:ulm,n_r_ic_max), workB(llm:ulm,n_r_ic_max)
      complex(cp) :: workA_1d(lm_max), workB_1d(lm_max)

      logical :: l_loop

      do n_movie=1,n_movies

         n_fields_oc =n_movie_fields(n_movie)
         n_fields_ic =n_movie_fields_ic(n_movie)
         n_fields    =n_fields_oc+n_fields_ic
         n_out       =n_movie_file(n_movie)
         n_const     =n_movie_const(n_movie)
         n_surface   =n_movie_surface(n_movie)
         const       =movie_const(n_movie)

         do n_field=n_fields_oc+1,n_fields

            n_field_type=n_movie_field_type(n_field,n_movie)
            n_store_last=n_movie_field_start(n_field,n_movie)-1

            !--- Axisymmetric dtFL or dtAB:
            if ( n_surface == 4 ) then

               do n_r=1,n_r_ic_max

                  if ( n_field_type == 21 ) then
                     call get_dtB(dtB,PadvLMIC_LMloc,llm,ulm,1,n_r_ic_max,n_r,.true.)
                  else if ( n_field_type == 22 ) then
                     call get_dtB(dtB,PdifLMIC_LMloc,llm,ulm,1,n_r_ic_max,n_r,.true.)
                  else if ( n_field_type == 25 ) then
                     call get_dtB(dtB,TadvLMIC_LMloc,llm,ulm,1,n_r_ic_max,n_r,.true.)
                  else if ( n_field_type == 26 ) then
                     call get_dtB(dtB,TdifLMIC_LMloc,llm,ulm,1,n_r_ic_max,n_r,.true.)
                  else
                     cycle
                  end if


                  if ( rank == 0 ) then
                     n_or=n_store_last+(n_r-1)*n_theta_max
                     !--- Store in frame field:
                     do n_theta_cal=1,n_theta_max
                        n_theta=n_theta_cal2ord(n_theta_cal)
                        frames(n_or+n_theta)=dtB(n_theta_cal)
                     end do
                  end if

               end do

            !--- dtBr:
            else if ( n_surface == 0 .or. n_surface == 1 ) then

               l_loop=.true.
               if ( n_surface == 1 .and. const >= r_icb ) l_loop= .false.

               if ( l_loop ) then
                  !------ Calculate needed radial derivatives:
                  if ( n_field_type == 36 ) then
                     call get_drNS_even(PadvLMIC_LMloc,workA,ulm-llm+1,1,ulm-llm+1, &
                          &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB,   &
                          &             chebt_ic,chebt_ic_even)
                  else if ( n_field_type == 37 ) then
                     call get_drNS_even(PdifLMIC_LMloc,workA,ulm-llm+1,1,ulm-llm+1, &
                          &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB,   &
                          &             chebt_ic,chebt_ic_even)
                  else if ( n_field_type == 40 ) then
                     call get_drNS_even(TadvLMIC_LMloc,workA,ulm-llm+1,1,ulm-llm+1, &
                          &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB,   &
                          &             chebt_ic,chebt_ic_even)
                  else if ( n_field_type == 41 ) then
                     call get_drNS_even(TdifLMIC_LMloc,workA,ulm-llm+1,1,ulm-llm+1, &
                          &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB,   &
                          &             chebt_ic,chebt_ic_even)
                  else
                     workA(:,:)=zero
                  end if

                  do n_r=1,n_r_ic_max
                     if ( n_surface == 0 ) then
                        rMov=r_ic(n_r)
                        n_or=n_store_last+(n_r-1)*n_theta_max*n_phi_max
                     else if ( n_surface == 1 ) then
                        rMov=const
                        n_or=n_store_last
                        if ( n_r /= n_const ) cycle  ! not desired radius
                     end if

                     if ( n_field_type == 35 .or. n_field_type == 43 .or. &
                     &    n_field_type == 44 .or. n_field_type == 27 .or. &
                     &    n_field_type == 31 .or. n_field_type == 32 .or. &
                     &    n_field_type == 38 .or. n_field_type == 49 .or. &
                     &    n_field_type == 39 ) then
                        dtBr(:,:)=0.0_cp
                        dtBt(:,:)=0.0_cp
                        dtBp(:,:)=0.0_cp
                     !--- Br:
                     else if ( n_field_type == 28 ) then
                        call gather_from_lo_to_rank0(PadvLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 29 ) then
                        call gather_from_lo_to_rank0(PdifLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                        !--- Jr:
                     else if ( n_field_type == 30 ) then
                        call gather_from_lo_to_rank0(aj_ic_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 33 ) then
                        call gather_from_lo_to_rank0(TadvLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 34 ) then
                        call gather_from_lo_to_rank0(TdifLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     !- Bz poloidal:
                     else if ( n_field_type == 13 ) then
                        call gather_from_lo_to_rank0(b_ic_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 36 ) then
                        call gather_from_lo_to_rank0(PadvLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 37 ) then
                        call gather_from_lo_to_rank0(PdifLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     !--- Jz poloidal:
                     else if ( n_field_type == 14 ) then
                        call gather_from_lo_to_rank0(aj_ic_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(dj_ic_LMloc(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 40 ) then
                        call gather_from_lo_to_rank0(TadvLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     else if ( n_field_type == 41 ) then
                        call gather_from_lo_to_rank0(TdifLMIC_LMloc(:,n_r), workA_1d)
                        call gather_from_lo_to_rank0(workA(:,n_r), workB_1d)
                        call get_Bpol(workA_1d,workB_1d,dtBr,dtBt,dtBp,rMov,.true.)
                     !--- Bphi toroidal:
                     else if ( n_field_type == 42 ) then
                        call gather_from_lo_to_rank0(aj_ic_LMloc(:,n_r), workA_1d)
                        call get_Btor(workA_1d,dtBt,dtBr,rMov,.true.)
                     else if ( n_field_type == 45 ) then
                        call gather_from_lo_to_rank0(TadvLMIC_LMloc(:,n_r), workA_1d)
                        call get_Btor(workA_1d,dtBt,dtBr,rMov,.true.)
                     else if ( n_field_type == 46 ) then
                        call gather_from_lo_to_rank0(TdifLMIC_LMloc(:,n_r), workA_1d)
                        call get_Btor(workA_1d,dtBt,dtBr,rMov,.true.)
                     !--- Btheta toroidal:
                     else if ( n_field_type == 50 ) then
                        call gather_from_lo_to_rank0(aj_ic_LMloc(:,n_r), workA_1d)
                        call get_Btor(workA_1d,dtBr,dtBt,rMov,.true.)
                     !--- Toroidal Potential
                     else if ( n_field_type == 51 ) then
                        call gather_from_lo_to_rank0(aj_ic_LMloc(:,n_r), workA_1d)
                        call lm2pt(workA_1d,dtBr,r_ic(n_r),.true.,.false.)
                     !--- Fieldlines for theta=const.
                     else if ( n_field_type == 52 ) then
                        call gather_from_lo_to_rank0(b_ic_LMloc(:,n_r), workA_1d)
                        call get_Btor(workA_1d,dtBr,dtBt,rMov,.true.)
                     else
                        cycle
                     end if

                     !--- Now rank=0 store the stuff for the inner core.
                     if ( rank == 0 ) then
                        do n_phi=1,n_phi_max
                           do n_theta_cal=1,n_theta_max
                              n_theta=n_theta_cal2ord(n_theta_cal)
                              n_o=n_or+(n_theta-1)*n_phi_max
                              if ( n_field_type == 13 .or. n_field_type == 36 .or. &
                              &    n_field_type == 37 .or. n_field_type == 14 .or. &
                              &    n_field_type == 40 .or. n_field_type == 41 ) then
                                 frames(n_phi+n_o) =                                 &
                                 &    cosTheta(n_theta_cal)*dtBr(n_theta_cal,n_phi) -&
                                 &    sinTheta(n_theta_cal)*dtBt(n_theta_cal,n_phi)
                              else
                                 frames(n_phi+n_o)=dtBr(n_theta_cal,n_phi)
                              end if ! Br or Bz ?
                           end do ! phi Loop
                        end do    ! theta loop
                     end if
                  end do          ! r loop

               end if ! l_loop ?


            end if ! n_surface ?

         end do  ! LOOP over fields in movie file

      end do ! Loop over movies

   end subroutine calc_dtB_frame_IC
!----------------------------------------------------------------------------
   subroutine get_dtB(dtB,dtBLM,DimB1,DimB2,DimB3,DimB4,n_r,l_ic)

      !-- Input variables:
      integer,     intent(in) :: n_r             ! No. of radial grid point
      logical,     intent(in) :: l_ic            ! =true if inner core field
      integer,     intent(in) :: DimB1,DimB2,DimB3,DimB4
      complex(cp), intent(in) :: dtBLM(DimB1:DimB2,DimB3:DimB4)

      !-- Output variables:
      real(cp),intent(out) ::  dtB(:)     ! Result Field with dim>=n_theta_block

      !-- Local variables:
      integer :: n_theta, l, lm
      real(cp) :: r_ratio, O_r
      real(cp) :: r_dep(l_max)  ! (r/r_ICB)**l / r_ICB
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)
      complex(cp) :: Tl_AX(l_max+1), Tl_AX_loc(l_max+1)

      !-- Calculate radial dependencies:
      !     for IC: (r/r_ICB)**l / r_ICB
      !     for OC: 1/r
      if ( l_ic ) then
         r_ratio =r_ic(n_r)/r_icb
         r_dep(1)=r_ratio/r_icb
         do l=2,l_max
            r_dep(l)=r_dep(l-1)*r_ratio
         end do
      else
         O_r=or1(n_r)
      end if

      Tl_AX_loc(:)=zero
      do l=1,l_max
         if ( l_ic ) then
            lm=lo_map%lm2(l,0)
            if ( lm >= DimB1 .and. lm <= DimB2 ) Tl_AX_loc(l+1)=r_dep(l)*dtBLM(lm,n_r)
         else
            lm=lm2(l,0)
            Tl_AX_loc(l+1)=O_r*dtBLM(lm,n_r)
         end if
      end do

      if ( l_ic ) then
#ifdef WITH_MPI
         call MPI_Reduce(Tl_AX_loc, Tl_AX, l_max+1, MPI_DEF_COMPLEX, MPI_SUM, 0, &
              &          MPI_COMM_WORLD, ierr)
#else
         Tl_AX(:)=Tl_AX_loc(:)
#endif
      else
         Tl_AX(:)=Tl_AX_loc(:)
      end if

      if ( l_ic ) then
         call toraxi_to_spat(Tl_AX, tmpt, tmpp, l_max)
      else
         call toraxi_to_spat(Tl_AX, tmpt, tmpp, l_R(n_r))
      end if

      do n_theta=1,n_theta_max
         dtB(n_theta)=O_sin_theta(n_theta)*tmpp(n_theta)
      end do

   end subroutine get_dtB
!-------------------------------------------------------------------------
   subroutine get_Bpol(PolLM,dPolLM,Br,Bt,Bp,rT,lIC)
      !
      !  Purpose of this subroutine is to calculate the components
      !  Br, Bt, and Bp of the poloidal magnetic field PolLM (l,m space)
      !  at the radial grid point r=rT and the block of theta grid points
      !  from n_theta=n_theta_start to n_theta=n_theta_start+n_theta_block-1
      !  and for all phis.
      !  For lIC=.true. the inner core field is calculated,
      !  to get the IC field for a conducting inner core PolLM has to be
      !  the poloidal field at the ICB.
      !

      !-- Input variables
      real(cp),    intent(in) :: rT             ! radius
      logical,     intent(in) :: lIC            ! true for inner core, special rDep !
      complex(cp), intent(in) :: PolLM(lm_max)  ! field in (l,m)-space for rT
      complex(cp), intent(in) :: dPolLM(lm_max) ! dr field in (l,m)-space for rT

      !-- Output variables
      real(cp), intent(out) :: Br(:,:),Bt(:,:),Bp(:,:)

      !-- Local variables
      integer :: lm,m,l
      real(cp) :: rRatio,rDep(0:l_max)
      real(cp) :: O_r_E_2
      complex(cp) :: cs1(lm_max),cs2(lm_max)
      complex(cp) :: zeros(lm_max)

      !-- Calculate radial dependencies: (r_ic(1)=r(n_r_max)=inner core boundary)
      !   Radial dependence = (r/r_ICB)**(l-1) / r_ICB**2
      if ( lIC ) then
         rRatio=rT/r_icb
         rDep(0)=one/(rT*r_icb)
         do l=1,l_max
            rDep(l)=rDep(l-1)*rRatio
         end do
      !--- NOTE: field for insulating inner core has same rDep but uses poloidal
      !          field at r=r_icb
      else
         O_r_E_2=one/(rT*rT)
      end if
      !--- NOTE: mantle potential field may be included by adding special rDep
      !          and using poloidal field at r_cmb

      !-- Get coeffs with radial dependence:
      if ( lIC ) then
         if ( l_cond_ic ) then
            do m=0,m_max,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  cs1(lm)=rDep(l)*PolLM(lm)
                  cs2(lm)=rDep(l)*( real(l+1,cp)*PolLM(lm) + rT*dPolLM(lm) )
               end do
            end do
         else
            do m=0,m_max,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  cs1(lm)=rDep(l)*PolLM(lm)
                  cs2(lm)=cs1(lm)/real(l,cp)
               end do
            end do
         end if
      else
         do m=0,m_max,minc
            do l=m,l_max
               lm=lm2(l,m)
               cs1(lm)=O_r_E_2*PolLM(lm)
               cs2(lm)=O_r_E_2*rT*dPolLM(lm)
            end do
         end do
      endif

      zeros(:)=zero

      call torpol_to_spat_single(cs1, cs2, zeros, Br, Bt, Bp, l_max)

   end subroutine get_Bpol
!-------------------------------------------------------------------------------------
   subroutine get_Btor(Tlm,Bt,Bp,rT,lIC)
      !
      !  Purpose of this subroutine is to calculate the components
      !  Bt and Bp of the toroidal magnetic field  Tlm (in l,m space)
      !  at the radial grid point r=rT and the block of theta grid points
      !  from n_theta=n_theta_start to n_theta=n_theta_start+n_theta_block-1
      !  and for all phis.
      !  For lIC=.true. the inner core field is calculated,
      !  to get the IC field for a conducting inner core Plm has to be
      !  the toroidal field at the ICB.
      !

      !-- Input variables:
      real(cp),    intent(in) :: rT            ! radius
      logical,     intent(in) :: lIC           ! true for inner core, special rDep !
      complex(cp), intent(in) :: Tlm(lm_max)   ! field in (l,m)-space for rT

      !-- Output variables:
      real(cp), intent(out) :: Bt(:,:),Bp(:,:)

      !-- Local variables:
      integer :: lm,m,l
      real(cp) :: rRatio,rDep(0:l_max)
      complex(cp) :: cs1(lm_max)
      complex(cp) :: zeros(lm_max)

      if ( lIC .and. .not. l_cond_ic )    return

      !-- Calculate radial dependencies: (r_ic(1)=r(n_r_max)=inner core boundary)
      !   Radial dependence = (r/r_ICB)**(l-1) / r_ICB**2
      if ( lIC ) then
         rRatio=rT/r_icb
         rDep(0)=one/r_icb
         do l=1,l_max
            rDep(l)=rDep(l-1)*rRatio
         end do
      else
         do l=0,l_max
            rDep(l)=one/rT
         end do
      end if
      !--- NOTE: mantle potential field may be included by adding special rDep
      !          and using poloidal field at r_cmb


      !-- Get coeffs with radial dependence:
      do m=0,m_max,minc
         do l=m,l_max
            lm=lm2(l,m)
            cs1(lm)=rDep(l)*Tlm(lm)
            zeros(lm)=zero
         end do
      end do

      call sphtor_to_spat(sht_l_single, zeros, cs1, Bt, Bp, l_max)

   end subroutine get_Btor
!-------------------------------------------------------------------------------------
   subroutine lm2pt(alm,aij,rT,lIC,lrComp)
      !
      ! Spherical harmonic transform from alm(l,m) to aij(phi,theta)
      ! Radial field components are calculated for lrComp=.true.
      ! Done within the range [n_theta_min,n_theta_min+n_theta_block-1]
      ! Used only for graphic output.
      !

      !-- Input variables:
      real(cp),    intent(in) ::  rT
      logical,     intent(in) :: lIC         ! true for inner core, extra factor !
      logical,     intent(in) :: lrComp      ! true for radial field components
      complex(cp), intent(in) :: alm(*)      ! field in (l,m)-space

      !-- Output variables:
      real(cp), intent(out) :: aij(:,:)  ! field in (theta,phi)-space

      !-- Local variables:
      integer :: l,lm     ! degree/order
      complex(cp) :: cs1(lm_max) ! help array
      real(cp) :: O_r_E_2,rRatio

      O_r_E_2=one/(rT*rT)
      rRatio=rT/r_icb

      !-- Multiplication with l(l+1)/r**2 for radial component:
      cs1(1)=zero
      do lm=2,lm_max
         l = lm2l(lm)
         cs1(lm)=alm(lm)
         if ( lrComp ) cs1(lm)=cs1(lm)*O_r_E_2*dLh(lm)
         if ( lIC )    cs1(lm)=rRatio**real(l+1,cp)*cs1(lm)
      end do

      call scal_to_spat(sht_l_single, cs1, aij, l_max)

   end subroutine lm2pt
!---------------------------------------------------------------------------------
end module out_dtB_frame
