module out_dtB_frame

   use truncation, only: n_r_max, lm_max, n_r_ic_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag, l_max, m_max, minc, n_theta_max,     &
       &                 n_r_tot, nlat_padded, n_phi_max, n_cheb_ic_max
   use precision_mod
   use radial_functions, only: r, or1, chebt_ic, r_ic, rscheme_oc, r_icb, &
       &                       dr_fac_ic, chebt_ic_even
   use blocking, only: lm2m, lm2l, lm2
   use horizontal_data, only: cosTheta, n_theta_cal2ord, sinTheta, &
       &                      dLh, O_sin_theta
   use dtB_mod, only: PstrLM, PadvLM, PdifLM, TstrLM, TadvLM, TdifLM, &
       &              PadvLMIC, PdifLMIC, TadvLMIC, TomeLM, TdifLMIC
   use movie_data, only: n_movie_type, n_movie_fields, n_movie_fields_ic, &
       &                 n_movie_file, n_movie_const, n_movie_surface,    &
       &                 movie_const, n_movie_field_type
   use logic, only: l_cond_ic
   use constants, only: zero, one, ci
   use radial_der_even, only: get_drNS_even
   use radial_der, only: get_dr
   use sht, only: torpol_to_spat_single, sphtor_to_spat, scal_to_spat, &
       &          toraxi_to_spat, sht_l_single

   implicit none

   private

   public :: write_dtB_frame

contains

   subroutine write_dtB_frame(n_movie,b,db,aj,dj,b_ic,db_ic,aj_ic,dj_ic)
      !
      !  Controls output of specific movie frames related to magnetic field
      !  production and diffusion.
      !
    
      !-- Input of variables:
      integer,     intent(in) :: n_movie
      complex(cp), intent(in) :: b(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: db(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: dj(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(lm_maxMag,n_r_ic_maxMag)
    
      !-- Local variables:
      integer :: n_type
      integer :: n_out
      integer :: n_fields_oc
      integer :: n_fields_ic
      integer :: n_fields,n_field
      integer :: n_const
    
      integer :: n_r,n_rC,n_r_loop_max,n
      integer :: n_theta,n_theta_cal,n_phi,n_pos,n_or
      integer :: lm,l,m
    
      integer :: n_surface
      integer :: n_field_type
      integer :: n_field_size
    
      real(cp) :: dtB(n_theta_max)
      real(cp) :: dtBframe(n_r_tot*n_theta_max)
      real(cp) :: dtBr(nlat_padded,n_phi_max),dtBt(nlat_padded,n_phi_max)
      real(cp) :: dtBp(nlat_padded,n_phi_max)
      real(cp) :: dtBrframe(n_r_tot*n_phi_max*n_theta_max)
      real(cp) :: const,rMov
    
      complex(cp) :: workA(lm_max,n_r_max)
      complex(cp) :: workB(lm_max,n_r_ic_max)
    
      logical :: l_loop
    
      n_type      =n_movie_type(n_movie)
      n_fields_oc =n_movie_fields(n_movie)
      n_fields_ic =n_movie_fields_ic(n_movie)
      n_fields    =max(n_fields_ic,n_fields_oc)
      n_out       =n_movie_file(n_movie)
      n_const     =n_movie_const(n_movie)
      n_surface   =n_movie_surface(n_movie)
      const       =movie_const(n_movie)
    
      !--- Axisymmetric dtFL or dtAB:
    
      if ( n_type == 31 .or. n_type == 32 .or. n_type == 33 .or. &
      &    n_type == 41 .or. n_type == 42 .or. n_type == 43 ) then
    
         n_field_size=n_theta_max*(n_r_max+n_r_ic_max-2)
    
         do n_field=1,n_fields
    
            n_field_type=n_movie_field_type(n_field,n_movie)
    
            !--- OC contribution:
            do n_r=1,n_r_max
    
               if ( n_field_type == 20 ) then
                  call get_dtB(dtB,PstrLM,lm_max,n_r_max,n_r,.false.)
               else if ( n_field_type == 21 ) then
                  call get_dtB(dtB,PadvLM,lm_max,n_r_max,n_r,.false.)
               else if ( n_field_type == 22 ) then
                  call get_dtB(dtB,PdifLM,lm_max,n_r_max,n_r,.false.)
               else if ( n_field_type == 23 ) then
                  call get_dtB(dtB,TstrLM,lm_max,n_r_max,n_r,.false.)
               else if ( n_field_type == 25 ) then
                  call get_dtB(dtB,TadvLM,lm_max,n_r_max,n_r,.false.)
               else if ( n_field_type == 26 ) then
                  call get_dtB(dtB,TdifLM,lm_max,n_r_max,n_r,.false.)
               end if
    
               !--- Store in frame field:
               do n_theta_cal=1,n_theta_max
                  n_theta=n_theta_cal2ord(n_theta_cal)
                  n_pos  =(n_r-1)*n_theta_max+n_theta
                  dtBframe(n_pos)=dtB(n_theta)
               end do
    
            end do
    
            !--- Now IC contribution:
            do n_r=2,n_r_ic_max-1
    
               if ( n_field_type == 21 ) then
                  call get_dtB(dtB,PadvLMIC,lm_max,n_r_ic_max,n_r,.true.)
               else if ( n_field_type == 22 ) then
                  call get_dtB(dtB,PdifLMIC,lm_max,n_r_ic_max,n_r,.true.)
               else if ( n_field_type == 25 ) then
                  call get_dtB(dtB,TadvLMIC,lm_max,n_r_ic_max,n_r,.true.)
               else if ( n_field_type == 26 ) then
                  call get_dtB(dtB,TdifLMIC,lm_max,n_r_ic_max,n_r,.true.)
               else
                  do n_theta=1,n_theta_max
                     dtB(n_theta)=0.0_cp
                  end do
               end if
    
               !--- Store in frame field:
               do n_theta_cal=1,n_theta_max
                  n_theta=n_theta_cal2ord(n_theta_cal)
                  n_pos  =(n_r_max+n_r-2)*n_theta_max+n_theta
                  dtBframe(n_pos)=dtB(n_theta)
               end do
    
            end do
    
            !--- Write frame field:
            write(n_out) (real(dtBframe(n),kind=outp),n=1,n_field_size)
    
         end do
    
         !--- dtBr:
      else ! non-axisymmetric fields
    
         do n_field=1,n_fields
    
            n_field_type=n_movie_field_type(n_field,n_movie)
    
            if ( n_field_type == 24 ) then
               ! This reduces omega effect to field production of axisymm. toroidal field:
               do n_r=1,n_r_max
                  do lm=1,lm_max
                     m=lm2m(lm)
                     if ( m == 0 ) then
                        workA(lm,n_r)=TomeLM(lm,n_r)
                     else
                        workA(lm,n_r)=zero
                     end if
                  end do
               end do
            else if ( n_field_type == 81 ) then
               ! This reduces poloidal field to the dipole contribution:
               do n_r=1,n_r_max
                  do lm=1,lm_max
                     l=lm2l(lm)
                     if ( l == 1 ) then
                        workA(lm,n_r)=b(lm,n_r)
                     else
                        workA(lm,n_r)=zero
                     end if
                  end do
               end do
            end if
    
            !------ Outer core contribution:
    
            !------ Calculate needed radial derivatives:
            if ( n_field_type == 35 ) then
               call get_dr(PstrLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            else if ( n_field_type == 36 ) then
               call get_dr(PadvLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            else if ( n_field_type == 37 ) then
               call get_dr(PdifLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            else if ( n_field_type == 38 ) then
               call get_dr(TstrLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            else if ( n_field_type == 39 ) then
               call get_dr(TomeLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            else if ( n_field_type == 40 ) then
               call get_dr(TadvLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            else if ( n_field_type == 41 ) then
               call get_dr(TdifLM,workA,lm_max,1,lm_max, &
                    &      n_r_max,rscheme_oc,nocopy=.true.)
            end if
    
            if ( n_surface == 0 ) then
               n_r_loop_max=n_r_max
               n_field_size=(n_r_max+n_r_ic_max-2)*n_theta_max*n_phi_max
            else if ( n_surface == 1 ) then
               n_r_loop_max=1
               n_field_size=n_theta_max*n_phi_max
            end if
    
            do n_rC=1,n_r_loop_max
               n_r=n_rC
               if ( n_surface == 0 ) then
                  rMov=r(n_r)
                  n_or=(n_r-1)*n_theta_max*n_phi_max
               else if ( n_surface == 1 ) then
                  n_r=n_const
                  rMov=const
                  n_or=0
               end if
    
               !-- Br:
               if ( n_field_type == 27 ) then
                  call get_Bpol(PstrLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,rMov, &
                       &        .false.)
               else if ( n_field_type == 28 ) then
                  call get_Bpol(PadvLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 29 ) then
                  call get_Bpol(PdifLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 81 ) then
                  !---------- get radial field diffusion and radial dipole field:
                  call get_Bpol(PdifLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
                  call get_Bpol(workA(:,n_r),workA(:,n_r),dtBt,dtBp,dtBp,  &
                       &        rMov,.false.)
                  !-- Jr:
               else if ( n_field_type == 30 ) then
                  call get_Bpol(aj(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,  &
                       &        rMov,.false.)
               else if ( n_field_type == 31 ) then
                  call get_Bpol(TstrLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 32 ) then
                  call get_Bpol(workA(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 33 ) then
                  call get_Bpol(TadvLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 34 ) then
                  call get_Bpol(TdifLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
                  !-- Bz poloidal
               else if ( n_field_type == 13 ) then
                  call get_Bpol(b(:,n_r),db(:,n_r),dtBr,dtBt,dtBp,rMov, &
                       &        .false.)
               else if ( n_field_type == 35 ) then
                  call get_Bpol(PstrLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 36 ) then
                  call get_Bpol(PadvLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 37 ) then
                  call get_Bpol(PdifLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
 
                  !-- Jz poloidal
               else if ( n_field_type == 14 ) then
                  call get_Bpol(aj(:,n_r),dj(:,n_r),dtBr,dtBt,dtBp,rMov, &
                       &        .false.)
               else if ( n_field_type == 38 ) then
                  call get_Bpol(TstrLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 39 ) then
                  call get_Bpol(workA(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 40 ) then
                  call get_Bpol(TadvLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
               else if ( n_field_type == 41 ) then
                  call get_Bpol(TdifLM(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                       &        rMov,.false.)
                  !--- Bp toroidal:
               else if ( n_field_type == 24 ) then
                  call get_Btor(workA(:,n_r),dtBt,dtBr,rMov,.false.)
               else if ( n_field_type == 42 ) then
                  call get_Btor(aj(:,n_r),dtBt,dtBr,rMov,.false.)
               else if ( n_field_type == 43 ) then
                  call get_Btor(TstrLM(:,n_r),dtBt,dtBr,rMov,.false.)
               else if ( n_field_type == 44 ) then
                  call get_Btor(workA(:,n_r),dtBt,dtBr,rMov,.false.)
               else if ( n_field_type == 45 ) then
                  call get_Btor(TadvLM(:,n_r),dtBt,dtBr,rMov,.false.)
               else if ( n_field_type == 46 ) then
                  call get_Btor(TdifLM(:,n_r),dtBt,dtBr,rMov,.false.)
               else if ( n_field_type == 49 ) then
                  call get_Btor(TomeLM(:,n_r),dtBt,dtBr,rMov,.false.)
                  !--- Bt toroidal
               else if ( n_field_type == 50 ) then
                  call get_Btor(aj(:,n_r),dtBr,dtBt,rMov,.false.)
                  !--- Toroidal Potential
               else if ( n_field_type == 51 ) then
                  call lm2pt(aj(:,n_r),dtBr,rMov,.false.,.false.)
                  !--- Fieldlines for theta=const.
               else if ( n_field_type == 52 ) then
                  call get_Btor(b(:,n_r),dtBr,dtBt,rMov,.false.)
               end if
 
               !--- Now store the stuff of theta block to dtBrframe:
               do n_phi=1,n_phi_max
                  do n_theta_cal=1,n_theta_max
                     n_theta=n_theta_cal2ord(n_theta_cal)
                     n_pos=n_or+n_phi+(n_theta-1)*n_phi_max
                     if ( n_field_type == 13 .or. n_field_type == 35 .or. &
                     &    n_field_type == 36 .or. n_field_type == 37 .or. &
                     &    n_field_type == 14 .or. n_field_type == 38 .or. &
                     &    n_field_type == 39 .or. n_field_type == 40 .or. &
                     &    n_field_type == 41 ) then
                        dtBrframe(n_pos) =                                   &
                        &    cosTheta(n_theta_cal)*dtBr(n_theta_cal,n_phi) - &
                        &    sinTheta(n_theta_cal)*dtBt(n_theta_cal,n_phi)
                     else if ( n_field_type == 81 ) then
                        if ( dtBr(n_theta_cal,n_phi) * &
                             dtBt(n_theta_cal,n_phi) < 0.0_cp ) then
                           dtBrframe(n_pos)=dtBr(n_theta_cal,n_phi)
                        else
                           dtBrframe(n_pos)=0.0_cp
                        end if
                     else
                        dtBrframe(n_pos)=dtBr(n_theta_cal,n_phi)
                     end if ! Br (Bp) or Bz ?
                  end do ! phi Loop
               end do    ! theta loop
    
            end do          ! r loop
    
            !------ Inner core contribution:
            l_loop=.true.
            if ( n_surface == 0 ) then
               n_r_loop_max=n_r_ic_max-1
            else if ( n_surface == 1 ) then
               n_r_loop_max=2
               if ( const >= r_icb ) l_loop= .false. 
            end if
            if ( l_loop ) then
    
               !------ Calculate needed radial derivatives:
               if ( n_field_type == 36 ) then
                  call get_drNS_even(PadvLMIC,workA,lm_max,1,lm_max,n_r_ic_max, &
                       &             n_cheb_ic_max,dr_fac_ic,workB,chebt_ic,    &
                       &             chebt_ic_even)
               else if ( n_field_type == 37 ) then
                  call get_drNS_even(PdifLMIC,workA,lm_max,1,lm_max,n_r_ic_max, &
                       &             n_cheb_ic_max,dr_fac_ic,workB,chebt_ic,    &
                       &             chebt_ic_even)
               else if ( n_field_type == 40 ) then
                  call get_drNS_even(TadvLMIC,workA,lm_max,1,lm_max,n_r_ic_max, &
                       &             n_cheb_ic_max,dr_fac_ic,workB,chebt_ic,    &
                       &             chebt_ic_even)
               else if ( n_field_type == 41 ) then
                  call get_drNS_even(TdifLMIC,workA,lm_max,1,lm_max,n_r_ic_max, &
                       &             n_cheb_ic_max,dr_fac_ic,workB,chebt_ic,    &
                       &             chebt_ic_even)
               end if
    
               do n_rC=2,n_r_loop_max
                  n_r=n_rC
                  if ( n_surface == 0 ) then
                     rMov=r_ic(n_r)
                     n_or=(n_r_max+n_r-2) * n_theta_max*n_phi_max
                  else if ( n_surface == 1 ) then
                     n_r=n_const
                     rMov=const
                     n_or=0
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
                     call get_Bpol(PadvLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,&
                          &        rMov,.true.)
                  else if ( n_field_type == 29 ) then
                     call get_Bpol(PdifLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,&
                          &        rMov,.true.)
                     !--- Jr:
                  else if ( n_field_type == 30 ) then
                     call get_Bpol(aj_ic(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.true.)
                  else if ( n_field_type == 33 ) then
                     call get_Bpol(TadvLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,&
                          &        rMov,.true.)
                  else if ( n_field_type == 34 ) then
                     call get_Bpol(TdifLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,&
                          &        rMov,.true.)
                     !- Bz poloidal:
                  else if ( n_field_type == 13 ) then
                     call get_Bpol(b_ic(:,n_r),db_ic(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.true.)
                  else if ( n_field_type == 36 ) then
                     call get_Bpol(PadvLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.true.)
                  else if ( n_field_type == 37 ) then
                     call get_Bpol(PdifLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.true.)
                     !--- Jz poloidal:
                  else if ( n_field_type == 14 ) then
                     call get_Bpol(aj_ic(:,n_r),dj_ic(:,n_r),dtBr,dtBt,dtBp, &
                          &        rMov,.true.)
                  else if ( n_field_type == 40 ) then
                     call get_Bpol(TadvLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,&
                          &        rMov,.true.)
                  else if ( n_field_type == 41 ) then
                     call get_Bpol(TdifLMIC(:,n_r),workA(:,n_r),dtBr,dtBt,dtBp,&
                          &        rMov,.true.)
                     !--- Bphi toroidal:
                  else if ( n_field_type == 42 ) then
                     call get_Btor(aj_ic(:,n_r),dtBt,dtBr,rMov,.true.)
                  else if ( n_field_type == 45 ) then
                     call get_Btor(TadvLMIC(:,n_r),dtBt,dtBr,rMov,.true.)
                  else if ( n_field_type == 46 ) then
                     call get_Btor(TdifLMIC(:,n_r),dtBt,dtBr,rMov,.true.)
                     !--- Btheta toroidal:
                  else if ( n_field_type == 50 ) then
                     call get_Btor(aj_ic(:,n_r),dtBr,dtBt,rMov,.true.)
                     !--- Toroidal Potential
                  else if ( n_field_type == 51 ) then
                     call lm2pt(aj_ic(:,n_r),dtBr,r_ic(n_r),.true.,.false.)
                     !--- Fieldlines for theta=const.
                  else if ( n_field_type == 52 ) then
                     call get_Btor(b_ic(:,n_r),dtBr,dtBt,rMov,.true.)
                  end if
 
                  !--- Now store the stuff:
                  do n_phi=1,n_phi_max
                     do n_theta_cal=1,n_theta_max
                        n_theta=n_theta_cal2ord(n_theta_cal)
                        n_pos=n_or+n_phi+(n_theta-1)*n_phi_max
                        if ( n_field_type == 13 .or. n_field_type == 36 .or. &
                        &    n_field_type == 37 .or. n_field_type == 14 .or. &
                        &    n_field_type == 40 .or. n_field_type == 41 ) then
                           dtBrframe(n_pos) =                                  &
                           &    cosTheta(n_theta_cal)*dtBr(n_theta_cal,n_phi) -&
                           &    sinTheta(n_theta_cal)*dtBt(n_theta_cal,n_phi)
                        else
                           dtBrframe(n_pos)=dtBr(n_theta_cal,n_phi)
                        end if ! Br or Bz ?
                     end do ! phi Loop
                  end do    ! theta loop
    
               end do          ! r loop
    
            end if ! l_loop ?
    
            !--- Write frame field:
            write(n_out) (real(dtBrframe(n),kind=outp),n=1,n_field_size)
    
         end do  ! LOOP over fields in movie file
    
      end if

   end subroutine write_dtB_frame
!----------------------------------------------------------------------------
   subroutine get_dtB(dtB,dtBLM,DimB1,DimB2,n_r,l_ic)

      !-- Input variables:
      integer,     intent(in) :: n_r             ! No. of radial grid point
      logical,     intent(in) :: l_ic            ! =true if inner core field
      integer,     intent(in) :: DimB1,DimB2
      complex(cp), intent(in) :: dtBLM(DimB1,DimB2)
    
      !-- Output variables:
      real(cp),intent(out) ::  dtB(:)     ! Result Field with dim>=n_theta_block
    
      !-- Local variables:
      integer :: n_theta         ! No. of theta
      integer :: l,lm            ! Degree, counter for degree/order combinations
      real(cp) :: r_ratio          ! r/r_ICB
      real(cp) :: O_r              ! 1/r
      real(cp) :: r_dep(l_max)     ! (r/r_ICB)**l / r_ICB
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)
      complex(cp) :: Tl_AX(1:l_max+1)
    
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
    
      Tl_AX(1)=zero
      do l=1,l_max
         lm=lm2(l,0)
         if ( l_ic ) then
            Tl_AX(l+1)=r_dep(l)*dtBLM(lm,n_r)
         else
            Tl_AX(l+1)=O_r*dtBLM(lm,n_r)
         end if
      end do

      call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:))

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

      do lm=1,lm_max
         zeros(lm)=zero
      end do

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
