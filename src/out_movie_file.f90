module out_movie
   !
   ! This module handles the storage of the relevant diagnostics in the public
   ! array frames(*) and then the writing of the corresponding movie files
   !

   use precision_mod
   use parallel_mod, only: rank
   use geos, only: cyl
   use truncation, only: n_phi_max, n_theta_max, minc, lm_max, l_max,    &
       &                 n_m_max, lm_maxMag, n_r_maxMag, n_r_ic_maxMag,  &
       &                 n_r_ic_max, n_r_max, nlat_padded
   use movie_data, only: frames, n_movie_fields, n_movies, n_movie_surface, &
       &                 n_movie_const, n_movie_field_type, lGeosField,     &
       &                 n_movie_field_start,n_movie_field_stop,            &
       &                 movie_const, movie_file, n_movie_file, lICField,   &
       &                 n_movie_fields_ic
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: orho1, orho2, or1, or2, or3, or4, beta,  &
       &                       r_surface, r_cmb, r, r_ic, temp0, l_R
   use physical_parameters, only: LFfac, radratio, ra, ek, pr, prmag, raxi, sc
   use num_param, only: vScale, tScale
   use blocking, only: lm2l, lm2, llmMag, ulmMag
   use horizontal_data, only: O_sin_theta, sinTheta, cosTheta,    &
       &                      n_theta_cal2ord, O_sin_theta_E2,    &
       &                      phi, theta_ord
   use fields, only: w_Rloc, b_Rloc, b_ic, bICB
   use sht, only: torpol_to_spat, toraxi_to_spat
   use logic, only: l_save_out, l_cond_ic, l_mag, l_full_sphere
   use constants, only: zero, half, one, two
   use output_data, only: runid, n_s_max
   use useful, only: abortRun

   implicit none

   private

   public :: store_movie_frame, write_movie_frame, get_fl

contains

   subroutine store_movie_frame(n_r,vr,vt,vp,br,bt,bp,sr,drSr,xir,phir,  &
              &                 dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt)
      !
      !  Controls output of movie frames.
      !  Usually called from radialLoop.
      !

      !-- Input variables:
      integer,     intent(in) :: n_r                ! radial grid point no.
      real(cp),    intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp),    intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp),    intent(in) :: sr(:,:),drSr(:,:)
      real(cp),    intent(in) :: xir(:,:),phir(:,:)
      real(cp),    intent(in) :: dvrdp(:,:),dvpdr(:,:)
      real(cp),    intent(in) :: dvtdr(:,:),dvrdt(:,:)
      real(cp),    intent(in) :: cvr(:,:)
      real(cp),    intent(in) :: cbr(:,:),cbt(:,:)

      !-- Local variables:
      integer :: n_movie        ! No. of movie
      integer :: n_field        ! No. of field
      integer :: n_surface      ! Surface (1=r,2=theta,3=phi)
      integer :: n_const        ! Gives surface
      integer :: n_field_type   ! Numbers field types
      integer :: n_store_last   ! Position i in frame(i) were field starts-1
      integer :: n_theta
      integer :: n_field_size
      integer :: n_fields
      logical :: lThetaFound
      complex(cp) :: bCMB(lm_max)

      if ( n_r==n_r_cmb .and. l_mag ) bCMB = b_Rloc(:,n_r_cmb)

      do n_movie=1,n_movies

         n_fields =n_movie_fields(n_movie)
         n_surface=n_movie_surface(n_movie)
         n_const  =n_movie_const(n_movie)

         select case ( n_surface )

            case(-1) ! Earth surface
               if ( n_r /= 1 ) cycle  ! not CMB radius
               do n_field=1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_sur(n_store_last,n_field_type,bCMB)
                  end if
               end do

            case(0) ! 3d
               do n_field=1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_3d(vr,vt,vp,br,bt,bp,sr,drSr,xir,phir,  &
                          &               dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
                          &               n_r,n_store_last,n_field_type)
                  end if
               end do

            case(1) ! Surface r=constant
               if ( n_r /= n_const ) cycle  ! not desired radius

               do n_field=1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_r(vr,vt,vp,br,bt,bp,sr,drSr,xir,phir, &
                          &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,n_r,    &
                          &              n_store_last,n_field_type)
                  end if
               end do

            case(2) ! Surface theta=constant
               !------ Test whether n_theta_movie is in the current theta block
               !       and find its position n_theta_movie_c:
               lThetaFound=.false.
               do n_theta=1,n_theta_max
                  if ( n_theta == n_const ) then
                     lThetaFound=.true.
                     exit
                  end if
               end do
               if ( .not. lThetaFound) cycle        ! Theta not found !

               do n_field=1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_t(vr,vt,vp,br,bt,bp,sr,drSr,xir,phir,   &
                          &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbt,n_r,  &
                          &              n_store_last,n_field_type,n_theta)
                    end if
               end do

            case(3)  ! Surface phi=const.
               do n_field=1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  n_field_size=(n_movie_field_stop(n_field,n_movie) - n_store_last)/2
                  if ( n_store_last >= 0 ) then
                     call store_fields_p(vr,vt,vp,br,bp,bt,sr,drSr,xir,phir,   &
                          &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt,  &
                          &              n_r,n_store_last,n_field_type,        &
                          &              n_const,n_field_size)
                  end if
               end do  ! Do loop over field for one movie

            case(4)  ! phi averages
               do n_field=1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_p_avg(vr,vt,vp,bp,sr,drSr,xir,phir,dvrdp,       &
                          &                  dvpdr,dvtdr, dvrdt,cvr,n_r,n_store_last,  &
                          &                  n_field_type)
                  end if
               end do  ! Do loop over field for one movie

         end select

      end do  ! Do loop over movies !

   end subroutine store_movie_frame
!----------------------------------------------------------------------------
   subroutine write_movie_frame(n_frame, time)
      !
      !  Writes different movie frames into respective output files.
      !  Called from rank 0 with full arrays in standard LM order.
      !

      !-- Input of variables:
      integer,  intent(in) :: n_frame
      real(cp), intent(in) :: time

      !-- Local variables:
      integer :: n_fields, n_fields_ic, n_fields_oc
      integer :: n_movie, n_surface, n_out, n_field, n, n_start, n_stop
      integer :: n_r_mov_tot, version
      real(cp) :: const
      real(cp) :: r_mov_tot(n_r_max+n_r_ic_max)

      do n_movie=1,n_movies

         n_surface  =n_movie_surface(n_movie)
         n_fields_oc=n_movie_fields(n_movie)
         n_fields_ic=n_movie_fields_ic(n_movie)
         n_fields   =max(n_fields_ic,n_fields_oc)
         n_out      =n_movie_file(n_movie)
         const      =movie_const(n_movie)
         if ( n_surface == 1 ) const=const/r_cmb

         !------ Open movie file:
         if ( l_save_out .and. rank == 0 ) then
            open(newunit=n_out, file=movie_file(n_movie), status='unknown', &
            &    form='unformatted', position='append')
         end if

         !------ Write header if this is the first frame:

         if ( n_frame == 1 .and. rank == 0 ) then

            !------ Start with info about movie type:
            version=3
            write(n_out) version
            write(n_out) n_surface, n_fields
            write(n_out) const
            write(n_out) (n_movie_field_type(n,n_movie),n=1,n_fields)

            !------ Combine OC and IC radial grid points:
            n_r_mov_tot=n_r_max
            r_mov_tot(:n_r_max)=r(:)
            if ( n_r_ic_max > 0 .and. lICField(n_movie) ) then
               n_r_mov_tot=n_r_mov_tot+n_r_ic_max
               r_mov_tot(n_r_max+1:)=r_ic(:)
            end if

            !------ Now other info about grid and parameters:
            write(n_out) runid          ! run identifyer (as set in namelist contrl)
            write(n_out) n_r_mov_tot, n_r_max, n_r_ic_max, n_theta_max, n_phi_max, &
            &            n_s_max, minc
            write(n_out) real(ra,outp), real(ek,outp), real(pr,outp), real(raxi,outp), &
            &            real(sc,outp), real(prmag,outp), real(radratio,outp),         &
            &            real(tScale,outp)

            !------ Write grid:
            if ( lGeosField(n_movie) ) &
            &   write(n_out) real(cyl(:)/r_cmb,kind=outp)
            write(n_out) real(r_mov_tot(:)/r_cmb,kind=outp)
            write(n_out) real(theta_ord(:),kind=outp)
            write(n_out) real(phi(:),kind=outp)

         end if  ! Write header ?

         if ( rank == 0 ) then
            !------ Write frame number, time and IC and MA rotation rates::
            write(n_out) real(time,outp)
         end if

         !------ Write frames:
         if ( rank == 0 ) then
            do n_field=1,n_fields
               n_start=n_movie_field_start(n_field,n_movie)
               if ( n_fields_oc > 0 ) then
                  n_stop =n_movie_field_stop(n_field,n_movie)
               end if
               if ( n_fields_ic > 0 ) then
                  n_stop=n_movie_field_stop(n_fields_oc + n_field,n_movie)
               end if
               write(n_out) (real(frames(n),kind=outp),n=n_start,n_stop)
            end do
         end if

         if ( l_save_out .and. rank==0 ) close(n_out)

      end do  ! Loop over movies

   end subroutine write_movie_frame
!----------------------------------------------------------------------------
   subroutine store_fields_sur(n_store_last,n_field_type,bCMB)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !--- Input variables:
      integer,     intent(in) :: n_store_last     ! Start position for storing -1
      integer,     intent(in) :: n_field_type     ! Defines field type
      complex(cp), intent(in) :: bCMB(lm_max)

      !--- Local variables:
      integer :: n_theta, n_theta_cal, n_phi, n_o

      !----- Magnetic field at surface (theta-blocks):
      real(cp) :: br_sur(nlat_padded,n_phi_max) ! Radial magnetic field in (phi,theta)-space
      real(cp) :: bt_sur(nlat_padded,n_phi_max) ! Latitudinal magnetic field
      real(cp) :: bp_sur(nlat_padded,n_phi_max) ! Azimuthal magnetic field.

      call get_B_surface(br_sur,bt_sur,bp_sur,bCMB)

      if ( n_field_type == 1 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=br_sur(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 2 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=bt_sur(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 3 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=bp_sur(n_theta_cal,n_phi)
            end do
         end do

      end if

   end subroutine store_fields_sur
!----------------------------------------------------------------------------
   subroutine store_fields_r(vr,vt,vp,br,bt,bp,sr,drSr,xir,phir,dvrdp,dvpdr,  &
              &              dvtdr,dvrdt,cvr,n_r,n_store_last,n_field_type)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
      real(cp), intent(in) :: xir(:,:),phir(:,:)
      real(cp), intent(in) :: dvrdp(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvrdt(:,:)
      real(cp), intent(in) :: cvr(:,:)

      integer,  intent(in) :: n_r
      integer,  intent(in) :: n_store_last     ! Start position in frame(*)-1
      integer,  intent(in) :: n_field_type     ! Defines field type

      !-- Local variables:
      integer :: n_theta, n_theta_cal, n_phi, n_o
      real(cp) ::  fac, fac_r, fac_t


      !--- Store data for all output thetas in the current block
      !    and all output phis:

      if ( n_field_type == 1 ) then

         fac=or2(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*br(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 2 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*bt(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 3 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*bp(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 4 ) then

         fac=or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*vr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 5 ) then

         fac_r=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*vt(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 6 ) then

         fac_r=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*vp(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 7 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=sr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 109 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=xir(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 112 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=phir(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 16 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * (                                    &
               &    cosTheta(n_theta_cal)*or1(n_r)*cvr(n_theta_cal,n_phi) - &
               &                        or2(n_r)*dvrdp(n_theta_cal,n_phi) + &
               &                                 dvpdr(n_theta_cal,n_phi) - &
               &                          beta(n_r)*vp(n_theta_cal,n_phi)   )
            end do
         end do

      else if ( n_field_type == 17 ) then ! Convective heat flux

         fac=or2(n_r)*temp0(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * vr(n_theta_cal,n_phi)*sr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 113 ) then ! Chemical heat flux

         fac=or2(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * vr(n_theta_cal,n_phi)*xir(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 91 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=drSr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 18 ) then

         !--- Helicity:
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=                                        &
               &    or4(n_r)*orho2(n_r)*vr(n_theta_cal,n_phi) *          &
               &                       cvr(n_theta_cal,n_phi) +          &
               &    or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)* (   &
               &                                 vt(n_theta_cal,n_phi) * &
               &                   ( or2(n_r)*dvrdp(n_theta_cal,n_phi) - &
               &                              dvpdr(n_theta_cal,n_phi) + &
               &                    beta(n_r)*vp(n_theta_cal,n_phi)  ) + &
               &                                 vp(n_theta_cal,n_phi) * &
               &                   (          dvtdr(n_theta_cal,n_phi) - &
               &                    beta(n_r)*vt(n_theta_cal,n_phi)    - &
               &                  or2(n_r)*dvrdt(n_theta_cal,n_phi) ) )
            end do
         end do

      else if ( n_field_type == 48 ) then

         !--- Radial component of vorticity:
         fac=or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*cvr(n_theta_cal,n_phi)
            end do
         end do

      else if (n_field_type == 108) then

         !-- Cylindrically radial component of magnetic field (Bs):
         fac_r = or2(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               fac_t=or1(n_r)*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac_r*br(n_theta_cal,n_phi)*sinTheta(n_theta_cal)+ &
               &                 fac_t*bt(n_theta_cal,n_phi)*cosTheta(n_theta_cal)
            end do
         end do

      end if

   end subroutine store_fields_r
!----------------------------------------------------------------------------
   subroutine store_fields_p(vr,vt,vp,br,bp,bt,sr,drSr,xir,phir,dvrdp,         &
              &              dvpdr,dvtdr, dvrdt,cvr,cbr,cbt,n_r,n_store_last,  &
              &              n_field_type,n_phi_const,n_field_size)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces phi=const. into array frames(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
      real(cp), intent(in) :: xir(:,:),phir(:,:)
      real(cp), intent(in) :: dvrdp(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvrdt(:,:)
      real(cp), intent(in) :: cvr(:,:)
      real(cp), intent(in) :: cbr(:,:),cbt(:,:)
      integer,  intent(in) :: n_r              ! No. of radial point
      integer,  intent(in) :: n_store_last     ! Start position in frame(*)-1
      integer,  intent(in) :: n_field_type     ! Defines field type
      integer,  intent(in) :: n_phi_const      ! No. of surface phi
      integer,  intent(in) :: n_field_size     ! Size of field

      !-- Local variables:
      integer :: n_phi_0,n_phi_180,n_theta,n_theta_cal,n_0,n_180
      real(cp) ::  fac,fac_r

      !--- Get phi no. for left and right halfspheres:
      n_phi_0=n_phi_const
      if ( mod(minc,2) == 1 ) then
         n_phi_180=n_phi_max/2+n_phi_0
      else
         n_phi_180=n_phi_0
      end if
      n_0=n_store_last+(n_r-1)*n_theta_max
      n_180=n_0+n_field_size

      if ( n_field_type == 1 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=one
         else
            fac=or2(n_r)
         end if
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*br(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*br(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 2 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=one
         else
            fac=or1(n_r)
         end if
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*bt(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*bt(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 3 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=one
         else
            fac=or1(n_r)
         end if
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*bp(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*bp(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 4 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*vScale
         else
            fac=or2(n_r)*orho1(n_r)*vScale
         end if
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac*vr(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*vr(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 5 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*vScale
         else
            fac=or1(n_r)*orho1(n_r)*vScale
         end if
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*vt(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*vt(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 6 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*vScale
         else
            fac=or1(n_r)*orho1(n_r)*vScale
         end if
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*vp(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*vp(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 7 ) then

         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=sr(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=sr(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 15 ) then ! vz

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac* ( cosTheta(n_theta_cal)*or1(n_r)* &
            &                          vr(n_theta_cal,n_phi_0) -       &
            &                          vt(n_theta_cal,n_phi_0) )
            frames(n_180+n_theta)=fac* ( cosTheta(n_theta_cal)*or1(n_r)* &
            &                            vr(n_theta_cal,n_phi_180) -     &
            &                            vt(n_theta_cal,n_phi_180) )
         end do

      else if ( n_field_type == 109 ) then

         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =xir(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=xir(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 112 ) then

         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =phir(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=phir(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 16 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac * (                                    &
            &    cosTheta(n_theta_cal)*or1(n_r)*cvr(n_theta_cal,n_phi_0) - &
            &                        or2(n_r)*dvrdp(n_theta_cal,n_phi_0) + &
            &                                 dvpdr(n_theta_cal,n_phi_0) - &
            &                          beta(n_r)*vp(n_theta_cal,n_phi_0)  )
            frames(n_180+n_theta)=fac * (                                    &
            &    cosTheta(n_theta_cal)*or1(n_r)*cvr(n_theta_cal,n_phi_180) - &
            &                        or2(n_r)*dvrdp(n_theta_cal,n_phi_180) + &
            &                                 dvpdr(n_theta_cal,n_phi_180) - &
            &                          beta(n_r)*vp(n_theta_cal,n_phi_180)  )
         end do

      else if ( n_field_type == 17 ) then ! Convective heat flux

         fac=or2(n_r)*temp0(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac*vr(n_theta_cal,n_phi_0)*sr(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*vr(n_theta_cal,n_phi_180) * &
            &                     sr(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 113 ) then ! Composition heat flux

         fac=or2(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac*vr(n_theta_cal,n_phi_0)*xir(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*vr(n_theta_cal,n_phi_180) * &
            &                     xir(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 91 ) then

         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=drSr(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=drSr(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 18 ) then

         !--- Helicity:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=                                      &
            &           or4(n_r)*orho2(n_r)*vr(n_theta_cal,n_phi_0) * &
            &                              cvr(n_theta_cal,n_phi_0) + &
            &    or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)* (   &
            &                               vt(n_theta_cal,n_phi_0) * &
            &                 ( or2(n_r)*dvrdp(n_theta_cal,n_phi_0) - &
            &                            dvpdr(n_theta_cal,n_phi_0) + &
            &                beta(n_r)*   vp(n_theta_cal,n_phi_0) ) + &
            &                               vp(n_theta_cal,n_phi_0) * &
            &                 (          dvtdr(n_theta_cal,n_phi_0) - &
            &                  beta(n_r)*   vt(n_theta_cal,n_phi_0) - &
            &                   or2(n_r)*dvrdt(n_theta_cal,n_phi_0) ) )
            frames(n_180+n_theta)=                                    &
            &         or4(n_r)*orho2(n_r)*vr(n_theta_cal,n_phi_180) * &
            &                            cvr(n_theta_cal,n_phi_180) + &
            &    or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)* (   &
            &                             vt(n_theta_cal,n_phi_180) * &
            &               ( or2(n_r)*dvrdp(n_theta_cal,n_phi_180) - &
            &                          dvpdr(n_theta_cal,n_phi_180) + &
            &              beta(n_r)*   vp(n_theta_cal,n_phi_180) ) + &
            &                             vp(n_theta_cal,n_phi_180) * &
            &               (          dvtdr(n_theta_cal,n_phi_180) - &
            &                 beta(n_r)*  vt(n_theta_cal,n_phi_180) - &
            &                 or2(n_r)*dvrdt(n_theta_cal,n_phi_180) ) )
         end do

      else if ( n_field_type == 47 ) then

         !--- Phi component of vorticity:
         fac=vScale*orho1(n_r)*or1(n_r)
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=                         &
            &            fac*O_sin_theta(n_theta_cal)*   &
            &    (          dvtdr(n_theta_cal,n_phi_0) - &
            &     beta(n_r)*   vt(n_theta_cal,n_phi_0) - &
            &      or2(n_r)*dvrdt(n_theta_cal,n_phi_0) )
            frames(n_180+n_theta)=                         &
            &              fac*O_sin_theta(n_theta_cal)*   &
            &    (          dvtdr(n_theta_cal,n_phi_180) - &
            &     beta(n_r)*   vt(n_theta_cal,n_phi_180) - &
            &      or2(n_r)*dvrdt(n_theta_cal,n_phi_180) )
         end do

         !--- Phi component of Lorentz-Force:

      else if ( n_field_type == 54 ) then

         fac_r=LFfac*or3(n_r)
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fac=fac_r*O_sin_theta(n_theta_cal)
            frames(n_0+n_theta)= fac *                                &
            &    ( cbr(n_theta_cal,n_phi_0)*bt(n_theta_cal,n_phi_0) - &
            &      cbt(n_theta_cal,n_phi_0)*br(n_theta_cal,n_phi_0) )
            frames(n_180+n_theta)= fac *                                  &
            &    ( cbr(n_theta_cal,n_phi_180)*bt(n_theta_cal,n_phi_180) - &
            &      cbt(n_theta_cal,n_phi_180)*br(n_theta_cal,n_phi_180) )
         end do

      end if

   end subroutine store_fields_p
!----------------------------------------------------------------------------
   subroutine store_fields_p_avg(vr,vt,vp,bp,sr,drSr,xir,phir,dvrdp,       &
              &                  dvpdr,dvtdr, dvrdt,cvr,n_r,n_store_last,  &
              &                  n_field_type)
      !
      !  Purpose of this subroutine is to store movie frames for phi-averaged
      !  surfaces into array frames(*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
      real(cp), intent(in) :: xir(:,:),phir(:,:)
      real(cp), intent(in) :: dvrdp(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvrdt(:,:)
      real(cp), intent(in) :: cvr(:,:)
      integer,  intent(in) :: n_r              ! No. of radial point
      integer,  intent(in) :: n_store_last     ! Start position in frame(*)-1
      integer,  intent(in) :: n_field_type     ! Defines field type

      !-- Local variables:
      integer :: n_theta,n_theta2,n_theta_cal,n_phi,n_0
      real(cp) ::  phi_norm

      real(cp) ::  fl(n_theta_max) ! Field for poloidal field lines

      !--- Get phi no. for left and right halfspheres:
      n_0=n_store_last+(n_r-1)*n_theta_max
      phi_norm=one/n_phi_max

      if ( n_field_type == 8 ) then

         !--- Field for axisymmetric poloidal field lines:
         call get_fl(fl,n_r,.false.)
         do n_theta_cal=1,n_theta_max,2

            n_theta    =n_theta_cal2ord(n_theta_cal)
            n_theta2   =n_theta_cal2ord(n_theta_cal+1)
            !call get_fl(fl,n_r,n_theta_cal,1,.false.)
            frames(n_0+n_theta) =fl(n_theta_cal)
            frames(n_0+n_theta2)=fl(n_theta_cal+1)
         end do

      else if ( n_field_type == 9 ) then

         !--- Axisymmetric B_phi:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+bp(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*or1(n_r)*O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 10 ) then

         !--- Field for axisymmetric velocity stream lines:
         call get_sl(fl,n_r)
         do n_theta_cal=1,n_theta_max,2

            n_theta =n_theta_cal2ord(n_theta_cal)
            n_theta2=n_theta_cal2ord(n_theta_cal+1)
            frames(n_0+n_theta) =fl(n_theta_cal)
            frames(n_0+n_theta2)=fl(n_theta_cal+1)
         end do

      else if ( n_field_type == 11 ) then

         !--- Axisymmetric v_phi:
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+orho1(n_r)*vp(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*or1(n_r)*O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 12 ) then

         !--- Axisymmetric T:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+sr(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 110 ) then

         !--- Axisymmetric C:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+xir(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 111 ) then

         !--- Axisymmetric phase field:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+phir(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 92 ) then

         !--- Axisymmetric dsdr:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+drSr(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 115 ) then

         !--- Axisymmetric convective heat flux:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+temp0(n_r)*or2(n_r)*vr(n_theta_cal,n_phi)* &
               &     sr(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 116 ) then

         !--- Axisymmetric flux of chemical composition
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+temp0(n_r)*or2(n_r)*vr(n_theta_cal,n_phi)* &
               &     xir(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 114 ) then

         !--- Kinetic energy
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+half*orho1(n_r)*or2(n_r)*(                         &
               &      or2(n_r)* vr(n_theta_cal,n_phi)*vr(n_theta_cal,n_phi) + &
               &     O_sin_theta_E2(n_theta_cal)*                             &
               &                vt(n_theta_cal,n_phi)*vt(n_theta_cal,n_phi) + &
               &     O_sin_theta_E2(n_theta_cal)*                             &
               &                vp(n_theta_cal,n_phi)*vp(n_theta_cal,n_phi) )
            end do
            frames(n_0+n_theta) =phi_norm*fl(1)
         end do

      else if ( n_field_type == 94 ) then

         !--- Axisymmetric v_s=sin(theta)*v_r+cos(theta)*v_theta
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+sinTheta(n_theta_cal)*or1(n_r)*  vr(n_theta_cal,n_phi)+ &
               &     cosTheta(n_theta_cal)*O_sin_theta(n_theta_cal)*               &
               &                                            vt(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*orho1(n_r)*or1(n_r)
         end do

      else if ( n_field_type == 95 ) then

         !--- Axisymmetric v_s*v_phi
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+or1(n_r)*  vr(n_theta_cal,n_phi)*vp(n_theta_cal,n_phi) + &
               &     cosTheta(n_theta_cal)*O_sin_theta_E2(n_theta_cal)*             &
               &                      vt(n_theta_cal,n_phi)*vp(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*orho2(n_r)*or2(n_r)
         end do

      else if ( n_field_type == 96 ) then

         !--- Axisymmetric v_z=cos(theta)*v_r-sin(theta)*v_theta
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+cosTheta(n_theta_cal)*or1(n_r)*  vr(n_theta_cal,n_phi)- &
               &                                            vt(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*orho1(n_r)*or1(n_r)
         end do

      else if ( n_field_type == 97 ) then

         !--- Axisymmetric v_z*v_phi
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+or1(n_r)*cosTheta(n_theta_cal)*O_sin_theta(n_theta_cal)*&
               &                      vr(n_theta_cal,n_phi)*vp(n_theta_cal,n_phi) -&
               & O_sin_theta(n_theta_cal)*vt(n_theta_cal,n_phi)*vp(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*orho2(n_r)*or2(n_r)
         end do

      else if ( n_field_type == 98 ) then

         !--- Axisymmetric v_s**2
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+sinTheta(n_theta_cal)*sinTheta(n_theta_cal)*or2(n_r)* &
               &                    vr(n_theta_cal,n_phi)*vr(n_theta_cal,n_phi)+ &
               &                    cosTheta(n_theta_cal)*cosTheta(n_theta_cal)* &
               &                                    O_sin_theta_E2(n_theta_cal)* &
               &                    vt(n_theta_cal,n_phi)*vt(n_theta_cal,n_phi)+ &
               &                             two*cosTheta(n_theta_cal)*or1(n_r)* &
               &                        vr(n_theta_cal,n_phi)*vt(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*orho2(n_r)*or2(n_r)
         end do


      else if ( n_field_type == 99 ) then

         !--- Axisymmetric v_z**2
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max   ! Average over phis
               fl(1)=fl(1)+cosTheta(n_theta_cal)*cosTheta(n_theta_cal)*or2(n_r)*  &
               &                    vr(n_theta_cal,n_phi)*vr(n_theta_cal,n_phi)+  &
               &                    vt(n_theta_cal,n_phi)*vt(n_theta_cal,n_phi)-  &
               &                             two*cosTheta(n_theta_cal)*or1(n_r)*  &
               &                        vr(n_theta_cal,n_phi)*vt(n_theta_cal,n_phi)
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)*orho2(n_r)*or2(n_r)
         end do

      else if ( n_field_type == 19 ) then

         !--- Axisymmetric helicity:
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            fl(1)=0.0_cp
            do n_phi=1,n_phi_max
               fl(1)=fl(1) +                                            &
               &            or4(n_r)*orho2(n_r)*vr(n_theta_cal,n_phi) * &
               &                               cvr(n_theta_cal,n_phi) + &
               &    or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)*(   &
               &                                vt(n_theta_cal,n_phi) * &
               &                  ( or2(n_r)*dvrdp(n_theta_cal,n_phi) - &
               &                             dvpdr(n_theta_cal,n_phi) + &
               &                 beta(n_r)*   vp(n_theta_cal,n_phi) ) + &
               &                                vp(n_theta_cal,n_phi) * &
               &                    (        dvtdr(n_theta_cal,n_phi) - &
               &                   beta(n_r)*   vt(n_theta_cal,n_phi) - &
               &                    or2(n_r)*dvrdt(n_theta_cal,n_phi) ) )
            end do
            frames(n_0+n_theta)=phi_norm*fl(1)
         end do

      end if

   end subroutine store_fields_p_avg
!----------------------------------------------------------------------------
   subroutine store_fields_t(vr,vt,vp,br,bt,bp,sr,drSr,xir,phir,dvrdp,    &
              &              dvpdr,dvtdr,dvrdt,cvr,cbt,n_r,n_store_last,  &
              &              n_field_type,n_theta)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
      real(cp), intent(in) :: xir(:,:),phir(:,:)
      real(cp), intent(in) :: dvrdp(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvrdt(:,:)
      real(cp), intent(in) :: cvr(:,:),cbt(:,:)
      integer,  intent(in) :: n_r              ! No. of radial grid point
      integer,  intent(in) :: n_store_last     ! Position in frame(*)-1
      integer,  intent(in) :: n_field_type     ! Defines field
      integer,  intent(in) :: n_theta          ! No. of theta in block

      !-- Local variables:
      integer :: n_phi, n_o
      real(cp) ::  fac


      n_o=n_store_last+(n_r-1)*n_phi_max

      if ( n_field_type == 1 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=one
         else
            fac=or2(n_r)
         end if
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*br(n_theta,n_phi)
         end do

      else if ( n_field_type == 2 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=O_sin_theta(n_theta)
         else
            fac=or1(n_r)*O_sin_theta(n_theta)
         end if
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*bt(n_theta,n_phi)
         end do

      else if ( n_field_type == 3 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=O_sin_theta(n_theta)
         else
            fac=or1(n_r)*O_sin_theta(n_theta)
         end if
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*bp(n_theta,n_phi)
         end do

      else if ( n_field_type == 4 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*vScale
         else
            fac=or2(n_r)*orho1(n_r)*vScale
         end if
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vr(n_theta,n_phi)
         end do

      else if ( n_field_type == 5 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*O_sin_theta(n_theta)*vScale
         else
            fac=or1(n_r)*orho1(n_r)*O_sin_theta(n_theta)*vScale
         end if
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vt(n_theta,n_phi)
         end do

      else if ( n_field_type == 6 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*O_sin_theta(n_theta)*vScale
         else
            fac=or1(n_r)*orho1(n_r)*O_sin_theta(n_theta)*vScale
         end if
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vp(n_theta,n_phi)
         end do

      else if ( n_field_type == 7 ) then

         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=sr(n_theta,n_phi)
         end do

      else if ( n_field_type == 109 ) then

         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=xir(n_theta,n_phi)
         end do

      else if ( n_field_type == 112 ) then

         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=phir(n_theta,n_phi)
         end do

      else if ( n_field_type == 13 ) then

         fac=-or1(n_r)*O_sin_theta(n_theta)
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*bt(n_theta,n_phi)
         end do

      else if ( n_field_type == 14 ) then

         fac=-or1(n_r)*O_sin_theta(n_theta)
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*cbt(n_theta,n_phi)
         end do

      else if ( n_field_type == 15 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac* ( cosTheta(n_theta)*or1(n_r)* &
            &                  vr(n_theta,n_phi) - vt(n_theta,n_phi) )
         end do

      else if ( n_field_type == 16 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            frames(n_phi+n_o)=fac * (                                  &
            &          cosTheta(n_theta)*or1(n_r)*cvr(n_theta,n_phi) - &
            &                          or2(n_r)*dvrdp(n_theta,n_phi) + &
            &                                   dvpdr(n_theta,n_phi) - &
            &                         beta(n_r)*   vp(n_theta,n_phi) )
         end do

      else if ( n_field_type == 17 ) then ! Convective heat flux

         fac=or2(n_r)*temp0(n_r)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vr(n_theta,n_phi)*sr(n_theta,n_phi)
         end do

      else if ( n_field_type == 113 ) then ! Chemical heat flux

         fac=or2(n_r)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vr(n_theta,n_phi)*xir(n_theta,n_phi)
         end do

      else if ( n_field_type == 91 ) then

         fac=vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=drSr(n_theta,n_phi)
         end do

      else if ( n_field_type == 18 ) then

         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=                                       &
            &                or4(n_r)*orho2(n_r)*vr(n_theta,n_phi) * &
            &                                   cvr(n_theta,n_phi) + &
            &          or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta)*( &
            &                                    vt(n_theta,n_phi) * &
            &                      ( or2(n_r)*dvrdp(n_theta,n_phi) - &
            &                                 dvpdr(n_theta,n_phi) + &
            &                     beta(n_r)*   vp(n_theta,n_phi) ) + &
            &                                    vp(n_theta,n_phi) * &
            &                      (          dvtdr(n_theta,n_phi) - &
            &                       beta(n_r)*   vt(n_theta,n_phi) - &
            &                        or2(n_r)*dvrdt(n_theta,n_phi) ) )
         end do

      end if

   end subroutine store_fields_t
!----------------------------------------------------------------------------
   subroutine store_fields_3d(vr,vt,vp,br,bt,bp,sr,drSr,xir,phir,dvrdp,  &
              &               dvpdr,dvtdr,dvrdt,cvr,cbr,cbt,n_r,         &
              &               n_store_last,n_field_type)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
      real(cp), intent(in) :: xir(:,:),phir(:,:)
      real(cp), intent(in) :: dvrdp(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvtdr(:,:),dvrdt(:,:)
      real(cp), intent(in) :: cvr(:,:)
      real(cp), intent(in) :: cbr(:,:),cbt(:,:)

      integer,  intent(in) :: n_r              ! No. of radial grid point
      integer,  intent(in) :: n_store_last     ! Position in frame(*)-1
      integer,  intent(in) :: n_field_type     ! Defines field

      !-- Local variables:
      integer :: n_phi,n_theta,n_theta_cal,n_o,n_or
      real(cp) ::  fac,fac_r


      n_or=n_store_last+(n_r-1)*n_theta_max*n_phi_max

      if ( n_field_type == 1 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=one
         else
            fac=or2(n_r)
         end if
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*br(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 2 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac_r=one
         else
            fac_r=or1(n_r)
         end if
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*bt(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 3 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac_r=one
         else
            fac_r=or1(n_r)
         end if
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*bp(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 4 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac=orho1(n_r)*vScale
         else
            fac=or2(n_r)*orho1(n_r)*vScale
         end if
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*vr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 5 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac_r=orho1(n_r)*vScale
         else
            fac_r=or1(n_r)*orho1(n_r)*vScale
         end if
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*vt(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 6 ) then

         if ( n_r == n_r_icb .and. l_full_sphere ) then
            fac_r=orho1(n_r)*vScale
         else
            fac_r=or1(n_r)*orho1(n_r)*vScale
         end if
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*vp(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 7 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=sr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 109 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=xir(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 112 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=phir(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 15 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * (                                   &
               &    cosTheta(n_theta_cal)*or1(n_r)*vr(n_theta_cal,n_phi) - &
               &    vt(n_theta_cal,n_phi) )
            end do
         end do

      else if ( n_field_type == 16 ) then

         !--- Z-component of helicity
         fac=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * (                                    &
               &    cosTheta(n_theta_cal)*or1(n_r)*cvr(n_theta_cal,n_phi) - &
               &                        or2(n_r)*dvrdp(n_theta_cal,n_phi) + &
               &                                 dvpdr(n_theta_cal,n_phi) - &
               &                       beta(n_r)*   vp(n_theta_cal,n_phi) )
            end do
         end do

      else if ( n_field_type == 17 ) then

         fac=or2(n_r)*temp0(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*vr(n_theta_cal,n_phi)*sr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 113 ) then

         fac=or2(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*vr(n_theta_cal,n_phi)*xir(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 91 ) then

         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=drSr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 18 ) then

         !--- Helicity
         fac=vScale*vScale*or2(n_r)*orho2(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * (                   &
               &          or2(n_r)*vr(n_theta_cal,n_phi) * &
               &                  cvr(n_theta_cal,n_phi) + &
               &          O_sin_theta_E2(n_theta_cal)* (   &
               &                   vt(n_theta_cal,n_phi) * &
               &     ( or2(n_r)*dvrdp(n_theta_cal,n_phi) - &
               &                dvpdr(n_theta_cal,n_phi) + &
               &    beta(n_r)*   vp(n_theta_cal,n_phi) ) + &
               &                   vp(n_theta_cal,n_phi) * &
               &     (          dvtdr(n_theta_cal,n_phi) - &
               &      beta(n_r)*   vt(n_theta_cal,n_phi) - &
               &       or2(n_r)*dvrdt(n_theta_cal,n_phi) ) ) )
            end do
         end do

      else if ( n_field_type == 47 ) then

         !--- Phi component of vorticity
         fac=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac *                               &
               &    O_sin_theta(n_theta_cal)*vp(n_theta_cal,n_phi) * &
               &               (          dvtdr(n_theta_cal,n_phi) - &
               &                beta(n_r)*   vt(n_theta_cal,n_phi) - &
               &                 or2(n_r)*dvrdt(n_theta_cal,n_phi) )
            end do
         end do

      else if ( n_field_type == 48 ) then

         !--- Radial component of vorticity:
         fac=or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*cvr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 53 ) then

         !--- Omega effect: Br*dvp/dr
         fac_r=or3(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=         fac*br(n_theta_cal,n_phi) * &
               &                         ( dvpdr(n_theta_cal,n_phi) - &
               &    (beta(n_r)+two*or1(n_r))*vp(n_theta_cal,n_phi) )
            end do
         end do

      else if ( n_field_type == 54 ) then

         !--- Phi component of Lorentz-Force:

         fac_r=LFfac*or3(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)= fac *                              &
               &    ( cbr(n_theta_cal,n_phi)*bt(n_theta_cal,n_phi) - &
               &      cbt(n_theta_cal,n_phi)*br(n_theta_cal,n_phi) )
            end do
         end do

      end if

   end subroutine store_fields_3d
!----------------------------------------------------------------------------
   subroutine get_sl(sl,n_r)
      !
      !  Return field sl whose contourlines are the stream lines
      !  of the axisymmetric poloidal velocity field.
      !
      !  .. math::
      !     s(r,\theta) = \dfrac{1}{r}\dfrac{\partial}{\partial \theta} u(r,\theta,m=0)
      !

      !-- Input variables:
      integer, intent(in) :: n_r             ! No. of radial grid point

      !-- Output variables:
      real(cp), intent(out) ::  sl(:)           ! Field for field lines

      !-- Local variables:
      integer :: n_theta         ! No. of theta
      integer :: l,lm            ! Degree, counter for degree/order combinations
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)
      complex(cp) :: Tl_AX(1:l_max+1)

      Tl_AX(1)=zero
      do l=1,l_max
         lm=lm2(l,0)
         Tl_AX(l+1)=or1(n_r)*w_Rloc(lm,n_r)
      end do

      call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:), l_R(n_r))

      do n_theta=1,n_theta_max
         sl(n_theta)=O_sin_theta(n_theta)*tmpp(n_theta)
      end do

   end subroutine get_sl
!----------------------------------------------------------------------------
   subroutine get_fl(fl,n_r,l_ic)
      !
      !  Return field fl whose contourlines are the fields lines
      !  of the axisymmetric poloidal mangetic field.
      !
      !  .. math::
      !     f(r,\theta) = \dfrac{1}{r}\dfrac{\partial}{\partial \theta} b(r,\theta,m=0)
      !
      !  This routine is called for l_ic=.true. only from rank 0 with full
      !  field b_ic in standard lm ordering available.
      !  The case l_ic=.false. is called from all ranks and uses b_Rloc.

      !-- Input variables:
      integer, intent(in) :: n_r             ! No. of radial grid point
      logical, intent(in) :: l_ic            ! =true if inner core field

      !-- Output variables:
      real(cp), intent(out) ::  fl(:)    ! Field for field lines

      !-- Local variables:
      integer :: n_theta         ! No. of theta
      integer :: l,lm            ! Degree, counter for degree/order combinations

      real(cp) :: r_ratio          ! r/r_ICB
      real(cp) :: r_dep(l_max)     ! (r/r_ICB)**l / r_ICB
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)
      complex(cp) :: Tl_AX(1:l_max+1)

      if ( l_ic ) then
         r_ratio =r_ic(n_r)/r_ic(1)
         r_dep(1)=r_ratio/r_ic(1)
         do l=2,l_max
            r_dep(l)=r_dep(l-1)*r_ratio
         end do
      end if

      Tl_AX(1)=zero
      do l=1,l_max
         lm=lm2(l,0)
         if ( l_ic ) then ! Inner Core
            if ( l_cond_ic ) then
               Tl_AX(l+1)=r_dep(l)*b_ic(lm,n_r)
            else
               Tl_AX(l+1)=r_dep(l)*bICB(lm)
            end if
         else             ! Outer Core
            Tl_AX(l+1)=or1(n_r)*b_Rloc(lm,n_r)
         end if
      end do

      if ( l_ic ) then
         call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:), l_max)
      else
         call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:), l_R(n_r))
      end if

      do n_theta=1,n_theta_max
         fl(n_theta)=O_sin_theta(n_theta)*tmpp(n_theta)
      end do

   end subroutine get_fl
!----------------------------------------------------------------------------
   subroutine get_B_surface(b_r,b_t,b_p,bCMB)
      !
      !  Upward continuation of laplacian field to Earths surface.
      !  Field is given by poloidal harmonic coefficients b at CMB.
      !  Spherical harmonic transforms of upward continued field
      !  to r/theta/phi vector components for all logitudes and
      !  latitude are returned in br/bt/bp.
      !  Note that this routine given the real components of the magnetic
      !  fields while other transforms in the code provide only:
      !  :math:`r^2 B_r`, :math:`r^2 \sin\theta B_\theta`,
      !  :math:`r^2 \sin\theta B_\phi`
      !

      !-- Input of variables:
      complex(cp), intent(in) :: bCMB(:)

      !-- Output:
      real(cp), intent(out) :: b_r(:,:) !Radial magnetic field in (theta,phi)-space
      real(cp), intent(out) :: b_t(:,:) !Latitudinal magnetic field
      real(cp), intent(out) :: b_p(:,:) !Azimuthal magnetic field.

      !-- Local variables:
      integer :: l,lm

      real(cp) :: r_ratio          ! r_cmb/r_surface
      real(cp) :: r_dep(l_max)     ! Radial dependence
      complex(cp) :: cs1(lm_max),cs2(lm_max) ! help arrays
      complex(cp) :: zerosc(lm_max)

      !-- Radial dependence:
      r_ratio=r_cmb/r_surface
      r_dep(1)=r_ratio/(r_surface*r_surface)  ! l=1 term
      do l=2,l_max
         r_dep(l)=r_dep(l-1)*r_ratio
      end do

      !-- Construct help arrays containing radial dependence
      !   and l dependence: dLh=l*(l+1)
      cs1(1)=zero
      cs2(1)=zero
      do lm=2,lm_max
         l = lm2l(lm)
         cs1(lm) = bCMB(lm)*r_dep(l) ! multiplication by l(l+1) in shtns.f90
         cs2(lm)= -bCMB(lm)*real(l,cp)*r_dep(l)
      end do

      zerosc(:)=zero
      call torpol_to_spat(cs1, cs2, zerosc, b_r, b_t, b_p, l_max)

   end subroutine get_B_surface
!----------------------------------------------------------------------------
end module out_movie
