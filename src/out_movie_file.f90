module out_movie

   use precision_mod
   use parallel_mod, only: coord_r, l_master_rank
   use communications, only: gather_all_from_mlo_to_master, gather_f, gather_Flm
   !@> TODO remove gather_f once finished
   use truncation, only: n_phi_max, n_theta_max, minc, lm_max, l_max,       &
       &                 n_m_max, lm_maxMag, n_r_maxMag, n_r_ic_maxMag,     &
       &                 n_r_ic_max, n_r_max, l_axi, n_r_icb, nThetaStart,  &
       &                 nThetaStop, n_lm_loc, n_mloMag_loc, n_r_cmb,       &
       &                 nlat_padded
   use movie_data, only: frames, n_movie_fields, n_movies, n_movie_surface, &
       &                 n_movie_const, n_movie_field_type,                 &
       &                 n_movie_field_start,n_movie_field_stop,            &
       &                 movieDipColat, movieDipLon, movieDipStrength,      &
       &                 movieDipStrengthGeo, t_movieS, n_movie_type,       &
       &                 lStoreMov, n_movie_file, n_movie_fields_ic,        &
       &                 movie_file, movie_const
   use radial_functions, only: orho1, orho2, or1, or2, or3, or4, beta,  &
       &                       r_surface, r_cmb, r, r_ic
   use physical_parameters, only: LFfac, radratio, ra, ek, pr, prmag
   use num_param, only: vScale, tScale
   use LMmapping, only: map_glbl_st
   use horizontal_data, only: O_sin_theta, sinTheta, cosTheta,    &
       &                      n_theta_cal2ord, O_sin_theta_E2,    &
       &                      osn1, phi, theta_ord
   use fields, only: b_ic, bICB, w_Rdist, b_Rdist
   use sht, only: torpol_to_spat_loc, toraxi_to_spat
   use logic, only: l_save_out, l_cond_ic, l_mag
   use constants, only: zero, one, two
   use out_dtB_frame, only: write_dtB_frame
   use output_data, only: runid
   use useful, only: abortRun

   implicit none

   private

   public :: store_movie_frame, write_movie_frame, get_fl

contains

   subroutine store_movie_frame(n_r,vr_loc,vt_loc,vp_loc,br_loc,bt_loc,     &
              &                 bp_loc,sr_loc,drSr_loc,dvrdp_loc,dvpdr_loc, &
              &                 dvtdr_loc,dvrdt_loc,cvr_loc,cbr_loc,cbt_loc)
      !
      !  Controls output of movie frames.
      !  Usually called from radialLoop.
      !

      !-- Input variables:
      integer,     intent(in) :: n_r                ! radial grid point no.
      real(cp),    intent(in) :: vr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: vt_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: vp_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: br_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: bt_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: bp_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: sr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: drSr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: dvrdp_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: dvpdr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: dvtdr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: dvrdt_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: cvr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: cbr_loc(n_phi_max,nThetaStart:nThetaStop)
      real(cp),    intent(in) :: cbt_loc(n_phi_max,nThetaStart:nThetaStop)

      !-- Local variables:
      integer :: n_movie        ! No. of movie
      integer :: n_field        ! No. of field
      integer :: n_surface      ! Surface (1=r,2=theta,3=phi)
      integer :: n_const        ! Gives surface
      integer :: n_field_type   ! Numbers field types
      integer :: n_store_last   ! Position i in frame(i) were field starts-1
      integer :: n_theta,n_field_size,n_fields
      logical :: lThetaFound

      !@> TODO remove useless local arrays:
      real(cp) :: vr(n_theta_max,n_phi_max), vt(n_theta_max,n_phi_max)
      real(cp) :: vp(n_theta_max,n_phi_max), sr(n_theta_max,n_phi_max)
      real(cp) :: cvr(n_theta_max,n_phi_max), br(n_theta_max,n_phi_max)
      real(cp) :: bt(n_theta_max,n_phi_max), bp(n_theta_max,n_phi_max)
      real(cp) :: cbr(n_theta_max,n_phi_max), cbt(n_theta_max,n_phi_max)
      real(cp) :: drSr(n_theta_max,n_phi_max), dvrdp(n_theta_max,n_phi_max)
      real(cp) :: dvpdr(n_theta_max,n_phi_max), dvtdr(n_theta_max,n_phi_max)
      real(cp) :: dvrdt(n_theta_max,n_phi_max)
      complex(cp) :: bCMB(lm_max)

      !@> TODO remove the gather_f from here at some point
      call gather_f(vr_loc, vr)
      call gather_f(vt_loc, vt)
      call gather_f(vp_loc, vp)
      call gather_f(br_loc, br)
      call gather_f(bt_loc, bt)
      call gather_f(bp_loc, bp)
      call gather_f(sr_loc, sr)
      call gather_f(drSr_loc, drSr)
      call gather_f(dvrdp_loc, dvrdp)
      call gather_f(dvpdr_loc, dvpdr)
      call gather_f(dvtdr_loc, dvtdr)
      call gather_f(dvrdt_loc, dvrdt)
      call gather_f(cvr_loc, cvr)
      call gather_f(cbr_loc, cbr)
      call gather_f(cbt_loc, cbt)
      if ( n_r == n_r_cmb ) call gather_Flm(b_Rdist(:,n_r_cmb), bCMB)
      !@> END suppression blocks once DONE...

      do n_movie=1,n_movies

         n_fields =n_movie_fields(n_movie)
         n_surface=n_movie_surface(n_movie)
         n_const  =n_movie_const(n_movie)

         if ( n_surface == -1 ) then ! Earth Surface

            if ( n_r /= 1 ) cycle  ! not CMB radius

            do n_field=1,n_fields
               n_field_type=n_movie_field_type(n_field,n_movie)
               n_store_last=n_movie_field_start(n_field,n_movie)-1
               if ( n_store_last >= 0 ) then
                  call store_fields_sur(n_store_last,n_field_type,bCMB)
               end if
            end do

         else if ( n_surface == 0 ) then ! 3d

            do n_field=1,n_fields
               n_field_type=n_movie_field_type(n_field,n_movie)
               n_store_last=n_movie_field_start(n_field,n_movie)-1
               if ( n_store_last >= 0 ) then
                  call store_fields_3d(vr,vt,vp,br,bt,bp,sr,drSr,           &
                       &               dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
                       &               n_r,n_store_last,n_field_type)
               end if
            end do

         else if ( n_surface == 1 ) then ! Surface r=constant

            if ( n_r /= n_const ) cycle  ! not desired radius

            do n_field=1,n_fields
               n_field_type=n_movie_field_type(n_field,n_movie)
               n_store_last=n_movie_field_start(n_field,n_movie)-1
               if ( n_store_last >= 0 ) then
                  call store_fields_r(vr,vt,vp,br,bt,bp,sr,drSr,     &
                       &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,   &
                       &              n_r,n_store_last,n_field_type)
               end if
            end do

         else if ( n_surface == 2 ) then ! Surface theta=constant

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
                  call store_fields_t(vr,vt,vp,br,bt,bp,sr,drSr,       &
                       &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbt, &
                       &              n_r,n_store_last,n_field_type,   &
                       &              n_theta)
                 end if
            end do

         else if ( abs(n_surface) == 3 ) then  ! Surface phi=const.

            do n_field=1,n_fields
               n_field_type=n_movie_field_type(n_field,n_movie)
               n_store_last=n_movie_field_start(n_field,n_movie)-1
               n_field_size=(n_movie_field_stop(n_field,n_movie) - n_store_last)/2
               if ( n_store_last >= 0 ) then
                  call store_fields_p(vr,vt,vp,br,bp,bt,sr,drSr,           &
                       &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
                       &              n_r,n_store_last,n_field_type,       &
                       &              n_const,n_field_size)
               end if
            end do  ! Do loop over field for one movie


         end if  ! Surface ?

      end do  ! Do loop over movies !

   end subroutine store_movie_frame
!----------------------------------------------------------------------------
   subroutine write_movie_frame(n_frame,time,b_LMdist,db_LMdist,aj_LMdist,   &
              &                 dj_LMdist,b_ic,db_ic,aj_ic,dj_ic,omega_ic, &
              &                 omega_ma)
      !
      !  Writes different movie frames into respective output files.
      !  Called from coord_r 0 with full arrays in standard LM order.
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time
      integer,     intent(in) :: n_frame
      real(cp),    intent(in) :: omega_ic,omega_ma
      complex(cp), intent(in) :: b_LMdist(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: db_LMdist(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: aj_LMdist(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: dj_LMdist(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(lm_maxMag,n_r_ic_maxMag)

      !-- Local variables:
      integer :: n_fields
      integer :: n_fields_ic
      integer :: n_fields_oc
      integer :: n_movie
      integer :: n_surface
      integer :: n_type
      integer :: n_out
      integer :: n_field,n,n_start,n_stop
      integer :: n_r,n_theta,n_phi
      integer :: n_r_mov_tot
      logical :: l_dtB_frame
      character(len=64) :: version
      real(cp) :: const
      real(cp) :: r_mov_tot(n_r_max+n_r_ic_max)
      real(outp) :: dumm(n_theta_max)
      complex(cp), allocatable :: b(:,:), aj(:,:), db(:,:), dj(:,:)

      l_dtB_frame = .false.

      do n_movie=1,n_movies
         n_type=n_movie_type(n_movie)
         if ( (.not. lStoreMov(n_movie)) .and. (n_type /= 99) ) then
            l_dtB_frame = .true.
         end if
      end do

      if ( l_dtB_frame ) then
         if ( l_master_rank ) then
            allocate( b(lm_maxMag,n_r_maxMag), aj(lm_maxMag,n_r_maxMag) )
            allocate( db(lm_maxMag,n_r_maxMag), dj(lm_maxMag,n_r_maxMag) )
         else
            allocate( b(1,1), aj(1,1), db(1,1), dj(1,1) )
         end if

         call gather_all_from_mlo_to_master(b_LMdist,b,n_r_max)
         call gather_all_from_mlo_to_master(db_LMdist,db,n_r_max)
         call gather_all_from_mlo_to_master(aj_LMdist,aj,n_r_max)
         call gather_all_from_mlo_to_master(dj_LMdist,dj,n_r_max)
      end if

      if ( l_master_rank ) then

         t_movieS(n_frame)=time

         do n_movie=1,n_movies

            n_type     =n_movie_type(n_movie)
            n_surface  =n_movie_surface(n_movie)
            n_fields_oc=n_movie_fields(n_movie)
            n_fields_ic=n_movie_fields_ic(n_movie)
            n_fields   =max(n_fields_ic,n_fields_oc)
            n_out      =n_movie_file(n_movie)
            const      =movie_const(n_movie)
            if ( n_surface == 1 ) const=const/r_cmb

            !------ Open movie file:
            if ( l_save_out ) then
               open(newunit=n_out, file=movie_file(n_movie), status='unknown', &
               &    form='unformatted', position='append')
            end if

            !------ Write header if this is the first frame:

            if ( n_frame == 1 ) then

               !------ Start with info about movie type:
               version='JW_Movie_Version_2'
               write(n_out) version
               write(n_out) real(n_type,kind=outp), real(n_surface,kind=outp), &
               &            real(const,kind=outp), real(n_fields,kind=outp)
               write(n_out) (real(n_movie_field_type(n,n_movie),kind=outp),n=1,n_fields)

               !------ Combine OC and IC radial grid points:
               n_r_mov_tot=n_r_max
               do n_r=1,n_r_max
                  r_mov_tot(n_r)=r(n_r)
               end do
               if ( n_r_ic_max > 0 ) then
                  n_r_mov_tot=n_r_mov_tot+n_r_ic_max-2
                  do n_r=1,n_r_ic_max-2
                     r_mov_tot(n_r_max+n_r)=r_ic(n_r+1)
                  end do
               end if

               !------ Now other info about grid and parameters:
               write(n_out) runid          ! run identifyer (as set in namelist contrl)
               dumm( 1)=real(n_r_mov_tot,kind=outp)
               dumm( 2)=real(n_r_max,kind=outp)
               dumm( 3)=real(n_theta_max,kind=outp) ! no. of theta points
               dumm( 4)=real(n_phi_max,kind=outp)   ! no. of phi points
               dumm( 5)=real(minc,kind=outp)        ! imposed symmetry
               dumm( 6)=real(ra,kind=outp)          ! control parameters
               dumm( 7)=real(ek,kind=outp)          ! (for information only)
               dumm( 8)=real(pr,kind=outp)          !      -"-
               dumm( 9)=real(prmag,kind=outp)       !      -"-
               dumm(10)=real(radratio,kind=outp)    ! ratio of inner / outer core
               dumm(11)=real(tScale,kind=outp)      ! timescale
               write(n_out) (dumm(n),n=1,11)

               !------ Write grid:
               write(n_out) (real(r_mov_tot(n_r)/r_cmb,kind=outp), n_r=1,n_r_mov_tot)
               write(n_out) (real(theta_ord(n_theta),kind=outp), n_theta=1,n_theta_max)
               write(n_out) (real(phi(n_phi),kind=outp), n_phi=1,n_phi_max)

            end if  ! Write header ?

            !------ Write frame number, time and IC and MA rotation rates::
            dumm(1)=real(n_frame,kind=outp)
            dumm(2)=real(t_movieS(n_frame),kind=outp)
            dumm(3)=real(omega_ic,kind=outp)
            dumm(4)=real(omega_ma,kind=outp)
            dumm(5)=real(movieDipColat,kind=outp)
            dumm(6)=real(movieDipLon,kind=outp)
            dumm(7)=real(movieDipStrength,kind=outp)
            dumm(8)=real(movieDipStrengthGeo,kind=outp)
            write(n_out) (dumm(n),n=1,8)

            !------ Write frames:
            if ( .not. lStoreMov(n_movie) ) then
               if ( n_type == 99 ) then
                  call abortRun('! Use TO output for Lorentz force!')
               else
                  call write_dtB_frame(n_movie,b,db,aj,dj,b_ic,db_ic,aj_ic,dj_ic)
               end if
            else
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

            if ( l_save_out ) close(n_out)

         end do  ! Loop over movies

      end if ! coord_r 0

      if ( l_dtB_frame ) deallocate( b, aj, db, dj )

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
   subroutine store_fields_r(vr,vt,vp,br,bt,bp,sr,drSr,dvrdp,dvpdr,dvtdr,  &
              &              dvrdt,cvr,n_r,n_store_last,n_field_type)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
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

      else if ( n_field_type == 17 ) then

         fac=-or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_store_last+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac * vr(n_theta_cal,n_phi)*drSr(n_theta_cal,n_phi)
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
   subroutine store_fields_p(vr,vt,vp,br,bp,bt,sr,drSr,dvrdp,dvpdr,dvtdr, &
              &              dvrdt,cvr,cbr,cbt,n_r,n_store_last,          &
              &              n_field_type,n_phi_const,n_field_size)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces phi=const. into array frames(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
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
      integer :: n_phi_0,n_phi_180,n_theta,n_theta2
      integer :: n_theta_cal,n_phi,n_0,n_180
      real(cp) ::  phi_norm,fac,fac_r

      real(cp) ::  fl(n_theta_max) ! Field for poloidal field lines

      !--- Get phi no. for left and right halfspheres:
      n_phi_0=n_phi_const
      if ( mod(minc,2) == 1 ) then
         n_phi_180=n_phi_max/2+n_phi_0
      else
         n_phi_180=n_phi_0
      end if
      n_0=n_store_last+(n_r-1)*n_theta_max
      n_180=n_0+n_field_size
      phi_norm=one/n_phi_max

      if ( n_field_type == 1 ) then

         fac=or2(n_r)
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*br(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*br(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 2 ) then

         fac=or1(n_r)
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*bt(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*bt(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 3 ) then

         fac=or1(n_r)
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*bp(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*bp(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 4 ) then

         fac=or2(n_r)*orho1(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac*vr(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*vr(n_theta_cal,n_phi_180)
         end do

      else if ( n_field_type == 5 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)  =fac*vt(n_theta_cal,n_phi_0)*   &
            &                     O_sin_theta(n_theta_cal)
            frames(n_180+n_theta)=fac*vt(n_theta_cal,n_phi_180)* &
            &                     O_sin_theta(n_theta_cal)
         end do

      else if ( n_field_type == 6 ) then

         fac=or1(n_r)*orho1(n_r)*vScale
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

      else if ( n_field_type == 8 ) then

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

      else if ( n_field_type == 17 ) then

         fac=-or2(n_r)*orho1(n_r)*vScale
         do n_theta_cal=1,n_theta_max
            n_theta=n_theta_cal2ord(n_theta_cal)
            frames(n_0+n_theta)=fac*vr(n_phi_0,n_theta)*drSr(n_theta_cal,n_phi_0)
            frames(n_180+n_theta)=fac*vr(n_phi_180,n_theta) * &
            &                     drSr(n_theta_cal,n_phi_180)
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
   subroutine store_fields_t(vr,vt,vp,br,bt,bp,sr,drSr,dvrdp,dvpdr,dvtdr, &
              &              dvrdt,cvr,cbt,n_r,n_store_last,n_field_type, &
              &              n_theta)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
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

         fac=or2(n_r)
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*br(n_theta,n_phi)
         end do

      else if ( n_field_type == 2 ) then

         fac=or1(n_r)*O_sin_theta(n_theta)
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*bt(n_theta,n_phi)
         end do

      else if ( n_field_type == 3 ) then

         fac=or1(n_r)*O_sin_theta(n_theta)
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*bp(n_theta,n_phi)
         end do

      else if ( n_field_type == 4 ) then

         fac=or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vr(n_theta,n_phi)
         end do

      else if ( n_field_type == 5 ) then

         fac=or1(n_r)*orho1(n_r)*O_sin_theta(n_theta)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vt(n_theta,n_phi)
         end do

      else if ( n_field_type == 6 ) then

         fac=or1(n_r)*orho1(n_r)*O_sin_theta(n_theta)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vp(n_theta,n_phi)
         end do

      else if ( n_field_type == 7 ) then

         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=sr(n_theta,n_phi)
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

      else if ( n_field_type == 17 ) then

         fac=-or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            frames(n_o+n_phi)=fac*vr(n_theta,n_phi)*drSr(n_theta,n_phi)
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
   subroutine store_fields_3d(vr,vt,vp,br,bt,bp,sr,drSr,dvrdp,dvpdr,dvtdr, &
              &               dvrdt,cvr,cbr,cbt,n_r,n_store_last,          &
              &               n_field_type)
      !
      !  Purpose of this subroutine is to store movie frames for
      !  surfaces r=const. into array frame(*,*)
      !

      !-- Input variables:
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),drSr(:,:)
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

         fac=or2(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*br(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 2 ) then

         fac_r=or1(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*bt(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 3 ) then

         fac_r=or1(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*bp(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 4 ) then

         fac=or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*vr(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 5 ) then

         fac_r=or1(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               frames(n_phi+n_o)=fac*vt(n_theta_cal,n_phi)
            end do
         end do

      else if ( n_field_type == 6 ) then

         fac_r=or1(n_r)*orho1(n_r)*vScale
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

         fac=-or2(n_r)*orho1(n_r)*vScale
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta=n_theta_cal2ord(n_theta_cal)
               n_o=n_or+(n_theta-1)*n_phi_max
               frames(n_phi+n_o)=fac*vr(n_theta_cal,n_phi)*drSr(n_theta_cal,n_phi)
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
      !  sl(r,theta)=d_theta v(r,theta,m=0)/r
      !

      !-- Input variables:
      integer, intent(in) :: n_r             ! No. of radial grid point

      !-- Output variables:
      real(cp), intent(out) ::  sl(:)           ! Field for field lines

      !-- Local variables:
      integer :: n_theta         ! No. of theta
      integer :: n_theta_nhs     ! Counter for thetas in north HS
      integer :: l,lm            ! Degree, counter for degree/order combinations

      real(cp) :: O_r              ! 1/r
      real(cp) :: O_sint           ! 1/sin(theta)
      complex(cp) :: w_Rloc(lm_max)
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)
      complex(cp) :: Tl_AX(1:l_max+1)

      !-- Calculate radial dependencies:
      O_r=or1(n_r)

      call gather_Flm(w_Rdist(:,n_r),w_Rloc) !@> TODO: one should find a better fix here!

      Tl_AX(1)=zero
      do l=1,l_max
         lm=map_glbl_st%lm2(l,0)
         Tl_AX(l+1)=O_r*w_Rloc(lm)
      end do

      call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:))

      do n_theta=1,n_theta_max,2 ! loop over thetas in Northern HS
         n_theta_nhs=(n_theta+1)/2
         O_sint=osn1(n_theta_nhs)
         sl(n_theta)  =O_sint*tmpp(n_theta)
         sl(n_theta+1)=O_sint*tmpp(n_theta+1)
      end do

   end subroutine get_sl
!----------------------------------------------------------------------------
   subroutine get_fl(fl,n_r,l_ic)
      !
      !  Return field fl whose contourlines are the fields lines
      !  of the axisymmetric poloidal mangetic field.
      !    fl(r,theta)=d_theta b(r,theta,m=0)/r
      !

      !  This routine is called for l_ic=.true. only from coord_r 0 with full
      !  field b_ic in standard lm ordering available.
      !  The case l_ic=.false. is called from all ranks and uses b_Rloc.

      !-- Input variables:
      integer, intent(in) :: n_r             ! No. of radial grid point
      logical, intent(in) :: l_ic            ! =true if inner core field

      !-- Output variables:
      real(cp), intent(out) ::  fl(:)    ! Field for field lines

      !-- Local variables:
      integer :: n_theta         ! No. of theta
      integer :: n_theta_nhs     ! Counter for thetas in north HS
      integer :: l,lm            ! Degree, counter for degree/order combinations

      real(cp) :: r_ratio          ! r/r_ICB
      real(cp) :: O_r              ! 1/r
      real(cp) :: O_sint           ! 1/sin(theta)
      real(cp) :: r_dep(l_max)     ! (r/r_ICB)**l / r_ICB
      complex(cp) :: b_Rloc(lm_max)
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)
      complex(cp) :: Tl_AX(1:l_max+1)

      if ( .not. l_ic ) call gather_Flm(b_Rdist(:,n_r),b_Rloc) !@> TODO: one should find a better fix here!

      if ( l_ic ) then
         r_ratio =r_ic(n_r)/r_ic(1)
         r_dep(1)=r_ratio/r_ic(1)
         do l=2,l_max
            r_dep(l)=r_dep(l-1)*r_ratio
         end do
      else
         O_r=or1(n_r)
      end if

      Tl_AX(1)=zero
      do l=1,l_max
         lm=map_glbl_st%lm2(l,0)
         if ( l_ic ) then ! Inner Core
            if ( l_cond_ic ) then
               Tl_AX(l+1)=r_dep(l)*b_ic(lm,n_r)
            else
               Tl_AX(l+1)=r_dep(l)*bICB(lm)
            end if
         else             ! Outer Core
            Tl_AX(l+1)=O_r*b_Rloc(lm)
         end if
      end do

      call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:))

      do n_theta=1,n_theta_max,2 ! loop over thetas in Northern HS
         n_theta_nhs=(n_theta+1)/2
         O_sint=osn1(n_theta_nhs)
         fl(n_theta)  =O_sint*tmpp(n_theta)
         fl(n_theta+1)=O_sint*tmpp(n_theta+1)
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
      !   r**2 br, r**2 sin(theta) bt, r**2 sin(theta) bp
      !

      !-- Input of variables:
      complex(cp), intent(in) :: bCMB(:)

      !-- Output:
      real(cp), intent(out) :: b_r(:,:) !Radial magnetic field in (phi,theta)-space
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
         l = map_glbl_st%lm2l(lm)
         cs1(lm) = bCMB(lm)*r_dep(l) ! multiplication by l(l+1) in shtns.f90
         cs2(lm)= -bCMB(lm)*real(l,cp)*r_dep(l)
      end do

      zerosc(:)=zero
      call torpol_to_spat_loc(cs1, cs2, zerosc, b_r, b_t, b_p, l_max)

   end subroutine get_B_surface
!----------------------------------------------------------------------------
end module out_movie
