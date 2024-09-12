module movie_data

   use iso_fortran_env, only: output_unit
   use parallel_mod
   use precision_mod
   use truncation, only: n_r_max, n_theta_max, n_phi_max, minc, n_r_ic_max, n_r_tot
   use logic, only: l_store_frame, l_save_out, l_movie, l_movie_oc, l_geosMovie, &
       &            l_movie_ic, l_HTmovie, l_dtBmovie, l_store_frame, l_save_out,&
       &            l_phaseMovie, l_dtphaseMovie
   use radial_data, only: nRstart,nRstop, n_r_icb, n_r_cmb, radial_balance
   use radial_functions, only: r_cmb, r_icb, r, r_ic
   use horizontal_data, only: theta_ord, phi, n_theta_ord2cal
   use output_data, only: n_log_file, log_file, tag
   use charmanip, only: capitalize, delete_string, dble2str
   use useful, only: logWrite, abortRun
   use constants, only: pi, one
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   real(cp), public :: movieDipColat,movieDipLon
   real(cp), public :: movieDipStrength,movieDipStrengthGeo

   !-- Info in movie type and were the frames are stored:
   integer, public, parameter :: n_movies_max=30  ! Max no. of different movies
   integer, public, parameter :: n_movie_fields_max=6 ! Max no. of fields per movie
   real(cp), public ::  movie_const(n_movies_max)
   character(len=80), public :: movie(n_movies_max)  ! Only for input
   character(len=72), public :: movie_file(n_movies_max)

   logical, public :: lStoreMov(n_movies_max),lICField(n_movies_max)
   logical, public :: lGeosField(n_movies_max), lPhaseField(n_movies_max)
   logical :: lAxiField(n_movies_max)
   integer, public :: n_movies
   integer, public :: n_movie_surface(n_movies_max)
   integer, public :: n_movie_const(n_movies_max)
   integer, public :: n_movie_fields(n_movies_max)
   integer, public :: n_movie_fields_ic(n_movies_max)
   integer, public :: n_movie_field_type(n_movie_fields_max,n_movies_max)
   integer, public :: n_movie_field_start(n_movie_fields_max,n_movies_max)
   integer, public :: n_movie_field_stop(n_movie_fields_max,n_movies_max)

   integer, public :: n_movie_file(n_movies_max)

   !-- Work arrays for storing movie frame:
   integer, public :: n_frame_work
   integer, public :: n_MD
   real(cp), public, allocatable :: frames(:)

   public :: initialize_movie_data, finalize_movie_data, &
   &         movie_gather_frames_to_rank0

contains

   subroutine initialize_movie_data
      !
      ! This routine is called during the initialization of the code.
      ! It allows to:
      !
      !    * Estimate the required memory imprint and allocate the arrays
      !      accordingly
      !    * Open the requested movie files
      !

      integer :: n

      movieDipColat      =0.0_cp
      movieDipLon        =0.0_cp
      movieDipStrength   =0.0_cp
      movieDipStrengthGeo=0.0_cp

      if ( .not. l_movie ) then
         l_movie_oc=.false.
         l_movie_ic=.false.
         l_HTmovie =.false.
         l_dtBmovie=.false.
         l_geosMovie=.false.
         l_store_frame=.false.
      else
         call get_movie_type()

         ! Allocate the required memory
         n_MD=maxval(n_movie_field_stop)
         n_frame_work=max(n_MD,1)
         allocate( frames(n_frame_work) )
         frames(:)=0.0_cp
         bytes_allocated = bytes_allocated+n_frame_work*SIZEOF_DEF_REAL

         if ( rank == 0 ) then
            !----- Open movie files on first processor only:
            if ( .not. l_save_out ) then
               do n=1,n_movies
                  open(newunit=n_movie_file(n), file=movie_file(n), &
                  &    status='new', form='unformatted')
               end do
            end if
         end if
      end if

   end subroutine initialize_movie_data
!----------------------------------------------------------------------------
   subroutine finalize_movie_data
      !
      ! Close movie files
      !
      integer :: n

      if ( rank == 0 .and. l_movie ) then
         if ( l_movie ) then
            do n=1,n_movies
               close(n_movie_file(n))
            end do
         end if
      end if

   end subroutine finalize_movie_data
!----------------------------------------------------------------------------
   subroutine get_movie_type
      !
      !  Purpose of this subroutine is to identify the different movie
      !  types from the input string movies(*).
      !  Note that generally blanks are not interpreted and that the
      !  interpretation is not case sensitive.
      !  In general two informations are needed:
      !
      !    1. A word FIELDINFO that identifies the field to be plotted
      !       (e.g. Br for radial magnetic field, see list below)
      !       Possible keywords are (optional text in brackets):
      !
      !         - B r[adial]     : radial magnetic field
      !         - B t[heta]      : theta component
      !         - B p[hi]        : azimuthal component
      !         - B h[orizontal] : the two horizontal components
      !         - B a[ll]        : all three components
      !         - FIELDLINE[S]   : field lines of axisymmetric
      !           or FL            poloidal field for phi=constant
      !         - AX[ISYMMETRIC] B
      !           or AB          : axisymmetric phi component of the
      !           magnetic field for phi=constant
      !         - V r[adial]     : radial velocity field
      !         - V t[heta]      : theta component
      !         - V p[hi]        : azimuthal component
      !         - V h[orizontal] : the two horizontal components
      !         - V a[ll]        : all three components
      !         - STREAMLINE[S]  : field lines of axisymmetric
      !           or SL          : poloidal field for phi=constant
      !         - AX[ISYMMETRIC] V
      !           or AV          : axisymmetric phi component of the
      !           velocity field for phi=constant
      !         - V z            : z component of velocity at equator
      !           and z component of the vorticity at
      !           the equator (closest point to equator)
      !         - Vo z          : z-component of vorticity
      !         - Vo r          : r-component of vorticity
      !         - Vo p          : phi-component of vorticity
      !         - T[emperature]  : sic
      !         - AX[ISYMMETRIC] T
      !           or AT          : axisymmetric T field for phi=constant
      !         - Heat t[ransport]: radial derivative of T
      !         - C[omposition]  : sic
      !         - AX[ISYMMETRIC] C
      !           or AC          : axisymmetric C field for phi=constant
      !         - FL Pro         : axisymmetric field line stretching
      !         - FL Adv         : axisymmetric field line advection
      !         - FL Dif         : axisymmetric field line diffusion
      !         - AB Pro         : axisymmetric (tor.) Bphi production
      !         - AB Dif         : axisymmetric (tor.) Bphi diffusion
      !         - Br Pro         : Br production
      !         - Br Adv         : Br advection
      !         - Br Dif         : Br diffusion
      !         - Jr             : Jr production
      !         - Jr Pro         : Jr production +  omega effects
      !         - Jr Adv         : Jr advection
      !         - Jr Dif         : Jr diffusion
      !         - Bz Pol         : poloidal Bz
      !         - Bz Pol Pro     : poloidal Bz production
      !         - Bz Pol Adv     : poloidal Bz advection
      !         - Bz Pol Dif     : poloidal Bz diffusion
      !         - Jz Tor         : poloidal Jz
      !         - Jz Tor Pro     : poloidal Jz production
      !         - Jz Tor Adv     : poloidal Jz advection
      !         - Jz Tor Dif     : poloidal Jz diffusion
      !         - Bp Tor         : toriodal Bphi
      !         - Bp Tor Pro     : toriodal Bphi production
      !         - Bp Tor Adv     : toriodal Bphi advection
      !         - Bp Tor Dif     : toriodal Bphi diffusion
      !         - HEL[ICITY]     : sic
      !         - AX[ISYMMETRIC HELICITY]  or
      !           AHEL           : axisymmetric helicity
      !         - Bt Tor         : toroidal Btheta
      !         - Pot Tor        : toroidal Potential
      !         - Pol Fieldlines : toroidal Potential
      !         - Br Shear       : azimuthal Shear of Br
      !         - Lorentz[force] : Lorentz force (only phi component)
      !         - Br Inv         : Inverse field apperance at CMB
      !
      !    2. A second information that identifies the coordinate
      !       to be kept constant (surface).
      !       E.g. r=number for surface r=constant with number given
      !       in units of the total core radius or
      !       theta/phi=number with number given in degrees
      !       Four keywords are also possible:
      !
      !         - CMB       : core mantle boundary
      !         - EQ[UATOR] : equatorial plane
      !         - SUR[FACE] : Earth surface (only magnetic field)
      !         - 3[D]      : 3D field throughout the OC [and IC for B]
      !
      !  On output the necessary information is coded into integers
      !  and is used in this form by further subroutines:
      !
      !     * n_movies = total number of movies
      !
      !     * n_movie_surface(n_movie) = defines surface
      !     * n_movie_surface =  1  : r=constant:
      !
      !                     - 2  : theta=constant
      !                     - 3  : phi=constant
      !                     - -1 : r=constant, Earth surface
      !                     - 0  : 3d volume
      !
      !     * n_movie_fields(n_movie) = no. of fields for outer core
      !     * n_movie_fields_ic(n_movie) = no. of fields for inner core
      !     * n_movie_field_type(n_field,n_movie) = defines field
      !     * n_movie_field_type:
      !
      !                      - = 1 : radial magnetic field
      !                      - = 2 : theta comp. of the magnetic field
      !                      - = 3 : azimuthal magnetic field
      !                      - = 4 : radial velocity field
      !                      - = 5 : theta comp. of the velocity field
      !                      - = 6 : azimuthal velocity field
      !                      - = 7 : temperature field
      !                      - = 8 : scalar field for field lines
      !                      - = 9 : axisymm. toroidal mag. field
      !                      - =10 : scalar field for stream lines
      !                      - =11 : axisymm. v_phi
      !                      - =12 : axisymm. T
      !                      - =13 : z-comp. of poloidal Bz
      !                      - =14 : z-comp. of poloidal Jz
      !                      - =15 : z-comp. of velocity
      !                      - =16 : z-comp. of vorticity
      !                      - =17 : convective heat flux T * vr
      !                      - =18 : helicity
      !                      - =19 : axisymmetric helicity
      !                      - =20 : axisymm field-line production
      !                      - =21 : axisymm field-line advection
      !                      - =22 : axisymm field-line diffusion
      !                      - =23 : axisymm Bphi production
      !                      - =24 : axisymm Bphi omega effect
      !                      - =25 : axisymm Bphi advection
      !                      - =26 : axisymm Bphi diffusion
      !                      - =27 : Br production
      !                      - =28 : Br advection
      !                      - =29 : Br diffusion
      !                      - =30 : Jr
      !                      - =31 : Jr production
      !                      - =32 : Jr omega effect
      !                      - =33 : Jr advection
      !                      - =34 : Jr diffusion
      !                      - =35 : poloidal Bz production
      !                      - =36 : poloidal Bz advection
      !                      - =37 : poloidal Bz diffusion
      !                      - =38 : poloidal Jz production
      !                      - =39 : poloidal Jz omega effect
      !                      - =40 : poloidal Jz advection
      !                      - =41 : poloidal Jz diffusion
      !                      - =42 : toroidal Bp
      !                      - =43 : toroidal Bp production
      !                      - =44 : toroidal Bp omega effect
      !                      - =45 : toroidal Bp advection
      !                      - =46 : toroidal Bp diffusion
      !                      - =47 : phi-comp. of vorticity
      !                      - =48 : r-comp. of vorticity
      !                      - =49 : toroidal Bp omega effect
      !                      - =50 : toroidal Bt
      !                      - =51 : toroidal Potential
      !                      - =52 : poloidal Fieldlines in theta=const
      !                      - =53 : Br dr ( vp/(r sin(theta))
      !                      - =54 : phi Lorentz force
      !                      - =61 : AS phi reynolds stress force
      !                      - =62 : AS phi advective stress force
      !                      - =63 : AS phi viscous stress force
      !                      - =64 : AS phi Lorentz force
      !                      - =66 : time derivative of axisym. v phi
      !                      - =67 : relative strength of axisym. v phi
      !                      - =81 : Br inverse appearence at CMB
      !                      - =91 : radial derivative of T or s
      !                      - =92 : axisymmetric radial derivative of T or s
      !                      - =94 : axisymmetric s component of velocity
      !                      - =95 : axisymmetric Reynolds stress correlation us*uphi
      !                      - =96 : axisymmetric z component of velocity
      !                      - =97 : axisymmetric Reynolds stress correlation uz*uphi
      !                      - =98 : axisymmetric square of s component of velocity
      !                      - =99 : axisymmetric square of z component of velocity
      !                      - =100: geostrophic Vs
      !                      - =101: geostrophic Vphi
      !                      - =102: geostrophic z vorticity
      !                      - =103: AS poloidal Br diffusion
      !                      - =104: AS toroidal Bp production
      !                      - =105: AS toroidal Bp dynamo term
      !                      - =106: AS toroidal Bp omega effect
      !                      - =107: AS toroidal Bp diffusion
      !                      - =108: Bs
      !                      - =109: composition field
      !                      - =110: axisymmetric chemical composition
      !                      - =111: axisymmetric phase field
      !                      - =112: phase field
      !                      - =113: radial derivative of chemical composition
      !                      - =114: axisymmetric kinetic energy
      !                      - =115: axisymmetric convective heat flux
      !                      - =116: axisymmetric radial derivative of composition
      !                      - =117: melting radius
      !                      - =118: temperature gradient along the solid/liquid interface
      !
      !     * n_movie_field_start(n_field,n_movie) = defines where first
      !       element of a field is stored in ``frames(*)``
      !     * n_movie_field_stop(n_field,n_movie) = defines where last
      !       element of a field is stored in ``frames(*)``
      !
      !     * The subroutine also defines appropriate file names for
      !       the movie files. These generally have the form TYPE_mov.TAG
      !

      !--- Local variables:
      character(len=80) :: word
      character(len=:), allocatable :: message, string, stringC
      character(len=:), allocatable :: file_name, typeStr
      real(cp) :: r_movie,theta_movie,phi_movie
      real(cp) :: phi_max, rad, const
      integer :: length
      integer :: i, n, n_ic, ns
      integer :: n_type, n_surface, n_const
      integer :: n_fields,n_fields_oc,n_fields_ic
      integer :: n_field_size, n_field_size_ic, n_field_start
      integer :: n_field_type(n_movie_fields_max)
      logical :: lStore, lIC, lAxi, lGeos, lPhase

      !--- Initialize first storage index:
      n_field_start=1

      !--- Converts from radiant to degree:
      rad=180.0_cp/pi
      n_type=0

      !--- Loop over max possible no of movies:
      l_movie_oc    =.false.
      l_movie_ic    =.false.
      l_HTmovie     =.false.
      l_dtBmovie    =.false.
      l_geosMovie   =.false.
      l_phaseMovie  =.false.
      l_dtphaseMovie=.false.
      l_store_frame =.false.
      n_field_type(:)=0

      do i=1,n_movies_max

         lStore=.true.
         lIC   =.false.
         lPhase=.false.
         lGeos =.false.
         lAxi  =.false.

         string=movie(i)

         if ( len_trim(string)  ==  0 ) cycle !blank string

         !--- Delete blanks, they are not interpreted
         call delete_string(string,' ',length)

         !--- Convert to capitals:
         call capitalize(string)

         !--- Identify movie type (fields to be stored):

         if ( index(string,'BR') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               typeStr=' radial magnetic field production '
               file_name='BrPro_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=27
               l_dtBmovie=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' radial magnetic field advection '
               file_name='BrAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=28
               l_dtBmovie=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' radial magnetic field diffusion '
               file_name='BrDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=29
               l_dtBmovie=.true.
            else if ( index(string,'INV') /= 0 ) then
               typeStr=' radial magnetic field diffusion '
               file_name='BrInv_'
               lIC=.false.
               n_fields=1
               n_field_type(1)=81
               l_dtBmovie=.true.
            else if ( index(string,'SHEAR') /= 0 ) then
               typeStr=' azimuthal shear of radial magnetic field'
               file_name='BrShear'
               lIC=.false.
               n_fields=1
               n_field_type(1)=53
               l_dtBmovie=.true.
            else
               typeStr=' radial magnetic field '
               n_fields=1
               file_name='Br_'
               lIC=.true.
               n_field_type(1)=1
            end if
         else if ( index(string,'BT') /= 0 ) then
            if ( index(string,'TOR') /= 0 ) then
               typeStr=' toroidal B theta'
               file_name='BtTor_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=50
               l_dtBmovie=.true.
            else
               typeStr=' theta comp. of magnetic field '
               n_fields=1
               file_name='Bt_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=2
            end if
         else if ( index(string,'BP' ) /= 0 .and. index(string,'TOR') == 0 .and. &
         &    index(string,'AX' ) == 0 .and. index(string,'AB' ) == 0 ) then
            typeStr=' azimuthal magnetic field '
            file_name='Bp_'
            lIC=.true.
            n_fields=1
            n_field_type(1)=3
            lAxi=.true.
         else if ( index(string,'BH') /= 0 .or. index(string,'BM') /= 0 ) then
            n_type=4 ! Horizontal field
            file_name='Bh_'
            lIC=.true.
         else if ( index(string,'BS') /= 0) then
            typeStr='cyl radial magnetic field'
            file_name='Bs_'
            lIC=.true.
            n_fields=1
            n_field_type(1)=108
         else if ( index(string,'BALL') /= 0 ) then
            typeStr=' all magnetic field components '
            n_fields=3
            file_name='B_'
            lIC=.true.
            n_fields=3
            n_field_type(1)=1
            n_field_type(2)=2
            n_field_type(3)=3
         else if ( index(string,'FIELDLINE') /= 0 .or. index(string,'FL') /= 0 ) then
            if ( index(string,'POL') /= 0 ) then
               typeStr=' pol. fieldlines'
               file_name='FLPol_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=52
               l_dtBmovie=.true.
               lAxi=.true.
            else if ( index(string,'PRO') /= 0 ) then
               typeStr=' axisymm. fieldline production '
               file_name='FLPro_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=20
               l_dtBmovie=.true.
               lAxi=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' axisymm. fieldline advection '
               file_name='FLAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=21
               l_dtBmovie=.true.
               lAxi=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' axisymm. fieldline diffusion '
               file_name='FLDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=22
               l_dtBmovie=.true.
               lAxi=.true.
            else
               typeStr=' axisymm. fieldlines '
               file_name='FL_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=8
               lAxi=.true.
            end if
         else if ( index (string,'LORENTZ') /= 0 .or. index(string,'LF') /= 0 ) then
            typeStr=' phi Lorentz force '
            file_name='LFp_'
            lIC=.false.
            n_fields=1
            n_field_type(1)=54
         else if ( ( index(string,'AX') /= 0 .and.  &
         &    index(string,'BP') /= 0 ) .or. index(string,'AB') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               typeStr=' axisymm. B phi production '
               file_name='ABPro_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=23
               l_dtBmovie=.true.
               lAxi=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' axisymm. B phi advection '
               file_name='ABAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=25
               l_dtBmovie=.true.
               lAxi=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' axisymm. B phi diffusion '
               file_name='ABDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=26
               l_dtBmovie=.true.
               lAxi=.true.
            else
               typeStr=' axisymm. B phi '
               file_name='AB_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=9
               lAxi=.true.
            end if
         else if ( index(string,'JR') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               typeStr=' radial current production '
               file_name='JrPro_'
               lIC=.true.
               n_fields=2
               n_field_type(1)=31
               n_field_type(2)=32
               l_dtBmovie=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' radial current advection '
               file_name='JrAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=33
               l_dtBmovie=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' radial current diffusion '
               file_name='JrDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=34
               l_dtBmovie=.true.
            else
               typeStr=' radial current '
               file_name='Jr_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=30
               l_dtBmovie=.true.
            end if
         else if ( index(string,'BZ') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               typeStr=' poloidal Bz production '
               file_name='BzPolPro_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=35
               l_dtBmovie=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' poloidal Bz advection '
               file_name='BzPolAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=36
               l_dtBmovie=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' poloidal Bz diffusion '
               file_name='BzPolDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=37
               l_dtBmovie=.true.
            else
               typeStr=' poloidal Bz '
               file_name='BzPol_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=13
               l_dtBmovie=.true.
            end if
         else if ( index(string,'BP') /= 0 .and. index(string,'TOR') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               typeStr=' toroidal Bphi production '
               file_name='BpTorPro_'
               lIC=.true.
               n_fields=2
               n_field_type(1)=43
               n_field_type(2)=49
               l_dtBmovie=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' toroidal Bphi advection '
               file_name='BpTorAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=45
               l_dtBmovie=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' toroidal Bphi diffusion '
               file_name='BpTorDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=46
               l_dtBmovie=.true.
            else
               typeStr=' toroidal Bphi '
               file_name='BpTor_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=42
               l_dtBmovie=.true.
            end if
         else if ( index(string,'JZ') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               typeStr=' poloidal Jz production '
               file_name='JzTorPro_'
               lIC=.true.
               n_fields=2
               n_field_type(1)=38
               n_field_type(2)=39
               l_dtBmovie=.true.
            else if ( index(string,'ADV') /= 0 ) then
               typeStr=' poloidal Jz advection '
               file_name='JzTorAdv_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=40
               l_dtBmovie=.true.
            else if ( index(string,'DIF') /= 0 ) then
               typeStr=' poloidal Jz diffusion '
               file_name='JzTorDif_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=41
               l_dtBmovie=.true.
            else
               typeStr=' poloidal Jz '
               file_name='JzTor_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=14
               l_dtBmovie=.true.
            end if
         ! Possible conflict with HEATTCONV r= (which contains VR)
         else if ( index(string,'VR') /= 0 .and. index(string,'CONV') == 0) then
            typeStr=' radial velocity field '
            file_name='Vr_'
            n_fields=1
            n_field_type(1)=4
         else if ( index(string,'VT') /= 0 ) then
            typeStr=' theta comp. of velocity field '
            file_name='Vt_'
            n_fields=1
            n_field_type(1)=5
         else if ( index(string,'VP') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. azimuthal velocity field '
               file_name='AV_'
               n_fields=1
               n_field_type(1)=11
               lAxi=.true.
            else if ( index(string,'GEOS') /= 0 ) then
               typeStr=' geos phi-component of velocity '
               file_name='geosVPHI_'
               n_fields=1
               lStore=.false.
               l_geosMovie=.true.
               n_field_type(1)=101
               lGeos=.true.
            else
               typeStr=' phi comp. of velocity field '
               file_name='Vp_'
               n_fields=1
               n_field_type(1)=6
            end if
         else if ( index(string,'RMELT') /= 0 ) then
            typeStr=' melting radius'
            file_name='rmelt_'
            n_fields=1
            lStore=.false.
            l_phaseMovie=.true.
            n_field_type(1)=117
            lPhase=.true.
         else if ( index(string,'DTRM') /= 0 .or. index(string,'DTPHASE') /= 0 ) then
            typeStr=' temperature gradient at melting radius'
            file_name='dt_rmelt_'
            n_fields=1
            lStore=.false.
            l_dtphaseMovie=.true.
            n_field_type(1)=118
            lPhase=.true.
         else if ( index(string,'VH') /= 0 .or. index(string,'VM') /= 0 ) then
            n_type=14
            file_name='Vh_'
         else if ( index(string,'VALL') /= 0 ) then
            typeStr=' all velocity components '
            file_name='V_'
            n_fields=3
            n_field_type(1)=4
            n_field_type(2)=5
            n_field_type(3)=6
         else if ( index(string,'STREAMLINE') /= 0 .or. index(string,'SL') /= 0 ) then
            typeStr=' axisym. meridional velocity streamlines '
            file_name='SL_'
            n_fields=1
            n_field_type(1)=10
            lAxi=.true.
         else if ( index(string,'VOR') /= 0 ) then
            if ( index(string,'Z') /= 0 ) then
               if ( index(string,'GEOS') /= 0 ) then
                  typeStr=' geos z-component of vorticity '
                  file_name='geosVorZ_'
                  n_fields=1
                  lStore=.false.
                  l_geosMovie=.true.
                  n_field_type(1)=102
                  lGeos=.true.
               else
                  typeStr=' z-component of vorticity '
                  file_name='VorZ_'
                  n_fields=1
                  n_field_type(1)=16
               end if
            else if ( index(string,'P') /= 0 ) then
               typeStr=' phi-component of vorticity '
               file_name='VorP_'
               n_fields=1
               n_field_type(1)=47
            else if ( index(string,'R') /= 0 ) then
               typeStr=' r-component of vorticity '
               file_name='VorR_'
               n_fields=1
               n_field_type(1)=48
            end if
         else if ( index(string,'VS2') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. vs**2 '
               file_name='AVS2_'
               n_fields=1
               n_field_type(1)=98
               lAxi=.true.
            end if
         else if ( index(string,'VS') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. s-component of velocity '
               file_name='AVS_'
               n_fields=1
               n_field_type(1)=94
               lAxi=.true.
            else if ( index(string,'GEOS') /= 0 ) then
               typeStr=' geos s-component of velocity '
               file_name='geosVS_'
               n_fields=1
               lStore=.false.
               l_geosMovie=.true.
               n_field_type(1)=100
               lGeos=.true.
            end if
         else if ( index(string,'REYS') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. vs*vphi '
               file_name='AReyS_'
               n_fields=1
               n_field_type(1)=95
               lAxi=.true.
            end if
         else if ( index(string,'REYZ') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. vz*vphi '
               file_name='AReyZ_'
               n_fields=1
               n_field_type(1)=97
               lAxi=.true.
            end if
         else if ( index(string,'VZ2') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. vz**2 '
               file_name='AVZ2_'
               n_fields=1
               n_field_type(1)=99
               lAxi=.true.
            end if
         else if ( index(string,'VZ') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr=' axisym. vz '
               file_name='AVZ_'
               n_fields=1
               n_field_type(1)=96
               lAxi=.true.
            else
               typeStr=' z-component of velocity '
               file_name='VZ_'
               n_fields=1
               n_field_type(1)=15
            end if
         else if ( ( index(string,'AX') /= 0 .and.  &
         &    index(string,'HEL' ) /= 0 ) .or. index(string,'AHEL') /= 0 ) then
            typeStr=' axisymmetric helicity '
            file_name='AH_'
            n_fields=1
            n_field_type(1)=19
            lAxi=.true.
         else if ( index(string,'HEL') /= 0 ) then
            typeStr=' helicity '
            file_name='HE_'
            n_fields=1
            n_field_type(1)=18
         else if ( index(string,'AX') /= 0 .and. &
         &    ( index(string,'TEM') /= 0 .or. index(string,'ENT') /= 0 ) ) then
            typeStr=' axisymmetric temp. '
            file_name='AT'
            n_fields=1
            n_field_type(1)=12
            lAxi=.true.
         else if ( index(string,'AX') /= 0 .and. &
         &    ( index(string,'COMP') /= 0 .or. index(string,'XI') /= 0 ) ) then
            typeStr=' axisymmetric comp. '
            file_name='AC'
            n_fields=1
            n_field_type(1)=110
            lAxi=.true.
         else if ( index(string,'AX') /= 0 .and. ( index(string,'PHASE') /= 0 ) ) then
            typeStr=' axisymmetric phase '
            file_name='APHI'
            n_fields=1
            n_field_type(1)=111
            lAxi=.true.
         else if ( index(string,'ENTROPY') /= 0 .or. index(string,'TEM') /= 0 ) then
            ns=index(string,'S')
            if ( ns > 0 ) then
               if ( string(ns:ns+2) == 'SUR' ) then
                  call abortRun('! No surface T field available !')
               end if
            end if
            typeStr=' temperature field '
            file_name='T_'
            n_fields=1
            n_field_type(1)=7
         else if ( index(string,'COMP') /= 0 .or. index(string,'XI') /= 0 ) then
            ns=index(string,'S')
            if ( ns > 0 ) then
               if ( string(ns:ns+2) == 'SUR' ) then
                  call abortRun('! No surface C field available !')
               end if
            end if
            typeStr=' Composition field '
            file_name='XI_'
            n_fields=1
            n_field_type(1)=109
         else if ( index(string,'PHASE') /= 0 ) then
            ns=index(string,'S')
            if ( ns > 0 ) then
               if ( string(ns:ns+2) == 'SUR' ) then
                  call abortRun('! No surface C field available !')
               end if
            end if
            typeStr=' phase field '
            file_name='PHI_'
            n_fields=1
            n_field_type(1)=112
         else if ( index(string,'HEATTCONV') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr='axisymmetric convective heat transport '
               file_name='AHT_'
               n_fields=1
               n_field_type(1)=115
               lAxi=.true.
            else
               ns=index(string,'S')
               if ( ns > 0 ) then
                  if ( string(ns:ns+2) == 'SUR' ) then
                     call abortRun('! No surface convective heat flux !')
                  end if
               end if
               typeStr=' radial convective heat transport '
               file_name='HT_'
               n_fields=1
               n_field_type(1)=17
            end if
         else if ( index(string,'HEATXCONV') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr='axisymmetric transport of composition '
               file_name='AHXi_'
               n_fields=1
               n_field_type(1)=116
               lAxi=.true.
            else
               ns=index(string,'S')
               if ( ns > 0 ) then
                  if ( string(ns:ns+2) == 'SUR' ) then
                     call abortRun('! No surface flux of composition !')
                  end if
               end if
               typeStr=' radial transport of composition '
               file_name='HXi_'
               n_fields=1
               n_field_type(1)=113
            end if
         else if ( index(string,'HEATF') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               typeStr='axisymmetric dSdr '
               file_name='AHF_'
               l_HTmovie=.true.
               n_fields=1
               n_field_type(1)=92
               lAxi=.true.
            else
               ns=index(string,'S')
               if ( ns > 0 ) then
                  if ( string(ns:ns+2) == 'SUR' ) then
                     call abortRun('! No surface T field available !')
                  end if
               end if
               typeStr='dSdr '
               file_name='HF_'
               l_HTmovie=.true.
               n_fields=1
               n_field_type(1)=91
            end if
         else if ( index(string,'AX') /= 0 .and. index(string,'EKIN') /= 0   ) then
            typeStr='axisymmetric kinetic energy '
            file_name='AEKIN_'
            n_fields=1
            n_field_type(1)=114
            lAxi=.true.
         else if ( index(string,'POT') /= 0 .and. index(string,'TOR') /= 0 ) then
            typeStr=' pot tor '
            file_name='PotTor_'
            lIC=.true.
            n_fields=1
            n_field_type(1)=51
            l_dtBmovie=.true.
         else
            message = 'Couldnt interpret movie field from string:'//string
            call abortRun(message)
         end if


         !--- Identify surface type:
         if ( n_field_type(1) == 81 ) then
            n_surface=1 !
            n_const=1   !
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=r_cmb
         else if ( lGeos ) then ! Geos average
            n_surface=-2 ! constant theta
            n_const=1   !
            n_field_size=n_phi_max*n_r_max
            n_field_size_ic=0
            const=r_cmb
         else if ( lPhase ) then ! Melting radius or temp gradient at rm
            n_surface=1  ! R=const.
            n_const=1   !
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=0
            const=r_cmb
         else if (   index(string,'AX') /= 0 .or. lAxi ) then
            !--- Axisymmetric stuff:
            n_surface=3  ! PHI=const.
            n_const=1
            n_field_size=n_r_max*n_theta_max
            n_field_size_ic=n_r_ic_max*n_theta_max
            const=0.0_cp
         else if ( index(string,'3D') /= 0 ) then
            n_surface=0  ! 3d
            n_const=0    ! Not needed
            file_name=file_name//'3D_'
            n_field_size=n_r_max*n_theta_max*n_phi_max
            n_field_size_ic=n_r_ic_max*n_theta_max*n_phi_max
         else if ( index(string,'CMB') /= 0 ) then
            n_surface=1 ! R=const. at CMB
            n_const=1
            file_name=file_name//'CMB_'
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=r_cmb
         else if ( index(string,'ICB') /= 0 ) then
            n_surface=1 ! R=const. at ICB
            n_const=n_r_max
            file_name=file_name//'ICB_'
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=r_icb
         else if ( index(string,'SUR') /= 0 ) then
            n_surface=-1
            n_const=1
            file_name=file_name//'SUR_'
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=one
         else if ( index(string,'R=') /= 0 .or. &
         &    index(string,'RAD=') /= 0 .or. index(string,'RADIUS=') /= 0 ) then
            n_surface=1  ! R=const.
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=0!n_field_size

            !------ Get Radius in fractions of outer core radius:
            if ( index(string,'R=') /= 0 ) then
               word=string(index(string,'R=')+2:length)
            else if ( index(string,'RAD=') /= 0 ) then
               word=string(index(string,'RAD=')+4:length)
            else if ( index(string,'RADIUS=') /= 0 ) then
               word=string(index(string,'RADIUS=')+7:length)
            end if
            read(word,*) r_movie


            !------ Choose closest radial grid point:
            if ( r_movie == 0 ) then
               n_const=n_r_icb
               const  =r_icb
            else if ( r_movie == 1 ) then
               n_const=n_r_cmb
               const  =r_cmb
            else if ( r_movie > 0 ) then
               r_movie=r_icb+r_movie
               n_const = minloc(abs(r - r_movie),1)
               const = r(n_const)
            else
               !------ Negative numbers signify inner core values in
               !       fractions of r_icb:
               if ( lIC ) then
                  n_field_size_ic=n_field_size
                  r_movie=-r_movie*r_icb
                  n_const = minloc(abs(r_ic - r_movie),1)
                  const = r_ic(n_const)
               end if
            end if

            call dble2str(r_movie,word)
            stringC='R='//trim(word)//'_'
            file_name=file_name//stringC

         else if ( index(string,'EQ') /= 0 ) then

            n_surface=2    ! Equator
            n_const=n_theta_max/2
            file_name=file_name//'EQU_'
            n_field_size=n_r_max*n_phi_max
            n_field_size_ic=n_r_ic_max*n_phi_max
            const=rad*theta_ord(n_const)
            n_const=n_theta_ord2cal(n_const) ! Scrambling if needed

         else if ( index(string,'T=') /= 0 .or. index(string,'THETA=') /= 0 ) then

            n_surface=2    ! Theta=const.
            n_field_size=n_r_max*n_phi_max
            n_field_size_ic=n_r_ic_max*n_phi_max

            !------ Get desired colatitude value:
            if ( index(string,'T=') /= 0 ) then
               word=string(index(string,'T=')+2:length)
            else if ( index(string,'THETA=') /= 0 ) then
               word=string(index(string,'THETA=')+6:length)
            end if
            read(word,*) theta_movie
            theta_movie=abs(theta_movie)
            theta_movie=theta_movie/rad

            !------ Choose closest colatitude grid point:
            n_const=minloc(abs(theta_ord - theta_movie),1)
            const=rad*theta_ord(n_const)
            n_const=n_theta_ord2cal(n_const)

            call dble2str(theta_movie,word)
            stringC='T='//trim(word)//'_'
            file_name=file_name//stringC

         else if ( index(string,'MER' ) /= 0 .or.  &
         &    index(string,'P='  ) /= 0  .or. index(string,'PHI=') /= 0 ) then

            n_surface=3  ! PHI=const.
            n_field_size=2*n_r_max*n_theta_max
            n_field_size_ic=2*n_r_ic_max*n_theta_max

            !------ Get desired longitude value:
            if ( index(string,'P=') /=0 ) then
               word=string(index(string,'P=')+2:length)
            else if ( index(string,'PHI=') /=0 ) then
               word=string(index(string,'PHI=')+4:length)
            end if
            read(word,*) phi_movie
            if ( phi_movie < 0.0_cp ) phi_movie=360.0_cp-phi_movie
            phi_max=360.0_cp/minc
            if ( minc > 1 ) then
               do n=minc-1,1,-1
                  if ( phi_movie > n*phi_max ) then
                     phi_movie=phi_movie-n*phi_max
                     exit
                  end if
               end do
            end if
            phi_movie=phi_movie/rad

            !------ Choose closest longitude grid point:
            n_const=minloc(abs(phi - phi_movie),1)
            const=rad*phi(n_const)

            call dble2str(phi_movie,word)
            stringC='P='//trim(word)//'_'
            file_name=file_name//stringC

         else
            message = 'Couldnt interpret movie surface from string:'//string
            call abortRun(message)
         end if

         !--- Now store the necessary information:
         !------ Increase number of movies:
         n_movies=n_movies+1
         lStoreMov(n_movies)=lStore
         lICField(n_movies)=lIC
         lGeosField(n_movies)=lGeos
         lAxiField(n_movies)=lAxi
         lPhaseField(n_movies)=lPhase

         !------ Translate horizontal movies:
         if ( n_type == 4 ) then
            if ( n_surface == 1 ) then
               typeStr=' theta and phi B components '
               n_fields=2
               n_field_type(1)=2
               n_field_type(2)=3
            else if ( n_surface == 2 ) then
               typeStr=' radial and phi B components '
               n_fields=2
               n_field_type(1)=1
               n_field_type(2)=3
            else if ( n_surface == 3 ) then
               typeStr=' radial and theta B components '
               n_fields=2
               n_field_type(1)=1
               n_field_type(2)=2
            end if
         else if ( n_type == 14 ) then
            if ( n_surface == 1 ) then
               typeStr=' theta and phi V components '
               n_fields=2
               n_field_type(1)=5
               n_field_type(2)=6
            else if ( n_surface == 2 ) then
               typeStr=' radial and phi V components '
               n_fields=2
               n_field_type(1)=4
               n_field_type(2)=6
            else if ( n_surface == 3 ) then
               typeStr=' radial and theta V components '
               n_fields=2
               n_field_type(1)=4
               n_field_type(2)=5
            end if
         end if

         !--- Inner core and/or outer core field?
         n_fields_oc=n_fields
         n_fields_ic=0
         if ( lIC ) then
            if ( n_surface == 1 .and. n_const < 0 ) then
               n_fields_oc=0
               n_fields_ic=n_fields
            else if ( n_surface == 0 .or. n_surface == 2 .or. n_surface == 3 ) then
               n_fields_ic=n_fields
            end if
         end if
         if ( n_fields_oc > 0 ) l_movie_oc= .true.
         if ( n_fields_ic > 0 ) l_movie_ic= .true.

         !------ Store name of movie file:
         file_name=file_name//'mov.'//tag
         call delete_string(file_name,' ',length)
         movie_file(n_movies)=file_name(1:length)

         !------ Store information about movie surface:
         n_movie_surface(n_movies)=n_surface
         n_movie_const(n_movies)  =n_const
         movie_const(n_movies)    =const

         !------ Store number of fields for this movie:
         n_movie_fields(n_movies)   =n_fields_oc
         n_movie_fields_ic(n_movies)=n_fields_ic

         !------ Store size of field and where it should be stored in
         !       the work array frames(*):
         do n=1,n_fields
            if ( n_fields_oc > 0 ) then
               n_movie_field_type(n,n_movies) =n_field_type(n)
               if ( lStore ) then
                  n_movie_field_start(n,n_movies)=n_field_start
                  n_field_start=n_field_start+n_field_size
                  n_movie_field_stop(n,n_movies) =n_field_start-1
                  l_store_frame=.true.
               else
                  n_movie_field_start(n,n_movies)=-1
                  n_movie_field_stop(n,n_movies) =-1
               end if
            end if
            if ( n_fields_ic > 0 ) then
               n_ic=n_fields_oc+n
               n_movie_field_type(n_ic,n_movies)=n_field_type(n)
               if ( lStore ) then
                  n_movie_field_start(n_ic,n_movies)=n_field_start
                  n_field_start=n_field_start+n_field_size_ic
                  n_movie_field_stop(n_ic,n_movies)=n_field_start-1
               else
                  n_movie_field_start(n_ic,n_movies)=-1
                  n_movie_field_stop(n_ic,n_movies) =-1
               end if
            end if
         end do

         !--- Write info about output files into log-file:
         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(newunit=n_log_file, file=log_file, status='unknown', &
               &    position='append')
            end if
            write(n_log_file,'(/,'' ! PRODUCING MOVIE-FILE :'',a64)') &
            &     movie_file(n_movies)
            write(n_log_file,*) '!    Containing ',typeStr

            if ( n_surface == -1 ) then
               write(n_log_file,*) '!    at the surface !'
            else if ( n_surface == -2 ) then
               write(n_log_file,*) '!    geos movie     !'
            else if ( n_surface == 0 ) then
               write(n_log_file,*) '!    in 3d !'
            else if ( n_surface == 1 .and. n_const == 1 ) then
               write(n_log_file,*) '!    at the core-mantle boundary !'
            else if ( n_surface == 1 .and. n_const > 1 ) then
               write(n_log_file,'('' !    at r='',f10.6)') r(n_const)
            else if ( n_surface == 1 .and. n_const < 0 ) then
               write(n_log_file,'('' !    at r='',f10.6)') r_ic(n_const)
            else if ( n_surface == 2 ) then
               write(n_log_file,'('' !    at theta='',f12.6)') const
            else if ( n_surface == 3 ) then
               write(n_log_file,'('' !    at phi='',f12.6)') rad*phi(n_const)
            end if
            if ( n_fields_ic > 0 ) &
                 write(n_log_file,'('' !    including inner core magnetic field.'')')

            if ( l_save_out ) close(n_log_file)
         end if

      end do     ! loop over all movies

   end subroutine get_movie_type
!----------------------------------------------------------------------------
   subroutine movie_gather_frames_to_rank0()
      !
      ! MPI communicators for movie files
      !

#ifdef WITH_MPI
      integer :: n_fields,n_surface,n_movie,n_const
      integer :: n_start,n_stop,n_field,n_field_type
      integer :: myTag, status(MPI_STATUS_SIZE)
      integer :: local_start,local_end,irank,sendcount
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(cp), allocatable :: field_frames_global(:)
      integer :: max_field_length,field_length

      max_field_length=0
      do n_movie=1,n_movies
         n_fields =n_movie_fields(n_movie)
         do n_field=1,n_fields
            n_start = n_movie_field_start(n_field,n_movie)
            n_stop  = n_movie_field_stop(n_field,n_movie)
            field_length=n_stop-n_start+1
            if (field_length > max_field_length) max_field_length=field_length
         end do
      end do
      if ( rank == 0 ) then
         allocate(field_frames_global(max_field_length))
      else
         ! This is only needed for debug runs with boundary check.
         allocate(field_frames_global(1))
      end if

      ! loop over all movies
      do n_movie=1,n_movies
         n_fields =n_movie_fields(n_movie)
         n_surface=n_movie_surface(n_movie)
         n_const  =n_movie_const(n_movie)

         select case(n_surface)

            case(-1) ! Earth Surface
               ! theta-phi surface for n_r=1 (CMB)
               ! frames is already existent on rank 0 with all
               ! needed values
               ! do nothing, pass to the next movie
               cycle

            case(0) ! 3d
               ! 3d, all grid points written to frames
               ! but only n_r=nRstart:nRstop on one rank,
               ! gather needed
               do n_field=1,n_fields
                  n_start = n_movie_field_start(n_field,n_movie)
                  n_stop  = n_movie_field_stop(n_field,n_movie)
                  field_length = n_stop-n_start+1

                  local_start=n_start+(nRstart-1)*n_phi_max*n_theta_max
                  local_end  =local_start+nR_per_rank*n_phi_max*n_theta_max-1
                  if (local_end > n_stop) then
                     call abortRun('local_end exceeds n_stop')
                  end if
                  do irank=0,n_procs-1
                     recvcounts(irank) = radial_balance(irank)%n_per_rank* &
                     &                   n_phi_max*n_theta_max
                  end do
                  displs(0)=0
                  do irank=1,n_procs-1
                     displs(irank)=displs(irank-1)+recvcounts(irank-1)
                  end do
                  sendcount=local_end-local_start+1

                  call MPI_Gatherv(frames(local_start),sendcount,MPI_DEF_REAL, &
                       &           field_frames_global,recvcounts,displs,      &
                       &           MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
                  if ( rank == 0 ) then
                     frames(n_start:n_stop)=field_frames_global(1:field_length)
                  end if
               end do

            case(1) ! Surface r=constant
               ! frames is set only for one rank, where n_r=n_const
               ! send to rank 0
               do n_field=1,n_fields
                  n_start=n_movie_field_start(n_field,n_movie)
                  n_stop =n_movie_field_stop(n_field,n_movie)
                  field_length = n_stop-n_start+1
                  myTag=7654+n_movie
                  if ( rank == 0 ) then
                     if ( (nRstart <= n_const) .and. (n_const <= nRstop) ) then
                        ! relevant frames already set on rank 0
                        ! do nothing
                     else
                        call MPI_Recv(frames(n_start),field_length,MPI_DEF_REAL,     &
                             &        MPI_ANY_SOURCE,mytag,MPI_COMM_WORLD,status,ierr)
                     end if
                  else
                     if ( (nRstart <= n_const) .and. (n_const <= nRstop) ) then
                        ! relevant frames are all on this rank  /= 0
                        ! send to rank 0
                        call MPI_Send(frames(n_start),field_length,MPI_DEF_REAL,&
                             &        0,mytag,MPI_COMM_WORLD,ierr)
                     end if
                  end if
               end do

            case(2) ! Surface theta=constant
               do n_field=1,n_fields
                  n_start = n_movie_field_start(n_field,n_movie)
                  n_stop  = n_movie_field_stop(n_field,n_movie)
                  field_length = n_stop-n_start+1

                  local_start=n_start+(nRstart-1)*n_phi_max
                  local_end  =local_start+nR_per_rank*n_phi_max-1
                  if ( local_end > n_stop ) then
                     call abortRun('local_end exceeds n_stop')
                  end if
                  do irank=0,n_procs-1
                     recvcounts(irank)=radial_balance(irank)%n_per_rank*n_phi_max
                  end do
                  displs(0)=0
                  do irank=1,n_procs-1
                     displs(irank)=displs(irank-1)+recvcounts(irank-1)
                  end do
                  sendcount=local_end-local_start+1

                  call MPI_Gatherv(frames(local_start),sendcount,MPI_DEF_REAL, &
                       &           field_frames_global,recvcounts,displs,      &
                       &           MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
                  if ( rank == 0 ) then
                     frames(n_start:n_stop)=field_frames_global(1:field_length)
                  end if
               end do  ! Do loop over field for one movie

            case(3)  ! Surface phi=const.
               ! all ranks have a part of the frames array for each movie
               ! we need to gather

               do n_field=1,n_fields
                  n_start      = n_movie_field_start(n_field,n_movie)
                  n_stop       = n_movie_field_stop(n_field,n_movie)
                  n_field_type = n_movie_field_type(n_field,n_movie)
                  field_length = n_stop-n_start+1

                  local_start=n_start+(nRstart-1)*n_theta_max
                  local_end  =local_start+nR_per_rank*n_theta_max-1
                  do irank=0,n_procs-1
                     recvcounts(irank)=radial_balance(irank)%n_per_rank*n_theta_max
                  end do
                  if ( local_end > n_stop ) then
                     call abortRun('local_end exceeds n_stop')
                  end if
                  displs(0)=0
                  do irank=1,n_procs-1
                     displs(irank)=displs(irank-1)+recvcounts(irank-1)
                  end do
                  sendcount=local_end-local_start+1

                  !-- Either only the axisymmetric or both slices
                  if ( lAxiField(n_movie) ) then
                     call MPI_Gatherv(frames(local_start),sendcount,MPI_DEF_REAL, &
                          &           field_frames_global,recvcounts,displs,      &
                          &           MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
                     if ( rank == 0 ) then
                        frames(n_start:n_stop)=field_frames_global(1:field_length)
                     end if

                  else ! Two phi slices

                     n_stop=n_start+field_length/2-1
                     call MPI_Gatherv(frames(local_start),sendcount,MPI_DEF_REAL, &
                          &           field_frames_global,recvcounts,displs,      &
                          &           MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
                     if ( rank == 0 ) then
                        frames(n_start:n_stop)=field_frames_global(1:field_length/2)
                     end if
                     n_start=n_stop+1
                     n_stop =n_start+field_length/2-1
                     local_start = local_start+field_length/2
                     local_end = local_end+field_length/2
                     call MPI_Gatherv(frames(local_start),sendcount,MPI_DEF_REAL, &
                          &           field_frames_global,recvcounts,displs,      &
                          &           MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
                     if ( rank == 0 ) then
                        frames(n_start:n_stop)=field_frames_global(1:field_length/2)
                     end if

                  end if
               end do  ! Do loop over field for one movie

         end select
      end do

      if ( rank == 0 ) deallocate(field_frames_global)
#endif

   end subroutine movie_gather_frames_to_rank0
!----------------------------------------------------------------------------
end module movie_data
