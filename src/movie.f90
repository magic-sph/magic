module movie_data

   use parallel_mod
   use precision_mod
   use truncation, only: n_r_max, n_theta_max, n_phi_max,      &
       &                 ldtBMem, minc, n_r_ic_max, lMovieMem, &
       &                 n_r_tot
   use logic, only:  l_store_frame, l_save_out, l_movie, &
       &             l_movie_oc, l_movie_ic, l_HTmovie,  &
       &             l_dtBmovie, l_store_frame, l_save_out
   use radial_data, only: nRstart,nRstop, n_r_icb, n_r_cmb, radial_balance
   use radial_functions, only: r_cmb, r_icb, r, r_ic
   use horizontal_data, only: theta, phi
   use output_data, only: n_log_file, log_file, tag
   use charmanip, only: capitalize,delete_string, str2dble,length_to_blank
   use useful, only: logWrite, abortRun
   use constants, only: pi, one
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   real(cp), public :: movieDipColat,movieDipLon
   real(cp), public :: movieDipStrength,movieDipStrengthGeo
   real(cp), public :: t_movieS(10000)

   !-- Info in movie type and were the frames are stored:
   integer, public, parameter :: n_movies_max=30  ! Max no. of different movies
   integer, public, parameter :: n_movie_fields_max=6 ! Max no. of fields per movie
   real(cp), public ::  movie_const(n_movies_max)
   character(len=80), public :: movie(n_movies_max)  ! Only for input
   character(len=72), public :: movie_file(n_movies_max)

   logical, public :: lStoreMov(n_movies_max),lICField(n_movies_max)
   integer, public :: n_movies
   integer, public :: n_movie_type(n_movies_max)
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

      if ( .not. l_movie ) then
         l_movie_oc=.false.
         l_movie_ic=.false.
         l_HTmovie =.false.
         l_dtBmovie=.false.
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
            do n=1,n_movies_max
               n_movie_file(n)=70+n
            end do
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
      !     * n_type(n_movie) =  movie type:
      !
      !                   - =  1  : Radial magnetic field
      !                   - =  2  : Theta component of magnetic field
      !                   - =  3  : Azimuthal magnetic field
      !                   - =  4  : Horizontal magnetic field
      !                   - =  5  : Total magnetic field (all compnents)
      !                   - =  8  : Axisymmetric azimuthal
      !                     magnetic field (phi=constant)
      !                   - =  9  : 3d magnetic field
      !                   - = 11  : Radial velocity field
      !                   - = 12  : Theta component of velocity field
      !                   - = 13  : Azimuthal velocity field
      !                   - = 14  : Horizontal velocity field
      !                   - = 15  : Total velocity field (all compnents)
      !                   - = 17  : Scalar field whose contours are the
      !                     stream lines of the axisymm. poloidal
      !                     velocity field (phi=constant)
      !                   - = 18  : Axisymmetric azimuthal
      !                     velocity field (phi=constant)
      !                   - = 19  : 3d velocity field
      !                   - = 20  : z component of vorticity
      !                   - = 21  : Temperature field
      !                   - = 22  : radial conv. heat transport
      !                   - = 23  : helicity
      !                   - = 24  : axisymmetric helicity
      !                   - = 25  : phi component of vorticity
      !                   - = 26  : radial component of vorticity
      !                   - = 28  : axisymmetric Temperature field
      !                     for phi=const.
      !                   - = 29  : 3d temperature field
      !                   - = 30  : Scalar field whose contours are the
      !                     fieldlines of the axisymm. poloidal
      !                     magnetic field (phi=constant)
      !                   - = 31  : field line production
      !                   - = 32  : field line advection
      !                   - = 33  : field line diffusion
      !                   - = 40  : Axisymmetric azimuthal
      !                     magnetic field (phi=constant)
      !                   - = 41  : Axis. Bphi production +  omega eff.
      !                   - = 42  : Axis. Bphi advection
      !                   - = 43  : Axis. Bphi diffusion
      !                   - = 44  : Axis. Bphi str.,dyn.,omega,diff.
      !                   - = 50  : Bz
      !                   - = 51  : Bz production
      !                   - = 52  : Bz advection
      !                   - = 53  : Bz diffusion
      !                   - = 60  : toroidal Bphi
      !                   - = 61  : toroidal Bphi production + omega eff.
      !                   - = 62  : toroidal Bphi advection
      !                   - = 63  : toroidal Bphi diffusion
      !                   - = 71  : Br production
      !                   - = 72  : Br advection
      !                   - = 73  : Br diffusion
      !                   - = 80  : Jr
      !                   - = 81  : Jr production
      !                   - = 82  : Jr advection
      !                   - = 83  : Jr diffusion
      !                   - = 90  : poloidal Jz pol.
      !                   - = 91  : poloidal Jz pol. production
      !                   - = 92  : poloidal Jz advection
      !                   - = 93  : poloidal Jz diffusion
      !                   - = 94  : z component of velovity
      !                   - = 95  : toroidal Btheta
      !                   - = 96  : toroidal Potential
      !                   - = 97  : Function for Poloidal Fieldlines
      !                   - = 98  : azimuthal shear of Br
      !                   - = 99  : phi component of Lorentz force
      !                   - =101  : Stress fields
      !                   - =102  : Force fields
      !                   - =103  : Br Inverse appearence at CMB
      !                   - =110  : radial heat flow
      !                   - =111  : Vz and Vorz north/south correlation
      !                   - =112  : axisymm dtB tersm for Br and Bp
      !                   - =113  : axisymm dSdr
      !                   - =114  : Cylindrically radial magnetic field
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
      !                      - =17 : radial derivative of T * vr
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
      !                      - =91 : radial derivative of T
      !                      - =92 : Vz north/south correlation
      !                      - =93 : Vorz north/south correlation
      !                      - =94 : Hel north/south correlation
      !                      - =101: AS poloidal Br production
      !                      - =102: AS poloidal Br dynamo term
      !                      - =103: AS poloidal Br diffusion
      !                      - =104: AS toroidal Bp production
      !                      - =105: AS toroidal Bp dynamo term
      !                      - =106: AS toroidal Bp omega effect
      !                      - =107: AS toroidal Bp diffusion
      !                      - =108: Bs
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
      logical :: lEquator
      integer :: length,length_fn,lengthC
      character(len=200) :: message
      character(len=80) :: string,word,stringC
      character(len=80) :: file_name
      character(len=50) :: typeStr
      integer :: n_theta,n_phi
      real(cp) :: r_movie,theta_movie,phi_movie
      real(cp) :: phi_max
      real(cp) :: rad
      real(cp) :: const
      integer :: i,n,n_ic
      integer :: ns
      integer :: n_type
      integer :: n_surface
      integer :: n_const
      integer :: n_fields,n_fields_oc,n_fields_ic
      integer :: n_field_size
      integer :: n_field_size_ic
      integer :: n_field_start
      integer :: n_field_type(n_movie_fields_max)
      integer :: n_rc,n_tc,n_pc

      logical :: lStore,lIC,foundGridPoint

      !--- Initialize first storage index:
      n_field_start=1

      !--- Converts from radiant to degree:
      rad=180.0_cp/pi

      !--- Loop over max possible no of movies:
      l_movie_oc   =.false.
      l_movie_ic   =.false.
      l_HTmovie    =.false.
      l_dtBmovie   =.false.
      l_store_frame=.false.

      n_rc=0
      n_tc=0
      n_pc=0

      do i=1,n_movies_max

         lStore=.true.
         lIC   =.false.

         string=movie(i)

         if ( len(trim(string))  ==  0 ) cycle !blank string

         !--- Delete blanks, they are not interpreted
         call delete_string(string,' ',length)

         !--- Convert to capitals:
         call capitalize(string)

         lEquator=.false.

         !--- Identify movie type (fields to be stored):

         if ( index(string,'BR') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               n_type=71
               typeStr=' radial magnetic field production '
               file_name='BrPro_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=27
            else if ( index(string,'ADV') /= 0 ) then
               n_type=72
               typeStr=' radial magnetic field advection '
               file_name='BrAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=28
            else if ( index(string,'DIF') /= 0 ) then
               n_type=73
               typeStr=' radial magnetic field diffusion '
               file_name='BrDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=29
            else if ( index(string,'INV') /= 0 ) then
               n_type=103
               typeStr=' radial magnetic field diffusion '
               file_name='BrInv_'
               lIC=.false.
               lStore=.false.
               n_fields=1
               n_field_type(1)=81
            else if ( index(string,'SHEAR') /= 0 ) then
               n_type=98
               typeStr=' azimuthal shear of radial magnetic field'
               file_name='BrShear'
               lIC=.false.
               lStore=.true.
               n_fields=1
               n_field_type(1)=53
            else
               n_type=1
               typeStr=' radial magnetic field '
               n_fields=1
               file_name='Br_'
               lIC=.true.
               n_field_type(1)=1
            end if
         else if ( index(string,'BT') /= 0 ) then
            if ( index(string,'TOR') /= 0 ) then
               n_type=95 ! Toroidal B theta
               typeStr=' toroidal B theta'
               file_name='BtTor_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=50
            else
               n_type=2
               typeStr=' theta comp. of magnetic field '
               n_fields=1
               file_name='Bt_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=2
            end if
         else if ( index(string,'BP' ) /= 0 .and. index(string,'TOR') == 0 .and. &
         &    index(string,'AX' ) == 0 .and. index(string,'AB' ) == 0 ) then
            n_type=3
            typeStr=' azimuthal magnetic field '
            file_name='Bp_'
            lIC=.true.
            n_fields=1
            n_field_type(1)=3
         else if ( index(string,'BH') /= 0 .or. index(string,'BM') /= 0 ) then
            n_type=4 ! Horizontal field
            file_name='Bh_'
            lIC=.true.
         else if( index(string,'BS') /= 0) then
            n_type=114
            typeStr='cyl radial magnetic field'
            file_name='Bs_'
            lIC=.true.
            n_fields=1
            n_field_type(1)=108
         else if ( index(string,'BALL') /= 0 ) then
            n_type=5 ! Total field
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
               n_type=97 ! Poloidal field lines
               typeStr=' pol. fieldlines'
               file_name='FLPol_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=52
            else if ( index(string,'PRO') /= 0 ) then
               n_type=31 ! Poloidal field lines production
               typeStr=' axisymm. fieldline production '
               file_name='FLPro_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=20
            else if ( index(string,'ADV') /= 0 ) then
               n_type=32 ! Poloidal field lines advection
               typeStr=' axisymm. fieldline advection '
               file_name='FLAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=21
            else if ( index(string,'DIF') /= 0 ) then
               n_type=33 ! Poloidal field lines diffusion
               typeStr=' axisymm. fieldline diffusion '
               file_name='FLDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=22
            else
               n_type=30
               typeStr=' axisymm. fieldlines '
               file_name='FL_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=8
            end if
         else if ( index (string,'LORENTZ') /= 0 .or. index(string,'LF') /= 0 ) then
            n_type=99
            typeStr=' phi Lorentz force '
            file_name='LFp_'
            lIC=.false.
            n_fields=1
            n_field_type(1)=54
         else if ( ( index(string,'AX') /= 0 .and.  &
         &    index(string,'BP') /= 0 ) .or. index(string,'AB') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               n_type=41
               typeStr=' axisymm. B phi production '
               file_name='ABPro_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=23
            else if ( index(string,'ADV') /= 0 ) then
               n_type=42
               typeStr=' axisymm. B phi advection '
               file_name='ABAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=25
            else if ( index(string,'DIF') /= 0 ) then
               n_type=43
               typeStr=' axisymm. B phi diffusion '
               file_name='ABDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=26
            else
               n_type=40
               typeStr=' axisymm. B phi '
               file_name='AB_'
               lIC=.true.
               n_fields=1
               n_field_type(1)=9
            end if
         else if ( index(string,'JR') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               n_type=81
               typeStr=' radial current production '
               file_name='JrPro_'
               lIC=.true.
               lStore=.false.
               n_fields=2
               n_field_type(1)=31
               n_field_type(2)=32
            else if ( index(string,'ADV') /= 0 ) then
               n_type=82
               typeStr=' radial current advection '
               file_name='JrAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=33
            else if ( index(string,'DIF') /= 0 ) then
               n_type=83
               typeStr=' radial current diffusion '
               file_name='JrDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=34
            else
               n_type=80
               typeStr=' radial current '
               file_name='Jr_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=30
            end if
         else if ( index(string,'BZ') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               n_type=51
               typeStr=' poloidal Bz production '
               file_name='BzPolPro_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=35
            else if ( index(string,'ADV') /= 0 ) then
               n_type=52
               typeStr=' poloidal Bz advection '
               file_name='BzPolAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=36
            else if ( index(string,'DIF') /= 0 ) then
               n_type=53
               typeStr=' poloidal Bz diffusion '
               file_name='BzPolDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=37
            else
               n_type=50
               typeStr=' poloidal Bz '
               file_name='BzPol_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=13
            end if
         else if ( index(string,'BP') /= 0 .and. index(string,'TOR') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               n_type=61
               typeStr=' toroidal Bphi production '
               file_name='BpTorPro_'
               lIC=.true.
               lStore=.false.
               n_fields=2
               n_field_type(1)=43
               n_field_type(2)=49
            else if ( index(string,'ADV') /= 0 ) then
               n_type=62
               typeStr=' toroidal Bphi advection '
               file_name='BpTorAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=45
            else if ( index(string,'DIF') /= 0 ) then
               n_type=63
               typeStr=' toroidal Bphi diffusion '
               file_name='BpTorDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=46
            else
               n_type=60
               typeStr=' toroidal Bphi '
               file_name='BpTor_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=42
            end if
         else if ( index(string,'JZ') /= 0 ) then
            if ( index(string,'PRO') /= 0 ) then
               n_type=91
               typeStr=' poloidal Jz production '
               file_name='JzTorPro_'
               lIC=.true.
               lStore=.false.
               n_fields=2
               n_field_type(1)=38
               n_field_type(2)=39
            else if ( index(string,'ADV') /= 0 ) then
               n_type=92
               typeStr=' poloidal Jz advection '
               file_name='JzTorAdv_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=40
            else if ( index(string,'DIF') /= 0 ) then
               n_type=93
               typeStr=' poloidal Jz diffusion '
               file_name='JzTorDif_'
               lIC=.true.
               lStore=.false.
               n_fields=1
               n_field_type(1)=41
            else
               n_type=90
               typeStr=' poloidal Jz '
               file_name='JzTor_'
               lIC=.true.
               lStore=.true.
               n_fields=1
               n_field_type(1)=14
            end if
         else if ( index(string,'VR') /= 0 ) then
            n_type=11
            typeStr=' radial velocity field '
            file_name='Vr_'
            n_fields=1
            n_field_type(1)=4
         else if ( index(string,'VT') /= 0 ) then
            n_type=12
            typeStr=' theta comp. of velocity field '
            file_name='Vt_'
            n_fields=1
            n_field_type(1)=5
         else if ( index(string,'VP') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=18
               typeStr=' axisym. azimuthal velocity field '
               file_name='AV_'
               n_fields=1
               n_field_type(1)=11
            else
               n_type=13
               typeStr=' phi comp. of velocity field '
               file_name='Vp_'
               n_fields=1
               n_field_type(1)=6
            end if
         else if ( index(string,'VH') /= 0 .or. index(string,'VM') /= 0 ) then
            n_type=14
            file_name='Vh_'
         else if ( index(string,'VA') /= 0 ) then
            n_type=15
            typeStr=' all velocity components '
            file_name='V_'
            n_fields=3
            n_field_type(1)=4
            n_field_type(2)=5
            n_field_type(3)=6
         else if ( index(string,'STREAMLINE') /= 0 .or. index(string,'SL') /= 0 ) then
            n_type=17
            typeStr=' axisym. meridional velocity streamlines '
            file_name='SL_'
            n_fields=1
            n_field_type(1)=10
         else if ( index(string,'VOR') /= 0 ) then
            if ( index(string,'Z') /= 0 ) then
               n_type=20
               typeStr=' z-component of vorticity '
               file_name='VorZ_'
               n_fields=1
               n_field_type(1)=16
            else if ( index(string,'P') /= 0 ) then
               n_type=20
               typeStr=' phi-component of vorticity '
               file_name='VorP_'
               n_fields=1
               n_field_type(1)=47
            else if ( index(string,'R') /= 0 ) then
               n_type=20
               typeStr=' r-component of vorticity '
               file_name='VorR_'
               n_fields=1
               n_field_type(1)=48
            end if
         else if ( index(string,'VS2') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=118
               typeStr=' axisym. vs**2 '
               file_name='AVS2_'
               n_fields=1
               n_field_type(1)=98
            end if
         else if ( index(string,'VS') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=115
               typeStr=' axisym. s-component of velocity '
               file_name='AVS_'
               n_fields=1
               n_field_type(1)=94
            end if
         else if ( index(string,'REYS') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=116
               typeStr=' axisym. vs*vphi '
               file_name='AReyS_'
               n_fields=1
               n_field_type(1)=95
            end if
         else if ( index(string,'REYZ') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=117
               typeStr=' axisym. vz*vphi '
               file_name='AReyZ_'
               n_fields=1
               n_field_type(1)=97
            end if
         else if ( index(string,'VZ2') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=119
               typeStr=' axisym. vz**2 '
               file_name='AVZ2_'
               n_fields=1
               n_field_type(1)=99
            end if
         else if ( index(string,'VZ') /= 0 ) then
            if ( index(string,'AX') /= 0 ) then
               n_type=120
               typeStr=' axisym. vz '
               file_name='AVZ_'
               n_fields=1
               n_field_type(1)=96
            else
               n_type=94
               typeStr=' z-component of velocity '
               file_name='VZ_'
               n_fields=1
               n_field_type(1)=15
            end if
         else if ( ( index(string,'AX'  ) /= 0 .and.  &
         &    index(string,'HEL' ) /= 0 ) .or. index(string,'AHEL') /= 0 ) then
            n_type=24
            typeStr=' axisymmetric helicity '
            file_name='AH_'
            n_fields=1
            n_field_type(1)=19
         else if ( index(string,'HEL') /= 0 ) then
            n_type=23
            typeStr=' helicity '
            file_name='HE_'
            n_fields=1
            n_field_type(1)=18
         else if ( index(string,'AX' ) /= 0 .and. &
         &    ( index(string,'TEM') /= 0 .or. index(string,'ENT') /= 0 ) ) then
            n_type=28
            typeStr=' axisymmetric temp. '
            file_name='AT'
            n_fields=1
            n_field_type(1)=12
         else if ( index(string,'ENTROPY') /= 0 .or. index(string,'TEM') /= 0 ) then
            ns=index(string,'S')
            if ( ns > 0 ) then
               if ( string(ns:ns+2) == 'SUR' ) then
                  call abortRun('! No surface T field available !')
               end if
            end if
            n_type=21
            typeStr=' temperature field '
            file_name='T_'
            n_fields=1
            n_field_type(1)=7
         else if ( ( index(string,'CONV' ) /= 0 .and.  &
         &    index(string,'HEAT' ) /= 0 ) .or. index(string,'HEATT') /= 0 ) then
            ns=index(string,'S')
            if ( ns > 0 ) then
               if ( string(ns:ns+2) == 'SUR' ) then
                  call abortRun('! No surface T field available !')
               end if
            end if
            n_type=22
            typeStr=' radial convective heat transport '
            file_name='HT_'
            l_HTmovie=.true.
            n_fields=1
            n_field_type(1)=17
         else if ( index(string,'AX' ) /= 0 .and. index(string,'HEATF') /= 0   ) then
            n_type=113
            typeStr='axisymmetric dSdr '
            file_name='AHF_'
            l_HTmovie=.true.
            n_fields=1
            n_field_type(1)=92
         else if ( index(string,'HEATF') /= 0 ) then
            ns=index(string,'S')
            if ( ns > 0 ) then
               if ( string(ns:ns+2) == 'SUR' ) then
                  call abortRun('! No surface T field available !')
               end if
            end if
            n_type=110
            typeStr=' radial heat transport '
            file_name='HF_'
            l_HTmovie=.true.
            n_fields=1
            n_field_type(1)=91
         else if ( index(string,'POT') /= 0 .and. index(string,'TOR') /= 0 ) then
            n_type=96
            typeStr=' radial convective heat transport '
            file_name='PotTor_'
            lIC=.true.
            lStore=.false.
            n_fields=1
            n_field_type(1)=51
         else
            message = 'Couldnt interpret movie field from string:'//string
            call abortRun(message)
         end if


         !--- Identify surface type:

         length_fn=len(trim(file_name))
         if ( n_type == 103 ) then
            n_surface=1 !
            n_const=1   !
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=r_cmb
         else if (   index(string,'AX') /= 0 .or.                     &
         &    file_name(1:2) == 'AV' .or. file_name(1:2) == 'AB' .or. &
         &    n_type == 30 .or. n_type == 31 .or. n_type == 32 .or.   &
         &    n_type == 33 .or. n_type == 40 .or. n_type == 41 .or.   &
         &    n_type == 42 .or. n_type == 43 .or. n_type == 17 ) then
            !--- Axisymmetric stuff:
            n_surface=3  ! PHI=const.
            n_const=1
            n_field_size=n_r_max*n_theta_max
            n_field_size_ic=n_r_ic_max*n_theta_max
            const=0.0_cp
         else if ( index(string,'3D') /= 0 ) then
            n_surface=0  ! 3d
            n_const=0    ! Not needed
            file_name=file_name(1:length_fn)//'3D_'
            n_field_size=n_r_max*n_theta_max*n_phi_max
            n_field_size_ic=n_r_ic_max*n_theta_max*n_phi_max
         else if ( index(string,'CMB') /= 0 ) then
            n_surface=1 ! R=const. at CMB
            n_const=1
            file_name=file_name(1:length_fn)//'CMB_'
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=r_cmb
         else if ( index(string,'ICB') /= 0 ) then
            n_surface=1 ! R=const. at ICB
            n_const=n_r_max
            file_name=file_name(1:length_fn)//'ICB_'
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=r_icb
         else if ( index(string,'SUR') /= 0 ) then
            n_surface=-1
            n_const=1
            file_name=file_name(1:length_fn)//'SUR_'
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size
            const=one
         else if ( index(string,'R=') /= 0 .or. &
         &    index(string,'RAD=') /= 0 .or. index(string,'RADIUS=') /= 0 ) then
            n_surface=1  ! R=const.
            n_rc=n_rc+1
            write(stringC,'(''R=C'',i1,''_'')') n_rc
            lengthC=length_to_blank(stringC)
            file_name=file_name(1:length_fn)//stringC(1:lengthC)
            n_field_size=n_phi_max*n_theta_max
            n_field_size_ic=n_field_size

            !------ Get Radius in fractions of outer core radius:
            if ( index(string,'R=') /= 0 ) then
               word=string(index(string,'R=')+2:length)
            else if ( index(string,'RAD=') /= 0 ) then
               word=string(index(string,'RAD=')+4:length)
            else if ( index(string,'RADIUS=') /= 0 ) then
               word=string(index(string,'RADIUS=')+7:length)
            end if
            call str2dble(word,r_movie)

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
               r_movie=-r_movie*r_icb
               n_const = minloc(abs(r_ic - r_movie),1)
               const = r_ic(n_const)
            end if

         else if ( index(string,'EQ') /= 0 .or. lEquator ) then

            n_surface=2    ! Equator
            n_const=n_theta_max
            file_name=file_name(1:length_fn)//'EQU_'
            n_field_size=n_r_max*n_phi_max
            n_field_size_ic=n_r_ic_max*n_phi_max
            const=rad*theta(n_const)

         else if ( index(string,'T=') /= 0 .or. index(string,'THETA=') /= 0 ) then

            n_surface=2    ! Theta=const.
            n_tc=n_tc+1
            write(stringC,'(''T=C'',i1,''_'')') n_tc
            lengthC=length_to_blank(stringC)
            file_name=file_name(1:length_fn)//stringC(1:lengthC)
            n_field_size=n_r_max*n_phi_max
            n_field_size_ic=n_r_ic_max*n_phi_max

            !------ Get desired colatitude value:
            if ( index(string,'T=') /= 0 ) then
               word=string(index(string,'T=')+2:length)
            else if ( index(string,'THETA=') /= 0 ) then
               word=string(index(string,'THETA=')+6:length)
            end if
            call str2dble(word,theta_movie)
            theta_movie=abs(theta_movie)
            theta_movie=theta_movie/rad

            !------ Choose closest colatitude grid point:
            foundGridPoint=.false.
            do n_theta=1,n_theta_max-1
               if ( theta(n_theta)  <= theta_movie .and. &
               &    theta(n_theta+1) >= theta_movie ) then
                  if ( theta(n_theta+1)-theta_movie < &
                  &    theta_movie-theta(n_theta) ) then
                     n_const=n_theta+1
                  else
                     n_const=n_theta
                  end if
                  foundGridPoint=.true.
                  exit
               end if
            end do
            if ( .not. foundGridPoint ) then
               if ( theta_movie-theta(n_theta_max) <= &
               &    theta(1)+180.0_cp/rad-theta_movie ) then
                  n_const=n_theta_max
               else
                  n_const=1
               end if
            end if
            const=rad*theta(n_const)

            !---------- Now switch to north/south order of thetas:
            if ( n_const < n_theta_max/2 ) then
               n_const=2*n_const-1
            else
               n_const=2*(n_theta_max-n_const+1)
            end if

         else if ( index(string,'MER' ) /= 0 .or.  &
              index(string,'P='  ) /= 0  .or. index(string,'PHI=') /= 0 ) then

            n_surface=3  ! PHI=const.
            n_pc=n_pc+1
            write(stringC,'(''P=C'',i1,''_'')') n_pc
            lengthC=length_to_blank(stringC)
            file_name=file_name(1:length_fn)//stringC(1:lengthC)
            n_field_size=2*n_r_max*n_theta_max
            n_field_size_ic=2*n_r_ic_max*n_theta_max

            !------ Get desired longitude value:
            if ( index(string,'P=') /=0 ) then
               word=string(index(string,'P=')+2:length)
            else if ( index(string,'PHI=') /=0 ) then
               word=string(index(string,'PHI=')+4:length)
            end if
            call str2dble(word,phi_movie)
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
            foundGridPoint=.false.
            do n_phi=1,n_phi_max-1
               if ( phi(n_phi)  <= phi_movie .and. phi(n_phi+1) >= phi_movie ) then
                  if ( phi(n_phi+1)-phi_movie < phi_movie-phi(n_phi) ) then
                     n_const=n_phi+1
                  else
                     n_const=n_phi
                  end if
                  foundGridPoint=.true.
                  exit
               end if
            end do
            if ( .not. foundGridPoint ) then
               if ( phi_movie-phi(n_phi_max) <= phi(1)+phi_max-phi_movie ) then
                  n_const=n_phi_max
               else
                  n_const=1
               end if
            end if

            const=rad*phi(n_const)

         else
            message = 'Couldnt interpret movie surface from string:'//string
            call abortRun(message)
         end if

         if ( n_field_type(1) == 54 .and.                    &
         &    ( ( n_surface /= 0 .and. n_surface /= 3 ) .or. &
         &    index(string,'AX') /= 0 ) ) then
            write(*,*) 'Sorry, can only prepare movie file for'
            write(*,*) 'phi component of LF in 3d or for phi-cut.'
            call abortRun('Stop run in movie')
         end if

         !--- Now store the necessary information:

         !------ Increase number of movies:
         n_movies=n_movies+1
         lStoreMov(n_movies)=lStore
         lICField(n_movies)=lIC
         if ( .not. lStore .and. n_field_type(1) /= 13 .and.          &
         &    n_field_type(1) /= 14 .and. n_field_type(1) /= 30 .and. &
         &    n_field_type(1) /= 42 .and. n_field_type(1) /= 50 .and. &
         &    n_field_type(1) /= 51 .and. n_field_type(1) /= 52 .and. &
         &    n_field_type(1) /= 54 ) l_dtBmovie= .true.

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
         length_fn=len(trim(file_name))
         file_name=file_name(1:length_fn)//'mov.'//tag
         call delete_string(file_name,' ',length)
         movie_file(n_movies)=file_name(1:length)

         !------ Store movie type:
         n_movie_type(n_movies)=n_type

         !------ Store information about movie surface:
         n_movie_surface(n_movies)=n_surface
         n_movie_const(n_movies)  =n_const
         movie_const(n_movies)    =const

         !------ Store number of fields for this movie:
         n_movie_fields(n_movies)   =n_fields_oc
         n_movie_fields_ic(n_movies)=n_fields_ic

         !------ Store size of field and where it should be stored in
         !       the work array frames(*) (see c_movie.f):
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
            else if ( n_surface == 0 ) then
               write(n_log_file,*) '!    in 3d !'
            else if ( n_surface == 1 .and. n_const == 1 ) then
               write(n_log_file,*) '!    at the core-mantle boundary !'
            else if ( n_surface == 1 .and. n_const > 1 ) then
               write(n_log_file,'('' !    at r='',f10.6)') r(n_const)
            else if ( n_surface == 1 .and. n_const < 0 ) then
               write(n_log_file,'('' !    at r='',f10.6)') r_ic(n_const)
            else if ( n_surface == 2 ) then
               write(n_log_file,'('' !    at theta='',f12.6)') rad*theta(n_const)
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
   subroutine movie_gather_frames_to_rank0
      !
      ! MPI communicators for movie files
      !

#ifdef WITH_MPI
      integer :: n_fields,n_surface,n_movie,n_const
      integer :: n_start,n_stop,n_field
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
      if (rank == 0) then
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
         if ( n_surface == -1 ) then ! Earth Surface
            ! theta-phi surface for n_r=1 (CMB)
            ! frames is already existent on rank 0 with all
            ! needed values
            ! do nothing, pass to the next movie

            cycle

         else if ( n_surface == 0 ) then ! 3d
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
               if (rank == 0) then
                  frames(n_start:n_stop)=field_frames_global(1:field_length)
               end if
            end do

         else if ( n_surface == 1 ) then ! Surface r=constant
            ! frames is set only for one rank, where n_r=n_const
            ! send to rank 0
            do n_field=1,n_fields
               n_start=n_movie_field_start(n_field,n_movie)
               n_stop =n_movie_field_stop(n_field,n_movie)
               field_length = n_stop-n_start+1
               myTag=7654+n_movie
               if (rank == 0) then
                  if ((nRstart <= n_const) .and. (n_const <= nRstop)) then
                     ! relevant frames already set on rank 0
                     ! do nothing
                  else
                     call MPI_Recv(frames(n_start),field_length,MPI_DEF_REAL,     &
                          &        MPI_ANY_SOURCE,mytag,MPI_COMM_WORLD,status,ierr)
                  end if
               else
                  if ((nRstart <= n_const) .and. (n_const <= nRstop)) then
                     ! relevant frames are all on this rank  /= 0
                     ! send to rank 0
                     call MPI_Send(frames(n_start),field_length,MPI_DEF_REAL,&
                          &        0,mytag,MPI_COMM_WORLD,ierr)
                  end if
               end if
            end do

         else if ( n_surface == 2 ) then ! Surface theta=constant
            do n_field=1,n_fields
               n_start = n_movie_field_start(n_field,n_movie)
               n_stop  = n_movie_field_stop(n_field,n_movie)
               field_length = n_stop-n_start+1

               local_start=n_start+(nRstart-1)*n_phi_max
               local_end  =local_start+nR_per_rank*n_phi_max-1
               if (local_end > n_stop) then
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
               if (rank == 0) then
                  frames(n_start:n_stop)=field_frames_global(1:field_length)
               end if
            end do  ! Do loop over field for one movie
         else if ( abs(n_surface) == 3 ) then  ! Surface phi=const.
            ! all ranks have a part of the frames array for each movie
            ! we need to gather

            do n_field=1,n_fields
               n_start = n_movie_field_start(n_field,n_movie)
               n_stop  = n_movie_field_stop(n_field,n_movie)
               field_length = n_stop-n_start+1

               local_start=n_start+(nRstart-1)*n_theta_max
               local_end  =local_start+nR_per_rank*n_theta_max-1
               if (local_end > n_stop) then
                  call abortRun('local_end exceeds n_stop')
               end if
               do irank=0,n_procs-1
                  recvcounts(irank)=radial_balance(irank)%n_per_rank*n_theta_max
               end do
               displs(0)=0
               do irank=1,n_procs-1
                  displs(irank)=displs(irank-1)+recvcounts(irank-1)
               end do
               sendcount=local_end-local_start+1
               call MPI_Gatherv(frames(local_start),sendcount,MPI_DEF_REAL, &
                    &           field_frames_global,recvcounts,displs,      &
                    &           MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
               if (rank == 0) then
                  frames(n_start:n_stop)=field_frames_global(1:field_length)
               end if
            end do  ! Do loop over field for one movie

         end if
      end do

      if (rank == 0) deallocate(field_frames_global)
#endif

   end subroutine movie_gather_frames_to_rank0
!----------------------------------------------------------------------------
end module movie_data
