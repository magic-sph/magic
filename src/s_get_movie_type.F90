!$Id$
!***********************************************************************
SUBROUTINE get_movie_type
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to identify the different movie    |
  !  |  types from the input string movies(*).                           |
  !  |  Note that generally blanks are not interpreted and that the      |
  !  |  interpretation is not case sensitive.                            |
  !  |  In general two informations are needed:                          |
  !  |    1) A word FIELDINFO that identifies the field to be plotted    |
  !  |         (e.g. Br for radial magnetic field, see list below)       |
  !  |      Possible keywords are [optional text in brackets             |
  !  |           B r[adial]     : radial magnetic field                  |
  !  |           B t[heta]      : theta component                        |
  !  |           B p[hi]        : azimuthal component                    |
  !  |           B h[orizontal] : the two horizontal components          |
  !  |           B a[ll]        : all three components                   |
  !  |           FIELDLINE[S]   : field lines of axisymmetric            |
  !  |           or FL            poloidal field for phi=constant        |
  !  |           AX[ISYMMETRIC] B                                        |
  !  |           or AB          : axisymmetric phi component of the      |
  !  |                            magnetic field for phi=constant        |
  !  |           V r[adial]     : radial velocity field                  |
  !  |           V t[heta]      : theta component                        |
  !  |           V p[hi]        : azimuthal component                    |
  !  |           V h[orizontal] : the two horizontal components          |
  !  |           V a[ll]        : all three components                   |
  !  |           STREAMLINE[S]  : field lines of axisymmetric            |
  !  |           or SL          : poloidal field for phi=constant        |
  !  |           AX[ISYMMETRIC] V                                        |
  !  |           or AV          : axisymmetric phi component of the      |
  !  |                            velocity field for phi=constant        |
  !  |           V z            : z component of velocity at equator     |
  !  |                            + z component of the vorticity at      |
  !  |                            the equator (closest point to equator) |
  !  |           Vo z          : z-component of vorticity                |
  !  |           Vo r          : r-component of vorticity                |
  !  |           Vo p          : phi-component of vorticity              |
  !  |           T[emperature]  : sic                                    |
  !  |           AX[ISYMMETRIC] T                                        |
  !  |           or AT          : axisymmetric T field for phi=constant  |
  !  |           Heat t[ransport]: radial derivative of T                |
  !  |                                                                   |
  !  |           FL Pro         : axisymmetric field line stretching     |
  !  |           FL Adv         : axisymmetric field line advection      |
  !  |           FL Dif         : axisymmetric field line diffusion      |
  !  |           AB Pro         : axisymmetric (tor.) Bphi production    |
  !  |           AB Dif         : axisymmetric (tor.) Bphi diffusion     |
  !  |           Br Pro         : Br production                          |
  !  |           Br Adv         : Br advection                           |
  !  |           Br Dif         : Br diffusion                           |
  !  |           Jr             : Jr production                          |
  !  |           Jr Pro         : Jr production +  omega effects         |
  !  |           Jr Adv         : Jr advection                           |
  !  |           Jr Dif         : Jr diffusion                           |
  !  |           Bz Pol         : poloidal Bz                            |
  !  |           Bz Pol Pro     : poloidal Bz production                 |
  !  |           Bz Pol Adv     : poloidal Bz advection                  |
  !  |           Bz Pol Dif     : poloidal Bz diffusion                  |
  !  |           Jz Tor         : poloidal Jz                            |
  !  |           Jz Tor Pro     : poloidal Jz production                 |
  !  |           Jz Tor Adv     : poloidal Jz advection                  |
  !  |           Jz Tor Dif     : poloidal Jz diffusion                  |
  !  |           Bp Tor         : toriodal Bphi                          |
  !  |           Bp Tor Pro     : toriodal Bphi production               |
  !  |           Bp Tor Adv     : toriodal Bphi advection                |
  !  |           Bp Tor Dif     : toriodal Bphi diffusion                |
  !  |           HEL[ICITY]     : sic                                    |
  !  |           AX[ISYMMETRIC HELICITY]  or                             |
  !  |           AHEL           : axisymmetric helicity                  |
  !  |           Bt Tor         : toroidal Btheta                        |
  !  |           Pot Tor        : toroidal Potential                     |
  !  |           Pol Fieldlines : toroidal Potential                     |
  !  |           Br Shear       : azimuthal Shear of Br                  |
  !  |           Lorentz[force] : Lorentz force (only phi component)     |
  !  |           Br Inv         : Inverse field apperance at CMB         |
  !  |                                                                   |
  !  |    2) A second information that identifies the coordinate         |
  !  |       to be kept constant (surface).                              |
  !  |       (e.g. r=number for surface r=constant with number given     |
  !  |        in units of the total core radius or                       |
  !  |        theta/phi=number with number given in degrees              |
  !  |       Three keyword are also possible:                            |
  !  |           CMB       : core mantle boundary                        |
  !  |           SUR[FACE] : Earth surface (only magnetic field)         |
  !  |           3[D]      : 3D field throughout the OC [and IC for B]   |
  !  |                                                                   |
  !  |  On output the necessary information is coded into integers       |
  !  |  and is used in this form by further subroutines:                 |
  !  |     n_movies = total number of movies                             |
  !  |     n_type(n_movie) =  movie type                                 |
  !  |                     =  1  : Radial magnetic field                 |
  !  |                     =  2  : Theta component of magnetic field     |
  !  |                     =  3  : Azimuthal magnetic field              |
  !  |                     =  4  : Horizontal magnetic field             |
  !  |                     =  5  : Total magnetic field (all compnents)  |
  !  |                     =  8  : Axisymmetric azimuthal                |
  !  |                             magnetic field (phi=constant)         |
  !  |                     =  9  : 3d magnetic field                     |
  !  |                     = 11  : Radial velocity field                 |
  !  |                     = 12  : Theta component of velocity field     |
  !  |                     = 13  : Azimuthal velocity field              |
  !  |                     = 14  : Horizontal velocity field             |
  !  |                     = 15  : Total velocity field (all compnents)  |
  !  |                     = 17  : Scalar field whose contours are the   |
  !  |                             stream lines of the axisymm. poloidal |
  !  |                             velocity field (phi=constant)         |
  !  |                     = 18  : Axisymmetric azimuthal                |
  !  |                             velocity field (phi=constant)         |
  !  |                     = 19  : 3d velocity field                     |
  !  |                     = 20  : z component of vorticity              |
  !  |                     = 21  : Temperature field                     |
  !  |                     = 22  : radial conv. heat transport           |
  !  |                     = 23  : helicity                              |
  !  |                     = 24  : axisymmetric helicity                 |
  !  |                     = 25  : phi component of vorticity            |
  !  |                     = 26  : radial component of vorticity         |
  !  |                     = 28  : axisymmetric Temperature field        |
  !  |                             for phi=const.                        |
  !  |                     = 29  : 3d temperature field                  |
  !  |
  !  |                     = 30  : Scalar field whose contours are the   |
  !  |                             fieldlines of the axisymm. poloidal   |
  !  |                             magnetic field (phi=constant)         |
  !  |                     = 31  : field line production                 |
  !  |                     = 32  : field line advection                  |
  !  |                     = 33  : field line diffusion                  |
  !  |                     = 40  : Axisymmetric azimuthal                |
  !  |                             magnetic field (phi=constant)         |
  !  |                     = 41  : Axis. Bphi production +  omega eff.   |
  !  |                     = 42  : Axis. Bphi advection                  |
  !  |                     = 43  : Axis. Bphi diffusion                  |
  !  |                     = 44  : Axis. Bphi str.,dyn.,omega,diff.      |
  !  |                     = 50  : Bz                                    |
  !  |                     = 51  : Bz production                         |
  !  |                     = 52  : Bz advection                          |
  !  |                     = 53  : Bz diffusion                          |
  !  |                     = 60  : toroidal Bphi                         |
  !  |                     = 61  : toroidal Bphi production + omega eff. |
  !  |                     = 62  : toroidal Bphi advection               |
  !  |                     = 63  : toroidal Bphi diffusion               |
  !  |                     = 71  : Br production                         |
  !  |                     = 72  : Br advection                          |
  !  |                     = 73  : Br diffusion                          |
  !  |                     = 80  : Jr                                    |
  !  |                     = 81  : Jr production                         |
  !  |                     = 82  : Jr advection                          |
  !  |                     = 83  : Jr diffusion                          |
  !  |                     = 90  : poloidal Jz pol.                      |
  !  |                     = 91  : poloidal Jz pol. production           |
  !  |                     = 92  : poloidal Jz advection                 |
  !  |                     = 93  : poloidal Jz diffusion                 |
  !  |                     = 94  : z component of velovity               |
  !  |                     = 95  : toroidal Btheta                       |
  !  |                     = 96  : toroidal Potential                    |
  !  |                     = 97  : Function for Poloidal Fieldlines      |
  !  |                     = 98  : azimuthal shear of Br                 |
  !  |                     = 99  : phi component of Lorentz force        |
  !  |                     =101  : Stress fields                         |
  !  |                     =102  : Force fields                          |
  !  |                     =103  : Br Inverse appearence at CMB          |
  !  |                     =110  : radial heat flow                      |
  !  |                     =111  : Vz and Vorz north/south correlation   |
  !  |                     =112  : axisymm dtB tersm for Br and Bp       |
  !  |                     =113  : axisymm dSdr                          |
  !  |                                                                   |
  !  |     n_movie_surface(n_movie) = defines surface                    |
  !  |     n_movie_surface =  1  : r=constant                            |
  !  |                        2  : theta=constant                        |
  !  |                        3  : phi=constant                          |
  !  |                       -1  : r=constant, Earth surface             |
  !  |                        0  : 3d volume                             |
  !  |     n_movie_fields(n_movie) = no. of fields for outer core        |
  !  |     n_movie_fields_ic(n_movie) = no. of fields for inner core     |
  !  |     n_movie_field_type(n_field,n_movie) = defines field           |
  !  |     n_movie_field_type = 1 : radial magnetic field                |
  !  |                        = 2 : theta comp. of the magnetic field    |
  !  |                        = 3 : azimuthal magnetic field             |
  !  |                        = 4 : radial velocity field                |
  !  |                        = 5 : theta comp. of the velocity field    |
  !  |                        = 6 : azimuthal velocity field             |
  !  |                        = 7 : temperature field                    |
  !  |                        = 8 : scalar field for field lines         |
  !  |                        = 9 : axisymm. toroidal mag. field         |
  !  |                        =10 : scalar field for stream lines        |
  !  |                        =11 : axisymm. v_phi                       |
  !  |                        =12 : axisymm. T                           |
  !  |                        =13 : z-comp. of poloidal Bz               |
  !  |                        =14 : z-comp. of poloidal Jz               |
  !  |                        =15 : z-comp. of velocity                  |
  !  |                        =16 : z-comp. of vorticity                 |
  !  |                        =17 : radial derivative of T * vr          |
  !  |                        =18 : helicity                             |
  !  |                        =19 : axisymmetric helicity                |
  !  |                        =20 : axisymm field-line production        |
  !  |                        =21 : axisymm field-line advection         |
  !  |                        =22 : axisymm field-line diffusion         |
  !  |                        =23 : axisymm Bphi production              |
  !  |                        =24 : axisymm Bphi omega effect            |
  !  |                        =25 : axisymm Bphi advection               |
  !  |                        =26 : axisymm Bphi diffusion               |
  !  |                        =27 : Br production                        |
  !  |                        =28 : Br advection                         |
  !  |                        =29 : Br diffusion                         |
  !  |                        =30 : Jr                                   |
  !  |                        =31 : Jr production                        |
  !  |                        =32 : Jr omega effect                      |
  !  |                        =33 : Jr advection                         |
  !  |                        =34 : Jr diffusion                         |
  !  |                        =35 : poloidal Bz production               |
  !  |                        =36 : poloidal Bz advection                |
  !  |                        =37 : poloidal Bz diffusion                |
  !  |                        =38 : poloidal Jz production               |
  !  |                        =39 : poloidal Jz omega effect             |
  !  |                        =40 : poloidal Jz advection                |
  !  |                        =41 : poloidal Jz diffusion                |
  !  |                        =42 : toroidal Bp                          |
  !  |                        =43 : toroidal Bp production               |
  !  |                        =44 : toroidal Bp omega effect             |
  !  |                        =45 : toroidal Bp advection                |
  !  |                        =46 : toroidal Bp diffusion                |
  !  |                        =47 : phi-comp. of vorticity               |
  !  |                        =48 : r-comp. of vorticity                 |
  !  |                        =49 : toroidal Bp omega effect             |
  !  |                        =50 : toroidal Bt                          |
  !  |                        =51 : toroidal Potential                   |
  !  |                        =52 : poloidal Fieldlines in theta=const   |
  !  |                        =53 : Br dr ( vp/(r sin(theta))            |
  !  |                        =54 : phi Lorentz force                    |
  !  |                        =61 : AS phi reynolds stress force         |
  !  |                        =62 : AS phi advective stress force        |
  !  |                        =63 : AS phi viscous stress force          |
  !  |                        =64 : AS phi Lorentz force                 |
  !  |                        =66 : time derivative of axisym. v phi     |
  !  |                        =67 : relative strength of axisym. v phi   |
  !  |                        =81 : Br inverse appearence at CMB         |
  !  |                        =91 : radial derivative of T               |
  !  |                        =92 : Vz north/south correlation           |
  !  |                        =93 : Vorz north/south correlation         |
  !  |                        =94 : Hel north/south correlation          |
  !  |                        =101: AS poloidal Br production            |
  !  |                        =102: AS poloidal Br dynamo term           |
  !  |                        =103: AS poloidal Br diffusion             |
  !  |                        =104: AS toroidal Bp production            |
  !  |                        =105: AS toroidal Bp dynamo term           |
  !  |                        =106: AS toroidal Bp omega effect          |
  !  |                        =107: AS toroidal Bp diffusion             |
  !  |     n_movie_field_start(n_field,n_movie) = defines where first    |
  !  |         element of a field is stored in frames(*)                 |
  !  |     n_movie_field_stop(n_field,n_movie) = defines where last      |
  !  |         element of a field is stored in frames(*)                 |
  !  |                 l_movie_oc : field for OC stored ?                |
  !  |                 l_movie_ic : fields for IC stored (only B-field)  |
  !  |                                                                   |
  !  |     The subroutine also defines appropriate file names for        |
  !  |     the movie files. These generally have the form:               |
  !  |                 TYPE_mov.TAG                                      |
  !  |     The TYPE corresponts to the input movie type,                 |
  !  |     for TOROIDAL AXISYMMETRIC fields we use TYPE=TAS.             |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE horizontal_data
  USE logic
  USE movie_data
  USE output_data
  USE charmanip, only: capitalize,delete_string, &
       str2dble,length_to_blank
  USE parallel_mod,only:rank
  IMPLICIT NONE

  !--- Local:
  logical :: lEquator
  integer :: length,length_fn,lengthC
  character(len=80) :: string,word,stringC
  character(len=80) :: file_name
  character(len=50) :: typeStr
  integer :: n_r,n_theta,n_phi
  real(kind=8) :: r_movie,theta_movie,phi_movie
  real(kind=8) :: phi_max
  real(kind=8) :: rad
  real(kind=8) :: const
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


  !--- End of Declaration
  !-----------------------------------------------------------------------


  !--- Initialize first storage index:
  n_field_start=1


  !--- Converts from radiant to degree:
  rad=45.d0/datan(1.d0)

  !--- Loop over max possible no of movies:

  l_movie_oc   =.false.
  l_movie_ic   =.false.
  l_HTmovie    =.false.
  l_dtBmovie   =.false.
  l_store_frame=.false.

  n_rc=0
  n_tc=0
  n_pc=0

  DO i=1,n_movies_max

     lStore=.true.
     lIC   =.false.

     string=movie(i)

     IF ( LEN(TRIM(string)) .EQ. 0 ) CYCLE !blank string

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
        endif
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
     else if ( index(string,'BP' ) /= 0 .AND. &
          index(string,'TOR') == 0 .AND. &
          index(string,'AX' ) == 0 .AND. &
          index(string,'AB' ) == 0 ) then
        n_type=3
        typeStr=' azimuthal magnetic field '
        file_name='Bp_'
        lIC=.true.
        n_fields=1
        n_field_type(1)=3
     else if ( index(string,'BH') /= 0 .OR. &
          index(string,'BM') /= 0 ) then
        n_type=4 ! Horizontal field
        file_name='Bh_'
        lIC=.true.
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
     else if ( index(string,'FIELDLINE') /= 0 .OR. &
          index(string,'FL') /= 0 ) then
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
     else if ( index (string,'LORENTZ') /= 0 .OR. &
          index(string,'LF') /= 0 ) then
        n_type=99
        typeStr=' phi Lorentz force '
        file_name='LFp_'
        lIC=.FALSE.
        n_fields=1
        n_field_type(1)=54
     else if ( ( index(string,'AX') /= 0 .AND.  &
          index(string,'BP') /= 0 ) .OR. &
          index(string,'AB') /= 0 ) then
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
        endif
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
        endif
     else if ( index(string,'BP') /= 0 .AND. &
          index(string,'TOR') /= 0 ) then
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
        endif
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
           lStore=.false.
           n_fields=1
           n_field_type(1)=14
        endif
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
        IF ( index(string,'AX') /= 0 ) THEN
           n_type=18
           typeStr=' axisym. azimuthal velocity field '
           file_name='AV_'
           n_fields=1
           n_field_type(1)=11
        ELSE
           n_type=13
           typeStr=' phi comp. of velocity field '
           file_name='Vp_'
           n_fields=1
           n_field_type(1)=6
        ENDIF
     else if ( index(string,'VH') /= 0 .OR. &
          index(string,'VM') /= 0 ) then
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
     else if ( index(string,'STREAMLINE') /= 0 .OR. &
          index(string,'SL') /= 0 ) then
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
     else if ( index(string,'VZ') /= 0 ) then
        n_type=94
        typeStr=' z-component of velocity '
        file_name='VZ_'
        n_fields=1
        n_field_type(1)=15
     else if ( ( index(string,'AX'  ) /= 0 .AND.  &
          index(string,'HEL' ) /= 0 ) .OR. &
          index(string,'AHEL') /= 0 ) then
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
     else if ( index(string,'AX' ) /= 0 .AND. &
          ( index(string,'TEM') /= 0 .OR.  &
          index(string,'ENT') /= 0 ) ) then
        n_type=28
        typeStr=' axisymmetric temp. '
        file_name='AT'
        n_fields=1
        n_field_type(1)=12
     else if ( index(string,'ENTROPY') /= 0 .OR. &
          index(string,'TEM') /= 0 ) then
        ns=index(string,'S')
        if ( string(ns:ns+2) == 'SUR' ) then
           write(*,*) '! No surface T field available !'
           stop
        end if
        n_type=21
        typeStr=' temperature field '
        file_name='T_'
        n_fields=1
        n_field_type(1)=7
     else if ( ( index(string,'CONV' ) /= 0 .AND.  &
          index(string,'HEAT' ) /= 0 ) .OR. &
          index(string,'HEATT') /= 0 ) then
        ns=index(string,'S')
        if ( string(ns:ns+2) == 'SUR' ) then
           write(*,*) '! No surface T field available !'
           stop
        end if
        n_type=22
        typeStr=' radial convective heat transport '
        file_name='HT_'
        l_HTmovie=.true.
        n_fields=1
        n_field_type(1)=17
     else if ( index(string,'AX' ) /= 0 .AND. &
               index(string,'HEATF') /= 0   ) then
        n_type=113
        typeStr='axisymmetric dSdr '
        file_name='AHF_'
        l_HTmovie=.true.
        n_fields=1
        n_field_type(1)=92
     else if ( index(string,'HEATF') /= 0 ) THEN
        ns=index(string,'S')
        if ( string(ns:ns+2) == 'SUR' ) then
           write(*,*) '! No surface T field available !'
           stop
        end if
        n_type=110
        typeStr=' radial heat transport '
        file_name='HF_'
        l_HTmovie=.true.
        n_fields=1
        n_field_type(1)=91
     else if ( index(string,'POT') /= 0 .AND. &
          index(string,'TOR') /= 0 ) then
        n_type=96
        typeStr=' radial convective heat transport '
        file_name='PotTor_'
        lIC=.true.
        lStore=.false.
        n_fields=1
        n_field_type(1)=51
     else
        write(*,*) 'Couldnt interpret movie field from'
        write(*,*) '! string:',string
        stop
     end if


     !--- Identify surface type:

     length_fn=len(trim(file_name))
     IF ( n_type == 103 ) THEN
        n_surface=1 !
        n_const=1   !
        n_field_size=n_phi_max*n_theta_max
        n_field_size_ic=n_field_size
        const=r_cmb
     else if (   index(string,'AX') /= 0 .OR. &
          file_name(1:2) == 'AV' .OR. &
          file_name(1:2) == 'AB' .OR. &
          n_type == 30 .OR. n_type == 31 .OR. &
          n_type == 32 .OR. n_type == 33 .OR. &
          n_type == 40 .OR. n_type == 41 .OR. &
          n_type == 42 .OR. n_type == 43 .OR. &
          n_type == 17 ) then
        !--- Axisymmetric stuff:
        n_surface=3  ! PHI=const.
        n_const=1
        n_field_size=n_r_max*n_theta_max
        n_field_size_ic=n_r_ic_max*n_theta_max
        const=0.d0
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
        const=1.d0
     else if ( index(string,'R='     ) /= 0 .OR. &
          index(string,'RAD='   ) /= 0 .OR. &
          index(string,'RADIUS=') /= 0 ) then
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
        IF ( r_movie == 0 ) THEN
           n_const=n_r_icb
           const  =r_icb
        ELSE IF ( r_movie == 1 ) THEN
           n_const=n_r_cmb
           const  =r_cmb
        ELSE IF ( r_movie > 0 ) THEN
           r_movie=r_icb+r_movie
           do n_r=2,n_r_max
              if ( r(n_r-1) > r_movie .AND. &
                   r(n_r)  <= r_movie ) then
                 if ( r(n_r-1)-r_movie < &
                      r_movie-r(n_r) ) then
                    n_const=n_r-1
                 else
                    n_const=n_r
                 end if
                 const=r(n_const)
                 exit
              end if
           end do
        ELSE
           !------ Negative numbers signify inner core values in
           !       fractions of r_icb:
           r_movie=-r_movie*r_icb
           do n_r=2,n_r_ic_max
              if ( r_ic(n_r-1) >= r_movie .AND. &
                   r_ic(n_r)  <= r_movie ) then
                 if ( r_ic(n_r-1)-r_movie < &
                      r_movie-r_ic(n_r) ) then
                    n_const=-(n_r-1)
                 else
                    n_const=-n_r
                 end if
                 const=r_ic(-n_const)
                 exit
              end if
           end do
        END IF

     ELSE IF ( INDEX(string,'EQ') /= 0 .OR. lEquator ) THEN

        n_surface=2    ! Equator
        n_const=n_theta_max
        file_name=file_name(1:length_fn)//'EQU_'
        n_field_size=n_r_max*n_phi_max
        n_field_size_ic=n_r_ic_max*n_phi_max
        const=rad*theta(n_const)

     else if ( index(string,'T='    ) /= 0 .OR. &
          index(string,'THETA=') /= 0 ) then

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
        theta_movie=dabs(theta_movie)
        theta_movie=theta_movie/rad

        !------ Choose closest colatitude grid point:
        foundGridPoint=.FALSE.
        do n_theta=1,n_theta_max-1
           if ( theta(n_theta)  <= theta_movie .AND. &
                theta(n_theta+1) >= theta_movie ) then
              if ( theta(n_theta+1)-theta_movie < &
                   theta_movie-theta(n_theta) ) then
                 n_const=n_theta+1
              else
                 n_const=n_theta
              end if
              foundGridPoint=.TRUE.
              exit
           end if
        end do
        if ( .not. foundGridPoint ) then
           if ( theta_movie-theta(n_theta_max) <= &
                theta(1)+180.d0/rad-theta_movie ) then
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

     else if ( index(string,'MER' ) /= 0 .OR.  &
          index(string,'P='  ) /= 0  .OR. &
          index(string,'PHI=') /= 0 ) then

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
        if ( phi_movie < 0.d0 ) phi_movie=360.d0-phi_movie
        phi_max=360.d0/minc
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
        foundGridPoint=.FALSE.
        do n_phi=1,n_phi_max-1
           if ( phi(n_phi)  <= phi_movie .AND. &
                phi(n_phi+1) >= phi_movie ) then
              if ( phi(n_phi+1)-phi_movie < &
                   phi_movie-phi(n_phi) ) then
                 n_const=n_phi+1
              else
                 n_const=n_phi
              end if
              foundGridPoint=.TRUE.
              exit
           end if
        end do
        if ( .not. foundGridPoint ) then
           if ( phi_movie-phi(n_phi_max) <= &
                phi(1)+phi_max-phi_movie ) then
              n_const=n_phi_max
           else
              n_const=1
           end if
        end if

        const=rad*phi(n_const)

     ELSE
        write(*,*) 'Couldnt interpret movie surface from'
        write(*,*) '! string:',string
        write(*,*) '! file name:',file_name
        stop
     END IF

     IF ( n_field_type(1) == 54 .AND.                &
          ( ( n_surface /= 0 .AND. n_surface /= 3 ) .OR. &
          index(string,'AX') /= 0 ) ) THEN
        WRITE(*,*) 'Sorry, can only prepare movie file for'
        WRITE(*,*) 'phi component of LF in 3d or for phi-cut.'
        STOP
     END IF

     !--- Now store the necessary information:

     !------ Increase number of movies:
     n_movies=n_movies+1
     lStoreMov(n_movies)=lStore
     lICField(n_movies)=lIC
     if ( .NOT. lStore .AND. n_field_type(1) /= 13 .AND. &
          n_field_type(1) /= 14 .AND. &
          n_field_type(1) /= 30 .AND. &
          n_field_type(1) /= 42 .AND. &
          n_field_type(1) /= 50 .AND. &
          n_field_type(1) /= 51 .AND. &
          n_field_type(1) /= 52 .AND. &
          n_field_type(1) /= 54 ) l_dtBmovie= .TRUE. 

     !------ Translate horizontal movies:
     IF ( n_type == 4 ) THEN
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
     END IF

     !--- Inner core and/or outer core field?
     n_fields_oc=n_fields
     n_fields_ic=0
     if ( lIC ) then
        if ( n_surface == 1 .AND. n_const < 0 ) then
           n_fields_oc=0
           n_fields_ic=n_fields
        else if ( n_surface == 0 .OR. n_surface == 2 .OR. &
             n_surface == 3 ) then
           n_fields_ic=n_fields
        end if
     end if
     if ( n_fields_oc > 0 ) l_movie_oc= .TRUE. 
     if ( n_fields_ic > 0 ) l_movie_ic= .TRUE. 

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
     DO n=1,n_fields
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
     END DO

     !--- Write info about output files into log-file:

     IF (rank.EQ.0) THEN
        IF ( l_save_out ) THEN
           OPEN(n_log_file,file=log_file,status='unknown', &
                POSITION='APPEND')
        END IF
        WRITE(n_log_file,                              &
             &     '(/,'' ! PRODUCING MOVIE-FILE :'',a64)') &
             movie_file(n_movies)
        WRITE(n_log_file,*) '!    Containing ',typeStr

        IF ( n_surface == -1 ) THEN
           WRITE(n_log_file,*) &
                &   '!    at the surface !'
        ELSE IF ( n_surface == 0 ) THEN
           WRITE(n_log_file,*) &
                &   '!    in 3d !'
        ELSE IF ( n_surface == 1 .AND. n_const == 1 ) THEN
           WRITE(n_log_file,*) &
                &   '!    at the core-mantle boundary !'
        ELSE IF ( n_surface == 1 .AND. n_const > 1 ) THEN
           WRITE(n_log_file, &
                &   '('' !    at r='',f10.6)') r(n_const)
        ELSE IF ( n_surface == 1 .AND. n_const < 0 ) THEN
           WRITE(n_log_file, &
                &   '('' !    at r='',f10.6)') r_ic(n_const)
        ELSE IF ( n_surface == 2 ) THEN
           WRITE(n_log_file, &
                &   '('' !    at theta='',f12.6)') rad*theta(n_const)
        ELSE IF ( n_surface == 3 ) THEN
           WRITE(n_log_file, &
                &   '('' !    at phi='',f12.6)') rad*phi(n_const)
        END IF
        IF ( n_fields_ic > 0 ) &
             WRITE(n_log_file, &
             &   '('' !    including inner core magnetic field.'')')

        IF ( l_save_out ) CLOSE(n_log_file)
     END IF

  END DO     ! loop over all movies


  !--- Test whether movie frame will fit into work arrays:
  if ( n_field_start > n_frame_work ) then
     write(*,*) '! Storage size of frames(*) does not suffice'
     write(*,*) '! to store all frames for one time step.'
     if ( lMovieMem == 0 ) then
        write(*,*) &
             &  '! For additional memory set lMovieMem=1 in truncation.f!'
     else
        write(*,*) '! Increase n_movie_work in truncation.f!'
     end if
     stop
  end if


  return
end subroutine get_movie_type
!----------------------------------------------------------------------------
