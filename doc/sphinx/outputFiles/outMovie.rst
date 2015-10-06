.. _secMovieFile:

Movie files ``*_mov.TAG``
=========================

.. note:: These files are written **only** when :ref:`l_movie = .true.
 <varl_movie>` or when a finite number of movie frames are asked for using the
 input parameters described in the :ref:`standard inputs section
 <secMovieStdin>` of the :ref:`output control namelist <secOutputNml>`.

These are unformatted fortran files containing time evolution of fields on
different surfaces - constant radius, colatitude or azimuth or on the full 3D
grid. The fields can be of various types like radial magnetic field or
velocity, entropy, helicity etc. The type of field and the type of surface can
be specified using a string that begins with the field name, followed by the
surface type (or 'full 3D', when a 3D movie is desired). One such example is as
follows:

    .. code-block:: fortran
     
     l_movie = .true.,
     n_movie_frames = 1000,
     movie(1) = "B r r=0.5",
     movie(2) = "V all 3D",
     movie(3) = "Hel Eq"

The code does not interpret any whitespaces and is not case-sensitive so
there's no difference between, say, ``B r cmb`` and ``brcmb``. For further
details and a list of keywords for different fields and surfaces, please refer
to the :ref:`movie <varmovie>` in the :ref:`output control namelist <secOutputNml>`.

These files are written by the subroutine :f:subr:`write_movie_frame <out_movie/write_movie_frame>`.

The movie files are suitably named to reflect the type of field and surface. Their names begin
with the keyword for the type of movie asked for, followed by the type of surface, followed by the word 'mov'.
Thus, a generic movie name looks like:
      
      ``Keyword_SurType_mov.TAG``

E.g: if one asks for the radial component of magnetic field on surface of CMB, the movie would be named as ``Br_CMB_mov.TAG``.

When asks multiple movies for same surface types but different surface levels, the surfaces are numbered with integers. Thus, for the following namelist input,

    .. code-block:: fortran

     l_movie = .true.,
     n_movie_frames = 1000,
     movie(1) = "B r r=0.5",
     movie(2) = "V p r=0.5",
     movie(3) = "V r r=0.8",

one would get the following movie files as output:
    
    .. code-block:: fortran

     Br_R=C1_mov.TAG
     Vp_R=C1_mov.TAG
     Vr_R=C2_mov.TAG

The structure of a generic movie file is as follows:

    .. code-block:: fortran

     !-------
     ! Line 1
     !-------

     version                       !Movie version: 'JW_Movie_Version_2'

     !-------
     ! Line 2
     !-------

     n_type, n_surface,            !Type of movie,
     const, n_fields               !Type of surface (r,theta,phi,CMB,Eq etc.)
    
     !-------
     ! Line 3
     !-------

     n_movie_field_type(1:n_fields, n_movie) !Type of fields (velocity,
                                             !mag field, vorticity etc.)
     !-------
     ! Line 4
     !-------
     
     runid
     
     !-------
     ! Line 5
     !-------
 
     n_r_mov_tot, n_r_max,          !Total number of
     n_theta_max, n_phi_max,        !radial grid points (including IC),
     minc, ra, ek, pr, prmag,       !grid data, physical parameters 
     radratio, tScale

     !-------
     ! Line 6
     !-------
 
     r_mov_tot(1:n_r_mov_tot)/r_cmb !All radii in terms of r_CMB
     
     !-------
     ! Line 7
     !-------
 
     theta(1:n_theta_max)           !All theta points
     
     !-------
     ! Line 8
     !-------
 
     phi(1:n_phi_max)               !All phi points

     !------------------------------------------------------------------- 
     
     !---------
     ! Frame N
     !---------

     !-----------
     ! Line 8 + N
     !-----------

     n_frame, t_movie(N), omega_ic, omega_ma, dipLat, dipLon, dipStr, dipStrGeo

     !---------------
     ! Line 8 + (N+1)
     !---------------

     frame_data(1:n_fields,n_start:n_stop)  !Desired field data on a
                                            !surface or 3D volume
                                            !n_start = start index of a field
                                            !n_stop  = last index of a field
     
     !-----------
     ! Frame N+1
     !-----------

     !---------------
     ! Line 8 + (N+2)
     !---------------

     n_frame, t_movie(N+1), omega_ic, omega_ma, dipLat, dipLon, dipStr, dipStrGeo

     !---------------
     ! Line 8 + (N+3)
     !---------------

     frame_data(1:n_fields,n_start:n_stop)  !Desired field data on a
                                            !surface or 3D volume
                                            !n_start = start index of a field
                                            !n_stop  = last index of a field

     ...
      
     !-----------
     ! Frame N+M                            !M is the desired number of movie frames
     !-----------

     !---------------
     ! Line 8 + (N+M)
     !---------------

     n_frame, t_movie(N+M), omega_ic, omega_ma, dipLat, dipLon, dipStr, dipStrGeo

     !---------------
     ! Line 8 + (N+M)
     !---------------

     frame_data(1:n_fields,n_start:n_stop)  !Desired field data on a
                                            !surface or 3D volume
                                            !n_start = start index of a field
                                            !n_stop  = last index of a field


The 2D movie files can be read and displayed using the python class :py:class:`Movie <magic.Movie>` as follows:

    >>> Movie()   #Lists out available movie files to choose from
    >>> M = Movie(file = 'Vr_R=C1_mov.TAG')

The 3D movie files can be read using the python class :py:class:`Movie3D <magic.Movie3D>`:
    
    >>> M = Movie3D(file = 'V_3D_mov.TAG')
