.. _secOutNmlStd:

Standard time-series outputs
----------------------------

The **log** outputs controls the output of all the default time series of the
file: kinetic and magnetic energies (:ref:`e_kin.TAG <secEkinFile>`,
:ref:`e_mag_oc.TAG <secEmagocFile>` and :ref:`e_mag_ic.TAG <secEmagicFile>`
files), dipole information (:ref:`dipole.TAG <secDipoleFile>` file), rotation
(:ref:`rot.TAG <secRotFile>`) parameters (:ref:`par.TAG <secParFile>`) and
various additional diagnostics (:ref:`heat.TAG <secHeatFile>`):

.. _varn_log_step:

* **n_log_step** (default :f:var:`n_log_step=50 <n_log_step>`) is an integer. This is the number of timesteps between two log outputs.

  .. warning:: Be careful: when using too small :f:var:`n_log_step`, the disk access will dramatically increases, thus decreasing the code performance.

.. _varn_logs:

* **n_logs** (default :f:var:`n_logs=0 <n_logs>`) is an integer. This is the number of log-information sets to be written.

.. _vart_log:

* **t_log**  (default  :f:var:`t_log=-1.0 -1.0 ... <t_log>`) is real array, which contains the times when log outputs are requested.

.. _vardt_log:

* **dt_log** (default :f:var:`dt_log=0.0 <dt_log>`) is a real, which defines the time interval between log outputs.

.. _vart_log_start:

* **t_log_start** (default :f:var:`t_log_start=0.0 <t_log_start>`) is a real, which defines the time to start writing log outputs.

.. _vart_log_stop:

* **t_log_stop** (default :f:var:`t_log_stop=0.0 <t_log_stop>`) is a real, which defines the time to stop writing log outputs.

.. _secRstIn:

Restart files
-------------

The **rst** outputs controls the output of restart files (:ref:`checkpoint_t_#.TAG <secRestartFile>`) (i.e. check points in time from which the code could be restarted):

.. _varn_rst_step:

* **n_rst_step** (default :f:var:`n_rst_step=0 <n_rst_step>`) is an integer. This is the number of timesteps between two restart files.

.. _varn_rsts:

* **n_rsts** (default :f:var:`n_rsts=1 <n_rsts>`) is an integer. This is the number of restart files to be written.

.. _vart_rst:

* **t_rst**  (default  :f:var:`t_rst=-1.0 -1.0 ... <t_rst>`) is real array, which contains the times when restart files are requested.

.. _vardt_rst:

* **dt_rst** (default :f:var:`dt_rst=0.0 <dt_rst>`) is a real, which defines the time interval between restart files.


* **t_rst_start** (default :f:var:`t_rst_start=0.0 <t_rst_start>`) is a real, which defines the time to start writing restart files.


* **t_rst_stop** (default :f:var:`t_rst_stop=0.0 <t_rst_stop>`) is a real, which defines the time to stop writing restart files.


* **n_stores** (default :f:var:`n_stores=0 <n_stores>`) is an integer. This is another way of requesting a certain number of restart files. However, instead of creating each time a new restart file, if ``n_stores > n_rsts``  the restart file is overwritten, which can possibly help saving some disk space.

.. warning:: The ``rst`` files can become quite big and writting them too frequently will slow down the code. Except for very special use, the default set up should be sufficient.

.. _secOutGraphFile:

Graphic files
-------------

The **graph** outputs controls the output of graphic files (:ref:`G_#.TAG <secGraphFile>`) which contain a snapshot the entropy, the velocity field and the magnetic fields:

.. _varn_graph_step:

* **n_graph_step** (default :f:var:`n_graph_step=0 <n_graph_step>`) is an integer. This is the number of timesteps between two graphic files.

.. _varn_graphs:

* **n_graphs** (default :f:var:`n_graphs=1 <n_graphs>`) is an integer. This is the number of graphic files to be written.

.. _vart_graph:

* **t_graph**  (default  :f:var:`t_graph=-1.0 -1.0 ... <t_graph>`) is real array, which contains the times when graphic files are requested.

.. _vardt_graph:

* **dt_graph** (default :f:var:`dt_graph=0.0 <dt_graph>`) is a real, which defines the time interval between graphic files.

.. _vart_graph_start:

* **t_graph_start** (default :f:var:`t_graph_start=0.0 <t_graph_start>`) is a real, which defines the time to start writing graphic files.

.. _vart_graph_stop:

* **t_graph_stop** (default :f:var:`t_graph_stop=0.0 <t_graph_stop>`) is a real, which defines the time to stop writing graphic files.



Spectra
-------

The **spec** outputs controls the output of spectra: kinetic energy spectra (:ref:`kin_spec_#.TAG <secKinSpecFile>`), magnetic energy spectra (:ref:`mag_spec_#.TAG <secMagSpecFile>`) and thermal spectra (:ref:`T_spec_#.TAG <secTSpecFile>`):

.. _varn_spec_step:

* **n_spec_step** (default :f:var:`n_spec_step=0 <n_spec_step>`) is an integer. This is the number of timesteps between two spectra.

.. _varn_specs:

* **n_specs** (default :f:var:`n_specs=0 <n_specs>`) is an integer. This is the number of spectra to be written.

.. _vart_spec:

* **t_spec**  (default  :f:var:`t_spec=-1.0 -1.0 ... <t_spec>`) is real array, which contains the times when spectra are requested.

.. _vardt_spec:

* **dt_spec** (default :f:var:`dt_spec=0.0 <dt_spec>`) is a real, which defines the time interval between spectra.

.. _vart_spec_start:

* **t_spec_start** (default :f:var:`t_spec_start=0.0 <t_spec_start>`) is a real, which defines the time to start writing spectra.

.. _vart_spec_stop:

* **t_spec_stop** (default :f:var:`t_spec_stop=0.0 <t_spec_stop>`) is a real, which defines the time to stop writing spectra.

.. _varl_2D_spectra:

* **l_2D_spectra** (default :f:var:`l_2D_spectra=.false. <l_2d_spectra>`) is a 
  logical. When set to ``.true.``, this logical enables the calculation of 2-D
  spectra in the :math:`(r,\ell)` and in the :math:`(r,m)` parameter spaces. 
  Those data are stored in the files named :ref:`2D_[mag|kin]_spec_#.TAG <sec2DSpectra>`.


Movie files
-----------

The **movie** outputs controls the output of movie files (:ref:`*_mov.TAG <secMovieFile>`). 


Specific inputs
+++++++++++++++

.. _varl_movie:

* **l_movie** (default :f:var:`l_movie=.false. <l_movie>`) is a logical. It needs to be turned on to get movie computed.

  Several movie-files can be produced during a run (it is now limited to 30 by
  the variable ``n_movies_max`` in the module :f:mod:`movie`). The movies are
  defined by a keyword determining the fields to be plotted and an expression
  that determines the nature of movie (:math:`r`-slice, :math:`\theta`-slice,
  :math:`\phi`-slice, etc.). The code searches this information in a
  character string  provided for each movie.  These strings are elements of the 
  array :ref:`movie <varmovie>`:

.. _varmovie:

* **movie** (default :f:var:`movie=' ', ' ', ... <movie>`) is a character string array. It contains the description of the movies one wants to compute.

  For example, to invoke a movie(file) that shows (stores) the radial magnetic
  component of the magnetic field at the CMB, you have to provide the line

    .. code-block:: fortran

        movie(1)="Br CMB",

  in the :ref:`&output <secOutputNml>` namelist. Here, ``Br`` is the keyword for 
  the radial component of the magnetic field and ``CMB`` is the expression that
  defines the movie surface. If, in addition, a movie of the temperature field 
  at the meridional slice ``phi=0`` and a movie of the :math:`z`-vorticity in 
  the equatorial plane are desired, the following line have to be added:

     .. code-block:: fortran

        movie(2)="Temp phi=0",
        movie(3)="Vortz eq",

  Note that the code does **not interpret spaces and ignores additional characters**
  that do not form a keyword or a surface definition. Thus, for example ``Br`` or ``B r``
  or ``Bradial`` are all interpreted as the same keyword. Furthermore, the
  interpretation is **not case-sensitive**. The following table gives the possible
  keywords for movie calculations and their corresponding physical meaning:


  .. tabularcolumns:: |p{3cm}|p{10cm}|

  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Keyword                                         | Fields stored in movie file                                                                                                                     |
  +=================================================+=================================================================================================================================================+
  | Br[radial]                                      | Radial component of the magnetic field :math:`B_r`.                                                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bt[heta]                                        | Latitudinal component of the magnetic field  :math:`B_\theta`.                                                                                  |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bp[hi]                                          | Azimuthal component of the magnetic field  :math:`B_\phi`.                                                                                      |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bh[orizontal]                                   | The two horizontal components of the magnetic field.                                                                                            |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bs                                              | Cylindrically radial component of the magnetic field :math:`B_s`.                                                                               |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Ba[ll]                                          | All magnetic field components.                                                                                                                  |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Fieldline[s] or FL                              | Axisymmetric poloidal field lines in a meridional cut.                                                                                          |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | AX[ISYMMETRIC] B or AB                          | Axisymmetric phi component of the magnetic field for :math:`\phi=cst.`                                                                          |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vr[adial]                                       | Radial component of the velocity field :math:`u_r`.                                                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vt[heta]                                        | Latitudinal component of the velocity field  :math:`u_\theta`.                                                                                  |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vp[hi]                                          | Azimuthal component of the velocity field  :math:`u_\phi`.                                                                                      |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vh[orizontal]                                   | Horizontal velocity field, two components depending on  the surface.                                                                            |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Va[ll]                                          | All velocity field components.                                                                                                                  |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Streamline[s] or SL                             | Field lines of axisymmetric poloidal field for :math:`\phi=cst.`                                                                                |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | AX[ISYMMETRIC] V or AV                          | Axisymmetric component of the velocity field for :math:`\phi=cst.`                                                                              |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vz                                              | Vertical component of the velocity :math:`u_z`.                                                                                                 |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vs                                              | Cylindrical radil component of the velocity :math:`u_s`.                                                                                        |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Voz                                             | Vertical component of the vorticity :math:`\omega_z`.                                                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vor                                             | Radial component of the vorticity  :math:`\omega_r`.                                                                                            |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vop                                             | Azimuthal component of vorticity  :math:`\omega_\phi`                                                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Tem[perature] or Entropy                        | Temperature/Entropy                                                                                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Entropy (or Tem[perature]) AX[ISYMMETRIC] or AT | Axisymmetric temperature/entropy field for :math:`\phi=cst.`                                                                                    |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Heat t[ransport]                                | Radial advection of temperature :math:`u_r\frac{\partial s}{\partial r}`                                                                        |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | HEATF AX[iSYMMETRIC]                            | Conducting heat flux :math:`\partial s /\partial r`                                                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Voz geos                                        | Vertical component of the vorticity :math:`\omega_z` averaged over the rotation axis.                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vs geos                                         | Cylindrical radial component of the velocity :math:`u_s` averaged over the rotation axis.                                                       |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Vp geos                                         | Azimuthal component of the velocity :math:`u_\phi` averaged over the rotation axis.                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | FL Pro                                          | Axisymmetric field line stretching.                                                                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | FL Adv                                          | Axisymmetric field line advection.                                                                                                              |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | FL Dif                                          | Axisymmetric field line diffusion.                                                                                                              |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | AB Pro                                          | Toroidal axisymmetric  field production.                                                                                                        |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | AB Dif                                          | Toroidal axisymmetric field diffusion.                                                                                                          |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Br Pro                                          | Production of radial magnetic field  :math:`B_r`.                                                                                               |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Br Adv                                          | Advection of radial magnetic field  :math:`B_r`.                                                                                                |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Br Dif                                          | Diffusion of radial magnetic field :math:`B_r`.                                                                                                 |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jr                                              | Radial component of the current :math:`j_r`.                                                                                                    |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jr Pro                                          | Production of radial current + :math:`\Omega`-effect.                                                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jr Adv                                          | Advection of the radial component of the current :math:`j_r`.                                                                                   |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jr Dif                                          | Diffusion of the radial component of the current :math:`j_r`.                                                                                   |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bz Pol                                          | Poloidal part of vertical component of the magnetic field  :math:`B_z`.                                                                         |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bz Pol Pro                                      | Production of the poloidal part of the vertical component of the magnetic field  :math:`B_z`.                                                   |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bz Pol Adv                                      | Advection  of the poloidal part of the vertical component of the magnetic field :math:`B_z`.                                                    |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bz Pol Dif                                      | Diffusion of the poloidal part of the vertical component of the magnetic field :math:`B_z`.                                                     |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jz Tor                                          | Toroidal part of the vertical component of the current (:math:`j_z`).                                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jz Tor Pro                                      | Production of the toroidal part of the vertical component of the current :math:`j_z`.                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jz Tor Adv                                      | Advection  of the toroidal part of the vertical component of the current :math:`j_z`.                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Jz Tor Dif                                      | Diffusion of the  toroidal part of the vertical component of the current :math:`j_z`.                                                           |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bp Tor                                          | Toroidal part of the azimuthal component of the magnetic field :math:`B_\phi`.                                                                  |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bp Tor Pro                                      | Production of the toroidal part of the azimuthal component of the magnetic field :math:`B_\phi`.                                                |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bp Tor Adv                                      | Advection of the toroidal part of the azimuthal component of the magnetic field :math:`B_\phi`.                                                 |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bp Tor Dif                                      | Diffusion of the toroidal part of the azimuthal component of the magnetic field :math:`B_\phi`.                                                 |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | HEL[ICITY]                                      | Kinetic helicity :math:`{\cal H}=\vec{u}\cdot(\vec{\nabla}\times\vec{u})`                                                                       |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | AX[ISYMMETRIC HELICITY] or AHEL                 | Axisymmetric component of the kinetic helicity.                                                                                                 |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Bt Tor                                          | Toroidal component of the latitudinal component of the magnetic field :math:`B_\theta`.                                                         |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Pot Tor                                         | Toroidal potential.                                                                                                                             |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Pol Fieldlines                                  | Poloidal fieldlines.                                                                                                                            |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Br Shear                                        | Azimuthal shear of the radial component of the magnetic field :math:`B_r`                                                                       |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Lorentz[force] or LF                            | Lorentz force (only :math:`\phi`-component).                                                                                                    |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
  | Br Inv                                          | Inverse field apperance at CMB.                                                                                                                 |
  +-------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------+

  The following table gives the possible surface expression for movie calculations 
  and their corresponding physical meaning:


  +--------------------+-------------------------------------------------+
  | Surface expression | Definition                                      |
  +====================+=================================================+
  | CMB                | Core-mantle boundary                            |
  +--------------------+-------------------------------------------------+
  | Surface            | Earth surface                                   |
  +--------------------+-------------------------------------------------+
  | EQ[uatot]          | Equatorial plane                                |
  +--------------------+-------------------------------------------------+
  | r=radius           | Radial cut at r=radius with radius given in     |
  |                    | units of the outer core radius.                 |
  +--------------------+-------------------------------------------------+
  | theta=colat        | Latitudinal cut at theta=colat given in degrees |
  +--------------------+-------------------------------------------------+
  | phi=phiSlice       | Azimuthal cut ath phi=phiSlice given in degrees.|
  +--------------------+-------------------------------------------------+
  | AX[isymmetric]     | Axisymmetric quantity in an azimuthal plane     |
  +--------------------+-------------------------------------------------+
  | 3D                 | 3D array                                        |
  +--------------------+-------------------------------------------------+


  Here is an additional example of the possible combinations to build your
  desired ``movie`` files.

  .. code-block:: fortran

     l_movie  = .true.,
     movie(1) = "Br CMB", 
     movie(2) = "Vr EQ",
     movie(3) = "Vortr r=0.8",
     movie(4) = "Bp theta=45",
     movie(5) = "Vp phi=10",
     movie(6) = "entropy AX",
     movie(7) = "vr 3D",

.. _secMovieStdin:
  
Standard inputs
+++++++++++++++

.. _varn_movie_step:

* **n_movie_step** (default :f:var:`n_movie_step=0 <n_movie_step>`) is an integer. This is the number of timesteps between two movie outputs.

* **n_movies** (default :f:var:`n_movies=1 <n_movies>`) is an integer. This is the number of movie outputs to be written.

.. _vart_movie:

* **t_movie**  (default  :f:var:`t_movie=-1.0 -1.0 ... <t_movie>`) is real array, which contains the times when movie outputs are requested.

.. _vardt_movie:

* **dt_movie** (default :f:var:`dt_movie=0.0 <dt_movie>`) is a real, which defines the time interval between movie outputs.


* **t_movie_start** (default :f:var:`t_movie_start=0.0 <t_movie_start>`) is a real, which defines the time to start writing movie outputs.


* **t_movie_stop** (default :f:var:`t_movie_stop=0.0 <t_movie_stop>`) is a real, which defines the time to stop writing movie outputs.


Field Averages
--------------

The code can perform on-the-fly time-averaging of entropy, velocity field and magnetic field. Respective graphic output and spectra are written into the corresponding files (with :ref:`G_ave.TAG <secGraphFile>`, :ref:`kin_spec_ave.TAG <secKinSpecAveFile>`,  :ref:`mag_spec_ave.TAG <secMagSpecAveFile>`). The time-averaged energies are written into the :ref:`log.TAG <secLogFile>` file.


.. _varl_average:

* **l_average** (default :f:var:`l_average=.false. <l_average>`) is a logical, which enables the time-averaging of fields when set to ``.true.``.

  .. warning:: Time-averaging has a large memory imprint as it requires the storage of 3-D arrays. Be careful, when using large truncations.

.. _varl_spec_avg:

* **l_spec_avg** (default :f:var:`l_spec_avg=.false. <l_spec_avg>`) is a logical, which enables the time-averaging of spectra when set to ``.true.``. It is always set to ``.true.``, if :ref:`l_average=.true. <varl_average>`.
