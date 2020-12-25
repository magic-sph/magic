
.. _secTOoutputFiles:

``TO`` outputs
==============

.. note:: These output files are **only** written when :ref:`l_TO=.true. <varl_TO>`

.. _secTayFile:

``Tay.TAG``
------------

This file contains the time series of the Taylorisation as well as some measures
of the relative geostrophic energy.

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with
the following options:

    >>> # To load the most recent 'Tay.TAG' file in a directory
    >>> ts = MagicTs(tag='Tay')


``TOnhs.TAG`` and ``TOshs.TAG``
--------------------------------

Those files correspond to the z-averaging of the axisymmetric phi component of the
Navier-Stokes equations. It contains the different cylindrical profiles of the
forces involved the zonal equation as well as some additional measures of the
Taylorization of the solution. shs corresponds to Southern Hemisphere (inside the
tangent cylinder), while nhs corresponds to Nothern Hemisphere).

Those files can be read using :py:class:`MagicTOHemi <magic.MagicTOHemi>` with
the following options:

    >>> # To load 'TOshs.test' and plot the time-averaged forces:
    >>> tos = MagicTOHemi(tag='test', hemi='s', iplot=True)

.. _secTO_movieFile:

``TO_mov.TAG`` files
--------------------

.. note:: This file is **only** written when :ref:`l_TOmovie=.true. <varl_TOmovie>`

This file contains the time evolution of the different forces that enter the
phi-average of the azimuthal component of the Navier-Stokes equation. This is a
special kind of :ref:`movie file <secMovieFile>` that contains seven different
azimuthally-averaged fields in a :math:`(r,\theta)` plane : the axisymmetric
zonal flow component , the azimuthal component of the Reynolds stresses, the
azimuthal component of advection, the azimuthal component of viscosity, the
azimuthal component of Lorentz force, the azimuthal component of Coriolis force
and the azimuthal component of the time-derivative. The structure of the file
is similar to a :ref:`movie file <secMovieFile>`, i.e. an unformatted fortran binary
file with a header that describes the type of the movie file. The detailed calculations
can be found in the subroutine :f:subr:`outTO <out_to_mod/outto>`.

On a whole, the structure of the file looks like follows:

   .. code-block:: fortran

      !----------
      ! Line 1
      !----------

      version

      !----------
      ! Line 2
      !----------

      n_type, n_surface, const, n_fields

      !----------
      ! Line 3
      !----------

      runid

      !----------
      ! Line 4
      !----------

      n_r_movie_max, n_r_max, n_theta_max, n_phi_tot, minc, ra, ek, pr, 
      prmag, radratio, tScale

      !----------
      ! Line 5
      !----------

      r(1), r(2), ..., r(n_r_movie_max)

      !----------
      ! Line 6
      !----------

      theta(1), theta(2), ..., theta(n_theta_max)

      !----------
      ! Line 7
      !----------

      phi(1), phi(2), ..., phi(n_theta_max)

      ...

      !----------
      ! Line 7+N
      !----------

      n_frame, t_movie(N), omega_ic, omega_ma, dipLat, dipLon, dipStr, dipStrGeo

      !----------
      ! Line 7+(N+1)
      !----------

      vphi(t=t_movie(N),phi=1,theta=1), 
      vphi(t=t_movie(N),phi=1,theta=2), 
      ..., 
      vphi(t=t_movie(N),phi=n_phi_max,theta=n_theta_max)

      !----------
      ! Line 7+(N+2)
      !----------

      rey(t=t_movie(N),phi=1,theta=1), 
      rey(t=t_movie(N),phi=1,theta=2), 
      ..., 
      rey(t=t_movie(N),phi=n_phi_max,theta=n_theta_max)

      ...

      !----------
      ! Line 7+(N+7)
      !----------

      dtVphi(t=t_movie(N),phi=1,theta=1), 
      dtVphi(t=t_movie(N),phi=1,theta=2), 
      ..., 
      dtVphi(t=t_movie(N),phi=n_phi_max,theta=n_theta_max)


This file can be read using :py:class:`TOMovie <magic.TOMovie>` with the following options:

    >>> # To load 'TO_mov.test' and time-average it:
    >>> to = TOMOvie(file='TO_mov.test', avg=True, levels=65, cm='seismic')

