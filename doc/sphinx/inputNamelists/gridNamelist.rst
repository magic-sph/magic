.. _secGridNml:

Grid namelist
=============

This namelist defines the resolution of the computations. Keep in mind that **MagIC** is a 3D pseudo-spectral spherical shell code using Chebyshev polynomial expansions in the radial and spherical harmonic expansions in the angular directions.

Outer Core
----------

.. _varn_r_max:

* **n_r_max** (default :f:var:`n_r_max=33 <n_r_max>`) is an integer which gives the number of grid points in the radial direction in the outer core (:math:`[r_i,r_o]`). It must be of the form ``4*n+1``, where ``n`` is an integer.

  .. note:: The possible values for :f:var:`n_r_max` are thus: 17, 21, 25, 33, 37, 41, 49, 61, 65, 73? 81, 97, 101, 121, 129, 145, 161, 257, 401, 513, ...

* **n_cheb_max** (default :f:var:`n_cheb_max=31 <n_cheb_max>`) is an integer which is the number of terms in the Chebyshev polynomial expansion to be used in the radial direction - the highest degree of Chebyshev polynomial used being ``n_cheb_max-1``. Note that ``n_cheb_max <= n_r_max``.

  .. note:: Adopting ``n_cheb_max=n_r_max-2`` is usually a good choice

* **n_phi_tot** (default :f:var:`n_phi_tot=192 <n_phi_tot>`) is an integer which gives the number of longitudinal/azimuthal grid points. It has the following constraints:
 
  - :f:var:`n_phi_tot` must be a multiple of :f:var:`minc` (see below)

  - ``n_phi_tot/minc`` must be a multiple of 4

  - :f:var:`n_phi_tot` must be a multiple of 16

  .. note:: The possible values for :f:var:`n_phi_max` are thus: 16, 32, 48, 64, 96, 128, 192, 256, 288? 320, 384, 400, 512, 576, 640, 768, 864, 1024, 1280, 1536, 1792, 2048, ...

Inner Core
----------

* **n_r_ic_max** (default :f:var:`n_r_ic_max=17 <n_r_ic_max>`) is an integer which gives the number of grid points in the radial direction in the inner core (:math:`[0,r_i]`). It too, must be of the form ``4*n+1``, where ``n`` is an integer.

* **n_cheb_ic_max** (default :f:var:`n_cheb_ic_max=15 <n_cheb_ic_max>`) is the number of terms in the Chebyshev polynomial expansion in the radial direction in the inner core. Only Chebyshev polynomials of even degrees are used in the expansion giving the highest degree used to be ``2*n_cheb_ic_max-2``. Note that here too, ``n_cheb_ic_max <= n_r_max``.

Symmetry and aliasing
---------------------

.. _varMinc:

* **minc** (default :f:var:`minc=1 <minc>`) is an integer which gives the longitudinal symmetry. e.g: ``minc=n`` would give an n-fold rotational symmetry in the azimuthal direction. One can use this to reduce computational costs when the symmetry of the solution is known. The orders of the spherical harmonic expansion (``m``) are multiples of :f:var:`minc`.

* **nalias** (default :f:var:`nalias=20 <nalias>`) is an integer which determines antialiasing used in the spherical harmonic representation. Note that ``20 <= nalias <= 30``.


The number of grid points in latitude :f:var:`n_theta_max = n_phi_tot/2 <n_theta_max>`. The
maximum degree (:f:var:`l_max`) and maximum order (:f:var:`m_max`) of the spherical
harmonic expansion are determined by :f:var:`nalias`:

  .. code-block:: fortran

	l_max = (nalias * n_theta_max)/30
