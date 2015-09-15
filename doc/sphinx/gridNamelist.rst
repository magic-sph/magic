.. _secGridNml:

Grid namelist
=============

This namelist defines the resolution of the computations. Keep in mind that **MagIC** is a 3D pseudo-spectral spherical shell code using Chebyshev polynomial expansions in the radial and spherical harmonic expansions in the angular directions.

Outer Core
----------

* **n_r_max** (default ``n_r_max=33``) is an integer which gives the number of grid points in the radial direction in the outer core (:math:`[r_i,r_o]`). It must be of the form ``4*n+1``, where ``n`` is an integer.

* **n_cheb_max** (default ``n_cheb_max=31``) is an integer which is the number of terms in the Chebyshev polynomial expansion to be used in the radial direction - the highest degree of Chebyshev polynomial used being ``n_cheb_max-1``. Note that ``n_cheb_max <= n_r_max``.

* **n_phi_tot** (default ``n_phi_tot=192``) is an integer which gives the number of longitudinal/azimuthal grid points. It has the following contraints:
 
  - ``n_phi_tot`` must be a multiple of ``minc`` (see below)

  - ``n_phi_tot/minc`` must be a multiple of 4

  - ``n_phi_tot`` must be a multiple of 16

Inner Core
----------

* **n_r_ic_max** (default ``n_r_ic_max=17``) is an integer which gives the number of grid points in the radial direction in the inner core (:math:`[0,r_i]`). It too, must be of the form ``4*n+1``, where ``n`` is an integer.

* **n_cheb_ic_max** (default ``n_cheb_ic_max=15``) is the number of terms in the Chebyshev polynomial expansion in the radial direction in the inner core. Only Chebyshev polynomials of even degrees are used in the expansion giving the highest degree used to be ``2*n_cheb_ic_max-2``. Note that here too, ``n_cheb_ic_max <= n_r_max``.

Symmetry and aliasing
---------------------

* **minc** (default ``minc=1``) is an integer which gives the longitudinal symmetry. e.g: ``minc=n`` would give an n-fold rotational symmetry in the azimuthal direction. One can use this to reduce computational costs when the symmetry of the solution is known. The orders of the spherical harmonic expansion (``m``) are multiples of ``minc``.

* **nalias** (default ``nalias=20``) is an integer which determines antialiasing used in the spherical harmonic representation. Note that ``20 <= nalias <= 30``.


The number of grid points in latitude ``n_theta_max = n_phi_tot/2``. The maximum degree (``l_max``) and maximum order (``m_max``) of the spherical harmonic expansion are determined by ``nalias``:

	``l_max = (nalias * n_theta_max)/30``
