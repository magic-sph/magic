.. _secBextnml:

External Magnetic Field Namelist
================================

The namelist :code:`&B_external`  provides options for imposing an external magnetic field.

Externally imposed magnetic field
---------------------------------

* **n_imp** (default :f:var:`n_imp = 0 <n_imp>`) is an integer controling the type of external field applied.

  .. tabularcolumns:: |c|p{12cm}|
  
  +---------+-------------------------------------------------------+
  | n_imp=0 | No external magnetic field                            |
  +---------+-------------------------------------------------------+
  | n_imp=1 | Follows idea of Uli Christensen of external field     |
  |         | compensating internal field such that radial component| 
  |         | of magnetic field vanishes at :math:`r/r_{cmb}=rrMP`  |
  |         | where ``rrMP`` is the 'magnetopause radius' input by  |
  |         | the user (see below)                                  |
  +---------+-------------------------------------------------------+
  | n_imp=2 | Uniform axisymmetric magnetic field of geometry given |
  |         | by ``l_imp`` (see below)                              |
  +---------+-------------------------------------------------------+
  | n_imp=3 | Uniform axisymmetric magnetic field which changes     |
  |         | direction according to the direction of the axial     |
  |         | dipole of the internal magnetic field                 |
  +---------+-------------------------------------------------------+
  | n_imp=4 | Same as ``n_imp=3`` but the amplitude of the external |
  |         | field is scaled to the amplitude of the axial dipole  |
  |         | of the internal field                                 |
  +---------+-------------------------------------------------------+
  | n_imp=7 | External field depends on internal axial dipole       |
  |         | through Special Heyner feedback functions             |
  +---------+-------------------------------------------------------+

* **rrMP** (default :f:var:`rrMP = 0.0 <rrmp>`) is a real which gives the value of 'magnetopause radius'. In other words, it gives the radius (as a fraction of ``r_cmb``) at which the radial component of the magnetic field vanishes due to cancelling out of external and internal magnetic field components. Used only when ``n_imp = 1``.

* **amp_imp** (default :f:var:`amp_imp = 0.0 <amp_imp>`) is a real which gives the amplitude of the external magnetic field.

* **expo_imp** (default :f:var:`expo_imp = 0.0 <expo_imp>`) is a real which gives the exponent of dependence of external magnetic field on the axial dipole of the internal magnetic field. Used for ``n_imp=7``.

* **bmax_imp** (default :f:var:`bmax_imp = 0.0 <bmax_imp>`) is a real which gives the location of the maximum of the ratio of the poloidal potentials :math:`g_{ext}/g_{int}`.

* **l_imp** (default :f:var:`l_imp = 1 <l_imp>`) is an integer which gives the geometry (degree of spherical harmonic) of the external magnetic field. The external field is always axisymmetric, hence ``m = 0`` always. This option is used when ``n_imp = 2,3`` or ``4``.

Current carrying loop
---------------------

To simulate experiments, an external current carrying loop, concentric to the sphere and in the equatorial plane, has been implemented in the code. It's radius is fixed at a distance :math:`a = r_{cmb}/0.8` to match conditions of the Maryland 3 metre experiment.

* **l_curr** (default :f:var:`l_curr = .false. <l_curr>`) is a logical that controls switching on or off of the current carrying loop.

* **amp_curr** (default :f:var:`amp_curr = 0.0 <amp_curr>`) is a real that gives the amplitude of magnetic field produced by the current carring loop.

.. warning::

 Note that an external magnetic field is incompatible with a region of low conductivity inside the spherical shell (i.e, if ``r_LCR < r_cmb``). Thus, while imposing an external magnetic field, make sure ``r_LCR > r_cmb`` (which is the default case). For details on ``r_LCR``, have a look at the section on :ref:`electrical conductivity <varnVarCond>` in the namelist for :ref:`physical parameters <secPhysNml>`.
