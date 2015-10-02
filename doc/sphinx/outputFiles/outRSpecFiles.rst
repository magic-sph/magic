
.. _secrBspecFiles:

Radial spectra ``rB[r|p]Spec.TAG``
----------------------------------

.. note:: This files are **only** written when :ref:`l_rMagSpec=.true. <varl_rMagSpec>`

Those files contain the time-evolution of the poloidal (**rBrSpec.TAG**) and
the toroidal (**rBpSpec.TAG**) magnetic energies for all radii including the
inner core and for spherical harmonic degrees from :math:`\ell=1` to
:math:`\ell=6`.  The calculations are done in the subroutines
:f:subr:`rBrSpec <radial_spectra/rbrspec()>` and :f:subr:`rBpSpec <radial_spectra/rbpspec()>`,
respectively. The outputs are stored as a fortran unformatted file which
follows the following structure for ``rBrSpec.TAG``:


   .. code-block:: fortran

       !-------------
       ! Line N
       !-------------

       time[N], 
       (real(e_p(l=1,n_r),kind=outp),n_r=1,n_r_tot-1), ! Poloidal energy for \ell=0
       (real(e_p(l=2,n_r),kind=outp),n_r=1,n_r_tot-1),
       ...
       (real(e_p(l=6,n_r),kind=outp),n_r=1,n_r_tot-1)  ! Poloidal energy for \ell=6

       !-------------
       ! Line N+1
       !-------------

       time[N], 
       (real(e_p_ax(l=1,n_r),kind=outp),n_r=1,n_r_tot-1), ! Poloidal energy for \ell=0 and m=0
       (real(e_p_ax(l=2,n_r),kind=outp),n_r=1,n_r_tot-1),
       ...
       (real(e_p_ax(l=6,n_r),kind=outp),n_r=1,n_r_tot-1)

       !-------------
       ! Line N+2
       !-------------

       time[N+1]
       ...

The ``rBpSpec.TAG`` files have exactly the same structure (just replacing the poloidal
energy by its toroidal counterpart).


.. warning:: Be careful that in this file, :f:var:`n_r_tot` is the **total** number of grid
             points (thus including the inner core).


Those files can be read using the python class :py:class:`MagicRSpec <magic.MagicRSpec>` with
the following options:

   >>> # Read the files BrSpec.testa, BrSpec.testb and BrSpec.testc and stack them
   >>> rsp = MagicRSpec(tag='test[a-c]', field='Br')
   >>> # Print time and the time evolution of e_pol(\ell=4) at the 10th radial grid point
   >>> print(rsp.time, rsp.e_pol[:, 10, 3])
