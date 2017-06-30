.. _secOutputFiles:

Output files
############

While some information of a run is written into ``STDOUT`` to monitor its
progress, most outputs are printed into dedicated files identified by the chosen
:ref:`TAG <varTag>` extension. These files can be parsed and analysed using the
:ref:`python classes <secPythonPostProc>`. The following pages describe the content
and the structure of the different type of output files:

   1. Most of the information found in ``STDOUT`` is
      also written to the **log-file** called :ref:`log.TAG <secLogFile>`. In addition,
      this file contains all input parameters, truncation, information on other
      output files, and some results like the time averaged energies (when
      :ref:`l_average=.true. <varl_average>`). 

   2. There are several ascii files that contain the **time-evolution of
      integrated quantities** (energies, heat fluxes, rotation rate, Reynolds
      numbers, etc.) that are systematically produced: 

         * Kinetic energies: :ref:`e_kin.TAG <secEkinFile>`, 
	 * Magnetic energies: :ref:`e_mag_oc.TAG <secEmagocFile>` and :ref:`e_mag_ic.TAG <secEmagicFile>`,
         * Rotation rates: :ref:`rot.TAG <secRotFile>`, 
	 * Informations about the dipolar component of the magnetic field: :ref:`dipole.TAG <secDipoleFile>`,
	 * Diagnostic parameters (Reynolds, Elsasser, etc.): :ref:`par.TAG <secParFile>`,
      
   3. There are **additional conditional time series** that contain the time-evolution of 
      other physical quantities that depend on the chosen 
      :ref:`input parameters <secOutputNml>`:

         * Angular momentum balance: :ref:`AM.TAG <secAMFile>`,
	 * Heat transport: :ref:`heat.TAG <secHeatFile>`,
	 * Helicity: :ref:`helicity.TAG <secHelicityFile>`,
	 * Power budget: :ref:`power.TAG <secpowerFile>` and :ref:`dtE.TAG <secdtEFile>`,
	 * Square velocities: :ref:`u_square.TAG <secu_squareFile>`,
	 * Drift rates: :ref:`drift[V|B][D|Q].TAG <secdriftFile>` and :ref:`iner[P|T].TAG <secinerFile>`,
	 * Torques: :ref:`SR[IC|MA].TAG <secSRFile>`,
	 * Geostrophy: :ref:`geos.TAG <secGeosFile>`,
	 * RMS calculations of the force balances: :ref:`dtVrms.TAG <secdtVrmsFile>` and :ref:`dtBrms.TAG <secdtBrmsFile>`,
	 * Kinetic energies perpendicular and parallel to the rotation axis: :ref:`perpPar.TAG <secperpParFile>`.

   4. **Time-averaged radial profiles**:

         * Kinetic energies: :ref:`eKinR.TAG <secEkinRFile>`,
         * Magnetic energies: :ref:`eMagR.TAG <secEmagRFile>`,
         * Diagnostic quantities: :ref:`parR.TAG <secParRfile>`,
         * Power budget: :ref:`powerR.TAG <secPowerRfile>`,
	 * Average temperature, entropy and pressure: :ref:`heatR.TAG <secHeatRfile>`,
	 * Heat fluxes: :ref:`fluxesR.TAG <secFluxesRfile>`,
	 * Temperature and horizontal velocities: :ref:`bLayersR.TAG <secBLayersRfile>`,
	 * Kinetic energies perpendicular and parallel to the rotation axis: :ref:`perpParR.TAG <secPerpParRfile>`.

   5. **Radial profiles of the transport properties** of the reference state (those files will
      only be produced when the appropriate input option is chosen):

         * Temperature, density and gravity: :ref:`anel.TAG <secAnelFile>`,
	 * Electrical conductivity: :ref:`varCond.TAG <secVarCondFile>`,
	 * Thermal conductivity: :ref:`varDiff.TAG <secVarDiffFile>`,
	 * Kinematic viscosity: :ref:`varVisc.TAG <secVarViscFile>`,
	 * Mapping of the Chebyshev grid: :ref:`rNM.TAG <secMappingFile>`.

   6. Kinetic energy, magnetic energy and temperature/entropy **spectra**:

         * Kinetic energy: :ref:`kin_spec_#.TAG <secKinSpecFile>`,
         * Magnetic energy: :ref:`kin_spec_#.TAG <secMagSpecFile>`,
         * Velocity square: :ref:`u2_spec_#.TAG <secu2SpecFile>`,
         * Temperature/entropy: :ref:`T_spec_#.TAG <secMagSpecFile>`,
         * Time-averaged kinetic energy: :ref:`kin_spec_ave.TAG <secKinSpecAveFile>`,
         * Time-averaged magnetic energy: :ref:`mag_spec_ave.TAG <secMagSpecAveFile>`,
         * Time-averaged temperature/entropy: :ref:`T_spec_ave.TAG <secTempSpecAveFile>`,
	 * 2-D ([r,\ell] and [r,m]) spectra: :ref:`2D_[mag|kin|u2]_spec_#.TAG <sec2DSpectra>`.

   7. Output snapshot that contains the 3-D components of the velocity field, the magnetic
      field and the temperature/entropy. Those files are named **graphic files**
      :ref:`G_#.TAG <secGraphFile>` (or :ref:`G_ave.TAG <secGraphFile>` for its time-averaged
      counterpart). 

   8. Time evolution of some chosen fields. Those files are named **movie files**:
      :ref:`*_mov.TAG <secMovieFile>`.

   9. Checkpoints outputs that will allow the code to restart. Those files are named
      **restart files**: :ref:`checkpoint_end.TAG <secRestartFile>`.

   10. **Time-evolution of the poloidal and toroidal coefficients** at diffent depths:

         * Time evolution of the poloidal magnetic field at the CMB: :ref:`B_coeff_cmb.TAG <secCmbFile>`,
	 * Time evolution of the potentials at several depths: :ref:`[V|T|B]_coeff_r#.TAG <secCoeffrFiles>`

   11. **Additional specific outputs**:

         * Torsional oscillations (see :ref:`here <secTOoutputFiles>`),
	 * Potential files: :ref:`V_lmr_#.TAG <secVpotFile>`, :ref:`B_lmr_#.TAG <secBpotFile>` and :ref:`T_lmr_#.TAG <secTpotFile>`,
	 * Potential vorticity files: ``PVZ.TAG`` and ``Vcy.TAG``,
	 * Magnetic spectra for various radii: :ref:`rB[r|p]Spec.TAG <secrBspecFiles>`.



.. toctree::
   :hidden:

   logFile.rst
   outTimeSeriesStd.rst
   outTimeSeriesSupp.rst
   outRadialFiles.rst
   outTransportProp.rst
   outSpecFiles.rst
   outGraph.rst
   outMovie.rst
   outRst.rst
   outCoeffFiles.rst
   outTOFiles.rst
   outRSpecFiles.rst
   outPotFiles.rst
