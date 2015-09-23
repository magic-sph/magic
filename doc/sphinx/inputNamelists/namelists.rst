.. _secNamelists:

Input parameters
################

True runtime input parameters are read from STDIN as namelists, a Fortran
feature. A namelist is identified by its unique name `&name`. The
name-statement is followed by the parameters that are part of the namelist in
the format `parameter=value,`. The namelist is closed by a backslash. The
subroutine `defaultNamelists` (in the module `Namelist.f90`) defines a default
value for each parameter. Only the parameters whose value should differ from
its default have to be stated in the namelist.

An example for the short namelist defining inner core parameters is

.. code-block:: fortran

   &inner_core
     sigma_ratio = 1.0,
     nRotIc      = 1

Comas can be used to seperate namelist entries since they are not interpreted by the code.

Magic uses eight namelists:

1. :ref:`&grid <secGridNml>` for resolution
2. :ref:`&control <secControlNml>` for control parameters and numerical parameters.
3. :ref:`&phys_param <secPhysNml>` for the physical parameters.
4. :ref:`&B_external <secBextnml>` for setting up an external field contribution
5. :ref:`&start_field <secStartNml>` to define the starting fields.
6. :ref:`&output_control <secOutputNml>` for defining the output.
7. :ref:`&mantle <secMantle>` for setting mantle parameters.
8. :ref:`&inner_core <secInnerCore>` for setting inner core parameters.

.. toctree::
   :hidden:
   :maxdepth: 1

   gridNamelist.rst
   controlNamelist.rst
   physNamelist.rst
   BextNamelist.rst
   startNamelist.rst
   outNamelist.rst
   mantle_icNamelist.rst
