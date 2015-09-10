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

.. code:: fortran

   &inner_core
     sigma_ratio = 1.0,
     nRotIc      = 1

Comas can be used to seperate namelist entries since they are not interpreted by the code.

Magic uses eight namelists:

1. ``&grid`` for resolution
2. ``&control`` for control parameters and numerical parameters.
3. ``&phys_param`` for the physical parameters.
4. ``&start_field`` to define the starting fields.
5. ``&output_control`` for defining the output.
6. ``&mantle`` for setting mantle parameters.
7. ``&inner_core`` for setting inner core parameters.
8. ``&B_external`` for setting up an external field contribution

.. toctree::
   :maxdepth: 1

   gridNamelist.rst
   physNamelist.rst
   controlNamelist.rst
   startNamelist.rst
   outputNamelist.rst
