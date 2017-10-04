**********
User Guide
**********

Numerical implementation
========================

Icepack is written in FORTRAN90 and runs on platforms using UNIX, LINUX,
and other operating systems. The code is not parallelized. (CHANGE IF OPENMP IS IMPLEMENTED)

Icepack consists of the sea ice column physics code, contained in the 
**columnphysics/** directory, and a **configuration/** directory that includes
a driver for testing the column physics and a set of scripts for configuring the tests.
Icepack is designed such that the column physics code may be used by a host sea ice
model without direct reference to the driver or scripts, although these may be consulted for 
guidance when coupling the column physics code to the host sea ice model 
(`CICE <https://github.com/CICE-Consortium/CICE>`_ may also be useful for this.)  Information
about the interface between the column physics and the driver or host sea ice model is
located in the :ref:`init` section.

.. _dirstructure:

Directory structure
-------------------

The present code distribution includes make files, several scripts and
some input files. One year of atmospheric forcing data is also available
from the code distribution web site (see the **README** file for
details).

**LICENSE.pdf**

**DistributionPolicy.pdf**
    license and policy for using and sharing the code

**README.md**
    basic information and pointers

**columnphysics/**
    the essential physics code

    **constants**  (CHECK - this will change)
      physical and numerical constants required for column package

      **cesm/icepack_constants.F90**

      **cice/icepack_constants.F90**

      **hadgem3/icepack_constants.F90**

    **icepack\_aerosol.F90**
        handles most work associated with the aerosol tracers

    **icepack\_age.F90**
        handles most work associated with the age tracer

    **icepack\_algae.F90**
        biogeochemistry

    **icepack\_atmo.F90**
        stability-based parameterization for calculation of turbulent ice–atmosphere fluxes

    **icepack\_brine.F90**
        evolves the brine height tracer

    **icepack\_firstyear.F90**
        handles most work associated with the first-year ice area tracer

    **icepack\_flux.F90**
        fluxes needed/produced by the model

    **icepack\_intfc.F90**
        interface routines for linking Icepack with a host sea ice model

    **icepack\_intfc\_shared.F90**
        interface routines for linking Icepack with a host sea ice model

    **icepack\_intfc\_tracers.F90**
        interface routines for linking Icepack with a host sea ice model

    **icepack\_itd.F90**
        utilities for managing ice thickness distribution

    **icepack\_kinds\_mod.F90**
        basic definitions of reals, integers, etc.

    **icepack\_mechred.F90**
        mechanical redistribution (ridging)

    **icepack\_meltpond\_cesm.F90**
        CESM melt pond parameterization

    **icepack\_meltpond\_lvl.F90**
        level-ice melt pond parameterization

    **icepack\_meltpond\_topo.F90**
        topo melt pond parameterization

    **icepack\_ocean.F90**  (CHECK THIS, not in directory now)
        mixed layer ocean model

    **icepack\_mushy\_physics.F90**
        physics routines for mushy thermodynamics

    **icepack\_orbital.F90**
        orbital parameters for Delta-Eddington shortwave parameterization

    **icepack\_shortwave.F90**
        shortwave and albedo parameterizations

    **icepack\_therm\_0layer.F90**
        zero-layer thermodynamics of :cite:`Semtner76`

    **icepack\_therm\_bl99.F90**
        multilayer thermodynamics of :cite:`BL99`

    **icepack\_therm\_itd.F90**
        thermodynamic changes mostly related to ice thickness distribution

    **icepack\_therm\_mushy.F90**
        mushy-theory thermodynamics of :cite:`THB13`

    **icepack\_therm\_shared.F90**
        code shared by all thermodynamics parameterizations

    **icepack\_therm\_vertical.F90**
        vertical growth rates and fluxes

    **icepack\_warnings.F90**
        utilities for writing warning and error messages

    **icepack\_zbgc.F90**
        driver for ice biogeochemistry and brine tracer motion

    **icepack\_zbgc\_shared.F90**
        parameters and shared code for biogeochemistry and brine height

    **icepack\_zsalinity.F90**
        vertical salinity parameterization of :cite:`JHE11`

**configuration/**
    drivers and scripts for testing Icepack in stand-alone mode
    
    **driver/**
        **icepack\_drv\_MAIN.F90**
            main program

        **icepack\_drv\_InitMod.F90**
            routines for initializing a run

        **icepack\_drv\_RunMod.F90**
            main driver routines for time stepping

        **icepack\_drv\_arrays\_column.F90**
            essential arrays to describe the state of the ice

        **icepack\_drv\_calendar.F90**
            keeps track of what time it is

        **icepack\_drv\_constants.F90**
            physical and numerical constants and parameters

        **icepack\_drv\_diagnostics.F90**
            miscellaneous diagnostic and debugging routines

        **icepack\_drv\_diagnostics\_bgc.F90**
            diagnostic routines for biogeochemistry

        **icepack\_drv\_domain\_size.F90**
            domain sizes

        **icepack\_drv\_flux.F90**
            fluxes needed/produced by the model

        **icepack\_drv\_forcing.F90**
            routines to read and interpolate forcing data for stand-alone model runs

        **icepack\_drv\_init.F90**
            general initialization routines

        **icepack\_drv\_init\_column.F90**
            initialization routines specific to the column physics

        **icepack\_drv\_restart.F90**
            driver for reading/writing restart files

        **icepack\_drv\_restart\_column.F90**  (CHECK: RENAME bgc)
            restart routines specific to the column physics

        **icepack\_drv\_restart\_shared.F90**
            code shared by all restart options

        **icepack\_drv\_state.F90**
            essential arrays to describe the state of the ice

        **icepack\_drv\_step\_mod.F90**
            routines for time stepping the major code components

    **scripts/**
        **Makefile**
            primary makefile

        **icepack.batch.csh**
            creates batch scripts for particular machines

        **icepack.build**
            compiles the code

        **icepack.launch.csh**
            creates script logic that runs the executable

        **icepack.run.setup.csh**
            sets up the run directory

        **icepack.settings**
            defines environment, model configuration and run settings

        **icepack.test.setup.csh**
            creates configurations for testing the model

        **icepack\_decomp.csh**
            defines the grid size

        **icepack\_in**
            namelist input data

        **machines/**
            macro definitions for the given computers

        **makdep.c**
            determines module dependencies

        **options/**
            other namelist configurations available from the icepack.create.case command line

        **parse\_namelist.sh**
            replaces namelist with command-line configuration

        **parse\_settings.sh**
            replaces settings with command-line configuration

        **tests/**
            scripts for configuring and running basic tests

**doc/**
    documentation

**icepack.create.case**
    main script for setting up a test case

A case (compile) directory is created upon initial execution of the script 
**icepack.create.case** at the user-specified location provided after the -c flag. 
Executing the command ``./icepack.create.case -h`` provides helpful information for 
this tool. Please refer to the `user guide <https://CICE-Consortium.github.io/Icepack/index.html>`_ 
for further information.

.. _grids:

Grid and boundary conditions 
-----------------------------------

The driver configures a collection of grid cells on which the column physics code 
will be run. This "horizontal" grid is a vector of length ``nx``, with a minimum length 
of 4.   
The grid vector is initialized with different sea ice conditions, such as open 
water, a uniform slab of ice, a multi-year ice thickness distribution with snow, 
and land. For simplicity, the same forcing values are applied to all grid cells. 

Icepack includes two vertical grids.  The basic vertical grid contains 
``nilyr`` equally spaced grid cells.  
History variables available for column output are ice and snow
temperature, ``Tinz`` and ``Tsnz``. These variables also include thickness
category as a fourth dimension.

In addition, there is a bio-grid that 
can be more finely resolved and includes additional nodes for boundary conditions.
It is used for solving the brine height variable :math:`h_b` and for
discretizing the vertical transport equations of biogeochemical tracers.
The bio-grid is a non-dimensional vertical grid which takes the value
zero at :math:`h_b` and one at the ice–ocean interface. The number of
grid levels is specified during compilation by setting
the variable ``NBGCLYR`` equal to an integer (:math:`n_b`) .

Ice tracers and microstructural properties defined on the bio-grid are
referenced in two ways: as ``bgrid`` :math:`=n_b+2` points and as
igrid\ :math:`=n_b+1` points. For both bgrid and igrid, the first and
last points reference :math:`h_b` and the ice–ocean interface,
respectively, and so take the values :math:`0` and :math:`1`,
respectively. For bgrid, the interior points :math:`[2, n_b+1]` are
spaced at :math:`1/n_b` intervals beginning with `bgrid(2)` = 
:math:`1/(2n_b)`. The ``igrid`` interior points :math:`[2, n_b]` are also
equidistant with the same spacing, but physically coincide with points
midway between those of ``bgrid``.


.. _testconfigs:

Test configurations
-------------------

*(CHECK) UPDATE with similar, correct information*

The column is located
near Barrow (71.35N, 156.5W). Options for choosing the column
configuration are given in **comp\_ice** (choose `RES col`) and in the
namelist file, **input\_templates/col/ice\_in**. Here, ``istep0`` and the
initial conditions are set such that the run begins September 1 with no
ice. 


.. _init:

Initialization and coupling
---------------------------

CHECK:  link to information about the column physics interface in section 2

Icepack’s parameters and variables are initialized in several
steps. Many constants and physical parameters are set in
**icepack\_constants.F90**. Namelist variables (:ref:`tabnamelist`),
whose values can be altered at run time, are handled in *input\_data*
and other initialization routines. These variables are given default
values in the code, which may then be changed when the input file
**icepack\_in** is read. Other physical constants, numerical parameters, and
variables are first set in initialization routines for each ice model
component or module. Then, if the ice model is being restarted from a
previous run, core variables are read and reinitialized in
*restartfile*, while tracer variables needed for specific configurations
are read in separate restart routines associated with each tracer or
specialized parameterization. Finally, albedo and other quantities
dependent on the initial ice state are set. Some of these parameters
will be described in more detail in the :ref:`tabnamelist`.

Two namelist variables control model initialization, ``ice_ic``
and ``restart``.  Setting ``ice_ic`` = 'default' causes the model to run using
constant forcing and initial values set in the code.  To start
from a file **filename**, set 
``restart`` = .true. and ``ice_ic`` = **filename**.  When restarting using the Icepack
driver, for simplicity the tracers are assumed to be set the same way (on/off) as in the
run that created the restart file; i.e. that the restart file contains exactly the 
information needed for the new run.  CICE is more flexible in this regard.

For stand-alone runs,
routines in **icepack\_drv\_forcing.F90** read and interpolate data from files,
and are intended merely for testing, although they can also provide guidance for 
the user to write his or her own routines. 


.. _parameters:

Choosing an appropriate time step
---------------------------------

Transport in thickness space imposes a restraint on the time
step, given by the ice growth/melt rate and the smallest range of
thickness among the categories,
:math:`\Delta t<\min(\Delta H)/2\max(f)`, where :math:`\Delta H` is the
distance between category boundaries and :math:`f` is the thermodynamic
growth rate. For the 5-category ice thickness distribution used as the
default in this distribution, this is not a stringent limitation:
:math:`\Delta t < 19.4` hr, assuming :math:`\max(f) = 40` cm/day.


Model output
------------

.. _history:

History output from Icepack is not currently supported in the Icepack driver, except
in restart files.
The sea ice model `CICE <https://github.com/CICE-Consortium/CICE>`_ provides extensive 
options for model output, including many derived output variables.

Diagnostic files
~~~~~~~~~~~~~~~~

Icepack writes diagnostic information for each grid cell as a separate file, 
**ice\_diag.\***, identified by the initial ice state of the grid cell (no ice, slab, land, etc).


Restart files
~~~~~~~~~~~~~

CHECK and CHANGE as needed re netCDF

CICE provides restart data in binary unformatted or netCDF formats, via
the ``IO_TYPE`` flag in **comp\_ice** and namelist variable
``restart_format``. 

The restart files created by the Icepack driver contain all of the variables needed
for a full, exact restart. The filename begins with the character string
‘iced.’, and the restart dump frequency is given by the namelist
variable ``dumpfreq``. The namelist variable ``ice_ic`` contains the
pointer to the filename from which the restart data is to be read.


Execution procedures
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

Scripts
-------------

Most of the scripts that configure, build and run Icepack tests are contained in 
the directory **configuration/scripts/**, except for **icepack.create.case**, which is
in the main directory.  

Users likely will need to create or edit some scripts for their computer's environment. 
Specific instructions for porting are provided below.

icepack.create.case generates a case. Use ``create.case -h`` for help with the tool.
  -c is the case name and location (required)

  -m is the machine name (required). Currently, there are working ports for NCAR yellowstone and cheyenne, AFRL thunder, NavyDSRC gordon and conrad, and LANL’s wolf machines.

  -s are comma separated optional env or namelist settings (default is 'null')

  -t is the test name and location (cannot be used with -c).

  -bd is used to specify the location of the baseline datasets (only used with -t)

  -bg is used to specify the icepack version name for generating baseline datasets (only used with -t)

  -bc is used to specify the icepack version name for comparison. I.e., the version name for the baseline dataset (only used with -t)

  -testid is used to specify a test ID (used only with -t or -ts)

  -ts is used to generate all test cases for a given test suite.


Several files are placed in the case directory

- **env.[machine]** defines the environment

- **icepack.settings** defines many variables associated with building and running the model

- **makdep.c** is a tool that will automatically generate the make dependencies

- **Macros.[machine]** defines the Makefile macros

- **Makefile** is the makefile used to build the model

- **icepack.build** is a script that builds and compiles the model

- **icepack\_in** is the namelist file

- **icepack.run** is a batch run script

- **icepack.submit** is a simple script that submits the icepack.run script

Once the case is created, all scripts and namelist are fully resolved.  Users can edit any
of the files in the case directory manually to change the model configuration.  The file
dependency is indicated in the above list.  For instance, if any of the files before
**icepack.build** in the list are edited, **icepack.build** should be rerun.

The **casescripts/** directory holds scripts used to create the case and can 
largely be ignored.  
In general, when **icepack.build** is executed, the model will build from scratch 
due to extensive preprocessing dependencies.  To change this behavior, edit the 
env variable ``ICE_CLEANBUILD`` in **icepack.settings**.  

The **icepack.submit** script simply submits the **icepack.run script**.  
You can also submit the **icepack.run** script on the command line.

To port, an **env.[machine]** and **Macros.[machine]** file have to be added to 
**configuration/scripts/machines/** and the 
**icepack.run.setup.csh** file needs to be modified.
 
- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files

- cd .. to **configuration/scripts/**

- Edit the **icepack.run.setup.csh** script to add a section for your machine 
  with batch settings and job launch settings

- Download and untar a forcing dataset to the location defined by 
  ``ICE_MACHINE_INPUTDATA`` in the env file

- Create a file in your home directory called **.cice\_proj** and add your preferred account name to the first line.


Directories
-------------

CHECK

The **icepack.create.case** script creates a case directory in the location specified 
by the ``-c`` or ``-t`` flags.  The **icepack.build** script 
creates the run directory defined by the env variable ``ICE_RUNDIR`` in 
**icepack.settings**, and it compiles the code there.  The run directory is further 
populated by the **icepack.run** script, which also runs the executable.  Specifying 
the test suite creates a directory containing subdirectories for each test.

Build and run logs will be copied from the run directory into the case **logs/** 
directory when complete.


Local modifications
--------------------------

Scripts and files can be changed in the case directory and then built from there, without 
changing them in your main directory.

You also can directly modify the namelist files (**icepack\_in**) in the run directory and
run the code by submitting the executable **icepack** directly.  Beware that any changes 
made in the run directory will be overwritten if scripts are later run from the case
directory.

Forcing data
------------

CHECK once we've settled on a forcing suite:

The code is currently configured to run in standalone mode on a 4-cell grid using 
atmospheric data, available as detailed on the `wiki <https://github.com/CICE-Consortium/Icepack/wiki/Testing-Icepack>`_.
These data files are designed only for testing the code, not for use in production 
runs or as observational data.  Please do not publish results based on these data
sets.  Module **configuration/driver/icepack\_drv\_forcing.F90**
can be modified to change the forcing data. 


Adding things
====================

We require that any changes made to the code be implemented in such a way that they can
be "turned off" through namelist flags.  In most cases, code run with such changes should 
be bit-for-bit identical with the unmodified code.  Occasionally, non-bit-for-bit changes
are necessary, e.g. associated with an unavoidable change in the order of operations. In
these cases, changes should be made in stages to isolate the non-bit-for-bit changes, 
so that those that should be bit-for-bit can be tested separately.

Tracers
--------------

.. _addtrcr:

Tracers added to Icepack will also require extensive modifications to the host
sea ice model, including initialization on the horizontal grid, namelist flags 
and restart capabilities.  Modifications to the Icepack driver should reflect
the modifications needed in the host model but are not expected to match completely.
We recommend that the logical namelist variable
``tr_[tracer]`` be used for all calls involving the new tracer outside of
**ice\_[tracer].F90**, in case other users do not want to use that
tracer.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, level-ice area and all of the
biogeochemistry tracers in this release are area-weighted tracers. In
the absence of sources and sinks, the total mass of a volume-weighted
tracer such as aerosol (kg) is conserved under transport in horizontal
and thickness space (the mass in a given grid cell will change), whereas
the aerosol concentration (kg/m) is unchanged following the motion, and
in particular, the concentration is unchanged when there is surface or
basal melting. The proper units for a volume-weighted mass tracer in the
tracer array are kg/m.

In several places in the code, tracer computations must be performed on
the conserved "tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies, and code must be included to account for new
dependencies (e.g., :math:`a_{lvl}` and :math:`a_{pnd}` in
**ice\_itd.F90** and **ice\_mechred.F90**).

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **icepack\_drv\_domain\_size.F90**: increase ``max_ntrcr`` (can also add option
   to **icepack.settings** and **icepack.build**)

#. **icepack\_drv\_state.F90**: declare ``nt_[tracer]`` and ``tr_[tracer]``

#. **icepack\_[tracer].F90**: create initialization, physics routines

#. **ice\_drv\_init.F90**: (some of this may be done in **ice\_[tracer].F90**
   instead)

   -  add new module and ``tr_[tracer]`` to list of used modules and
      variables

   -  add logical namelist variable ``tr_[tracer]``

   -  initialize namelist variable

   -  print namelist variable to diagnostic output file

   -  increment number of tracers in use based on namelist input (``ntrcr``)

   -  define tracer types (``trcr_depend`` = 0 for ice area tracers, 1 for
      ice volume, 2 for snow volume, 2+``nt_``[tracer] for dependence on
      other tracers)

#. **icepack\_itd.F90**, **icepack\_mechred.F90**: Account for new dependencies
   if needed.

#. **icepack\_drv\_InitMod.F90**: initialize tracer (includes reading restart
   file)

#. **icepack\_drv\_RunMod.F90**, **icepack\_drv\_step\_mod.F90**:

   -  call routine to write tracer restart data

   -  call physics routines in **icepack\_[tracer].F90** (often called from
      **icepack\_drv\_step\_mod.F90**)

#. **icepack\_drv\_restart.F90**: define restart variables

#. **icepack\_in**: add namelist variables to *tracer\_nml* and
   *icefields\_nml*

#. If strict conservation is necessary, add diagnostics as noted for
   topo ponds in Section :ref:`ponds`.


Troubleshooting 
================

Check the FAQ: https://github.com/CICE-Consortium/Icepack/wiki

.. _setup:

Initial setup
-------------

If there are problems, you can manually edit 
the env, Macros, and **icepack.run** files in the case directory until things are 
working properly.  Then you can copy the env and Macros files back to 
**configuration/scripts/machines**.  

- Changes made directly in the run directory, e.g. to the namelist file, will be overwritten
  if scripts in the case directory are run again later.

- If changes are needed in the **icepack.run.setup.csh** script, it must be manually modified.

.. _restarttrouble:

Restarts
--------

- Manual restart tests require the path to the restart file be included in ``ice_in`` in the 
  namelist file.

- Ensure that ``kcatbound`` is the same as that used to create the restart file.  
  Other configuration parameters, such as ``NICELYR``, must also be consistent between runs.

.. _testtrouble:

Testing
--------

- Tests using a debug flag that traps underflows will fail unless a "flush-to-zero" flag 
  is set in the Macros file.  This is due to very small exponential values in the delta-Eddington
  radiation scheme.

Debugging hints
---------------

CHECK write utility in column physics interface, for checking parameter values

A printing utility is available in the driver that can be helpful when debugging the
code. Not all of these will work everywhere in the code, due to possible
conflicts in module dependencies.

*debug\_icepack* (**configuration/driver/ice\_diagnostics.F90**)
    A wrapper for *print\_state* that is easily called from numerous
    points during initialization and the timestepping loop

*print\_state* (**configuration/driver/ice\_diagnostics.F90**)
    Print the ice state and forcing fields for a given grid cell.

Known bugs
----------

-   With the old CCSM radiative scheme (``shortwave`` = ‘default’ or
    ‘ccsm3’), a sizable fraction (more than 10%) of the total shortwave
    radiation is absorbed at the surface but should be penetrating into
    the ice interior instead. This is due to use of the aggregated,
    effective albedo rather than the bare ice albedo 
    when ``snowpatch`` < 1.

Interpretation of albedos
-------------------------

The snow-and-ice albedo, ``albsni``, and diagnostic albedos ``albice``, ``albsno``,
and ``albpnd`` are merged over categories but not scaled (divided) by the
total ice area. (This is a change from CICE v4.1 for ``albsni``.) The latter
three history variables represent completely bare or completely snow- or
melt-pond-covered ice; that is, they do not take into account the snow
or melt pond fraction (``albsni`` does, as does the code itself during
thermodyamic computations). This is to facilitate comparison with
typical values in measurements or other albedo parameterizations. The
melt pond albedo ``albpnd`` is only computed for the Delta-Eddington
shortwave case.

With the Delta-Eddington parameterization, the albedo depends on the
cosine of the zenith angle (:math:`\cos\varphi`, ``coszen``) and is zero if
the sun is below the horizon (:math:`\cos\varphi < 0`). Therefore
time-averaged albedo fields would be low if a diurnal solar cycle is
used, because zero values would be included in the average for half of
each 24-hour period. To rectify this, a separate counter is used for the
averaging that is incremented only when :math:`\cos\varphi > 0`. The
albedos will still be zero in the dark, polar winter hemisphere.

Proliferating subprocess parameterizations
------------------------------------------

With the addition of several alternative parameterizations for sea ice
processes, a number of subprocesses now appear in multiple parts of the
code with differing descriptions. For instance, sea ice porosity and
permeability, along with associated flushing and flooding, are
calculated separately for mushy thermodynamics, topo and level-ice melt
ponds, and for the brine height tracer, each employing its own
equations. Likewise, the BL99 and mushy thermodynamics compute freeboard
and snow–ice formation differently, and the topo and level-ice melt pond
schemes both allow fresh ice to grow atop melt ponds, using slightly
different formulations for Stefan freezing. These various process
parameterizations will be compared and their subprocess descriptions
possibly unified in the future.


.. _testing:

Testing Icepack
================

.. _basic:

Individual tests and test suites
--------------------------------

The Icepack scripts support both setup of individual tests as well as test suites.  Individual
tests are run from the command line,

  ``> ./icepack.create.case -t smoke -m wolf -s diag1,debug -testid myid``

where -m designates a specific machine.  Test suites are multiple tests that are specified in 
an input file and are started on the command line,

  ``> ./icepack.create.case -ts base_suite -m wolf -testid myid``

Invoking **icepack.create.case** with -t or -ts requires a testid to uniquely name test directories.  The format
of the case directory name for a test will always be 
${machine}_${test}_${grid}_${pes}_${soptions}.${testid}

To build and run a test, the process is the same as a case,
  cd into the test directory,
  
  run icepack.build
  
  run icepack.submit

The test results will be generated in a local file called **test_output**.

When running a test suite, the **icepack.create.case** command line automatically generates all the tests
under a directory named ${test_suite}.${testid}.  It then automatically builds and submits all
tests.  When the tests are complete, run the **results.csh** script to see the results from all the
tests.

Tests are defined under **configuration/scripts/tests/**.  The tests currently supported are:
  smoke   - Runs the model for default length.  The length and options can
            be set with the -s command line option.  The test passes if the
            model completes successfully.
  restart - Runs the model for 14 months, writing a restart file at month 3 and
            again at the end of the run.  Runs the model a second time starting from the
            month 3 restart and writing a restart at month 12 of the model run.
            The test passes if both runs complete and
            if the restart files at month 12 from both runs are bit-for-bit identical.

Please run ``./icepack.create.case -h`` for additional details.

.. _additional:

Additional testing options
--------------------------

There are several additional options on the ``icepack.create.case`` command line for testing that
provide the ability to regression test and compare tests to each other.

  ``-bd`` defines a baseline directory where tests can be stored for regression testing
  
  ``-bg`` defines a version name that where the current tests can be saved for regression testing
  
  ``-bc`` defines a version name that the current tests should be compared to for regression testing
  
  ``-td`` provides a way to compare tests with each other

To use ``-bg``,
  ``> icepack.create.case -ts base_suite -m wolf -testid v1 -bg version1 -bd $SCRATCH/ICEPACK_BASELINES``
  will copy all the results from the test suite to ``$SCRATCH/ICEPACK_BASELINES/version1``.

To use ``-bc``,
  ``> icepack.create.case -ts base_suite -m wolf -testid v2 -bc version1 -bd $SCRATCH/ICEPACK_BASELINES``
  will compare all the results from this test suite to results saved before in $SCRATCH/ICEPACK_BASELINES/version1``.

``-bc`` and ``-bg`` can be combined,
  ``>icepack.create.case -ts base_suite -m wolf -testid v2 -bg version2 -bc version1 -bd $SCRATCH/ICEPACK_BASELINES``
  will save the current results to ``$SCRATCH/ICEPACK_BASELINES/version2`` and compare the current results to
  results save before in ``$SCRATCH/ICEPACK_BASELINES/version1``.

``-bg``, ``-bc``, and ``-bd`` are used for regression testing.  There is a default ``-bd`` on each machine.

``-td`` allows a user to compare one test result to another.  For instance,

CHECK provide example suitable for Icepack. This one doesn't work because it relies on MPI

  ``> icepack.create.case -t smoke -m wolf -s run5day -testid t01``

  ``> icepack.create.case -t smoke -m wolf -s run5day -testid t01 -td smoke_gx3_8x2_run5day``

  An additional check will be done for the second test (because of the ``-td`` argument), and it will compare
  the output from the first test "smoke_gx3_8x2_run5day" to the output from its test "smoke_gx3_4x2_run5day"
  and generate a result for that.  It's important that the first test complete before the second test is done.  Also, the ``-td`` option works only if the testid and the machine are the same for the baseline run and the current run.

.. _format:

Test suite format
-----------------

The format for the test suite file is relatively simple.  It is a text file with white space delimited 
columns, e.g. **base\_suite.ts**

.. _tab-test:

.. csv-table:: Table 7
   :header: "#Test", "Grid", "PEs", "Sets", "BFB-compare"
   :widths: 7, 7, 7, 15, 15

   "smoke", "col", "1x1", "diag1,run1year", ""
   "smoke", "col", "1x1", "debug,run1year", ""
   "restart", "col", "1x1", "debug", ""
   "restart", "col", "1x1", "diag1", ""
   "restart", "col", "1x1", "pondcesm", ""
   "restart", "col", "1x1", "pondlvl", ""
   "restart", "col", "1x1", "pondtopo", ""


The first column is the test name, the second the grid, the third the pe count, the fourth column is
the ``-s`` options and the fifth column is the ``-td`` argument. (The grid and PEs columns are provided for compatibility with the similar CICE scripts.)  The fourth and fifth columns are optional.
The argument to ``-ts`` defines which filename to choose and that argument can contain a path.  ``icepack.create.case`` 
will also look for the filename in **configuration/scripts/tests/** where some preset test suites are defined.

Example Tests (Quickstart)
--------------------------

To generate a baseline dataset for a test case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``./icepack.create.case -t smoke -m wolf -bg icepackv6.0.0 -testid t00``

``cd wolf_smoke_col_1x1.t00``

``./icepack.build``

``./icepack.submit``

After job finishes, check output

``cat test_output``


To run a test case and compare to a baseline dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``./icepack.create.case -t smoke -m wolf -bc icepackv6.0.0 -testid t01``

``cd wolf_smoke_col_1x1.t01``

``./icepack.build``

``./icepack.submit``

After job finishes, check output

``cat test_output``


To run a test suite to generate baseline data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``./icepack.create.case -m wolf -ts base_suite -testid t02 -bg icepackv6.0.0bs``

``cd base_suite.t02``

Once all jobs finish, concatenate all output

``./results.csh``   All tests results will be stored in results.log

To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"

``./timeseries.csh <directory>``

``ls \*.png``


To run a test suite to compare to baseline data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``./icepack.create.case -m wolf -ts base_suite -testid t03 -bc icepackv6.0.0bs``

``cd base_suite.t03``

Once all jobs finish, concatenate all output

``./results.csh``   All tests results will be stored in results.log

To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"

``./timeseries.csh <directory>``

``ls \*.png``


To compare to another test
~~~~~~~~~~~~~~~~~~~~~~~~~~

CHECK needs a different example for Icepack

`First:`

./icepack.create.case -m wolf -t smoke -testid t01 -p 8x2

cd wolf_smoke_gx3_8x2.t01

./icepack.build

./icepack.submit

# After job finishes, check output

cat test_output

`Then, do the comparison:` 

./icepack.create.case -m wolf -t smoke -testid t01 -td smoke_gx3_8x2 -s thread -p 4x1

cd wolf_smoke_gx3_4x1_thread.t01

./icepack.build

./icepack.submit

# After job finishes, check output

cat test_output


Additional Details
------------------
- In general, the baseline generation, baseline compare, and test diff are independent.
- Use the ``-bd`` flag to specify the location where you want the baseline dataset
    to be written.  Without specifying ``-bd``, the baseline dataset will be written
    to the default baseline directory found in the **env.<machine>** file (``ICE_MACHINE_BASELINE``).
- If ``-bd`` is not passed, the scripts will look for baseline datasets in the default 
    baseline directory found in the **env.<machine>** file (``ICE_MACHINE_BASELINE``).
    If the ``-bd`` option is passed, the scripts will look for baseline datasets in the
    location passed to the ``-bd`` argument.
- To generate a baseline dataset for a specific version (for regression testing),
    use ``-bg <version_name>``.  The scripts will then place the baseline dataset
    in ``$ICE_MACHINE_BASELINE/<version_name>/``
- The ``-testid`` flag allows users to specify a testing id that will be added to the
    end of the case directory.  For example, 
    ``./icepack.create.case -m wolf -t smoke -testid t12``
    creates the directory **wolf_smoke_col_1x1.t12**.  This flag is REQUIRED if using ``-t`` or ``-ts``.


.. _tabnamelist:

Table of namelist options
=========================

CHECK

.. _tab-namelist:

.. csv-table:: Table 7
   :header: "variable", "options/format", "description", "recommended value"
   :widths: 15, 15, 30, 15 

   "*setup_nml*", "", "", ""
   "", "", "*Time, Diagnostics*", ""
   "``days_per_year``", "``360`` or ``365``", "number of days in a model year", "365"
   "``use_leap_years``", "true/false", "if true, include leap days", ""
   "``year_init``", "yyyy", "the initial year, if not using restart", ""
   "``istep0``", "integer", "initial time step number", "0"
   "``dt``", "seconds", "thermodynamics time step length", "3600."
   "``npt``", "integer", "total number of time steps to take", ""
   "``ndtd``", "integer", "number of dynamics/advection/ridging/steps per thermo timestep", "1"
   "", "", "*Initialization/Restarting*", ""
   "``runtype``", "``initial``", "start from ``ice_ic``", ""
   "", "``continue``", "restart using ``pointer_file``", ""
   "``ice_ic``", "``default``", "latitude and sst dependent", "default"
   "", "``none``", "no ice", ""
   "", "path/file", "restart file name", ""
   "``restart``", "true/false", "initialize using restart file", "``.true.``"
   "``use_restart_time``", "true/false", "set initial date using restart file", "``.true.``"
   "``restart_format``", "nc", "read/write  restart files (use with PIO)", ""
   "", "bin", "read/write binary restart files", ""
   "``lcdf64``", "true/false", "if true, use 64-bit  format", ""
   "``restart_dir``", "path/", "path to restart directory", ""
   "``restart_ext``", "true/false", "read/write halo cells in restart files", ""
   "``restart_file``", "filename prefix", "output file for restart dump", "‘iced’"
   "``pointer_file``", "pointer filename", "contains restart filename", ""
   "``dumpfreq``", "``y``", "write restart every ``dumpfreq_n`` years", "y"
   "", "``m``", "write restart every ``dumpfreq_n`` months", ""
   "", "``d``", "write restart every ``dumpfreq_n`` days", ""
   "``dumpfreq_n``", "integer", "frequency restart data is written", "1"
   "``dump_last``", "true/false", "if true, write restart on last time step of simulation", ""
   "", "", "*Model Output*", ""
   "``bfbflag``", "true/false", "for bit-for-bit diagnostic output", ""
   "``diagfreq``", "integer", "frequency of diagnostic output in ``dt``", "24"
   "", "*e.g.*, 10", "once every 10 time steps", ""
   "``diag_type``", "``stdout``", "write diagnostic output to stdout", ""
   "", "``file``", "write diagnostic output to file", ""
   "``diag_file``", "filename", "diagnostic output file (script may reset)", ""
   "``print_global``", "true/false", "print diagnostic data, global sums", "``.false.``"
   "``print_points``", "true/false", "print diagnostic data for two grid points", "``.false.``"
   "``latpnt``", "real", "latitude of (2) diagnostic points", "" 
   "``lonpnt``", "real", "longitude of (2) diagnostic points", ""
   "``dbug``", "true/false", "if true, write extra diagnostics", "``.false.``"
   "``histfreq``", "string array", "defines output frequencies", ""
   "", "``y``", "write history every ``histfreq_n`` years", ""
   "", "``m``", "write history every ``histfreq_n`` months", ""
   "", "``d``", "write history every ``histfreq_n`` days", ""
   "", "``h``", "write history every ``histfreq_n`` hours", ""
   "", "``1``", "write history every time step", ""
   "", "``x``", "unused frequency stream (not written)", ""
   "``histfreq_n``", "integer array", "frequency history output is written", ""
   "", "0", "do not write to history", ""
   "``hist_avg``", "true", "write time-averaged data", "``.true.``"
   "", "false", "write snapshots of data", ""
   "``history\_dir``", "path/", "path to history output directory", ""
   "``history\_file``", "filename prefix", "output file for history", "‘iceh’"
   "``write\_ic``", "true/false", "write initial condition", ""
   "``incond\_dir``", "path/", "path to initial condition directory", ""
   "``incond\_file``", "filename prefix", "output file for initial condition", "‘iceh’"
   "``runid``", "string", "label for run (currently CESM only)", ""
   "", "", "", ""
   "*grid_nml*", "", "", ""
   "", "", "*Grid*", ""
   "``grid_format``", "``nc``", "read  grid and kmt files", "‘bin’"
   "", "``bin``", "read direct access, binary file", ""
   "``grid_type``", "``rectangular``", "defined in *rectgrid*", ""
   "", "``displaced_pole``", "read from file in *popgrid*", ""
   "", "``tripole``", "read from file in *popgrid*", ""
   "", "``regional``", "read from file in *popgrid*", ""
   "``grid_file``", "filename", "name of grid file to be read", "‘grid’"
   "``kmt_file``", "filename", "name of land mask file to be read", "‘kmt’"
   "``gridcpl_file``", "filename", "input file for coupling grid info", ""
   "``kcatbound``", "``0``", "original category boundary formula", "0"
   "", "``1``", "new formula with round numbers", ""
   "", "``2``", "WMO standard categories", ""
   "", "``-1``", "one category", ""
   "", "", "", ""
   "*domain_nml*", "", "", ""
   "", "", "*Domain*", ""
   "``nprocs``", "integer", "number of processors to use", ""
   "``processor_shape``", "``slenderX1``", "1 processor in the y direction (tall, thin)", ""
   "", "``slenderX2``", "2 processors in the y direction (thin)", ""
   "", "``square-ice``", "more processors in x than y, :math:`\sim` square", ""
   "", "``square-pop``", "more processors in y than x, :math:`\sim` square", ""
   "``distribution_type``", "``cartesian``", "distribute blocks in 2D Cartesian array", ""
   "", "``roundrobin``", "1 block per proc until blocks are used", ""
   "", "``sectcart``", "blocks distributed to domain quadrants", ""
   "", "``sectrobin``", "several blocks per proc until used", ""
   "", "``rake``", "redistribute blocks among neighbors", ""
   "", "``spacecurve``", "distribute blocks via space-filling curves", ""
   "``distribution_weight``", "``block``", "full block size sets ``work_per_block``", ""
   "", "``latitude``", "latitude/ocean sets ``work_per_block``", ""
   "``ew_boundary_type``", "``cyclic``", "periodic boundary conditions in x-direction", ""
   "", "``open``", "Dirichlet boundary conditions in x", ""
   "``ns_boundary_type``", "``cyclic``", "periodic boundary conditions in y-direction", ""
   "", "``open``", "Dirichlet boundary conditions in y", ""
   "", "``tripole``", "U-fold tripole boundary conditions in y", ""
   "", "``tripoleT``", "T-fold tripole boundary conditions in y", ""
   "``maskhalo_dyn``", "true/false", "mask unused halo cells for dynamics", ""
   "``maskhalo_remap``", "true/false", "mask unused halo cells for transport", ""
   "``maskhalo_bound``", "true/false", "mask unused halo cells for boundary updates", ""
   "", "", "", ""
   "*tracer_nml*", "", "", ""
   "", "", "*Tracers*", ""
   "``tr_iage``", "true/false", "ice age", ""
   "``restart_age``", "true/false", "restart tracer values from file", ""
   "``tr_FY``", "true/false", "first-year ice area", ""
   "``restart_FY``", "true/false", "restart tracer values from file", ""
   "``tr_lvl``", "true/false", "level ice area and volume", ""
   "``restart_lvl``", "true/false", "restart tracer values from file", ""
   "``tr_pond_cesm``", "true/false", "CESM melt ponds", ""
   "``restart_pond_cesm``", "true/false", "restart tracer values from file", ""
   "``tr_pond_topo``", "true/false", "topo melt ponds", ""
   "``restart_pond_topo``", "true/false", "restart tracer values from file", ""
   "``tr_pond_lvl``", "true/false", "level-ice melt ponds", ""
   "``restart_pond_lvl``", "true/false", "restart tracer values from file", ""
   "``tr_aero``", "true/false", "aerosols", ""
   "``restart_aero``", "true/false", "restart tracer values from file", ""
   "*thermo_nml*", "", "", ""
   "", "", "*Thermodynamics*", ""
   "``kitd``", "``0``", "delta function ITD approximation", "1"
   "", "``1``", "linear remapping ITD approximation", ""
   "``ktherm``", "``0``", "zero-layer thermodynamic model", ""
   "", "``1``", "Bitz and Lipscomb thermodynamic model", ""
   "", "``2``", "mushy-layer thermodynamic model", ""
   "``conduct``", "``MU71``", "conductivity :cite:`MU71`", ""
   "", "``bubbly``", "conductivity :cite:`PETB07`", ""
   "``a_rapid_mode``", "real", "brine channel diameter", "0.5x10 :math:`^{-3}` m"
   "``Rac_rapid_mode``", "real", "critical Rayleigh number", "10"
   "``aspect_rapid_mode``", "real", "brine convection aspect ratio", "1"
   "``dSdt_slow_mode``", "real", "drainage strength parameter", "-1.5x10 :math:`^{-7}` m/s/K"
   "``phi_c_slow_mode``", ":math:`0<\phi_c < 1`", "critical liquid fraction", "0.05"
   "``phi_i_mushy``", ":math:`0<\phi_i < 1`", "solid fraction at lower boundary", "0.85"
   "", "", "", ""
   "*dynamics_nml*", "", "", ""
   "", "", "*Dynamics*", ""
   "``kdyn``", "``0``", "dynamics OFF", "1"
   "", "``1``", "EVP dynamics", ""
   "", "``2``", "EAP dynamics", ""
   "``revised_evp``", "true/false", "use revised EVP formulation", ""
   "``ndte``", "integer", "number of EVP subcycles", "120"
   "``advection``", "``remap``", "linear remapping advection", "‘remap’"
   "", "``upwind``", "donor cell advection", ""
   "``kstrength``", "``0``", "ice strength formulation :cite:`Hibler79`", "1"
   "", "``1``", "ice strength formulation :cite:`Rothrock75`", ""
   "``krdg_partic``", "``0``", "old ridging participation function", "1"
   "", "``1``", "new ridging participation function", ""
   "``krdg_redist``", "``0``", "old ridging redistribution function", "1"
   "", "``1``", "new ridging redistribution function", ""
   "``mu_rdg``", "real", "e-folding scale of ridged ice", ""
   "``Cf``", "real", "ratio of ridging work to PE change in ridging", "17."
   "", "", "", ""
   "*shortwave_nml*", "", "", ""
   "", "", "*Shortwave*", ""
   "``shortwave``", "``default``", "NCAR CCSM3 distribution method", ""
   "", "``dEdd``", "Delta-Eddington method", ""
   "``albedo_type``", "``default``", "NCAR CCSM3 albedos", "‘default’"
   "", "``constant``", "four constant albedos", ""
   "``albicev``", ":math:`0<\alpha <1`", "visible ice albedo for thicker ice", ""
   "``albicei``", ":math:`0<\alpha <1`", "near infrared ice albedo for thicker ice", ""
   "``albsnowv``", ":math:`0<\alpha <1`", "visible, cold snow albedo", ""
   "``albsnowi``", ":math:`0<\alpha <1`", "near infrared, cold snow albedo", ""
   "``ahmax``", "real", "albedo is constant above this thickness", "0.3 m"
   "``R_ice``", "real", "tuning parameter for sea ice albedo from Delta-Eddington shortwave", ""
   "``R_pnd``", "real", "... for ponded sea ice albedo …", ""
   "``R_snw``", "real", "... for snow (broadband albedo) …", ""
   "``dT_mlt``", "real", ":math:`\Delta` temperature per :math:`\Delta` snow grain radius", ""
   "``rsnw_mlt``", "real", "maximum melting snow grain radius", ""
   "``kalg``", "real", "absorption coefficient for algae", ""
   "", "", "", ""
   "*ponds_nml*", "", "", ""
   "", "", "*Melt Ponds*", ""
   "``hp1``", "real", "critical ice lid thickness for topo ponds", "0.01 m"
   "``hs0``", "real", "snow depth of transition to bare sea ice", "0.03 m"
   "``hs1``", "real", "snow depth of transition to pond ice", "0.03 m"
   "``dpscale``", "real", "time scale for flushing in permeable ice", ":math:`1\times 10^{-3}`"
   "``frzpnd``", "``hlid``", "Stefan refreezing with pond ice thickness", "‘hlid’"
   "", "``cesm``", "CESM refreezing empirical formula", ""
   "``rfracmin``", ":math:`0 \le r_{min} \le 1`", "minimum melt water added to ponds", "0.15"
   "``rfracmax``", ":math:`0 \le r_{max} \le 1`", "maximum melt water added to ponds", "1.0"
   "``pndaspect``", "real", "aspect ratio of pond changes (depth:area)", "0.8"
   "", "", "", ""
   "*zbgc_nml*", "", "", ""
   "", "", "*Biogeochemistry*", ""
   "``tr_brine``", "true/false", "brine height tracer", ""
   "``restart_hbrine``", "true/false", "restart tracer values from file", ""
   "``skl_bgc``", "true/false", "biogeochemistry", ""
   "``bgc_flux_type``", "``Jin2006``", "ice–ocean flux velocity of :cite:`JDWSTWLG06`", ""
   "", "``constant``", "constant ice–ocean flux velocity", ""
   "``restart_bgc``", "true/false", "restart tracer values from file", ""
   "``restore_bgc``", "true/false", "restore nitrate/silicate to data", ""
   "``bgc_data_dir``", "path/", "data directory for bgc", ""
   "``sil_data_type``", "``default``", "default forcing value for silicate", ""
   "", "``clim``", "silicate forcing from ocean climatology :cite:`GLBA06`", ""
   "``nit_data_type``", "``default``", "default forcing value for nitrate", ""
   "", "``clim``", "nitrate forcing from ocean climatology :cite:`GLBA06`", ""
   "", "``sss``", "nitrate forcing equals salinity", ""
   "``tr_bgc_C_sk``", "true/false", "algal carbon tracer", ""
   "``tr_bgc_chl_sk``", "true/false", "algal chlorophyll tracer", ""
   "``tr_bgc_Am_sk``", "true/false", "ammonium tracer", ""
   "``tr_bgc_Sil_sk``", "true/false", "silicate tracer", ""
   "``tr_bgc_DMSPp_sk``", "true/false", "particulate DMSP tracer", ""
   "``tr_bgc_DMSPd_sk``", "true/false", "dissolved DMSP tracer", ""
   "``tr_bgc_DMS_sk``", "true/false", "DMS tracer", ""
   "``phi_snow``", "real", "snow porosity for brine height tracer", ""
   "", "", "", ""
   "*forcing_nml*", "", "", ""
   "", "", "*Forcing*", ""
   "``formdrag``", "true/false", "calculate form drag", ""
   "``atmbndy``", "``default``", "stability-based boundary layer", "‘default’"
   "", "``constant``", "bulk transfer coefficients", ""
   "``fyear_init``", "yyyy", "first year of atmospheric forcing data", ""
   "``ycycle``", "integer", "number of years in forcing data cycle", ""
   "``atm_data_format``", "``nc``", "read  atmo forcing files", ""
   "", "``bin``", "read direct access, binary files", ""
   "``atm_data_type``", "``default``", "constant values defined in the code", ""
   "", "``LYq``", "AOMIP/Large-Yeager forcing data", ""
   "", "``monthly``", "monthly forcing data", ""
   "", "``ncar``", "NCAR bulk forcing data", ""
   "", "``oned``", "column forcing data", ""
   "``atm_data_dir``", "path/", "path to atmospheric forcing data directory", ""
   "``calc_strair``", "true", "calculate wind stress and speed", ""
   "", "false", "read wind stress and speed from files", ""
   "``highfreq``", "true/false", "high-frequency atmo coupling", ""
   "``natmiter``", "integer", "number of atmo boundary layer iterations", ""
   "``calc_Tsfc``", "true/false", "calculate surface temperature", "``.true.``"
   "``precip_units``", "``mks``", "liquid precipitation data units", ""
   "", "``mm_per_month``", "", ""
   "", "``mm_per_sec``", "(same as MKS units)", ""
   "``tfrz_option``", "``minus1p8``", "constant ocean freezing temperature (:math:`-1.8^\circ C`)", ""
   "", "``linear_salt``", "linear function of salinity (ktherm=1)", ""
   "", "``mushy_layer``", "matches mushy-layer thermo (ktherm=2)", ""
   "``ustar_min``", "real", "minimum value of ocean friction velocity", "0.0005 m/s"
   "``fbot_xfer_type``", "``constant``", "constant ocean heat transfer coefficient", ""
   "", "``Cdn\_ocn``", "variable ocean heat transfer coefficient", ""
   "``update_ocn_f``", "true", "include frazil water/salt fluxes in ocn fluxes", ""
   "", "false", "do not include (when coupling with POP)", ""
   "``l_mpond_fresh``", "true", "retain (topo) pond water until ponds drain", ""
   "", "false", "release (topo) pond water immediately to ocean", ""
   "``oceanmixed_ice``", "true/false", "active ocean mixed layer calculation", "``.true.`` (if uncoupled)"
   "``ocn_data_format``", "``nc``", "read  ocean forcing files", ""
   "", "``bin``", "read direct access, binary files", ""
   "``sss_data_type``", "``default``", "constant values defined in the code", ""
   "", "``clim``", "climatological data", ""
   "", "``near``", "POP ocean forcing data", ""
   "``sst_data_type``", "``default``", "constant values defined in the code", ""
   "", "``clim``", "climatological data", ""
   "", "``ncar``", "POP ocean forcing data", ""
   "``ocn_data_dir``", "path/", "path to oceanic forcing data directory", ""
   "``oceanmixed_file``", "filename", "data file containing ocean forcing data", ""
   "``restore_sst``", "true/false", "restore sst to data", ""
   "``trestore``", "integer", "sst restoring time scale (days)", ""
   "``restore_ice``", "true/false", "restore ice state along lateral boundaries", ""
   "", "", "", ""
   "*icefields_tracer_nml*", "", "", ""
   "", "", "*History Fields*", ""
   "``f_<var>``", "string", "frequency units for writing ``<var>`` to history", ""
   "", "``y``", "write history every ``histfreq_n`` years", ""
   "", "``m``", "write history every ``histfreq_n`` months", ""
   "", "``d``", "write history every ``histfreq_n`` days", ""
   "", "``h``", "write history every ``histfreq_n`` hours", ""
   "", "``1``", "write history every time step", ""
   "", "``x``", "do not write ``<var>`` to history", ""
   "", "``md``", "*e.g.,* write both monthly and daily files", ""
   "``f_<var>_ai``", "", "grid cell average of ``<var>`` (:math:`\times a_i`)", ""

