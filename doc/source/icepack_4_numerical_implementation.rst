**********
User Guide
**********

Numerical implementation
========================

Icepack is written in FORTRAN90 and runs on platforms using UNIX, LINUX,
and other operating systems. The code is not parallelized. (CHANGE IF OPENMP IS IMPLEMENTED)

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
RM? **README\_v1**
    basic information and pointers

**columnphysics/**
    the essential physics code

CHANGE    **constants**

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
        mechanical redistribution component (ridging)

    **icepack\_meltpond\_cesm.F90**
        CESM melt pond parameterization

    **icepack\_meltpond\_lvl.F90**
        level-ice melt pond parameterization

    **icepack\_meltpond\_topo.F90**
        topo melt pond parameterization

CHECK THIS    **icepack\_ocean.F90**
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
        vertical salinity parameterization of (CITE Jeffery)

**configuration/**
    drivers and scripts for testing Icepack in stand-alone mode
    **drivers/**
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

CHECK THIS        **icepack\_drv\_constants.F90**
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

RENAME bgc        **icepack\_drv\_restart\_column.F90**
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
CHECK            namelist input data (data paths depend on particular system)

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
ADD LINK            scripts for configuring and running basic tests

**doc/**
    documentation

**icepack.create.case**
    main script for setting up a test case

ADD LINK 
A case (compile) directory is created upon initial execution of the script 
**icepack.create.case** at the user-specified location provided after the -c flag. 
Executing the command ``./icepack.create.case -h`` provides helpful information for 
this tool. Please refer to the (LINK)user guide for further information.


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
spaced at :math:`1/n_b` intervals beginning with `bgrid(2)` :math:` =
1/(2n_b)`. The ``igrid`` interior points :math:`[2, n_b]` are also
equidistant with the same spacing, but physically coincide with points
midway between those of ``bgrid``.


.. _testconfigs:

Test configurations
-------------------

*UPDATE*

The column is located
near Barrow (71.35N, 156.5W). Options for choosing the column
configuration are given in **comp\_ice** (choose `RES col`) and in the
namelist file, **input\_templates/col/ice\_in**. Here, ``istep0`` and the
initial conditions are set such that the run begins September 1 with no
ice. 


.. _init:

Initialization and coupling
---------------------------

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
will be described in more detail in :ref:`tabnamelist`.

Two namelist variables control model initialization, ``ice\_ic``
and ``restart``.  Setting ``ice\_ic`` = default causes the model to run using
constant forcing and initial values set in the code.  To start
from a file **filename**, set 
``restart`` = true and ``ice\_ic`` = **filename**.  When restarting using the Icepack
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

(LINK)
CICE provides extensive options for model output, including many derived output variables.

Diagnostic files
~~~~~~~~~~~~~~~~

Icepack writes diagnostic information for each grid cell as a separate file, 
**ice\_diag.\***.


Restart files
~~~~~~~~~~~~~

CHANGE as needed re netCDF

CICE provides restart data in binary unformatted or netCDF formats, via
the `IO\_TYPE` flag in **comp\_ice** and namelist variable
`restart\_format`. 

The restart files created by the Icepack driver contain all of the variables needed
for a full, exact restart. The filename begins with the character string
‘iced.’, and the restart dump frequency is given by the namelist
variable `dumpfreq`. The namelist variable `ice\_ic` contains the
pointer to the filename from which the restart data is to be read.


Execution procedures
====================

Quick-start instructions are provided in section :ref:`quickstart`

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

  -s are comma separated optional env or namelist settings (default is "null")

  -t is the test name and location (cannot be used with -c).

  -bd is used to specify the location of the baseline datasets (only used with -t)

  -bg is used to specify the icepack version name for generating baseline datasets (only used with -t)

  -bc is used to specify the icepack version name for comparison. I.e., the version name for the baseline dataset (only used with -t)

  -testid is used to specify a test ID (used only with -t or -ts)

  -ts is used to generate all test cases for a given test suite.


Several files are placed in the case directory

 - **env.**[machine] defines the environment

 - **icepack.settings** defines many variables associated with building and running the model

 - **makdep.c** is a tool that will automatically generate the make dependencies

 - **Macros.**[machine] defines the Makefile macros

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
env variable ``ICE\_CLEANBUILD`` in **icepack.settings**.  

The **icepack.submit** script simply submits the **icepack.run script**.  
You can also submit the **icepack.run** script on the command line.

To port, an **env.**[machine] and **Macros.**[machine] file have to be added to 
**configuration/scripts/machines/** and the 
**icepack.run.setup.csh** file needs to be modified.
 - cd to **configuration/scripts/machines/**
 - Copy existing env and Macros files to new names for your new machine
 - Edit the env and Macros files
 - cd .. to **configuration/scripts/**
 - Edit the **icepack.run.setup.csh** script to add a section for your machine 
with batch settings and job launch settings
 - Download and untar a forcing dataset to the location defined by 
``ICE\_MACHINE\_INPUTDATA`` in the env file
 - Create a file in your home directory called .cice\_proj and add your preferred account name to the first line.


Directories
-------------

CHECK

The **icepack.create.case** script creates a case directory in the location specified 
by the ``-c`` or ``-t`` flags.  The **icepack.build** (or equivalent test suite) script 
creates the run directory defined by the env variable ``ICE\_RUNDIR`` in 
**icepack.settings**, and it compiles the code there.  The run directory is further 
populated by the **icepack.run** script, which also runs the executable.

Build and run logs will be copied from the run directory into the case **logs/** 
directory when complete.


Local modifications
--------------------------

If there are problems, you can manually edit 
the env, Macros, and **icepack.run** files in the case directory until things are 
working properly.  Then you can copy the env and Macros files back to 
**configuration/scripts/machines**.  You will have to manually modify the 
**icepack.run.setup.csh** script if there any changes needed there.

You can also directly modify the namelist files (**icepack\_in**) in the run directory and
run the code by submitting the executable **icepack** directly.  Beware that any changes 
make in the run directory will be overwritten if scripts are later run from the case
directory.

Forcing data
------------

FINISH:

The code is currently configured to run in standalone mode on a 4-cell grid using 
atmospheric data, available as detailed on the `wiki <https://github.com/CICE-Consortium/Icepack/wiki/Testing-Icepack>`_.
These data files are designed only for testing the code, not for use in production 
runs or as observational data.  Please do not publish results based on these data
sets.  Module **configuration/driver/icepack\_drv\_forcing.F90**
can be modified to change the forcing data. 



Adding Tracers
====================

.. _addtrcr:

Tracers added to Icepack will also require extensive modifications to the host
sea ice model, including initialization on the horizontal grid, namelist flags 
and restart capabilities.  Modifications to the Icepack driver should reflect
the modifications needed in the host model but are not expected to match completely.
We recommend that the logical namelist variable
``tr\_[tracer]`` be used for all calls involving the new tracer outside of
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
the conserved “tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies, and code must be included to account for new
dependencies (e.g., :math:`a_{lvl}` and :math:`a_{pnd}` in
**ice\_itd.F90** and **ice\_mechred.F90**).

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **icepack\_drv\_domain\_size.F90**: increase ``max\_ntrcr`` (can also add option
   to **icepack.settings** and **icepack.build**)

#. **icepack\_drv\_state.F90**: declare `nt\_[tracer]` and `tr\_[tracer]`

#. **icepack\_[tracer].F90**: create initialization, physics routines

#. **ice\_drv\_init.F90**: (some of this may be done in **ice\_[tracer].F90**
   instead)

   -  add new module and ``tr\_[tracer]`` to list of used modules and
      variables

   -  add logical namelist variable ``tr\_[tracer]``

   -  initialize namelist variable

   -  print namelist variable to diagnostic output file

   -  increment number of tracers in use based on namelist input (``ntrcr``)

   -  define tracer types (``trcr\_depend`` = 0 for ice area tracers, 1 for
      ice volume, 2 for snow volume, 2+``nt\_``[tracer] for dependence on
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

