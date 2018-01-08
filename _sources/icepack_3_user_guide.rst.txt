:tocdepth: 4

**********
User Guide
**********

Implementation
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

The present code distribution includes source code for the column physics,
source code for the driver, and the scripts.  Forcing data is available from the ftp site.
The directory structure of Icepack is as follows.  All columnphysics filename have a prefix
of icepack\_ and all driver files are prefixed with icedrv\_*.

**LICENSE.pdf**
  license and policy for using and sharing the code

**DistributionPolicy.pdf**
  license and policy for using and sharing the code

**README.md**
  basic information and pointers

**columnphysics/**
   directory for columnphysics source code, see :ref:`dev_colphys`

**configuration/scripts/**
   directory of support scripts, see :ref:`dev_scripts`

**configuration/driver/**
   directory of icepack driver code, see :ref:`dev_driver`

**doc/**
    documentation

**icepack.create.case**
  main icepack script for creating cases

A case (compile) directory is created upon initial execution of the script 
**icepack.create.case** at the user-specified location provided after the -c flag. 
Executing the command ``./icepack.create.case -h`` provides helpful information for 
this tool.

.. _grids:

Grid and boundary conditions 
----------------------------

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

Initialization and Forcing
--------------------------

Icepack’s parameters and variables are initialized in several
steps. Many constants and physical parameters are set in
**icepack\_constants.F90**. In the current driver implementation,
a namelist file is read to setup the model.
Namelist values are given default
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
initial values set in the code.  To start
from a file **filename**, set 
``restart`` = .true. and ``ice_ic`` = **filename**.  When restarting using the Icepack
driver, for simplicity the tracers are assumed to be set the same way (on/off) as in the
run that created the restart file; i.e. that the restart file contains exactly the 
information needed for the new run.  CICE is more flexible in this regard.

For stand-alone runs,
routines in **icedrv\_forcing.F90** read and interpolate data from files,
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


.. _history:

Model output
------------

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


.. _bgc-hist:

Biogeochemistry History Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note:  History output is not provided with Icepack.  This documentation
indicates what is available for output and is implemented in CICE.

CHECK and FIX!

The biogeochemical history fields specified in icefields\_bgc\_nml are
written when ‘x’ is replaced with a time interval: step (‘1’), daily
(‘d’), monthly (‘m’), or yearly (‘y’). Several of these flags turn on
multiple history variables according to the particular ecosystem
prescribed in **icepack\_in**. For example, biogeochemical fluxes from the
ice to ocean will be saved monthly in the history output if

::

    f_fbio = 'm'

However, only the biogeochemical tracers which are active will be saved.
This includes at most fNit nitrate, fAm ammonium, fN algal nitrogen,
fDOC dissolved organic carbon, fDON dissolved organic nitrogen, fFep
particulate iron, fFed dissolved iron, fSil silicate, fhum humic matter,
fPON passive mobile tracer, fDMS DMS, fDMSPd dissolved DMSP and fDMSPp
particulate DMSP.

:ref:`tab-bio-history` lists the
biogeochemical tracer history flags along with a short description and
the variable or variables saved. Not listed are flags appended with
\_ai, i.e. f\_fbio\_ai. These fields are identical to their counterpart.
i.e. f\_fbio, except they are averaged by ice area.

:ref:`tab-bio-history` :*Biogeochemical History variables*

.. _tab-bio-history:

.. csv-table:: Table 5
   :header: "History Flag", "Definition", "Variable(s)", "Units"
   :widths: 10, 25, 20, 10

   "f\_faero\_atm", "atmospheric aerosol deposition flux", "faero\_atm", "kg m\ :math:`^{-2}` s\ :math:`^{-1}`"
   "f\_faero\_ocn", "aerosol flux from ice to ocean", "faero\_ocn", "kg m\ :math:`^{-2}` s\ :math:`^{-1}`"
   "f\_aero", "aerosol mass (snow and ice ssl and int)", "aerosnossl, aerosnoint,aeroicessl, aeroiceint", "kg/kg"
   "f\_fbio", "biological ice to ocean flux", "fN, fDOC, fNit, fAm,fDON,fFep\ :math:`^a`, fFed\ :math:`^a`, fSil,fhum, fPON, fDMSPd,fDMS, fDMSPp, fzaero", "mmol m\ :math:`^{-2}` s\ :math:`^{-1}`"
   "f\_zaero", "bulk z-aerosol mass fraction", "zaero", "kg/kg"
   "f\_bgc\_S", "bulk z-salinity", "bgc\_S", "ppt"
   "f\_bgc\_N", "bulk algal N concentration", "bgc\_N", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_C", "bulk algal C concentration", "bgc\_C", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DOC", "bulk DOC concentration", "bgc\_DOC", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DON", "bulk DON concentration", "bgc\_DON", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DIC", "bulk DIC concentration", "bgc\_DIC", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_chl", "bulk algal chlorophyll concentration", "bgc\_chl", "mg chl m\ :math:`^{-3}`"
   "f\_bgc\_Nit", "bulk nitrate concentration", "bgc\_Nit", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_Am", "bulk ammonium concentration", "bgc\_Am", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_Sil", "bulk silicate concentration", "bgc\_Sil", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DMSPp", "bulk particulate DMSP concentration", "bgc\_DMSPp", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DMSPd", "bulk dissolved DMSP concentration", "bgc\_DMSPd", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DMS", "bulk DMS concentration", "bgc\_DMS", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_Fe", "bulk dissolved and particulate iron conc.", "bgc\_Fed, bgc\_Fep", ":math:`\mu\,`\ mol m\ :math:`^{-3}`"
   "f\_bgc\_hum", "bulk humic matter concentration", "bgc\_hum", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_PON", "bulk passive mobile tracer conc.", "bgc\_PON", "mmol m\ :math:`^{-3}`"
   "f\_upNO", "Total algal :math:`{\mbox{NO$_3$}}` uptake rate", "upNO", "mmol m\ :math:`^{-2}` d\ :math:`^{-1}`"
   "f\_upNH", "Total algal :math:`{\mbox{NH$_4$}}` uptake rate", "upNH", "mmol m\ :math:`^{-2}` d\ :math:`^{-1}`"
   "f\_bgc\_ml", "upper ocean tracer concentrations", "ml\_N, ml\_DOC, ml\_Nit,ml\_Am, ml\_DON, ml\_Fep\ :math:`^b`,ml\_Fed\ :math:`^b`, ml\_Sil, ml\_hum, ml\_PON,ml\_DMS, ml\_DMSPd, ml\_DMSPp", "mmol m\ :math:`^{-3}`"
   "f\_bTin", "ice temperature on the bio grid", "bTizn", ":math:`^o`\ C"
   "f\_bphi", "ice porosity on the bio grid", "bphizn", "%"
   "f\_iDin", "brine eddy diffusivity on the interface bio grid", "iDin", "m\ :math:`^{2}` d\ :math:`^{-1}`"
   "f\_iki", "ice permeability on the interface bio grid", "ikin", "mm\ :math:`^{2}`"
   "f\_fbri", "ratio of brine tracer height to ice thickness", "fbrine", "1"
   "f\_hbri", "brine tracer height", "hbrine", "m"
   "f\_zfswin", "internal ice PAR on the interface bio grid", "zfswin", "W m\ :math:`^{-2}`"
   "f\_bionet", "brine height integrated tracer concentration", "algalN\_net, algalC\_net, chl\_net, pFe\ :math:`^c`\ \_net, dFe\ :math:`^c`\ \_net, Sil\_net, Nit\_net, Am\_net, hum\_net, PON\_net, DMS\_net, DMSPd\_net, DMSPp\_net, DOC\_net, zaero\_net, DON\_net", "mmol m\ :math:`^{-2}`"
   "f\_biosnow", snow integrated tracer concentration", "algalN\_snow, algalC\_snow,chl\_snow, pFe\ :math:`^c`\ \_snow, dFe\ :math:`^c`\ \_snow,Sil\_snow, Nit\_snow, Am\_snow, hum\_snow, PON\_snow, DMS\_snow, DMSPd\_snow, DMSPp\_snow, DOC\_snow, zaero\_snow, DON\_snow", "mmol m\ :math:`^{-2}`"
   "f\_grownet", "Net specific algal growth rate", "grow\_net", "m d\ :math:`^{-1}`"
   "f\_PPnet", "Net primary production", "PP\_net", "mgC m\ :math:`^{-2}` d\ :math:`^{-1}`"
   "f\_algalpeak", "interface bio grid level of peak chla", "peak\_loc", "1"
   "f\_zbgc\_frac", "mobile fraction of tracer", "algalN\_frac, chl\_frac, pFe\_frac,dFe\_frac, Sil\_frac, Nit\_frac,Am\_frac, hum\_frac, PON\_frac,DMS\_frac, DMSPd\_frac, DMSPp\_frac,DOC\_frac, zaero\_frac, DON\_frac", "1"


:math:`^a` units are :math:`\mu`\ mol m\ :math:`^{-2}` s\ :math:`^{-1}`

:math:`^b` units are :math:`\mu`\ mol m\ :math:`^{-3}`

:math:`^c` units are :math:`\mu`\ mol m\ :math:`^{-2}`


Running Icepack
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

.. _scripts:

Scripts
-------

The icepack scripts are written to allow quick setup of cases and tests.  Once a case is 
generated, users can manually modify the namelist and other files to custom configure
the case.  Several settings are available via scripts as well.

Most of the scripts that configure, build and run Icepack are contained in 
the directory **configuration/scripts/**, except for **icepack.create.case**, which is
in the main directory.  **icepack.create.case** is the main script that generates a case. 

Users may need to port the scripts to their local machine.
Specific instructions for porting are provided in :ref:`porting`.

``icepack.create.case -h`` will provide the latest information about how to use the tool.
There are three usage modes,

* ``-c`` creates individual stand alone cases.
* ``-t`` creates individual tests.  Tests are just cases that have some extra automation in order to carry out particular tests such as exact restart.
* ``-ts`` creates a test suite.  Test suites are predefined sets of tests and ``-ts`` provides the ability to quick setup, build, and run a full suite of tests.

All modes will require use of ``-m`` to specify the machine and case and test modes 
can use ``-s`` to define specific options.  ``-t`` and ``-ts`` will require ``-testid`` to be set 
and both of the test modes can use ``-bd``, ``-bg``, ``-bc``, and ``-td`` to compare with other results.
Testing will be described in greater detail in the :ref:`testing` section.

Again, ``icepack.create.case -h`` will show the latest usage information including 
the available ``-s`` options, the current ported machines, and the test choices.

To create a case, run **icepack.create.case**::

  icepack.create.case -c mycase -m machine
  cd mycase

Once a case/test is created, several files are placed in the case directory

- **env.[machine]** defines the environment
- **icepack.settings** defines many variables associated with building and running the model
- **makdep.c** is a tool that will automatically generate the make dependencies
- **Macros.[machine]** defines the Makefile macros
- **Makefile** is the makefile used to build the model
- **icepack.build** is a script that builds and compiles the model
- **icepack\_in** is the namelist input file
- **icepack.run** is a batch run script
- **icepack.submit** is a simple script that submits the icepack.run script

All scripts and namelist are fully resolved in the case.  Users can edit any
of the files in the case directory manually to change the model configuration,
build options, or batch settings.  The file
dependency is indicated in the above list.  For instance, if any of the files before
**icepack.build** in the list are edited, **icepack.build** should be rerun.

The **casescripts/** directory holds scripts used to create the case and can 
largely be ignored.  Once a case is created, the **icepack.build** script should be run
interactively and then the case should be submitted by executing the 
**icepack.submit** script interactively.  The **icepack.submit** script
simply submits the **icepack.run script**.  
You can also submit the **icepack.run** script on the command line.

Some hints:

- To change namelist, manually edit the **icepack_in** file
- To change batch settings, manually edit the top of the **icepack.run** file
- To turn on the debug compiler flags, set ``ICE_BLDDEBUG`` in **icepack.setttings** to true
- To change compiler options, manually edit the Macros file
- To clean the build before each compile, set ``ICE_CLEANBUILD`` in **icepack.settings** to true
  To not clean before the build, set ``ICE_CLEANBUILD`` in **icepack.settings** to false

To build and run::

  ./icepack.build
  ./icepack.submit

The build and run log files will be copied into the logs directory in the case directory.
Other model output will be in the run directory.  The run directory is set in **icepack.settings**
via the **ICE_RUNDIR** variable.  To modify the case setup, changes should be made in the
case directory, NOT the run directory.

.. _porting:

Porting
-------

To port, an **env.[machine]** and **Macros.[machine]** file have to be added to the
**configuration/scripts/machines/** directory and the 
**configuration/scripts/icepack.run.setup.csh** file needs to be modified.
 
- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files

- cd .. to **configuration/scripts/**

- Edit the **icepack.run.setup.csh** script to add a section for your machine 
  with batch settings and job launch settings

- Download and untar a forcing dataset to the location defined by 
  ``ICE_MACHINE_INPUTDATA`` in the env file

In fact, this process almost certainly will require some iteration.  The easiest way 
to carry this out is to create an initial set of changes as described above, then 
create a case and manually modify the **env.[machine]** file and **Macros.[machine]** 
file until the case can build and run.  Then copy the files from the case 
directory back to **configuration/scripts/machines/** and update 
the **configuration/scripts/icepack.run.setup.csh** file, retest, 
and then add and commit the updated machine files to the repository.

Machine Account Settings
~~~~~~~~~~~~~~~~~~~~~~~~

The machine account default is specified by the variable ``ICE_MACHINE_ACCT`` in 
the **env.[machine]** file.  The easiest way to change a user's default is to 
create a file in your home directory called **.cice\_proj** and add your 
preferred account name to the first line.  
There is also an option (``-a``) in **icepack.create.case** to define the account number.  
The order of precedent is **icepack.create.case** command line option, 
**.cice\_proj** setting, and then value in the **env.[machine]** file.

Forcing data
------------

CHECK once we've settled on a forcing suite:

The code is currently configured to run in standalone mode on a 4-cell grid using 
atmospheric data, available as detailed on the `wiki <https://github.com/CICE-Consortium/Icepack/wiki/Testing-Icepack>`_.
These data files are designed only for testing the code, not for use in production 
runs or as observational data.  Please do not publish results based on these data
sets.  Module **configuration/driver/icedrv\_forcing.F90**
can be modified to change the forcing data. 

The input data space is defined on a per machine basis by the ``ICE_MACHINE_INPUTDATA`` 
variable in the **env.[machine]** file.  That file space is often shared among multiple 
users, and it can be desirable to consider using a common file space with group read 
and write permissions such that a set of users can update the inputdata area as 
new datasets are available.


Run Directories
---------------

The **icepack.create.case** script creates a case directory.  However, the model 
is actually built and run under the ``ICE_OBJDIR`` and ``ICE_RUNDIR`` directories
as defined in the **icepack.settings** file.

Build and run logs will be copied from the run directory into the case **logs/** 
directory when complete.


Local modifications
-------------------

Scripts and other case settings can be changed manually in the case directory and
used.  Source code can be modified in the main sandbox.  When changes are made, the code
should be rebuilt before being resubmitted.  It is always recommended that users
modify the scripts and input settings in the case directory, NOT the run directory.
In general, files in the run directory are overwritten by versions in the case
directory when the model is built, submitted, and run.


.. _testing:

Testing Icepack
================

The section documents primarily how to use the Icepack scripts to carry 
out icepack testing.  Exactly what to test is a separate question and
depends on the kinds of code changes being made.  Prior to merging
changes to the CICE Consortium master, changes will be reviewed and
developers will need to provide a summary of the tests carried out.

There is a base suite of tests provided by default with Icepack and this
may be a good starting point for testing.


.. _indtests:

Individual Tests
----------------

The Icepack scripts support both setup of individual tests as well as test suites.  Individual
tests are run from the command line::

  ./icepack.create.case -t smoke -m wolf -s diag1,debug -testid myid -a P0000000

where ``-m`` designates a specific machine, ``-a`` designates the account number 
for the queue manager, and testid is a user defined string that allows
test cases to be uniquely identified.
The format of the case directory name for a test will always be 
``[machine]_[test]_[grid]_[pes]_[soptions].[testid]``

To build and run a test, the process is the same as a case.  cd to the 
test directory, run the build script, and run the submit script::

 cd [test_case]
 ./icepack.build
 ./icepack.submit

The test results will be copied into a local file called **test_output**.
To check those results::

 cat test_output

Tests are defined under **configuration/scripts/tests/**.  The tests currently supported are:

-  smoke   - Runs the model for default length.  The length and options can
            be set with the ``-s`` command line option.  The test passes if the
            model completes successfully.
-  restart - Runs the model for 14 months, writing a restart file at month 3 and
            again at the end of the run.  Runs the model a second time starting from the
            month 3 restart and writing a restart at month 12 of the model run.
            The test passes if both runs complete and
            if the restart files at month 12 from both runs are bit-for-bit identical.

Please run ``./icepack.create.case -h`` for the latest information.

.. _testsuites:

Test suites
-----------

Test suites are multiple tests that are specified via
an input file.  When invoking the test suite option (``-ts``) with **icepack.create.case**,
all tests will be created, built, and submitted automatically under
a directory called [suite_name].[testid]::

  ./icepack.create.case -ts base_suite -m wolf -testid myid -a P00000

Like an individual test, the ``-testid`` option must be specified and can be any 
string.  Once the tests are complete, results can be checked by running the
results.csh script in the [suite_name].[testid]::

  cd base_suite.[testid]
  ./results.csh

Please run ``./icepack.create.case -h`` for additional details.

The predefined test suites are defined under **configuration/scripts/tests** and the files defining 
the suites
have a suffix of .ts in that directory.  The format for the test suite file is relatively simple.  
It is a text file with white space delimited 
columns, e.g. **base_suite.ts**

.. _tab-test:

.. csv-table:: Table 6
   :header: "Test", "Grid", "PEs", "Sets", "BFB-compare"
   :widths: 7, 7, 7, 15, 15

   "smoke", "col", "1x1", "diag1", ""
   "smoke", "col", "1x1", "diag1,run1year", "smoke_col_1x1_diag1_run1year"
   "smoke", "col", "1x1", "debug,run1year", ""
   "restart", "col", "1x1", "debug", ""
   "restart", "col", "1x1", "diag1", ""
   "restart", "col", "1x1", "pondcesm", ""
   "restart", "col", "1x1", "pondlvl", ""
   "restart", "col", "1x1", "pondtopo", ""

The first column is the test name, the second the grid, the third the pe count, the fourth column is
the ``-s`` options and the fifth column is the ``-td`` argument. (The grid and PEs columns are provided 
for compatibility with the similar CICE scripts.)  The fourth and fifth columns are optional.
The argument to ``-ts`` defines which filename to choose and that argument can contain a path.  
**icepack.create.case** 
will look for the filename in the local directory, in **configuration/scripts/tests/**, or in the path defined
by the ``-ts`` option.

.. _regtesting:

Regression testing
------------------

The **icepack.create.case** options ``-bg``, ``-bc``, and ``-bd`` are used for regression testing.
There are several additional options on the **icepack.create.case** command line for testing that
provide the ability to regression test and compare tests to each other.  These options only
work in test (``-t``) or test suite (``-ts``) mode, not in case (``-c``) mode.

  ``-bd`` defines a top level directory where tests can be stored for regression testing.  The
  default is defined by ``ICE_MACHINE_BASELINE`` defined in ``env.[machine]``.
  
  ``-bg`` defines a directory name where the current tests can be saved for regression testing.  
  It's handy for this name to be related to the model version.  This directory will be created
  below the directory associated with ``-bd``.
  
  ``-bc`` defines the directory name that the current tests should be compared to for regression 
  testing.  This directory will be added to the directory associated with ``-bd``.
  
To create a baseline, use ``-bg``::

  icepack.create.case -ts base_suite -m wolf -testid v1 -bg version1 -bd $SCRATCH/ICEPACK_BASELINES -a P000000

will copy all the results from the test suite to ``$SCRATCH/ICEPACK_BASELINES/version1``.

To compare to a prior result, use ``-bc``::

  icepack.create.case -ts base_suite -m wolf -testid v2 -bc version1 -bd $SCRATCH/ICEPACK_BASELINES -a P000000

will compare all the results from this test suite to results saved before in $SCRATCH/ICEPACK_BASELINES/version1``.

To both create and compare, ``-bc`` and ``-bg`` can be combined::

  icepack.create.case -ts base_suite -m wolf -testid v2 -bg version2 -bc version1 -bd $SCRATCH/ICEPACK_BASELINES -a P000000

will save the current results to ``$SCRATCH/ICEPACK_BASELINES/version2`` and compare the current results to
results save before in ``$SCRATCH/ICEPACK_BASELINES/version1``.

In summary, 

- an individual test will have a case name like 
  ``[machine]_[test]_[grid]_[pes]_[soptions].[testid]``.
- A test suite will generate the individual tests under a directory called ``[suite_name].[testid]``.
- ``-bg`` will copy test results to the ``[bd_directory]/[bg_directory]/[test_name]``.
- ``-bc`` will compare results from  ``[bd_directory]/[bc_directory]/[test_name]``.

.. _comptesting:

Comparison testing
------------------

This feature is primarily used in test suites and has limited use in icepack, but is being
described for completeness.  If modifications to the column physics modules in
Icepack code generate differences (i.e. results are not bit-for-bit), then full 
comparisons tests will be necessary in CICE, comparing the modified column 
physics with the current version.

``-td`` provides a way to compare tests with each other.  The test is always compared relative to
the current case directory.  For instance::

  icepack.create.case -t smoke -m wolf -testid t01

creates a test case named wolf_smoke_col_1x1.t01::

  icepack.create.case -t smoke -m wolf -s run1year -testid t01 -td smoke_col_1x1

will create a test case named wolf_smoke_col_1x1_run1year.t01.  
An additional check will be done for the second test (because of the ``-td`` argument), and it will compare
the output from the first test "smoke_col_1x1" to the output from its test "smoke_col_1x1_run1year"
and generate a result for that.  It's important that the first test complete before the second test is done
and that the tests are created in parallel directories.
The ``-td`` option works only if the testid and the machine are the same for the baseline run and the 
current run, a basic feature associated with test suites.

Test Reporting
----------------------

The Icepack testing scripts have the capability of posting the test results
to an online dashboard, located `on CDash <http://my.cdash.org/index.php?project=myICEPACK>`_.

To post test suite results to CDash, add the ``-report`` option to **icepack.create.case**.
The base_suite will attempt to post the test results on CDash when the suite is complete.

If the results cannot be posted to CDash, the following information will be displayed::

 CTest submission failed.  To try the submission again run 
    ./run_ctest.csh -submit
 If you wish to submit the test results from another server, copy the 
 icepack_ctest.tgz file to another server and run 
    ./run_ctest.csh -submit

Examples
---------

To generate a baseline dataset for a test case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  ./icepack.create.case -t smoke -m wolf -bg icepackv6.0.0 -testid t00
  cd wolf_smoke_col_1x1.t00
  ./icepack.build
  ./icepack.submit

After job finishes, check output::

  cat test_output

To run a test case and compare to a baseline dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  ./icepack.create.case -t smoke -m wolf -bc icepackv6.0.0 -testid t01
  cd wolf_smoke_col_1x1.t01
  ./icepack.build
  ./icepack.submit

After job finishes, check output::

  cat test_output

To run a test suite to generate baseline data, review results, plot timeseries, and report results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  ./icepack.create.case -m wolf -ts base_suite -testid t02 -bg icepackv6.0.0bs -report

Once all jobs finish, concatenate all output and manually report results::

  cd base_suite.t02
  cat results.log

To plot a timeseries of "total ice extent", "total ice area", and "total ice volume"::

  ./timeseries.csh <directory>
  ls *.png

To run a test suite, compare to baseline data, generate a new baseline, and report the results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  ./icepack.create.case -m wolf -ts base_suite -testid t03 -bc icepackv6.0.0bs -bg icepackv6.0.0new -report

Case Settings
=====================

There are two important files that define the case, **icepack.settings** and 
**icepack_in**.  **icepack.settings** is a list of env variables that define many
values used to setup, build and run the case.  **icepack_in** is the input namelist file
for the icepack driver.  Variables in both files are described below.

.. _tabsettings:

Table of icepack settings
--------------------------

The **icepack.settings** file is reasonably well self documented.  Several of
the variables defined in the file are not used in Icepack.  They exist
to support the CICE model.

.. csv-table:: Table 7
   :header: "variable", "options/format", "description", "recommended value"
   :widths: 15, 15, 25, 20

   "ICE_MACHINE", " ", "machine name", "set by icepack.create.case"
   "ICE_CASENAME", " ", "case name", "set by icepack.create.case"
   "ICE_SANDBOX", " ", "sandbox directory", "set by icepack.create.case"
   "ICE_SCRIPTS", " ", "scripts directory", "set by icepack.create.case"
   "ICE_CASEDIR", " ", "case directory", "set by icepack.create.case"
   "ICE_RUNDIR", " ", "run directory", "set by icepack.create.case"
   "ICE_OBJDIR", " ", "compile directory", "${ICE_RUNDIR}/compile"
   "ICE_RSTDIR", " ", "unused", "${ICE_RUNDIR}/restart"
   "ICE_HSTDIR", " ", "unused", "${ICE_RUNDIR}/history"
   "ICE_LOGDIR", " ", "log directory", "${ICE_CASEDIR}/logs"
   "ICE_RSTPFILE", " ", "unused", "undefined"
   "ICE_DRVOPT", " ", "unused", "icepack"
   "ICE_CONSTOPT", " ", "unused", "cice"
   "ICE_IOTYPE", " ", "unused", "none"
   "ICE_CLEANBUILD", "true,false", "automatically clean before building", "true"
   "ICE_GRID", "col", "grid", "col"
   "ICE_NXGLOB", "4", "number of gridcells", "4"
   "ICE_NTASKS", "1", "number of tasks, must be set to 1", "1"
   "ICE_NTHRDS", "1", "number of threads per task, must be set to 1", "1"
   "ICE_TEST", " ", "test setting if using a test", "set by icepack.create.case"
   "ICE_TESTNAME", " ", "test name if using a test", "set by icepack.create.case"
   "ICE_BASELINE", " ", "baseline directory name, associated with icepack.create.case -bd", "set by icepack.create.case"
   "ICE_BASEGEN", " ", "baseline directory name for regression generation, associated with icepack.create.case -bg ", "set by icepack.create.case"
   "ICE_BASECOM", " ", "baseline directory name for regression comparison, associated with icepack.create.case -bc ", "set by icepack.create.case"
   "ICE_BFBCOMP", " ", "location of case for comparison, associated with icepack.create.case -td", "set by icepack.create.case"
   "ICE_SPVAL", " ", "unused", "UnDeFiNeD"
   "ICE_RUNLENGTH", " ", "batch run length default", "  00:10:00"
   "ICE_THREADED", "true,false", "force threading in compile, will always compile threaded if NTHRDS is gt 1", "false"
   "NICELYR", " ", "number of vertical layers in the ice", "7"
   "NSNWLYR", " ", "number of vertical layers in the snow", "1"
   "NICECAT", " ", "number of ice thickness categories", "5"
   "TRAGE", "0,1", "ice age tracer", "1"
   "TRFY", "0,1", "first-year ice area tracer", "1"
   "TRLVL", "0,1", "deformed ice tracer", "1"
   "TRPND", "0,1", "melt pond tracer", "1"
   "NTRAERO", " ", "number of aerosol tracers", "1"
   "TRBRI", "0,1", "brine height tracer", "0"
   "TRZS", "0,1", "zsalinity tracer, needs TRBRI=1", "0"
   "TRBGCS", "0,1", "skeletal layer tracer, needs TRBGCZ=0", "0"
   "TRBGCZ", "0,1", "zbgc tracers, needs TRBGCS=0 and TRBRI=1", "0"
   "NBGCLYR", " ", "number of zbgc layers", "7"
   "TRZAERO", "0-6", "number of z aerosol tracers", "0"
   "TRALG", "0,1,2,3", "number of algal tracers", "0"
   "TRDOC", "0,1,2,3", "number of dissolved organic carbon", "0"
   "TRDIC", "0,1", "number of dissolved inorganic carbon", "0"
   "TRDON", "0,1", "number of dissolved organic nitrogen", "0"
   "TRFEP", "0,1,2", "number of particulate iron tracers", "0"
   "TRFED", "0,1,2", "number of dissolved iron tracers", "0"
   "CAM_ICE", " ", "unused", "no"
   "DITTO", " ", "unused", "no"
   "BARRIERS", " ", "unused", "no"
   "ICE_BLDDEBUG", "true,false", "turn on compile debug flags", "false"
   "NUMIN", " ", "unused", "11"
   "NUMAX", " ", "unused", "99"


.. _tabnamelist:

Table of namelist inputs
--------------------------

CHECK

Namelist is part of the icepack driver code and is used to setup testing of the
column physics.

.. _tab-namelist:

.. csv-table:: Table 8
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
   "``ice_ic``", "``default``", "latitude and sst dependent", "default"
   "", "``none``", "no ice", ""
   "", "path/file", "restart file name", ""
   "``restart_dir``", "path/", "path to restart directory", ""
   "``dumpfreq``", "``y``", "write restart every ``dumpfreq_n`` years", "y"
   "", "``m``", "write restart every ``dumpfreq_n`` months", ""
   "", "``d``", "write restart every ``dumpfreq_n`` days", ""
   "", "", "*Model Output*", ""
   "``diagfreq``", "integer", "frequency of diagnostic output in ``dt``", "24"
   "", "*e.g.*, 10", "once every 10 time steps", ""
   "``diag_file``", "filename", "diagnostic output file (script may reset)", ""
   "", "", "", ""
   "*grid_nml*", "", "", ""
   "", "", "*Grid*", ""
   "``kcatbound``", "``0``", "original category boundary formula", "0"
   "", "``1``", "new formula with round numbers", ""
   "", "``2``", "WMO standard categories", ""
   "", "``-1``", "one category", ""
   "", "", "", ""
   "*tracer_nml*", "", "", ""
   "", "", "*Tracers*", ""
   "``tr_iage``", "true/false", "ice age", ""
   "``tr_FY``", "true/false", "first-year ice area", ""
   "``tr_lvl``", "true/false", "level ice area and volume", ""
   "``tr_pond_cesm``", "true/false", "CESM melt ponds", ""
   "``tr_pond_topo``", "true/false", "topo melt ponds", ""
   "``tr_pond_lvl``", "true/false", "level-ice melt ponds", ""
   "``tr_aero``", "true/false", "aerosols", ""
   "", "", "", ""
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
   "``shortwave``", "``ccsm3``", "NCAR CCSM3 distribution method", "'dEdd'"
   "", "``dEdd``", "Delta-Eddington method", ""
   "``albedo_type``", "``ccsm3``", "NCAR CCSM3 albedos", "‘ccsm3’"
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
   "*forcing_nml*", "", "", ""
   "", "", "*Forcing*", ""
   "``formdrag``", "true/false", "calculate form drag", ""
   "``atmbndy``", "``default``", "stability-based boundary layer", "‘default’"
   "", "``constant``", "bulk transfer coefficients", ""
   "``fyear_init``", "yyyy", "first year of atmospheric forcing data", ""
   "``ycycle``", "integer", "number of years in forcing data cycle", ""
   "``atm_data_type``", "``default``", "constant values defined in the code", ""
   "", "``clim``", "monthly climatology", ""
   "", "``CFS``", "CFS model output", ""
   "", "``ISPOL``", "ISPOL experiment data", ""
   "``data_dir``", "path/", "path to forcing data directory", ""
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
   "", "``mushy``", "matches mushy-layer thermo (ktherm=2)", ""
   "``ustar_min``", "real", "minimum value of ocean friction velocity", "0.0005 m/s"
   "``fbot_xfer_type``", "``constant``", "constant ocean heat transfer coefficient", ""
   "", "``Cdn_ocn``", "variable ocean heat transfer coefficient", ""
   "``update_ocn_f``", "true", "include frazil water/salt fluxes in ocn fluxes", ""
   "", "false", "do not include (when coupling with POP)", ""
   "``l_mpond_fresh``", "true", "retain (topo) pond water until ponds drain", ""
   "", "false", "release (topo) pond water immediately to ocean", ""
   "``oceanmixed_ice``", "true/false", "active ocean mixed layer calculation", "``.true.`` (if uncoupled)"
   "``ocn_data_type``", "``default``", "constant values defined in the code", ""
   "", "``ISPOL``", "ISPOL experiment data", ""
   "``bgc_data_type``", "``default``", "constant values defined in the code", ""
   "", "``ISPOL``", "ISPOL experiment data", ""
   "``oceanmixed_file``", "filename", "data file containing ocean forcing data", ""
   "``restore_ocn``", "true/false", "restore sst to data", ""
   "``trestore``", "integer", "sst restoring time scale (days)", ""
   "", "", "", ""
   "*zbgc_nml*", "", "", ""
   "", "", "*Biogeochemistry*", ""
   "``tr_brine``", "true/false", "brine height tracer (needs TRBRI 1 in comp_ice)", ".true."
   "``restart_hbrine``", "true/false", "restart the brine height tracer (automatically turned on if restart = .true.)", ""
   "``tr_zaero``", "true/false", "turns on black carbon and dust aerosols", ""
   "``modal_aero``", "true/false", "turns on a modal aerosol option", ""
   "``skl_bgc``", "true/false", "turns on a single bottom layer biogeochemistry. z_tracers and solve_zbgc must be false", ""
   "``z_tracers``", "true/false", "turns on a vertically resolved transport", ""
   "``dEdd_algae``", "true/false", "Include radiative impact of algae and aerosols in the delta-Eddington shortwave scheme. Requires shortwave = 'dEdd'.", ""
   "``solve_zbgc``", "true/false", "turns on the biochemistry using z_tracers (specify algal numbers in comp_ice TRALG)", ""
   "``bgc_flux_type``", "``Jin2006``", "ice–ocean flux type for bottom layer tracers only :cite:`JDWSTWLG06`", ""
   "``restore_bgc``", "true/false", "restores upper ocean concentration fields to data values for nitrate and silicate", ""
   "``restart_bgc``", "true/false", "restarts biogeochemical tracers (automatically turned on if restart = .true.)", ""
   "``scale_bgc``", "true/false", "Initialize biogeochemical profiles to scale with prognosed salinity profile", ""
   "``solve_zsal``", "true/false", "prognostic salinity tracer used with ktherm = 1", ""
   "``restart_zsal``", "true/false", "restarts zsalinity", ""
   "``bgc_data_dir``", "``/nitrate_and_silicate/forcing_directory/``", "", ""
   "``sil_data_type``", "``default``", "fixed, spatially homogeneous value for silicate. 'clim' data file (see ice_forcing_bgc.F90) :cite:`GLBA06`", ""
   "``nit_data_type``", "``default``", "fixed, spatially homogeneous value for nitrate. 'clim' data file (see ice_forcing_bgc.F90) :cite:`GLBA06`", ""
   "``fe_data_type``", "``default``", "fixed, spatially homogeneous value for iron", ""
   "``tr_bgc_Nit``", "true/false", "nitrate tracer", ""
   "``tr_bgc_C``", "true/false", "dissolved organic carbon tracers and dissolved inorganic carbon tracers (not yet implemented)", ""
   "``tr_bgc_chl``", "true/false", "dummy variable for now. Chl is simply fixed ratio of algal Nitrogen", ""
   "``tr_bgc_Am``", "true/false", "Ammonium", ""
   "``tr_bgc_Sil``", "true/false", "Silicate", ""
   "``tr_bgc_DMS``", "true/false", "Three tracers: DMS dimethyl sulfide, DMSPp (particulate, assumed to be a fixed ratio of sulfur to algal nitrogen) and DMSPd (dissolved)", ""
   "``tr_bgc_PON``", "true/false", "passive purely mobile ice tracer with ocean concentration equivalent to nitrate", ""
   "``tr_bgc_hum``", "true/false", "refractory DOC or DON (units depend on the ocean source)", ""
   "``tr_bgc_DON``", "true/false", "dissolved organic nitrogen", ""
   "``tr_bgc_Fe``", "true/false", "dissolved iron and particulate iron", ""
   "``grid_o``", "real", "ice-ocean surface layer thickness (bgc transport scheme)", "0.006"
   "``grid_o_t``", "real", "ice-atmosphere surface layer thickness (bgc transport scheme)", "0.006"
   "``l_sk``", "real", "length scale in gravity drainage parameterization (bgc transport scheme)", "0.024"
   "``grid_oS``", "real", "ice-ocean surface layer thickness (zsalinity transport scheme)", "0.0"
   "``l_skS``", "real", "ice-atmosphere surface layer thickness (zsalinity transport scheme)", "0.028"
   "``phi_snow``", "real", "snow porosity at the ice-snow interface. if :math:`<0` then phi_snow is computed from snow density", "-0.3"
   "``initbio_frac``", "real", "for each bgc tracer, specifies the fraction of the ocean concentration that is retained in the ice during initial new ice formation", "0.8"
   "``frazil_scav``", "real", "for each bgc tracer, specifies the fraction or multiple of the ocean concentration that is retained in the ice from frazil ice formation", "0.8"
   "``ratio_si2N_diatoms``", "real", "algal Si to N (:math:`mol/mol`) for diatoms", "1.8"
   "``ratio_si2N_sp``", "real", "algal Si to N (:math:`mol/mol`) for small phytoplankton", "0.0"
   "``ratio_si2N_phaeo``", "real", "algal Si to N (:math:`mol/mol`) for phaeocystis", "0.0"
   "``ratio_S2N_diatoms``", "real", "algal S to N (:math:`mol/mol`) for diatoms", "0.03"
   "``ratio_S2N_sp``", "real", "algal S to N (:math:`mol/mol`) for small phytoplankton", "0.03"
   "``ratio_S2N_phaeo``", "real", "algal S to N (:math:`mol/mol`) for phaeocystis", "0.03"
   "``ratio_Fe2C_diatoms``", "real", "algal Fe to C (:math:`\mu mol/mol`) for diatoms", "0.0033"
   "``ratio_Fe2C_sp``", "real", "algal Fe to C (:math:`\mu mol/mol`) for small phytoplankton", "0.0033"
   "``ratio_Fe2C_phaeo``", "real", "algal Fe to C (:math:`\mu mol/mol`) for phaeocystis", "0.1"
   "``ratio_Fe2N_diatoms``", "real", "algal Fe to N (:math:`\mu mol/mol`) for diatoms", "0.023"
   "``ratio_Fe2N_sp``", "real", "algal Fe to N (:math:`\mu mol/mol`) for small phytoplankton", "0.023"
   "``ratio_Fe2N_phaeo``", "real", "algal Fe to N (:math:`\mu mol/mol`) for phaeocystis", "0.7"
   "``ratio_Fe2DON``", "real", "Fe to N of DON (:math:`nmol/mol`)", "0.023"
   "``ratio_Fe2DOC_s``", "real", "Fe to C of DOC for saccharids (:math:`nmol/mol`)", "0.1"
   "``ratio_Fe2DOC_l``", "real", "Fe to C of DOC for lipids (:math:`nmol/mol`)", "0.033"
   "``fr_resp``", "real", "fraction of algal growth lost due to respiration", "0.05"
   "``tau_min``", "real", "rapid mobile to stationary exchanges (:math:`s`)", "5200.0"
   "``tau_max``", "real", "long time mobile to stationary exchanges (:math:`s`)", "1.73e5"
   "``algal_vel``", "real", "0.5 :math:`cm/day (m/s)`", "1.11e-8"
   "``R_dFe2dust``", "real", "g/g (3.5% content)", "0.035"
   "``dustFe_sol``", "real", "solubility fraction", "0.005" 
   "``chlabs_diatoms``", "real", "diatoms chl absorption (:math:`1/m/(mg/m^3)`)", "0.03"
   "``chlabs_sp``", "real", "small phytoplankton chl absorption (:math:`1/m/(mg/m^3)`)", "0.01"
   "``chlabs_phaeo``", "real", "phaeocystis chl absorption (:math:`1/m/(mg/m^3)`)", "0.05"
   "``alpha2max_low_diatoms``", "real", "diatoms light limitation (:math:`(W/m^2)^{-1}`)", "0.8"    
   "``alpha2max_low_sp``", "real", "small phytoplankton light limitation (:math:`(W/m^2)^{-1}`)", "0.67"    
   "``alpha2max_low_phaeo``", "real", "phaeocystis light limitation (:math:`(W/m^2)^{-1}`)", "0.67"    
   "``beta2max_diatoms``", "real", "diatoms light inhibition (:math:`(W/m^2)^{-1}`)", "0.018"    
   "``beta2max_sp``", "real", "small phytoplankton light inhibition (:math:`(W/m^2)^{-1}`)", "0.0025"    
   "``beta2max_phaeo``", "real", "phaeocystis light inhibition (:math:`(W/m^2)^{-1}`)", "0.01" 
   "``mu_max_diatoms``", "real", "diatoms maximum growth rate (:math:`day^{-1}`)", "1.2" 
   "``mu_max_sp``", "real", "small phytoplankton maximum growth rate (:math:`day^{-1}`)", "0.851"  
   "``mu_max_phaeo``", "real", "phaeocystis maximum growth rate (:math:`day^{-1}`)", "0.851" 
   "``grow_Tdep_diatoms``", "real", "diatoms Temperature dependence of growth (:math:`^o`\ C\ :math:`^{-1})`", "0.06" 
   "``grow_Tdep_sp``", "real", "small phytoplankton Temperature dependence of growth :math:`^o`\ C\ :math:`^{-1}`", "0.06"  
   "``grow_Tdep_phaeo``", "real", "phaeocystis Temperature dependence of growth :math:`^o`\ C\ :math:`^{-1}`", "0.06"  
   "``fr_graze_diatoms``", "real", "diatoms fraction grazed", "0.01" 
   "``fr_graze_sp``", "real", "small phytoplankton fraction grazed", "0.1"  
   "``fr_graze_phaeo``", "real", "phaeocystis fraction grazed", "0.1" 
   "``mort_pre_diatoms``", "real", "diatoms mortality (:math:`day^{-1}`)", "0.007" 
   "``mort_pre_sp``", "real", "small phytoplankton mortality (:math:`day^{-1}`)", "0.007"  
   "``mort_pre_phaeo``", "real", "phaeocystis mortality (:math:`day^{-1}`)", "0.007" 
   "``mort_Tdep_diatoms``", "real", "diatoms temperature dependence of mortality :math:`^o`\ C\ :math:`^{-1}`", "0.03" 
   "``mort_Tdep_sp``", "real", "small phytoplankton temperature dependence of mortality (:math:`^o`\ C\ :math:`^{-1}`)", "0.03"  
   "``mort_Tdep_phaeo``", "real", "phaeocystis temperature dependence of mortality (:math:`^o`\ C\ :math:`^{-1}`)", "0.03" 
   "``k_exude_diatoms``", "real", "diatoms algal exudation (:math:`day^{-1}`)", "0.0" 
   "``k_exude_sp``", "real", "small phytoplankton algal exudation (:math:`day^{-1}`)", "0.0"  
   "``k_exude_phaeo``", "real", "phaeocystis algal exudation (:math:`day^{-1}`)", "0.0"           
   "``K_Nit_diatoms``", "real", "datoms nitrate half saturation (:math:`mmol/m^3`)", "1.0" 
   "``K_Nit_sp``", "real", "small phytoplankton nitrate half saturation (:math:`mmol/m^3`)", "1.0"  
   "``K_Nit_phaeo``", "real", "phaeocystis nitrate half saturation (:math:`mmol/m^3`)", "1.0"           
   "``K_Am_diatoms``", "real", "diatoms ammonium half saturation (:math:`mmol/m^3`)", "0.3" 
   "``K_Am_sp``", "real", "small phytoplankton ammonium half saturation (:math:`mmol/m^3`)", "0.3"  
   "``K_Am_phaeo``", "real", "phaeocystis ammonium half saturation (:math:`mmol/m^3`)", "0.3"   
   "``K_Sil_diatoms``", "real", "diatoms silicate half saturation (:math:`mmol/m^3`)", "4.0" 
   "``K_Sil_sp``", "real", "small phytoplankton silicate half saturation (:math:`mmol/m^3`)", "0.0"  
   "``K_Sil_phaeo``", "real", "phaeocystis silicate half saturation (:math:`mmol/m^3`)", "0.0" 
   "``K_Fe_diatoms``", "real", "diatoms iron half saturation (:math:`nM`)", "1.0" 
   "``K_Fe_sp``", "real", "small phytoplankton iron half saturation (:math:`nM`)", "0.2"  
   "``K_Fe_phaeo``", "real", "phaeocystis iron half saturation (:math:`nM`)", "0.1"   
   "``f_don_protein``", "real", "fraction of spilled grazing to proteins", "0.6" 
   "``kn_bac_protein``", "real", "Bacterial degredation of DON (:math:`day^{-1}`)", "0.03"                
   "``f_don_Am_protein``", "real", "fraction of remineralized DON to ammonium", "0.25"
   "``f_doc_s``", "real", "fraction of mortality to DOC saccharids", "0.4"
   "``f_doc_l``", "real", "fraction of mortality to DOC lipids", "0.4"  
   "``f_exude_s``", "real", "fraction of exudation to DOC saccharids", "1.0"
   "``f_exude_l``", "real", "fraction of exudation to DOC lipids", "1.0"  
   "``k_bac_s``", "real", "bacterial degredation of DOC (:math:`day^{-1}`) saccharids", "0.03"
   "``k_bac_l``", "real", "bacterial degredation of DOC (:math:`day^{-1}`) lipids", "0.03"  
   "``T_max``", "real", "maximum temperature (:math:`^o`\ C)", "0.0"
   "``fsal``", "real", "Salinity limitation (ppt)", "1.0"
   "``op_dep_min``", "real", "Light attenuates for optical depths exceeding min", "0.1"
   "``fr_graze_s``", "real", "fraction of grazing spilled or slopped", "0.5"
   "``fr_graze_e``", "real", "fraction of assimilation excreted", "0.5"
   "``fr_mort2min``", "real", "fractionation of mortality to Am", "0.5"
   "``fr_dFe``", "real", "fraction of remineralized nitrogen (algal iron)", "0.3"
   "``k_nitrif``", "real", "nitrification rate (:math:`day^{-1}`)", "0.0"
   "``t_iron_conv``", "real", "desorption loss pFe to dFe (day)", "3065.0"
   "``max_loss``", "real", "restrict uptake to % of remaining value", "0.9"
   "``max_dfe_doc1``", "real", "max ratio of dFe to saccharides in the ice (:math:`nM Fe/\mu M C`)", "0.2"
   "``fr_resp_s``", "real", "DMSPd fraction of respiration loss as DMSPd", "0.75"
   "``y_sk_DMS``", "real", "fraction conversion given high yield", "0.5"
   "``t_sk_conv``", "real", "Stefels conversion time (:math:`day`)", "3.0"
   "``t_sk_ox``", "real", "DMS oxidation time (:math:`day`)", "10.0"
   "``algaltype_diatoms``", "real", "mobility type between stationary <--> mobile for diatoms", "0.0"
   "``algaltype_sp``", "real", "mobility type between stationary <--> mobile for small phytoplankton", "0.5"
   "``algaltype_phaeo``", "real", "mobility type between stationary <--> mobile for phaeocystis", "0.5"
   "``nitratetype``", "real", "mobility type between stationary <--> mobile for nitrate", "-1.0"
   "``ammoniumtype``", "real", "mobility type between stationary <--> mobile for ammonium", "1.0"
   "``silicatetype``", "real", "mobility type between stationary <--> mobile for silicate", "-1.0"
   "``dmspptype``", "real", "mobility type between stationary <--> mobile for DMSP particulate", "0.5"
   "``dmspdtype``", "real", "mobility type between stationary <--> mobile for DMSP dissolved", "-1.0"
   "``humtype``", "real", "mobility type between stationary <--> mobile for humic matter", "1.0"
   "``doctype_s``", "real", "mobility type between stationary <--> mobile for DOC saccharids", "0.5"
   "``doctype_l``", "real", "mobility type between stationary <--> mobile for DOC lipids", "0.5"
   "``dontype_protein``", "real", "mobility type between stationary <--> mobile for proteins", "0.5"
   "``fedtype_1``", "real", "mobility type between stationary <--> mobile for FeD", "0.5"
   "``feptype_1``", "real", "mobility type between stationary <--> mobile for FeP", "0.5"
   "``zaerotype_bc1``", "real","mobility type between stationary <--> mobile for zaerotype_bc1",  "1.0"
   "``zaerotype_bc2``", "real", "mobility type between stationary <--> mobile for zaerotype_bc2", "1.0"
   "``zaerotype_dust1``", "real", "mobility type between stationary <--> mobile for dust1", "1.0"
   "``zaerotype_dust2``", "real", "mobility type between stationary <--> mobile for dust2", "1.0"
   "``zaerotype_dust3``", "real", "mobility type between stationary <--> mobile for dust3", "1.0"
   "``zaerotype_dust4``", "real", "mobility type between stationary <--> mobile for dust4", "1.0"
   "``ratio_C2N_diatoms``", "real", "diatom algal C to N ratio (:math:`mol/mol`)", "7.0"
   "``ratio_C2N_sp``", "real", "small phytoplankton algal C to N ratio (:math:`mol/mol`)", "7.0"
   "``ratio_C2N_phaeo``", "real", "phaeocystis algal C to N ratio (:math:`mol/mol`)", "7.0"
   "``ratio_chl2N_diatoms``", "real", "diatom algal chlorophyll to N ratio (:math:`mg/mmol`)", "2.1"
   "``ratio_chl2N_sp``", "real", "small phytoplankton algal chlorophyll to N ratio (:math:`mg/mmol`)", "1.1"
   "``ratio_chl2N_phaeo``", "real", "phaeocystis algal chlorophyll to N ratio (:math:`mg/mmol`)", "0.84"
   "``F_abs_chl_diatoms``", "real", "diatom scales absorbed radiation for dEdd", "2.0"
   "``F_abs_chl_sp``", "real", "small phytoplankton scales absorbed radiation for dEdd", "4.0"
   "``F_abs_chl_phaeo``", "real", "phaeocystis scales absorbed radiation for dEdd", "5.0"       
   "``ratio_C2N_proteins``", "real", "ratio of C to N in proteins (:math:`mol/mol`)", "7.0"
   "", "", "", ""

.. commented out below
..   "``dbug``", "true/false", "if true, write extra diagnostics", "``.false.``"
..   "``atm_data_format``", "``nc``", "read  atmo forcing files", ""
..   "", "``bin``", "read direct access, binary files", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "", "``NICE``", "N-ICE experiment data", ""

  
.. _tuning:

BGC Tuning Parameters
~~~~~~~~~~~~~~~~~~~~~

Biogeochemical tuning parameters are specified as namelist options in
**icepack\_in**. Table :ref:`tab-bio-tracers2` provides a list of parameters
used in the reaction equations, their representation in the code, a
short description of each and the default values. Please keep in mind
that there has only been minimal tuning of the model.

:ref:`tab-bio-tracers2` :*Biogeochemical Reaction Parameters*

.. _tab-bio-tracers2:

.. csv-table:: Table 9
   :header: "Text Variable", "Variable in code", "Description", "Value", "units"
   :widths: 7, 20, 15, 15, 15

   ":math:`f_{graze}`", "fr\_graze(1:3)", "fraction of growth grazed", "0, 0.1, 0.1", "1"
   ":math:`f_{res}`", "fr\_resp", "fraction of growth respired", "0.05", "1"
   ":math:`l_{max}`", "max\_loss", "maximum tracer loss fraction", "0.9", "1"
   ":math:`m_{pre}`", "mort\_pre(1:3)", "maximum mortality rate", "0.007, 0.007, 0.007", "day\ :math:`^{-1}`"
   ":math:`m_{T}`", "mort\_Tdep(1:3)", "mortality temperature decay", "0.03, 0.03, 0.03", ":math:`^o`\ C\ :math:`^{-1}`"
   ":math:`T_{max}`", "T\_max", "maximum brine temperature", "0", ":math:`^o`\ C"
   ":math:`k_{nitr}`", "k\_nitrif", "nitrification rate", "0", "day\ :math:`^{-1}`"
   ":math:`f_{ng}`", "fr\_graze\_e", "fraction of grazing excreted", "0.5", "1"
   ":math:`f_{gs}`", "fr\_graze\_s", "fraction of grazing spilled", "0.5", "1"
   ":math:`f_{nm}`", "fr\_mort2min", "fraction of mortality to :math:`{\mbox{NH$_4$}}`", "0.5", "1"
   ":math:`f_{dg}`", "f\_don", "frac. spilled grazing to :math:`{\mbox{DON}}`", "0.6", "1"
   ":math:`k_{nb}`", "kn\_bac :math:`^a`", "bacterial degradation of :math:`{\mbox{DON}}`", "0.03", "day\ :math:`^{-1}`"
   ":math:`f_{cg}`", "f\_doc(1:3)", "fraction of mortality to :math:`{\mbox{DOC}}`", "0.4, 0.4, 0.2 ", "1"
   ":math:`R_{c:n}^c`", "R\_C2N(1:3)", "algal carbon to nitrogen ratio", "7.0, 7.0, 7.0", "mol/mol"
   ":math:`k_{cb}`", "k\_bac1:3\ :math:`^a`", "bacterial degradation of DOC", "0.03, 0.03, 0.03", "day\ :math:`^{-1}`"
   ":math:`\tau_{fe}`", "t\_iron\_conv", "conversion time pFe :math:`\leftrightarrow` dFe", "3065.0 ", "day"
   ":math:`r^{max}_{fed:doc}`", "max\_dfe\_doc1", "max ratio of dFe to saccharids", "0.1852", "nM Fe\ :math:`/\mu`\ M C"
   ":math:`f_{fa}`", "fr\_dFe  ", "fraction of remin. N to dFe", "0.3", "1"
   ":math:`R_{fe:n}`", "R\_Fe2N(1:3)", "algal Fe to N ratio", "0.023, 0.023, 0.7", "mmol/mol"
   ":math:`R_{s:n}`", "R\_S2N(1:3)", "algal S to N ratio", "0.03, 0.03, 0.03", "mol/mol"
   ":math:`f_{sr}`", "fr\_resp\_s", "resp. loss as DMSPd", "0.75", "1"
   ":math:`\tau_{dmsp}`", "t\_sk\_conv", "Stefels rate", "3.0", "day"
   ":math:`\tau_{dms}`", "t\_sk\_ox", "DMS oxidation rate", "10.0", "day"
   ":math:`y_{dms}`", "y\_sk\_DMS", "yield for DMS conversion", "0.5", "1"
   ":math:`K_{{\mbox{NO$_3$}}}`", "K\_Nit(1:3)", ":math:`{\mbox{NO$_3$}}` half saturation constant", "1,1,1", "mmol/m\ :math:`^{3}`"
   ":math:`K_{{\mbox{NH$_4$}}}`", "K\_Am(1:3)", ":math:`{\mbox{NH$_4$}}` half saturation constant", "0.3, 0.3, 0.3", "mmol/m\ :math:`^{-3}`"
   ":math:`K_{{\mbox{SiO$_3$}}}`", "K\_Sil(1:3)", "silicate half saturation constant", "4.0, 0, 0", "mmol/m\ :math:`^{-3}`"
   ":math:`K_{{\mbox{fed}}}`", "K\_Fe(1:3)", "iron half saturation constant", "1.0, 0.2, 0.1", ":math:`\mu`\ mol/m\ :math:`^{-3}`"
   ":math:`op_{min}`", "op\_dep\_min", "boundary for light attenuation", "0.1", "1"
   ":math:`chlabs`", "chlabs(1:3)", "light absorption length per chla conc.", "0.03, 0.01, 0.05", "1\ :math:`/`\ m\ :math:`/`\ (mg\ :math:`/`\ m\ :math:`^{3}`)"
   ":math:`\alpha`", "alpha2max\_low(1:3)", "light limitation factor", "0.25, 0.25, 0.25", "m\ :math:`^2`/W"
   ":math:`\beta`", "beta2max(1:3)", "light inhibition factor", "0.018, 0.0025, 0.01", "m\ :math:`^2`/W"
   ":math:`\mu_{max}`", "mu\_max(1:3)", "maximum algal growth rate", "1.44, 0.851, 0.851", "day\ :math:`^{-1}`"
   ":math:`\mu_T`", "grow\_Tdep(1:3)", "temperature growth factor", "0.06, 0.06, 0.06", "day\ :math:`^{-1}`"
   ":math:`f_{sal}`", "fsal", "salinity growth factor", "1", "1"
   ":math:`R_{si:n}`", "R\_Si2N(1:3)", "algal silicate to nitrogen", "1.8, 0, 0", "mol/mol"

:math:`^a` only (1:2) of DOC and DOC parameters have physical meaning


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

Underflows
-----------

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

Known bugs and other issues
---------------------------

-   With the old CCSM radiative scheme (``shortwave`` = ‘default’ or
    ‘ccsm3’), a sizable fraction (more than 10%) of the total shortwave
    radiation is absorbed at the surface but should be penetrating into
    the ice interior instead. This is due to use of the aggregated,
    effective albedo rather than the bare ice albedo 
    when ``snowpatch`` < 1.

-   The linear remapping algorithm for thickness is not monotonic for tracers.

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
