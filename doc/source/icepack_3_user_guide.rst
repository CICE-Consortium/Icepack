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

Initialization and Forcing
---------------------------

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


Running Icepack
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

.. _scripts:

Scripts
-------------

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
-------------

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

Input Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The input data space is defined on a per machine basis by the ``ICE_MACHINE_INPUTDATA`` 
variable in the **env.[machine]** file.  That file space is often shared among multiple 
users, and it can be desirable to consider using a common file space with group read 
and write permissions such that a set of users can update the inputdata area as 
new datasets are available.

The latest input data will be available thru the CICE Consortium github wiki.

Machine Account Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The machine account default is specified by the variable ``ICE_MACHINE_ACCT`` in 
the **env.[machine]** file.  The easiest way to change a user's default is to 
create a file in your home directory called **.cice\_proj** and add your 
preferred account name to the first line.  
There is also an option (``-a``) in **icepack.create.case** to define the account number.  
The order of precedent is **icepack.create.case** command line option, 
**.cice\_proj** setting, and then value in the **env.[machine]** file.


Run Directories
-----------------

The **icepack.create.case** script creates a case directory.  However, the model 
is actually built and run under the ``ICE_OBJDIR`` and ``ICE_RUNDIR`` directories
as defined in the **icepack.settings** file.

Build and run logs will be copied from the run directory into the case **logs/** 
directory when complete.


Local modifications
--------------------------

Scripts and other case settings can be changed manually in the case directory and
used.  Source code can be modified in the main sandbox.  When changes are made, the code
should be rebuilt before being resubmitted.  It is always recommended that users
modify the scripts and input settings in the case directory, NOT the run directory.
In general, files in the run directory are overwritten by versions in the case
directory when the model is built, submitted, and run.


Forcing data
------------

CHECK once we've settled on a forcing suite:

The code is currently configured to run in standalone mode on a 4-cell grid using 
atmospheric data, available as detailed on the `wiki <https://github.com/CICE-Consortium/Icepack/wiki/Testing-Icepack>`_.
These data files are designed only for testing the code, not for use in production 
runs or as observational data.  Please do not publish results based on these data
sets.  Module **configuration/driver/icedrv\_forcing.F90**
can be modified to change the forcing data. 



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
--------------------------------

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
--------------------------------

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

.. csv-table:: Table 7
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
----------------------------------

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
----------------------------------

.. THIS IS NOT REALLY USED IN ICEPACK, should we leave this out of the icepack documentation

This feature is primarily used in test suites and has limited use in icepack, but is being
described for completeness.

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

Icepack Test Reporting
------------------------------------

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
--------------------------

To generate a baseline dataset for a test case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  ./icepack.create.case -t smoke -m wolf -bg icepackv6.0.0 -testid t00
  cd wolf_smoke_col_1x1.t00
  ./icepack.build
  ./icepack.submit

After job finishes, check output::

  cat test_output

To run a test case and compare to a baseline dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. csv-table:: Table 9
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
   "*zbgc_nml*", "", "", ""
   "", "", "*Biogeochemistry*", ""
   "``tr_brine``", "true/false", "brine height tracer", ""
   "``skl_bgc``", "true/false", "biogeochemistry", ""
   "``bgc_flux_type``", "``Jin2006``", "ice–ocean flux velocity of :cite:`JDWSTWLG06`", ""
   "", "``constant``", "constant ice–ocean flux velocity", ""
   "``restore_bgc``", "true/false", "restore nitrate/silicate to data", ""
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


.. commented out below
..   "``dbug``", "true/false", "if true, write extra diagnostics", "``.false.``"
..   "``atm_data_format``", "``nc``", "read  atmo forcing files", ""
..   "", "``bin``", "read direct access, binary files", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "", "``NICE``", "N-ICE experiment data", ""

Adding things
====================

We require that any changes made to the code be implemented in such a way that they can
be "turned off" through namelist flags.  In most cases, code run with such changes should 
be bit-for-bit identical with the unmodified code.  Occasionally, non-bit-for-bit changes
are necessary, e.g. associated with an unavoidable change in the order of operations. In
these cases, changes should be made in stages to isolate the non-bit-for-bit changes, 
so that those that should be bit-for-bit can be tested separately.

.. _addtrcr:

Tracers
--------------

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

#. **icedrv\_domain\_size.F90**: increase ``max_ntrcr`` (can also add option
   to **icepack.settings** and **icepack.build**)

#. **icedrv\_state.F90**: declare ``nt_[tracer]`` and ``tr_[tracer]``

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

#. **icedrv\_InitMod.F90**: initialize tracer (includes reading restart
   file)

#. **icedrv\_RunMod.F90**, **icedrv\_step\_mod.F90**:

   -  call routine to write tracer restart data

   -  call physics routines in **icepack\_[tracer].F90** (often called from
      **icedrv\_step\_mod.F90**)

#. **icedrv\_restart.F90**: define restart variables

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
