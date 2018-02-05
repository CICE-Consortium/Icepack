:tocdepth: 3

.. _running_icepack:

Running Icepack
====================

Quick-start instructions are provided in the :ref:`quickstart` section.

.. _scripts:

Scripts
-------

The icepack scripts are written to allow quick setup of cases and tests.  Once a case is 
generated, users can manually modify the namelist and other files to custom configure
the case.  Several settings are available via scripts as well.

Overview
~~~~~~~~

Most of the scripts that configure, build and run Icepack are contained in 
the directory **configuration/scripts/**, except for **icepack.setup**, which is
in the main directory.  **icepack.setup** is the main script that generates a case. 

Users may need to port the scripts to their local machine.
Specific instructions for porting are provided in :ref:`porting`.

``icepack.setup -h`` will provide the latest information about how to use the tool.
``icepack.setup --help`` will provide an extended version of the help.
There are three usage modes,

* ``--case`` or ``-c`` creates individual stand alone cases.
* ``--test`` creates individual tests.  Tests are just cases that have some extra automation in order to carry out particular tests such as exact restart.
* ``--suite`` creates a test suite.  Test suites are predefined sets of tests and ``--suite`` provides the ability to quick setup, build, and run a full suite of tests.

All modes will require use of ``--mach`` or ``-m`` to specify the machine and case and test modes 
can use ``--set`` or ``-s`` to define specific options.  ``--test`` and ``--suite`` will require ``--testid`` to be set 
and both of the test modes can use ``--bdir``, ``--bgen``, ``--bcmp``, and ``--diff`` to generate (save) results and compare results with prior results.
Testing will be described in greater detail in the :ref:`testing` section.

Again, ``icepack.setup --help`` will show the latest usage information including 
the available ``--set`` options, the current ported machines, and the test choices.

To create a case, run **icepack.setup**::

  icepack.setup -c mycase -m machine
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
- To change batch settings, manually edit the top of the **icepack.run** or **icepack.test** (if running a test) file
- To turn on the debug compiler flags, set ``ICE_BLDDEBUG`` in **icepack.setttings** to true
- To change compiler options, manually edit the Macros file
- To clean the build before each compile, set ``ICE_CLEANBUILD`` in **icepack.settings** to true.  To not clean before the build, set ``ICE_CLEANBUILD`` in **icepack.settings** to false

To build and run::

  ./icepack.build
  ./icepack.submit

The build and run log files will be copied into the logs directory in the case directory.
Other model output will be in the run directory.  The run directory is set in **icepack.settings**
via the **ICE_RUNDIR** variable.  To modify the case setup, changes should be made in the
case directory, NOT the run directory.

.. _case_options:

Command Line Options
~~~~~~~~~~~~~~~~~~~~

``icepack.setup -h`` provides a summary of the command line options.  There are three different modes, ``--case``, ``--test``, and ``--suite``.  This section provides details about the relevant options for setting up cases with examples.
Testing will be described in greater detail in the :ref:`testing` section.

``--help``, ``-h`` 
  prints ``icepack.setup`` help information to the terminal and exits.

``--version``
  prints the Icepack version to the terminal and exits.

``--case``, ``-c`` CASE
  specifies the case name.  This can be either a relative path of an absolute path.  This cannot be used with --test or --suite.  Either ``--case``, ``--test``, or ``--suite`` is required.

``--mach``, ``-m`` MACHINE
  specifies the machine name.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  This is required in all modes.

``--env``,  ``-e`` ENVIRONMENT1,ENVIRONMENT2,ENVIRONMENT3
  specifies the environment or compiler associated with the machine.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  Each machine can have multiple supported environments including support for different compilers or other system setups.  When used with ``--suite`` or ``--test``, the ENVIRONMENT can be a set of comma deliminated values with no spaces and the tests will then be run for all of those environments.  With ``--case``, only one ENVIRONMENT should be specified. (default is intel)
  
``--pes``,  ``-p`` MxN
  specifies the number of tasks and threads the case should be run on.  This only works with ``--case``.  The format is tasks x threads or "M"x"N" where M is tasks and N is threads and both are integers. The current icepack driver is purely serial so setting multiple tasks or multiple threads will have no impact.  (default is 1x1)

``--acct``  ACCOUNT
  specifies a batch account number.  This is optional.  See :ref:`account` for more information.

``--grid``, ``-g`` GRID
  specifies the grid.  This is a string and for the current icepack driver, only col is supported. (default = col)

``--set``,  ``-s`` SET1,SET2,SET3
  specifies the optional settings for the case.  This is only used with ``--case`` or ``--test``.  The settings for ``--suite`` are defined in the suite file.  Multiple settings can be specified by providing a comma deliminated set of values without spaces between settings.  The available settings are in **configurations/scripts/options** and ``icepack.setup --help`` will also list them.  These settings files can change either the namelist values or overall case settings (such as the debug flag).

For Icepack, when setting up cases, the ``--case`` and ``--mach`` must be specified.  
It's also recommended that ``--env`` be set explicitly as well.  
At the present time, ``--pes`` and ``--grid`` cannot vary from 1x1 and col respectively
which are the defaults.  ``--acct`` is not normally used.  A more convenient method 
is to use the **~/cice\_proj** file, see :ref:`account`.  The ``--set`` option can be 
extremely handy.  The ``--set`` options are documented in :ref:`settings`.

.. _settings:

Preset Options
~~~~~~~~~~~~~~

There are several preset options.  These are hardwired in 
**configurations/scripts/options** and are specfied for a case or test by 
the ``--set`` command line option.  You can see the full list of settings 
by doing ``icepack.setup --help``.  

The default icepack namelist and icepack settings are specified in the 
files **configuration/scripts/icepack_in** and 
**configuration/scripts/icepack.settings** respectively.  When picking a 
preset setting (option), the set_env.setting and set_nml.setting will be used to 
change the defaults.  This is done as part of the ``icepack.setup`` and the
modifications are resolved in the **icepack.settings** and **icepack_in** file placed in 
the case directory.  If multiple options are chosen and then conflict, then the last
option chosen takes precedent.  Not all options are compatible with each other.

Some of the options are

``debug`` which turns on the compiler debug flags

``short``, ``medium``, ``long`` which change the batch time limit

``diag1`` which turns on diagnostics each timestep

``leap`` which turns on the leap year

``pondcesm``, ``pondlvl``, ``pondtopo`` which turn on the various pond schemes

``run10day``, ``run1year``, etc which specifies a run length

``swccsm3`` which turns on the ccsm3 shortwave and albedo computation

``thermo1`` which on turns on the Bitz-Lipscomb thermodynamics model (default is mushy-layer)

``bgc*`` which turns of various bgc configurations

and there are others.  To add a new option, just add the appropriate file in **configuration/scripts/options**.  Some of the options settings like ``smoke`` and ``restart`` are specifically geared toward setting up tests.  For more information, see :ref:`dev_options`

Examples
~~~~~~~~~

The simplest case is just to setup a default configurations specifying the
case name, machine, and environment::

  icepack.setup --case mycase1 --mach spirit --env intel

To add some optional settings, one might do::

  icepack.setup --case mycase2 --mach spirit --env intel --set debug,diag1,run1year,pondtopo

Once the cases are created, users are free to modify the icepack.settings and icepack_in namelist to further modify their setup.

.. _porting:

Porting
-------

To port, an **env.[machine]_[environment]** and **Macros.[machine]_[environment}** file have to be added to the
**configuration/scripts/machines/** directory and the 
**configuration/scripts/icepack.batch.csh** file needs to be modified.
In general, the machine is specified in ``icepack.setup`` with ``--mach``
and the environment (compiler) is specified with ``--env``.
 
- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files

- cd .. to **configuration/scripts/**

- Edit the **icepack.batch.csh** script to add a section for your machine 
  with batch settings and job launch settings

- Download and untar a forcing dataset to the location defined by 
  ``ICE_MACHINE_INPUTDATA`` in the env file

In fact, this process almost certainly will require some iteration.  The easiest way 
to carry this out is to create an initial set of changes as described above, then 
create a case and manually modify the **env.[machine]** file and **Macros.[machine]** 
file until the case can build and run.  Then copy the files from the case 
directory back to **configuration/scripts/machines/** and update 
the **configuration/scripts/icepack.batch.csh** file, retest, 
and then add and commit the updated machine files to the repository.

.. _account:

Machine Account Settings
~~~~~~~~~~~~~~~~~~~~~~~~

The machine account default is specified by the variable ``ICE_MACHINE_ACCT`` in 
the **env.[machine]** file.  The easiest way to change a user's default is to 
create a file in your home directory called **.cice\_proj** and add your 
preferred account name to the first line.  
There is also an option (``--acct``) in **icepack.setup** to define the account number.  
The order of precedent is **icepack.setup** command line option, 
**.cice\_proj** setting, and then value in the **env.[machine]** file.

Forcing data
------------

The code is currently configured to run in standalone mode on a 4-cell grid using 
atmospheric data, available as detailed in :ref:`testforce` and on the `wiki <https://github.com/CICE-Consortium/Icepack/wiki/Testing-Icepack>`_.
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

The **icepack.setup** script creates a case directory.  However, the model 
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
