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

.. _overview:

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
* ``--suite`` creates a test suite.  Test suites are predefined sets of tests and ``--suite`` provides the ability to quickly setup, build, and run a full suite of tests.

All modes will require use of ``--mach`` or ``-m`` to specify the machine.  Use of ``--env`` is also recommended to specify the compilation environment.  ``--case`` and ``--test`` modes can use ``--set`` or ``-s`` which will turn on various model options.  ``--test`` and ``--suite`` will require ``--testid`` to be set and can use ``--bdir``, ``--bgen``, ``--bcmp``, and ``--diff`` to generate (save) results for regression testing (comparison with prior results). ``--tdir`` will specify the location of the test directory.
Testing will be described in greater detail in the :ref:`testing` section.

Again, ``icepack.setup --help`` will show the latest usage information including 
the available ``--set`` options, the current ported machines, and the test choices.

To create a case, run **icepack.setup**::

  icepack.setup -c mycase -m machine -e intel
  cd mycase

Once a case/test is created, several files are placed in the case directory

- **env.[machine]_[env]** defines the machine environment
- **icepack.settings** defines many variables associated with building and running the model
- **makdep.c** is a tool that will automatically generate the make dependencies
- **Macros.[machine]_[env]** defines the Makefile macros
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
submits the **icepack.run** or **icepack.test** script.  

Some hints:

- To change namelist, manually edit the **icepack_in** file
- To change batch settings, manually edit the top of the **icepack.run** or **icepack.test** (if running a test) file
- When the run scripts are submitted, the current **icepack_in**, **icepack.settings**, and **env.[machine]** files are copied from the case directory into the run directory.  Users should generally not edit files in the run directory as these are overwritten when following the standard workflow.  **icepack.settings** can be sourced to establish the case values in the login shell.  An alias like the following can be established to quickly switch between case and run directories::

    alias  cdrun 'cd `\grep "setenv ICE_RUNDIR"  icepack.settings | awk "{print "\$"NF}"`'
    alias cdcase 'cd `\grep "setenv ICE_CASEDIR" icepack.settings | awk "{print "\$"NF}"`'

- To turn on the debug compiler flags, set ``ICE_BLDDEBUG`` in **icepack.setttings** to true
- To change compiler options, manually edit the Macros file.  To add user defined preprocessor macros, modify ``ICE_CPPDEFS`` in **icepack.settings** using the syntax ``-DCICE_MACRO``.
- To clean the build before each compile, set ``ICE_CLEANBUILD`` in **icepack.settings** to true.  To not clean before the build, set ``ICE_CLEANBUILD`` in **icepack.settings** to false

To build and run::

  ./icepack.build
  ./icepack.submit

The build and run log files will be copied into the logs subdirectory in the case directory.
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

``--setvers``
  Updates the stored value of the Icepack version in the sandbox and exits  See :ref:`version` for more information.

``--docintfc``
  Runs a script that updates the public interfaces in the documentation.  This script parses the source code directly.  See :ref:`docintfc` for more information.

``--case``, ``-c`` CASE
  specifies the case name.  This can be either a relative path of an absolute path.  This cannot be used with --test or --suite.  Either ``--case``, ``--test``, or ``--suite`` is required.

``--mach``, ``-m`` MACHINE
  specifies the machine name.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  This is required in all modes and is paired with ``--env`` to define the compilation environment.

``--env``,  ``-e`` ENVIRONMENT1,ENVIRONMENT2,ENVIRONMENT3
  specifies the compilation environment associated with the machine.  This should be consistent with the name defined in the Macros and env files in **configurations/scripts/machines**.  Each machine can have multiple supported environments including support for different compilers, different compiler versions, different mpi libraries, or other system settigs.  When used with ``--suite`` or ``--test``, the ENVIRONMENT can be a set of comma deliminated values with no spaces and the tests will then be run for all of those environments.  With ``--case``, only one ENVIRONMENT should be specified. (default is intel)
  
``--pes``,  ``-p`` MxN
  specifies the number of tasks and threads the case should be run on.  This only works with ``--case``.  The format is tasks x threads or "M"x"N" where M is tasks and N is threads and both are integers. The current icepack driver is purely serial so setting multiple tasks or multiple threads will have no impact.  (default is 1x1)

``--acct``  ACCOUNT
  specifies a batch account number.  This is optional.  See :ref:`account` for more information.

``--queue`` QUEUE
  specifies a batch queue name.  This is optional.  See :ref:`queue` for more information.

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
the case directory.  If multiple options are chosen that conflict, then the last
option chosen takes precedence.  Not all options are compatible with each other.

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

Once the cases are created, users are free to modify the **icepack.settings** and **icepack_in** namelist to further modify their setup.


.. _cicecpps:

C Preprocessor (CPP) Macros
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a few C Preprocessor Macros supported in the Icepack model.  These
support certain coding features to be excluded or included during the compile.  They
exist in part to support the CICE model and other applications that use Icepack.

For standalone Icepack, The CPPs are defined by the `CPPDEFS` variable in the Icepack
Makefile.  They are defined
by passing the -D[CPP] to the C and Fortran compilers (ie. -DNO_I8) and this
is what needs to be set in the ``CPPDEFS`` variable.  The value of ``ICE_CPPDEFS`` in
**icepack.settings** is copied into the Makefile ``CPPDEFS`` variable as are settings
hardwired into the **Macros.[machine]_[environment]** file.

A list of available CPPs can be found in :ref:`tabcpps`.

.. _version:

Model Version Control
~~~~~~~~~~~~~~~~~~~~~~~~

Managing the internal representation of the model version is handled through the
**icepack.setup** script.  The ``--version`` option displays the version value
on the terminal.  The ``--setvers`` option updates the version defined in the 
sandbox.  It is highly recommended that any changes to the version name be done
through this interface to make sure it's done correctly and comprehensively.
The version name should just include the string associated with the major, minor,
and similar.  For instance,::

  icepack.setup --version

returns

  ./icepack.setup: This is ICEPACK_v1.0.0.d0003

and::

  icepack.setup --setvers v1.0.0.d0004

would update the version.  Always check the string by doing
``icepack.setup --version`` after invoking ``icepack.setup --setvers``.

The version is not updated in the repository unless the code changes associated
with the new version are pushed to the repository.

.. _otherscripts:

Other Scripts Tools
~~~~~~~~~~~~~~~~~~~~~~~~

There are other scripts that come with icepack.  These include

- setup_run_dirs.csh.  This scripts is added to the case directory.  Invoking it creates all the run directories manually.  This script is automatically called as part of the run script, but sometimes it's useful to create these directories before submitting in order to stage custom input files or other data.

.. _porting:

Porting
-------

To port, an **env.[machine]_[environment]** and **Macros.[machine]_[environment]** file have to be added to the
**configuration/scripts/machines/** directory and the 
**configuration/scripts/icepack.batch.csh** file needs to be modified.
In addition **configuration/scripts/icepack.launch.csh** may need to
be modified if simply running the binary directly will not work.
In general, the machine is specified in ``icepack.setup`` with ``--mach``
and the environment (compiler) is specified with ``--env``.  mach and env 
in combination define the compiler, compiler version, supporting libaries,
and batch information.  Multiple compilation environments can be created for
a single machine by choosing unique env names.

- cd to **configuration/scripts/machines/**

- Copy an existing env and a Macros file to new names for your new machine

- Edit your env and Macros files, update as needed

- cd .. to **configuration/scripts/**

- Edit the **icepack.batch.csh** script to add a section for your machine 
  with batch settings and job launch settings

- Edit the **icepack.launch.csh** script to add a section for your machine 
  if executing the binary directly is not supported

- Download and untar a forcing dataset to the location defined by 
  ``ICE_MACHINE_INPUTDATA`` in the env file

In fact, this process almost certainly will require some iteration.  The easiest way 
to carry this out is to create an initial set of changes as described above, then 
create a case and manually modify the **env.[machine]** file and **Macros.[machine]** 
file until the case can build and run.  Then copy the files from the case 
directory back to **configuration/scripts/machines/** and update 
the **configuration/scripts/icepack.batch.csh** file, retest, 
and then add and commit the updated machine files to the repository.

.. _machvars: 

Machine variables
~~~~~~~~~~~~~~~~~~~~~

There are several machine specific variables defined in the **env.$[machine]**.  These
variables are used to generate working cases for a given machine, compiler, and batch
system.  Some variables are optional.

.. csv-table:: *Machine Settings*
   :header: "variable", "format", "description"
   :widths: 15, 15, 25

   "ICE_MACHINE_MACHNAME", "string", "machine name"
   "ICE_MACHINE_MACHINFO", "string", "machine information"
   "ICE_MACHINE_ENVNAME", "string", "env/compiler name"
   "ICE_MACHINE_ENVINFO", "string", "env/compiler information"
   "ICE_MACHINE_MAKE", "string", "make command"
   "ICE_MACHINE_WKDIR", "string", "root work directory"
   "ICE_MACHINE_INPUTDATA", "string", "root input data directory"
   "ICE_MACHINE_BASELINE", "string", "root regression baseline directory"
   "ICE_MACHINE_SUBMIT", "string", "batch job submission command"
   "ICE_MACHINE_TPNODE", "integer", "machine maximum MPI tasks per node"
   "ICE_MACHINE_ACCT", "string", "batch default account"
   "ICE_MACHINE_QUEUE", "string", "batch default queue"
   "ICE_MACHINE_BLDTHRDS", "integer", "number of threads used during build"
   "ICE_MACHINE_QSTAT", "string", "batch job status command (optional)"
   "ICE_MACHINE_QUIETMODE", "true/false", "flag to reduce build output (optional)"

.. _cross_compiling:

Cross-compiling
~~~~~~~~~~~~~~~~~~~~~~~~

It can happen that the model must be built on a platform and run on another, for example when the run environment is only available in a batch queue. The program **makdep** (see :ref:`overview`), however, is both compiled and run as part of the build process.

In order to support this, the Makefile uses a variable ``CFLAGS_HOST`` that can hold compiler flags specfic to the build machine for the compilation of makdep. If this feature is needed, add the variable ``CFLAGS_HOST`` to the **Macros.[machine]_[environment]** file. For example : ::

  CFLAGS_HOST = -xHost

.. _account:

Machine Account Settings
~~~~~~~~~~~~~~~~~~~~~~~~

The machine account default is specified by the variable ``ICE_MACHINE_ACCT`` in 
the **env.[machine]** file.  The easiest way to change a user's default is to 
create a file in your home directory called **.cice\_proj** and add your 
preferred account name to the first line.  
There is also an option (``--acct``) in **icepack.setup** to define the account number.  
The order of precedence is **icepack.setup** command line option, 
**.cice\_proj** setting, and then value in the **env.[machine]** file.

.. _queue:

Machine Queue Settings
~~~~~~~~~~~~~~~~~~~~~~~~

The machine queue default is specified by the variable ``ICE_MACHINE_QUEUE`` in 
the **env.[machine]** file.  The easiest way to change a user's default is to 
create a file in your home directory called **.cice\_queue** and add your 
preferred account name to the first line.  
There is also an option (``--queue``) in **icepack.setup** to define the queue name on a case basis.
The order of precedence is **icepack.setup** command line option, 
**.cice\_queue** setting, and then value in the **env.[machine]** file.

.. _laptops:

Porting to Laptop or Personal Computers
-----------------------------------------
To get the required software necessary to build and run Icepack, a `conda <https://docs.conda.io/en/latest/>`_ environment file is available at :

``configuration/scripts/machines/environment.yml``.

This configuration is supported by the Consortium on a best-effort basis on macOS and GNU/Linux. It is untested under Windows, but might work using the `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_.

Once you have installed Miniconda and created the ``icepack`` conda environment by following the procedures in this section, Icepack should run on your machine without having to go through the formal :ref:`porting` process outlined above.

.. _install_miniconda:

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

We recommend the use of the `Miniconda distribution <https://docs.conda.io/en/latest/miniconda.html>`_ to create a self-contained conda environment from the ``environment.yml`` file.
This process has to be done only once.
If you do not have Miniconda or Anaconda installed, you can install Miniconda by following the `official instructions  <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_, or with these steps:

On macOS:

.. code-block:: bash

  # Download the Miniconda installer to ~/Downloads/miniconda.sh
  curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/Downloads/miniconda.sh
  # Install Miniconda
  bash ~/Downloads/miniconda.sh
  
  # Follow the prompts
  
  # Close and reopen your shell


On GNU/Linux:

.. code-block:: bash

  # Download the Miniconda installer to ~/miniconda.sh
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
  # Install Miniconda
  bash ~/miniconda.sh
  
  # Follow the prompts
  
  # Close and reopen your shell
  

Note: on some Linux distributions (including Ubuntu and its derivatives), the csh shell that comes with the system is not compatible with conda.
You will need to install the tcsh shell (which is backwards compatible with csh), and configure your system to use tcsh as csh:

.. code-block:: bash

  # Install tcsh
  sudo apt-get install tcsh
  # Configure your system to use tcsh as csh
  sudo update-alternatives --set csh /bin/tcsh

.. _init_shell:

Initializing your shell for use with conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend initializing your default shell to use conda.
This process has to be done only once.

The Miniconda installer should ask you if you want to do that as part of the installation procedure.
If you did not answer "yes", you can use one of the following procedures depending on your default shell.
Bash should be your default shell if you are on macOS (10.14 and older) or GNU/Linux.

Note: answering "yes" during the Miniconda installation procedure will only initialize the Bash shell for use with conda.

If your Mac has macOS 10.15 or higher, your default shell is Zsh. 

These instructions make sure that the ``conda`` command is available when you start your shell by modifying your shell's startup file.
Also, they make sure not to activate the "base" conda environment when you start your shell.
This conda environment is created during the Miniconda installation but is not used for Icepack. 

For Bash:

.. code-block:: bash

  # Install miniconda as indicated above, then initialize your shell to use conda:
  source $HOME/miniconda3/bin/activate
  conda init bash
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For Zsh (Z shell):

.. code-block:: bash

  # Initialize Zsh to use conda
  source $HOME/miniconda3/bin/activate
  conda init zsh
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For tcsh:

.. code-block:: bash
  
  # Install miniconda as indicated above, then initialize your shell to use conda:
  source $HOME/miniconda3/etc/profile.d/conda.csh
  conda init tcsh
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For fish:

.. code-block:: bash
  
  # Install miniconda as indicated above, then initialize your shell to use conda:
  source $HOME/miniconda3/etc/fish/conf.d/conda.fish
  conda init fish
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

For xonsh:

.. code-block:: bash

  # Install miniconda as indicated above, then initialize your shell to use conda:
  source-bash $HOME/miniconda3/bin/activate
  conda init xonsh
  
  # Don't activate the "base" conda environment on shell startup
  conda config --set auto_activate_base false
  
  # Close and reopen your shell

.. _init_shell_manually:

Initializing your shell for conda manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you prefer not to modify your shell startup files, you will need to run the appropriate ``source`` command below (depending on your default shell) before using any conda command, and before compiling and running Icepack.
These instructions make sure the ``conda`` command is available for the duration of your shell session.

For Bash and Zsh:

.. code-block:: bash

  # Initialize your shell session to use conda:
  source $HOME/miniconda3/bin/activate

For tcsh:

.. code-block:: bash
  
  # Initialize your shell session to use conda:
  source $HOME/miniconda3/etc/profile.d/conda.csh


For fish:

.. code-block:: bash
  
  # Initialize your shell session to use conda:
  source $HOME/miniconda3/etc/fish/conf.d/conda.fish

For xonsh:

.. code-block:: bash

  # Initialize your shell session to use conda:
  source-bash $HOME/miniconda3/bin/activate


.. _create_conda_env:

Creating Icepack directories and the conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The conda configuration expects some directories and files to be present at ``$HOME/icepack-dirs``:

.. code-block:: bash

  cd $HOME
  mkdir -p icepack-dirs/runs icepack-dirs/baseline icepack-dirs/input
  # Download the required forcing from https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Input-Data
  # and untar it at $HOME/icepack-dirs/input

This step needs to be done only once.

If you prefer that some or all of the Icepack directories be located somewhere else, you can create a symlink from your home to another location:

.. code-block:: bash

  
  # Create the Icepack directories at your preferred location
  cd ${somewhere}
  mkdir -p icepack-dirs/runs icepack-dirs/baseline icepack-dirs/input
  # Download the required forcing from https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Input-Data
  # and untar it at icepack-dirs/input
  
  # Create a symlink to icepack-dirs in your $HOME
  cd $HOME
  ln -s ${somewhere}/icepack-dirs icepack-dirs

Note: if you wish, you can also create a complete machine port for your computer by leveraging the conda configuration as a starting point. See :ref:`porting`.

Next, create the "icepack" conda environment from the ``environment.yml`` file in the Icepack source code repository.  You will need to clone Icepack to run the following command:

.. code-block:: bash

  conda env create -f configuration/scripts/machines/environment.yml

This step needs to be done only once.  If you ever need to update the conda environment
because the required packages change or packages are out of date, do

.. code-block:: bash

  conda env update -f configuration/scripts/machines/environment.yml

.. _using_conda_env:

Using the conda configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the general instructions in :ref:`overview`, using the ``conda`` machine name and ``macos`` or ``linux`` as compiler names.

On macOS:

.. code-block:: bash

  ./icepack.setup -m conda -e macos -c ~/icepack-dirs/cases/case1
  cd ~/icepack-dirs/cases/case1
  ./icepack.build
  ./icepack.run

On GNU/Linux:

.. code-block:: bash

  ./icepack.setup -m conda -e linux -c ~/icepack-dirs/cases/case1
  cd ~/icepack-dirs/cases/case1
  ./icepack.build
  ./icepack.run

A few notes about the conda configuration:

- This configuration always runs the model interactively, such that ``./icepack.run`` and ``./icepack.submit`` are the same.
- You should not update the packages in the ``icepack`` conda environment, nor install additional packages.
- It is not recommeded to run other test suites than ``quick_suite`` or ``travis_suite`` on a personal computer.
- The conda environment is automatically activated when compiling or running the model using the ``./icepack.build`` and ``./icepack.run`` scripts in the case directory. These scripts source the file ``env.conda_{linux.macos}``, which calls ``conda activate icepack``.
- The environment also contains the Sphinx package necessesary to build the HTML documentation. For this use case you must manually activate the environment:

  .. code-block:: bash
  
    cd doc
    conda activate icepack
    make html
    # Open build/html/index.html in your browser
    conda deactivate  # to deactivate the environment

.. _force:

Forcing data
------------

The input data space is defined on a per machine basis by the ``ICE_MACHINE_INPUTDATA`` 
variable in the **env.[machine]** file.  That file space is often shared among multiple 
users, and it can be desirable to consider using a common file space with group read 
and write permissions such that a set of users can update the inputdata area as 
new datasets are available.

The code is currently configured to run in standalone mode on a 4-cell grid using 
atmospheric data, available as detailed on the 
`wiki <https://github.com/CICE-Consortium/Icepack/wiki/Testing-Icepack>`_.
These data files are designed only for testing the code, not for use in production 
runs or as observational data.  Please do not publish results based on these data
sets.  Module **configuration/driver/icedrv\_forcing.F90**
can be modified to change the forcing data. 

Icepack requires near surface atmospheric data at a single point which are set
in ``forcing_nml`` with the ``atm_data_type`` in the namelist (see :ref:`tabsettings`).
The required fields to force icepack include: downwelling long wave and shortwave 
radiative fluxes, latent and sensible heat fluxes, precipitation rate, and near 
surface potential temperature and specific humidity.  The filenames ``atm_data_file``,
``ocn_data_file``, ``ice_data_file``, and ``bgc_data_file``
must also be provided for options other than the default and climatological forcing
cases.  Current filenames can be found in the options scripts in
**configuration/scripts/options** and in the forcing data directories.


1) **Climate Forecast System (CFS)**

   Hourly atmospheric forcing from the National Centers for Environmental Prediction's (NCEP) 
   Climate Forecast System, version 2 (CFSv2) :cite:`Saha14` were utilized to generate
   a one-year time series for Icepack testing. These data were used to create the annual cycle at a 
   point in the Beaufort Sea (70N, 220W) for the period of January 1 00:00UTC - December 31 23:00UTC, 2015. 
   Additional locations can be provided for both hemispheres for the period of 1999-2015 for 
   future testing. This dataset can be used to run for several years to reach equilibrium of the annual
   cycle. 

   Atmospheric forcing fields consist of 2-m air temperature (K), specific humidity (kg/kg),
   10-m wind velocity in the x and y directions (m/s), downward solar radiation (:math:`W/m^2`), 
   downward longwave radiation (:math:`W/m^2`), and precipitation (:math:`kg/m^2/s`). 
   Icepack's boundary layer calculation is used to derive sensible and latent heat fluxes.
   In the namelist, set ``atm_data_type = CFS`` to use CFS atmospheric forcing.


2) **Field campaign derived**

   a) **Norwegian Young Sea Ice cruise (N-ICE)**

    Atmospheric, oceanic, and biogeochemical forcing are available from the 2015 Norwegian Young Sea Ice Cruise 
    (N-ICE) :cite:`Duarte17`. These data are available daily, except for incoming atmospheric radiative forcing,
    which are available 6-hourly. The data correspond to the Arctic Ocean north of Svalbard along the N-ICE drift 
    track (83N, 16E to 80N, 5E) from April 24, 2015 to June 6, 2015.

    Atmospheric forcing fields from :cite:`Duarte17` consist of 2-m air temperature (K), 2-m specific humidity (kg/kg), 
    10-m wind velocity in the x and y directions (m/s), downward solar radiation (:math:`W/m^2`), and precipitation
    (:math:`kg/m^2/s`). Icepack's boundary layer calculation is used to derive sensible and latent heat fluxes. 
    In the namelist, set ``atm_data_type = NICE`` to use N-ICE atmospheric forcing.

    Oceanic forcing fields are available from a Parallel Ocean Program (POP) 1-degree (gx1v3) simulation :cite:`Collins06`.
    These fields consist of sea surface temperature (K), sea surface salinity (ppt), boundary layer depth (m),
    ocean velocity in the x and y direction (m/s), and deep ocean heat flux (:math:`W/m^2`). 
    In the namelist, set ``ocn_data_type = NICE`` to use N-ICE oceanic forcing.

    Biogeochemical forcing fields are available from the World Ocean Atlas :cite:`WOA13`. The biogeochemical fields provided
    are nitrate concentration (:math:`mmol/m^3`) and silicate concentration (:math:`mmol/m^3`). In the namelist, set
    ``bgc_data_type = NICE`` to use N-ICE biogeochemical forcing.

   b) **Ice Station Polarstern (ISPOL)**

    Atmospheric, oceanic, and biogeochemical forcing are available from the 2004 Ice Station Polarstern
    (ISPOL) :cite:`Jeffery14`. These data can be used with both :cite:`Bitz99` and mushy layer thermodynamics. 
    These data are available daily, except for incoming atmospheric radiative forcing,
    which are available 6-hourly. The data correspond to the Weddell Sea (67.9S, 54W) from June 16, 2004 
    to December 31, 2004.

    Atmospheric forcing fields from :cite:`Jeffery14` consist of 2-m air temperature (K), 2-m specific humidity (kg/kg), 10-m wind 
    velocity in the x and y directions (m/s), downward solar radiation (:math:`W/m^2`), and precipitation
    (:math:`kg/m^2/s`). Icepack's boundary layer calculation is used to derive sensible and latent heat fluxes. 
    In the namelist, set ``atm_data_type = ISPOL`` to use ISPOL atmospheric forcing.

    Oceanic forcing fields are available from :cite:`Jeffery14` derived from a POP 1-degree (gx1v3 simulation) :cite:`Collins06`. 
    These consist of sea surface temperature (K), sea surface salinity (ppt), boundary layer depth (m), 
    ocean velocity in the x and y direction (m/s), and deep ocean heat flux (:math:`W/m^2`). 
    In the namelist, set ``ocn_data_type = ISPOL`` to use ISPOL oceanic forcing.

    Biogeochemical forcing fields are available from the World Ocean Atlas :cite:`WOA13`. The biogeochemical fields provided
    are nitrate concentration (:math:`mmol/m^3`) and silicate concentration (:math:`mmol/m^3`). In the namelist, set
    ``bgc_data_type = ISPOL`` to use ISPOL biogeochemical forcing.

   c) **Surface HEat Budget of the Arctic (SHEBA)**

    The ice opening and closing rates (1/s) are derived from the SHEBA data and have been used 
    previously in Cecilia Bitz's column model. For additional information see the following websites:

    - https://atmos.washington.edu/~bitz/column_model/
    - https://atmos.washington.edu/~bitz/column_model/notes_forcing_data

    At present, only the opening and closing rates (1/s) are used from the forcing data. 
    In the namelist, set ``ocn_data_type = SHEBA`` to use this forcing in Icepack.

3) **Climatological** - Maykut and Untersteiner 1971 :cite:`Maykut71`

   The climatological forcing consists of a monthly climatology of downward radiative fluxes, air temperature, 
   relative humidity and wind speed compiled from Arctic ice station observations shown in Table 1 from
   :cite:`Lindsay98`. Icepack's boundary layer calculation is used to derive sensible and latent heat fluxes.  
   The snowfall follows the idealized specification used by :cite:`Semtner76` . 
   To adjust the ice thickness a fixed heating of 6 :math:`W/m^2` is applied to the bottom of the ice.
   This may be thought of as containing about 2 :math:`W/m^2` of ocean heating and an adjustment of 
   about 4 :math:`W/m^2` for biases in the forcings or the model. In the namelist, set ``atm_data_type = clim`` 
   to use climatological atmospheric forcing.


Run Directories
---------------

The **icepack.setup** script creates a case directory.  However, the model 
is actually built and run under the ``ICE_OBJDIR`` and ``ICE_RUNDIR`` directories
as defined in the **icepack.settings** file.  It's important to note that when the
run script is submitted, the current **icepack_in**, **icepack.settings**, and **env.[machine]**
files are copied from the case directory into the run directory.  Users should 
generally not edit files in the run directory as these are overwritten when following
the standard workflow.

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
