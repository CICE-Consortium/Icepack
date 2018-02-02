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
