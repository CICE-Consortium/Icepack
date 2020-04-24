:tocdepth: 3

.. _dev_scripts:

Scripts Implementation
========================

The scripts are the third part of the icepack package.  They support setting up
cases, building, and running the icepack stand-alone model.

File List
--------------

The directory structure under configure/scripts is as follows.

| **configuration/scripts/**
|        **Makefile**              primary makefile
|        **icepack.batch.csh**     creates batch scripts for particular machines
|        **icepack.build**         compiles the code
|        **icepack.launch.csh**    creates script logic that runs the executable
|        **icepack.run.setup.csh** sets up the run scripts
|        **icepack.run.suite.csh** sets up the test suite
|        **icepack.settings**      defines environment, model configuration and run settings
|        **icepack.test.setup.csh**   creates configurations for testing the model
|        **icepack_decomp.csh**    defines the grid size
|        **icepack_in**            namelist input data
|        **machines/**             machine specific files to set env and Macros
|        **makdep.c**              determines module dependencies
|        **options/**              other namelist configurations available from the icepack.setup command line
|        **parse_namelist.sh**     replaces namelist with command-line configuration
|        **parse_namelist_from_settings.sh**   replaces namelist with values from icepack.settings
|        **parse_settings.sh**     replaces settings with command-line configuration
|        **tests/**                scripts for configuring and running basic tests

.. _dev_strategy:

Strategy
-----------

The icepack scripts are implemented such that everything is resolved after
**icepack.setup** is called.  This is done by both copying specific files
into the case directory and running scripts as part of the **icepack.setup**
command line to setup various files.

**icepack.setup** drives the case setup.  It is written in csh.  All supporting
scripts are relatively simple csh or sh scripts.

The file **icepack.settings** specifies a set of env defaults for the case.  The file
**icepack_in** defines the namelist input for the icepack driver.

.. _dev_options:

Preset Case Options
---------------------


``icepack.setup -s`` option allows the user to choose some predetermined icepack
settings and namelist.  Those options are defined in **configurations/scripts/options/**
and the files are prefixed by either set_env, set_nml, or test_nml.  When **icepack.setup**
is executed, the appropriate files are read from **configurations/scripts/options/**
and the **icepack.settings** and/or **icepack_in** files are updated in the case directory
based on the values in those files.

The filename suffix determines the name of the -s option.  So, for instance, 

  ``icepack.setup -s diag1,debug,bgcISPOL``

will search for option files with suffixes of diag1, debug, and bgcISPOL and then
apply those settings.  

**parse_namelist.sh**, **parse_settings.sh**, and **parse_namelist_from_settings.sh** 
are the three scripts that modify **icepack_in** and **icepack.settings**.

To add new options, just add new files to the **configurations/scripts/options/** directory
with appropriate names and syntax.  The set_nml file syntax is the same as namelist
syntax and the set_env files are consistent with csh setenv syntax.  See other files for
examples of the syntax.

.. _dev_machines:

Machines
-----------

Machine specific information is contained in **configuration/scripts/machines**.  That
directory contains a Macros file and an env file for each supported machine.
For more information on porting to a new machine, see :ref:`porting`.  

.. _dev_testing:

Test scripts
-------------

Under **configuration/scripts/tests** are several files including the scripts to 
setup the various tests, such as smoke and restart tests (**test_smoke.script**, 
**test_restart.script**).
and the files that describe which options files are needed for each test 
(ie. **test_smoke.files**, **test_restart.files**).
A baseline test script (**baseline.script**) is also there to setup the general regression
and comparison testing.  That directory also contains the preset test suites 
(ie. **base_suite.ts**) and a file that supports post-processing on the model
output (**timeseries.csh**).  There is also a script **report_results.csh** that 
pushes results from test suites back to the CICE-Consortium test results wiki page.

To add a new test (for example newtest), several files may be needed,

- **configuration/scripts/tests/test_newtest.script** defines how to run the test.  This chunk
  of script will be incorporated into the case test script
- **configuration/scripts/tests/test_newtest.files** list the set of options files found in
  **configuration/scripts/options/** needed to
  run this test.  Those files will be copied into the test directory when the test is invoked
  so they are available for the **test_newtest.script** to use.
- some new files may be needed in **configuration/scripts/options/**.  These could be 
  relatively generic **set_nml** or **set_env** files, or they could be test specific files 
  typically carrying a prefix of **test_nml**.

Generating a new test, particularly the **test_newtest.script** usually takes some iteration 
before it's working properly.

