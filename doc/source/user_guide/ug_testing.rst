:tocdepth: 3

.. _testing:

Testing Icepack
================

This section documents primarily how to use the Icepack scripts to carry 
out icepack testing.  Exactly what to test is a separate question and
depends on the kinds of code changes being made.  Prior to merging
changes to the CICE Consortium master, changes will be reviewed and
developers will need to provide a summary of the tests carried out.

There is a base suite of tests provided by default with Icepack and this
may be a good starting point for testing.

The testing scripts support several features
 - Ability to test individual (via ``--test``)or multiple tests (via ``--suite``)
   using an input file to define the suite or suites
 - Ability to use test suite defined in the package or test suites defined by the user
 - Ability to store test results for regresssion testing (``--bgen``)
 - Ability to compare results to prior baselines to verify bit-for-bit (``--bcmp``)
 - Ability to define where baseline tests are stored
 - Ability to compare tests against each other (``--diff``)

.. _indtests:

Individual Tests
----------------

The Icepack scripts support both setup of individual tests as well as test suites.  Individual
tests are run from the command line::

  ./icepack.setup --test smoke --mach conrad --env cray --set diag1,debug --testid myid 

Tests are just like cases but have some additional scripting around them.  Individual
tests can be created and manually modified just like cases.
Many of the command line arguments for individual tests
are similar to :ref:`case_options` for ``--case``.  
For individual tests, the following command line options can be set

``--test`` TESTNAME
     specifies the test type.  This is probably either smoke or restart but see `icepack.setup --help` for the latest.  This is required instead of ``--case``.

``--testid`` ID
     specifies the testid.  This is required for every use of ``--test`` and ``--suite``.  This is a user defined string that will allow each test to have a unique case and run directory name.  This is also required.

``--tdir`` PATH
     specifies the test directory.  Testcases will be created in this directory.  (default is .)

``--mach`` MACHINE (see :ref:`case_options`)

``--env`` ENVIRONMENT1 (see :ref:`case_options`)

``--set`` SET1,SET2,SET3 (see :ref:`case_options`)

``--acct`` ACCOUNT (see :ref:`case_options`)

``--grid`` GRID (see :ref:`case_options`)

``--pes`` MxN (see :ref:`case_options`)

Like ``--case``, ``--grid`` and ``--pes`` are not particularly
useful right now within Icepack since the model can only run serially and only
with the col grid setting.  
There are several additional options that come with ``--test`` that are not available
with ``--case`` for regression and comparision testing,

``--bdir`` DIR
     specifies the top level location of the baseline results.  This is used in conjuction with ``--bgen`` and ``--bcmp``.  The default is set by ICE_MACHINE_BASELINE in the env.[machine]_[environment] file.

``--bgen`` DIR
     specifies the name of the directory under [bdir] where test results will be stored.  When this flag is set, it automatically creates that directory and stores results from the test under that directory.  If DIR is set to ``default``, then the scripts will automatically generate a directory name based on the Icepack hash and the date and time.  This can be useful for tracking the baselines by hash.

``--bcmp`` DIR
     specifies the name of the directory under [bdir] that the current tests will be compared to.  When this flag is set, it automatically invokes regression testing and compares results from the current test to those prior results.  If DIR is set to ``default``, then the script will automatically generate the last directory name in the [bdir] directory.  This can be useful for automated regression testing.

``--diff`` LONG_TESTNAME
     invokes a comparison against another local test.  This allows different tests to be compared to each other.  The restrictions are that the test has to already be completed and the testid has to match.

The format of the case directory name for a test will always be 
``[machine]_[env]_[test]_[grid]_[pes]_[sets].[testid]``
The [sets] will always be sorted alphabetically by the script so ``--set debug,diag1`` and
``--set diag1,debug`` produces the same testname and test with _debug_diag1 in that order.

To build and run a test, the process is the same as a case.  cd to the 
test directory, run the build script, and run the submit script::

 cd [test_case]
 ./icepack.build
 ./icepack.submit

The test results will be generated in a local file called **test_output**.
To check those results::

 cat test_output

Tests are defined under **configuration/scripts/tests/**.  The tests currently supported are:

-  smoke   - Runs the model for default length.  The length and options can
            be set with the ``--set`` command line option.  The test passes if the
            model completes successfully.
-  restart - Runs the model for 14 months, writing a restart file at month 3 and
            again at the end of the run.  Runs the model a second time starting from the
            month 3 restart and writing a restart at month 12 of the model run.
            The test passes if both runs complete and
            if the restart files at month 12 from both runs are bit-for-bit identical.

Please run ``./icepack.setup --help`` for the latest information.


.. _examplediff:

Individual Test Examples
~~~~~~~~~~~~~~~~~~~~~~~~

 1) **Basic default single test**
     
    Define the test, mach, env, and testid.
    ::

      ./icepack.setup --test smoke --mach wolf --env gnu --testid t00
      cd wolf_gnu_smoke_col_1x1.t00
      ./icepack.build
      ./icepack.submit
      ./cat test_output


 2) **Simple test with some options**

    Add ``--set``
    ::

      ./icepack.setup --test smoke --mach wolf --env gnu --set diag1,debug --testid t00
      cd wolf_gnu_smoke_col_1x1_debug_diag1.t00
      ./icepack.build
      ./icepack.submit
      ./cat test_output


 3) **Single test, generate baseline dataset**

    Add ``--bgen``
    ::

      ./icepack.setup --test smoke --mach wolf -env gnu --bgen icepack.v01 --testid t00 --set diag1
      cd wolf_gnu_smoke_col_1x1_diag1.t00
      ./icepack.build
      ./icepack.submit
      ./cat test_output


 4) **Single test, compare results to a prior baseline**

    Add ``--bcmp``.  For this to work,
    the prior baseline must exist and have the exact same base testname 
    [machine]_[env]_[test]_[grid]_[pes]_[sets] 
    ::

      ./icepack.setup --test smoke --mach wolf -env gnu --bcmp icepack.v01 --testid t01 --set diag1
      cd wolf_gnu_smoke_col_1x1_diag1.t01
      ./icepack.build
      ./icepack.submit
      ./cat test_output


 5) **Simple test, generate a baseline dataset and compare to a prior baseline**

    Use ``--bgen`` and ``--bcmp``.  The prior baseline must exist already.
    ::

      ./icepack.setup --test smoke --mach wolf -env gnu --bgen icepack.v02 --bcmp icepack.v01 --testid t02 --set diag1
      cd wolf_gnu_smoke_col_1x1_diag1.t02
      ./icepack.build
      ./icepack.submit
      ./cat test_output


 6) **Simple test, comparison against another test**

    Use ``--diff``.  This feature is primarily used in test suites and has 
    limited use in icepack, but is being described for completeness.

    ``--diff`` provides a way to compare tests with each other.  
    For this to work, the tests have to be run in a specific order and
    the testids need to match.  The test 
    is always compared relative to the current case directory.

    To run the first test,
    ::

      ./icepack.setup --test smoke --mach wolf -env gnu --testid tx01 --set debug
      cd wolf_gnu_smoke_col_1x1_debug.tx01
      ./icepack.build
      ./icepack.submit
      ./cat test_output

    Then to run the second test and compare to the results from the first test
    ::

      ./icepack.setup --test smoke --mach wolf -env gnu --testid tx01 --diff smoke_col_1x1_debug
      cd wolf_gnu_smoke_col_1x1.tx01
      ./icepack.build
      ./icepack.submit
      ./cat test_output

    The scripts will add a [machine]_[environment] to the beginning of the diff 
    argument and the same testid to the end of the diff argument.  Then the runs 
    will be compared for bit-for-bit and a result will be produced in test_output.  
    This is really more useful in CICE and for test suites right now.  For example, 
    CICE uses this feature to compare results from different pe counts or 
    decompositions, single threaded vs multi-threaded, and so forth.

.. _testsuites:

Test suites
------------

Test suites support running multiple tests specified via
an input file or files.  When invoking the test suite option (``--suite``) with **icepack.setup**,
all tests will be created, built, and submitted automatically under
a directory called testsuite.[testid].[$date] as part of involing the suite.
Because the tests are built and submitted automatically, 
this feature does not allow for customization of cases or tests like
individual cases and tests do::

  ./icepack.setup --suite base_suite --mach wolf --env gnu --testid myid

Like an individual test, the ``--testid`` option must be specified and can be any 
string.  

If using the ``--tdir`` option, that directory must not exist before the script is run.  The tdir directory will be
created by the script and it will be populated by all tests as well as scripts that support the test suite::

  ./icepack.setup --suite base_suite --mach wolf --env gnu --testid myid --tdir /scratch/$user/testsuite.myid

Once the tests are complete, results can be checked by running the
results.csh script in the [suite_name].[testid]::

  cd testsuite.[testid]
  ./results.csh

The predefined test suites are defined under **configuration/scripts/tests** and 
the files defining the suites
have a suffix of .ts in that directory.  The format for the test suite file 
is relatively simple.  
It is a text file with white space delimited 
columns that define a handful of values in a specific order.  
The first column is the test name, the second the grid, the third the pe count, 
the fourth column is
the ``--set`` options and the fifth column is the ``--diff`` argument. 
(The grid and PEs columns are provided 
for compatibility with the similar CICE scripts.)  The fourth and fifth columns are 
optional.
Lines that begin with # or are blank are ignored.  For example,
::

   #Test   Grid  PEs  Sets                Diff
    smoke   col  1x1  diag1  
    smoke   col  1x1  diag1,run1year  smoke_col_1x1_diag1
    smoke   col  1x1  debug,run1year  
   restart  col  1x1  debug  
   restart  col  1x1  diag1  
   restart  col  1x1  pondcesm  
   restart  col  1x1  pondlvl  
   restart  col  1x1  pondtopo  

The argument to ``--suite`` defines the test suite (.ts) filename or filenames and that argument 
can contain a path.  
**icepack.setup** 
will look for the filename in the local directory, in **configuration/scripts/tests/**, 
or in the path defined by the ``--suite`` option.

Because many of the command line options are specified in the input file, ONLY the
following options are valid for suites,

``--suite`` suitename1,suitename2
  required, input filename with comma delimited list of suite or suites

``--mach`` MACHINE
  required

``--env`` ENVIRONMENT1,ENVIRONMENT2
  strongly recommended

``--acct`` ACCOUNT
  optional

``--tdir`` PATH
  optional

``--testid`` ID
  required

``--bdir`` DIR
  optional, top level baselines directory and defined by default by ICE_MACHINE_BASELINE in **env.[machine]_[environment]**.

``--bgen`` DIR
  recommended, test output is copied to this directory under [bdir]

``--bcmp`` DIR
  recommended, test output are compared to prior results in this directory under [bdir]

``--report``
  This is only used by ``--suite`` and when set, invokes a script that sends the test results to the results page when all tests are complete.  Please see :ref:`testreporting` for more information.

Please see :ref:`case_options` and :ref:`indtests` for more details about how these options are used.


Test Suite Examples
~~~~~~~~~~~~~~~~~~~~~~~~

 1) **Basic test suite**
     
    Specify suite, mach, env, testid.
    ::

      ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a
      cd base_suite.v01a
      #wait for runs to complete
      ./results.csh

 2) **Basic test suite with user defined test directory**
     
    Specify suite, mach, env, testid.
    ::

      ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a --tdir /scratch/$user/ts.v01a
      cd /scratch/$user/ts.v01a
      #wait for runs to complete
      ./results.csh


 3) **Multiple test suites on multiple environments**

      Specify multiple envs.
      ::

        ./icepack.setup --suite base_suite,quick_suite --mach conrad --env cray,pgi,intel,gnu --testid v01a
        cd testsuite.v01a
        #wait for runs to complete
        ./results.csh

      The interface supports both multiple suites and multiple environments from a single
      command line invokation.  Each env or suite can also be run as a separate invokation 
      of `icepack.setup` but if that approach is taken, it is recommended that different testids be used.


 4) **Basic test suite, store baselines in user defined name**

      Add ``--bgen``
      ::

        ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a --bgen icepack.v01a
        cd testsuite.v01a
        #wait for runs to complete
        ./results.csh

      This will store the results in the default [bdir] directory under the subdirectory icepack.v01a.


 5) **Basic test suite, store baselines in user defined top level directory**

      Add ``--bgen`` and ``--bdir``
      ::

        ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a --bgen icepack.v01a --bdir /tmp/user/ICEPACK_BASELINES
        cd testsuite.v01a
        #wait for runs to complete
        ./results.csh

      This will store the results in /tmp/user/ICEPACK_BASELINES/icepack.v01a.


 6) **Basic test suite, store baselines in auto-generated directory**

      Add ``--bgen default``
      ::

        ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a --bgen default
        cd testsuite.v01a
        #wait for runs to complete
        ./results.csh

      This will store the results in the default [bdir] directory under a directory name generated by the script that includes the hash and date.


 7) **Basic test suite, compare to prior baselines**

      Add ``--bcmp``
      ::

        ./icepack.setup --suite base_suite --mach conrad --env cray --testid v02a --bcmp icepack.v01a
        cd testsuite.v02a
        #wait for runs to complete
        ./results.csh

      This will compare to results saved in the baseline [bdir] directory under
      the subdirectory icepack.v01a.  You can use other regression options as well
      (``--bdir`` and ``--bgen``)


 8) **Basic test suite, use of default string in regression testing**

      default is a special argument to ``--bgen`` and ``--bcmp``.  When used, the
      scripts will automate generation of the directories.  In the case of ``--bgen``,
      a unique directory name consisting of the hash and a date will be created.
      In the case of ``--bcmp``, the latest directory in [bdir] will automatically
      be specified.  This provides a number of useful features

       - the ``--bgen`` directory will be named after the hash automatically
       - the ``--bcmp`` will always find the most recent set of baselines
       - the ``--bcmp`` reporting will include information about the comparison directory
         name which will include hash information
       - automation can be invoked easily, especially if ``--bdir`` is used to separate
         results

      Imagine the case where the default settings are used and ``--bdir`` is used to 
      create a unique location.  You could easily carry out regular builds automatically via,
      ::

        set mydate = `date -u "+%Y%m%d"`
        git clone https://github.com/myfork/icepack icepack.$mydate
        cd icepack.$mydate
        ./icepack.setup --suite base_suite --mach conrad --env cray,gnu,intel,pgi --testid $mydate --bcmp default --bgen default --bdir /tmp/work/user/ICEPACK_BASELINES_MASTER

      When this is invoked, a new set of baselines will be generated and compared to the prior
      results each time without having to change the arguments.


 9) **Create and test a custom suite**

      Create your own input text file consisting of 5 columns of data,
       - Test
       - Grid
       - pes
       - sets (optional)
       - diff test (optional)

      such as
      ::

         > cat mysuite
         smoke    col  1x1  diag1,debug
         restart  col  1x1
         restart  col  1x1  diag1,debug    restart_col_1x1
         restart  col  1x1  mynewoption,diag1,debug

      then use that input file, mysuite
      ::

        ./icepack.setup --suite mysuite --mach conrad --env cray --testid v01a --bgen default
        cd mysuite.v01a
        #wait for runs to complete
        ./results.csh

      You can use all the standard regression testing options (``--bgen``, ``--bcmp``, 
      ``--bdir``).  Make sure any "diff" testing that goes on is on tests that
      are created earlier in the test list, as early as possible.  Unfortunately,
      there is still no absolute guarantee the tests will be completed in the correct 
      sequence.


.. _testreporting:

Test Reporting
---------------

The Icepack testing scripts have the capability to post test results
to the official `wiki page <https://github.com/CICE-Consortium/Test-Results/wiki>`_.
You may need write permission on the wiki.  If you are interested in using the
wiki, please contact the consortium.

To post results, once a test suite is complete, run ``results.csh`` and
``report_results.csh`` from the suite directory,
::

  ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a
  cd testsuite.v01a
  #wait for runs to complete
  ./results.csh
  ./report_results.csh

The reporting can also be automated by adding ``--report``
::

  ./icepack.setup --suite base_suite --mach conrad --env cray --testid v01a --report

With ``--report``, the suite will create all the tests, build and submit them,
wait for all runs to be complete, and run the results and report_results scripts.

.. _testplotting:

Test Plotting
----------------

The Icepack scripts include a script (``timeseries.csh``) that will generate a timeseries 
figure from the diagnostic output file.  
When running a test suite, the ``timeseries.csh`` script is automatically copied to the suite directory.  
If the ``timeseries.csh`` script is to be used on a test / case that is not a part of a test suite, 
users will need to run the ``timeseries.csh`` script from the tests directory 
(``./configuration/scripts/tests/timeseries.csh``), or copy it to a local directory and run it 
locally (``cp configuration/scripts/tests/timeseries.csh .`` followed by 
``./timeseries.csh /path/to/ice_diag.full_ITD``. The plotting script can be run
on any of the output files - icefree, slab, full_ITD, land).  To generate the figure, 
run the ``timeseries.csh`` script and pass the full path to the ice_diag file as an argument.  

For example:

Run the test suite. ::

$ ./icepack.setup -m conrad -e intel --suite base_suite -acct <account_number> --testid t00

Wait for suite to finish then go to the directory. ::

$ cd testsuite.t00

Run the timeseries script on the desired case. ::

$ ./timeseries.csh /p/work1/turner/ICEPACK_RUNS/conrad_intel_smoke_col_1x1_diag1_run1year.t00/ice_diag.full_ITD
    
The output figures are placed in the directory where the ice_diag file is located.

This plotting script can be used to plot the following variables:

  - area fraction
  - average ice thickness (m)
  - average snow depth (m)
  - air temperature (C)
  - shortwave radiation (:math:`W/m^2`)
  - longwave radiation (:math:`W/m^2`)
  - snowfall
  - average salinity (ppt)
  - surface temperature (C)
  - outward longwave flux (:math:`W/m^2`)
  - sensible heat flux (:math:`W/m^2`)
  - latent heat flux (:math:`W/m^2`)
  - top melt (m)
  - bottom melt (m)
  - lateral melt (m)
  - new ice (m)
  - congelation (m)
  - snow-ice (m)
  - initial energy change (:math:`W/m^2`)
