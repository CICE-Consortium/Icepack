:tocdepth: 3

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

.. csv-table:: *Tests*
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
