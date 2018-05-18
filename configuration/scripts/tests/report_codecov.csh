#!/bin/csh -f
set use_curl = 1
if (`where curl` == "") then
  echo "use_curl = 0"
  set use_curl = 0
  if (`where wget` == "" ) then
    echo "curl/wget not available"
    (>&2 echo "ERROR: Code coverage reporting (--codecov) needs 'curl' or 'wget' to upload results")
    exit(1)
  endif
endif

# token from https://codecov.io/gh/CICE-Consortium/Icepack/settings
setenv CODECOV_TOKEN "df12b574-8dce-439d-8d3b-ed7428d7598a"

# set name for the report, visible on the codecov user interface
set report_name = "`git rev-parse HEAD`:`git rev-parse --abbrev-ref HEAD` on `hostname`"

# The test-coverage files (*.gcno,*.gcda) must reside next to the source code
# for the coverage reporting to work. However, the coverage files are created
# for each test. For that reason, this script will copy the coverage files over
# to the source directory and report the results, one test at a time. Codecov.io
# is clever enough to report the cumulative results.
echo "Looping over test cases and uploading test coverage"
set testdirs=`ls -d ${ICE_MACHINE_WKDIR}/*`
foreach dir ($testdirs)
  echo "## Submitting results from ${dir}"
  set test_suite_id = "`printf ${dir} | sed 's/_.*\././'`"  # gives <suite>.<id>
  cp $dir/compile/*.{gcno,gcda} ${ICE_SANDBOX}/columnphysics/
  if ( $status == 0 ) then
      echo "Uploading coverage results to codecov.io"
      if ( $use_curl == 1 ) then
          bash -c "bash <(curl -s https://codecov.io/bash) -N '${report_name} ${test_suite_id}'"
      else
          bash -c "bash <(wget -O - https://codecov.io/bash) -N '${report_name} ${test_suite_id}'"
      endif
  else
      echo "No coverage files found for this test"
  endif
end
