#!/bin/csh -f

source ./icepack.settings
source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = icepack.test
set subfile = icepack.submit

set nthrds = ${ICE_NTHRDS}

#==========================================

# Print information about this test to stdout
echo ""
echo "Test name    : ${ICE_TESTNAME}"
echo "Test case dir: ${ICE_CASEDIR}"
echo "Test         : ${ICE_TEST}"
echo "BaseGen      : ${ICE_BASEGEN}"
echo "BaseCom      : ${ICE_BASECOM}"

# Create test script that runs icepack.run, and validates
#==========================================

# Write the batch code into the job file
${ICE_SCRIPTS}/icepack.batch.csh ${jobfile}
if ($status != 0) then
  exit -1
endif

cat >> ${jobfile} << EOF2

cd ${ICE_CASEDIR}
source ./icepack.settings || exit 2
source ./env.\${ICE_MACHCOMP} || exit 2

# Check to see if executable exists in ICE_RUNDIR
if ( ! -f ${ICE_RUNDIR}/icepack ) then
  echo "icepack executable does not exist in ${ICE_RUNDIR}  "
  echo "Please run icepack.build before this test."
  exit 99
endif

# Initial test results and Reset test results for rerun
mv -f ${ICE_CASEDIR}/test_output ${ICE_CASEDIR}/test_output.prev
echo "#---" >! ${ICE_CASEDIR}/test_output
cat ${ICE_CASEDIR}/test_output.prev | grep -i "${ICE_TESTNAME} build" >> ${ICE_CASEDIR}/test_output
echo "PEND ${ICE_TESTNAME} run" >> ${ICE_CASEDIR}/test_output
rm -f ${ICE_CASEDIR}/test_output.prev

EOF2

if ( -f ${ICE_SCRIPTS}/tests/test_${ICE_TEST}.script) then
  echo "${0:t} using test_${ICE_TEST}.script"
  cat >> ${jobfile} < ${ICE_SCRIPTS}/tests/test_${ICE_TEST}.script
else
  echo "${0:t} ERROR: ${ICE_SCRIPTS}tests/test_${ICE_TEST}.script not found"
  exit -1
endif
cat >> ${jobfile} < ${ICE_SCRIPTS}/tests/baseline.script

chmod +x ${jobfile}

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${ICE_MACHINE_SUBMIT} ./${jobfile}
echo "\`date\` \${0}: ${ICE_CASENAME} job submitted"  >> ${ICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
exit 0
