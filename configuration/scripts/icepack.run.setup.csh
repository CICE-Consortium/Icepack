#!/bin/csh -f

#echo ${0}
echo "running icepack.run.setup.csh"

source ./icepack.settings
source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = icepack.run
set subfile = icepack.submit

set nthrds = ${ICE_NTHRDS}

#==========================================

# Write the batch code into the job file
${ICE_SCRIPTS}/icepack.batch.csh ${jobfile}
if ($status != 0) then
  echo "${0}: ERROR icepack.batch.csh aborted"
  exit -1
endif

#==========================================

cat >> ${jobfile} << EOF1

#--------------------------------------------

cd ${ICE_CASEDIR}
source ./icepack.settings || exit 2
source ./env.\${ICE_MACHCOMP} || exit 2

echo " "
echo "\${0}:"

set  stamp   = \`date '+%y%m%d-%H%M%S'\`
set ICE_RUNLOG_FILE = "icepack.runlog.\${stamp}"

#--------------------------------------------

./setup_run_dirs.csh

#--------------------------------------------
cd \${ICE_RUNDIR}

setenv OMP_NUM_THREADS ${nthrds}

cp -f \${ICE_CASEDIR}/icepack_in \${ICE_RUNDIR}
cp -f \${ICE_CASEDIR}/env.\${ICE_MACHCOMP} \${ICE_RUNDIR}
cp -f \${ICE_CASEDIR}/icepack.settings \${ICE_RUNDIR}
echo " "
echo "ICEPACK rundir is \${ICE_RUNDIR}"
echo "ICEPACK log file is \${ICE_RUNLOG_FILE}"
echo "ICEPACK run started : \`date\`"

EOF1

#==========================================

# Write the job launching logic into the job file
${ICE_SCRIPTS}/icepack.launch.csh ${jobfile}
if ($status != 0) then
  echo "${0}: ERROR icepack.launch.csh aborted"
  exit -1
endif

#==========================================

cat >> ${jobfile} << EOFE
echo "ICEPACK run finished: \`date\`"
echo " "

#--------------------------------------------

if !(-d \${ICE_LOGDIR}) mkdir -p \${ICE_LOGDIR}
cp -p \${ICE_RUNLOG_FILE} \${ICE_LOGDIR}
foreach file (ice_diag.*)
  cp -p \${file} \${ICE_LOGDIR}/\${file}.\${stamp}
end

grep ' ICEPACK COMPLETED SUCCESSFULLY' \${ICE_RUNLOG_FILE}
if ( \$status != 0 ) then
  echo "ICEPACK run did not complete - see \${ICE_LOGDIR}/\${ICE_RUNLOG_FILE}"
  echo "\`date\` \${0}: \${ICE_CASENAME} run did NOT complete \${ICE_RUNLOG_FILE}"  >> \${ICE_CASEDIR}/README.case
  exit -1
endif

echo "\`date\` \${0}: \${ICE_CASENAME} run completed \${ICE_RUNLOG_FILE}"  >> \${ICE_CASEDIR}/README.case
echo "done \${0}"

EOFE

#==========================================

chmod +x ${jobfile}

#==========================================

cat >! ${subfile} << EOFS
#!/bin/csh -f 

${ICE_MACHINE_SUBMIT} ./${jobfile}
echo "\`date\` \${0}: ${ICE_CASENAME} job submitted"  >> ${ICE_CASEDIR}/README.case

EOFS

chmod +x ${subfile}
exit 0
