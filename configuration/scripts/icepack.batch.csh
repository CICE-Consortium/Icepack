#! /bin/csh -f

if ( $1 != "" ) then
  echo ${0:t} ${1}
else
  echo ${0:t}
endif

source ./icepack.settings
source ${ICE_CASEDIR}/env.${ICE_MACHCOMP} || exit 2

set jobfile = $1

set ntasks = ${ICE_NTASKS}
set nthrds = ${ICE_NTHRDS}
set maxtpn = ${ICE_MACHINE_TPNODE}
set acct   = ${ICE_ACCOUNT}

@ ncores = ${ntasks} * ${nthrds}
@ taskpernode = ${maxtpn} / $nthrds
@ nnodes = ${ntasks} / ${taskpernode}
if (${nnodes} * ${taskpernode} < ${ntasks}) @ nnodes = $nnodes + 1
set taskpernodelimit = ${taskpernode}
if (${taskpernodelimit} > ${ntasks}) set taskpernodelimit = ${ntasks}
@ corespernode = ${taskpernodelimit} * ${nthrds}

set ptile = $taskpernode
if ($ptile > ${maxtpn} / 2) @ ptile = ${maxtpn} / 2

set shortcase = `echo ${ICE_CASENAME} | cut -c1-15`

#==========================================

cat >! ${jobfile} << EOF0
#!/bin/csh -f 
EOF0

#==========================================

if (${ICE_MACHINE} =~ cheyenne*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
###PBS -m ae 
#PBS -V
#PBS -q ${ICE_MACHINE_QUEUE}
#PBS -N ${ICE_CASENAME}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${corespernode}:mpiprocs=${taskpernodelimit}:ompthreads=${nthrds}
#PBS -l walltime=${ICE_RUNLENGTH}
EOFB

else if (${ICE_MACHINE} =~ hobart*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
###PBS -m ae 
#PBS -V
#PBS -q short
#PBS -N ${ICE_CASENAME}
#PBS -l nodes=1:ppn=24
EOFB

else if (${ICE_MACHINE} =~ izumi*) then
cat >> ${jobfile} << EOFB
#PBS -j oe 
###PBS -m ae 
#PBS -V
#PBS -q ${ICE_MACHINE_QUEUE}
#PBS -N ${ICE_CASENAME}
#PBS -l nodes=1:ppn=48
#PBS -l walltime=${ICE_RUNLENGTH}
EOFB

else if (${ICE_MACHINE} =~ thunder* || ${ICE_MACHINE} =~ gordon* || ${ICE_MACHINE} =~ conrad* || ${ICE_MACHINE} =~ gaffney* || ${ICE_MACHINE} =~ koehr*) then
cat >> ${jobfile} << EOFB
#PBS -N ${shortcase}
#PBS -q ${ICE_MACHINE_QUEUE}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=${ICE_RUNLENGTH}
#PBS -j oe
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ onyx*) then
cat >> ${jobfile} << EOFB
#PBS -N ${ICE_CASENAME}
#PBS -q ${ICE_MACHINE_QUEUE}
#PBS -A ${acct}
#PBS -l select=${nnodes}:ncpus=${maxtpn}:mpiprocs=${taskpernode}
#PBS -l walltime=${ICE_RUNLENGTH}
#PBS -j oe
###PBS -M username@domain.com
###PBS -m be
EOFB

else if (${ICE_MACHINE} =~ cori*) then
@ nthrds2 = ${nthrds} * 2
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
###SBATCH -A ${acct}
#SBATCH --qos shared
#SBATCH --ntasks ${ncores}
#SBATCH --time ${ICE_RUNLENGTH}
#SBATCH --cpus-per-task ${nthrds2}
#SBATCH --constraint haswell
###SBATCH -e filename
###SBATCH -o filename
###SBATCH --mail-type FAIL
###SBATCH --mail-user username@domain.com
EOFB

else if (${ICE_MACHINE} =~ badger*) then
cat >> ${jobfile} << EOFB
#SBATCH -J ${ICE_CASENAME}
#SBATCH -t ${ICE_RUNLENGTH}
#SBATCH -A ${acct}
#SBATCH -N ${nnodes}
#SBATCH -e slurm%j.err
#SBATCH -o slurm%j.out
###SBATCH --mail-type END,FAIL
###SBATCH --mail-user=eclare@lanl.gov
#SBATCH --qos=standby
EOFB

else if (${ICE_MACHINE} =~ loft*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ high_Sierra*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ testmachine*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else if (${ICE_MACHINE} =~ travisCI*) then
cat >> ${jobfile} << EOFB
# nothing to do
EOFB

else
  echo "${0} ERROR: ${ICE_MACHINE} unknown"
  exit -1
endif

exit 0
