#!/bin/csh -f

#echo ${0}
echo "running icepack.launch.csh"

source ./icepack.settings

set jobfile = $1

#==========================================

if (${ICE_MACHINE} =~ gordon* || ${ICE_MACHINE} =~ conrad* || ${ICE_MACHINE} =~ onyx*  || ${ICE_MACHINE} =~ daley* || ${ICE_MACHINE} =~ banting* ) then
cat >> ${jobfile} << EOFR
aprun -n 1 -N 1 -d 1 ./icepack >&! \$ICE_RUNLOG_FILE
EOFR

#==========================================

else if (${ICE_MACHINE} =~ perlmutter) then
cat >> ${jobfile} << EOFR
srun --cpu-bind=cores ./icepack >&! \$ICE_RUNLOG_FILE
EOFR

#==========================================
else
cat >> ${jobfile} << EOFR
./icepack >&! \$ICE_RUNLOG_FILE
EOFR

endif

exit 0
