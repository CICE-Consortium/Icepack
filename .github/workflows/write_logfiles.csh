#!/bin/csh 

#echo "hello"

foreach logfile (case*/logs/icepack.runlog* testsuite.*/*/logs/icepack.runlog*)
  echo "### ${logfile} ###"
  tail -20 $logfile
  echo " "
end

