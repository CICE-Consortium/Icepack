#!/bin/csh -f

#--- inputs ---

#echo "${0:t}  input ICE_DECOMP_GRID  = $ICE_DECOMP_GRID"

set grid = $ICE_DECOMP_GRID

#--- computation ---

if (${grid} == 'col') then
  set nxglob = 4
else
  echo "${0:t}: ERROR unknown grid ${grid}"
  exit -9
endif

#--- outputs ---

setenv ICE_DECOMP_NXGLOB $nxglob

#echo "${0:t} output ICE_DECOMP_NXGLOB   = $ICE_DECOMP_NXGLOB"

exit 0
