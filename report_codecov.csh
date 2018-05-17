#!/bin/csh -f
if (`where curl` == "") then
  (>&2 echo "ERROR: Code coverage reporting (--codecov) needs 'curl' to upload results")
  exit(1)
endif

echo "Looping over test cases and uploading test coverage"
foreach dir ("${ICE_MACHINE_WKDIR}"/*/)
  echo "## Submitting results from ${dir}"
  cp "$dir"/compile/*.{gcno,gcda} "${ICE_SANDBOX}"/columnphysics/ &&
                     bash <(curl -s https://codecov.io/bash)
done
