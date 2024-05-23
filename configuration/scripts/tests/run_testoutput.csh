#/bin/csh -f

./icepack.setup -m conda -e macos --case basecase -s diag1,run1year
cd basecase
./icepack.build
./icepack.submit
cd logs
foreach file (ice*[0-9])
  cp ${file} ${file}.txt
end
conda activate icepack
../../configuration/scripts/tests/timeseries.csh ice_diag.full_ITD.*[0-9]

echo "run_testoutput DONE"
