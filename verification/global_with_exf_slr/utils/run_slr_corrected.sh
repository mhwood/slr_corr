cd ../
mkdir run_slr_corrected
rm -r run_slr_corrected/mnc_0001
rm -r run_slr_corrected/mnc_0002
rm run_slr_corrected/*
cd run_slr_corrected
ln -s ../input_slr_corrected/* .
ln -s ../data/* .
ln -s ../run_spinup/pickup.0000014600* .
ln -s ../run_spinup/pickup_cd.0000014600* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv
