cd ../
mkdir run_slr_uncorrected
rm -r run_slr_uncorrected/mnc_0001
rm -r run_slr_uncorrected/mnc_0002
rm run_slr_uncorrected/*
cd run_slr_uncorrected
ln -s ../input_slr_uncorrected/* .
ln -s ../data/* .
ln -s ../run_spinup/pickup.0000014600* .
ln -s ../run_spinup/pickup_cd.0000014600* .
ln -s ../build/mitgcmuv .
mpirun -np 1 ./mitgcmuv
