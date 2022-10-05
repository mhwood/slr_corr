cd ../
mkdir run_1_yr_balanced_precip
rm -r run_1_yr_balanced_precip/mnc_0001
rm -r run_1_yr_balanced_precip/mnc_0002
rm run_1_yr_balanced_precip/*
cd run_1_yr_balanced_precip
ln -s ../input_1_yr_balanced_precip/* .
ln -s ../data/* .
ln -s ../run_spinup/pickup.0000014600* .
ln -s ../run_spinup/pickup_cd.0000014600* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv
