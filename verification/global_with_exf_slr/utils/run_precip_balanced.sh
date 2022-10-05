cd ../
mkdir run_precip_balanced
rm -r run_precip_balanced/mnc_0001
rm -r run_precip_balanced/mnc_0002
rm run_precip_balanced/*
cd run_precip_balanced
ln -s ../input_precip_balanced/* .
ln -s ../data/* .
ln -s ../run_spinup/pickup.0000014600* .
ln -s ../run_spinup/pickup_cd.0000014600* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv
