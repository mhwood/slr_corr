# global_with_exf_slr

This verification experiment was designed to test the `slr_corr` package. The experiment contains two configurations - an "uncorrected" experiement which doesn't use the package, and a "correctd" experiement which uses package. In the "uncorrected" experiment, the model is run with a transient sea level rise trend which differents from the provided "observations". In the "corrected" experiment, the sea level rise signal is adjusted online with the `slr_corr` package so that the mean sea level is consistent with "observations".

The instructions below provide step-by-step instructions to implement the experiments in a fresh clone of MITgcm. There are three main steps and one optional step:
- Build the model with `slr_corr` in a fresh MITgcm clone
- Run the "uncorrected" experiment
- Run the "corrected" experiment
- Create a comparison plot (optional)

## Building the model with `slr_corr` in a fresh MITgcm clone
Begin by cloning a fresh copy of [MITgcm](https://github.com/MITgcm/MITgcm):
```
git clone https://github.com/MITgcm/MITgcm.git
```
Next, create a directory for the `slr_corr` package and add all of the files from `slr_corr/pkg`:
```
cd MITgcm/pkg
mkdir slr_corr
cp [this_repo]/pkg/* slr_corr
cd ../
```
Next, create a directory for the verification experiment:
```
mkdir configurations
mkdir configurations/slr_corr_test
cd configurations/slr_corr_test
cp -r [this_repo]/verification/global_with_exf_slr .
```
Now, add files specific to your cloned version of MITgcm and edit them for the `slr_corr` package:
```
cd global_with_exf_slr
cp ../../../model/inc/PARAMS.h code
cp ../../../model/src/packages_boot.F code
cp ../../../model/src/packages_init_fixed.F code
cp ../../../model/src/packages_readparms.F code
cp ../../../pkg/exf/exf_getffields.F code
```
Instructions for implementing the modifications to these files are provided on the [mods](https://github.com/mhwood/slr_corr/tree/main/mods) page.

As a final step before builing, choose whether you would like to use mpi or not:
```
cp code/SIZE.h_2_proc code/SIZE.h          # with mpi
cp code/SIZE.h_no_mpi code/SIZE.h          # without mpi
```

Once these changes have been made to the modified code files, the model can be built. Building is OS-dependent but generally follows the following steps:
```
mkdir build
cd build
../../../../tools/genmake2 -mpi -mods ../code -optfile ../../../../tools/build_options/darwin_amd64_gfortran
make depend
make
cd ..
```

## Running the "uncorrected" experiment
Now that the model is built, the "uncorrected" experiment can be prepared and run. First, gather all of the pertinent files from the existing verification experiment:
```
cp ../../../verification/tutorial_global_oce_latlon/input/*.bin data
```
Note that the external forcing conditions from the `global_with_exf` experiment are not necessry to copy - these files are replaced with data provided with this experiment. Now, the model can be run for a 1 year (or longer) simulation:
```
mkdir run_uncorrected
cd run_uncorrected
ln -s ../input/* .
ln -s ../namelist_uncorrected/* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv
```



