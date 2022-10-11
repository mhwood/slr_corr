# global_with_exf_slr

This verification experiment was designed to test the `slr_corr` package. The experiment contains two configurations, with and without the package. In the "uncorrected" experiment, the model is run with a transient sea level rise trend. In the "corrected" experiment, the sea level rise trend is adjusted online with the `slr_corr` package so that the mean sea level is consistent with observations.

The instructions below provide step-by-step instructions to implement the experiments in a fresh clone of MITgcm.

## Adding `slr_corr` to MITgcm
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
cp ../../../pkg/exf/exf_getffields.F
```
Instructions for implementing the modifications to these files is provided on the [mods](https://github.com/mhwood/slr_corr/tree/main/mods) page

Once these changes have been made to the modified code files, the model can be built. Building is OS-dependent but generally follows the following steps:
```
mkdir build
cd build
../../../../tools/genmake2 -mpi -mods ../code -optfile ../../../../tools/build_options/darwin_amd64_gfortran
make depend
make
cd ..
```

