# global_with_exf_slr

This verification experiment was designed to test the `slr_corr` package. The experiment contains two configurations - an "uncorrected" experiement which doesn't use the package, and a "corrected" experiement which uses package. In the "uncorrected" experiment, the model is run with a transient sea level rise trend which differs from the provided "observations". In the "corrected" experiment, the sea level rise signal is adjusted online with the `slr_corr` package so that the mean sea level is consistent with "observations".

The instructions below provide step-by-step instructions to implement the experiments in a fresh clone of MITgcm. The instructions below outline how to:
- Build the model with `slr_corr` in a fresh MITgcm clone
- Run the "uncorrected" experiment (optional)
- Run the "corrected" experiment using `slr_corr`
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

As a final step before builing, choose whether you would like to use mpi or not. In this example, we will use mpi but the code will work without it (using 2 threads). To configure with mpi, use the SIZE.h file with 2 procs:
```
cp code/SIZE.h_2_proc code/SIZE.h
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
cp ../../../verification/tutorial_global_oce_latlon/input/*.bin input
```
Note that the external forcing conditions from the `global_with_exf` experiment are not necessry to copy - these files are replaced with data provided with this experiment. In particular, the precipitation and evaporation fields were created using the mean seasonal cycle from ECCOv4r4 data with an added transient sea level induced by tuning the precipitation field.

Now, the model can be run for a 1 year (or longer) simulation starting in 1992:
```
mkdir run_uncorrected
cd run_uncorrected
ln -s ../input/* .
ln -s ../namelist_uncorrected/* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv
cd ..
```


## Running the "corrected" experiment
The steps to run the "corrected" experiment are similar to the "uncorrected" experiment. If you did not run the "uncorrected" experiment, be sure to grab the input binaries from the tutorial as described above. Running the "corrected" experiement can be done with:
```
mkdir run_corrected
cd run_corrected
ln -s ../input/* .
ln -s ../namelist_corrected/* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv
cd ..
```

## Plotting a comparison of mean sea level
To plot a comparison of the "uncorrected" and "corrected" model results along with the "observations", there is a convenient Python script provided in the `utils` directory:
```
cd utils
python3 plot_mean_EtaN_comparison.py -d ../
```
This script will generate the plot shown on the front page of this repo.
