# global_with_exf_slr

This verification experiment was designed to test the `slr_corr` package. The experiment contains two configurations, with and without the package. In the "uncorrected" experiment, the model is run with a transient sea level rise trend. In the "corrected" experiment, the sea level rise trend is adjusted online with the `slr_corr` package so that the mean sea level is consistent with observations.

The instructions below provide step-by-step instructions to implement the experiments in a fresh clone of MITgcm.

## Adding `slr_corr` to MITgcm
Begin by cloning a fresh copy of [MITgcm](https://github.com/MITgcm/MITgcm):
```
git clone https://github.com/MITgcm/MITgcm.git
```
Next, create a directory for the `slr_corr` package and add all of the files from `slr_corr/pkg':
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


