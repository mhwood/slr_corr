# slr_corr: An MITgcm package for online Sea Level Rise CORRections
This MITgcm package was designed to provide online adjustments to model runs such that modeled mean sea level is consistent with a prescribed set of observations.

![Demo of the slr-corr package](https://github.com/mhwood/slr_corr/blob/main/slr_corr_demo.png)

## Repository Overview
This repository contains the following sub-directories:

#### pkg
The pkg directory contains all of the pertinent files required to compile MITgcm with `slr_corr`. These files should be placed in the pkg directory of MITgcm in a subdirectory titled `slr_corr`. For a complete description of the package, see the README file in the pkg directory.

#### mods
The mods directory contains all of the existing MITgcm src or inc files which were modified to implement the `slr_corr` package. A summary of the changes to each file has been added in the README so that similar changes can be implemented in various versions of MITgcm.

#### verification
The verification directory contains experiments that demonstrate the use of `slr_corr`. Instructions for compiling and running each experiment are provided in their respective READMEs.

## Using this package
As for most MITgcm packages, there are both compile-time and run-time considerations for `slr_corr`.

#### Compile-time Steps
1. Ensure that `slr_corr` is added in the MITgcm/pkg directory.
2. Edit default `src` and `inc` files to include `slr_corr` (see instructions in the [mods](https://github.com/mhwood/slr_corr/tree/main/mods) directory).
2. Add `slr_corr` to the packages.conf file in the code directory.
3. Add the SLR_CORR_SIZE.h header file to the code directory and edit to match the observations timeseries size. 

#### Run-time Steps
1. Add `useSlr_corr = .TRUE.` to the the data.pkg file.
2. Provide a `data.slr_corr` file for run time.
3. Provide a timeseries of mean sea level for run time.

See the [verification](https://github.com/mhwood/slr_corr/tree/main/verification) experiments for a demo of these steps.
