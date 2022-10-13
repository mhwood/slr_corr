# Verification Experiments

At this time, there is only one verification experiment. If more are developed, a brief summary will be added here.

## global_with_exf_slr
This verification experiment follows the [global_with_exf](https://github.com/MITgcm/MITgcm/tree/master/verification/global_with_exf) experiment provided with the main branch of MITgcm. The only difference is that the precipitation and evaporation files are edited to induce a global mean sea level rise. The global signal has a trend and a seasonal cycle, similar to global observations. The verification experiment has two configurations - with and without `slr_corr` - to verify the package works as desired.
