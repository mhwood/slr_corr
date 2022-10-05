### Summary of MITgcm Changes for `slr_corr` (

1. `PARAMS.h`: This file is inside the `model/inc` directory. Here, two lines are added for the `slr_corr` package. First `useSlr_corr` is defined as a LOGICAL and then put it in the namelist:
```
1075      LOGICAL useSlr_corr                             
```
and
```
1089     &        useRunClock, useEMBED_FILES, useSlr_corr,           
```

2. `packages_boot.F`: This file is inside the `model/src` directory. For this file, the `slr_corr` package is added in the same location where the `RunClock` package is included. This occurs in three places:
```
95     &          useSlr_corr,        
```
and
```
160      useSlr_corr  = .FALSE.                      
```
and
```
406 #ifdef ALLOW_EMBED_FILES                                         
407       CALL PACKAGES_PRINT_MSG( useSlr_corr,'SLR_CORR', ' ' ) 
408 #endif
```

3. `packages_init_fixed.F`: This file is inside the `model/src` directory. For this file, the `slr_corr` package is added near the end of the file before the default `ALLOW_MYPACKAGE` block:
```
622 #ifdef ALLOW_SLR_CORR
623       IF (useSlr_corr) THEN
624 # ifdef ALLOW_DEBUG
625         IF (debugMode) CALL DEBUG_CALL('SLR_CORR_INIT_FIXED',myThid)
626 # endif
627         CALL SLR_CORR_INIT_FIXED(myThid)
628       ENDIF
629 #endif
```

4. `packages_readparms.F`: This file is inside the `model/src` directory. For this file, the `slr_corr` package is added after the block for the `diagnostics` package:
```
363 #ifdef ALLOW_SLR_CORR
364 C--   if useSlr_corr=T, set slr_corr choices
365 C      otherwise, just set pkgStatus=-1 and return
366        IF (useSlr_corr) CALL SLR_CORR_READPARMS( myThid )
367 #endif /* ALLOW_SLR_CORR */
```

5. `packages_init_variables.F`: This file is inside the `model/src` directory. For this file, the `slr_corr` package is added before the default `ALLOW_MYPACKAGE` block:
```
472 #ifdef ALLOW_SLR_CORR
473       IF ( useSlr_corr ) THEN
474         CALL SLR_CORR_INIT_VARIA( myThid )
475       ENDIF
476 #endif /* ALLOW_SLR_CORR */
```

6. `exf_getffields.F`: This file is inside the `pkg/exf` directory. For this file, a line has been added after the precip fields are read:
```
301       CALL SLR_CORR_ADJUST_PRECIP( myTime, myIter, myThid, precip)
```
This line will most likely need to be moved further down in the script to accomodate additional ctrl adjustments.
