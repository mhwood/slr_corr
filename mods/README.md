# Summary of MITgcm modifications for `slr_corr`
The files in this directory and this summary are provided because these files undergo routine changes as MITgcm evolves. By implementing these changes manually, one can ensure that the version of MITgcm being used will accomodate the `slr_corr` package and compile without error.

## PARAMS.h
This file is inside the `model/inc` directory. Here, two lines are added for the `slr_corr` package. First `useSlr_corr` is defined as a LOGICAL and then put it in the namelist:
```
1075      LOGICAL useSlr_corr                             
```
and
```
1089     &        useRunClock, useEMBED_FILES, useSlr_corr,           
```

## packages_boot.F
This file is inside the `model/src` directory. For this file, the `slr_corr` package is added in the same location where the `RunClock` package is included. This occurs in three places:
```
95     &          useSlr_corr,        
```
and
```
160      useSlr_corr  = .FALSE.                      
```
and
```
406 #ifdef ALLOW_SLR_CORR                                        
407       CALL PACKAGES_PRINT_MSG( useSlr_corr,'SLR_CORR', ' ' ) 
408 #endif
```

## packages_init_fixed.F
This file is inside the `model/src` directory. For this file, the `slr_corr` package is added near the end of the file before the default `ALLOW_MYPACKAGE` block:
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

## packages_readparms.F
This file is inside the `model/src` directory. For this file, the `slr_corr` package is added after the block for the `diagnostics` package:
```
363 #ifdef ALLOW_SLR_CORR
364 C--   if useSlr_corr=T, set slr_corr choices
365 C      otherwise, just set pkgStatus=-1 and return
366        IF (useSlr_corr) CALL SLR_CORR_READPARMS( myThid )
367 #endif /* ALLOW_SLR_CORR */
```

## exf_getffields.F
This file is inside the `pkg/exf` directory. For this file, a line has been added after the precip fields are read:
```
301 #ifdef ALLOW_SLR_CORR
302       IF (useSlr_corr) THEN
303       CALL SLR_CORR_ADJUST_PRECIP( myTime, myIter, myThid, precip)
304       ENDIF
305 #endif
```
This line will most likely need to be moved further down in the script to accomodate additional ctrl adjustments.
