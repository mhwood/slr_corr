C CPP options file for SLR_CORR package
C Use this file for selecting options within the SLR_CORR package

#ifndef SLR_CORR_OPTIONS_H
#define SLR_CORR_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_SLR_CORR

C balance sea level rise (each time step by default)
#define ALLOW_SLR_CORR_BALANCE

#endif /* ALLOW_SLR_CORR */
#endif /* ALLOW_SLR_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***