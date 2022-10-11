C------------------------------------------------------------------------------|
C                           SLR_CORR_FIELDS.h
C------------------------------------------------------------------------------|

#ifdef ALLOW_SLR_CORR

#include "SLR_CORR_SIZE.h"

CBOP
C     !ROUTINE: SLRC_CORR.h
C     !INTERFACE:
C     #include "SLR_CORR.h"

C     !DESCRIPTION:
C     *==========================================================*
C     | SLR_CORR_FIELDS.h
C     | o Header file containing fields to provide adjustments
C     |   to the precip field to balance EtaN online
C     *==========================================================*
CEOP


C------------------------------------------------------------------------------|
C     Create COMMON blocks for the diagnostics_vec variables
C------------------------------------------------------------------------------|

      COMMON /SLR_CORR_FIELDS_R/
     & slrc_obs_timeseries, slr_average_etans,
     & volumes_above_zero
      _RL slrc_obs_timeseries(slrc_n_obs)
      _RL slr_average_etans(slrc_max_average_recs)
      _RL volumes_above_zero(slrc_est_order+1)


#ifdef ALLOW_USE_MPI
      COMMON /SLR_CORR_MPI_FIELDS/
     & procs_volume_above_zero,
     & procs_precip_volume_flux,
     & procs_evap_volume_flux,
     & procs_wet_area
      _RL procs_volume_above_zero(nPx*nPy)
      _RL procs_precip_volume_flux(nPx*nPy)
      _RL procs_evap_volume_flux(nPx*nPy)
      _RL procs_wet_area(nPx*nPy)
#endif


#endif /* ALLOW_SLR_CORR */
