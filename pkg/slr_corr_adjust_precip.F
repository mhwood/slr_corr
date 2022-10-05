C----&------------------------------------------------------------------xxxxxxx|
C                           slr_corr_adjust_precip
C----&------------------------------------------------------------------xxxxxxx|

#include "SLR_CORR_OPTIONS.h"
#include "EXF_OPTIONS.h"

C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: SLR_CORR_ADJUST_PRECIP

C !INTERFACE:
      SUBROUTINE SLR_CORR_ADJUST_PRECIP( myTime, myIter, myThid,
     &                                   precipArr)

C     !DESCRIPTION:
C     Initialize SLR_CORR parameters and constants.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

#include "EXF_FIELDS.h"

#include "SLR_CORR_PARAM.h"
#include "SLR_CORR_FIELDS.h"

C     !INPUT PARAMETERS:
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      _RL precipArr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
      INTEGER, PARAMETER :: debug = 1
      INTEGER n_average_points
      _RL max_obs_time
      _RL mean_EtaN_target
      INTEGER valid_pts
      INTEGER count0, count1
      INTEGER i, j, bi, bj
      _RL volume_above_zero_target
      _RL volume_above_zero
      _RL volume_above_zero_difference
      _RL precip_volume_flux
      _RL evap_volume_flux
      _RL wet_area
      _RL precip_adjustment

C     TMP LOCAL VARS FOR NOW TO CALL THE MON FUNCTION FOR ETA MEAN
      _RL theMin, theMax, theMean, theSD, theDel2, theVol

CEOP
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOC

#ifdef ALLOW_SLR_CORR

C----&------------------------------------------------------------------xxxxxxx|
C     Define the slr_corr variables
C----&------------------------------------------------------------------xxxxxxx|
    
      if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  BEGIN SLR_CORR_ADJUST_PRECIP DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

      _BEGIN_MASTER(myThid)

C----&------------------------------------------------------------------xxxxxxx|
C   Calculate the mean etan target
C----&------------------------------------------------------------------xxxxxxx|

C     Approach for now is to average available points in the future
C     Sort of like a forward difference
C     Should consider centered or backward averaging later

C     First get the lower and upper bound steps in the obs model
      max_obs_time = slrc_obs_period*slrc_n_obs
      count0=1
      count1=slrc_n_obs
      DO i=1,slrc_n_obs
         if (
     &   (myTime.ge.slrc_obs_start_time+(i-1)*slrc_obs_period)
     &   .and.(i.ge.count0)) then
         count0 = i
         endif
      ENDDO
      count1 = count0 +1

C     Next get the number of steps which are used for the smoothing
      n_average_points = CEILING(slrc_balancePeriod/slrc_obs_period)

C     Now, find the mean EtaN in the average period
C     This will be the target EtaN to adjust the precip toward
      mean_EtaN_target = 0.0
      valid_pts = 0
      DO i=1,n_average_points
            if ((count0+i-1).le.slrc_n_obs) then
              mean_EtaN_target = mean_EtaN_target
     &                         + slrc_obs_timeseries(count0+i-1)
              valid_pts = valid_pts +1
            endif
      ENDDO
      mean_EtaN_target = mean_EtaN_target/valid_pts
C      mean_EtaN_target = slrc_obs_timeseries(count0+i-1)

      PRINT *, 'slrc_obs_timeseries',slrc_obs_timeseries
      PRINT *, 'myTime', myTime
      PRINT *, 'lTime',slrc_obs_start_time+(count0-1)*slrc_obs_period
      PRINT *, 'uTime',slrc_obs_start_time+(count1-1)*slrc_obs_period
      PRINT *, 'count0',count0
      PRINT *, 'count1',count1
      PRINT *, 'n_average_points',n_average_points
      PRINT *, 'mean_EtaN_target',mean_EtaN_target

C----&------------------------------------------------------------------xxxxxxx|
C     Calculate the volume flux deficit
C----&------------------------------------------------------------------xxxxxxx|

C     This is all done in only 1 proc
C     Need to update this to do exchanges for use in MPI

C     Get the volume above 0
      volume_above_zero = 0.0
      precip_volume_flux = 0.0
      evap_volume_flux = 0.0
      wet_area = 0.0
      DO bi=1,nSx
      DO bj=1,nSy
      DO i=1,sNx
      DO j=1,sNy
            volume_above_zero = volume_above_zero
     &       + EtaN(i,j,bi,bj) * rA(i,j,bi,bj) * maskInC(i,j,bi,bj)
            precip_volume_flux = precip_volume_flux
     &       + precip(i,j,bi,bj) * rA(i,j,bi,bj) * maskInC(i,j,bi,bj)
            evap_volume_flux = evap_volume_flux
     &       + evap(i,j,bi,bj) * rA(i,j,bi,bj) * maskInC(i,j,bi,bj)
            wet_area = wet_area + rA(i,j,bi,bj) * maskInC(i,j,bi,bj)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO

      volume_above_zero_target = mean_EtaN_target*wet_area

C     Calculate the precip adjustment (m3/s)
      precip_adjustment = 
     & (volume_above_zero_target-volume_above_zero)/deltaT
     & - precip_volume_flux + evap_volume_flux

C     Convert precip adjustment to a mean per area cell
      precip_adjustment = precip_adjustment/wet_area

      PRINT *, 'wet_area',wet_area
      PRINT *, 'precip_volume_flux',precip_volume_flux
      PRINT *, 'evap_volume_flux',precip_volume_flux
      PRINT *, 'volume_above_zero',volume_above_zero
      PRINT *, 'volume_above_zero_target',volume_above_zero_target
      PRINT *, 'precip_adjustment',precip_adjustment
      PRINT *, 'deltaT',deltaT

      DO bi=1,nSx
      DO bj=1,nSy
      DO i=1,sNx
      DO j=1,sNy
            precipArr(i,j,bi,bj) = precip(i,j,bi,bj)
     C       + precip_adjustment * maskInC(i,j,bi,bj) 
      ENDDO 
      ENDDO
      ENDDO
      ENDDO

C----&------------------------------------------------------------------xxxxxxx|
C   Modify the precip field with new tendency
C----&------------------------------------------------------------------xxxxxxx|

C     Scale the precip to meet this volume flux

C----&------------------------------------------------------------------xxxxxxx|
C     Wrap up the code
C----&------------------------------------------------------------------xxxxxxx|

      _END_MASTER(myThid)

C--   Everyone must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_SLR_CORR */

      if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  END SLR_CORR_ADJUST_PRECIP DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

      RETURN

      END
CEOC
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|