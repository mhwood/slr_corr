C----&------------------------------------------------------------------xxxxxxx|
C                           slr_corr_adjust_precip
C----&------------------------------------------------------------------xxxxxxx|

#include "SLR_CORR_OPTIONS.h"
#include "EXF_OPTIONS.h"
#include "ECCO_OPTIONS.h"

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
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EXF_FIELDS.h"
#include "ECCO_SIZE.h"
#include "ECCO.h"

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
      INTEGER i, j, bi, bj, n
      INTEGER kSrf
      _RL volume_above_zero_target
      _RL volume_above_zero
      _RL volume_above_zero_difference
      _RL precip_volume_flux
      _RL evap_volume_flux
      _RL wet_area, corr_area
      _RL next_volume_above_zero
      _RL precip_adjustment
      _RL alpha

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
C   Calculate the mean EtaN target from provided observations
C----&------------------------------------------------------------------xxxxxxx|

C     First get the lower and upper bound steps in the obs
C     timeseries
      max_obs_time = slrc_obs_period*slrc_n_obs
      count0=1
      count1=slrc_n_obs
      DO i=1,slrc_n_obs
         if (
     &   ((myTime+deltaT) .ge. 
     &     slrc_obs_start_time+(i-1)*slrc_obs_period)
     &   .and.(i.ge.count0)) then
         count0 = i
         endif
      ENDDO
      count1 = count0 +1

C     Now, interpolate to get the target EtaN
C     This is the target EtaN to adjust the precip toward
      alpha = (myTime+deltaT
     &       -(slrc_obs_start_time+(count0-1)*slrc_obs_period))
     &       /  (slrc_obs_period)
      mean_EtaN_target = (1-alpha) * slrc_obs_timeseries(count0) 
     &                   + alpha * slrc_obs_timeseries(count1)
     &                  + slrc_obs_ref

C      PRINT *, 'slrc_obs_timeseries',slrc_obs_timeseries
C      PRINT *, 'myTime', myTime
C      PRINT *, 'lTime',slrc_obs_start_time+(count0-1)*slrc_obs_period
C      PRINT *, 'uTime',slrc_obs_start_time+(count1-1)*slrc_obs_period
C      PRINT *, 'count0',count0
C      PRINT *, 'count1',count1
C      PRINT *, 'alpha', alpha
C      PRINT *, 'lEtaN',slrc_obs_timeseries(count0)
C      PRINT *, 'uEtaN',slrc_obs_timeseries(count1)
C      PRINT *, 'mean_EtaN_target',mean_EtaN_target

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,I10)') "  myIter ", myIter
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3)') "  myTime ", myTime
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I5)') "  count0 ", count0
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I5)') "  count1 ", count1
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3)') "  alpha ", alpha
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3)') "  mean_EtaN_target ", mean_EtaN_target
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C----&------------------------------------------------------------------xxxxxxx|
C     Calculate the volume flux deficit/excess
C----&------------------------------------------------------------------xxxxxxx|

      volume_above_zero = 0.0
      precip_volume_flux = 0.0
      evap_volume_flux = 0.0
      wet_area = 0.0
      corr_area = 0.0

      IF ( fluidIsAir ) THEN
       kSrf = 0
      ELSEIF ( usingPCoords ) THEN
       kSrf = Nr
      ELSE
       kSrf = 1
      ENDIF

#ifndef ALLOW_USE_MPI
C     Without MPI, we can just calculate these values with what we know
      DO bi=1,nSx
      DO bj=1,nSy
      DO i=1,sNx
      DO j=1,sNy
            volume_above_zero = volume_above_zero
#ifdef ALLOW_ECCO
#if (defined ATMOSPHERIC_LOADING && defined ALLOW_IB_CORR)
     &       + m_eta_dyn(i,j,bi,bj)
#else
     &       + m_eta(i,j,bi,bj)
#endif
#else
     &       + EtaN(i,j,bi,bj)
#endif
     &       * rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            precip_volume_flux = precip_volume_flux
     &       + precip(i,j,bi,bj) * rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            evap_volume_flux = evap_volume_flux
     &       + evap(i,j,bi,bj) * rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            wet_area = wet_area + rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
#ifdef MODIDY_POSITIVE_PRECIP_ONLY
            IF(precip(i,j,bi,bj).GT.0. _d 0)THEN
             corr_area = corr_area + rA(i,j,bi,bj)*maskC(i,j,kSrf,bi,bj)
            ENDIF
#endif
      ENDDO 
      ENDDO
      ENDDO
      ENDDO
#else
C     When using MPI, proc 0 keeps track of everyone's values
      CALL SEND_FLUXES_TO_PROC_0(myThid)

C      PRINT *,'procs_volume_above_zero',procs_volume_above_zero
C      PRINT *,'procs_precip_volume_flux',procs_precip_volume_flux
C      PRINT *,'procs_evap_volume_flux',procs_evap_volume_flux
C      PRINT *,'procs_wet_area',procs_wet_area

      DO i=1,nPx*nPy
            volume_above_zero = volume_above_zero
     &       + procs_volume_above_zero(i)
            precip_volume_flux = precip_volume_flux
     &       + procs_precip_volume_flux(i)
            evap_volume_flux = evap_volume_flux
     &       + procs_evap_volume_flux(i)
            wet_area = wet_area + procs_wet_area(i)
            corr_area = corr_area + procs_corr_area(i)
      ENDDO
#endif

      volume_above_zero_target = mean_EtaN_target*wet_area

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,F20.3,A)') "  volume_above_zero", 
     &                          volume_above_zero, " m^3"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3,A)') "  volume_above_zero_target", 
     &                          volume_above_zero_target, " m^3"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3,A)') "  precip_volume_flux", 
     &                          precip_volume_flux, " m^3/s"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3,A)') "  evap_volume_flux", 
     &                          evap_volume_flux, " m^3/s"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3,A)') "  wet_area", 
     &                          wet_area, " m^3"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3,A)') "  corr_area",
     &                          corr_area, " m^2"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C----&------------------------------------------------------------------xxxxxxx|
C     Estimate what the next volume above 0 will be
C----&------------------------------------------------------------------xxxxxxx|

C     Deactivated this section for now
C     This idea is to estimate the next target volume based on the 
C     past few values
C     This code could be used for improvements later

C     Here, we will estimate V at the next timestep (V(i+1))
C     based on the current timestep (V(i))
C     i.e. V(i+1) = V(i) + dV/dt(i) * deltaT
C     dV/dt can be calculated with a backward difference or higher 
C     order versions of this 
C     the variable slrc_est_order is used to set the order of the approximation
C     the default is set to 2 in SLR_CORR_SIZE (linear approximation)

C      DO n=1,slrc_est_order
C      volumes_above_zero(n+1) = volumes_above_zero(n)
C      ENDDO
C      volumes_above_zero(1) = volume_above_zero

C      IF (myIter.eq.nIter0) then
C      DO n=1,slrc_est_order
C      volumes_above_zero(n+1) = volume_above_zero
C      ENDDO
C      ENDIF

C     Linear approximation (backwards difference)
C      next_volume_above_zero = volume_above_zero
C     & + volumes_above_zero(1)-volumes_above_zero(2)

C      PRINT *, 'next_volume_above_zero',next_volume_above_zero
C      PRINT *, 'volumes_above_zero(1)',volumes_above_zero(1)
C      PRINT *, 'volumes_above_zero(2)',volumes_above_zero(2)
C      PRINT *, 'next/current',next_volume_above_zero/volume_above_zero

C----&------------------------------------------------------------------xxxxxxx|
C     Calculate the precip adjustment
C----&------------------------------------------------------------------xxxxxxx|

C     Calculate the precip adjustment (m3/s)
      precip_adjustment = 
     & (volume_above_zero_target-volume_above_zero)/deltaT
C     & - precip_volume_flux + evap_volume_flux

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,F20.3,A)') "  precip_adjustment", 
     &                          precip_adjustment, " m^3/s"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C     Convert precip adjustment to a mean per area cell
C     Total area (denominator) is wet_area by default and corr_area if 
C     MODIDY_POSITIVE_PRECIP_ONLY is defined in EXF_OPTIONS.h
#ifdef MODIDY_POSITIVE_PRECIP_ONLY
      precip_adjustment = precip_adjustment/corr_area
#else
      precip_adjustment = precip_adjustment/wet_area
#endif

#ifndef ALLOW_USE_MPI
      slrc_precip_adjustment = precip_adjustment
#else
C     With MPI, only proc 0 has added up all the precips from
C     the other additions
      if (mpiMyId.eq.0) then 
      slrc_precip_adjustment = precip_adjustment
      else
      slrc_precip_adjustment = 0.0
      endif
#endif

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,F20.3,A)') "  precip_adjustment", 
     &                          precip_adjustment, " m/s"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3,A)') "  deltaT", 
     &                          deltaT, " s"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C----&------------------------------------------------------------------xxxxxxx|
C     Apply the precip adjustment to the precip grid
C----&------------------------------------------------------------------xxxxxxx|

#ifdef ALLOW_USE_MPI
C     When using MPI, proc 0 has to send out the precip adjustment
C     to everyone else before its applied
      CALL SEND_ADJUSTMENT_TO_ALL_PROCS(myThid)
#endif

C      print *,'slrc_precip_adjustment',slrc_precip_adjustment

      DO bi=1,nSx
      DO bj=1,nSy
      DO i=1,sNx
      DO j=1,sNy
            precipArr(i,j,bi,bj) = precip(i,j,bi,bj)
     &       + slrc_precip_adjustment * maskC(i,j,kSrf,bi,bj) 
      ENDDO 
      ENDDO
      ENDDO
      ENDDO


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











#ifdef ALLOW_USE_MPI
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: SUBROUTINE SEND_FLUXES_TO_PROC_0

C !INTERFACE:
      SUBROUTINE SEND_FLUXES_TO_PROC_0(myThid)

C     !DESCRIPTION:
C     Sending surface fluxes from processes > 0 to master process 0 
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EXF_FIELDS.h"
#include "ECCO_SIZE.h"
#include "ECCO.h"
#include "SLR_CORR_PARAM.h"
#include "SLR_CORR_FIELDS.h"

C     !LOCAL VARIABLES:
      INTEGER vec_id, myThid
      INTEGER pid, status(MPI_STATUS_SIZE), ierror
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER, PARAMETER :: debug = 0
      INTEGER i, j, bi, bj
      INTEGER kSrf

      _RL proc_volume_above_zero
      _RL proc_precip_volume_flux
      _RL proc_evap_volume_flux
      _RL proc_wet_area
      _RL proc_corr_area

      _RL vol_tmp
      _RL precip_tmp
      _RL evap_tmp
      _RL area_tmp
      _RL corr_area_tmp
CEOP
C----&------------------------------------------------------------------xxxxxxx|
CBOC 

C     /==================================================================\
C     | Note: MPI_SEND params: data_to_send, send_count, send_type,      |
C     |                        destination_ID, tag, comm, ierror         |
C     \==================================================================/

C     /==================================================================\
C     | Note: MPI_RECV params: received_data, receive_count,             |
C     |                        receive_type, sender_ID, tag, comm,       |
C     |                        status, ierror                            |
C     \==================================================================/

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SEND_FLUXES_TO_PROC_0 DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      IF ( fluidIsAir ) THEN
       kSrf = 0
      ELSEIF ( usingPCoords ) THEN
       kSrf = Nr
      ELSE
       kSrf = 1
      ENDIF

      proc_volume_above_zero = 0.0
      proc_precip_volume_flux = 0.0
      proc_evap_volume_flux = 0.0
      proc_wet_area = 0.0
      proc_corr_area = 0.0
      DO bi=1,nSx
      DO bj=1,nSy
      DO i=1,sNx
      DO j=1,sNy
            proc_volume_above_zero = proc_volume_above_zero
#ifdef ALLOW_ECCO
#if (defined ATMOSPHERIC_LOADING && defined ALLOW_IB_CORR)
     &       + m_eta_dyn(i,j,bi,bj)
#else
     &       + m_eta(i,j,bi,bj)
#endif
#else
     &       + EtaN(i,j,bi,bj)
#endif
     &       * rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            proc_precip_volume_flux = proc_precip_volume_flux
     &       + precip(i,j,bi,bj) * rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            proc_evap_volume_flux = proc_evap_volume_flux
     &       + evap(i,j,bi,bj) * rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            proc_wet_area = proc_wet_area
     &        + rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
#ifdef MODIDY_POSITIVE_PRECIP_ONLY
            IF(precip(i,j,bi,bj).GT.0. _d 0)THEN
             proc_corr_area = proc_corr_area
     &         + rA(i,j,bi,bj) * maskC(i,j,kSrf,bi,bj)
            ENDIF
#endif
      ENDDO 
      ENDDO
      ENDDO
      ENDDO

      procs_volume_above_zero(mpiMyId+1) = proc_volume_above_zero
      procs_precip_volume_flux(mpiMyId+1) = proc_precip_volume_flux
      procs_evap_volume_flux(mpiMyId+1) = proc_evap_volume_flux
      procs_wet_area(mpiMyId+1) = proc_wet_area
      procs_corr_area(mpiMyId+1) = proc_corr_area


C ---- Here all procs besides 0 send their volume sums to
C      proc 0 in this order:
C      1. volume above 0 (calculated from EtaN)
C      2. precip volume flux
C      3. eval volume flux
C      4. wet_area
C      5. corr_area
      

        IF (mpiMyId .gt. 0) then
          if (debug .eq. 1) then
          WRITE(msgBuf,'(A,I5,A)') 'mpiMyId = ',mpiMyId,
     &     '>0  => Calling MPI_SEND'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

C      PRINT *,'proc_volume_above_zero',proc_volume_above_zero
C      PRINT *,'proc_precip_volume_flux',proc_precip_volume_flux
C      PRINT *,'proc_evap_volume_flux',proc_evap_volume_flux
C      PRINT *,'proc_wet_area',proc_wet_area

          endif
          call MPI_SEND(proc_volume_above_zero,
     &             1, MPI_DOUBLE, 0, 1, MPI_COMM_MODEL, ierror)
          call MPI_SEND(proc_precip_volume_flux,
     &             1, MPI_DOUBLE, 0, 2, MPI_COMM_MODEL, ierror)
          call MPI_SEND(proc_evap_volume_flux,
     &             1, MPI_DOUBLE, 0, 3, MPI_COMM_MODEL, ierror)
          call MPI_SEND(proc_wet_area,
     &             1, MPI_DOUBLE, 0, 4, MPI_COMM_MODEL, ierror)
          call MPI_SEND(proc_corr_area,
     &             1, MPI_DOUBLE, 0, 5, MPI_COMM_MODEL, ierror)
        ENDIF

C ---- Here proc 0 receives the volume sums in this order:
C      1. volume above 0 (calculated from EtaN)
C      2. precip volume flux
C      3. eval volume flux
C      4. wet_area
C      5. corr_area

        IF (mpiMyId .eq. 0) then
          if (debug .eq. 1) then
          WRITE(msgBuf,'(A)')  'mpiMyId=0  => Calling MPI_RECV'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

          DO pid=1, nPx*nPy-1

              call MPI_RECV(vol_tmp, 1, MPI_DOUBLE,
     &               pid, 1 ,MPI_COMM_MODEL, status, ierror)
              call MPI_RECV(precip_tmp, 1, MPI_DOUBLE,
     &               pid, 2 ,MPI_COMM_MODEL, status, ierror)
              call MPI_RECV(evap_tmp, 1, MPI_DOUBLE,
     &               pid, 3 ,MPI_COMM_MODEL, status, ierror)
              call MPI_RECV(area_tmp, 1, MPI_DOUBLE,
     &               pid, 4 ,MPI_COMM_MODEL, status, ierror)
              call MPI_RECV(corr_area_tmp, 1, MPI_DOUBLE,
     &               pid, 5 ,MPI_COMM_MODEL, status, ierror)

C            PRINT *, 'pid',pid
C            PRINT *,  'vol_tmp', vol_tmp
C            PRINT *,  'precip_tmp', precip_tmp
C            PRINT *,  'evap_tmp', evap_tmp
C            PRINT *,  'area_tmp', area_tmp

      procs_volume_above_zero(pid+1) = vol_tmp
      procs_precip_volume_flux(pid+1) = precip_tmp
      procs_evap_volume_flux(pid+1) = evap_tmp
      procs_wet_area(pid+1) = area_tmp
      procs_corr_area(pid+1) = corr_area_tmp

          ENDDO
         ENDIF


      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SEND_FLUXES_TO_PROC_0 DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif


       END

CEOC
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
#endif /* ALLOW_USE_MPI */



#ifdef ALLOW_USE_MPI
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: SUBROUTINE SEND_ADJUSTMENT_TO_ALL_PROCS

C !INTERFACE:
      SUBROUTINE SEND_ADJUSTMENT_TO_ALL_PROCS(myThid)

C     !DESCRIPTION:
C     Sending precip adjustment from process 0 to all other processes
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "SLR_CORR_PARAM.h"
#include "SLR_CORR_FIELDS.h"

C     !LOCAL VARIABLES:
      INTEGER vec_id, myThid
      INTEGER pid, status(MPI_STATUS_SIZE), ierror
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER, PARAMETER :: debug = 0
      INTEGER i, j, bi, bj

      _RL adjustment_tmp
CEOP
C----&------------------------------------------------------------------xxxxxxx|
CBOC 

C     /==================================================================\
C     | Note: MPI_SEND params: data_to_send, send_count, send_type,      |
C     |                        destination_ID, tag, comm, ierror         |
C     \==================================================================/

C     /==================================================================\
C     | Note: MPI_RECV params: received_data, receive_count,             |
C     |                        receive_type, sender_ID, tag, comm,       |
C     |                        status, ierror                            |
C     \==================================================================/

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SEND_ADJUSTMENT_TO_ALL_PROCS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif


C ---- Here proc 0 sends the adjustment to all other procs

        IF (mpiMyId .eq. 0) then
          if (debug .eq. 1) then
          WRITE(msgBuf,'(A,I5,A)') 'mpiMyId = ',mpiMyId,
     &     '>0  => Calling MPI_SEND'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

          DO pid=1, nPx*nPy-1
          call MPI_SEND(slrc_precip_adjustment,
     &             1, MPI_DOUBLE, pid, 1, MPI_COMM_MODEL, ierror)
          ENDDO
          
        ENDIF

C ---- Here all proc besides 0 receive the adjutment info

        IF (mpiMyId .gt. 0) then
          if (debug .eq. 1) then
          WRITE(msgBuf,'(A)')  'mpiMyId=0  => Calling MPI_RECV'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif


          call MPI_RECV(adjustment_tmp, 1, MPI_DOUBLE,
     &               0, 1 ,MPI_COMM_MODEL, status, ierror)

         slrc_precip_adjustment = adjustment_tmp
        ENDIF


      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SEND_ADJUSTMENT_TO_ALL_PROCS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif


       END

CEOC
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
#endif /* ALLOW_USE_MPI */

