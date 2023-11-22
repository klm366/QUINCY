!>
!! @brief Flow simulation of reservoir cascade
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-04-20)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hd_reservoir_cascade
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,                   ONLY: wp

  USE mo_time_base,              ONLY: rdaylen

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reservoir_cascade

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_reservoir_cascade'

CONTAINS

  !------------------------------------------------------------------------------------------------
  SUBROUTINE reservoir_cascade(inflow, nc, steplen, retention, nres, res_content, outflow)
  !------------------------------------------------------------------------------------------------
  ! routine based on jsbach3 routine kasglob by Stefan Hagemann
  !   cascade of n equal linear reservoirs.
  !   linear reservoir approach: The outflow from a reservoir is proportional to its content.
  !------------------------------------------------------------------------------------------------

    REAL(wp),  INTENT(in)    :: inflow(:)        !< inflow to the cascade [m3/s]
    REAL(wp),  INTENT(in)    :: steplen          !< time step [s]
    INTEGER,   INTENT(in)    :: nc               !< current block size
    REAL(wp),  INTENT(in)    :: retention(:)     !< retention parameter (residence time of water in the reservoir) [day]
    REAL(wp),  INTENT(in)    :: nres(:)          !< number of reservoirs
    REAL(wp),  INTENT(inout) :: res_content(:,:) !< intermediate content of reservoir cascade [m3]
    REAL(wp),  INTENT(out)   :: outflow(:)       !< outflow at the end of the cascade [m3/s]

    ! local variables
    ! ---------------

    REAL(wp) :: &
      & amod,              &
      & akdiv(nc),         &
      & inflow_cascade(nc)   !< inflow into each level of the cascade

    INTEGER :: res, nres_max
    INTEGER :: ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':reservoir_cascade'

    !$ACC DATA CREATE(akdiv, inflow_cascade)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(amod)
    DO ic = 1, nc
      IF (nres(ic) > 0._wp) THEN ! land - without internal drainage cells
        ! unit conversions: [day] -> [seconds]
        amod  = retention(ic) * rdaylen
        akdiv(ic) = 1._wp / (amod + steplen) * steplen
      ELSE                       ! ocean and lakes
        akdiv(ic) = 0._wp
      END IF
      ! conversions: [m3/s] -> [m3]
      inflow_cascade(ic) = inflow(ic) * steplen
    END DO
    !$ACC END PARALLEL LOOP

    nres_max = SIZE(res_content(1,:))    ! maximum number of reservoirs = 'vertical' dimension of reservoir cascade

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    !$ACC LOOP SEQ
    DO res = 1, nres_max
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        IF (res <= NINT(nres(ic))) THEN
          res_content(ic,res) = res_content(ic,res) + inflow_cascade(ic)
          outflow(ic) = res_content(ic,res) * akdiv(ic)
          res_content(ic,res) = res_content(ic,res) - outflow(ic)
        ELSE
          outflow(ic) = inflow_cascade(ic)
        END IF
        inflow_cascade(ic) = outflow(ic)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP

    ! conversions: [m3] -> [m3/s]
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      outflow(ic) = outflow(ic) / steplen
    END DO
    !$ACC END LOOP

    !$ACC END PARALLEL

    !$ACC END DATA
    !$ACC WAIT(1)

  END SUBROUTINE reservoir_cascade

#endif
END MODULE mo_hd_reservoir_cascade
