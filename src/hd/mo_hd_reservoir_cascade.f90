!>
!! @brief Flow simulation of reservoir cascade
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-04-20)
!!
MODULE mo_hd_reservoir_cascade
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text

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

    REAL(wp),  INTENT(in)    :: inflow(:)        !< inflow to the cascade
    REAL(wp),  INTENT(in)    :: steplen          !< time steps per day
    INTEGER,   INTENT(in)    :: nc               !< current block size
    REAL(wp),  INTENT(in)    :: retention(:)     !< retention parameter (residence time of water in the reservoir)
    REAL(wp),  INTENT(in)    :: nres(:)          !< number of reservoirs
    REAL(wp),  INTENT(inout) :: res_content(:,:) !< intermediate content of reservoir cascade
    REAL(wp),  INTENT(out)   :: outflow(:)       !< outflow at the end of the cascade

    ! local variables
    ! ---------------

    REAL(wp), ALLOCATABLE    :: amod(:)       !vg wuerde ich gern beides abschaffen
    REAL(wp), ALLOCATABLE    :: akdiv(:)
    REAL(wp), ALLOCATABLE    :: inflow_cascade(:) !< inflow into each level of the cascade

    INTEGER ::  res, nres_max


    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reservoir_cascade'

    ALLOCATE(amod(nc)) 
    ALLOCATE(akdiv(nc))
    ALLOCATE(inflow_cascade(nc))

    amod(:) = retention(:) * nres(:) / AINT(nres(:))    !! vg: lieber im inifile anpassen!!
    akdiv = 1._wp / (amod(:) + steplen) * steplen       !! vg: soll das so??

    inflow_cascade(:) = inflow(:)

    nres_max = SIZE(res_content(1,:))    ! maximum number of reservoirs = 'vertical' dimension of reservoir cascade
    DO res = 1, nres_max
      WHERE (res <= NINT(nres(:)))
        res_content(:,res) = res_content(:,res) + inflow_cascade(:)
        outflow(:) = res_content(:,res) * akdiv
        res_content(:,res) = res_content(:,res) - outflow(:)
      ELSEWHERE
        outflow(:) = inflow_cascade(:)
      ENDWHERE
      inflow_cascade(:) = outflow(:)
    END DO

    DEALLOCATE(amod, akdiv)

  END SUBROUTINE reservoir_cascade

#endif
END MODULE mo_hd_reservoir_cascade
