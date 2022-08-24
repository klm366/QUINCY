!>
!! @brief Contains the routines for the mix processes
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-08-24)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_mix_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish, message_text, em_warn

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_runoff_transport

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! @brief   Compute multi-scale intertile exchange processes
  !!
  PURE ELEMENTAL SUBROUTINE calc_runoff_transport(            &
    & dtime, source_fract, target_fract, runoff_source,       &
    & wtr_latflow_source, wtr_latflow_target, weq_pond_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime                 !< timestep length
    REAL(wp), INTENT(in)    :: source_fract          !< area fraction of water source tile
    REAL(wp), INTENT(in)    :: target_fract          !< area fraction of water target tile
    REAL(wp), INTENT(inout) :: runoff_source,      & !< surface runoff on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        weq_pond_target       !< surface water storage on target tile

    !> Update target ponds with surface runoff from source cell
    !!   Surface runoff is scaled with the tile fractions and moved into the target tile's
    !!   pond storage
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    IF (target_fract > 1.0e-10_wp) THEN
      weq_pond_target    = weq_pond_target + runoff_source * dtime / rhoh2o * (source_fract / target_fract)
      wtr_latflow_target = wtr_latflow_target + runoff_source * (source_fract / target_fract)
      wtr_latflow_source = wtr_latflow_source - runoff_source
      runoff_source      = 0._wp
    END IF

  END SUBROUTINE calc_runoff_transport

#endif
END MODULE mo_mix_process
