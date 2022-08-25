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

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_runoff_transport, calc_drainage_transport

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

    REAL(wp), INTENT(in)    :: dtime,              & !< timestep length
      &                        source_fract,       & !< area fraction of water source tile
      &                        target_fract          !< area fraction of water target tile
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

  ! ===============================================================================================================================
  !>
  !! @brief   Compute multi-scale intertile exchange processes
  !!
  PURE ELEMENTAL SUBROUTINE calc_drainage_transport(            &
    & dtime, source_fract, target_fract, field_cap_target,      &
    & latdrain_source, drainage_source, wtr_latflow_source,     &
    & wtr_latflow_target, wtr_soil_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,              & !< timestep length
      &                        source_fract,       & !< area fraction of water source tile
      &                        target_fract,       & !< area fraction of water target tile
      &                        field_cap_target      !< soil field capacity of target cell
    REAL(wp), INTENT(inout) :: latdrain_source,    & !< lateral drainage on source tile
      &                        drainage_source,    & !< bottom drainage on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        wtr_soil_target       !< soil water storage on target tile

    ! Local variables
    REAL(wp)                :: latflow_act           !< Actual lateral flow during time step
                                                     !! considering the target storage limits

    !> Update target soil with lateral drainage from source cell
    !!   Lateral drainage is scaled with the tile fractions and moved into the target tile's
    !!   soil storage. It flows into the same layer as it originates from on the source tile.
    !!   Excess drainage is added to the bottom drainage of the source tile.
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    IF (target_fract > 1.0e-10_wp) THEN
      ! Scale and add all lateral drainage to target tile soil layer
      latflow_act        = latdrain_source * (source_fract / target_fract)
      wtr_soil_target    = wtr_soil_target + latflow_act * dtime / rhoh2o
      ! Reduce lateral flow by upper soil water limit and correct for any excess soil moisture
      latflow_act        = latflow_act - MAX(wtr_soil_target - field_cap_target, 0._wp) * rhoh2o / dtime
      wtr_soil_target    = MIN(wtr_soil_target, field_cap_target)
      wtr_latflow_target = wtr_latflow_target + latflow_act
      ! Rescale the actually transported water back to the source tile
      latflow_act        = latflow_act * (target_fract / source_fract)
      wtr_latflow_source = wtr_latflow_source - latflow_act
      ! Add any residual lateral flow to the bottom drainage of the source tile
      drainage_source    = drainage_source + MAX(latdrain_source - latflow_act, 0._wp)
      latdrain_source    = 0._wp
    END IF

  END SUBROUTINE calc_drainage_transport

#endif
END MODULE mo_mix_process
