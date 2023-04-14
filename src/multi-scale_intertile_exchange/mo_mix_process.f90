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

  PUBLIC :: calc_runoff_transport, calc_drainage_transport, calc_intertile_diffusion

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of surface runoff
  !!
  PURE ELEMENTAL SUBROUTINE calc_runoff_transport(             &
    & l_lag, dtime, retention_time, wtr_pond_max_source,       &
    & source_fract, target_fract, runoff_source,               &
    & wtr_latflow_source, wtr_latflow_target, wtr_pond_source, &
    & wtr_pond_target, runoff_retention_source                 &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    LOGICAL,  INTENT(in)    :: l_lag                      !< T: apply lag to overland flow
    REAL(wp), INTENT(in)    :: dtime,                   & !< timestep length                        [s]
      &                        retention_time,          & !< retention time for lateral transport   [s]
      &                        wtr_pond_max_source,     & !< maximum surface water storage on source tile [m]
      &                        source_fract,            & !< area fraction of water source tile
      &                        target_fract               !< area fraction of water target tile
    REAL(wp), INTENT(inout) :: runoff_source,           & !< surface runoff on source tile          [kg m-2 s-1]
      &                        wtr_latflow_source,      & !< lateral water flow from source tile    [kg m-2 s-1]
      &                        wtr_latflow_target,      & !< lateral water flow towards target tile [kg m-2 s-1]
      &                        wtr_pond_source,         & !< surface water storage on source tile   [m]
      &                        wtr_pond_target,         & !< surface water storage on target tile   [m]
      &                        runoff_retention_source    !< intermediate storage for overland flow [m]

    REAL(wp)                :: overland_flx_source,     & !< surface flow leaving source tile       [m]
      &                        overland_flx_target        !< surface flow arriving at target tile   [m]

    !> Update target ponds with surface runoff/pond overflow from source tile
    !!   Surface runoff is scaled with the tile fractions and moved into the target tile's
    !!   pond storage.
    !!   In case of long transport routes (upland to lowland) a retention storage is used.
    !! @todo the retention store lag time needs to be physics based, not fixed to some arbitrary value

    IF (target_fract > 1.0e-10_wp) THEN
      ! Compute pond overflow on source tile and update source tile pond storage
      overland_flx_source = MAX(0._wp, wtr_pond_source + runoff_source * dtime / rhoh2o - wtr_pond_max_source)
      wtr_pond_source     = MIN(wtr_pond_max_source, wtr_pond_source + runoff_source * dtime / rhoh2o)
      ! Set source tile runoff to zero as it goes down the hill
      runoff_source       = 0._wp

      ! Compute flows for and update runoff retention storage
      IF (l_lag) THEN
        runoff_retention_source = runoff_retention_source + overland_flx_source
        overland_flx_source     = runoff_retention_source * MIN(1.0_wp, dtime / retention_time)
        runoff_retention_source = MAX(0._wp, runoff_retention_source - overland_flx_source)
      END IF

      ! Scale overland flow to account for tile fraction differences
      overland_flx_target = overland_flx_source * (source_fract / target_fract)
      ! Update target pond storage
      wtr_pond_target     = wtr_pond_target + overland_flx_target
      ! Update lateral flows
      wtr_latflow_source  = wtr_latflow_source - overland_flx_source * rhoh2o / dtime
      wtr_latflow_target  = wtr_latflow_target + overland_flx_target * rhoh2o / dtime

    END IF

  END SUBROUTINE calc_runoff_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of lateral drainage
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
      &                        field_cap_target      !< soil field capacity of target tile
    REAL(wp), INTENT(inout) :: latdrain_source,    & !< lateral drainage on source tile
      &                        drainage_source,    & !< bottom drainage on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        wtr_soil_target       !< soil water storage on target tile

    ! Local variables
    REAL(wp)                :: latflow_act           !< Actual lateral flow during time step
                                                     !! considering the target storage limits

    !> Update target soil with lateral drainage from source tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Lateral drainage is scaled with the tile fractions and moved into the target tile's
    !!   soil storage. It flows into the same layer as it originates from on the source tile.
    !!   Excess drainage is added to the total drainage of the source tile.
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    IF (target_fract > 1.0e-10_wp) THEN
      ! Substract diagnosed lateral drainage from total drainage
      drainage_source    = drainage_source - latdrain_source
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
      ! Add any residual lateral flow to the total drainage of the source tile
      drainage_source    = drainage_source + MAX(latdrain_source - latflow_act, 0._wp)
      latdrain_source    = 0._wp
    END IF

  END SUBROUTINE calc_drainage_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles
  !!
  PURE ELEMENTAL SUBROUTINE calc_intertile_diffusion(           &
    & dtime, source_fract, target_fract, field_cap_source,      &
    & field_cap_target, wtr_latflow_source, wtr_latflow_target, &
    & wtr_soil_source, wtr_soil_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,               & !< timestep length
      &                        source_fract,        & !< area fraction of water source tile
      &                        target_fract,        & !< area fraction of water target tile
      &                        field_cap_source,    & !< soil field capacity of source tile
      &                        field_cap_target       !< soil field capacity of target tile
    REAL(wp), INTENT(inout) :: wtr_latflow_source,  & !< lateral water flow from source tile
      &                        wtr_latflow_target,  & !< lateral water flow towards target tile
      &                        wtr_soil_source,     & !< soil water storage on source tile
      &                        wtr_soil_target        !< soil water storage on target tile

    ! Local variables
    REAL(wp)                :: diffus_source,       & !< Diffusion flux at source tile per time step
      &                        diffus_target,       & !< Diffusion flux at target tile per time step
      &                        wtr_soil_rel_source, & !< Relative soil water content in source tile
      &                        wtr_soil_rel_target, & !< Relative soil water content in target tile
      &                        wtr_soil_rel_mean      !< Mean relative soil water content of both tiles

    !> Compute lateral diffusion flux between soil water storages of source and target tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Lateral diffusion is implemented such that is compensates the the differences
    !!   in relative soil water content between both tiles. Thus, the flow is not necessarily
    !!   from the source to the target tile but follows the moisture gradient.
    !! @todo This dummy implementation sets the lateral diffusion flux to 10% of what would
    !!   be needed to fully compensate the soil moisture difference. Later on, soil hydraulic
    !!   conductivity should be used.

    IF (source_fract > 1.0e-10_wp .AND. target_fract > 1.0e-10_wp .AND. &
        field_cap_source > 1.0e-10_wp .AND. field_cap_target > 1.0e-10_wp) THEN
      ! Compute relative soil water content as well as the fraction-weighted relative soil water content mean between both tiles
      wtr_soil_rel_source = wtr_soil_source / field_cap_source
      wtr_soil_rel_target = wtr_soil_target / field_cap_target
      wtr_soil_rel_mean  = 1.0 / (source_fract + target_fract) * &
        & (wtr_soil_rel_source * source_fract + wtr_soil_rel_target * target_fract)
      ! Compute diffusive flow between both tiles
      diffus_target = (wtr_soil_rel_mean - wtr_soil_rel_target) * field_cap_target * 0.1
      diffus_source = (wtr_soil_rel_mean - wtr_soil_rel_source) * field_cap_source * 0.1
      ! Add diffusive flow to soil storages
      wtr_soil_target = wtr_soil_target + diffus_target
      wtr_soil_source = wtr_soil_source + diffus_source
      ! Add diffusive flow to overal lateral flow
      wtr_latflow_target = wtr_latflow_target + diffus_target * rhoh2o / dtime
      wtr_latflow_source = wtr_latflow_source + diffus_source * rhoh2o / dtime
    END IF

  END SUBROUTINE calc_intertile_diffusion

#endif
END MODULE mo_mix_process
