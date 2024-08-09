!> Contains the utilities for the soil and snow energy process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_sse_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  !USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reconstruct_soil_column_states

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_util'

CONTAINS

  SUBROUTINE reconstruct_soil_column_states(                                                   &
    & delta_soilwtr_sl, delta_soilice_sl, wtr_soil_sl, ice_soil_sl, weq_wsat_soil_sl, snow_fract_old, &
    wtr_soil_snowfree_sl, ice_soil_snowfree_sl, wtr_soil_withsnow_sl, ice_soil_withsnow_sl     &
    )

    REAL(wp), INTENT(in) ::   &
      & delta_soilwtr_sl,     &      !< Water content difference between snow covered and free soil column [m]
      & delta_soilice_sl,     &      !< Ice content difference between snow covered and free soil column   [m]
      & wtr_soil_sl,          &      !< Tile mean soil water content                                       [m]
      & ice_soil_sl,          &      !< Tile mean soil ice content                                         [m]
      & weq_wsat_soil_sl,     &      !< Maximum soil moisture capacity                                     [m]
      & snow_fract_old               !< Snow fraction during last time step                                [-]

    REAL(wp), INTENT(out) ::  &
      & wtr_soil_snowfree_sl, &      !< Water content of snow free soil column    [m]
      & ice_soil_snowfree_sl, &      !< Ice content of snow free soil column      [m]
      & wtr_soil_withsnow_sl, &      !< Water content of snow covered soil column [m]
      & ice_soil_withsnow_sl         !< Ice content of snow covered soil column   [m]

    REAL(wp) ::               &
      & weq_soil_snowfree,    &      !< Total soil moisture for snow free soil column     [m]
      & weq_soil_withsnow,    &      !< Total soil moisture for snow covered soil column  [m]
      & excess,               &      !< Amount of soil moisture exceeding its upper limit [m]
      & scale_wtr,            &      !< Scaling factor for soil water                     [-]
      & scale_ice                    !< Scaling factor for soil ice                       [-]

    CHARACTER(len=*), PARAMETER :: routine = modname//':reconstruct_soil_column_states'

    IF ((ABS(delta_soilwtr_sl) > 1.0e-6_wp .OR. ABS(delta_soilice_sl) > 1.0e-6_wp) .AND. &
      & snow_fract_old > 1.0e-6_wp) THEN
      ! Initialize water and ice content for both soil columns, from last time steps difference
      ! and current snow fraction
      wtr_soil_snowfree_sl = wtr_soil_sl - snow_fract_old * delta_soilwtr_sl
      ice_soil_snowfree_sl = ice_soil_sl - snow_fract_old * delta_soilice_sl
      wtr_soil_withsnow_sl = wtr_soil_snowfree_sl + delta_soilwtr_sl
      ice_soil_withsnow_sl = ice_soil_snowfree_sl + delta_soilice_sl

      ! Check whether soil water/ice state is still within bounds and adjust it if necessary.
      ! This is done scaling the water and ice differences such that no soil bounds are violated.
      weq_soil_snowfree = MAX(1.0e-10_wp, wtr_soil_snowfree_sl) + MAX(1.0e-10_wp, ice_soil_snowfree_sl)
      weq_soil_withsnow = MAX(1.0e-10_wp, wtr_soil_withsnow_sl) + MAX(1.0e-10_wp, ice_soil_withsnow_sl)

      IF (wtr_soil_snowfree_sl < 0._wp .OR. wtr_soil_withsnow_sl < 0._wp .OR. &
        & ice_soil_snowfree_sl < 0._wp .OR. ice_soil_withsnow_sl < 0._wp .OR. &
        & weq_soil_snowfree > weq_wsat_soil_sl + 1.0e+10_wp .OR.              &
        & weq_soil_withsnow > weq_wsat_soil_sl + 1.0e+10_wp) THEN

        ! Calculate scale factors for the differences in soil water and soil ice between the fractions
        scale_wtr = 1._wp
        scale_ice = 1._wp

        IF (wtr_soil_sl < 1.0e-10_wp) THEN
          ! Set scaling to zero in case the average soil water is already close to zero
          scale_wtr = 0._wp
        ELSE
          ! Ensure scaled soil water doesn't fall below zero
          IF (wtr_soil_snowfree_sl < 1.0e-10_wp) THEN
              scale_wtr = min(scale_wtr, &
                & 1.0_wp - ABS(wtr_soil_snowfree_sl - 1.0e-10_wp) / (wtr_soil_sl - wtr_soil_snowfree_sl))
          END IF
          IF (wtr_soil_withsnow_sl < 1.0e-10_wp) THEN
              scale_wtr = min(scale_wtr, &
                & 1.0_wp - ABS(wtr_soil_withsnow_sl - 1.0e-10_wp) / (wtr_soil_sl - wtr_soil_withsnow_sl))
          END IF
        END IF

        IF (ice_soil_sl < 1.0e-10_wp) THEN
          ! Set scaling to zero in case the average soil ice is already close to zero
          scale_ice = 0._wp
        ELSE
          ! Ensure water and ice components don't drop below zero
          IF (ice_soil_snowfree_sl < 1.0e-10_wp) THEN
              scale_ice = min(scale_ice, &
                & 1.0_wp - ABS(ice_soil_snowfree_sl - 1.0e-10_wp) / (ice_soil_sl - ice_soil_snowfree_sl))
          END IF
          IF (ice_soil_withsnow_sl < 1.0e-10_wp) THEN
              scale_ice = min(scale_ice, &
                & 1.0_wp - ABS(ice_soil_withsnow_sl - 1.0e-10_wp) / (ice_soil_sl - ice_soil_withsnow_sl))
          END IF
        END IF

        ! Ensure total soil column moisture doesn't exceed maximum capacity
        IF (weq_soil_snowfree > weq_wsat_soil_sl + 1.0e-10_wp) THEN
            excess = weq_soil_snowfree - weq_wsat_soil_sl
            scale_wtr = min(scale_wtr, max(0._wp, &
              & 1._wp - excess / (weq_soil_snowfree - (wtr_soil_sl + ice_soil_sl))))
            scale_ice = min(scale_ice, max(0._wp, &
              & 1._wp - excess / (weq_soil_snowfree) - (wtr_soil_sl + ice_soil_sl)))
        END IF
        IF (weq_soil_withsnow > weq_wsat_soil_sl + 1.0e-10_wp) THEN
            excess = weq_soil_withsnow - weq_wsat_soil_sl
            scale_wtr = min(scale_wtr, max(0._wp, &
              & 1._wp - excess / (weq_soil_withsnow - (wtr_soil_sl + ice_soil_sl))))
            scale_ice = min(scale_ice, max(0._wp, &
              & 1._wp - excess / (weq_soil_withsnow - (wtr_soil_sl + ice_soil_sl))))
        END IF

        ! Recalculate with bounds-adjusted differences
        wtr_soil_snowfree_sl = wtr_soil_sl - snow_fract_old * delta_soilwtr_sl * scale_wtr
        ice_soil_snowfree_sl = ice_soil_sl - snow_fract_old * delta_soilice_sl * scale_ice
        wtr_soil_withsnow_sl = wtr_soil_snowfree_sl + delta_soilwtr_sl * scale_wtr
        ice_soil_withsnow_sl = ice_soil_snowfree_sl + delta_soilice_sl * scale_ice
      END IF
    ELSE
      ! Initialize water and ice content for both soil columns equally
      ! if snow fraction or water/ice differences are very low
      wtr_soil_snowfree_sl = wtr_soil_sl
      ice_soil_snowfree_sl = ice_soil_sl
      wtr_soil_withsnow_sl = wtr_soil_sl
      ice_soil_withsnow_sl = ice_soil_sl
    END IF

  END SUBROUTINE reconstruct_soil_column_states

#endif
END MODULE mo_sse_util
