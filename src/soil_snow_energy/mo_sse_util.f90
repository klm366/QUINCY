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
  USE mo_exception, ONLY: finish, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reconstruct_soil_column_states, rescale_soil_column_states

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

  SUBROUTINE rescale_soil_column_states(      &
    & snow_fract, snow_fract_old,             &
    & state_withsnow, state_snowfree          &
    )

    REAL(wp), INTENT(in) ::    &
      & snow_fract,            &      !< Snow fraction during this time step     [-]
      & snow_fract_old                !< Snow fraction during last time step     [-]

    REAL(wp), INTENT(inout) :: &
      & state_withsnow,        &      !< Soil state for snow covered soil column [m]/[K]
      & state_snowfree                !< Soil state for snow free soil column    [m]/[K]

    REAL(wp) ::                &
      & state_withsnow_old,    &      !< Old soil state for snow covered soil column     [m]/[K]
      & state_snowfree_old,    &      !< Old soil state for snow free soil column        [m]/[K]
      & state_withsnow_new,    &      !< Updated soil state for snow covered soil column [m]/[K]
      & state_snowfree_new,    &      !< Updated Soil state for snow free soil column    [m]/[K]
      & state_mean,            &      !< Tile mean soil state                            [m]/[K]
      & delta_snow_fract              !< Change in snow fraction                         [-]

    LOGICAL :: do_scaling             !< T: rescale soil column states

    CHARACTER(len=*), PARAMETER :: routine = modname//':rescale_soil_column_states'
    CHARACTER(len=4096) :: message_text_long

    ! Initialize old states
    state_withsnow_old = state_withsnow
    state_snowfree_old = state_snowfree

    ! Disable scaling if snow fraction almost zero or one
    do_scaling = (snow_fract >= 1.0e-6 .AND. snow_fract <= 1._wp - 1.0e-6)

    ! Calculate the tile mean state and snow fraction change
    state_mean = snow_fract_old * state_withsnow_old + (1._wp - snow_fract_old) * state_snowfree_old
    delta_snow_fract = snow_fract - snow_fract_old

    ! Compute the state change in the increasing fraction (snow free or covered) and deduce
    ! the change in the decreasing fraction assuming a constant mean state
    IF (delta_snow_fract > 0._wp .AND. do_scaling) THEN
      ! Snow covered fraction increasing
      state_withsnow_new = MAX(0._wp, 1._wp / snow_fract * &
        & (state_withsnow_old * snow_fract_old + state_snowfree_old * delta_snow_fract))
      state_snowfree_new = MAX(0._wp, (state_mean - state_withsnow_new * snow_fract) / (1._wp - snow_fract))
    ELSE IF (delta_snow_fract < 0._wp .AND. do_scaling) THEN
      ! Snow free fraction increasing
      state_snowfree_new = MAX(0._wp, 1._wp / (1._wp - snow_fract) * &
        & (state_snowfree_old * (1._wp - snow_fract_old) + state_withsnow_old * ABS(delta_snow_fract)))
      state_withsnow_new = MAX(0._wp, (state_mean - state_snowfree_new * (1._wp - snow_fract)) / snow_fract)
    ELSE
      ! No scaling needed
      do_scaling         = .FALSE.
      state_withsnow_new = state_withsnow_old
      state_snowfree_new = state_snowfree_old
    END IF

    IF (do_scaling) THEN
      ! Check that no physical bounds are violated
      IF (state_snowfree_new < 0._wp .OR. state_withsnow_new < 0._wp) THEN
        ! WRITE (message_text,*) 'New states are outside of physical bounds'
        WRITE (message_text_long,*) 'Outside of physical bounds ', NEW_LINE('a'),   &
          &  'Old snow fraction:  ', snow_fract_old,     NEW_LINE('a'),          &
          &  'New snow fraction:  ', snow_fract,         NEW_LINE('a'),          &
          &  'Old state snowfree: ', state_snowfree_old, NEW_LINE('a'),          &
          &  'Old state withsnow: ', state_withsnow_old, NEW_LINE('a'),          &
          &  'New state snowfree: ', state_snowfree_new, NEW_LINE('a'),          &
          &  'New state withsnow: ', state_withsnow_new, NEW_LINE('a'),          &
          &  'Old mean state:     ', state_mean,         NEW_LINE('a'),          &
          &  'New mean state:     ', snow_fract * state_withsnow_new + (1._wp - snow_fract) * state_snowfree_new
        CALL finish(TRIM(routine), message_text_long)
      END IF

      ! Check that rescaling worked as expected
      IF (ABS(snow_fract * state_withsnow_new + (1._wp - snow_fract) * state_snowfree_new - state_mean) > 1.0E-10) THEN
        ! WRITE (message_text,*) 'State rescaling failed with a difference of ', &
        !   & snow_fract * state_withsnow_new + (1._wp - snow_fract) * state_snowfree_new - state_mean
        WRITE (message_text_long,*) 'State rescaling failed! ', NEW_LINE('a'),   &
          &  'Old snow fraction:  ', snow_fract_old,     NEW_LINE('a'),          &
          &  'New snow fraction:  ', snow_fract,         NEW_LINE('a'),          &
          &  'Old state snowfree: ', state_snowfree_old, NEW_LINE('a'),          &
          &  'Old state withsnow: ', state_withsnow_old, NEW_LINE('a'),          &
          &  'New state snowfree: ', state_snowfree_new, NEW_LINE('a'),          &
          &  'New state withsnow: ', state_withsnow_new, NEW_LINE('a'),          &
          &  'Old mean state:     ', state_mean,         NEW_LINE('a'),          &
          &  'New mean state:     ', snow_fract * state_withsnow_new + (1._wp - snow_fract) * state_snowfree_new
        CALL finish(TRIM(routine), message_text_long)
      END IF
    END IF

    ! Update soil states
    state_withsnow = state_withsnow_new
    state_snowfree = state_snowfree_new

  END SUBROUTINE rescale_soil_column_states

#endif
END MODULE mo_sse_util
