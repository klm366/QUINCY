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
  USE mo_exception, ONLY: finish, message, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reconstruct_soil_column_states, rescale_soil_column_states

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_util'

CONTAINS

  SUBROUTINE reconstruct_soil_column_states(                                               &
    & ice_frac_snowfree_sl, wtr_soil_sl, ice_soil_sl, weq_wsat_soil_sl, snow_fract_old,    &
    wtr_soil_snowfree_sl, ice_soil_snowfree_sl, wtr_soil_withsnow_sl, ice_soil_withsnow_sl &
    )

    REAL(wp), INTENT(in) ::   &
      & ice_frac_snowfree_sl, &      !< Fraction of soil ice in snow free soil column [-]
      & wtr_soil_sl,          &      !< Tile mean soil water content                  [m]
      & ice_soil_sl,          &      !< Tile mean soil ice content                    [m]
      & weq_wsat_soil_sl,     &      !< Maximum soil moisture capacity                [m]
      & snow_fract_old               !< Snow fraction during last time step           [-]

    REAL(wp), INTENT(out) ::  &
      & wtr_soil_snowfree_sl, &      !< Water content of snow free soil column    [m]
      & ice_soil_snowfree_sl, &      !< Ice content of snow free soil column      [m]
      & wtr_soil_withsnow_sl, &      !< Water content of snow covered soil column [m]
      & ice_soil_withsnow_sl         !< Ice content of snow covered soil column   [m]

    REAL(wp) ::               &
      & weq_soil_sl,          &      !< Total soil moisture for snow free soil column     [m]
      & column_mean                  !< New soil column mean state                        [m]

    CHARACTER(len=*), PARAMETER :: routine = modname//':reconstruct_soil_column_states'

    ! Re-construct the soil ice and water fractions for the individual columns. While the
    ! mean state will have changed during the soil hydrology computation, the relation
    ! between soil ice in the two columns should have remained similar.
    ! The reconstruction assumes that the total soil moisture content is identical in both
    ! columns and only the respective fraction of water and ice differ.

    weq_soil_sl = wtr_soil_sl + ice_soil_sl

    IF (snow_fract_old > 1.0E-10_wp) THEN
      ice_soil_snowfree_sl = MAX(0._wp, MIN(weq_soil_sl, ice_frac_snowfree_sl * ice_soil_sl))
      ice_soil_withsnow_sl = (ice_soil_sl - ice_soil_snowfree_sl * (1._wp - snow_fract_old)) / snow_fract_old
      wtr_soil_snowfree_sl = MAX(0._wp, weq_soil_sl - ice_soil_snowfree_sl)
      wtr_soil_withsnow_sl = MAX(0._wp, weq_soil_sl - ice_soil_withsnow_sl)
    ELSE
      ! Initialize water and ice content for both soil columns equally
      ! if snow fraction is very low
      ice_soil_snowfree_sl = ice_soil_sl
      ice_soil_withsnow_sl = ice_soil_sl
      wtr_soil_snowfree_sl = wtr_soil_sl
      wtr_soil_withsnow_sl = wtr_soil_sl
    END IF

    ! Sanity checks for ...
    ! ... water conservation
    column_mean = wtr_soil_withsnow_sl * snow_fract_old + wtr_soil_snowfree_sl * (1._wp - snow_fract_old)
    IF (ABS(wtr_soil_sl - column_mean) > 1.0E-13_wp) THEN
      WRITE (message_text,*) 'Soil water content not conserved by ', wtr_soil_sl - column_mean, ' m.', NEW_LINE('a'), &
        & 'Snow fraction is ', snow_fract_old
      CALL finish(TRIM(routine), message_text)
    END IF

    ! ... ice conservation
    column_mean = ice_soil_withsnow_sl * snow_fract_old + ice_soil_snowfree_sl * (1._wp - snow_fract_old)
    IF (ABS(ice_soil_sl - column_mean) > EPSILON(1.0_wp)) THEN
      WRITE (message_text,*) 'Soil ice content not conserved by ', ice_soil_sl - column_mean, ' m.'
      CALL finish(TRIM(routine), message_text)
    END IF

    ! Sanity check for upper soil moisture constrains
    IF (wtr_soil_withsnow_sl + ice_soil_withsnow_sl > weq_wsat_soil_sl + EPSILON(1._wp) .OR. &
      & wtr_soil_snowfree_sl + ice_soil_snowfree_sl > weq_wsat_soil_sl + EPSILON(1._wp) .OR. &
      & wtr_soil_withsnow_sl + ice_soil_withsnow_sl < 0._wp                             .OR. &
      & wtr_soil_snowfree_sl + ice_soil_snowfree_sl < 0._wp) THEN
      WRITE (message_text,*) 'Soil column moisture exceeds bounds:', NEW_LINE('a'),   &
        & 'snow:          ', snow_fract_old, NEW_LINE('a'), &
        & 'ice_fract:     ', ice_frac_snowfree_sl, NEW_LINE('a'), &
        & 'withsnow:      ', wtr_soil_withsnow_sl, ice_soil_withsnow_sl, NEW_LINE('a'), &
        & 'snowfree:      ', wtr_soil_snowfree_sl, ice_soil_snowfree_sl, NEW_LINE('a'), &
        & 'both:          ', wtr_soil_sl, ice_soil_sl, weq_wsat_soil_sl
      CALL finish(TRIM(routine), message_text)
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

    LOGICAL :: partial_snow_cover     !< T: partial snow cover exists

    CHARACTER(len=*), PARAMETER :: routine = modname//':rescale_soil_column_states'
    CHARACTER(len=4096) :: message_text_long

    ! Initialize old states
    state_withsnow_old = state_withsnow
    state_snowfree_old = state_snowfree

    ! Calculate the tile mean state and snow fraction change
    state_mean = snow_fract_old * state_withsnow_old + (1._wp - snow_fract_old) * state_snowfree_old
    delta_snow_fract = snow_fract - snow_fract_old

    ! Disable scaling if actual snow fraction is almost zero or one
    partial_snow_cover = (snow_fract >= 1.0e-10_wp .AND. snow_fract <= 1._wp - 1.0e-10_wp)

    ! Compute the state change in the increasing fraction (snow free or covered) and deduce
    ! the change in the decreasing fraction assuming a constant mean state
    IF (delta_snow_fract > 0._wp .AND. partial_snow_cover) THEN
      ! Snow covered fraction increasing
      state_withsnow_new = MAX(0._wp, 1._wp / snow_fract * &
        & (state_withsnow_old * snow_fract_old + state_snowfree_old * delta_snow_fract))
      state_snowfree_new = MAX(0._wp, (state_mean - state_withsnow_new * snow_fract) / (1._wp - snow_fract))
    ELSE IF (delta_snow_fract < 0._wp .AND. partial_snow_cover) THEN
      ! Snow free fraction increasing
      state_snowfree_new = MAX(0._wp, 1._wp / (1._wp - snow_fract) * &
        & (state_snowfree_old * (1._wp - snow_fract_old) + state_withsnow_old * ABS(delta_snow_fract)))
      state_withsnow_new = MAX(0._wp, (state_mean - state_snowfree_new * (1._wp - snow_fract)) / snow_fract)
    ELSE IF (partial_snow_cover) THEN
      ! No scaling needed due to constant (partial) snow fraction
      state_withsnow_new = state_withsnow_old
      state_snowfree_new = state_snowfree_old
    ELSE
      ! Only one column effectively exists because snow fraction is either zero or one, use mean state
      state_withsnow_new = state_mean
      state_snowfree_new = state_mean
    END IF

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
    IF (ABS(snow_fract * state_withsnow_new + (1._wp - snow_fract) * state_snowfree_new - state_mean) > 1.0E-13_wp) THEN
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

    ! Update soil states
    state_withsnow = state_withsnow_new
    state_snowfree = state_snowfree_new

  END SUBROUTINE rescale_soil_column_states

#endif
END MODULE mo_sse_util
