!> Contains the routines for the disturb processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>

!NEC$ options "-finline-functions"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"

MODULE mo_disturb_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control, ONLY: acc_stream
  USE mo_kind,      ONLY: lp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: burned_fract_jsbach, broken_woody_fract_jsbach

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_process'

CONTAINS

  SUBROUTINE burned_fract_jsbach( &
    & fire_rel_hum_threshold,                    & ! in
    & fire_litter_threshold,                     & ! in
    & fire_minimum,                              & ! in
    & fire_tau,                                  & ! in
    & q_rel_air_climbuf,                         & ! in
    & fuel,                                      & ! in
    & burned_fract                               & ! inout
    & )

    ! Input Arguments
    REAL(lp), INTENT(in) ::     &
      & fire_rel_hum_threshold, &  ! Maximum relative humidity for fire
      & fire_litter_threshold,  &  ! Minimum amount of litter necessary for fire [mol(C)/m^2(grid box)]
      & fire_minimum,           &  ! Minimum fraction woody PFTs burning every year
      & fire_tau,               &  ! Return period of fire at 0% relative humidity [year]
      & q_rel_air_climbuf(:),   &  ! Relative humidity (smoothed) [%]
      & fuel(:)                    ! Amount of fuel [mol(C)/m^2(grid box)]

    ! Output Arguments
    REAL(lp), INTENT(inout) :: burned_fract(:)

    ! Locals
    INTEGER  :: ic, nc

    REAL(lp), PARAMETER :: delta_time_yr = 1._lp / 365._lp

    nc = SIZE(fuel)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ! Minimum fraction burnt every day (needed with natural LCC)
      burned_fract(ic) = fire_minimum * delta_time_yr

      ! If conditions for fire are given ...
      IF (q_rel_air_climbuf(ic) < fire_rel_hum_threshold .AND. fuel(ic) > fire_litter_threshold) THEN
        burned_fract(ic) = burned_fract(ic)                                                                        &
          &               + (fire_rel_hum_threshold - q_rel_air_climbuf(ic)) / (fire_tau * fire_rel_hum_threshold) &
          &                 * delta_time_yr
      END IF

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE burned_fract_jsbach

  ! Calculate windbreak for woody types
  SUBROUTINE broken_woody_fract_jsbach( &
    & wind_threshold,                                  & ! in
    & wind_damage_scale,                               & ! in
    & cover_fract_pot,                                 & ! in
    & max_wind_yDay,                                   & ! in
    & max_wind_clim,                                   & ! in
    & damaged_fract                                    & ! inout
    & )

    ! Arguments
    REAL(lp), INTENT(in) ::       &
      & wind_threshold,           &
      & wind_damage_scale,        &
      & cover_fract_pot(:),       &
      & max_wind_yDay(:),         &
      & max_wind_clim(:)

    REAL(lp), INTENT(inout) :: &
      & damaged_fract(:)

    ! Locals
    INTEGER  :: ic, nc

    REAL(lp), PARAMETER :: delta_time_yr = 1._lp / 365._lp

    nc = SIZE(max_wind_clim)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      IF (cover_fract_pot(ic) > EPSILON(1._lp) .AND. max_wind_yDay(ic) > max_wind_clim(ic) * wind_threshold) THEN
        damaged_fract(ic) = MIN( 1._lp - EPSILON(1._lp) / cover_fract_pot(ic),                                          &
        &                        delta_time_yr * wind_damage_scale * max_wind_yDay(ic) ** 3._lp / max_wind_clim(ic) )
      ENDIF
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE broken_woody_fract_jsbach

#endif
END MODULE mo_disturb_process
