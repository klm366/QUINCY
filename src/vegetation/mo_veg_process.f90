!>
!! @brief Contains the routines for the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: canopy_conductance, water_stress_factor, foliage_projected_cover, canopy_albedo

  INTERFACE canopy_conductance
    MODULE PROCEDURE canopy_conductance_unstressed_echam5
    MODULE PROCEDURE canopy_conductance_stressed
  END INTERFACE canopy_conductance

  INTERFACE canopy_albedo
    MODULE PROCEDURE canopy_albedo_echam5
  END INTERFACE canopy_albedo

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_process'

CONTAINS

  REAL(wp) ELEMENTAL FUNCTION canopy_conductance_unstressed_echam5(lai, par) RESULT(conductance)

    USE mo_veg_constants, ONLY: k => conductance_k, a => conductance_a, b => conductance_b, c => conductance_c

    REAL(wp), INTENT(in) :: &
      & lai, &
      & par

    ! Local variables
    !
!!$    REAL(wp) :: d(SIZE(par))
!!$    REAL(wp) :: zpar(SIZE(par))
    REAL(wp) :: d, zpar

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':canopy_conductance_unstressed_echam5'

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    zpar = MAX(1.E-10_wp, par)

!!$    WHERE (lai > EPSILON(1._wp))
    IF (lai > EPSILON(1._wp)) THEN
      d = (a + b*c) / (c * zpar)
      conductance = ( LOG((d * EXP(k*lai) + 1._wp) / (d + 1._wp)) * b / (d * zpar) - &
                    & LOG((d + EXP(-k*lai)) / (d + 1._wp))                           &
                    ) / (k * c)
    ELSE
      conductance = EPSILON(1._wp)
    END IF

  END FUNCTION canopy_conductance_unstressed_echam5

  REAL(wp) ELEMENTAL FUNCTION canopy_conductance_stressed(cond_unstressed, water_stress, air_is_saturated) RESULT(conductance)

    REAL(wp), INTENT(in) :: &
      & cond_unstressed, &
      & water_stress
    LOGICAL,  INTENT(in) :: air_is_saturated

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':canopy_conductance_stressed'

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    IF (air_is_saturated) THEN
      conductance = EPSILON(1._wp)
    ELSE
      conductance = cond_unstressed * water_stress
    END IF

  END FUNCTION canopy_conductance_stressed

  REAL(wp) ELEMENTAL FUNCTION water_stress_factor (ws, ws_fc, ws_pwp)

    ! TBD: Maybe change later so that it uses geographically dependend critical values (like wilting point)

    USE mo_soil_constants, ONLY: ws_crit_fract

    REAL(wp), INTENT(in) :: ws      !< Soil water content
    REAL(wp), INTENT(in) :: ws_fc   !< Soil water content at field capacity
    REAL(wp), INTENT(in) :: ws_pwp  !< Soil water content at permanent wilting point

    REAL(wp) :: ws_crit             !< Soil water content at critical value

    ws_crit = ws_fc * ws_crit_fract

    IF (ws_crit - ws_pwp < EPSILON(1._wp)) THEN  ! This should only happen if there's no water in soil (ws=ws_fc=ws_pwp=0)
      water_stress_factor = 0._wp
    ELSE
      water_stress_factor = MAX (0._wp, &
        MIN (1._wp, (ws - ws_pwp) / (ws_crit - ws_pwp) ))
    END IF

  END FUNCTION water_stress_factor

  REAL(wp) ELEMENTAL FUNCTION foliage_projected_cover(fpc_max, lai, clumpiness)

    REAL(wp), INTENT(in) :: &
      & fpc_max,            &
      & lai,                &
      & clumpiness

    foliage_projected_cover = fpc_max * (1._wp - EXP(-lai/clumpiness))

  END FUNCTION foliage_projected_cover

  ! TBD: At the moment, this computes the albedo for canopy and ground. Should be split between the veg and soil processes
  REAL(wp) ELEMENTAL FUNCTION canopy_albedo_echam5( surface_temp, forest_fract, grnd_snow_fract, canopy_snow_fract, lai, &
                                                  &  albedo_snow_min, albedo_snow_max, albedo_background)

    USE mo_physical_constants, ONLY: tmelt
    USE mo_veg_constants,      ONLY: &
      & albedo_snow_can,             & !> Albedo of snow on canopy
      & sky_view_factor                !> Sky view factor

    REAL(wp), INTENT(in) :: &
      surface_temp,         &
      forest_fract,         &
      grnd_snow_fract,      &
      canopy_snow_fract,    &
      lai,                  &
      albedo_snow_min,      &
      albedo_snow_max,      &
      albedo_background

    REAL(wp) ::               &
      & sky_view_fract,       &
      & forest_view_fract,    &
      & snow_albedo_grnd,     &  !< Albedo of snow on ground
      & min_temp_snow_albedo, &  !< Temperature below which albedo of snow is at maximum
      & albedo_grnd,          &  !< Albedo on ground
      & albedo_can               !< Albedo on canopy

    min_temp_snow_albedo = tmelt - 5._wp

    IF (surface_temp >= tmelt) THEN
      snow_albedo_grnd = albedo_snow_min
    ELSE IF (surface_temp < min_temp_snow_albedo) THEN
      snow_albedo_grnd = albedo_snow_max
    ELSE
      snow_albedo_grnd = albedo_snow_min &
                       & + (albedo_snow_max - albedo_snow_min) * (tmelt - surface_temp) / (tmelt - min_temp_snow_albedo)
    END IF

    albedo_grnd = grnd_snow_fract * snow_albedo_grnd + (1._wp - grnd_snow_fract) * albedo_background

    sky_view_fract = EXP(-sky_view_factor * MAX(lai, 2._wp))  ! TBD: the value 2.0 should be replaced by stem_area
    forest_view_fract = forest_fract * (1._wp - sky_view_fract)

    albedo_can = canopy_snow_fract  * albedo_snow_can + (1._wp - canopy_snow_fract) * albedo_background

    canopy_albedo_echam5 = MAX(albedo_background, &
                         &     (1._wp - forest_view_fract) * albedo_grnd + forest_view_fract * albedo_can)

  END FUNCTION canopy_albedo_echam5

#endif
END MODULE mo_veg_process
