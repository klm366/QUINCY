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

  PUBLIC :: get_canopy_conductance, get_water_stress_factor, get_foliage_projected_cover

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_process'

CONTAINS

  REAL(wp) ELEMENTAL FUNCTION get_canopy_cond_unstressed_simple(lai, par) RESULT(conductance)

    USE mo_veg_constants, ONLY: k => conductance_k, a => conductance_a, b => conductance_b, c => conductance_c

    REAL(wp), INTENT(in) ::     &
                           lai, &
                           par

    ! Local variables
    !
!!$    REAL(wp) :: d(SIZE(par))
!!$    REAL(wp) :: zpar(SIZE(par))
    REAL(wp) :: d, zpar

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_canopy_cond_unstressed_simple'

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

  END FUNCTION get_canopy_cond_unstressed_simple

  
  REAL(wp) ELEMENTAL FUNCTION get_canopy_cond_stressed_simple(cond_unstressed, water_stress, air_is_saturated) RESULT(conductance)

    REAL(wp), INTENT(in) :: &
      & cond_unstressed, &
      & water_stress
    LOGICAL,  INTENT(in) :: air_is_saturated

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_canopy_cond_stressed_simple'

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    IF (air_is_saturated) THEN
      conductance = EPSILON(1._wp)
    ELSE
      conductance = cond_unstressed * water_stress
    END IF

  END FUNCTION get_canopy_cond_stressed_simple

  REAL(wp) ELEMENTAL FUNCTION get_water_stress_factor (w_soil_sl, w_soil_fc_sl, w_soil_pwp_sl)

    ! TBD: Maybe change later so that it uses geographically dependend critical values (like wilting point)

    USE mo_soil_constants, ONLY: ws_crit_fract

    REAL(wp), INTENT(in) :: w_soil_sl  !< Soil water content
    REAL(wp), INTENT(in) :: w_soil_fc_sl   !< Soil water content at field capacity
    REAL(wp), INTENT(in) :: w_soil_pwp_sl  !< Soil water content at permanent wilting point

    REAL(wp) :: ws_crit             !< Soil water content at critical value

    ws_crit = w_soil_fc_sl * ws_crit_fract

    IF (ws_crit - w_soil_pwp_sl < EPSILON(1._wp)) THEN  ! This should only happen if there's no water in soil (ws=w_soil_fc_sl=w_soil_pwp_sl=0)
      get_water_stress_factor = 0._wp
    ELSE
      get_water_stress_factor = MAX (0._wp, &
        MIN (1._wp, (w_soil_sl - w_soil_pwp_sl) / (ws_crit - w_soil_pwp_sl) ))
    END IF

  END FUNCTION get_water_stress_factor

  REAL(wp) ELEMENTAL FUNCTION get_foliage_projected_cover(fpc_max, lai, clumpiness)

    REAL(wp), INTENT(in) :: &
      & fpc_max,            &
      & lai,                &
      & clumpiness

    get_foliage_projected_cover = fpc_max * (1._wp - EXP(-lai/clumpiness))

  END FUNCTION get_foliage_projected_cover


#endif
END MODULE mo_veg_process
