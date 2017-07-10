!>
!! @brief Contains the routines for the plant respiration 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-02-17)
!!
MODULE mo_veg_respiration

  USE mo_kind


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_maintenance_respiration, temperature_response_respiration

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_respiration'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculate temperature response of maintenance respiration 
  !! following the temperature function of Lloyd and Taylor, 1996, Functional Ecology
  !!
  !! Input: temperature 
  !!
  !! Output: temperature rate modifier
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION temperature_response_respiration(temp, temp_mavg) RESULT(temp_response)

    USE mo_veg_constants,           ONLY: t_maint_k1, t_maint_k2, t_maint_k3, f_resp_acclim, t_acclim_zero

    IMPLICIT NONE

    ! 0.1 InOut
    REAL(wp), INTENT(in) :: temp                      !< instantaneous temperature (K)
    REAL(wp), INTENT(in) :: temp_mavg                 !< average temperature for acclimation (K) 
    REAL(wp)             :: temp_response             !< temperature response of respiration (1 at 10°C)
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':temperature_response_respiration'


    ! IF flag_t_resp_acclimation is false: temp_mavg == t_acclim_zero
    temp_response = (10._wp**(f_resp_acclim * (temp_mavg - t_acclim_zero))) * &         ! this term is 1 if flag_t_resp_acclimation is false
                     exp(t_maint_k1 * (1._wp/t_maint_k2 - 1._wp / (temp - t_maint_k3))) 

  END FUNCTION temperature_response_respiration

  !-----------------------------------------------------------------------------------------------------
  !> Calculate maintenance respiration
  !! following the base rate as defined by Sprugel, 1996, modified by Temperature according to 
  !!   Lloyd and Taylor, 1996. 
  !!
  !! Input: tissue N and temperature (leaf respiration calculated in photosynthesis and NOT modified
  !!        here)
  !!
  !! Output: maintenance respiration 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_maintenance_respiration( t_stem                , &
                                                        t_root                , &
                                                        temp_mavg              , &
                                                        sap_wood_nitrogen     , & 
                                                        fine_root_nitrogen    , &
                                                        coarse_root_nitrogen)   & 
                                                        RESULT(maint_respiration)

    USE mo_veg_constants, ONLY: fmaint_rate_base

    IMPLICIT NONE

    ! 0.1 InOut
    REAL(wp), INTENT(in) :: t_stem                   !< temperature of the stem (K) 
    REAL(wp), INTENT(in) :: t_root                   !< temperature of the roots (K)
    REAL(wp), INTENT(in) :: temp_mavg                 !< average temperature for acclimation (K) 
    REAL(wp), INTENT(in) :: sap_wood_nitrogen        !< nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: coarse_root_nitrogen     !< nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: fine_root_nitrogen       !< nitrogen content of sapwood (mol N / m2)
                                                     ! (calculated in photosynthesis; µmol C / m2 / s)
    REAL(wp)             :: maint_respiration        !< maintenance respiration (µmol C / m2 / s)
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_maintenance_respiration'


    maint_respiration = temperature_response_respiration(t_stem, temp_mavg) * & 
                        fmaint_rate_base * sap_wood_nitrogen + &
                        temperature_response_respiration(t_root, temp_mavg) * & 
                        fmaint_rate_base * ( coarse_root_nitrogen + fine_root_nitrogen ) 
     
  END FUNCTION calc_maintenance_respiration

END MODULE mo_veg_respiration
