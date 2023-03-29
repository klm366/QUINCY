!>
!! @brief Contains some functions for physcial schemes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_phy_schemes
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_thermo_iface, ONLY: specific_humidity, sat_pres_water, sat_pres_ice, potential_temperature

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: specific_humidity, qsat_water, qsat_ice, sat_pres_water, &
    & q_effective, surface_dry_static_energy, heat_transfer_coef, &
    & thermal_radiation, lwnet_from_lwdown, potential_temperature, &
    & register_exchange_coefficients_procedure, exchange_coefficients

  ABSTRACT INTERFACE
    PURE SUBROUTINE i_exchange_coefficients_procedure(           &
      & dz,                                                 &
      & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
      & km, kh                                              &
      & )
      IMPORT :: wp

      REAL(wp), INTENT(in) :: &
        dz,        &
        thetam1,   &
        pqm1 ,     &
        mwind,     &
        rough_m,   &
        theta_sfc, &
        qsat_sfc
        !
      REAL(wp), INTENT(out) :: &
        kh,      &
        km

    END SUBROUTINE
  END INTERFACE

  PROCEDURE(i_exchange_coefficients_procedure), POINTER :: registered_exchange_coefficients_procedure => NULL()

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  SUBROUTINE register_exchange_coefficients_procedure(exchange_coefficients_procedure)

    PROCEDURE(i_exchange_coefficients_procedure) :: exchange_coefficients_procedure

    registered_exchange_coefficients_procedure => exchange_coefficients_procedure

  END SUBROUTINE register_exchange_coefficients_procedure

  SUBROUTINE exchange_coefficients(                         &
    & dz,                                                 &
    & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
    & km, kh                                              &
    & )

    REAL(wp), DIMENSION(:), INTENT(in) :: &
      dz,        &
      thetam1,   &
      pqm1 ,     &
      mwind,     &
      rough_m,   &
      theta_sfc, &
      qsat_sfc
      !
    REAL(wp), DIMENSION(:), INTENT(out) :: &
      kh,      &
      km

    INTEGER :: ic

    DO ic=1,SIZE(dz)
      CALL registered_exchange_coefficients_procedure(          &
        & dz(ic),                                                 &
        & pqm1(ic), thetam1(ic), mwind(ic), rough_m(ic), theta_sfc(ic), qsat_sfc(ic), &
        & km(ic), kh(ic)                                              &
        & )
    END DO

  END SUBROUTINE exchange_coefficients

  REAL(wp) ELEMENTAL FUNCTION qsat_water(temperature, pressure)

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]

    qsat_water = specific_humidity(sat_pres_water(temperature), pressure)

  END FUNCTION qsat_water

  REAL(wp) ELEMENTAL FUNCTION qsat_ice(temperature, pressure)

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]

    qsat_ice = specific_humidity(sat_pres_ice(temperature), pressure)

  END FUNCTION qsat_ice

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION q_effective
  !! 
  !! out: q_effective
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) ELEMENTAL FUNCTION q_effective(qsat, qair, fsat, fair)

    REAL(wp), INTENT(in) :: &
      & qsat,               & !< Surface saturation specific humidity
      & qair,               & !< Air specific humidity
      & fsat, fair            !< Weighing factors for qsat and qair accounting for only partially wet surface.

    !$claw sca routine

    q_effective = fsat * qsat + (1._wp -fair) * qair

  END FUNCTION q_effective


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION surface_dry_static_energy
  !! 
  !! out: surface_dry_static_energy
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) ELEMENTAL FUNCTION surface_dry_static_energy(t_srf, qsat_srf)

    USE mo_jsb_physical_constants, ONLY: &
      & cvd,       & !< Specific heat of dry air at constant volume
      & cpd,       & !< Specific heat of dry air at constant pressure
      & cpvd1        !< cpv/cpd-1

    REAL(wp), INTENT(in) :: &
      & t_srf,                      & !< Surface temperature
      & qsat_srf                      !< Surface saturation specific humidity

      !$claw sca routine

      ! surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf)
      surface_dry_static_energy = t_srf * cpd
      ! surface_dry_static_energy = t_srf * cvd

  END FUNCTION surface_dry_static_energy


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION thermal_radiation
  !! 
  !! out: thermal_radiation
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) ELEMENTAL FUNCTION thermal_radiation(t_srf)

    USE mo_jsb_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: t_srf

    thermal_radiation = stbo * zemiss_def * t_srf**4._wp

  END FUNCTION thermal_radiation


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION lwnet_from_lwdown
  !! 
  !! out: lwnet_from_lwdown
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) ELEMENTAL FUNCTION lwnet_from_lwdown(lwdown, t_srf)

    USE mo_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: lwdown   ! downward longwave radiation
    REAL(wp), INTENT(in) :: t_srf    ! surface temperature

    !$claw sca routine

    lwnet_from_lwdown = zemiss_def * (lwdown - stbo * t_srf**4._wp)

  END FUNCTION lwnet_from_lwdown


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION heat_transfer_coef
  !! 
  !! out: heat_transfer_coef
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) ELEMENTAL FUNCTION heat_transfer_coef(drag, steplen)

    USE mo_jsb_physical_constants, ONLY: grav, cvdifts

    REAL(wp), INTENT(in) :: &
      & drag,      &
      & steplen

    !$claw sca routine

    heat_transfer_coef = drag / (cvdifts * grav * steplen)

  END FUNCTION heat_transfer_coef

#endif
END MODULE mo_phy_schemes
