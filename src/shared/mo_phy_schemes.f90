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
  USE mo_jsb_thermo_iface, ONLY: sat_specific_humidity, sat_pres_water

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sat_specific_humidity, sat_pres_water, q_effective, surface_dry_static_energy, heat_transfer_coef, &
    thermal_radiation, lwnet_from_lwdown

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

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
      & cvd, cvv !< Specific heats of dry and moist air at constant volume

    REAL(wp), INTENT(in) :: &
      & t_srf,                      & !< Surface temperature
      & qsat_srf                      !< Surface saturation specific humidity

      !$claw sca routine

      ! surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf)
      surface_dry_static_energy = t_srf * ( cvd + (cvv - cvd) * qsat_srf)

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
