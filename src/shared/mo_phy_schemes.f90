!>
!! @brief Contains some functions for physcial schemes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Adapted for SLM purposes (by Sönke Zaehle (2016-12-22)
!!
MODULE mo_phy_schemes
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  !PUBLIC :: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef, &
  PUBLIC :: q_effective, surface_dry_static_energy, heat_transfer_coef, &
    thermal_radiation, lwnet_from_lwdown

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

!  !>
!  ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
!  ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
!  ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
!  !
!  REAL(wp) ELEMENTAL FUNCTION qsat(temperature, pressure)
!
!    USE mo_echam_convect_tables, ONLY: tlucua, &                       ! Table for water vapor pressure e_s multiplied by R_d/R_v
!                                 jptlucu1, jptlucu2, lookuperror
!    USE mo_jsb_physical_constants, ONLY: rvd1                    ! = R_v/R_d - 1
!
!    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
!    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]
!
!    INTEGER :: it
!    REAL(wp) :: tluc
!
!    it = NINT(temperature*1000._wp)
!!!$    it = MIN(MAX(NINT(temperature*1000._wp), 50000),400000)    ! Temporary fix
!    tluc = tlucua(it) !<- TODO catch seg fault for unplausible temperatures? (e.g. see mo_jsb4_forcing)
!
!    IF (it >= jptlucu1 .AND. it <= jptlucu2) THEN
!      qsat = tluc / (pressure - rvd1*tluc)
!    ELSE
!      qsat = HUGE(1._wp)
!    END IF
!
!  END FUNCTION qsat

  REAL(wp) ELEMENTAL FUNCTION q_effective(qsat, qair, fsat, fair)

    REAL(wp), INTENT(in) :: &
      & qsat,               & !< Surface saturation specific humidity
      & qair,               & !< Air specific humidity
      & fsat, fair            !< Weighing factors for qsat and qair accounting for only partially wet surface.

    q_effective = fsat * qsat + (1._wp -fair) * qair

  END FUNCTION q_effective

  REAL(wp) ELEMENTAL FUNCTION surface_dry_static_energy(t_srf, qsat_srf)

!    USE mo_jsb_physical_constants, ONLY: &
    USE mo_shared_constants, ONLY: &
      & cpd,       & !< Specific heat of dry air at constant pressure
      & cpvd1        !< cpv/cpd-1

    REAL(wp), INTENT(in) :: &
      & t_srf,                      & !< Surface temperature
      & qsat_srf                      !< Surface saturation specific humidity

      surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf)

  END FUNCTION surface_dry_static_energy

  REAL(wp) ELEMENTAL FUNCTION thermal_radiation(t_srf)

!    USE mo_jsb_physical_constants, ONLY: &
    USE mo_shared_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: t_srf

    thermal_radiation = stbo * zemiss_def * t_srf**4

  END FUNCTION thermal_radiation

  REAL(wp) ELEMENTAL FUNCTION lwnet_from_lwdown(lwdown, t_srf)

!    USE mo_physical_constants, ONLY: &
    USE mo_shared_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: lwdown   ! downward longwave radiation
    REAL(wp), INTENT(in) :: t_srf    ! surface temperature

    lwnet_from_lwdown = zemiss_def * (lwdown - stbo * t_srf**4)

  END FUNCTION lwnet_from_lwdown

  REAL(wp) ELEMENTAL FUNCTION heat_transfer_coef(drag, steplen)

!    USE mo_jsb_physical_constants, ONLY: &
    USE mo_shared_constants, ONLY: &
      grav, cvdifts

    REAL(wp), INTENT(in) :: &
      & drag,      &
      & steplen

    heat_transfer_coef = drag / (cvdifts * grav * steplen)

  END FUNCTION heat_transfer_coef

#endif
END MODULE mo_phy_schemes
