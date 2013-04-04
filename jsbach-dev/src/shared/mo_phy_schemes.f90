!>
!! @brief Contains some functions for physcial schemes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_phy_schemes

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: qsat, surface_dry_static_energy, thermal_radiation, heat_transfer_coef

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  !>
  ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
  ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
  ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
  !
  REAL(wp) ELEMENTAL FUNCTION qsat(temperature, pressure)

    USE mo_convect_tables, ONLY: tlucua, &                       ! Table for water vapor pressure e_s multiplied by R_d/R_v
                                 jptlucu1, jptlucu2, lookuperror
    USE mo_physical_constants, ONLY: vtmpc1                      ! = R_v/R_d - 1

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]

    INTEGER :: it
    REAL(wp) :: tluc

    it = NINT(temperature*1000._wp)
    tluc = tlucua(it)

    IF (it >= jptlucu1 .AND. it <= jptlucu2) THEN
      qsat = tluc / (pressure - vtmpc1*tluc)
    ELSE
      qsat = HUGE(1._wp)
    END IF

  END FUNCTION qsat

  REAL(wp) ELEMENTAL FUNCTION  surface_dry_static_energy(t_srf, qsat_srf, q_air, fact_qsat_srf, fact_q_air)

    USE mo_physical_constants, ONLY: &
      & cpd,       & !< Specific heat of dry air at constant pressure
      & vtmpc2       !< cpv/cpd-1

    REAL(wp), INTENT(in) :: &
      & t_srf,         &
      & qsat_srf,      &
      & q_air,         &
      & fact_qsat_srf, &
      & fact_q_air

    surface_dry_static_energy = t_srf * cpd * ( 1._wp + &
                                                & vtmpc2 * ( fact_qsat_srf * qsat_srf + (1._wp - fact_q_air) * q_air) &
                                              )

  END FUNCTION surface_dry_static_energy

  REAL(wp) ELEMENTAL FUNCTION thermal_radiation(t_srf)

    USE mo_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: t_srf

    thermal_radiation = stbo * zemiss_def * t_srf

  END FUNCTION thermal_radiation

  REAL(wp) ELEMENTAL FUNCTION heat_transfer_coef(drag, steplen)

    USE mo_physical_constants, ONLY: grav
    USE mo_echam_vdiff_params, ONLY: cvdifts  ! TBD: change to 1. for standalone JSBACH

    REAL(wp), INTENT(in) :: &
      & drag,      &
      & steplen

    heat_transfer_coef = drag / (cvdifts * grav * steplen)

  END FUNCTION heat_transfer_coef

END MODULE mo_phy_schemes
