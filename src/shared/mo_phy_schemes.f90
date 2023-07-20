!>
!! @brief Contains some functions for physcial schemes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Adapted for QUINCY purposes (by Sönke Zaehle (2016-12-22)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_phy_schemes
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qsat
  PUBLIC :: q_effective, surface_dry_static_energy, heat_transfer_coef, &
    &       thermal_radiation, lwnet_from_lwdown
  PUBLIC :: calc_peaked_arrhenius_function

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  !>
  ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
  ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
  ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
  !
  REAL(wp) ELEMENTAL FUNCTION qsat(temperature, pressure)

    USE mo_jsb_convect_tables_iface, ONLY: tlucua, &                       ! Table for water vapor pressure e_s multiplied by R_d/R_v
      &                                    jptlucu1, jptlucu2
    USE mo_jsb_physical_constants, ONLY: rvd1                    ! = R_v/R_d - 1

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]

    INTEGER :: it
    REAL(wp) :: tluc

    !$claw sca routine

    it = NINT(temperature*1000._wp)
!!$    it = MIN(MAX(NINT(temperature*1000._wp), 50000),400000)    ! Temporary fix
    tluc = tlucua(it) !<- TODO catch seg fault for unplausible temperatures? (e.g. see mo_jsb4_forcing)

    IF (it >= jptlucu1 .AND. it <= jptlucu2) THEN
      qsat = tluc / (pressure - rvd1*tluc)
    ELSE
      qsat = HUGE(1._wp)
    END IF

  END FUNCTION qsat


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
      & cpd,       & !< Specific heat of dry air at constant pressure
      & cpvd1        !< cpv/cpd-1

    REAL(wp), INTENT(in) :: &
      & t_srf,                      & !< Surface temperature
      & qsat_srf                      !< Surface saturation specific humidity

      !$claw sca routine

      surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf)

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

    USE mo_jsb_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: lwdown   ! downward longwave radiation
    REAL(wp), INTENT(in) :: t_srf    ! surface temperature

    !$claw sca routine

    lwnet_from_lwdown = zemiss_def * (lwdown - stbo * t_srf**4._wp)

  END FUNCTION lwnet_from_lwdown


  ! ======================================================================================================= !
  !> calc heat_transfer_coef
  !>
  REAL(wp) ELEMENTAL FUNCTION heat_transfer_coef(drag, steplen, alpha)

    USE mo_jsb_physical_constants,      ONLY: grav

    REAL(wp),           INTENT(in) :: drag      !< ..
    REAL(wp),           INTENT(in) :: steplen   !< ..
    REAL(wp), OPTIONAL, INTENT(in) :: alpha     !< differs for standalone & coupled simulations

    !$claw sca routine
    heat_transfer_coef = drag / (alpha * grav * steplen)
  END FUNCTION heat_transfer_coef

  ! ======================================================================================================= !
  !>Calculates temperature response factor according to the peaked Arrhenius equation as defined in
  !>  Medlyn et al. 2002, PCE, eq. 18
  !>  NB if sensitivity of zero is specified, 1 will be returned
  !>
  !>  Input: temperature, (de-)activiation energy, optimum temperature
  !>
  !>  Output: rate modifier (unitless) 
  !>
  ELEMENTAL FUNCTION calc_peaked_arrhenius_function(temp,Ea,Ed,t_opt) RESULT(rate_modifier)

    USE mo_jsb_physical_constants,  ONLY: r_gas
    USE mo_jsb_math_constants,      ONLY: eps8
    !------------------------------------------------------------------------------------------------------ !
    REAL(wp), INTENT(in)  :: temp          !< temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp), INTENT(in)  :: Ed            !< deactivation energy [J mol-1]
    REAL(wp), INTENT(in)  :: t_opt         !< optimum temperature [K] 
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless] 
    !------------------------------------------------------------------------------------------------------ !
    REAL(wp)              :: hlp1, hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_peaked_arrhenius_function'

    IF(Ea > eps8) THEN
      hlp1 = temp - t_opt
      hlp2 = temp * t_opt * r_gas
      rate_modifier = (Ed * EXP(Ea * hlp1 / hlp2)) / (Ed - Ea * (1._wp - EXP(Ed * hlp1 / hlp2))) 
    ELSE
      rate_modifier = 1._wp
    ENDIF   

  END FUNCTION calc_peaked_arrhenius_function 

#endif
END MODULE mo_phy_schemes
