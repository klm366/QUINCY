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
  USE mo_aes_thermo, ONLY: potential_temperature

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef, &
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
  ! PROCEDURE(i_exchange_coefficients_procedure), POINTER :: exchange_coefficients => NULL()

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  SUBROUTINE register_exchange_coefficients_procedure(exchange_coefficients_procedure)

    PROCEDURE(i_exchange_coefficients_procedure) :: exchange_coefficients_procedure

    registered_exchange_coefficients_procedure => exchange_coefficients_procedure
    ! exchange_coefficients => exchange_coefficients_procedure

  END SUBROUTINE register_exchange_coefficients_procedure

! #ifndef _OPENACC
!   ELEMENTAL &
! #endif
!   SUBROUTINE exchange_coefficients(                         &
!     & dz,                                                 &
!     & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
!     & km, kh                                              &
!     & )

!     REAL(wp), INTENT(in) :: &
!       dz,        &
!       thetam1,   &
!       pqm1 ,     &
!       mwind,     &
!       rough_m,   &
!       theta_sfc, &
!       qsat_sfc
!       !
!     REAL(wp), INTENT(out) :: &
!       kh,      &
!       km

!     CALL registered_exchange_coefficients_procedure(          &
!       & dz,                                                 &
!       & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
!       & km, kh                                              &
!       & )

!   END SUBROUTINE exchange_coefficients
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
    ! print*,'AAA thetam1 ',MINVAL(thetam1),MAXVAL(thetam1)
    ! print*,'BBB theta_sfc ',MINVAL(theta_sfc),MAXVAL(theta_sfc)

  END SUBROUTINE exchange_coefficients
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

      ! surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf)
      surface_dry_static_energy = t_srf * cpd

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
