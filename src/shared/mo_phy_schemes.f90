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
  USE mo_jsb_convect_tables_iface, ONLY: init_convect_tables, prepare_ua_index_spline, lookup_ua_spline

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_convect_tables
  PUBLIC :: calc_qsat_from_esat, calc_qsat, qsat, q_effective, surface_dry_static_energy, heat_transfer_coef, &
    thermal_radiation, lwnet_from_lwdown

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  !>
  !! @brief Returns saturation specific humidity for given air pressure and saturated vapor pressure
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: specific_humidity    by Julia Nabel
  !
  ! === specific_humidity() =======================================================================================================
  !
  ! Specific humidity is computed here by assuming that each component of a mixture of dry air and vapor can be considered as an
  ! ideal gas. So, let m_d and m_v the masses of the two components in a fixed volume V. Then the specific humidity "q" is
  ! defined as
  ! (1) q = m_v/(m_v+m_d)
  ! Closely related is the moisture mixing ratio "r"
  ! (2) r = m_v/m_d
  ! Hence,
  ! (3) q = r/(1+r)
  ! Next it is assumed that both, dry air and vapor behave like ideal gases:
  ! (4) p_d*V = m_d*R_d*T and p_v*V = m_v*R_v*T,
  ! where p_d is the partial pressure of the dry air and p_v the vapor pressure, R_d and R_v the gas constants and T the
  ! temperature. Noting that air pressure "p" is the sum of the two partial pressures (p=p_d+p_v) the division of the two gas
  ! equation by each other gives by using (2):
  ! (5) r = eps*p_v/(p-p_v),
  ! where eps=R_d/R_v is the ratio of gas constants. Entering ths into (3) one finally obtains for the specific humidity:
  !                 p_v
  ! (6) q = eps -------------
  !             p-(1-eps)*p_v
  !
  PURE ELEMENTAL FUNCTION calc_qsat_from_esat(vapor_pressure, air_pressure)
  
    USE mo_jsb_physical_constants, ONLY: rdv ! = Rd/Rv

    REAL(wp),INTENT(in) :: vapor_pressure         !! sat. vapor pressure multiplied by Rd/Rv [N/m^2]
    REAL(wp),INTENT(in) :: air_pressure           !! air pressure at bottom [N/m^2]
    REAL(wp)            :: calc_qsat_from_esat    !! []

    !$claw sca routine

    calc_qsat_from_esat = rdv * vapor_pressure / (air_pressure - (1._wp- rdv) * vapor_pressure)
  
  END FUNCTION calc_qsat_from_esat

  !>
  ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
  !
  ! If l_esat_only present and = .TRUE., qsat will actually not contain the sat. specific humidity on return but
  ! the sat. vapor pressure multiplied by R_d/R_v
  !
  SUBROUTINE calc_qsat(caller, nc, temperature, pressure, qsat, iblk, blk_size, l_esat_only)

    USE mo_jsb_physical_constants,   ONLY: rvd1                    ! = R_v/R_d - 1

    CHARACTER(LEN=*), INTENT(in) :: caller
    INTEGER,  INTENT(in)  :: nc
    REAL(wp), INTENT(in)  :: temperature(nc)          ! Air temperature [K]
    REAL(wp), INTENT(in)  :: pressure(nc)             ! Pressure [Pa]
    REAL(wp), INTENT(out) :: qsat(nc)                 ! Saturation specific humidity
    INTEGER,  INTENT(in), OPTIONAL  :: iblk, blk_size
    LOGICAL,  INTENT(in), OPTIONAL  :: l_esat_only

    INTEGER  :: ic
    INTEGER  :: idx(nc)
    REAL(wp) :: za(nc), zua(nc), zes
    LOGICAL  :: l_esat

    CHARACTER(len=*), PARAMETER :: routine = 'mo_phy_schemes:calc_qsat'

    l_esat = .FALSE.
    IF (PRESENT(l_esat_only)) l_esat = l_esat_only

    !$ACC DATA &
    !$ACC PRESENT(temperature, pressure, qsat) &
    !$ACC CREATE(idx, za, zua)

    IF (PRESENT(iblk) .AND. PRESENT(blk_size)) THEN
      CALL prepare_ua_index_spline(caller, 1, nc, temperature(1:nc), idx(1:nc), za(1:nc), &
        &                          klev=0, kblock=iblk, kblock_size=blk_size)
    ELSE
      CALL prepare_ua_index_spline(caller, 1, nc, temperature(1:nc), idx(1:nc), za(1:nc))
    END IF
    CALL lookup_ua_spline(1,nc,idx(1:nc), za(1:nc), zua(1:nc))

    IF (l_esat) THEN
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        qsat(ic) = zua(ic)
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        zes = zua(ic) / pressure(ic)            ! sat. vapor pressure * Rd/Rv / pressure
        qsat(ic) = zes / (1._wp - rvd1 * zes)
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC END DATA

  END SUBROUTINE calc_qsat

  REAL(wp) ELEMENTAL FUNCTION qsat(temperature, pressure)

    USE mo_jsb_convect_tables_iface, ONLY: tlucua, &                       ! Table for water vapor pressure e_s multiplied by R_d/R_v
                                 jptlucu1, jptlucu2
    USE mo_jsb_physical_constants, ONLY: rvd1                    ! = R_v/R_d - 1

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]

    INTEGER :: it
    REAL(wp) :: tluc

    !$claw sca routine

    it = NINT(temperature*1000._wp)
    tluc = tlucua(it)

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
