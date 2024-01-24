!> Summary: QUINCY routines for physcial schemes
!>
!>#### ICON-Land
!>
!> ------------------------------------------------------------------------
!> Copyright (C) 2004-2024, MPI-M, MPI-BGC, DWD, DKRZ, KIT, ETH, MeteoSwiss  
!> Contact information: <https://icon-model.org>
!>
!> See AUTHORS.TXT for a list of authors  
!> See LICENSES/ for license information  
!> SPDX-License-Identifier: BSD-3-Clause
!> ------------------------------------------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### defines routines for physcial schemes adapted from JSBACH for QUINCY
!>
MODULE mo_phy_schemes
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: qsat, update_drag
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

    REAL(wp), INTENT(in) :: drag      !< ..
    REAL(wp), INTENT(in) :: steplen   !< ..
    REAL(wp), INTENT(in) :: alpha     !< differs for standalone & coupled simulations

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

  ! ======================================================================================================= !
  !>
  !> calculates surface drag and exchange coefficients for jsbach standalone simulation
  !>
  SUBROUTINE update_drag(nc, time_step_len, &
      & temp_air, pressure, qair, wind, &
      & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, &
      & height_wind, height_humidity, coef_ril_tm1, coef_ril_t, coef_ril_tp1, &
      & drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, &
      & t_srf_upd, zril_old)

    USE mo_jsb_physical_constants, ONLY: grav, rd, cpd, rvd1, cpvd1, von_karman, rd_o_cpd

    INTEGER,  INTENT(IN)  :: nc
    REAL(wp), INTENT(IN)  :: time_step_len
    REAL(wp), INTENT(IN)  :: temp_air(:)
    REAL(wp), INTENT(IN)  :: pressure(:)
    REAL(wp), INTENT(IN)  :: qair(:)
    REAL(wp), INTENT(IN)  :: wind(:)
    REAL(wp), INTENT(IN)  :: t_srf_proc(:)             !< Surface temperature
    REAL(wp), INTENT(IN)  :: fact_q_air_proc(:)
    REAL(wp), INTENT(IN)  :: fact_qsat_srf_proc(:)
    REAL(wp), INTENT(IN)  :: rough_h_srf_proc(:)
    REAL(wp), INTENT(IN)  :: rough_m_srf_proc(:)
    REAL(wp), INTENT(in)  :: height_wind               !< Defines lowest layer height-> where wind measurements are taken
    REAL(wp), INTENT(in)  :: height_humidity           !< Defines lowest layer height-> where humidity measurements are taken
    REAL(wp), INTENT(in)  :: coef_ril_tm1              !< Weighting factor for richardson numbers at different steps
    REAL(wp), INTENT(in)  :: coef_ril_t                !<   that are used to calculate a drag coef. that approximates
    REAL(wp), INTENT(in)  :: coef_ril_tp1              !<   the drag coef. at time t, but helps to maintain stability.

    REAL(wp), INTENT(OUT) :: drag_srf(:)               !< Surface drag
    REAL(wp), INTENT(OUT) :: t_acoef (:)
    REAL(wp), INTENT(OUT) :: t_bcoef (:)
    REAL(wp), INTENT(OUT) :: q_acoef (:)
    REAL(wp), INTENT(OUT) :: q_bcoef (:)
    REAL(wp), INTENT(OUT) :: pch     (:)

    ! Optional variables for the computation of a filtered Richardson number
    REAL(wp), OPTIONAL, INTENT(IN)    :: t_srf_upd(:)  !< Updated surface temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: zril_old(:)   !< Richardson number at previous time step

    INTEGER  :: ic
    LOGICAL  :: l_zril_filter
    REAL(wp) :: zcons9, zcons11, zcons12, zsigma
    REAL(wp) :: air_pressure, zdu2, ztvd, ztvir, qsat_surf(nc), qsat_surf_upd(nc)
    REAL(wp) :: ztvs_act, ztvs_upd, zg, zgh, zril, zril_act, zril_upd
    REAL(wp) :: zcons, zchnl, zcfnchl, zcfhl, zscfl, zucfhl

    REAL(wp), PARAMETER :: cb = 5._wp
    REAL(wp), PARAMETER :: cc = 5._wp
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_drag'

    ! in jsb3: rd = GasConstantDryAir; cpd = SpecificHeatDryAirConstPressure; grav = Gravity
    zcons9  = 3._wp * cb
    zcons11 = 3._wp * cb * cc
    zcons12 = time_step_len * grav / rd
    zsigma  = 0.99615_wp        ! corresponds to lowest level of 47 layer echam

    ! Initializations
    drag_srf(:) = 0._wp
    t_acoef (:) = 0._wp
    t_bcoef (:) = 0._wp
    q_acoef (:) = 0._wp
    q_bcoef (:) = 0._wp
    pch     (:) = 0._wp

    IF (PRESENT(t_srf_upd) .AND. PRESENT(zril_old)) THEN
      l_zril_filter = .TRUE.
    ELSE
      l_zril_filter = .FALSE.
    END IF

    qsat_surf(:)     = qsat(t_srf_proc(:), pressure(:))
    IF (l_zril_filter) THEN
      qsat_surf_upd(:) = qsat(t_srf_upd(:),  pressure(:))
    END IF

    DO ic = 1, nc

      !------------------------------------------------------------------------------------
      ! Approximation of cdrag
      !------------------------------------------------------------------------------------
      ! squared wind shear ! minimum wind speed square from echam/mo_surface_land.f90:precalc_land
      zdu2 = MAX(wind(ic)**2, 1._wp)

      ! virtual potential air temperature (see mo_surface_boundary.f90)
      ! according to Saucier, WJ Principles of Meteoroligical Analyses
      ! tv = t * (1 + 0.61 * q )    ! virtual temperature
      ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
      ! tvd = tair * (100000/p_pa)^rd_o_cps * (1 + rvd1 * qair) ! virtual potential temperature
      air_pressure = zsigma * pressure(ic)
      ztvd = ( temp_air(ic) * ( 100000._wp / air_pressure)**rd_o_cpd ) * ( 1._wp + rvd1 * qair(ic))
      ztvir = temp_air(ic) * ( 1._wp + rvd1 * qair(ic) )

      ! virtual potential surface temperature for actual and updated states
      ztvs_act = t_srf_proc(ic) * ( 100000._wp / pressure(ic))**rd_o_cpd * &
        & ( 1._wp + rvd1 * ( fact_qsat_srf_proc(ic) * qsat_surf(ic) + ( 1._wp - fact_q_air_proc(ic) ) &
        & * qair(ic)))
      IF (l_zril_filter) THEN
        ztvs_upd = t_srf_upd(ic) * ( 100000._wp / pressure(ic))**rd_o_cpd * &
          & ( 1._wp + rvd1 * ( fact_qsat_srf_proc(ic) * qsat_surf_upd(ic) + ( 1._wp - fact_q_air_proc(ic) ) &
          & * qair(ic)))
      END IF

      ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
      ! adapted according to jsb3 offline modifications by Marvin, Philipp et al. Jan. 24, 2017 (r8893)
      IF(height_wind > 0._wp .AND. height_humidity > 0._wp) THEN
        zg  = height_wind * grav
        zgh = height_humidity * grav ! jsb3 "HeightHumidity and HeightTemperature have to be equal"
      ELSE
        zg = ztvir * rd * LOG(1._wp / zsigma)
        zgh = zg
      ENDIF

      ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
      ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
      IF (l_zril_filter) THEN
        zril_act = zg * ( ztvd - ztvs_act ) / ( zdu2 * (ztvd + ztvs_act) / 2._wp )
        zril_upd = zg * ( ztvd - ztvs_upd ) / ( zdu2 * (ztvd + ztvs_upd) / 2._wp )

        ! Richardson number needs to be filtered for offline simulations in order to maintain stability
        zril = coef_ril_tm1 * zril_old(ic) &
          &  + coef_ril_t   * zril_act     &
          &  + coef_ril_tp1 * zril_upd
        zril_old(ic) = zril
      ELSE
        zril = zg * ( ztvd - ztvs_act ) / ( zdu2 * (ztvd + ztvs_act) / 2._wp )
      END IF

      ! Neutral drag coefficient for momentum and heat
      zchnl = von_karman**2 / (LOG(1._wp + zg / (grav * rough_m_srf_proc(ic) )) &
                * LOG( ( grav * rough_m_srf_proc(ic) + zgh ) / (grav * rough_h_srf_proc(ic) )))

      ! account for stable/unstable case: helper variables
      zscfl = SQRT ( 1._wp + 5._wp * ABS(zril))
      zucfhl = 1._wp / (1._wp + zcons11 * zchnl * SQRT(ABS(zril) * (1._wp &
        & + zgh / (grav * rough_h_srf_proc(ic)))))

      ! ignoring cloud water correction (see mo_surface_land.f90)
      zcons = zcons12 * pressure(ic) / ( temp_air(ic) * (1._wp + rvd1 * qair(ic)))
      zcfnchl  = zcons * SQRT(zdu2) * zchnl

      ! Stable / Unstable case
      IF ( zril > 0._wp ) THEN
        zcfhl = zcfnchl / (1._wp + zcons9 * zril * zscfl)
      ELSE
        zcfhl = zcfnchl * (1._wp - zcons9 * zril * zucfhl)
      END IF

      drag_srf(ic) = zcfhl
      pch(ic) = zcfhl / zcfnchl * zchnl

      !---------------------------------------------------------------------------------------------------------------
      ! Computation of Richtmeyr-Morton Coefficients
      ! This follows now Jan Polcher's explicit solution, i.e. atmospheric conditions at t+1 are assumed to be valid
      !---------------------------------------------------------------------------------------------------------------
      t_acoef(ic) = 0.0_wp
      t_bcoef(ic) = cpd * (1 + cpvd1 * qair(ic)) * temp_air(ic) + zgh
      q_acoef(ic) = 0.0_wp
      q_bcoef(ic) = qair(ic)

    END DO

  END SUBROUTINE update_drag

#endif
END MODULE mo_phy_schemes
