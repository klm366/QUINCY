!> Routines for physical schemes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains some routines for physical schemes
!>
MODULE mo_phy_schemes
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: lp
  USE mo_jsb_thermo_iface, ONLY: specific_humidity, sat_pres_water, sat_pres_ice, sat_pres_mixed, potential_temperature
#ifndef __NO_AES__
  USE mo_jsb_surface_exchange_iface, ONLY: sfc_exchange_coefficients
#endif
  USE mo_jsb_control, ONLY: acc_stream

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: specific_humidity, qsat_water, qsat_ice, qsat_mixed, sat_pres_water, sat_pres_ice, sat_pres_mixed, &
    &       update_drag, &
    &       q_effective, surface_dry_static_energy, heat_transfer_coef, &
    &       thermal_radiation, lwnet_from_lwdown, potential_temperature, &
    &       calc_peaked_arrhenius_function
#ifndef __NO_AES__
    PUBLIC :: sfc_exchange_coefficients
#endif

    ! Note: nvhpc doesn't support function pointers
    ! & register_exchange_coefficients_procedure, exchange_coefficients, registered_exchange_coefficients_procedure

  ! ABSTRACT INTERFACE
  !   PURE SUBROUTINE i_exchange_coefficients_procedure(           &
  !     & dz,                                                 &
  !     & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
  !     & km, kh, interp_fac_2m, interp_fac_10m               &
  !     & )
  !     IMPORT :: lp

  !     REAL(lp), INTENT(in) :: &
  !       dz,        &
  !       thetam1,   &
  !       pqm1 ,     &
  !       mwind,     &
  !       rough_m,   &
  !       theta_sfc, &
  !       qsat_sfc
  !       !
  !     REAL(lp), INTENT(out) :: &
  !       km,         &
  !       kh,         &
  !       interp_fac_2m, &
  !       interp_fac_10m

  !   END SUBROUTINE
  ! END INTERFACE

  ! PROCEDURE(i_exchange_coefficients_procedure), POINTER :: registered_exchange_coefficients_procedure => NULL()

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  ! SUBROUTINE register_exchange_coefficients_procedure(exchange_coefficients_procedure)

  !   PROCEDURE(i_exchange_coefficients_procedure) :: exchange_coefficients_procedure

  !   registered_exchange_coefficients_procedure => exchange_coefficients_procedure

  ! END SUBROUTINE register_exchange_coefficients_procedure

  ! SUBROUTINE exchange_coefficients(                       &
  !   & dz,                                                 &
  !   & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
  !   & km, kh, interp_fac_2m, interp_fac_10m               &
  !   & )

  !   REAL(lp), DIMENSION(:), INTENT(in) :: &
  !     dz,        &
  !     thetam1,   &
  !     pqm1 ,     &
  !     mwind,     &
  !     rough_m,   &
  !     theta_sfc, &
  !     qsat_sfc
  !     !
  !   REAL(lp), DIMENSION(:), INTENT(out) :: &
  !     km,         &
  !     kh,         &
  !     interp_fac_2m, &
  !     interp_fac_10m

  !   INTEGER :: ic, nc

  !   nc = SIZE(dz)

  !   !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
  !   DO ic = 1, nc
  !     CALL registered_exchange_coefficients_procedure(                                &
  !       & dz(ic),                                                                     &
  !       & pqm1(ic), thetam1(ic), mwind(ic), rough_m(ic), theta_sfc(ic), qsat_sfc(ic), &
  !       & km(ic), kh(ic), interp_fac_2m(ic), interp_fac_10m(ic)                       &
  !       & )
  !   END DO
  !   !$ACC END PARALLEL LOOP

  ! END SUBROUTINE exchange_coefficients

  ! ====================================================================================================== !
  !>
  !> Functions to calculate saturation specific humidity
  !>
  !> Input: temperature [K], pressure [Pa]
  !> Output: saturation specific humidity [kg/kg]
  !>
#ifndef _OPENACC
ELEMENTAL &
#endif
  REAL(lp) FUNCTION qsat_water(temperature, pressure)

    !$ACC ROUTINE SEQ

    REAL(lp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(lp), INTENT(in) :: pressure         ! Pressure [Pa]

    qsat_water = specific_humidity(sat_pres_water(temperature), pressure)

  END FUNCTION qsat_water

#ifndef _OPENACC
ELEMENTAL &
#endif
  REAL(lp) FUNCTION qsat_ice(temperature, pressure)

    !$ACC ROUTINE SEQ

    REAL(lp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(lp), INTENT(in) :: pressure         ! Pressure [Pa]

    qsat_ice = specific_humidity(sat_pres_ice(temperature), pressure)

  END FUNCTION qsat_ice

#ifndef _OPENACC
ELEMENTAL &
#endif
  REAL(lp) FUNCTION qsat_mixed(temperature, pressure)

    !$ACC ROUTINE SEQ

    REAL(lp), INTENT(in) :: temperature      !< Air temperature [K]
    REAL(lp), INTENT(in) :: pressure         !< Pressure [Pa]

    qsat_mixed = specific_humidity(sat_pres_mixed(temperature), pressure)

  END FUNCTION qsat_mixed

  !> ====================================================================================================== !

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION q_effective
  !!
  !! out: q_effective
  !-----------------------------------------------------------------------------------------------------
  REAL(lp) FUNCTION q_effective(qsat, qair, fsat, fair)
!NEC$ always_inline

    !$ACC ROUTINE SEQ

    REAL(lp), INTENT(in) :: &
      & qsat,               & !< Surface saturation specific humidity
      & qair,               & !< Air specific humidity
      & fsat, fair            !< Weighing factors for qsat and qair accounting for only partially wet surface.

    q_effective = fsat * qsat + (1._lp -fair) * qair

  END FUNCTION q_effective

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION surface_dry_static_energy
  !!
  !! out: surface_dry_static_energy
  !-----------------------------------------------------------------------------------------------------
  REAL(lp) FUNCTION surface_dry_static_energy(t_srf, qsat_srf, cpd_or_cvd, jsb_standalone)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: &
      & cpvd1        !< cpv/cpd-1

    REAL(lp), INTENT(in) :: &
      & t_srf,                      & !< Surface temperature
      & qsat_srf,                   & !< Surface saturation specific humidity
      & cpd_or_cvd

    LOGICAL, INTENT(in) :: &
      & jsb_standalone

    IF (jsb_standalone) THEN
      ! Todo: Check if this formulation is valid. We re-implement it here, as it was used up to now in
      !       jsbach:dev, and we do not want to change the model behavior with this merge.
      surface_dry_static_energy = t_srf * cpd_or_cvd * ( 1._lp + cpvd1 * qsat_srf)
    ELSE
      surface_dry_static_energy = t_srf * cpd_or_cvd
    END IF

  END FUNCTION surface_dry_static_energy


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION thermal_radiation
  !!
  !! out: thermal_radiation
  !-----------------------------------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(lp) FUNCTION thermal_radiation(t_srf)

    USE mo_jsb_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(lp), INTENT(in) :: t_srf

    thermal_radiation = stbo * zemiss_def * t_srf**4._lp

  END FUNCTION thermal_radiation


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION lwnet_from_lwdown
  !!
  !! out: lwnet_from_lwdown
  !-----------------------------------------------------------------------------------------------------
  REAL(lp) FUNCTION lwnet_from_lwdown(lwdown, t_srf)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(lp), INTENT(in) :: lwdown   ! downward longwave radiation
    REAL(lp), INTENT(in) :: t_srf    ! surface temperature

    lwnet_from_lwdown = zemiss_def * (lwdown - stbo * t_srf**4._lp)

  END FUNCTION lwnet_from_lwdown


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION heat_transfer_coef
  !!
  !! out: heat_transfer_coef
  !-----------------------------------------------------------------------------------------------------
  REAL(lp) FUNCTION heat_transfer_coef(drag, dtime, alpha)
!NEC$ always_inline

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: grav

    REAL(lp), INTENT(in) :: &
      & drag,      &
      & dtime,   &
      & alpha

    heat_transfer_coef = drag / (alpha * grav * dtime)

  END FUNCTION heat_transfer_coef

  ! ======================================================================================================= !
  !> Calculates temperature response factor for biological processes
  !>   according to the peaked Arrhenius equation as defined in:
  !>   Medlyn et al. 2002, PCE, eq. 18
  !>
  !>   NB if sensitivity of zero is specified, 1 will be returned
  !>
  !>  Input: temperature, (de-)activiation energy, optimum temperature
  !>
  !>  Output: rate modifier (unitless)
  !>
  ELEMENTAL FUNCTION calc_peaked_arrhenius_function(temp, Ea, Ed, temp_opt) RESULT(rate_modifier)

    USE mo_jsb_physical_constants,  ONLY: r_gas, Tzero
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: temp_ref_tresponse, temp_freeze_thres_bio, temp_freeze_window_bio

    !------------------------------------------------------------------------------------------------------ !
    REAL(lp), INTENT(in)  :: temp           !< temperature [K]
    REAL(lp), INTENT(in)  :: Ea             !< activation energy [J mol-1]
    REAL(lp), INTENT(in)  :: Ed             !< deactivation energy [J mol-1]
    REAL(lp), INTENT(in)  :: temp_opt       !< optimum temperature [K]
    REAL(lp)              :: rate_modifier  !< temperature rate modifier [unitless]
    REAL(lp)              :: frost_modifier !< reduction of rate modifier when approaching -2 degC [-]
                                            !!   only for biological processes
    !------------------------------------------------------------------------------------------------------ !
    REAL(lp)              :: hlp1, hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_peaked_arrhenius_function'

    IF (Ea > eps8) THEN
      hlp1            = temp - temp_opt
      hlp2            = temp * temp_opt * r_gas
      frost_modifier  = MIN(MAX((temp - (Tzero - temp_freeze_thres_bio)) / temp_freeze_window_bio, 0.0_lp), 1.0_lp)
      rate_modifier   = (Ed * EXP(Ea * hlp1 / hlp2)) / (Ed - Ea * (1._lp - EXP(Ed * hlp1 / hlp2)))
      rate_modifier   = rate_modifier * frost_modifier
    ELSE
      rate_modifier = 1._lp
    END IF

  END FUNCTION calc_peaked_arrhenius_function

  ! ======================================================================================================= !
  !>
  !> calculates surface drag and exchange coefficients for jsbach standalone simulation
  !>
  SUBROUTINE update_drag(nc, time_step_len, &
      & config_model_scheme, &
      & temp_air, pressure, qair, wind, &
      & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, &
      & height_wind, height_humidity, coef_ril_tm1, coef_ril_t, coef_ril_tp1, &
      & drag_srf, ddrag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, &
      & veg_height, t_srf_upd, zril_old)

    USE mo_jsb_model_class,        ONLY: MODEL_QUINCY, MODEL_JSBACH
    USE mo_jsb_physical_constants, ONLY: grav, rd, cpd, rvd1, cpvd1, von_karman, rd_o_cpd

    INTEGER,  INTENT(IN)  :: nc
    REAL(lp), INTENT(IN)  :: time_step_len
    INTEGER,  INTENT(IN)  :: config_model_scheme
    REAL(lp), INTENT(IN)  :: temp_air(:)
    REAL(lp), INTENT(IN)  :: pressure(:)
    REAL(lp), INTENT(IN)  :: qair(:)
    REAL(lp), INTENT(IN)  :: wind(:)
    REAL(lp), INTENT(IN)  :: t_srf_proc(:)             !< Surface temperature
    REAL(lp), INTENT(IN)  :: fact_q_air_proc(:)
    REAL(lp), INTENT(IN)  :: fact_qsat_srf_proc(:)
    REAL(lp), INTENT(IN)  :: rough_h_srf_proc(:)
    REAL(lp), INTENT(IN)  :: rough_m_srf_proc(:)
    REAL(lp), INTENT(in)  :: height_wind               !< Defines lowest layer height-> where wind measurements are taken
    REAL(lp), INTENT(in)  :: height_humidity           !< Defines lowest layer height-> where humidity measurements are taken
    REAL(lp), INTENT(in)  :: coef_ril_tm1              !< Weighting factor for richardson numbers at different steps
    REAL(lp), INTENT(in)  :: coef_ril_t                !<   that are used to calculate a drag coef. that approximates
    REAL(lp), INTENT(in)  :: coef_ril_tp1              !<   the drag coef. at time t, but helps to maintain stability.

    REAL(lp), INTENT(OUT) :: drag_srf(:)               !< Surface drag
    REAL(lp), INTENT(OUT) :: ddrag_srf(:)              !< d(Surface drag)/d(t_srf_proc)
    REAL(lp), INTENT(OUT) :: t_acoef (:)
    REAL(lp), INTENT(OUT) :: t_bcoef (:)
    REAL(lp), INTENT(OUT) :: q_acoef (:)
    REAL(lp), INTENT(OUT) :: q_bcoef (:)
    REAL(lp), INTENT(OUT) :: pch     (:)

    ! Optional variable for QUINCY
    REAL(lp), OPTIONAL, INTENT(IN)    :: veg_height(:)     !< vegetation height
    ! Optional variables for the computation of a filtered Richardson number
    REAL(lp), OPTIONAL, INTENT(IN)    :: t_srf_upd(:)  !< Updated surface temperature
    REAL(lp), OPTIONAL, INTENT(inout) :: zril_old(:)   !< Richardson number at previous time step

    REAL(lp) :: height(nc)
    INTEGER  :: ic
    LOGICAL  :: l_zril_filter
    REAL(lp) :: zcons9, zcons11, zcons12, zsigma
    REAL(lp) :: air_pressure, zdu2, ztvd, ztvir, qsat_surf, qsat_surf_upd, r_exner_surf, qsurf_eff
    REAL(lp) :: ztvs_act, ztvs_upd, zg, zgh, zril, zril_act, zril_upd
    REAL(lp) :: zcons, zchnl, zcfnchl, zcfhl, sfact, zscfl, ufact, zucfhl
    !REAL(lp) :: zcons8, zcdnl, zucfl, zcfncl, zcfml   ! quincy specific

    REAL(lp) :: dqsat_surf, dqsurf_eff, dztvs_act, dzril, dzscfl, dzucfhl, dzcfhl

    REAL(lp), PARAMETER :: cb = 5._lp
    REAL(lp), PARAMETER :: cc = 5._lp
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_drag'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(height)

    ! in jsb3: rd = GasConstantDryAir; cpd = SpecificHeatDryAirConstPressure; grav = Gravity
    ! zcons8  = 2._lp * cb  ! quincy specific
    zcons9  = 3._lp * cb
    zcons11 = 3._lp * cb * cc
    zcons12 = time_step_len * grav / rd
    zsigma  = 0.99615_lp        ! corresponds to lowest level of 47 layer echam

    IF (PRESENT(t_srf_upd) .AND. PRESENT(zril_old)) THEN
      l_zril_filter = .TRUE.
    ELSE
      l_zril_filter = .FALSE.
    END IF
    IF (PRESENT(veg_height)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        height(ic) = veg_height(ic)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        height(ic) = 0._lp
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(qsat_surf, dqsat_surf, qsat_surf_upd, zdu2, air_pressure, ztvd, ztvir)       &
    !$ACC   PRIVATE(r_exner_surf, qsurf_eff, dqsurf_eff, ztvs_act, dztvs_act, ztvs_upd, zg, zgh) &
    !$ACC   PRIVATE(zril, dzril, zril_act, zril_upd, zchnl, sfact, zscfl, dzscfl, ufact, zucfhl) &
    !$ACC   PRIVATE(dzucfhl, zcons, zcfnchl, zcfhl, dzcfhl)
    DO ic = 1, nc
      ! NOTE: a segmentation fault for unplausible surface temperatures (outside [50,400] Kelvin) would occur in "calc_qsat()"
      qsat_surf = qsat_mixed(t_srf_proc(ic), pressure(ic))
      dqsat_surf = 1e3_lp * (qsat_mixed(t_srf_proc(ic) + 1e-3_lp, pressure(ic)) - qsat_surf)
      IF (l_zril_filter) THEN
        qsat_surf_upd = qsat_mixed(t_srf_upd(ic), pressure(ic))
      END IF

      !------------------------------------------------------------------------------------
      ! Approximation of cdrag
      !------------------------------------------------------------------------------------
      ! squared wind shear ! minimum wind speed square from echam/mo_surface_land.f90:precalc_land
      zdu2 = MAX(wind(ic) ** 2, 1._lp)

      ! virtual potential air temperature (see mo_surface_boundary.f90)
      ! according to Saucier, WJ Principles of Meteoroligical Analyses
      ! tv = t * (1 + 0.61 * q )    ! virtual temperature
      ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
      ! tvd = tair * (100000/p_pa)^rd_o_cps * (1 + rvd1 * qair) ! virtual potential temperature
      air_pressure = zsigma * pressure(ic)
      ztvd = ( temp_air(ic) * ( 100000._lp / air_pressure)**rd_o_cpd ) * ( 1._lp + rvd1 * qair(ic))
      ztvir = temp_air(ic) * ( 1._lp + rvd1 * qair(ic) )

      r_exner_surf = ( 100000._lp / pressure(ic))**rd_o_cpd
      qsurf_eff = fact_qsat_srf_proc(ic) * qsat_surf + ( 1._lp - fact_q_air_proc(ic) ) * qair(ic)
      dqsurf_eff = fact_qsat_srf_proc(ic) * dqsat_surf

      ! virtual potential surface temperature for actual and updated states
      ztvs_act = t_srf_proc(ic) * r_exner_surf * (1._lp + rvd1 * qsurf_eff)
      dztvs_act = r_exner_surf * (1._lp + rvd1 * (qsurf_eff + t_srf_proc(ic) * dqsurf_eff))
      IF (l_zril_filter) THEN
        ztvs_upd = t_srf_upd(ic) * ( 100000._lp / pressure(ic))**rd_o_cpd * &
          & ( 1._lp + rvd1 * ( fact_qsat_srf_proc(ic) * qsat_surf_upd + ( 1._lp - fact_q_air_proc(ic) ) &
          & * qair(ic)))
      END IF

      ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
      ! adapted according to jsb3 offline modifications by Marvin, Philipp et al. Jan. 24, 2017 (r8893)
      IF(height_wind > 0._lp .AND. height_humidity > 0._lp) THEN
        SELECT CASE (config_model_scheme)
        CASE (MODEL_JSBACH)
          zg  = height_wind * grav
          zgh = height_humidity * grav ! jsb3 "HeightHumidity and HeightTemperature have to be equal"
        CASE (MODEL_QUINCY)
          zg  = (height_wind + 3._lp / 10._lp * height(ic)) * grav
          zgh = zg
        END SELECT
      ELSE
        zg = ztvir * rd * LOG(1._lp / zsigma)
        zgh = zg
      ENDIF

      ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
      ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
      IF (l_zril_filter) THEN
        zril_act = zg * ( ztvd - ztvs_act ) / ( zdu2 * (ztvd + ztvs_act) / 2._lp )
        zril_upd = zg * ( ztvd - ztvs_upd ) / ( zdu2 * (ztvd + ztvs_upd) / 2._lp )

        ! Richardson number needs to be filtered for offline simulations in order to maintain stability
        zril = coef_ril_tm1 * zril_old(ic) &
          &  + coef_ril_t   * zril_act     &
          &  + coef_ril_tp1 * zril_upd
        zril_old(ic) = zril
        dzril = -4._lp * zg / zdu2 * ztvd / (ztvd + ztvs_act)**2 * dztvs_act * coef_ril_t
      ELSE
        SELECT CASE (config_model_scheme)
        CASE (MODEL_JSBACH)
          zril = zg * ( ztvd - ztvs_act ) / ( zdu2 * (ztvd + ztvs_act) / 2._lp )
        CASE (MODEL_QUINCY)
          zril = MIN(zg * ( ztvd - ztvs_act ) / ( zdu2 * (ztvd + ztvs_act) / 2._lp ), 0.5_lp)
        END SELECT
      dzril = -4._lp * zg / zdu2 * ztvd / (ztvd + ztvs_act)**2 * dztvs_act
      END IF

      ! Neutral drag coefficient for
      ! momentum
      ! zcdnl = (von_karman / LOG(1._lp + zg / (grav * rough_m_srf_proc(ic) ))) ** 2
      ! heat
      zchnl = von_karman ** 2 / (LOG(1._lp + zg / (grav * rough_m_srf_proc(ic) )) &
                * LOG( ( grav * rough_m_srf_proc(ic) + zgh ) / (grav * rough_h_srf_proc(ic) )))

      ! account for stable/unstable case: helper variables
      ! The stability factors are multiplied by Ri to cancel a singularity in the derivative of zucfhl.
      sfact = SQRT(1._lp + 5._lp * ABS(zril))
      zscfl = zril * sfact
      dzscfl = dzril * (sfact + 2.5_lp * zril / sfact)

      ufact = zcons11 * zchnl * SQRT(ABS(zril) * (1._lp + zgh / (grav * rough_h_srf_proc(ic))))
      ! momentum
      ! zucfl  = 1._lp / (1._lp + zcons11 * zcdnl * SQRT(ABS(zril) * (1._lp + zg  / (grav * rough_m_srf_proc(ic)))))
      ! heat
      zucfhl = - zril / (1._lp + ufact)
      dzucfhl = - dzril * (1._lp + 0.5_lp * ufact) / (1._lp + ufact)**2

      ! ignoring cloud water correction (see mo_surface_land.f90)
      zcons = zcons12 * pressure(ic) / ztvir
      ! momentum
      ! zcfncl   = zcons * SQRT(zdu2) * zcdnl
      ! heat
      zcfnchl  = zcons * SQRT(zdu2) * zchnl

      ! Stable / Unstable case
      IF ( zril > 0._lp ) THEN
        ! momentum
        ! zcfml = zcfncl  / (1._lp + zcons8 * zril / zscfl)
        ! heat
        zcfhl = zcfnchl / (1._lp + zcons9 * zscfl)
        dzcfhl = - zcfnchl * zcons9 * dzscfl / (1._lp + zcons9 * zscfl)**2
      ELSE
        ! zcfml = zcfncl  * (1._lp - zcons8 * zril * zucfl)
        zcfhl = zcfnchl * (1._lp + zcons9 * zucfhl)
        dzcfhl = zcfnchl * zcons9 * dzucfhl
      END IF

      drag_srf(ic) = zcfhl
      ddrag_srf(ic) = dzcfhl
      pch(ic) = zcfhl / zcfnchl * zchnl

      !---------------------------------------------------------------------------------------------------------------
      ! Computation of Richtmeyr-Morton Coefficients
      ! This follows now Jan Polcher's explicit solution, i.e. atmospheric conditions at t+1 are assumed to be valid
      !---------------------------------------------------------------------------------------------------------------
      t_acoef(ic) = 0.0_lp
      t_bcoef(ic) = cpd * (1._lp + cpvd1 * qair(ic)) * temp_air(ic) + zgh
      q_acoef(ic) = 0.0_lp
      q_bcoef(ic) = qair(ic)

    END DO
    !$ACC END PARALLEL LOOP
    !$ACC END DATA

  END SUBROUTINE update_drag

#endif
END MODULE mo_phy_schemes
