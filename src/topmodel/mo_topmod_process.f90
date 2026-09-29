!>
!> Contains the routines for the topmodel wetland diagnostics
!>
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>
!> author
!>  Stiig Wilkenskjeld, MPI-M Hamburg
!>
!> Revision History
!> Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!>

MODULE mo_top_process
#ifndef __NO_JSBACH__

  USE mo_kind,     ONLY: wp
  USE mo_top_util, ONLY: gamm_p, gamm_inv_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_topmodel_runmean, update_peat_runmean
  PUBLIC :: update_water_table, update_water_balance_peatland, update_inun_frac, update_aet
  PUBLIC :: get_climatological_inun

  CHARACTER(len=*), PARAMETER :: modname = 'mo_top_process'

CONTAINS

  !>
  !      Updates equilibrium evapotranspiration (aet)
  !
  !>
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE update_aet( &
    & tmelt, &
    & temp_air, net_radiation, &
    & aet &
  )
  !$acc routine seq

    ! Interface parameters
    REAL(wp), INTENT(in)  :: tmelt                  ! Melting point of water [K]
    REAL(wp), INTENT(in)  :: temp_air               ! Air temperature [K]
    REAL(wp), INTENT(in)  :: net_radiation          ! Net surface radiation [W m-2]
    REAL(wp), INTENT(out) :: aet                    ! Equilibrium evapotranspiration flux [kg m-2 s-1]

    ! Locals
    REAL(wp) :: t_c, t_ref, gamma, lambda, s

    aet = 0._wp;
    IF (net_radiation > 0._wp) THEN
      t_c   = temp_air - tmelt                      ! convert temperature to Celsius
      t_ref = t_c + 237.3_wp

      gamma = (65.05_wp + t_c * 0.064_wp)           ! psychrometer constant [Pa/K], adjusted for temperature
      lambda = (2.495e6_wp - t_c * 2380._wp)        ! latent heat of vaporization of water [J/kg], adjusted for temperature

      s = 2.503e6_wp * EXP(17.269_wp * t_c/t_ref)/(t_ref*t_ref) ! saturated vapour pressure

      aet = ((s/(s+gamma))/lambda) * net_radiation  ! equilibrium evapotranspiration

    ENDIF

  END SUBROUTINE update_aet

  !>
  !      Updates the water balance for individual peat tile
  !
  !>
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE update_water_balance_peatland(           &
    & tspd, delta_time,                                    &
    & peat_runoff_fac, peat_runoff_exp,                    &
    & precipitation, frac_liquid_water, frac_snow,         &
    & melt, aet,                                           &
    & peat_water_volume, peat_water_table,                 &
    & peat_water_balance                                   &
  )
  !$acc routine seq

    ! Interface parameters
    INTEGER,  INTENT(in)    :: tspd                    ! Time steps per day
    REAL(wp), INTENT(in)    :: delta_time              ! Time step [s]
    REAL(wp), INTENT(in)    :: peat_runoff_fac
    REAL(wp), INTENT(in)    :: peat_runoff_exp
    REAL(wp), INTENT(in)    :: precipitation           ! Total precipitation [m/s]
    REAL(wp), INTENT(in)    :: frac_liquid_water       ! Fraction of total water liquid (from last timestep)
    REAL(wp), INTENT(in)    :: frac_snow               ! Fraction of tile covered with snow
    REAL(wp), INTENT(in)    :: melt                    ! Flux of melting snow [kg m-2 s-1]
    REAL(wp), INTENT(in)    :: aet                     !
    REAL(wp), INTENT(inout) :: peat_water_volume       ! Volume of water in peat [m]
    REAL(wp), INTENT(inout) :: peat_water_table        ! Height of water table [m]
    REAL(wp), INTENT(inout) :: peat_water_balance      ! Rate change of peat water [kg m-2 s-1]

    ! Parameters
    REAL(wp), PARAMETER :: acro_por = 0.9_wp - 0.08_wp     ! acrotelm_por - Fgas
    REAL(wp), PARAMETER :: minvtot = 0.25_wp               ! minimum water content at surface [mm3 mm-3]
    REAL(wp), PARAMETER :: max_h_acro = 0.3_wp             ! Wania's z1
    REAL(wp), PARAMETER :: zmin = 0.1_wp
    REAL(wp), PARAMETER :: max_wt_height = 0.02_wp         ! Wania's maxh (previously 0.1_wp)
    REAL(wp), PARAMETER :: az = (acro_por - minvtot)/zmin
    REAL(wp), PARAMETER :: eff_porespace = 0.9_wp - minvtot
!     REAL(wp), PARAMETER :: loc_peat_runoff_exp = 10._wp
!     REAL(wp), PARAMETER :: loc_peat_runoff_fac = 1.e-4_wp

    ! Locals
    REAL(wp) :: runoff                                 ! runoff in m/timestep
    REAL(wp) :: aet_modifier

    runoff = 0._wp

    IF (peat_water_volume > max_h_acro * eff_porespace) THEN
      peat_water_table = peat_water_volume - max_h_acro * eff_porespace ! wt = Vol above surface
      IF (peat_water_table > max_wt_height) THEN
        runoff = runoff + peat_water_table - max_wt_height
      ENDIF
    ELSE
      peat_water_table = -(max_h_acro*eff_porespace - peat_water_volume)/eff_porespace
    ENDIF

! StW: Are we using the right peat_runoff_{exp,fac}'s here
!     runoff = runoff + loc_peat_runoff_fac / tspd * EXP(peat_water_volume * loc_peat_runoff_exp)
!     IF (frac_snow < 0.8_wp) THEN
!       runoff = runoff + peat_runoff_fac / tspd * EXP(peat_water_table * peat_runoff_exp)
!     ENDIF
    runoff = runoff + frac_liquid_water * peat_runoff_fac / tspd * EXP(peat_water_table * peat_runoff_exp)

    ! determine water balance: use modified aet for evap, assume than only non-snow-covered ground is in connection with atmosphere
!     aet_modifier = 1.02_wp + 0.75_wp * peat_water_table
    aet_modifier = MAX(0._wp,(1.02_wp + 0.15_wp * TAN(peat_water_table)/(0.7_wp*max_h_acro)))

    peat_water_balance = (1._wp - frac_snow) * (precipitation - aet * aet_modifier) + melt - runoff * 1000._wp / delta_time
    peat_water_volume  = MAX(0._wp,peat_water_volume + peat_water_balance / 1000._wp * delta_time)

  END SUBROUTINE update_water_balance_peatland

  !      Updates the water table for individual tile
  !
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE update_water_table(     &
    &  with_permafrost, tmelt,            &
    &  soil_depth, soil_depth_min,        &
    &  dz_sl, lower_bound_sl,             &
    &  saturation_threshold,              &
    &  field_capacity,                    &
    &  t_soil_sl,                         &
    &  wtr_soil_sl, ice_soil_sl,          &
    &  water_table,                       &
    &  frac_liquid_top                    &
    & )
  !$acc routine seq

    ! Parameters
    LOGICAL,  INTENT(in)    :: with_permafrost      ! Soil may contain frozen water?
    REAL(wp), INTENT(in)    :: tmelt                ! Melting point of water [K]
    REAL(wp), INTENT(in)    :: soil_depth           ! Total depth of soil [m]
    REAL(wp), INTENT(in)    :: soil_depth_min       ! Minimum total depth of soil for point to be considered [m]
    REAL(wp), INTENT(in)    :: dz_sl(:)             ! Thickness of individual soil layers [m]
    REAL(wp), INTENT(in)    :: lower_bound_sl(:)    ! Depth of lower boundary of individual soil layers [m]
    REAL(wp), INTENT(in)    :: saturation_threshold ! Threshold for fully saturated soil layer []
    REAL(wp), INTENT(in)    :: field_capacity       ! Field capacity of soil []
    ! Input variables
    REAL(wp), INTENT(in)    :: t_soil_sl(:)         ! Temperature of each soil layer [K]
    REAL(wp), INTENT(in)    :: wtr_soil_sl(:)       ! Column of liquid water in each soil layer [m]
    REAL(wp), INTENT(in)    :: ice_soil_sl(:)       ! Column of frozen water in each soil layer [m]
    ! Output variables
    REAL(wp), INTENT(out)   :: water_table          ! Height of water table [m]
    REAL(wp), INTENT(out)   :: frac_liquid_top      ! fraction of liquid water in top <=0.5m (full layers) []

    ! Local variables
    REAL(wp) :: weq
    REAL(wp) :: saturation
    REAL(wp) :: dz                          ! Thickness of current layer
    REAL(wp) :: local_depth
    REAL(wp) :: frac_liquid
    LOGICAL :: layer_below_saturated, layer_below_unfrozen
    INTEGER :: nsoil, isoil

    nsoil = SIZE(dz_sl)
    layer_below_saturated =  .TRUE.
    layer_below_unfrozen  = .FALSE.
    local_depth = 0._wp
    water_table = -MIN(lower_bound_sl(nsoil),soil_depth)
    frac_liquid_top    = 0._wp
    DO isoil = nsoil,1,-1
      ! StW: Cmp with JSB3 version: glaciers and lakes should not appear on this tile so we can just start calc.

      IF (lower_bound_sl(MAX(1,isoil-1)) > soil_depth) CYCLE ! No soil in this layer

      ! Thickness of soil in current layer
      dz = dz_sl(isoil)
      IF (lower_bound_sl(isoil) > soil_depth) dz = soil_depth - lower_bound_sl(MAX(1,isoil-1))

      ! Absolute and relative water content (liquid and frozen) in layer isoil
      frac_liquid = 1._wp
      weq = wtr_soil_sl(isoil)
      IF (with_permafrost) weq = weq + ice_soil_sl(isoil)
      IF (weq > 0._wp) THEN
        IF (with_permafrost) THEN
          frac_liquid = wtr_soil_sl(isoil) / weq
        ELSE
          ! Piecewise linear function, 0 for T < ~0.5 degC, increasing to 1 between ~-0.5 degC and 0.5 degC and 1 afterwards
          frac_liquid = MAX(0._wp, MIN(1._wp, (t_soil_sl(isoil) - (tmelt-0.5_wp))))
        ENDIF
      ENDIF

      ! Sum for avg. liquid frac in upper 0.5 m
      IF (lower_bound_sl(isoil) <= 0.5_wp) THEN
        frac_liquid_top = frac_liquid_top + frac_liquid * dz
        local_depth = local_depth + dz
      ENDIF

      ! Water table for grids with sufficient soil_depth
      IF (soil_depth > soil_depth_min) THEN
        saturation = MAX(0._wp, MIN(1._wp, weq / dz / field_capacity))
        IF (layer_below_saturated) THEN
          IF ((t_soil_sl(isoil) < tmelt) .AND. .NOT. layer_below_unfrozen) THEN
            water_table = water_table + dz
          ELSE
            layer_below_unfrozen = .TRUE.
            IF (saturation > saturation_threshold) THEN
              water_table = water_table + dz
            ELSE
              layer_below_saturated = .FALSE.
              water_table = water_table + saturation * dz
            ENDIF
          ENDIF
        ENDIF
      ENDIF

    ENDDO ! isoil

    IF (local_depth > 0) THEN
      frac_liquid_top = frac_liquid_top / local_depth
    ELSE
      frac_liquid_top = 0._wp
    ENDIF

  END SUBROUTINE update_water_table

  !>
  !      Updates inundated area
  !
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE update_inun_frac(                   &
    & soil_depth_min,                            &
    & use_cti_basin, use_cti_flood, use_cti_max, &
    & use_wt_max, use_openwater,                 &
    & cti_min, cti_max, cti_min_flooding,        &
    & efold, wt_min, wt_max,                     &
    & soil_depth,                                &
    & cti_mean, cti_basin, cti_openwater,        &
    & cti_phi, cti_mu, cti_chi,                  &
    & mask, water_table, frac_liquid_water,      &
    & frac_inun, frac_inun_uncorr, wt_inun       &
  )
  !$acc routine seq

    ! Interface parameters
    !  Constants
    REAL(wp), INTENT(in)    :: soil_depth_min
    LOGICAL,  INTENT(in)    :: use_cti_basin
    LOGICAL,  INTENT(in)    :: use_cti_flood
    LOGICAL,  INTENT(in)    :: use_cti_max
    LOGICAL,  INTENT(in)    :: use_wt_max
    LOGICAL,  INTENT(in)    :: use_openwater
    REAL(wp), INTENT(in)    :: cti_min
    REAL(wp), INTENT(in)    :: cti_max
    REAL(wp), INTENT(in)    :: cti_min_flooding
    REAL(wp), INTENT(in)    :: efold
    REAL(wp), INTENT(in)    :: wt_min
    REAL(wp), INTENT(in)    :: wt_max
    REAL(wp), INTENT(in)    :: soil_depth
    REAL(wp), INTENT(in)    :: cti_mean
    REAL(wp), INTENT(in)    :: cti_basin
    REAL(wp), INTENT(in)    :: cti_openwater
    REAL(wp), INTENT(in)    :: cti_phi
    REAL(wp), INTENT(in)    :: cti_mu
    REAL(wp), INTENT(in)    :: cti_chi
    !  Topmodel variables
    LOGICAL,  INTENT(in)    :: mask
    REAL(wp), INTENT(in)    :: water_table
    REAL(wp), INTENT(in)    :: frac_liquid_water
    REAL(wp), INTENT(inout) :: frac_inun
    REAL(wp), INTENT(inout) :: frac_inun_uncorr
    REAL(wp), INTENT(inout) :: wt_inun

    ! Parameters
    REAL(wp), PARAMETER :: water_table_threshold = 0.95_wp

    ! Locals
    REAL(wp) :: cti_min_eff, cti_max_eff, cti_ref, cti_loc, cti_tmp
    REAL(wp) :: frac_wetland, wt_inun_loc
    REAL(wp) :: p_loc, p_min, p_max

    frac_inun        = 0._wp
    frac_inun_uncorr = 0._wp
    wt_inun          = 0._wp

    IF ((cti_mean > cti_min) .AND. (water_table > -water_table_threshold * soil_depth) .AND. mask) THEN

      ! Determine inundated area
      cti_ref = cti_mean
      IF (use_cti_basin) THEN
        cti_ref = cti_basin
      ENDIF

      cti_min_eff = efold * (wt_min - water_table) + cti_ref
      cti_max_eff = cti_max
      IF (use_wt_max) THEN
        cti_max_eff = efold * (wt_max - water_table) + cti_ref
      ENDIF

      IF (use_cti_flood) THEN
        cti_min_eff = MAX(cti_min_eff,cti_min_flooding)
      ENDIF

      p_min = 1._wp
      p_max = 1._wp
      IF (cti_phi /= cti_phi) THEN ! IsNaN(cti_phi)
        frac_wetland = 0._wp
      ELSE
        IF (cti_min_eff >= cti_mu) THEN
          p_min = gamm_p(cti_phi, (cti_min_eff-cti_mu)/cti_chi)
          IF (use_wt_max .OR. use_cti_max) THEN
            p_max = gamm_p(cti_phi, (cti_max_eff-cti_mu)/cti_chi)
          ENDIF
        ENDIF
      ENDIF

      IF (use_openwater) THEN
        p_max = 1._wp - cti_openwater
      ENDIF

      frac_wetland = MAX(0._wp,p_max - p_min)

      wt_inun_loc = 0._wp
      IF (frac_wetland > 0._wp) THEN
        p_loc = p_min + 0.5_wp * frac_wetland
        cti_tmp = gamm_inv_p(p_loc, cti_phi)
        cti_loc = cti_tmp * cti_chi + cti_mu
        wt_inun_loc = water_table + 1._wp/efold * (cti_loc - cti_ref)
      ENDIF

      IF (soil_depth > soil_depth_min) THEN
        frac_inun        = frac_wetland * frac_liquid_water
        frac_inun_uncorr = frac_wetland
        wt_inun          = wt_inun_loc
      ENDIF

    ENDIF

  END SUBROUTINE update_inun_frac

  !>
  !      Updates necessary averages for the topmodel
  !
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE update_topmodel_runmean(                       &
    & lexp, lfirst_year, lyear, lmonth, lday,               &
    & mlen, ylen, tspd,                                     &
    & summer_state,                                         &
    & nyear_lt,                                             &
    & larctic,                                              &
    & frac_inun, frac_inun_sum_month, frac_inun_mean_month, &
    & frac_inun_min_year, frac_inun_min_lt,                 &
    & frac_inun_summ_year, frac_inun_summ_lt,               &
    & frac_permanent_wetland,                               &
    & wt_inun, wt_inun_sum_month, wt_inun_mean_month,       &
    & wt_inun_min_year, wt_inun_min_lt                      &
    & )
  !$acc routine seq

    USE mo_top_util, ONLY: time_relaxation

    LOGICAL,  INTENT(in)    :: lexp, lyear, lmonth, lday ! .TRUE. if first step of a(n) Experiment/Year/Month/Day
    LOGICAL,  INTENT(in)    :: lfirst_year                     ! .TRUE. if first year of experiment (finished)
    INTEGER,  INTENT(in)    :: mlen                      ! Length of month of previous timestep [days] (by lmonth only)
    INTEGER,  INTENT(in)    :: ylen                      ! Length of year of previous timestep [days] (by lyear only)
    INTEGER,  INTENT(in)    :: tspd                      ! Time steps per day
    INTEGER,  INTENT(in)    :: summer_state              ! Length of summer [months] (>0) when avg. is to be performed
                                                         ! 0: summer month, -1 : non-summer month
    REAL(wp), INTENT(in)    :: nyear_lt                  ! efolding time for long term averages [years]
    LOGICAL,  INTENT(in)    :: larctic                   ! This is an "Arctic" point: i.e. lat > lat_const_wetland_type
    ! Various sums/averages/mimimums derived from the instant inundated fraction. All unit []
    REAL(wp), INTENT(in)    :: frac_inun                 ! Instant inundated fraction
    REAL(wp), INTENT(inout) :: frac_inun_sum_month       ! Intramonthly sum
    REAL(wp), INTENT(inout) :: frac_inun_mean_month      ! Monthly average
    REAL(wp), INTENT(inout) :: frac_inun_min_year        ! Minimum of monthly averages in current year
    REAL(wp), INTENT(inout) :: frac_inun_min_lt          ! Multiyear average of minimum monthly average with each year
    REAL(wp), INTENT(inout) :: frac_inun_summ_year       ! Sum/average of summer months monthly averages
    REAL(wp), INTENT(inout) :: frac_inun_summ_lt         ! Multiyear average of summer average of monthly averages
    REAL(wp), INTENT(inout) :: frac_permanent_wetland    ! Fraction considered permanently inundated
    ! Various sum/averages/minimums derived from the instant water table of the inundated fraction. All unit [m]
    REAL(wp), INTENT(in)    :: wt_inun                   ! Instant water table of the inundated fraction
    REAL(wp), INTENT(inout) :: wt_inun_sum_month         ! Intramonthly sum
    REAL(wp), INTENT(inout) :: wt_inun_mean_month        ! Monthly average
    REAL(wp), INTENT(inout) :: wt_inun_min_year          ! Minimum of monthly averages in current year
    REAL(wp), INTENT(inout) :: wt_inun_min_lt            ! Multiyear average of minimum of monthly averages

    IF (lday) THEN

      IF (lmonth) THEN

        IF (lyear) THEN

          ! Reset yearly values
          frac_inun_summ_year =    0._wp
          frac_inun_min_year  = 1000._wp
          wt_inun_min_year    = 1000._wp

        ENDIF ! lyear

        ! frac_inun: Monthly average of inundated fraction
        frac_inun_mean_month = frac_inun_sum_month / (mlen * tspd)
        frac_inun_sum_month  = 0._wp

        ! wt_inun: Monthly average water table height of inundated fraction
        wt_inun_mean_month = wt_inun_sum_month / (mlen * tspd)
        wt_inun_sum_month  = 0._wp

        ! frac_inun_min_year: Minimum monthly inundated fraction within a year
        frac_inun_min_year = MIN(frac_inun_min_year,frac_inun_mean_month)

        ! wt_inun_min_of_year: Minimum monthly water table in inundated fraction
        wt_inun_min_year = MIN(wt_inun_min_year,wt_inun_mean_month)

        ! frac_inun_summ_year
        IF (summer_state >= 0) THEN
          frac_inun_summ_year = frac_inun_summ_year + frac_inun_mean_month
        ENDIF
        IF (summer_state > 0) THEN
          frac_inun_summ_year = frac_inun_summ_year / summer_state
        ENDIF

        IF (lyear) THEN

          ! Update long term means
          CALL time_relaxation(  wt_inun_min_year,   wt_inun_min_lt, nstep=nyear_lt,linit=lexp)
          CALL time_relaxation(frac_inun_min_year, frac_inun_min_lt ,nstep=nyear_lt,linit=lexp)
          CALL time_relaxation(frac_inun_summ_year,frac_inun_summ_lt,nstep=nyear_lt,linit=lexp)

          ! Permanent wetland
          IF (larctic) THEN
            frac_permanent_wetland = frac_inun_summ_lt
          ELSE
            frac_permanent_wetland = frac_inun_min_lt
          ENDIF

        ENDIF ! lyear

      ENDIF ! lmonth

    ENDIF ! lday

    ! Any time step: Just add to sums
    frac_inun_sum_month = frac_inun_sum_month + frac_inun
    wt_inun_sum_month   = wt_inun_sum_month + wt_inun

  END SUBROUTINE update_topmodel_runmean

  !>
  !      Updates necessary averages for peatland calculations
  !
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE update_peat_runmean(                           &
    & lexp, lfirst_year, lyear, lmonth, lday,               &
    & mlen, ylen, tspd,                                     &
    & nyear_lt, nstep_relaxation,                           &
    & rain, snow,                                           &
    & precip_sum_day, precip_mean_day, precip_relax,        &
    & aet,                                                  &
    & aet_sum_day,  aet_mean_day, aet_relax,                &
    & poaet_mean_day,                                       &
    & poaet_sum_year, poaet_mean_year,                      &
    & poaet_lt                                              &
  )
  !$acc routine seq
    LOGICAL,  INTENT(in)    :: lexp, lyear, lmonth, lday       ! .TRUE. if first step of a(n) Experiment/Year/Month/Day
    LOGICAL,  INTENT(in)    :: lfirst_year                     ! .TRUE. if first year of experiment (finished)
    INTEGER,  INTENT(in)    :: mlen                            ! Length of month of previous timestep [days] (by lmonth only)
    INTEGER,  INTENT(in)    :: ylen                            ! Length of year of previous timestep [days] (by lyear only)
    INTEGER,  INTENT(in)    :: tspd                            ! Time steps per day
    REAL(wp), INTENT(in)    :: nyear_lt                        ! efolding time for long term averages [years]
    REAL(wp), INTENT(in)    :: nstep_relaxation                ! #efolding steps for precipitation and aet averages
    REAL(wp), INTENT(in)    :: rain, snow                      ! Rain and snow [kg m-2 s-1] (from A2L_)
    REAL(wp), INTENT(inout) :: precip_sum_day, precip_mean_day ! Daily sum/mean of total precipitation [kg m-2 s-1]
    REAL(wp), INTENT(inout) :: precip_relax                    ! Efolded average of precipitation [kg m-2 s-1]
    REAL(wp), INTENT(in)    :: aet                             ! Potential evapotranspiration [kg m-2 s-1]
    REAL(wp), INTENT(inout) :: aet_sum_day,  aet_mean_day      ! Daily sum/mean of potential evapotranspiration [kg m-2 s-1]
    REAL(wp), INTENT(inout) :: aet_relax                       ! Efolded average of aet [kg m-2 s-1]
    REAL(wp), INTENT(inout) :: poaet_mean_day                  ! Daily mean of precip/potential evapotranspiration []
    REAL(wp), INTENT(inout) :: poaet_sum_year, poaet_mean_year ! Yearly sum/mean of daily mean precip/pot evapotranspiration []
    REAL(wp), INTENT(inout) :: poaet_lt                        ! Long term mean of year mean of poaet []

    IF (lday) THEN
      precip_mean_day = precip_sum_day / tspd
      precip_sum_day  = 0._wp

      aet_mean_day = aet_sum_day / tspd
      aet_sum_day  = 0._wp

      ! Set threshold: minimum aet_mean allowed is 0.01 -> to fix problems at very cold temperatures
      poaet_mean_day = precip_mean_day / MAX(aet_mean_day,1.e-8_wp)

      poaet_sum_year = poaet_sum_year + poaet_mean_day

      IF (lmonth) THEN

        IF (lyear) THEN

          poaet_mean_year = poaet_sum_year / ylen
          poaet_sum_year  = 0._wp

          ! Sensible initial values
          IF (lfirst_year) THEN
            poaet_lt = poaet_mean_year
          ELSE
            poaet_lt = (poaet_lt * (nyear_lt - 1._wp) + poaet_mean_year) / nyear_lt
          ENDIF

        ENDIF ! lyear
      ENDIF ! lmonth
    ENDIF ! lday

    ! Any time step: Add to sums
    precip_sum_day = precip_sum_day + rain + snow
    aet_sum_day    = aet_sum_day + aet

    ! Update efolded averages
    aet_relax    = (   aet_relax * (nstep_relaxation - 1._wp) + aet        ) / nstep_relaxation
    precip_relax = (precip_relax * (nstep_relaxation - 1._wp) + rain + snow) / nstep_relaxation

  END SUBROUTINE update_peat_runmean

  !>
  !      Updates inundated fraction and it's watertable from climatology
  !
  !>
  !>
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE get_climatological_inun(      &
    & frac_inun_clim,                      & ! Input
    & frac_inun, wt_inun                   & ! Output
  )
  !$acc routine seq

    ! Interface parameters
    REAL(wp), INTENT(in)  :: frac_inun_clim            ! Climatological frac_inun
    REAL(wp), INTENT(out) :: frac_inun                 ! Actual frac_inun
    REAL(wp), INTENT(out) :: wt_inun                   ! Actual wt_inun

    frac_inun = frac_inun_clim
    IF (frac_inun_clim > 0._wp) THEN
      wt_inun = 0.25_wp
    ELSE
      wt_inun = 0.00_wp
    ENDIF

  END SUBROUTINE get_climatological_inun

#endif
END MODULE mo_top_process
