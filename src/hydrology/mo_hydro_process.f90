!>
!! @brief Contains the routines for the hydro processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish, message_text, em_warn

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_surface_hydrology_land, calc_surface_hydrology_glacier, calc_soil_hydrology, &
    & calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare,         &
    & get_canopy_conductance, get_water_stress_factor

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  REAL(wp), PARAMETER ::  zwdmin = 0.05_wp
  REAL(wp), PARAMETER ::  zwdcri = 0.9_wp
  REAL(wp), PARAMETER ::  zdrmin = 0.001_wp / (3600._wp*1000._wp)  ! factor for slow drainage (~ 10 hours)
  REAL(wp), PARAMETER ::  zdrmax = 0.1_wp   / (3600._wp*1000._wp)  ! factor for fast drainage (~ 4 days)
  REAL(wp), PARAMETER ::  zdrexp = 1.5_wp
  ! todo: get rid of zfcmin. The parameter seems to be needed only for numerical reasons.
  REAL(wp), PARAMETER ::  zfcmin = 1.e-10_wp              ! minimum field capacity (or rather epsilon?)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! @brief   Compute surface hydrology on glacier-free land
  !!
  !! @par Revision History
  !!   Routine based on the former JSBACH3 routine update_surf_down
  !!   Rewritten in Mai 2016 by V. Gayler
  !!   Adapted for JSBACH4 in February 2017 by R. Schnur
  !!
  !! @param [in]     lstart                T: Start of experiment
  !! @param [in]     dtime                 Time step
  !! @param [in]     l_dynsnow             T: Compute density of snow dynamically
  !! @param [in]     snow_depth_max        Snow depth limit [m water equivalent]; -1. for no limit
  !! @param [in]     t_unfilt              Surface tempemperature [K] (unfiltered)
  !! @param [in]     wind_10m              wind speed at 10m height [m/s]
  !! @param [in]     t_air                 lowest layer atmosphere temperature [K]
  !! @param [in]     skinres_canopy_max    Capacity of the canopy skin reservoir (also used to limit snow on canopy) [m]
  !! @param [in]     skinres_max           Total capacity of the skin reservoirs, i.e. soil and canopy [m]
  !! @param [in]     fract_snow            surface snow fraction []
  !! @param [in]     fract_wet             skin reservoir fraction []
  !! @param [in]     hcap_grnd             heat capacity of the uppermost soil layer [J m-2 K-1]
  !! @param [in]     evapotrans            evapotranspiration (including sublimation) [kg m-2 s-1]
  !! @param [in]     evapopot              potential evaporation (if there was enough water/ice) [kg m-2 s-1]
  !! @param [in]     transpiration         transpiration [kg m-2 s-1]
  !! @param [in]     rain                  liquid precipitation [kg m-2 s-1]
  !! @param [in]     snow                  solid precipitation [kg m-2 s-1]
  !! @param [in,out] wtr_skin              water content of the skin reservoir (vegetation and bare soil) [m]
  !! @param [in,out] weq_snow_soil         snow depth at the ground [m water equivalent]
  !! @param [in,out] weq_snow_can          snow depth on the canopy [m water equivalent]
  !! @param [in,out] snow_soil_dens        snow density on soil at non-glacier points [m water equivalent]
  !! @param [out]    q_snocpymlt           Heating due to snow melt on canopy [W m-2]
  !! @param [out]    snow_accum            snow budget change within time step [m water equivalent]
  !! @param [out]    snowmelt_soil         snow/ice melt at land points (excluding canopy) [kg m-2 s-1]
  !! @param [out]    evapotrans_soil       evapotranspiration from soil w/o skin and snow reservoirs [kg m-2 s-1]
  !! @param [out]    evapo_skin             evaporation from skin reservoir [kg m-2 s-1]
  !! @param [out]    evapo_snow             evaporation from snow [kg m-2 s-1]
  !! @param [out]    water_to_soil         Water available for infiltration into the soil [m water equivalent]
  !! @param [out]    evapo_deficit          Evaporation deficit flux due to inconsistent treatment of snow evap. [m]
  !
  PURE ELEMENTAL SUBROUTINE calc_surface_hydrology_land (           &
    & lstart, dtime, l_dynsnow, snow_depth_max,                     &
    & t_unfilt, wind_10m, t_air,                                    &
    & skinres_canopy_max, skinres_max, fract_snow, fract_wet,       &
    & hcap_grnd, evapotrans, evapopot, transpiration, rain, snow,   &
    & wtr_skin, weq_snow_soil, weq_snow_can, snow_soil_dens,        &
    & q_snocpymlt, snow_accum, snowmelt_soil,                       &
    & evapotrans_soil, evapo_skin, evapo_snow, water_to_soil, evapo_deficit)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf, dens_snow_min, dens_snow_max, dens_snow
    USE mo_hydro_constants,        ONLY: InterceptionEfficiency
    USE mo_sse_constants,          ONLY: snow_depth_min

    LOGICAL,  INTENT(in) :: lstart, l_dynsnow
    REAL(wp), INTENT(in) :: dtime, snow_depth_max
    !$claw model-data
    REAL(wp), INTENT(in) ::                                                         &
      & t_unfilt, wind_10m, t_air, skinres_canopy_max, skinres_max,                 &
      & fract_snow, fract_wet, hcap_grnd, evapotrans, evapopot, transpiration, rain, snow
    REAL(wp), INTENT(inout) :: &
      & wtr_skin, weq_snow_soil, weq_snow_can, snow_soil_dens
    REAL(wp), INTENT(out)   ::                   &
      & q_snocpymlt, snow_accum, snowmelt_soil,  &
      & evapotrans_soil, evapo_skin, evapo_snow, water_to_soil, evapo_deficit
    !$claw end model-data
    !
    !  local variables
    !
    REAL(wp) ::               &
      & rain_in_m,            & !< amount of rainfall within time step [m]
      & new_snow,             & !< amount of snowfall within time step [m]
      & weq_snow_soil_old,    & !< amount of snow prior to update [m]
      & evapotrans_in_m,      & !< amount of evapotranspiration within time step [m]
      & evapotrans_soil_in_m, & !< amount of evapotranspiration from soil (w/o snow and skin) within time step [m]
      & evapo_soil_in_m,      & !< amount of ground evaporation within time step [m]
      & evapo_skin_in_m,      & !< amount of evaporation from skin reservoir within time step [m]
      & transpiration_in_m,   & !< amount of transpiration within time step [m]
      & evapo_snow_in_m,      & !< amount of sublimation from snow within time step [m]
      & snowmelt_can,         & !< snow melt from canopy [m]
      & snowmelt_pot_in_m,    & !< potential amount of snow melt [m]
      & new_snow_can,         & !< amount of snowfall on canopy within time step [m]
      & new_snow_soil,        & !< amount of snowfall to soil within time step [m]
      & exp_t, exp_w,         & !< exponents needed for unloading of snow due to temperature / wind
      & snow_blown,           & !< amount of snow blown from canopy to the ground
      & canopy_pre_snow,      & !< snow depth on the canopy prior to its update
      & evapotrans_no_snow,   & !< evapotranspiration without snow evaporation [m]
      & evapo_snow_pot_in_m,  & !< potential snow evaporation
      & evapo_skin_pot_in_m,  & !< potential evaporation from wet surface (skin reservoir)
      & rain_to_skinres,      & !< amount of rain going into the skin reservoir
      & wtr_skin_pre_rain       !< amount of water in the skin reservoir prior to the update [m]


    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp ! [K]
    REAL(wp), PARAMETER :: zc2=1.87E5_wp       ! [Ks]
    REAL(wp), PARAMETER :: zc3=1.56E5_wp       ! [m]

    !$claw sca

    snow_accum = 0._wp
    weq_snow_soil_old = weq_snow_soil
    snowmelt_soil = 0._wp
    evapotrans_soil = 0._wp
    evapo_skin = 0._wp
    evapo_snow = 0._wp
    water_to_soil = 0._wp
    evapo_deficit = 0._wp

    !-----------------------------------------------------------------
    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m           = rain          * dtime/rhoh2o
    new_snow            = snow          * dtime/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m     = evapotrans    * dtime/rhoh2o
    transpiration_in_m  = transpiration * dtime/rhoh2o
    evapo_snow_pot_in_m = fract_snow    * evapopot * dtime/rhoh2o
    evapotrans_no_snow  = evapotrans_in_m - evapo_snow_pot_in_m
    evapo_skin_pot_in_m = (1._wp - fract_snow) * fract_wet * evapopot * dtime/rhoh2o

    !---------------------------------------------
    !  Budgets of snow (canopy, ground)
    !---------------------------------------------

    ! amount of snow remaining on the canopy
    new_snow_can = MIN(new_snow * InterceptionEfficiency, MAX(skinres_canopy_max - weq_snow_can, 0._wp))

    ! remaining snow falls on the ground
    new_snow_soil = new_snow - new_snow_can

    ! update snow on the canopy
    ! Note: evaporation happens from canopy, as long as there is enough snow
    canopy_pre_snow = weq_snow_can
    weq_snow_can = MIN(MAX(0._wp, canopy_pre_snow + new_snow_can + evapo_snow_pot_in_m), skinres_canopy_max)
    evapo_snow_in_m = evapo_snow_pot_in_m - (weq_snow_can - canopy_pre_snow - new_snow_can)
    ! unloading of snow from the canopy due to melting
    exp_t = MAX(0._wp, t_air-zc1)/zc2 * dtime
    snowmelt_can = weq_snow_can * (1._wp-EXP(-exp_t))
    weq_snow_can = weq_snow_can - snowmelt_can

    ! unloading of snow from the canopy due to wind
    exp_w = wind_10m/zc3 * dtime
    snow_blown = weq_snow_can * (1._wp-EXP(-exp_w))
    weq_snow_can  = weq_snow_can  - snow_blown
    new_snow_soil = new_snow_soil + snow_blown

    ! Heating due to snow melt on canopy
    ! @Todo Why is heating only computed for canopy?
    q_snocpymlt = snowmelt_can * rhoh2o * alf / dtime

    !  Snowfall and sublimation
    !-----------------------------------------------
    weq_snow_soil = weq_snow_soil + new_snow_soil + evapo_snow_in_m

    ! Correction if there was too much snow evaporation
    IF (weq_snow_soil < 0._wp) THEN
      evapotrans_no_snow = evapotrans_no_snow + weq_snow_soil
      evapo_deficit = weq_snow_soil
      weq_snow_soil = 0._wp
    ELSE
      evapo_deficit = 0._wp
    END IF

    !  Snow melt
    !------------------------

    IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
      IF (t_unfilt > tmelt) THEN
        snowmelt_pot_in_m = hcap_grnd * (t_unfilt-tmelt)/(alf*rhoh2o)     !  potential snow and ice melt [m]
        IF (weq_snow_soil > 0._wp) THEN
          snowmelt_soil = MAX(MIN(snowmelt_pot_in_m, weq_snow_soil), 0._wp)   ! snow melt limited by actual snow depth [m]
          weq_snow_soil = weq_snow_soil - snowmelt_soil                       ! reduce snow depth according to melting
        ! surface temperature correction for snow/ice melt is applied in task "snowmelt_correction"
        END IF
      END IF
    END IF

    !  Snow budget and meltwater
    !-----------------------------------------------------
    wtr_skin = wtr_skin + snowmelt_can                                    ! Add melt water from canopy to skin reservoir
    water_to_soil = snowmelt_soil + MAX(0._wp, wtr_skin-skinres_max)      ! Excess water enters the soil, and
    wtr_skin = MIN(skinres_max, wtr_skin)                                 ! skin reservoir is limited to the maximum
    snow_accum = new_snow + evapo_snow_pot_in_m - snowmelt_soil - snowmelt_can

    IF (snow_depth_max >= 0._wp) THEN
      ! Limit snow depth to avoid grid cells with infinitely growing snow depths in cooler climates
      water_to_soil = water_to_soil + MAX(0._wp, weq_snow_soil-snow_depth_max)
      weq_snow_soil = MIN(snow_depth_max, weq_snow_soil)
    END IF

    !-----------------------------------------------------------
    !  Budget of water in skin reservoir (canopy and ground)
    !-----------------------------------------------------------
    ! Interception of rain
    rain_to_skinres = MIN(rain_in_m * InterceptionEfficiency, MAX(skinres_max - wtr_skin, 0._wp)) ! Part of rain goes to skin
                                                                                                  ! reservoir
    ! Note: the MAX in following line is only to account for apparent precision error (small negative values)
    water_to_soil = water_to_soil + MAX(rain_in_m - rain_to_skinres, 0._wp)        ! Rest of rain goes to soil
    wtr_skin_pre_rain = wtr_skin
    wtr_skin = MIN(MAX(0._wp, wtr_skin_pre_rain + rain_to_skinres + evapo_skin_pot_in_m), skinres_max)

    ! Note: at this point ignoring the amount of dew that exceeds the skin reservoir
    !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
    evapo_skin_in_m = wtr_skin - (wtr_skin_pre_rain + rain_to_skinres)
    evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m
    IF (evapo_skin_pot_in_m < 0._wp) THEN
      evapo_deficit = evapo_deficit + evapo_skin_pot_in_m - evapo_skin_in_m
    END IF
    evapo_soil_in_m = evapotrans_soil_in_m - transpiration_in_m  ! Evaporation from ground

    ! positive values of evaporation and transpiration (dew) are added to water_to_soil
    ! negative fluxes change soil moisture later in calc_soil_hydrology
    IF (evapo_soil_in_m > 0._wp) water_to_soil = water_to_soil + evapo_soil_in_m
    IF (transpiration_in_m > 0._wp) water_to_soil = water_to_soil + transpiration_in_m

    ! Compute snow density
    IF (l_dynsnow) THEN
      IF (weq_snow_soil_old * rhoh2o / snow_soil_dens > snow_depth_min) THEN
        ! Update density of old snow following eq. 31 of Verseghy, 1991
        snow_soil_dens = dens_snow_max - (dens_snow_max - snow_soil_dens) * EXP(-0.01_wp * dtime/3600._wp)
        ! Calculate weighted mean density from old and fresh snow
        IF (weq_snow_soil > weq_snow_soil_old) THEN
          snow_soil_dens =  (snow_soil_dens * weq_snow_soil_old + dens_snow_min * (weq_snow_soil - weq_snow_soil_old))  &
            &              / weq_snow_soil
        END IF
      ELSE
        snow_soil_dens = dens_snow_min  ! fresh snow
      END IF
    ELSE
      snow_soil_dens = dens_snow
    END IF

    ! Transform fluxes in output from [m] to [kg m-2 s-1]
    snowmelt_soil = snowmelt_soil * rhoh2o / dtime

    evapotrans_soil = evapotrans_soil_in_m * rhoh2o / dtime
    evapo_skin = evapo_skin_in_m  * rhoh2o / dtime
    evapo_snow = fract_snow * evapopot

  END SUBROUTINE calc_surface_hydrology_land

  ! ===============================================================================================================================
  !>
  !! @brief   calc surface hydrology on glaciers
  !!
  !! @par Revision History
  !!   Routine based on the former JSBACH3 routine update_surf_down
  !!   Rewritten in Mai 2016 by V. Gayler
  !!   Adapted for JSBACH4 in February 2017 by R. Schnur
  !!
  !! @param [in]     lstart                T: Start of experiment
  !! @param [in]     dtime                 Time step
  !! @param [in]     t_unfilt              surface tempemperature [K] (unfiltered)
  !! @param [in]     fract_snow            surface snow fraction []
  !! @param [in]     hcap_grnd             heat capacity of the uppermost soil layer [J m-2 K-1]
  !! @param [in]     evapotrans            evapotranspiration (including sublimation) [kg m-2 s-1]
  !! @param [in]     evapopot              potential evaporation [kg m-2 s-1]
  !! @param [in]     rain                  liquid precipitation [kg m-2 s-1]
  !! @param [in]     snow                  solid precipitation [kg m-2 s-1]
  !! @param [in,out] weq_glac              glacier depth (snow and ice) [m water equivalent]
  !! @param [out]    q_snocpymlt           Heating due to snow melt on canopy [W m-2]
  !! @param [out]    snowmelt              Snow/ice melt at glacier points [kg m-2 s-1]
  !! @param [out]    runoff_glac           glacier runoff (rain+snow/ice melt, no calving) [kg m-2 s-1]
  !! @param [out]    pme_glacier           precipitation minus sublimation on glacier [kg m-2 s-1]
  !
  PURE ELEMENTAL SUBROUTINE calc_surface_hydrology_glacier ( &
    & lstart, dtime,                                         &
    & t_unfilt,                                              &
    & fract_snow, hcap_grnd,                                 &
    & evapotrans, evapopot, rain, snow,                      &
    & weq_glac, q_snocpymlt,                                 &
    & snowmelt, runoff_glac,                                 &
    & pme_glacier)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf

    LOGICAL,  INTENT(in) :: lstart
    REAL(wp), INTENT(IN) :: dtime
    !$claw model-data
    REAL(wp), INTENT(in) ::       &
      & t_unfilt,                 &
      & fract_snow, hcap_grnd, evapotrans, evapopot, rain, snow
    REAL(wp), INTENT(inout) :: &
      & weq_glac
    REAL(wp), INTENT(out)   :: &
      & q_snocpymlt, snowmelt, runoff_glac, pme_glacier
    !$claw end model-data

    !
    !  local variables
    !
    REAL(wp) ::              &
      & rain_in_m,           & !< amount of rainfall within time step [m]
      & new_snow,            & !< amount of snowfall within time step [m]
      & evapotrans_in_m,     & !< amount of evapotranspiration within time step [m]
      & snowmelt_pot_in_m,   & !< potential amount of snow melt [m]
      & evapo_snow_pot_in_m, & !< potential snow evaporation
      & pme_glacier_in_m,    & !< P-E in m
      & runoff_glac_in_m,    & !< Runoff in m
      & snowmelt_in_m          !< Snow/ice melt in m

    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp !< [K]
    REAL(wp), PARAMETER :: zc2=1.87E5_wp       !< [Ks]
    REAL(wp), PARAMETER :: zc3=1.56E5_wp       !< [m]

    !$claw sca

    !----------------------------------------------------------------------------------------------

    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m          = rain       * dtime/rhoh2o
    new_snow           = snow       * dtime/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m    = evapotrans * dtime/rhoh2o
    evapo_snow_pot_in_m = fract_snow * evapopot * dtime/rhoh2o

    !  Snowfall and sublimation on glaciers
    !---------------------------------------

    weq_glac         = weq_glac  + new_snow + evapo_snow_pot_in_m   ! glacier depth [m]
    pme_glacier_in_m = rain_in_m + new_snow + evapotrans_in_m       ! P-E on glaciers [m]
    runoff_glac_in_m = rain_in_m                                    ! there is no infiltration on glaciers


    !  Snow and glacier melt
    !------------------------

    snowmelt_in_m = 0._wp
    IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
      IF (t_unfilt > tmelt) THEN
        snowmelt_pot_in_m = hcap_grnd * (t_unfilt-tmelt)/(alf*rhoh2o) !  potential snow and ice melt [m]
        snowmelt_in_m     = snowmelt_pot_in_m                       ! there is an unlimited amount of snow
        weq_glac          = weq_glac - snowmelt_in_m                ! reduce glacier depth according to melting
        runoff_glac_in_m  = runoff_glac_in_m + snowmelt_in_m        ! add melt water to the runoff
        ! surface temperature correction for snow/ice melt is applied in task "snowmelt_correction"
      END IF
    END IF

    q_snocpymlt = 0._wp

    snowmelt    = snowmelt_in_m    * rhoh2o / dtime
    runoff_glac = runoff_glac_in_m * rhoh2o / dtime
    pme_glacier = pme_glacier_in_m * rhoh2o / dtime

  END SUBROUTINE calc_surface_hydrology_glacier

  !
  ! @brief      Calculates runoff, drainage and soil moisture
  !
  ! @param [in]      nc                  Vector lenth
  ! @param [in]      nsoil               Number of soil layers
  ! @param [in]      nlat                (Effective) number of latitudes
  ! @param [in]      dtime               Time step [s]
  ! @param [in]      l_infil_subzero     Allow infiltration of supercolled water
  ! @param [in]      hydrotiles_type    Property defining whether the actual tile is source or sink for lateral flow
  ! @param [in]      soil_depth_sl       Thicknesses of soil layers until bedrocj [m]
  ! @param [in]      root_depth_sl       Thicknesses of soil layers until rooting depth [m]
  ! @param [in]      fract_org_sl        Fractions per soil layer of organic comonent []
  ! @param [in]      hyd_cond_sat_sl     Hydraulic conductivity of saturated soil [m/s]
  ! @param [in]      matrix_pot_sl       Soil matrix potential [m]
  ! @param [in]      bclapp_sl           Exponent B in Clapp and Hornberger
  ! @param [in]      pore_size_index_sl  Soil pore size distribution index used in Van Genuchten method
  ! @param [in]      vol_porosity_sl     Soil porosity                      - ice not considered [m/m]
  ! @param [in]      vol_field_cap_sl    Volumetric field capacity          - ice not considered [m/m]
  ! @param [in]      vol_p_wilt_sl       Volumetric permanent wilting point - ice not considered [m/m]
  ! @param [in]      ws_max_mineral      Maximum water holding capacity of non-organic soil [m]
  ! @param [in]      oro_stddev          Sub-grid standard deviation of orography [m]
  ! @param [in]      t_soil_sl1          Temperature of the upper soil layer [K]
  ! @param [in]      water_to_soil       Amount of water going into the soil [m]
  ! @param [in]      evapotrans_soil     Evapotranspiration from the soil (w/o snow or canopy) [kg m-2 s-1]
  ! @param [in]      transpiration       Evapotranspiration [kg m-2 s-1]
  ! @param [in,out]  wtr_soil_sl         Content un-frozen water of soil layers in rooting zone [m]
  ! @param [in,out]  ice_soil_sl         Content of frozen water of soil layers in rooting zone [m]
  ! @param [out]     wtr_rootzone        Liquid water content of the root zone [m]
  ! @param [out]     wsat_sl             Maximum possible amount of liquid water in soil layer
  !                                        - reduced by amount of ice [m]
  ! @param [out]     field_cap_sl        Liquid water content at field capacity in soil layer [m]
  ! @param [out]     p_wilt_sl           Water content at permanent wilting point in soil layers
  !                                        - reduced relative to ice fraction) [m]
  ! @param [out]     runoff              Runoff
  ! @param [out]     drainage            Drainage
  !
  SUBROUTINE calc_soil_hydrology(                                                        &
    ! in
    & nc, nsoil, nlat, dtime, ltpe_closed, ltpe_open,                                    &
    & l_infil_subzero, hydrotiles_type,                                                  &
    & soil_depth_sl, root_depth_sl, fract_org_sl,                                        &
    & hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, vol_porosity_sl,    &
    & vol_field_cap_sl, vol_p_wilt_sl, ws_max_mineral, oro_stddev, t_soil_sl1,           &
    & water_to_soil_tile, evapotrans_soil, transpiration,                                &
    ! inout
    & wtr_soil_sl, ice_soil_sl, tpe_overflow, w_pond,                                    &
    ! out
    & wtr_rootzone, wsat_sl, field_cap_sl, p_wilt_sl, runoff, drainage                   &
    & )

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o
    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max, &
      & UPLAND_, LOWLAND_
    USE mo_hydro_util,             ONLY: get_amount_in_rootzone

    INTEGER,  INTENT(in) :: &
      & nc, nsoil, nlat
    REAL(wp), INTENT(in) :: &
      & dtime
    LOGICAL, INTENT(in) :: ltpe_closed
    LOGICAL, INTENT(in) :: ltpe_open
    LOGICAL, INTENT(in) :: l_infil_subzero
    INTEGER, INTENT(in) :: hydrotiles_type
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & soil_depth_sl, root_depth_sl, fract_org_sl,    &
      & hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl
    REAL(wp), INTENT(in), DIMENSION(:) ::               &
      & ws_max_mineral, oro_stddev, t_soil_sl1, water_to_soil_tile, &
      & evapotrans_soil, transpiration
    REAL(wp), INTENT(inout), DIMENSION(:) :: &
      & wtr_rootzone, tpe_overflow, w_pond
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & wtr_soil_sl, ice_soil_sl
    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & runoff, drainage
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & wsat_sl, field_cap_sl, p_wilt_sl

    ! Local variables
    REAL(wp) ::                       &
      ! variables for infiltration calculation
      & sigma_0,                      & !< minimum value
      & sigma_max,                    & !< resolution-dependent maximum
      & steepness,                    & !< sub-grid scale steepness of the orography
      & zbws, zb1, zbm, zconw1, zvol, &
      ! Other variables
      & evapo_soil_in_m(nc),          & !< evaporation from the soil (without transpiration) [m]
      & transpiration_in_m(nc),       & !< transpiration within time step [m]
      & ws_rel(nc),                   & !< relative soil moisture (0: wilting point, 1: field capacity)
      & infilt(nc),                   & !< infiltration [m]
      & water_to_soil(nc)               !< water available at soil surface [m]

    REAL(wp) ::              &
      & ws_max(nc),          & ! water holding capacity of the soil [m]
      & ws_min_drain(nc),    & ! minimum soil moisture for drainage
      & ws_incl_ice(nc),     & ! soil moisture including ice [m]
      & ice_rootzone(nc)       ! amount of ice in the rootzone [m]
    
    REAL(wp) ::              &
      & hlp_2d(nc,nsoil)       ! helper array needed with openacc

    INTEGER :: jl, i

    INTEGER, PARAMETER :: ilog = 0 ! Switch for debugging output
    INTEGER :: jllog = 0


    !$claw acc data pcreate(evapo_soil_in_m,transpiration_in_m,infilt,ws_max,ws_min_drain,ws_incl_ice,ice_rootzone)
    !$claw acc data present(transpiration_in_m,transpiration,wsat_sl,vol_porosity_sl,soil_depth_sl,ice_soil_sl, &
    !$claw acc              field_cap_sl,vol_field_cap_sl,p_wilt_sl,vol_p_wilt_sl)
    !$claw acc parallel
    !$claw acc loop gang
    !
    !> Calculation of soil saturation, field capacity and wilting point taking into account soil ice
    !
    DO i = 1, nc
      !$claw acc loop seq
      DO jl = 1, nsoil
        ! Note: here we are calculating soil saturation (this does not correspond to weq_rootzone_max, which was
        !    calculated from field capacity)
        wsat_sl     (i,jl) = MAX(vol_porosity_sl (i,jl) * soil_depth_sl(i,jl) - ice_soil_sl(i,jl), 0._wp)  ! -> diffusion
        field_cap_sl(i,jl) = MAX(vol_field_cap_sl(i,jl) * soil_depth_sl(i,jl) - ice_soil_sl(i,jl), 0._wp)  ! -> percolation
        ! todo: to be discussed: Does it make sense, to reduce water content at wilting point relative to the ice fraction?
        IF (vol_porosity_sl(i,jl) > 0._wp) THEN
          p_wilt_sl(i,jl) = wsat_sl(i,jl) * vol_p_wilt_sl(i,jl) / vol_porosity_sl(i,jl)
        ELSE
          p_wilt_sl(i,jl) = 0._wp
        END IF
      END DO
    END DO
    !$claw acc end parallel
    !$claw acc end data

    ! Note: ws_max calculated here (based on field capacity) is different from wsat calculated
    !       below (based on saturation). ws_max does not account for ice.
    hlp_2d(:,:) = vol_field_cap_sl(:,:) * soil_depth_sl(:,:)
    CALL get_amount_in_rootzone(hlp_2d(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), ws_max(:))

    !$claw acc data present(fract_org_sl,ws_max,ws_max_mineral)
    !$claw acc parallel loop gang vector
    ! Careful: ws_max just calculated from vol_field_cap_sl (depending on organic matter) is
    !          overwritten with weq_rootzone_max from the initial file (not depending on organic matter)
    !          in grid cells with no significant organic fraction. This might lead to unexpected
    !          inconsistencies!
    DO i = 1, nc
      IF (SUM(fract_org_sl(i,:)) <= EPSILON(1._wp)) THEN
        ws_max(i) = ws_max_mineral(i)
      END IF
    END DO
    !$claw acc end parallel
    !$claw acc end data

    CALL get_amount_in_rootzone(ice_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    !$claw acc present(ws_min_drain,ws_max,oro_stddev,runoff,infilt)    &
    !$claw acc end parallel

    !>
    !> 1. Infiltration
    !>
    !
    ! f(w) = 1 - (1 - w/w_max)**b
    !
    ! b: shape parameter, defining the sub-gid scale steepness of the orographie:
    !
    ! b = (sigma_z - sigma_0) / (sigma_z + sigma_max)
    !
    ! sigma_z: standard deviation of topographic height
    ! sigma_0: minimum value (100 m); below b = 0.01
    ! sigma_max: resolution dependent maximum

    sigma_0   = oro_var_min
    sigma_max = oro_var_max * 64._wp / REAL(nlat, wp) ! @todo: site level

    !$claw acc parallel
    !$claw acc loop gang vector
    DO i=1,nc
      steepness = MAX(0._wp, oro_stddev(i)-sigma_0) / (oro_stddev(i)+sigma_max)

      ! todo: parameters should be defined as parameters
      zbws   = MAX(MIN(steepness,0.5_wp), 0.01_wp)
      zb1    = 1._wp+zbws
      zbm    = 1._wp/zb1
      zconw1 = ws_max(i)*zb1

      ! Add lateral flow from upland tiles
      IF (hydrotiles_type == LOWLAND_) THEN
          water_to_soil(i) = water_to_soil_tile(i) + 0.9_wp * w_pond(i)
          w_pond(i) = 0.1_wp * w_pond(i)
      ELSE
          water_to_soil(i) = water_to_soil_tile(i)
      END IF

      !>  1.1 Surface runoff and infiltration
      !   -----------------------------------
      ! todo: The amount of runoff and infiltration is calculated regarding the complete root zone,
      !    rather than the uppermost soil layer (or the soil until bedrock, as this has nothing to
      !    do with plants).
      !    Using the uppermost soil layer only might need sub-time-steps...

      ! No drainage is possible below a certain minimum amount of water (zwdmin). This corresponds
      ! to the layered bucket scheme, and has nothing to do with field capacity or wilting point.
      ws_min_drain(i) = zwdmin * ws_max(i)    ! minimum amount of root zone soil water for drainage

      IF (ltpe_closed) THEN
        ! Terra planet setup without runoff
        runoff(i) = 0._wp
        infilt(i) = water_to_soil(i)
      ELSE
        IF (.NOT. l_infil_subzero .AND. (t_soil_sl1(i) < tmelt)) THEN
          !pdv: Temperature dependence should not be included when using the 5-layer snow scheme because:
          !pdv  when tsurf == 0.C; T_soil < 0.C -- and all snow melt is going to surfcae runoff
          runoff(i) = water_to_soil(i)       !  no infiltration as uppermost soil layer is frozen -> runoff
          infilt(i) = 0._wp
        ELSE
          ws_incl_ice(i) = wtr_rootzone(i) + ice_rootzone(i)            ! soil water including ice
          IF (water_to_soil(i) > 0._wp .AND. ws_incl_ice(i) > ws_min_drain(i)) THEN ! soilwater content above limit for drainage
            ! compare jsbach3 documentation section 2.3.2.2
            ws_rel(i) = MIN(1._wp, ws_incl_ice(i) / ws_max(i))             ! relative soil moisture
            zvol   = (1._wp - ws_rel(i))**zbm - water_to_soil(i)/zconw1    ! factor taking account subgrid scale
                                                                           ! inhomogeneity and steepness

            runoff(i) = water_to_soil(i) - (ws_max(i)-ws_incl_ice(i)) ! if >0: amount of water exceeding water holding capacity
            IF (zvol > 0._wp) runoff(i) = runoff(i) + ws_max(i)*zvol**zb1
            runoff(i) = MAX(MIN(runoff(i), water_to_soil(i)), 0._wp)   ! no negative runoff, runoff limited by available water
            infilt(i) = water_to_soil(i) - runoff(i)                   ! water not going into runoff is infiltrated
          ELSE
            runoff(i) = 0._wp                                            ! no runoff as soil is too dry
            infilt(i) = water_to_soil(i)                                ! all water infiltrated
          END IF
        END IF
      END IF

      ! Use 50% of surface runoff as lateral flow towards lowland tiles
      IF (hydrotiles_type == UPLAND_) THEN
          w_pond(i) = w_pond(i) + 0.5_wp * runoff(i)
          runoff(i) = 0.5_wp * runoff(i) 
      END IF

    END DO

    !
    !> 1.2 Update soil water with respect to evapotransration fluxes
    !
    ! Note: Even though evapotranspiration was calculated earlier, the corresponding water is taken
    !       from soil only here. Moving digest_evaporation to the beginnig of this routine would
    !       lead to wetter soils, as more water could infiltrate.
    !$claw acc data present(evapo_soil_in_m,evapotrans_soil,transpiration_in_m,transpiration)
    !$claw acc parallel
    !$claw acc loop gang
    DO i = 1, nc
      transpiration_in_m(i) = transpiration(i)  * dtime / rhoh2o
      ! Note: at this point ignoring the amount of dew that exceeds the skin reservoir
      !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
      evapo_soil_in_m(i) = evapotrans_soil(i) * dtime / rhoh2o - transpiration_in_m(i)
    END DO

    CALL digest_evapotrans(nc, nsoil, root_depth_sl, soil_depth_sl,  &
      &     field_cap_sl, transpiration_in_m, evapo_soil_in_m, wtr_soil_sl, ice_soil_sl)

    !$claw acc end parallel
    !---------------------
    !> 2. drainage and runoff
    ! ----------------------------------------------------------------------------------

    !$claw acc parallel
    !$claw acc loop gang
    DO i=1,nc
      ! initilalization of drainage
      drainage(i) = 0._wp
    END DO
    !$claw acc end parallel

    !$claw acc end data

    ! For numerical reasons, half of the water is infiltrated
    ! before soilhyd and the other half of it afterwards

    !$claw acc data present(infilt)
    !$claw acc parallel loop gang vector
    DO i = 1, nc
      infilt(i) = infilt(i) / 2._wp
    END DO
    !$claw acc end data

    ! update soil moisture according to infiltration
    !    - first half of the water
    !$claw sca forward update savepoint=ser_jsb_hydro_calc_soil_hydrology_digest_infilt_0
    CALL digest_infilt(nc, nsoil, ltpe_closed, field_cap_sl, infilt, &
      &     wtr_soil_sl, drainage, tpe_overflow)

    ! Drainage calculation with multi-layer soil scheme
    CALL soilhyd(nc, nsoil, dtime, ilog, jllog, wtr_soil_sl, soil_depth_sl,     &
      &     wsat_sl, field_cap_sl, hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, &
      &     matrix_pot_sl, drainage, pore_size_index_sl, p_wilt_sl,             &
      &     ice_soil_sl, tpe_overflow, ltpe_closed, ltpe_open)

    ! update soil moisture according to infiltration
    !    - second half of the water
    !$claw sca forward update savepoint=ser_jsb_hydro_calc_soil_hydrology_digest_infilt_1
    CALL digest_infilt(nc, nsoil, ltpe_closed, field_cap_sl, infilt, &
         wtr_soil_sl, drainage, tpe_overflow)

    ! recalculate soil moisture from the soil moisture of the individual levels
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    !$claw acc data present(runoff,drainage)
    !$claw acc parallel loop gang vector
    DO i = 1, nc
      runoff(i)   = runoff  (i) * rhoh2o / dtime
      drainage(i) = drainage(i) * rhoh2o / dtime
    END DO
    !$claw acc end data
    !$claw acc end data
  END SUBROUTINE calc_soil_hydrology

  !------------------------------------------------------------------------------------------------
  SUBROUTINE soilhyd(nc, nsoil, dtime, ilog, jllog, wtr_soil_sl, dsoil,                   &
    & wsat_sl, field_cap_sl, hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, matrix_pot_sl,  &
    & drainage, pore_size_index_sl, p_wilt_sl, ice_soil_sl, tpe_overflow, ltpe_closed, ltpe_open)

    !----------------------------------------------------------------------------------------------
    !
    !  SOIL HYDROLOGY - Routine to be build into REMO to calculate soil hydrology for 5 soil layers
    !
    !  Calculation of percolation and vertical diffusion of water
    !
    !  Programming and development: Stefan Hagemann (MPI-M)
    !----------------------------------------------------------------------------------------------

    INTEGER,  INTENT(IN)            :: nc                         ! vector length
    INTEGER,  INTENT(IN)            :: nsoil                      ! number of soil layers
    REAL(wp), INTENT(IN)            :: dtime                      ! model time step legth [s]
    INTEGER,  INTENT(IN)            :: ilog                       !
    INTEGER,  INTENT(IN)            :: jllog                      !
    REAL(wp), INTENT(IN)            :: dsoil(:,:)                 ! soil depth until bedrock within each layer [m]
    REAL(wp), INTENT(IN)            :: wsat_sl(:,:)               ! Maximum storage in each soil layer (porosity)
                                                                  !    (reduced by ice) [m]
    REAL(wp), INTENT(IN)            :: field_cap_sl(:,:)          ! field capacity of the soil layer (reduced by ice) [m]
    REAL(wp), INTENT(IN)            :: hyd_cond_sat_sl(:,:)       ! hydraulic conductivity of saturated soil [m/s]
    REAL(wp), INTENT(IN)            :: vol_porosity_sl(:,:)       ! volumetric soil porosity [m/m]
    REAL(wp), INTENT(IN)            :: bclapp_sl(:,:)             ! exponent B in Clapp and Hornberger
    REAL(wp), INTENT(IN)            :: matrix_pot_sl(:,:)         ! matrix potential [m]
    REAL(wp), INTENT(IN)            :: pore_size_index_sl(:,:)    ! soil pore size distribution index used in Van Genuchten
    REAL(wp), INTENT(IN)            :: p_wilt_sl(:,:)             ! wilting point - reduced relative to soil ice fraction
                                                                  !    todo: discuss if ice reduction is wanted here
    REAL(wp), INTENT(INOUT)         :: wtr_soil_sl(:,:)           ! soil moisture of each layer [m]
    REAL(wp), INTENT(INOUT)         :: drainage(:)                ! amount of drainage within timestep [m]
    REAL(wp), INTENT(INOUT)         :: tpe_overflow(:)
    REAL(wp), INTENT(IN), OPTIONAL  :: ice_soil_sl(:,:)           ! amount of ice in the soil [m]

    LOGICAL, INTENT(IN) :: ltpe_closed
    LOGICAL, INTENT(IN) :: ltpe_open

    INTEGER  :: i, j                        ! looping index
    REAL(wp) :: percolation(nc,0:nsoil)     ! downward flux from one soil layer to the layer below [m/s]
                                            !     and later amount of perculation per time step [m]
    REAL(wp) :: diffus_sl(nc,nsoil)         ! diffusivity at mid-depth of the layer
    REAL(wp) :: ws_min_drain(nc)            ! minimum soil moisture for drainage
    REAL(wp) :: ws_rel(nc,nsoil)            ! normalized (~relative) soil mositure
    REAL(wp) :: ws_crit(nc)                 ! soil moisture above which drainage strongly increases [m]
    REAL(wp) :: drain_crit(nc)              ! drainage above a soil moisture of ws_crit [m/s]
    REAL(wp) :: ws_tmp(nc)                  ! temporary array for soil moisture
    REAL(wp) :: drain_sl(nc)                ! drainage from the current layer
    LOGICAL  :: soil_above(nc)              ! is there a soil layer above?
    LOGICAL  :: soil_below(nc)              ! does the layer below have a soil fraction?
    REAL(wp) :: zda(nc,nsoil)               ! diffusion coefficients
    REAL(wp) :: zdb(nc,nsoil)               !     ''
    REAL(wp) :: zdc(nc,nsoil)               !     ''
    REAL(wp) :: ztri(nc,nsoil)              !     ''
    REAL(wp) :: diffusivity(nc,nsoil-1)     ! diffusivity between two layers, i.e. at the bottom of each layer
    REAL(wp) :: zdiflog(10), wslog(nsoil)   ! needed for debugging output
    REAL(wp) :: remoist(nc)                 ! for TPE

    REAL(wp) :: zmvg_sl(nc,nsoil)           ! exponent M in Van Genuchten scheme = pore_size_index / (pore_size_index + 1)


    !$claw acc data &
    !$claw acc present(dsoil, wsat_sl, field_cap_sl, hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, matrix_pot_sl, &
    !$claw acc         pore_size_index_sl, p_wilt_sl, wtr_soil_sl, drainage, tpe_overflow, ice_soil_sl)          &
    !$claw acc pcreate(diffus_sl, ws_min_drain, ws_rel, ws_crit, drain_crit, ws_tmp, drain_sl, soil_above,       &
    !$claw acc         soil_below, zda, zdb, zdc, ztri, diffusivity, zkredu, zdiflog, remoist, zmvg_sl,          &
    !$claw acc         percolation, wslog)

    !$claw acc parallel loop gang vector collapse(2)
    DO j=1,nc
      DO i=1,nsoil
        ws_rel(j,i) = 0._wp
      END DO
    END DO

    !$claw acc parallel loop gang vector collapse(2)
    DO j = 1, nc
      DO i = 1, nsoil
        zmvg_sl(j,i) = pore_size_index_sl(j,i) / (pore_size_index_sl(j,i)+1._wp)
      END DO
    END DO
    !----------------------------------------------------------------------------------------------
    ! Percolation = Gravitational drainage
    !
    ! de Vrese: percolation scheme using Van Genuchten-formulation and Euler upstream
    ! ----------------------------------------------------------------------------------------------
    ! Van Genuchten-formulation
    !
    ! percolation
    !     = relative hydraulic conductivity (K_r) * hydraulic conductivity at staturation (hydr_cond_sat)
    !
    !    K_r = SQRT(THETA) * (1-(1-THETA**1/m)**m)**2
    !
    ! with
    !  THETA = (theta - theta_r)/(theta_s - theta_r)
    !
    !  THETA:   normalised soil moisture     -> ws_rel
    !  theta:   current soil moisture        -> wtr_soil_sl
    !  theta_s: soil moisture at stauration  -> field_cap
    !  theta_r: residual soil moisture       -> p_wilt

    !$claw acc parallel default(none)
    !$claw acc loop gang vector
    DO j = 1, nc

      ! Percolation of the lowest soil level is zero
      percolation(j,nsoil) = 0._wp

      ! Percolation of the remaining soil levels
      !$claw acc loop seq
      DO i = 1, nsoil-1
        IF (field_cap_sl(j,i) > zfcmin .AND. field_cap_sl(j,i) > p_wilt_sl(j,i)) THEN
          ws_rel(j,i) = MAX(wtr_soil_sl(j,i) - p_wilt_sl(j,i), 0._wp) / (field_cap_sl(j,i) - p_wilt_sl(j,i))
          ws_rel(j,i) = MIN(ws_rel(j,i), 1._wp)
          percolation(j,i) =  hyd_cond_sat_sl(j,i) * SQRT(ws_rel(j,i))                                      &
            &                  * (1._wp - (1._wp - ws_rel(j,i)**(1._wp/zmvg_sl(j,i)))**zmvg_sl(j,i))**2._wp &
            &                  * dtime     ! [m/s] -> [m]
        ELSE
          percolation(j,i) = 0._wp
        END IF
      END DO  ! soil levels (until nsoil-1)

      !
      !> Drainage
      !
      ! Drainage of the lowest layer above bedrock. It replaces percolation calculated above
      !$claw acc loop seq
      DO i = 1, nsoil

        ! is there a layer above or a layer with soil below?
        soil_above(j) = (i > 1)
        IF (i == nsoil) THEN
          soil_below(j) = .false.
        ELSE
          soil_below(j) = dsoil(j,i+1) > 0._wp
        END IF

        ! Note: field_cap_sl must not be mixed up with field capacity in a hydrological sense. It is
        !    rather the maximum filling of a multi level bucket scheme.

        ! Above a critical point, fast drainage occurs.
        ws_crit(j) = field_cap_sl(j,i) * zwdcri     ! soil moisture above which drainage strongly increases

        ! calculate drainage of the lowest layer above bedrock
        IF (dsoil(j,i) > 0._wp .AND. .NOT. soil_below(j)) THEN       ! lowest layer above bedrock

          IF (wtr_soil_sl(j,i) <= p_wilt_sl(j,i)) THEN
            drain_sl(j) = 0._wp              ! no drainage without mobile water
          ELSE

            ! Slow drainage of mobile water

            IF (field_cap_sl(j,i) > p_wilt_sl(j,i)) THEN
              ! slow drainage = slow drainage factor * relative soil moisture
              drain_sl(j) = zdrmin * (wtr_soil_sl(j,i) - p_wilt_sl(j,i)) / (field_cap_sl(j,i) - p_wilt_sl(j,i))
            ELSE
              drain_sl(j) = 0._wp
            END IF

            ! Fast drainage of mobile water

            ! Note: here we use wsat_sl (from porosity, reduced by amount of soil ice) and no longer
            !    field_cap_sl as in the above drainage calculations.

            IF (wtr_soil_sl(j,i) > ws_crit(j)) THEN  ! enough liquid water for fast drainage
              drain_crit(j) = (zdrmax-zdrmin) * ((wtr_soil_sl(j,i)-ws_crit(j)) / (wsat_sl(j,i) - ws_crit(j)))**zdrexp
              drain_sl(j) = drain_sl(j) + MIN(drain_crit(j), (wtr_soil_sl(j,i)-ws_crit(j))/dtime)
            END IF
            ! drainage is limited by the available water of the layer
            drain_sl(j) = MIN(drain_sl(j), (wtr_soil_sl(j,i)-p_wilt_sl(j,i))/dtime)
          END IF
          percolation(j,i) = drain_sl(j) * dtime  ! [m/s] -> [m]

        ELSE IF (dsoil(j,i) > 0._wp) THEN  ! there are soil layers below

          ! Lateral Drainage (i.e. drainage from upper soil layers)
          !      The idea is that drainage occurs through cracks from an upper soil layer directly
          !      down to the bedrock layer.
          ! Todo: This concept needs to be discussed. Should lateral drainage occur at all?

          ! Drainage is calculated with percolation removed from wtr_soil_sl
          !   de Vrese --- wtr_soil_sl corrected for incomming percolation

          ! Remove percolating water from the corresponding soil layer and add percolation from the above layer
          ! Note: In contrast to the below 'real' update of soil water w.r.t. percolation,
          !       ws_min_drain and field_cap_sl are not accounted for, here. Besides, water from
          !       the above layer is added immediately.
          IF (.NOT. soil_above(j)) THEN       ! surface layer
            ws_tmp(j) = wtr_soil_sl(j,i) - percolation(j,i)
          ELSE
            ws_tmp(j) = wtr_soil_sl(j,i) - percolation(j,i) + percolation(j,i-1)
          END IF

          IF (ws_tmp(j) > p_wilt_sl(j,i)) THEN          ! drainage occurs only with mobile water
            ! Slow drainage
            IF (field_cap_sl(j,i) > p_wilt_sl(j,i)) THEN
              drain_sl(j) = zdrmin * (ws_tmp(j) - p_wilt_sl(j,i)) / (field_cap_sl(j,i)-p_wilt_sl(j,i))
            ELSE
              drain_sl(j) = 0._wp
            END IF
            ! Fast drainage
            IF (ws_tmp(j) > ws_crit(j)) THEN
              drain_crit(j) = (zdrmax-zdrmin) * ((ws_tmp(j) - ws_crit(j)) / (wsat_sl(j,i) - ws_crit(j)))**zdrexp
              drain_sl(j) = drain_sl(j) + MIN(drain_crit(j), (ws_tmp(j) - ws_crit(j))/dtime)
            END IF
            ! drainage is limited by the soil moisture within the layer
            drain_sl(j) = MIN(drain_sl(j) * dtime, (ws_tmp(j) - p_wilt_sl(j,i)))
          END IF
        END IF

        IF (.NOT. ltpe_closed) THEN
          IF (dsoil(j,i) > 0._wp .AND. soil_below(j) .AND. ws_tmp(j) > p_wilt_sl(j,i)) THEN
            wtr_soil_sl(j,i) = wtr_soil_sl(j,i) - drain_sl(j) ! remove drainage of the layer from soil moisture
            drainage(j) = drainage(j) + drain_sl(j)           ! and add it to the total drainage
          END IF
        END IF                                       ! no drainage for TPE closed case

      END DO  ! soil level

      !----------------------------------------------------------------------------------------------
      !  water diffusion
      ! -----------------

      ! numerical recipes/Richtmyer & -Morton diffusion

      ! calculating the diffusivity at mid depth of layer i
      !$claw acc loop seq
      DO i = 1, nsoil

        ws_rel(j,i) = 0._wp
        IF (wtr_soil_sl(j,i) > 0._wp) THEN
          ws_rel(j,i) =  wtr_soil_sl(j,i) / dsoil(j,i)  &
            &           / (vol_porosity_sl(j,i) - (ice_soil_sl(j,i) / dsoil(j,i)))
        END IF

        diffus_sl(j,i) = 0._wp
        IF (ws_rel(j,i) > zfcmin) THEN
          diffus_sl(j,i) =  bclapp_sl(j,i) * hyd_cond_sat_sl(j,i) * matrix_pot_sl(j,i)  &
            &               / (wtr_soil_sl(j,i) / dsoil(j,i)) *  ws_rel(j,i) ** (bclapp_sl(j,i)+3._wp)

        END IF
      END DO

      ! calculating the diffusivity at the bottom of layer i
      !$claw acc loop seq
      DO i = 1, nsoil-1
        IF (diffus_sl(j,i+1) == 0._wp .OR. diffus_sl(j,i) == 0._wp) THEN
          diffusivity(j,i) = 0._wp
        ELSE
          diffusivity(j,i) = diffus_sl(j,i+1) * diffus_sl(j,i) &
            &                 / ((diffus_sl(j,i+1)*dsoil(j,i) + diffus_sl(j,i)*dsoil(j,i+1))/(dsoil(j,i+1)+dsoil(j,i)))
        END IF
      END DO

      ! calculation of diffusion coefficients

      ! deepest layer
      zda(j,nsoil) = 0._wp
      IF (dsoil(j,nsoil) > 0._wp) THEN
        zdb(j,nsoil) = diffusivity(j,nsoil-1) * dtime / dsoil(j,nsoil) / (dsoil(j,nsoil)+dsoil(j,nsoil-1)) * 2._wp
      ELSE
        zdb(j,nsoil) = 0._wp
      END IF

      ! middle layers
      !$claw acc loop seq
      DO i = nsoil-1, 2, -1
        IF (dsoil(j,i) > 0._wp) THEN
          IF (dsoil(j,i+1) > 0._wp) THEN
            zda(j,i) = diffusivity(j,i) * dtime / dsoil(j,i) / (dsoil(j,i)+dsoil(j,i+1)) * 2._wp
          ELSE
            zda(j,i) = 0._wp
          END IF
          zdb(j,i) = diffusivity(j,i-1) * dtime / dsoil(j,i) / (dsoil(j,i)+dsoil(j,i-1)) * 2._wp
        ELSE
          zda(j,i) = 0._wp
          zdb(j,i) = 0._wp
        END IF
      END DO

      ! uppermost layer
      zdb(j, 1) = 0._wp
      IF (dsoil(j,1) > 0) THEN
        IF (dsoil(j,2) > 0) THEN
          zda(j,1) = diffusivity(j,1) * dtime / dsoil(j,1) / (dsoil(j,1)+dsoil(j,2)) * 2._wp
        ELSE
          zda(j,1) = 0._wp
        END IF
      ELSE
        zda(j,1) = 0._wp
      END IF

      ! temporary storage for diagnostic printing
      IF (ilog == 1 .OR. ilog == 3) THEN 
        DO i = 1, nsoil
          wslog(i) = wtr_soil_sl(jllog, i)
        END DO
      END IF


      ! ---------------------
      !  update soil wetness  - with respect to diffusion
      ! ---------------------
      ! routine TRIDIAG from Numerical Recipes, p. 43
      !   -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, wtr_soil_sl(T)=RI, wtr_soil_sl(T+1)=UI

      IF (dsoil(j,1) > 0._wp) THEN
        zdc(j,1) = zda(j,1) + zdb(j,1) + 1._wp
        wtr_soil_sl(j,1) = wtr_soil_sl(j,1) / zdc(j,1)
      END IF

      ! decomposition and forward substitution
      !$claw acc loop seq
      DO i = 2, nsoil
        IF (dsoil(j,i) > 0._wp) THEN
          ztri(j,i) = -zda(j,i-1) / zdc(j,i-1)
          zdc(j,i)  = zda(j,i) + zdb(j,i) + 1._wp + zdb(j,i)*ztri(j,i)
          wtr_soil_sl(j,i) = (wtr_soil_sl(j,i) / dsoil(j,i) + zdb(j,i) * wtr_soil_sl(j,i-1)/dsoil(j,i-1)) &
           &                 / zdc(j,i) * dsoil(j,i)
        END IF
      END DO

      ! backsubstitution
      !$claw acc loop seq
      DO i = nsoil-1, 1, -1
        IF (dsoil(j,i+1) > 0._wp) THEN
          wtr_soil_sl(j,i) = wtr_soil_sl(j,i) - ztri(j,i+1) * wtr_soil_sl(j,i+1) * (dsoil(j,i)/dsoil(j,i+1))
        END IF
      END DO

      ! temporary storage for diagnostic printing
      IF (ilog == 1 .OR. ilog == 3) THEN
        !!$claw acc loop seq ! This raises an error by the compiler
        DO i = 1, nsoil
          zdiflog(i) = wslog(i) - wtr_soil_sl(jllog,i)
        END DO
      END IF


      ! treatment of layer overflow due to diffusion
      IF (.NOT. ltpe_closed) THEN
        ! limit soil moisture to field capacity and add overflow to drainage
        !$claw acc loop seq
        DO i=1,nsoil
          IF (dsoil(j,i) > 0._wp .AND. wtr_soil_sl(j,i) > field_cap_sl(j,i)) THEN
            drainage(j) = drainage(j) + wtr_soil_sl(j,i)-field_cap_sl(j,i)
            wtr_soil_sl(j,i) = field_cap_sl(j,i)
          END IF
        END DO
      ELSE
        ! Modification for TPE (closed)
        !$claw acc loop seq
        DO i=1,nsoil
          IF (dsoil(j,i) > 0._wp) THEN
            IF (wtr_soil_sl(j,i) > field_cap_sl(j,i)) THEN
              ! redirect overflow from drainage to overflow pool
              tpe_overflow(j) = tpe_overflow(j) + wtr_soil_sl(j,i) - field_cap_sl(j,i)
              wtr_soil_sl(j,i) = field_cap_sl(j,i)
            ELSE ! wtr_soil_sl <= field_cap_sl
              ! refill soil moisture from overflow pool
              remoist(j) = MIN(tpe_overflow(j), field_cap_sl(j,i) - wtr_soil_sl(j,i))
              wtr_soil_sl(j,i) = wtr_soil_sl(j,i) + remoist(j)
              tpe_overflow(j) = tpe_overflow(j) - remoist(j)
            END IF
          END IF
        END DO
      END IF

      ! Modification for TPE (open)
      IF (ltpe_open) THEN
        !$claw acc loop seq
        DO i = 1, nsoil
          !IF (i == nsoil - 2 .OR. i == nsoil - 1 .OR. i == nsoil) THEN   ! special adaptation of Shirisha
          IF (i == nsoil-1 .OR. i == nsoil) THEN
            IF (dsoil(j,i) > 0._wp) THEN
             ! keep the two lowest soil layers always very wet
             wtr_soil_sl(j,i) = MAX(wtr_soil_sl(j,i),field_cap_sl(j,i) * 0.9_wp)
            END IF
          END IF
        END DO
      END IF

#ifndef _OPENACC
      !----------------------------------------------------------------------------------------------
      ! debugging output for grid cell jllog
      IF (ilog == 2 .OR. ilog == 3) THEN
        write(*,'(A10,1X, 11(E16.9E2, 1X) )')                                                            &
          & "SOILHYD1: ", hyd_cond_sat_sl(jllog,1:nsoil)*1000._wp, wtr_soil_sl(jllog,1:nsoil)*1000._wp,  &
          & field_cap_sl(jllog,1:nsoil)*1000._wp
      END IF
     !----------------------------------------------------------------------------------------------
#endif

      ! ---------------------
      !  update soil wetness  - with respect to percolation
      ! ---------------------
      ! todo: to be discussed: does the below method make sense? Wouldn't it be more natural (and
      !    clearer) to give up the bucket-like concept of saturated soil water piling up from below?

      ! 1. Take away drained soil water from the upper soil layer

      ! percolation per time step must be smaller than the usable soil layer water content
      !$claw acc loop seq
      DO i = 1, nsoil
        ws_min_drain(j) = 0._wp
        IF (dsoil(j,i) > 0._wp .AND. percolation(j,i) > 0._wp) THEN
          ! Note: this calculation of ws_min_drain is different from the calcuation in
          ! calc_soil_hydrology (with respect to the root whole zone), as field_cap_sl (used here)
          ! is reduced by the amount of ice, whereas ws_max (used above) does not account for ice.
          ws_min_drain(j) = zwdmin * field_cap_sl(j,i)              ! minimum soil moisture for drainage
          IF (wtr_soil_sl(j,i)-percolation(j,i) < ws_min_drain(j)) THEN   ! more percolation than water available
            IF (wtr_soil_sl(j,i) > ws_min_drain(j)) THEN                  !   there is water for perculation
              percolation(j,i) = wtr_soil_sl(j,i) - ws_min_drain(j)       !   all available water perculates
              wtr_soil_sl(j,i) = ws_min_drain(j)                          !   soil moisture is limited to minimum
            ELSE
              percolation(j,i) = 0._wp
            END IF
          ELSE                                                ! there is enough soil water
            wtr_soil_sl(j,i) = wtr_soil_sl(j,i) - percolation(j,i)
          END IF
        END IF
      END DO

      ! 2. Add drained soil water to the soil layer below

      ! Percolation is added to downward layer or, if lowest layer or bedrock is reached, to
      ! drainage (not TPE closed) resp. overflow pool (TPE closed)
      !$claw acc loop seq
      DO i = 1, nsoil-1
        IF (dsoil(j,i+1) > 0._wp) THEN
          wtr_soil_sl(j,i+1) = wtr_soil_sl(j,i+1) + percolation(j,i)
        END IF
        IF (.NOT. ltpe_closed) THEN
          IF (dsoil(j,i+1) <= 0._wp) THEN      ! bedrock below
            drainage(j) = drainage(j) + percolation(j,i)
            percolation(j,i) = 0._wp
          END IF
        ELSE  ! TPE closed case
          IF (dsoil(j,i+1) <= 0._wp) THEN      ! bedrock below
            tpe_overflow(j) = tpe_overflow(j) + percolation(j,i)
            percolation(j,i) = 0._wp
          END IF
        END IF
      END DO

      ! 3. Correct percolation if maximum moisture is exceeded, and correct soil moisture of the
      !    above layer accordingly

      ! Percolation overflow is not allowed from layers below the surface layer.
      ! Instead it is assumed that water piles upwards from saturated layers.
      !$claw acc loop seq
      DO i = nsoil, 2, -1
        IF (dsoil(j,i) > 0._wp .AND. wtr_soil_sl(j,i) > field_cap_sl(j,i)) THEN
          ! soil moisture above field capacity:
          !   reduce percolation from the above layer by the amount that does not fit into the below layer
          percolation(j,i-1) = MAX(percolation(j,i-1) - wtr_soil_sl(j,i) + field_cap_sl(j,i), 0._wp)
          wtr_soil_sl(j,i-1) = wtr_soil_sl(j,i-1) + wtr_soil_sl(j,i) - field_cap_sl(j,i)
          wtr_soil_sl(j,i) = field_cap_sl(j,i)
        END IF
      END DO

      IF (.NOT. ltpe_closed) THEN
        drainage(j) = drainage(j) + percolation(j,nsoil)
      ELSE  ! TPE closed case
        ! redirect percolation from drainage to overflow pool
        tpe_overflow(j) = tpe_overflow(j) +  percolation(j,nsoil)
      END IF

      ! 4. Add overflow of the surface layer to drainage

      ! Overflow from the surface layer enters drainage, resp. the overflow pool for TPE closed.
      ! This should not happen, this is just for security.
      IF (.NOT. ltpe_closed) THEN
        IF (wtr_soil_sl(j,1) > field_cap_sl(j,1)) THEN
          drainage(j) = drainage(j) + wtr_soil_sl(j,1) - field_cap_sl(j,1)
          wtr_soil_sl(j,1) = field_cap_sl(j,1)
        END IF
      ELSE      ! TPE closed case
        IF (wtr_soil_sl(j,1) > field_cap_sl(j,1)) THEN
          ! redirect overflow to overflow pool
          tpe_overflow(j) = tpe_overflow(j) + wtr_soil_sl(j,1) - field_cap_sl(j,1)
          wtr_soil_sl(j,1) = field_cap_sl(j,1)
        END IF
      END IF
    END DO
    !$claw acc end parallel
    !$claw acc end data

#ifndef _OPENACC
      !----------------------------------------------------------------------------------------------
      ! debugging output for grid cell jllog
      IF (ilog == 1 .OR. ilog == 3) THEN
        DO i=1, nsoil
          zdiflog(nsoil+i) = percolation(jllog,i)
        END DO
        write(*,'(A10,1X,11(E15.6E3,1X))') "SOILFLUX",zdiflog(1:10)*1000._wp, drainage(jllog)*1000._wp
        write(*,'(A10,1X,10(E16.9E2,1X))') "SOILHYD2: ", wtr_soil_sl(jllog,1:nsoil)*1000._wp, dsoil(jllog,1:nsoil)*1000._wp
      END IF
#endif

  END SUBROUTINE soilhyd

  SUBROUTINE digest_evapotrans( &
    & nc,                       &
    & nsoil,                    &
    & droot,                    &
    & dsoil,                    &
    & field_cap_sl,             &
    & transpiration_in_m,       &
    & evapo_soil_in_m,          &
    & wtr_soil_sl, ice_soil_sl  &
    & )

    !----------------------------------------------------------------------------------------------
    !  The routine is based on routine soilchange by Stefan Hagemann. It replaces the part for
    !  evapotranspiration changes (isch=1)
    !
    !  Update soil moisture due to evaporation and transpiration fluxes. The fluxes had been
    !  calculated with the old bucket scheme, and are not necessarily consistent to soil moisture
    !  of the multi-layer scheme.
    !  A potential water deficit is distributed between the soil layers.
    !----------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nc                      ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    REAL(wp), INTENT(in)     :: droot(:,:)              ! root depth within layer [m]
    REAL(wp), INTENT(in)     :: dsoil(:,:)              ! soil depth until bedrock within layer [m]
    REAL(wp), INTENT(in)     :: field_cap_sl(:,:)       ! field capacity of the layer (reduced by ice)
    REAL(wp), INTENT(in)     :: transpiration_in_m(:)   ! transpiration [m]
    REAL(wp), INTENT(in)     :: evapo_soil_in_m(:)      ! bare soil evaporation [m]

    REAL(wp), INTENT(inout)  :: wtr_soil_sl(:,:)        ! water content of the soil layer [m]
    REAL(wp), INTENT(inout)  :: ice_soil_sl(:,:)        ! ice content of the soil layer [m]

    ! local variables

    INTEGER  :: jk, i
    REAL(wp) :: deficit_evapo(nc)  ! water deficit due to evaporation
    REAL(wp) :: deficit_trans(nc)  ! water deficit due to tranpiration
    REAL(wp) :: deficit_sl(nc)     ! water deficit within a soil layer
    REAL(wp) :: remaining(nc)      ! water remaining in the soil layer
    REAL(wp) :: rootfract(nc)      ! fraction of the soil layer within the root zone
    REAL(wp) :: fixed(nc)          ! amount of water below the wilting point, not available for plant
    REAL(wp) :: trans_layer(nc)    ! transpiration from the layer
    REAL(wp) :: ws(nc)             ! water and ice within all soil layers
    REAL(wp) :: root_depth(nc)     ! Total depth of rootzone
    REAL(wp) :: deficit_evapotrans(nc) ! water deficit due to evapotranspiration

    INTEGER :: max2d(2)

    ! parameters
    ! @todo - take from mo_hydro_constants.f90
    REAL(wp), PARAMETER :: zwilt  = 0.35_wp

    !----------------------------------------------------------------------------------------------
    !  update soil moisture due to evaporation and transpiration fluxes
    !----------------------------------------------------------------------------------------------

    !$claw acc data present(droot, dsoil, field_cap_sl, transpiration_in_m, evapo_soil_in_m, &
    !$claw acc wtr_soil_sl, ice_soil_sl) pcreate(deficit_evapo, deficit_trans, deficit_sl, &
    !$claw acc remaining, rootfract, fixed, trans_layer, ws, root_depth, deficit_evapotrans)

    !  bare soil evaporation
    ! -----------------------
    !$claw acc parallel
    !$claw acc loop gang vector
    DO i = 1, nc
      deficit_evapo(i) = 0._wp   ! water missing for evaporation due to limited soil water
      IF (evapo_soil_in_m(i) < 0._wp) THEN          ! pos. evaporation flux (not dew)
        wtr_soil_sl(i,1) = wtr_soil_sl(i,1) + evapo_soil_in_m(i) ! take water from uppermost soil layer
        IF (wtr_soil_sl(i,1) < 0._wp) THEN                    ! if there is not enough soil water
          deficit_evapo(i) = wtr_soil_sl(i,1)                 !    this is the amount of water missing
          wtr_soil_sl(i,1) = 0._wp                            !    soil moisture cannot be negative
        END IF
      END IF

      !  transpiration
      ! ---------------

      root_depth(i) = SUM(droot(i,:))

      ! sealed grid cells (rooting_depth zero) and positive transpiration
      !    => no water reachable for transpiration, add flux to deficit term
      ! @todo: what does rooting depth=zero for bare soil tile?
      IF(root_depth(i) <= 0._wp .AND. transpiration_in_m(i) < 0._wp) THEN
        deficit_trans(i) = transpiration_in_m(i)
      ELSE
        deficit_trans(i) = 0._wp
      END IF

      !$claw acc loop seq
      DO jk =1, nsoil

        deficit_sl(i)   = 0._wp
        rootfract(i)    = 0._wp
        fixed(i)        = 0._wp
        trans_layer(i)  = 0._wp
        remaining(i)    = 0._wp

        ! grid cell not sealed and positive transpiration (not dew)
        !    => take water equally from all over the root zone
        IF (root_depth(i) > 0._wp .AND. transpiration_in_m(i) < 0._wp) THEN

          ! water needed for transpiration from the current layer
          trans_layer(i) = transpiration_in_m(i) * droot(i,jk)/root_depth(i)

          IF (droot(i,jk) > 0._wp) THEN

            rootfract(i) = droot(i,jk) / dsoil(i,jk)             ! fraction of the soil layer that is within the root zone
            ! todo: should we use here p_wilt_sl (instead of filed_cap_sl*zwilt)?
            fixed(i) = field_cap_sl(i,jk) * zwilt * rootfract(i) ! amount of water in the root zone, not available for plants

            IF (wtr_soil_sl(i,jk) * rootfract(i) >= fixed(i)) THEN      ! there is water in the root zone available for plants
              remaining(i) = wtr_soil_sl(i,jk) * rootfract(i) + trans_layer(i)  ! water that would remain after transpiration
              IF (remaining(i) < fixed(i)) THEN                           ! transpiration exceeds available water
                deficit_sl(i) = remaining(i) - fixed(i)
                wtr_soil_sl(i,jk) = fixed(i)     + wtr_soil_sl(i,jk) * (1._wp-rootfract(i))
              ELSE                                                        ! enough soil water
                deficit_sl(i) = 0._wp
                wtr_soil_sl(i,jk) = remaining(i) + wtr_soil_sl(i,jk) * (1._wp-rootfract(i))
              END IF
            ELSE                                            ! too dry: no water available for transpiration
              deficit_sl(i) = trans_layer(i)
            END IF

            ! sum up deficits from the different layers
            deficit_trans(i) = deficit_trans(i) + deficit_sl(i)

          END IF
        END IF
      END DO

      ! -----------------------------------------
      !  Total deficit from evapotranspiration
      ! -----------------------------------------
      deficit_evapotrans(i) = deficit_evapo(i) + deficit_trans(i)

      ! --------------------------
      !  soil moisture correction
      ! --------------------------
      ! In case of a water deficit due to evapotranspiration, soil moisture must be
      ! changed to close the land surface water balance.

      ws(i) = SUM(wtr_soil_sl(i,:)) + SUM(ice_soil_sl(i,:))  ! water and ice within the column
      DO jk=1, nsoil
        IF (deficit_evapotrans(i) < 0._wp .AND. ws(i) > 0._wp) THEN
          ! reduce soil moisture and ice relative to the water content of each layer
          wtr_soil_sl(i,jk) = MAX(0._wp, wtr_soil_sl(i,jk) + deficit_evapotrans(i) * wtr_soil_sl(i,jk)/ws(i))
          ice_soil_sl(i,jk) = MAX(0._wp, ice_soil_sl(i,jk) + deficit_evapotrans(i) * ice_soil_sl(i,jk)/ws(i))
        END IF
      END DO
    END DO ! nc loop

    !$claw acc end parallel
    !$claw acc end data

#ifndef _OPENACC
    ! make sure there is no negative soil ice
    IF (ANY(ice_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil ice: ', MINVAL(ice_soil_sl), ' at ', MINLOC(ice_soil_sl)
      CALL finish ('digest_evapotrans', message_text)
      ! CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      ! ice_soil_sl(:,:) = MAX(ice_soil_sl(:,:), 0._wp)
    END IF

    ! make sure there is no negative soil moisture and soil moisture does not exceed field capacity
    IF (ANY(wtr_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture: ', MINVAL(wtr_soil_sl), ' at ', MINLOC(wtr_soil_sl)
      CALL finish ('digest_evapotrans', message_text)
      ! CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      ! wtr_soil_sl(:,:) = MAX(wtr_soil_sl(:,:), 0._wp)
    END IF
    IF (ANY(wtr_soil_sl(:,:) > field_cap_sl(:,:) + zfcmin)) THEN
      max2d(:) = MAXLOC(wtr_soil_sl(:,:) - field_cap_sl(:,:))
      WRITE (message_text,*) 'liquid water content above field capacity at ', max2d(:), &
        &   ': wtr_soil_sl: ', wtr_soil_sl(max2d(1),max2d(2)),                          &
        &   ' field_cap_sl: ', field_cap_sl(max2d(1),max2d(2))
      CALL finish ('digest_evapotrans', message_text)
      ! CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      ! wtr_soil_sl(:,:) = MIN(wtr_soil_sl(:,:), field_cap_sl(:,:))
    END IF
#endif

  END SUBROUTINE digest_evapotrans
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE digest_infilt ( &
    & nc,                    &
    & nsoil,                 &
    & ltpe_closed,           &
    & field_cap_sl,          &
    & infilt,                &
    & wtr_soil_sl,           &
    & drainage,              &
    & tpe_overflow           &
    & )

  !------------------------------------------------------------------------------------------------
  !  The routine is based on routine soilchange by Stefan Hagemann. It replaces the part for
  !  infiltration changes (isch=2)
  !------------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nc                      ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    LOGICAL,  INTENT(in)     :: ltpe_closed             ! TPE closed case
    REAL(wp), INTENT(in)     :: field_cap_sl(:,:)       ! maximum amount of liquid water in layer
    REAL(wp), INTENT(in)     :: infilt(:)               ! infiltration in [m/s]

    REAL(wp), INTENT(inout)  :: wtr_soil_sl(:,:)        ! water content of the soil layer [m]
    REAL(wp), INTENT(inout)  :: drainage(:)             ! drainage [m]
    REAL(wp), INTENT(inout)  :: tpe_overflow(:)         ! For TPE

    ! local variables

    INTEGER  :: jk, i
    REAL(wp) :: excess(nc)         ! excess water after infiltration

    !------------------------------------------------------------------------------------------------
    !  update soil moisture according to infiltration
    !------------------------------------------------------------------------------------------------

    !$claw acc data &
    !$claw acc present(field_cap_sl, infilt, wtr_soil_sl, drainage, tpe_overflow) &
    !$claw acc create(excess)
    !$claw acc parallel
    !$claw acc loop gang vector
    DO i = 1, nc
      excess(i) = infilt(i)                   ! amount of water for infiltration at the surface [m]
      !$claw acc loop seq
      DO  jk = 1, nsoil
        wtr_soil_sl(i,jk) = wtr_soil_sl(i,jk) + excess(i)
        IF (.NOT. ltpe_closed) THEN
          IF (wtr_soil_sl(i,jk) > field_cap_sl(i,jk)) THEN       ! soil moisture exceeds field capacity
            excess(i) = wtr_soil_sl(i,jk) - field_cap_sl(i,jk)   !   excess water -> goes into level below
            wtr_soil_sl(i,jk) = field_cap_sl(i,jk)               !   limit soil moisture to field capacity
          ELSE
            excess(i) = 0._wp
          END IF
        ELSE  ! TPE closed case
          IF (wtr_soil_sl(i,jk) > field_cap_sl(i,jk)) THEN
            tpe_overflow(i) = tpe_overflow(i) + wtr_soil_sl(i,jk) - field_cap_sl(i,jk)
            wtr_soil_sl(i,jk) = field_cap_sl(i,jk)
          END IF
          excess(i) = 0._wp
        END IF
      END DO

      ! if all soil layers are at field capacity and there is still excess water, it is added to drainage
      IF (excess(i) > 0._wp) THEN
        drainage(i) = drainage(i) + excess(i)
      END IF
    END DO
    !$claw acc end parallel
    !$claw acc end data

    ! make sure there is no negative soil moisture
!TODO deal with that
#ifndef _OPENACC
    IF (ANY(wtr_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture: ', MINVAL(wtr_soil_sl), ' at ', MINLOC(wtr_soil_sl)
      IF (ANY(wtr_soil_sl(:,:) < -EPSILON(1._wp))) THEN
        CALL finish ('digest_infilt', message_text)
      ELSE
        CALL message ('digest_infilt', message_text, all_print=.TRUE., level=em_warn)
        wtr_soil_sl(:,:) = MAX(wtr_soil_sl(:,:), 0._wp)
      END IF
    END IF
#endif

  END SUBROUTINE digest_infilt

  !--------------  !!
  ! @brief      Calculates the snow fraction on lake ice.
  !
  ! @param [in]     weq_snow_lice    Snow depth on lake ice [m water equivalent]
  ! @param [out]    fract_snow_lice  Fraction of snow on lake ice
  !!
  ELEMENTAL SUBROUTINE calc_wskin_fractions_lice( &
    & weq_snow_lice,                              & ! in
    & fract_snow_lice                             & ! out
    & )

    !$claw model-data
    REAL(wp),  INTENT(in)    :: weq_snow_lice
    REAL(wp),  INTENT(out)   :: fract_snow_lice
    !$claw end model-data

    !$claw sca

    fract_snow_lice = TANH(weq_snow_lice * 100._wp)

  END SUBROUTINE calc_wskin_fractions_lice

  ELEMENTAL SUBROUTINE calc_wskin_fractions_veg( &
    & dtime,                                     & ! in
    & steplen,                                   & ! in
    & skinres_max,                               & ! in
    & oro_stddev,                                & ! in
    & t_srf_old,                                 & ! in
    & press_srf,                                 & ! in
    & drag_srf,                                  & ! in
    & q_air,                                     & ! in
    & wtr_skin,                                  & ! in
    & weq_snow_soil,                             & ! in
    & weq_snow_can,                              & ! in
    & fract_snow_can,                            & ! inout
    & fract_wet    ,                             & ! out
    & fract_snow_soil                            & ! out
    & )

    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in) :: &
      & dtime,              &
      & steplen

    !$claw model-data
    REAL(wp), INTENT(in) :: &
      & oro_stddev,         &
      & t_srf_old,          &
      & press_srf,          &
      & drag_srf,           &
      & q_air,              &
      & skinres_max,        &
      & wtr_skin,           &
      & weq_snow_soil,      &
      & weq_snow_can

    REAL(wp), INTENT(out) :: &
      & fract_snow_can,      &
      & fract_wet    ,       &
      & fract_snow_soil
    !$claw end model-data

    REAL(wp) ::       &
      & qsat_srf_old, &
      & evapo_pot

    !$claw sca

    ! Snow fraction on canopy
    ! Wet skin reservoir fraction on soil and canopy (between 0 and 1)
    IF (skinres_max > EPSILON(1._wp)) THEN
      fract_snow_can = MIN(1._wp, weq_snow_can / skinres_max)
      fract_wet      = MIN(1._wp, wtr_skin     / skinres_max)
    ELSE
      fract_snow_can = 0._wp
      fract_wet      = 0._wp
    END IF

    ! Snow fraction on soil
    !$claw sca forward
    fract_snow_soil = Get_snow_fract_noforest(weq_snow_soil, oro_stddev) ! snow on soil below forest
    fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, &
                               fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))


    ! Potential evaporation using old values of air and surface humidity
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evapo_pot =  -1._wp * heat_transfer_coef(drag_srf, steplen) * (q_air - qsat_srf_old)  ! Positive upwards

    ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy; same for skin reservoir
    IF (fract_snow_soil > 0._wp) THEN
      ! @todo Shouldn't one take rhoice here insteady of rhoh2o?
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evapo_pot * dtime         &
        &                                            / (rhoh2o * (weq_snow_soil + weq_snow_can)) &
        &                                    )
    END IF
    IF (fract_wet > 0._wp ) THEN
      fract_wet = fract_wet / MAX(1._wp, (1._wp - fract_snow_soil) * evapo_pot * dtime   &
        &                                    / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin))  &
        &                                )
    END IF

  END SUBROUTINE calc_wskin_fractions_veg

  ELEMENTAL SUBROUTINE calc_wskin_fractions_bare( &
    & dtime,                                      & ! in
    & steplen,                                    & ! in
    & skinres_max,                                & ! in
    & oro_stddev,                                 & ! in
    & t_srf,                                      & ! in
    & press_srf,                                  & ! in
    & drag_srf,                                   & ! in
    & q_air,                                      & ! in
    & wtr_skin,                                   & ! in
    & weq_snow_soil,                              & ! in
    & fract_wet,                                  & ! out
    & fract_snow_soil                             & ! out
    & )

    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in) :: &
      & dtime,              &
      & steplen

    !$claw model-data
    REAL(wp), INTENT(in) :: &
      & oro_stddev,         &
      & t_srf,              &
      & press_srf,          &
      & drag_srf,           &
      & q_air,              &
      & skinres_max,        &
      & wtr_skin,           &
      & weq_snow_soil

    REAL(wp), INTENT(out) :: &
      & fract_wet,           &
      & fract_snow_soil
    !$claw end model-data

    REAL(wp) ::       &
      & t_srf_old,    &
      & qsat_srf_old, &
      & evapo_pot

    !$claw sca

    IF (skinres_max > EPSILON(1._wp)) THEN
      fract_wet = MIN(1._wp, wtr_skin / skinres_max)
    ELSE
      fract_wet = 0._wp
    END IF

    ! Snow fraction on soil
    !$claw sca forward
    fract_snow_soil = Get_snow_fract_noforest(weq_snow_soil, oro_stddev) ! snow on soil below forest

    ! Potential evaporation using old values of air and surface humidity
    t_srf_old = t_srf
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evapo_pot =  heat_transfer_coef(drag_srf, steplen) * (qsat_srf_old - q_air) ! Positive upwards

    ! Reduce snow cover on soil if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy together
    IF (fract_snow_soil > 0._wp .AND. weq_snow_soil > EPSILON(1._wp)) THEN
      fract_snow_soil = fract_snow_soil &
        & / MAX(1._wp, fract_snow_soil * evapo_pot * dtime / (rhoh2o * MAX(EPSILON(1._wp), weq_snow_soil)))
    END IF
    ! Same for water cover on soil
    IF (fract_wet > 0._wp .AND. wtr_skin > EPSILON(1._wp)) THEN
      fract_wet = fract_wet &
        & / MAX(1._wp, (1._wp - fract_snow_soil) * evapo_pot * dtime / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin)))
    END IF

  END SUBROUTINE calc_wskin_fractions_bare

  !< Roesch et al. 2002, Climate Dynamics
  !
  REAL(wp) ELEMENTAL FUNCTION Get_snow_fract_noforest(snow, orodev)

    USE mo_hydro_constants, ONLY: wsn2fract_const, wsn2fract_eps, wsn2fract_sigfac

    REAL(wp), INTENT(in) :: &
      & snow,   &
      & orodev

    !$claw sca routine

    Get_snow_fract_noforest = wsn2fract_const * TANH(snow * 100._wp) &
                               & * SQRT(snow * 1000._wp / (snow * 1000._wp + wsn2fract_eps + wsn2fract_sigfac * orodev))

  END FUNCTION Get_snow_fract_noforest

  REAL(wp) ELEMENTAL FUNCTION get_canopy_cond_unstressed_simple(lai, par) RESULT(conductance)

    USE mo_hydro_constants, ONLY: k => conductance_k, a => conductance_a, b => conductance_b, c => conductance_c

    REAL(wp), INTENT(in) ::     &
                           lai, &
                           par

    ! Local variables
    !
!!$    REAL(wp) :: d(SIZE(par))
!!$    REAL(wp) :: zpar(SIZE(par))
    REAL(wp) :: d, zpar

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_canopy_cond_unstressed_simple'

    !$claw sca routine

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    zpar = MAX(1.E-10_wp, par)

!!$    WHERE (lai > EPSILON(1._wp))
    IF (lai > EPSILON(1._wp)) THEN
      d = (a + b*c) / (c * zpar)
      conductance = ( LOG((d * EXP(k*lai) + 1._wp) / (d + 1._wp)) * b / (d * zpar) - &
                    & LOG((d + EXP(-k*lai)) / (d + 1._wp))                           &
                    ) / (k * c)
    ELSE
      conductance = EPSILON(1._wp)
    END IF

  END FUNCTION get_canopy_cond_unstressed_simple

  REAL(wp) ELEMENTAL FUNCTION get_canopy_cond_stressed_simple(cond_unstressed, water_stress, air_is_saturated) RESULT(conductance)

    REAL(wp), INTENT(in) :: &
      & cond_unstressed, &
      & water_stress
    LOGICAL,  INTENT(in) :: air_is_saturated

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_canopy_cond_stressed_simple'

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    !$claw sca routine

    IF (air_is_saturated) THEN
      conductance = EPSILON(1._wp)
    ELSE
      conductance = cond_unstressed * water_stress
    END IF

  END FUNCTION get_canopy_cond_stressed_simple

  !
  ! @brief      Gets the water stress factor.
  !
  ! @param      w_soil               The w soil sl
  ! @param      w_soil_max           The w soil fc sl
  ! @param      w_soil_crit_fract    The w soil crit fract
  ! @param      w_soil_wilt_fract    The w soil pwp sl
  ! @param      water_stress_factor  The water stress factor
  !
  ! @return     The water stress factor. !
  !
  ELEMENTAL FUNCTION get_water_stress_factor ( &
    & w_soil, w_soil_max, w_soil_crit_fract, w_soil_wilt_fract) RESULT(water_stress_factor)

    ! TODO
    ! TBD: Maybe change later so that it uses geographically dependend critical values (like wilting point)

    REAL(wp), INTENT(in) :: w_soil            !< Soil water content
    REAL(wp), INTENT(in) :: w_soil_max        !< Soil water content at field capacity
    REAL(wp), INTENT(in) :: w_soil_crit_fract !< Fraction of max. soil water content at critical point
    REAL(wp), INTENT(in) :: w_soil_wilt_fract !< Fraction of max. soil water content at wilting point
    REAL(wp)             :: water_stress_factor

    REAL(wp) :: w_crit, w_wilt                       !< Soil water content at critical/wilting point

    !$claw sca routine

    w_crit = w_soil_max * w_soil_crit_fract
    w_wilt = w_soil_max * w_soil_wilt_fract

    IF (w_crit - w_wilt > 0._wp) THEN
      water_stress_factor = MAX (0._wp, &
        &                        MIN (1._wp, (w_soil - w_wilt) / (w_crit - w_wilt) ))
    ELSE
      water_stress_factor = 0._wp
    END IF

  END FUNCTION get_water_stress_factor

#endif
END MODULE mo_hydro_process
