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
  USE mo_exception, ONLY: message, finish, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_surface_hydrology_land, calc_surface_hydrology_glacier, calc_surface_runoff, &
    & calc_soil_hydrology, calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare, &
    & get_canopy_conductance, get_water_stress_factor, calc_orographic_features

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  REAL(wp), PARAMETER ::  zwdmin = 0.05_wp

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

  ! ===============================================================================================================================
  !>
  !! @brief   Calculates surface runoff
  
  SUBROUTINE calc_surface_runoff(                       &
    ! in
    & nc, ltpe_closed, l_infil_subzero,                 &
    & steepness, t_soil_sl1, water_to_soil,             &
    & soil_depth_sl, root_depth_sl, vol_field_cap_sl,   &
    & wtr_soil_sl, ice_soil_sl,                         &
    & infilt, runoff                                    &
    & )

    USE mo_hydro_util,             ONLY: get_amount_in_rootzone

    INTEGER,  INTENT(in) :: nc
    LOGICAL, INTENT(in)  :: ltpe_closed
    LOGICAL, INTENT(in)  :: l_infil_subzero

    REAL(wp), INTENT(in), DIMENSION(:) ::               &
      & steepness,           & !> Shape parameter describing orography    [/]
      & t_soil_sl1,          & !> Temperature of the uppermost soil layer [K]
      & water_to_soil          !> Amount of water at the soil surface     [m]

    REAL(wp), INTENT(in), DIMENSION(:,:) ::             &
      & soil_depth_sl, root_depth_sl, vol_field_cap_sl, &
      & wtr_soil_sl, ice_soil_sl

    REAL(wp), INTENT(out), DIMENSION(:) ::              &
      & infilt,              & !> Infiltration flux into the soil         [m]
      & runoff                 !> Surface runoff                          [m]

    ! Local variables
    REAL(wp) ::              &
      & weq_rootzone_max(nc),& ! maximum water holding capacity of the rootzone (liquid water + ice) [m]
      & wtr_rootzone(nc),    & ! amount of liquid water in the rootzone [m]
      & ice_rootzone(nc)       ! amount of ice in the rootzone [m]

    !---------------------
    !> Infiltration and surface runoff
    !---------------------
    !> Runoff is either set to zero for Terraplanet setup or is computed via the
    !! ARNO Scheme (default). Note that both, infiltration and runoff, are further
    !! modified in the soil hydrology routine.

    IF (ltpe_closed) THEN
      ! Terra planet setup without runoff
      runoff(:) = 0._wp
      infilt(:) = water_to_soil(:)

    ELSE
      CALL get_amount_in_rootzone(vol_field_cap_sl(:,:) * soil_depth_sl(:,:), &
        &                         soil_depth_sl(:,:), root_depth_sl(:,:), weq_rootzone_max(:))
      CALL get_amount_in_rootzone(ice_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

      CALL arno_scheme(nc, l_infil_subzero, t_soil_sl1,                          &
        &              water_to_soil(:), wtr_rootzone(:) + ice_rootzone(:),      &
        &              steepness(:), weq_rootzone_max(:),                        &
        &              runoff(:), infilt(:))

    END IF

  END SUBROUTINE calc_surface_runoff

  ! ===============================================================================================================================
  !>
  !! @brief   Calculates vertical water movement through the soil column
  !
  ! @param [in]      nc                  Vector lenth
  ! @param [in]      nsoil               Number of soil layers
  ! @param [in]      dtime               Time step [s]
  ! @param [in]      l_infil_subzero     Allow infiltration at soil temperatures below 0degC
  ! @param [in]      soil_depth_sl       Thicknesses of soil layers until bedrocj [m]
  ! @param [in]      root_depth_sl       Thicknesses of soil layers until rooting depth [m]
  ! @param [in]      hyd_cond_sat_sl     Hydraulic conductivity of saturated soil [m/s]
  ! @param [in]      matrix_pot_sl       Soil matrix potential [m]
  ! @param [in]      bclapp_sl           Exponent B in Clapp and Hornberger
  ! @param [in]      pore_size_index_sl  Soil pore size distribution index used in Van Genuchten method
  ! @param [in]      vol_porosity_sl     Soil porosity                      - ice not considered [m/m]
  ! @param [in]      vol_field_cap_sl    Volumetric field capacity          - ice not considered [m/m]
  ! @param [in]      vol_p_wilt_sl       Volumetric permanent wilting point - ice not considered [m/m]
  ! @param [in]      evapotrans_soil     Evapotranspiration from the soil (w/o snow or canopy) [kg m-2 s-1]
  ! @param [in]      transpiration       Evapotranspiration [kg m-2 s-1]
  ! @param [in]      wtr_soil_pot_scool_sl Amount of potentially supercooled water [m]
  ! @param [in,out]  wtr_soil_sl         Content un-frozen water of soil layers in rooting zone  [m]
  ! @param [in,out]  ice_soil_sl         Content of frozen water of soil layers in rooting zone  [m]
  ! @param [out]     wsat_sl             Maximum possible amount of water in soil layer          [m]
  ! @param [out]     field_cap_sl        Water content at field capacity in soil layer           [m]
  ! @param [out]     p_wilt_sl           Water content at permanent wilting point in soil layers [m]
  ! @param [out]     runoff              Runoff
  ! @param [out]     drainage            Drainage
  !
  SUBROUTINE calc_soil_hydrology(                                                        &
    ! in
    & nc, l_fract, lat, lon, nsoil, dtime, ltpe_closed, ltpe_open,                       &
    & soilhydmodel, interpol_mean,                                      &
    & w_soil_wilt_fract, soil_depth_sl, root_depth_sl,                                   &
    & hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, vol_porosity_sl,    &
    & vol_field_cap_sl, vol_p_wilt_sl, wtr_soil_pot_scool_sl,     &
    & evapotrans_soil, transpiration,                                     &
    ! inout
    & infilt, runoff, wtr_soil_sl, ice_soil_sl, tpe_overflow, wtr_soilhyd_res,           &
    ! out
    & wtr_wsat_sl, wtr_field_cap_sl, wtr_p_wilt_sl,                                      &
    & infilt_over, drainage, drain_rock, drain_lat                                       &
    & )

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_util,             ONLY: get_amount_in_rootzone

    INTEGER,  INTENT(in) :: &
      & nc, nsoil
    INTEGER, INTENT(in) :: &
      & soilhydmodel, interpol_mean
    REAL(wp), INTENT(in) :: &
      & dtime, w_soil_wilt_fract
    LOGICAL, INTENT(in) :: ltpe_closed
    LOGICAL, INTENT(in) :: ltpe_open
    LOGICAL, INTENT(in), DIMENSION(:) :: l_fract      ! Tile fraction > 0
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & soil_depth_sl, root_depth_sl, hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, &
      & vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl, wtr_soil_pot_scool_sl
    REAL(wp), INTENT(in), DIMENSION(:) ::               &
      & lat, lon,                              &
      & evapotrans_soil, transpiration
    REAL(wp), INTENT(inout), DIMENSION(:) :: &
      & tpe_overflow, infilt, runoff
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & wtr_soil_sl, ice_soil_sl
    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & infilt_over, drainage,             &
      & drain_rock, drain_lat, wtr_soilhyd_res
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & wtr_wsat_sl, wtr_field_cap_sl, wtr_p_wilt_sl

    ! Local variables
    REAL(wp) ::                       &
      & evapo_soil_in_m(nc),          & !< evaporation from the soil (without transpiration) [m]
      & transpiration_in_m(nc)          !< transpiration within time step [m]

    REAL(wp) ::              &
      & weq_rootzone_max(nc),& ! maximum water holding capacity of the rootzone (liquid water + ice) [m]
      & wtr_rootzone(nc),    & ! amount of liquid water in the rootzone [m]
      & ice_rootzone(nc)       ! amount of ice in the rootzone [m]

    REAL(wp) ::                 &
      & drain_lat_sl(nc,nsoil)

    REAL(wp) ::              &
      & hlp_2d(nc,nsoil)       ! helper array needed with openacc

    INTEGER :: jl, i

    INTEGER, PARAMETER :: ilog = 0 ! Switch for debugging output

    !$claw acc data pcreate(evapo_soil_in_m,transpiration_in_m,infilt,weq_rootzone_max,ws_min_drain,ws_incl_ice,ice_rootzone)
    !$claw acc data present(transpiration_in_m,transpiration,wtr_wsat_sl,wtr_vol_porosity_sl,soil_depth_sl,ice_soil_sl, &
    !$claw acc              wtr_field_cap_sl,vol_field_cap_sl,wtr_p_wilt_sl,vol_p_wilt_sl)
    !$claw acc parallel
    !$claw acc loop gang
    !
    !> Calculation of soil saturation, field capacity and wilting point
    !
    DO i = 1, nc
      !$claw acc loop seq
      DO jl = 1, nsoil
        ! Note: here we are calculating soil saturation (this does not correspond to weq_rootzone_max, which was
        !    calculated from field capacity)
        wtr_wsat_sl     (i,jl) = MAX(vol_porosity_sl (i,jl) * soil_depth_sl(i,jl) - ice_soil_sl(i,jl), 0._wp)  ! -> diffusion
        wtr_field_cap_sl(i,jl) = MAX(vol_field_cap_sl(i,jl) * soil_depth_sl(i,jl) - ice_soil_sl(i,jl), 0._wp)  ! -> percolation
        ! @todo: to be discussed: Does it make sense, to reduce water content at wilting point relative to the ice fraction?
        IF (vol_porosity_sl(i,jl) > 0._wp) THEN
          wtr_p_wilt_sl(i,jl) = wtr_wsat_sl(i,jl) * vol_p_wilt_sl(i,jl) / vol_porosity_sl(i,jl)
        ELSE
          wtr_p_wilt_sl(i,jl) = 0._wp
        END IF
      END DO
    END DO
    !$claw acc end parallel
    !$claw acc end data

    ! Note: weq_rootzone_max calculated here (based on field capacity) is different from wsat calculated
    !       below (based on saturation).
    hlp_2d(:,:) = vol_field_cap_sl(:,:) * soil_depth_sl(:,:)
    CALL get_amount_in_rootzone(hlp_2d(:,:),      soil_depth_sl(:,:), root_depth_sl(:,:), weq_rootzone_max(:))
    CALL get_amount_in_rootzone(ice_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    !$claw acc present(ws_min_drain,weq_rootzone_max,steepness,runoff,infilt)    &
    !$claw acc end parallel

    ! Initialize diagnostic runoff and drainage variables
    infilt_over(:)    = 0._wp
    drain_rock(:)     = 0._wp
    drain_lat_sl(:,:) = 0._wp

    ! Check for unrealistic water to soil fluxes
    DO i = 1, nc
      IF (infilt(i) > 1.0e+3_wp .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'Infiltration flux is far outside of realistic bounds at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Infiltration: ', infilt(i)
        CALL finish ('calc_soil_hydrology', message_text)
      END IF
    END DO

    !
    !$claw acc data present(evapo_soil_in_m,evapotrans_soil,transpiration_in_m,transpiration)
    !$claw acc parallel
    !$claw acc loop gang
    DO i = 1, nc
      transpiration_in_m(i) = transpiration(i)  * dtime / rhoh2o
      ! Note: at this point ignoring the amount of dew that exceeds the skin reservoir
      !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
      evapo_soil_in_m(i) = evapotrans_soil(i) * dtime / rhoh2o - transpiration_in_m(i)
    END DO
    !$claw acc end parallel

    !---------------------
    !> 2. Vertical water transport within the soil column
    !---------------------

    !$claw acc parallel
    !$claw acc loop gang
    DO i = 1, nc
      drainage(i) = 0._wp
    END DO
    !$claw acc end parallel
    !$claw acc end data

    wtr_soilhyd_res(:) = 0._wp

    CALL soilhyd(nc, l_fract, lat, lon, nsoil, dtime, soilhydmodel, interpol_mean,   &
      &     w_soil_wilt_fract, soil_depth_sl, wtr_field_cap_sl, wtr_wsat_sl,         &
      &     hyd_cond_sat_sl, vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl,       &
      &     bclapp_sl, matrix_pot_sl, pore_size_index_sl,                            &
      &     wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl,                         &
      &     transpiration_in_m, evapo_soil_in_m, root_depth_sl,                      &
      &     infilt, infilt_over, drain_rock, drain_lat_sl,                           &
      &     tpe_overflow, wtr_soilhyd_res, ltpe_closed, ltpe_open)

    ! Aggregate the different runoff flux components
    drain_lat(:) = SUM(drain_lat_sl(:,:),  DIM=2)
    drainage(:)  = drainage(:) + drain_rock(:) + drain_lat(:)
    runoff(:)    = runoff(:) + infilt_over(:)


    ! recalculate soil moisture from the soil moisture of the individual levels
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    !$claw acc data present(runoff,drainage)
    !$claw acc parallel loop gang vector
    DO i = 1, nc
      infilt(i)      = infilt     (i) * rhoh2o / dtime
      runoff(i)      = runoff     (i) * rhoh2o / dtime
      infilt_over(i) = infilt_over(i) * rhoh2o / dtime
      drainage(i)    = drainage   (i) * rhoh2o / dtime
      drain_rock(i)  = drain_rock (i) * rhoh2o / dtime
      drain_lat(i)   = drain_lat  (i) * rhoh2o / dtime
    END DO
    !$claw acc end data
    !$claw acc end data

  END SUBROUTINE calc_soil_hydrology


  !------------------------------------------------------------------------------------------------
  !! @brief   Calculation of subgrid orographic features
  !!
  !! @par Revision History
  !!   Written by Tobias Stacke (MPI-M) in January 2023 based on the
  !!   ARNO Scheme infiltration routine

  SUBROUTINE calc_orographic_features(nc, nlat, oro_stddev, steepness)

    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max

    INTEGER,  INTENT(in) :: &
      & nc, nlat

    REAL(wp), INTENT(in)  :: oro_stddev(:)  !< standard deviation of orography                   [m]
    REAL(wp), INTENT(out) :: steepness(:)   !< parameter defining the subgrid slope distribution [/]

    REAL(wp) ::                       &
      & sigma_0,                      & !< minimum value (100 m); below b = 0.01
      & sigma_max                       !< resolution-dependent maximum
    !
    ! steepness b: shape parameter, defining the sub-gid scale steepness of the orographie:
    !
    ! b = (oro_stddev - sigma_0) / (oro_stddev + sigma_max)
    !
    sigma_0   = oro_var_min
    sigma_max = oro_var_max * 64._wp / REAL(nlat, wp) ! @todo: site level

    ! Compute steepness parameter
    steepness(:) = MAX(0._wp, oro_stddev(:) - sigma_0) / (oro_stddev(:) + sigma_max)
    ! Limit parameter to realistic bounds
    steepness(:) = MAX(MIN(steepness(:), 0.5_wp), 0.01_wp)

  END SUBROUTINE calc_orographic_features

  !------------------------------------------------------------------------------------------------
  !! @brief   Calculation of surface runoff based on the ARNO Scheme
  !! (E. Todini, The ARNO rainfall-runoff model, doi:10.1016/S0022-1694(96)80016-3.)
  !!
  !! @par Revision History
  !!   Revised by Tobias Stacke (MPI-M) in January 2023

  SUBROUTINE arno_scheme(nc, l_infil_subzero, t_soil_sl1, water_to_soil, weq_rootzone, &
    &                    steepness, weq_rootzone_max,                                  &
    &                    runoff, infilt)

    USE mo_jsb_physical_constants, ONLY: tmelt

    INTEGER,  INTENT(in) :: nc

    LOGICAL, INTENT(in)  :: l_infil_subzero

    REAL(wp), INTENT(in) ::           &
      & t_soil_sl1(:),                & !< temperature of uppermost soil layer               [K]
      & water_to_soil(:),             & !< amount of water going into the soil               [m]
      & weq_rootzone(:),              & !< liquid water + ice within the root zone           [m]
      & steepness(:),                 & !< parameter defining the subgrid slope distribution [/]
      & weq_rootzone_max(:)             !< maximum water holding capacity (liquid + ice) of the root zone [m]

    REAL(wp), INTENT(out) ::          &
      & runoff(:),                    & !< surface runoff                                    [m]
      & infilt(:)                       !< infiltration of water into the soil               [m]

    ! Local variables
    INTEGER  ::                       &
      & i                               !< grid cell index

    REAL(wp) ::                       &
      & ws_min_drain,                 & !< minimum amount of root zone soil water for drainage
      & ws_rel,                       & !< relative rootzone soil moisture (water + ice)
      & zb1, zbm, zconw1, zvol

    !$claw acc parallel
    !$claw acc loop gang vector
    DO i = 1, nc

      ! todo: parameters should be defined as parameters
      zb1    = 1._wp + steepness(i)
      zbm    = 1._wp / zb1
      zconw1 = weq_rootzone_max(i) * zb1

      !>  Surface runoff and infiltration
      !   -----------------------------------
      !   f(w) = 1 - (1 - w/w_max)**b

      ! No drainage is possible below a certain minimum amount of water (zwdmin). This corresponds
      ! to the layered bucket scheme, and has nothing to do with field capacity or wilting point.
      ws_min_drain = zwdmin * weq_rootzone_max(i)

      ! Temperature dependence of infiltration should not be included when using the 5-layer snow scheme because:
      !     when tsurf == 0.C and T_soil < 0.C --> all snow melt is going to surface runoff
      IF (.NOT. l_infil_subzero .AND. (t_soil_sl1(i) < tmelt)) THEN
        runoff(i) = water_to_soil(i)       !  no infiltration as uppermost soil layer is frozen -> runoff
        infilt(i) = 0._wp

      ELSE
        IF (water_to_soil(i) > 0._wp .AND. weq_rootzone(i) > ws_min_drain) THEN   ! soilwater content above limit for drainage
          ! compare jsbach3 documentation section 2.3.2.2
          ws_rel = MIN(1._wp, weq_rootzone(i) / weq_rootzone_max(i))              ! relative rootzone soil moisture
          zvol   = (1._wp - ws_rel)**zbm - water_to_soil(i) / zconw1              ! factor taking account subgrid scale
                                                                                  ! inhomogeneity and steepness

          runoff(i) = water_to_soil(i) - (weq_rootzone_max(i) - weq_rootzone(i))  ! > 0: amount of water exceeding soil capacity
          IF (zvol > 0._wp) runoff(i) = runoff(i) + weq_rootzone_max(i) * zvol**zb1
          runoff(i) = MAX(MIN(runoff(i), water_to_soil(i)), 0._wp)                ! no runoff < 0, limited by available water
          infilt(i) = water_to_soil(i) - runoff(i)                                ! water not going into runoff is infiltrated

        ELSE
          runoff(i) = 0._wp                                                       ! no runoff as soil is too dry
          infilt(i) = water_to_soil(i)                                            ! all water infiltrated

        END IF
      END IF

    END DO

  END SUBROUTINE arno_scheme

  !------------------------------------------------------------------------------------------------
  !! @brief   Calculation of vertical soil water movement for multiple soil layers
  !!
  !! @par Revision History
  !!   Routine based on the soil hydrology proposed by Christian Reick (MPI-M) and written
  !!   by Philipp de Vrese (MPI-M)
  !!   Revised by Tobias Stacke (MPI-M) in September 2022

  SUBROUTINE soilhyd(nc, l_fract, lat, lon, nsoil, dtime,          &
    & soilhydmodel, interpol_mean,                                 &
    & w_soil_wilt_fract, dsoil,                                    &
    & wtr_field_cap_sl, wtr_wsat_sl, hyd_cond_sat_sl,              &
    & vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl,            &
    & bclapp_sl, matrix_pot_sl, pore_size_index_sl,                &
    & wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl,             &
    & transpiration_in_m, evapo_soil_in_m, root_depth_sl,          &
    & infilt, infilt_over, drain_rock, drain_lat,                  &
    & tpe_overflow, wtr_soilhyd_res,                               &
    & ltpe_closed, ltpe_open)

    INTEGER,  INTENT(IN)            :: nc                         ! vector length
    LOGICAL,  INTENT(IN)            :: l_fract(:)                 ! Tile fraction > 0
    REAL(wp), INTENT(IN)            :: lat(:)                     ! Grid cell latitudes
    REAL(wp), INTENT(IN)            :: lon(:)                     ! Grid cell longitudes
    INTEGER,  INTENT(IN)            :: nsoil                      ! number of soil layers
    INTEGER,  INTENT(in)            :: soilhydmodel               ! hydrology model scheme
    INTEGER,  INTENT(in)            :: interpol_mean              ! lateral interpolation scheme
    REAL(wp), INTENT(IN)            :: dtime                      ! model time step legth [s]
    REAL(wp), INTENT(in)            :: w_soil_wilt_fract       ! fractional plant wilting point [/]
    REAL(wp), INTENT(IN)            :: dsoil(:,:)                 ! soil depth until bedrock within each layer [m]
    REAL(wp), INTENT(IN)            :: wtr_field_cap_sl(:,:)      ! field capacity of the soil layer (reduced by ice) [m]
    REAL(wp), INTENT(IN)            :: wtr_wsat_sl(:,:)           ! saturation capacity of the soil layer (reduced by ice) [m]
    REAL(wp), INTENT(IN)            :: hyd_cond_sat_sl(:,:)       ! hydraulic conductivity of saturated soil [m/s]
    REAL(wp), INTENT(IN)            :: vol_porosity_sl(:,:)       ! Soil porosity                      - ice not considered [m/m]
    REAL(wp), INTENT(IN)            :: vol_field_cap_sl(:,:)      ! Volumetric field capacity          - ice not considered [m/m]
    REAL(wp), INTENT(IN)            :: vol_p_wilt_sl(:,:)         ! Volumetric permanent wilting point - ice not considered [m/m]
    REAL(wp), INTENT(IN)            :: bclapp_sl(:,:)             ! exponent B in Clapp and Hornberger
    REAL(wp), INTENT(IN)            :: matrix_pot_sl(:,:)         ! matrix potential [m]
    REAL(wp), INTENT(IN)            :: pore_size_index_sl(:,:)    ! soil pore size distribution index used in Van Genuchten
    REAL(wp), INTENT(IN)            :: wtr_soil_pot_scool_sl(:,:) ! potentially supercooled water [m]
    REAL(wp), INTENT(IN)            :: root_depth_sl(:,:)         ! thicknesses of soil layers until rooting depth [m]
    REAL(wp), INTENT(IN)            :: transpiration_in_m(:)      ! amount of transpiration within timestep [m]
    REAL(wp), INTENT(IN)            :: evapo_soil_in_m(:)         ! amount of soil evaporation within timestep [m]
    REAL(wp), INTENT(INOUT)         :: infilt(:)                  ! amount of infiltration within timestep [m]
    REAL(wp), INTENT(INOUT)         :: ice_soil_sl(:,:)           ! amount of ice in the soil [m]
    REAL(wp), INTENT(INOUT)         :: wtr_soil_sl(:,:)           ! soil moisture of each layer [m]
    REAL(wp), INTENT(INOUT)         :: infilt_over(:)             ! amount of infiltration overflow (Horton + Dunne) [m]
    REAL(wp), INTENT(INOUT)         :: drain_rock(:)              ! drainage from lowest soil layer [m]
    REAL(wp), INTENT(INOUT)         :: drain_lat(:,:)             ! overflow from subsurface soil layers [m]
    REAL(wp), INTENT(INOUT)         :: tpe_overflow(:)
    REAL(wp), INTENT(OUT)           :: wtr_soilhyd_res(:)         ! Residual error of the vertical transport scheme

    LOGICAL, INTENT(IN) :: ltpe_closed
    LOGICAL, INTENT(IN) :: ltpe_open

    ! Local variables
    INTEGER  :: i, j                         ! looping index (cells index and soil layer))
    INTEGER  :: nsoil_act(nc)                ! deepest soil layer index (above the bedrock boundary)

    REAL(wp) :: weq_soil(nc)                 ! total column soil moisture (water + ice)
    REAL(wp) :: ws_inter_sl(nc,nsoil)        ! soil moisture with infiltration added to first layer
    REAL(wp) :: ws_vol_sl(nc,nsoil)          ! volumetric soil mositure
    REAL(wp) :: ice_impeda_sl(nc,nsoil)      ! ice impedance factor on all soil layers
    REAL(wp) :: hyd_cond_sl(nc,nsoil)        ! hydraulic conductivity of a layer [m/s]
    REAL(wp) :: diffus_sl(nc,nsoil)          ! diffusivity at mid-depth of the layer[m2/s]
    REAL(wp) :: hyd_cond_li(nc,nsoil)        ! hydraulic conductivity at upper layer interface
    REAL(wp) :: diffus_li(nc,nsoil)          ! diffusivity at upper layer interface
    REAL(wp) :: drain_sl(nc)                 ! drainage from the current layer
    REAL(wp) :: evapo_soil_sl1(nc)           ! bare soil evaporation extracted from top layer
    REAL(wp) :: transpiration_sl(nc,nsoil)   ! transpiration extracted from every soil layer
    REAL(wp) :: deficit_evapotrans(nc)       ! water that could not be extracted from the soil [m]
    REAL(wp) :: remoist(nc)                  ! for TPE
    REAL(wp) :: zmvg_sl(nc,nsoil)            ! exponent M in Van Genuchten scheme = pore_size_index / (pore_size_index + 1)
    REAL(wp) :: wtr_rsdl_sl(nc,nsoil)        ! residual liquid water that can not move through the soil (supercooled or < pwilt_sl)
    REAL(wp) :: weq_wsat_sl(nc,nsoil)        ! Maximum storage in each soil layer (porosity) (not reduced by ice) [m]
    REAL(wp) :: weq_field_cap_sl(nc,nsoil)   ! field capacity of the soil layer (not reduced by ice) [m]
    REAL(wp) :: weq_p_wilt_sl(nc,nsoil)      ! wilting point - not reduced relative to soil ice fraction [m]
    REAL(wp) :: wtr_transp_residual(nc)      ! numerical error caused by the vertical transport scheme [m]
    REAL(wp) :: wtr_transp_corr(nc)          ! uncorrected part of the transport residual [m]
    REAL(wp) :: wtr_exinfcap(nc)             ! surface runoff due to exceeding the soil infiltration capacity [m]
    REAL(wp) :: wtr_exwsat(nc)               ! surface runoff due to exceeding the soil moisture saturation [m]

    !$claw acc data &
    !$claw acc present(dsoil, wtr_field_cap_sl, hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, matrix_pot_sl, &
    !$claw acc         pore_size_index_sl, wtr_soil_sl, tpe_overflow, ice_soil_sl)          &
    !$claw acc pcreate(diffus_sl, ws_tmp, drain_sl, soil_above, wtr_rsdl_sl,      &
    !$claw acc         remoist, zmvg_sl, wtr_transp_residual, wtr_transp_corr, wtr_exinfcap, wtr_exwsat)

    !$claw acc parallel loop gang vector collapse(2)

    DO i = 1, nc
      IF (ANY(wtr_soil_sl(i,:) < 0._wp) .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'negative soil moisture (routine start) at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Soil moisture:    ', wtr_soil_sl(i,:)
        CALL finish ('soilhyd', message_text)
      END IF
    END DO

    !> Calculation of soil saturation, field capacity and wilting point taking into account soil ice
    !
    weq_wsat_sl     (:,:) = vol_porosity_sl (:,:) * dsoil(:,:)
    weq_field_cap_sl(:,:) = vol_field_cap_sl(:,:) * dsoil(:,:)
    weq_p_wilt_sl   (:,:) = vol_p_wilt_sl(:,:)    * dsoil(:,:)

    !> Initialize surface runoff components
    wtr_exinfcap(:) = infilt(:)
    wtr_exwsat(:)   = 0._wp

    DO i = 1, nc
      nsoil_act(i) = 0
      DO j = 1, nsoil
        IF (dsoil(i,j) > 0._wp) nsoil_act(i) = nsoil_act(i) + 1
        zmvg_sl(i,j) = pore_size_index_sl(i,j) / (pore_size_index_sl(i,j)+1._wp)
        ! supercooled water in the soil should not be able to move vertically
        ! @todo at some point p_wilt has to be renamed to something more appropriate, as it
        !       basically describes immobile adhesion water which is not really the same
        !       as the plant wilting point.
        wtr_rsdl_sl(i,j) = MIN(MAX(wtr_soil_pot_scool_sl(i,j), weq_p_wilt_sl(i,j)), wtr_soil_sl(i,j))
      END DO
    END DO

    ws_inter_sl(:,:) = wtr_soil_sl(:,:)
    ws_inter_sl(:,1) = wtr_soil_sl(:,1) + infilt(:)

    ! Determining hydraulic conductivity and diffusivity on layers and at layer interfaces
    !   Note that the indexing uses the upper interface -> interface Ii = I(i,i-1)
    !   For numerical reasons, soil loops always cover all soil layers, even if they are not active (bedrock)
    CALL get_soilhyd_cond_properties(                                       &
          & soilhydmodel, interpol_mean,                                    &
          & dtime,                                                          &
          & nc,nsoil,                                                       &
          & nsoil_act(:), dsoil(:,:),                                       &
          & ws_inter_sl(:,:), ice_soil_sl(:,:), wtr_soil_pot_scool_sl(:,:), &
          & weq_wsat_sl(:,:), weq_field_cap_sl(:,:), weq_p_wilt_sl(:,:),    &
          & hyd_cond_sat_sl(:,:), matrix_pot_sl(:,:),                       &
          & bclapp_sl(:,:), zmvg_sl(:,:),                                   &
          & ice_impeda_sl(:,:),                                             &
          & hyd_cond_sl(:,:), diffus_sl(:,:),                               &
          & hyd_cond_li(:,:), diffus_li(:,:)                                &
          )

    !< Diagnose evapotranspiration fluxes from the different soil layer based on the available
    !    mobile water
    CALL diagnose_evapotrans(nc, l_fract, lat, lon, nsoil, w_soil_wilt_fract, root_depth_sl, &
      &                      dsoil, wtr_field_cap_sl, transpiration_in_m, evapo_soil_in_m,   &
      &                      wtr_soil_sl - wtr_rsdl_sl, evapo_soil_sl1, transpiration_sl,    &
      &                      deficit_evapotrans)

    DO i = 1, nc
      ! Diagnose infiltration from the hydraulic conductivity of the top soil layer [m s-1 -->] [m]
      ! @todo this is implicitely done already within the ARNO scheme, thus it cannot be re-applied here.
      !       As soon as the Arno Scheme is retired and the different soil properties associated with
      !       orography are explizitly resolved, infiltration could be computed based on the conductivity instead.
      !       However, that with the current model configuration this approach either results in too high
      !       infiltration (if top layer conductivity is not constraint with the amout of water + infiltration
      !       of the top layer) or too low if the soil is not fully saturated. Smart ideas are needed here!
      ! infilt(i)       = MIN(infilt(i), hyd_cond_sl(i,1) * dtime)

      ! Any water exceeding the infiltration capacity is given to the surface runoff (Horton runoff)
      !   Currently zero, as it is already included in the surface runoff from the ARNO scheme
      wtr_exinfcap(i) = MAX(wtr_exinfcap(i) - infilt(i), 0._wp)

      ! Prescribe drainage from hydraulic conductivity of the lowest layer [m s-1] --> [m]
      drain_rock(i) = hyd_cond_sl(i,nsoil_act(i)) * dtime

      !< store actual soil water state, infiltration (upper boundary condition) and evaporation 
      !  to compute the residual transport error
      wtr_transp_residual(i) = SUM(wtr_soil_sl(i,:)) + infilt(i) + evapo_soil_sl1(i) + SUM(transpiration_sl(i,:))

      ! Determination of Vertical movement of mobile water through the soil.
      !   Note that normalized storage values (volumetric soil moisture and storage change) are
      !   needed for the transport and immobile water is temporarily removed from the soil.
      DO j = 1, nsoil
        IF (dsoil(i,j) > 0._wp) THEN
          ws_vol_sl(i,j) = (wtr_soil_sl(i,j) - wtr_rsdl_sl(i,j)) / dsoil(i,j)
          transpiration_sl(i,j) = transpiration_sl(i,j) / dsoil(i,j)
        ELSE
          ! without soil there is no change within the soil and the values do not matter
          ws_vol_sl(i,j) = wtr_soil_sl(i,j) - wtr_rsdl_sl(i,j)
        END IF
      END DO

      ! CALL calc_vertical_transport(nc, n, d, alpha, dzf, top_bound, bot_bound, S, C, P, Q, X_old)
      !   *-bound, S are fluxes i.e. x/sec
      ws_vol_sl(i,1:nsoil_act(i))  = calc_vertical_transport(           &
                         & nsoil_act(i), dtime, 1._wp,                  &
                         & dsoil(i,1:nsoil_act(i)),                     & ! Layer thickness
                         & (infilt(i) + evapo_soil_sl1(i)) / dtime,     & ! Upper boundary condition
                         & drain_rock(i) / dtime,                       & ! Lower boundary condition
                         & transpiration_sl(i, 1:nsoil_act(i)) / dtime, & ! Normalized storage change term (Transpiration)
                         & SPREAD(1._wp, DIM=1, ncopies=nsoil_act(i)),  &
                         & diffus_li(i,1:nsoil_act(i)),                 & ! Diffusivity
                         & hyd_cond_li(i,1:nsoil_act(i)),               & ! Hydraulic conductivity
                         & ws_vol_sl(i,1:nsoil_act(i))                  & ! Normalized soil moisture
                         & )

      ! Converting back to absolute water content
      DO j = 1, nsoil
        IF (dsoil(i,j) > 0._wp) THEN
          wtr_soil_sl(i,j) = ws_vol_sl(i,j) * dsoil(i,j) + wtr_rsdl_sl(i,j)
        ELSE
          wtr_soil_sl(i,j) = ws_vol_sl(i,j) + wtr_rsdl_sl(i,j)
        END IF
      END DO

      !< If lowest soil layer becomes negative, reduce bottom drainage
      !    In rare cases - when the bottom layer is very thin - also above layers might get negative.
      !    Generally, this only happens right after initialization.
      IF (wtr_soil_sl(i, nsoil_act(i)) < 0._wp .AND. l_fract(i)) THEN
        DO j = nsoil_act(i), 1, -1
          IF (wtr_soil_sl(i,j) < 0._wp) THEN
            WRITE (message_text,*) 'Soil moisture correction needed at ', &
              &  lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
              & 'Soil moisture:   ', wtr_soil_sl(i,:), NEW_LINE('a'), &
              & 'Bottom drainage: ', drain_rock(i)
            CALL message('soilhyd', message_text, all_print=.TRUE.)
            drain_rock(i) = MAX((drain_rock(i) + wtr_soil_sl(i,j)), 0._wp)
            wtr_soil_sl(i,j) = 0._wp
            WRITE (message_text,*) 'After correction: wtr_soil_sl(i,:): ', &
              & lat(i),'N and ',lon(i), 'E: ', NEW_LINE('a'), &
              & 'Soil moisture:   ', wtr_soil_sl(i,:), NEW_LINE('a'), &
              & 'Bottom drainage: ', drain_rock(i)
            CALL message('soilhyd', message_text, all_print=.TRUE.)
          END IF
        END DO
      END IF

      !< subtract new soil state and bottom drainage (lower boundary condition) from water balance
      wtr_transp_residual(i) = wtr_transp_residual(i) - (SUM(wtr_soil_sl(i,:)) + drain_rock(i))
      wtr_soilhyd_res(i) = wtr_transp_residual(i)

    END DO

    !< Correct for the vertical transport error.
    DO i = 1, nc
      !< 1. Complain if the overall transport error is too large
      IF (ABS(wtr_transp_residual(i)) > 1.0e-10_wp .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'Soil water transport residual too large at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Transport residual: ', wtr_transp_residual(i), NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(i), NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(i), NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_rock(i), NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(i,:), NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(i,:)
        CALL finish ('soilhyd', message_text)
      END IF

      !< 2. Modify lower boundary condition (drainage) as much as possible
      IF (ABS(wtr_transp_residual(i)) > 0._wp) THEN
        wtr_transp_corr(i)     = MIN(wtr_transp_residual(i) + drain_rock(i), 0._wp)
        drain_rock(i)          = MAX(drain_rock(i) + wtr_transp_residual(i), 0._wp)
        wtr_transp_residual(i) = wtr_transp_corr(i)
      END IF

      !< 3. Modify soil layers if necessary
      DO j = nsoil, 1, -1
        IF (ABS(wtr_transp_residual(i)) > 0._wp) THEN
          ! Make sure to correct the residual and not to put negative soil moisture onto the residual
          IF (wtr_soil_sl(i,j) >= 0._wp .OR. wtr_soil_sl(i,j) + wtr_transp_residual(i) > 0._wp) THEN
            wtr_transp_corr(i)     = MIN(wtr_transp_residual(i) + wtr_soil_sl(i,j), 0._wp)
            wtr_soil_sl(i,j)       = MAX(wtr_soil_sl(i,j) + wtr_transp_residual(i), 0._wp)
            wtr_transp_residual(i) = wtr_transp_corr(i)
          END IF
        END IF
      END DO

      !> 4. Complain if correction is not sufficient
      IF (ABS(wtr_transp_residual(i)) > 1.0e-10_wp .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'Cannot correct soil water transport residual at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Transport residual: ', wtr_transp_residual(i), NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(i), NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(i), NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_rock(i), NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(i,:), NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(i,:)
        CALL finish ('soilhyd', message_text)
      END IF
    END DO

    ! Correct remaining evapotranspiration deficit
    ! @todo this is an unfortunate necessity as sometimes the evaporative demand is not fully
    !       backed by the water available in the land surface reservoirs.
    !       Changing the soil ice content here obviously violates the energy balance. A correction
    !       term is needed for this, if no better solution can be found to avoid the evaporation deficit.
    DO i = 1, nc
      weq_soil(i) = SUM(wtr_soil_sl(i,:)) + SUM(ice_soil_sl(i,:))  ! water and ice within the column
      DO j = 1, nsoil
        IF (deficit_evapotrans(i) < 0._wp .AND. weq_soil(i) > 0._wp) THEN
          ! reduce soil moisture and ice relative to the water content of each layer
          wtr_soil_sl(i,j) = MAX(0._wp, wtr_soil_sl(i,j) + deficit_evapotrans(i) * wtr_soil_sl(i,j) / weq_soil(i))
          ice_soil_sl(i,j) = MAX(0._wp, ice_soil_sl(i,j) + deficit_evapotrans(i) * ice_soil_sl(i,j) / weq_soil(i))
        END IF
      END DO
      ! @todo even after the correction, there are still cells with zero soil moisture and a negative
      !       deficit_evapotrans left. However, these do not seem to violate the overall water balance.
      !       Are these cells which should not be included in the soil hydrology computation?
    END DO

    ! Additional corrections for excess soil moisture
    DO i = 1, nc
      ! In case of over-saturation water would not have infiltrated in the first place --> pile
      ! water upwards
      DO j = nsoil, 1, -1
        IF (j > 1) THEN
          ! If saturation moisture is exceeded, add the excess water to next higher soil layer
          IF (wtr_soil_sl(i,j) > wtr_wsat_sl(i,j)) THEN
            wtr_soil_sl(i,j-1) = wtr_soil_sl(i,j-1) + (wtr_soil_sl(i,j) - wtr_wsat_sl(i,j))
            wtr_soil_sl(i,j)   = wtr_wsat_sl(i,j)
          END IF
        ELSE IF (j == 1) THEN
          ! Top layer saturation excess is added to the surface runoff (Dunne runoff)
          IF (wtr_soil_sl(i,j) > wtr_wsat_sl(i,j)) THEN
            wtr_exwsat(i)      = wtr_soil_sl(i,j) - wtr_wsat_sl(i,j)
            wtr_soil_sl(i,j)   = wtr_wsat_sl(i,j)
          END IF
        END IF
      END DO
      !
      ! Lateral Drainage (i.e. drainage from upper soil layers)
      ! ---- In supersaturated (or frozen) soils the water above the field capacity is assumed to 
      !      move lateraly until it reaches the stream (or the neighbouring topo-tile).
      !      In this implementation all water above field capacity is added to the drainage. However,
      !      this is a very simplistic assumption. For e.g. high resolution simulations it might
      !      be worthwile to consider lateral drainage as a reservoir outflow with a predefined
      !      lagtime to allow for oversaturation of soils.
      ! @todo In the current model configuration lateral drainage seems to deteriorate river
      !      discharge while it has little effect on evapotranspiration. Therefore it is 
      !      disabled, but this should be revised, e.g. for high resolution simulations.
      ! DO j = nsoil, 1, -1
      !   ! If maximum moisture is exceeded, add the excess water to lateral drainage
      !   IF (wtr_soil_sl(i,j) > wtr_field_cap_sl(i,j)) THEN
      !     drain_sl(i)      = wtr_soil_sl(i,j) - wtr_field_cap_sl(i,j) ! Put excess water in lateral drainage
      !     wtr_soil_sl(i,j) = wtr_soil_sl(i,j) - drain_sl(i)           ! remove lateral drainage of the layer from soil moisture
      !     IF (.NOT. ltpe_closed) THEN
      !       ! Add excess water to lateral drainage
      !       drain_lat(i,j)   = drain_lat(i,j)  + drain_sl(i)
      !     ELSE
      !       ! no drainage for TPE closed case --> redirect overflow from drainage to overflow pool
      !       tpe_overflow(i) = tpe_overflow(i) + drain_sl(i)
      !     END IF
      !   ELSE IF (ltpe_closed) THEN
      !     ! refill soil moisture from overflow pool
      !     remoist(i)       = MIN(tpe_overflow(i), wtr_field_cap_sl(i,j) - wtr_soil_sl(i,j))
      !     wtr_soil_sl(i,j) = wtr_soil_sl(i,j) + remoist(i)
      !     tpe_overflow(i)  = tpe_overflow(i) - remoist(i)
      !   END IF
      ! END DO

      ! Modification for TPE (open)
      IF (ltpe_open) THEN
        !$claw acc loop seq
        DO j = 1, nsoil
          !IF (i == nsoil - 2 .OR. i == nsoil - 1 .OR. i == nsoil) THEN   ! special adaptation of Shirisha
          IF (i == nsoil-1 .OR. i == nsoil) THEN
            IF (dsoil(i,j) > 0._wp) THEN
             ! keep the two lowest soil layers always very wet
             wtr_soil_sl(i,j) = MAX(wtr_soil_sl(i,j),wtr_field_cap_sl(i,j) * 0.9_wp)
            END IF
          END IF
        END DO
      END IF

    END DO

    ! Add surface runoff components as infiltration overflow
    infilt_over(:) = wtr_exinfcap(:) + wtr_exwsat(:)

    ! Check value range for soil moisture
    DO i = 1, nc
      IF (ANY(wtr_soil_sl(i,:) < -EPSILON(1.0_wp)) .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'negative soil moisture (routine end) at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Infiltration:     ', infilt(i), NEW_LINE('a'), &
          & 'Soil evaporation: ', evapo_soil_sl1(i), NEW_LINE('a'), &
          & 'Bottom drainage:  ', drain_rock(i), NEW_LINE('a'), &
          & 'Transpiration:    ', transpiration_sl(i,:), NEW_LINE('a'), &
          & 'Soil moisture:    ', wtr_soil_sl(i,:)
        CALL finish ('soilhyd', message_text)
      ELSE IF (ANY(wtr_soil_sl(i,:) < 0.0_wp) .AND. l_fract(i)) THEN
        ! WRITE (message_text,*) 'setting negative soil moisture (routine end) to zero at ', &
        !   & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
        !   & 'Soil moisture:    ', wtr_soil_sl(i,:)
        ! Enable for debugging purpose, otherwise it fills up the logfile
        ! CALL message ('soilhyd', message_text)
        wtr_soil_sl(i,:) = MAX(wtr_soil_sl(i,:), 0.0_wp)
      END IF

      IF (ANY(wtr_soil_sl(i,:) > wtr_wsat_sl(i,:) + EPSILON(1.0_wp)) .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'soil moisture exceeds saturation capacity(routine end) at ', &
          &  lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(i), NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(i), NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_rock(i), NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(i,:), NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(i,:), NEW_LINE('a'), &
          & 'Saturation capcity: ', wtr_wsat_sl(i,:)
        CALL finish ('soilhyd', message_text)
      END IF
    END DO

    !$claw acc end parallel
    !$claw acc end data

  END SUBROUTINE soilhyd

  !----------------------------------------------------------------------------------------------

  SUBROUTINE diagnose_evapotrans( &
    & nc,                         &
    & l_fract,                    &
    & lat,                        &
    & lon,                        &
    & nsoil,                      &
    & w_soil_wilt_fract,          &
    & droot,                      &
    & dsoil,                      &
    & wtr_field_cap_sl,           &
    & transpiration_in_m,         &
    & evapo_soil_in_m,            &
    & wtr_soil_sl,                &
    & evapo_soil_sl1,             &
    & transpiration_sl,           &
    & deficit_evapotrans          &
    & )

    !----------------------------------------------------------------------------------------------
    !< This routine distributes the evaporative demand of soil evaporation and transpiration
    !! between the soil layers of the root zone. Excess evaporation that cannot be satisfied
    !! from the root zone is subtracted from deeper layers.
    !! Note, that the routine only diagnoses the layered evaporation fluxes but the actual 
    !! soil moisture update is done in the vertical transport routine. Also, it only consideres
    !! negative evapotranspiration fluxes (actual evapotranspiration) as positive fluxes (dew)
    !! are already considered in calc_surface_hydrology.
    !! The routine is based on routine soilchange by Stefan Hagemann.
    !----------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nc                      ! vector length
    LOGICAL,  INTENT(IN)     :: l_fract(:)              ! Tile fraction > 0
    REAL(wp), INTENT(IN)     :: lat(:)                  ! Grid cell latitudes
    REAL(wp), INTENT(IN)     :: lon(:)                  ! Grid cell longitudes
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    REAL(wp), INTENT(in)     :: w_soil_wilt_fract       ! fractional plant wilting point [/]
    REAL(wp), INTENT(in)     :: droot(:,:)              ! root depth within layer [m]
    REAL(wp), INTENT(in)     :: dsoil(:,:)              ! soil depth until bedrock within layer [m]
    REAL(wp), INTENT(in)     :: wtr_field_cap_sl(:,:)   ! field capacity of the layer (reduced by ice)
    REAL(wp), INTENT(in)     :: transpiration_in_m(:)   ! transpiration [m]
    REAL(wp), INTENT(in)     :: evapo_soil_in_m(:)      ! bare soil evaporation [m]

    REAL(wp), INTENT(in)     :: wtr_soil_sl(:,:)        ! water content of the soil layer [m]

    REAL(wp), INTENT(OUT)    :: evapo_soil_sl1(:)       ! Soil evaporation from the top soil layer [m]
    REAL(wp), INTENT(OUT)    :: transpiration_sl(:,:)   ! Transpiration from soil layers [m]
    REAL(wp), INTENT(OUT)    :: deficit_evapotrans(:)   ! Water deficit due to evapotranspiration [m]

    ! local variables

    INTEGER  :: i, j
    REAL(wp) :: dummy_wtr_soil(nc, nsoil) ! dummy soil moisture storage [m]
    REAL(wp) :: deficit_evapo(nc)         ! water deficit due to evaporation [m]
    REAL(wp) :: deficit_trans(nc)         ! water deficit due to tranpiration [m]
    REAL(wp) :: deficit(nc)               ! water deficit within the actual soil layer [m]
    REAL(wp) :: remaining(nc)             ! water remaining in the soil layer [m]
    REAL(wp) :: rootfract(nc)             ! fraction of the soil layer within the root zone [/]
    REAL(wp) :: fixed(nc)                 ! amount of water below the wilting point, not available for plants [m]
    REAL(wp) :: root_depth(nc)            ! Total depth of rootzone [m]

    !----------------------------------------------------------------------------------------------
    !  diagnose evaporation and transpiration fluxes
    !----------------------------------------------------------------------------------------------

    !$claw acc data present(droot, dsoil, wtr_field_cap_sl, transpiration_in_m, evapo_soil_in_m, &
    !$claw acc wtr_soil_sl, evapo_soil_sl1, transpiration_sl, deficit_evapotrans) pcreate(deficit_evapo, &
    !$claw acc deficit_trans, deficit_layer, remaining, rootfract, fixed, root_depth, dummy_wtr_soil)

    !  bare soil evaporation
    ! -----------------------
    !$claw acc parallel
    !$claw acc loop gang vector


    DO i = 1, nc
      IF (ANY(wtr_soil_sl(i,:) < 0._wp) .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'negative soil moisture (routine start) at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Soil moisture: ', wtr_soil_sl(i,:)
        CALL finish ('diagnose_evapotrans', message_text)
      END IF
    END DO

    ! Initialize dummy soil moisture storage, because the real soil moisture storage
    !   is only to be changed during the vertical transport
    dummy_wtr_soil(:,:)   = wtr_soil_sl(:,:)
    evapo_soil_sl1(:)     = 0._wp
    transpiration_sl(:,:) = 0._wp

    DO i = 1, nc

      ! 1. Soil evaporation
      ! ---------------

      deficit_evapo(i) = 0._wp                                         ! water missing for evaporation due to limited soil water
      IF (evapo_soil_in_m(i) < 0._wp) THEN                             ! actual evaporation flux (not dew)
        dummy_wtr_soil(i,1) = dummy_wtr_soil(i,1) + evapo_soil_in_m(i) ! take water from uppermost soil layer
        IF (dummy_wtr_soil(i,1) < 0._wp) THEN                          ! if there is not enough soil water
          deficit_evapo(i) = dummy_wtr_soil(i,1)                       ! this is the amount of water missing
          dummy_wtr_soil(i,1) = 0._wp                                  ! soil moisture cannot be negative
        END IF
        evapo_soil_sl1(i) = evapo_soil_in_m(i) - deficit_evapo(i)      ! reduce soil evaporation if top layer soil is too dry
      END IF

      !  2. Transpiration
      ! ---------------

      root_depth(i) = SUM(droot(i,:))

      ! sealed grid cells (rooting_depth zero) and actual transpiration
      !    => no water reachable for transpiration, add flux to deficit term
      ! @todo what does rooting depth=zero for bare soil tile?
      IF(root_depth(i) <= 0._wp .AND. transpiration_in_m(i) < 0._wp) THEN
        deficit_trans(i) = transpiration_in_m(i)                       ! add transpiration to deficit term
      ELSE
        deficit_trans(i) = 0._wp
      END IF

      !$claw acc loop seq
      DO j =1, nsoil

        deficit(i)            = 0._wp
        rootfract(i)          = 0._wp
        fixed(i)              = 0._wp
        transpiration_sl(i,j) = 0._wp
        remaining(i)          = 0._wp

        ! grid cell not sealed and actual transpiration (not dew) takes place
        !    => take water equally from all over the root zone
        IF (root_depth(i) > 0._wp .AND. transpiration_in_m(i) < 0._wp) THEN

          ! water needed for transpiration from the current layer
          transpiration_sl(i,j) = transpiration_in_m(i) * droot(i,j) / root_depth(i)

          IF (droot(i,j) > 0._wp) THEN

            rootfract(i) = droot(i,j) / dsoil(i,j) ! fraction of the soil layer that is within the root zone
            ! @todo at some point wtr_field_cap * w_soil_wilt_fract needs to be changed to the real wilting point moisture of the
            !       specific plant growing on this tile
            fixed(i) = wtr_field_cap_sl(i,j) * w_soil_wilt_fract * rootfract(i)         ! amount of unavailable root zone water

            IF (dummy_wtr_soil(i,j) * rootfract(i) >= fixed(i)) THEN                    ! there is water in the root zone available
              remaining(i) = dummy_wtr_soil(i,j) * rootfract(i) + transpiration_sl(i,j) ! remaining water after transpiration
              IF (remaining(i) < fixed(i)) THEN                                         ! transpiration exceeds available water
                deficit(i) = remaining(i) - fixed(i)
                dummy_wtr_soil(i,j) = fixed(i) + dummy_wtr_soil(i,j) * (1._wp-rootfract(i))
              ELSE                                                                      ! enough soil water
                deficit(i) = 0._wp
                dummy_wtr_soil(i,j) = remaining(i) + dummy_wtr_soil(i,j) * (1._wp-rootfract(i))
              END IF
            ELSE                                                                        ! not enough water for transpiration
              deficit(i) = transpiration_sl(i,j)
            END IF
            transpiration_sl(i,j) = transpiration_sl(i,j) - deficit(i)                  ! reduce transpiration if necessary

            ! sum up deficits from the different layers
            deficit_trans(i) = deficit_trans(i) + deficit(i)

          END IF
        END IF
      END DO

      ! -----------------------------------------
      !  Total deficit from evapotranspiration
      ! -----------------------------------------
      deficit_evapotrans(i) = deficit_evapo(i) + deficit_trans(i)

      ! If the deficit can be reduced by redistributing transpiration extraction from different layers,
      !   it should be done already here. The deficit should be taken from layers starting at the
      !   surface to accelerate the feedback on evapotranspiration
      DO j = 1, nsoil
        IF (deficit_evapotrans(i) < 0._wp .AND. dummy_wtr_soil(i,j) > wtr_field_cap_sl(i,j) * w_soil_wilt_fract) THEN
          ! add water that can be still transpired to the transpiration to reduce the overall deficit
          remaining(i)          = dummy_wtr_soil(i,j) - &
            & MAX(dummy_wtr_soil(i,j) + deficit_evapotrans(i), wtr_field_cap_sl(i,j) * w_soil_wilt_fract)
          deficit_evapotrans(i) = deficit_evapotrans(i) + remaining(i)
          dummy_wtr_soil(i,j)   = dummy_wtr_soil(i,j)  - remaining(i)
          transpiration_sl(i,j) = transpiration_sl(i,j) - remaining(i)
        END IF
      END DO

    END DO ! nc loop

    !$claw acc end parallel
    !$claw acc end data

#ifndef _OPENACC
    ! make sure there is no negative soil moisture after applying ET
    DO i = 1, nc
      IF (ANY(dummy_wtr_soil(i,:) < 0._wp) .AND. l_fract(i)) THEN
        WRITE (message_text,*) 'negative soil moisture (routine end) at ', &
          & lat(i),'N and ',lon(i),'E:', NEW_LINE('a'), &
          & 'Soil evaporation: ', evapo_soil_sl1(i), NEW_LINE('a'), &
          & 'Transpiration:    ', transpiration_sl(i,:), NEW_LINE('a'), &
          & 'ET deficit:       ', deficit_evapotrans(i), NEW_LINE('a'), &
          & 'Soil moisture:    ', dummy_wtr_soil(i,:)
        CALL finish ('diagnose_evapotrans', message_text)
      END IF
    END DO
#endif

  END SUBROUTINE diagnose_evapotrans

  !----------------------------------------------------------------------------------------------
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

  FUNCTION calc_vertical_transport(n, dt, alpha, dzf, top_bound, bot_bound, S, C, P, Q, X_old) &
    &      RESULT(X_new)

  ! INPUT
  INTEGER, INTENT(in) :: n                          ! number of active soil levels
  REAL(wp), Intent(in):: dt, alpha,            &    ! time step
                       & top_bound, bot_bound, &    ! fluxes at upper and lower boundary
                       & dzf(n),               &    ! soil layer thickness 
                       & C(n),                 &    ! heat capacity if used for heat transport,
                                                    ! set to 1 for water transport
                       & S(n),                 &    ! prescribed water/heat change flux per layer
                       & Q(n),                 &    ! Diffusivity at (upper) layer boundaries
                       & P(n),                 &    ! Conductivity at (upper) layer boundaries
                       & X_old(n)                   ! Normalized soil state for every layer

  ! OUTPUT
  REAL(wp)           :: X_new(n)

  ! LOCAL
  INTEGER            :: j
  REAL(wp)           :: dzh(n), matrix(n,4), &
                      & bet, gam(n)

  dzh(:) = 0.0                             ! thickness of half levels
  DO j=2,n
    dzh(j) = (dzf(j-1) + dzf(j)) / 2._wp
  END DO

  ! Note that P & Q refer to the inter face above a given lvl not below.
  ! Hence for j == 1, there is no P(j) and Q(j)
  ! Make sure that this is consistent with the definition of input params.

  DO j = 1, n
    IF(j .eq. 1) THEN
      ! Top soil layer
      matrix(j,1) =  0.0_wp
      matrix(j,2) =  P(j+1) / dzh(j+1) + dzf(j) * C(j) / (alpha * dt)
      matrix(j,3) = -P(j+1) / dzh(j+1)
      matrix(j,4) = -Q(j+1) + dzf(j) * S(j) + dzf(j) * C(j) * X_old(j) / (alpha * dt) + top_bound
    ELSE IF (j .eq. n) THEN
      ! Bottom soil layers
      matrix(j,1) = -P(j) / dzh(j)
      matrix(j,2) =  P(j) / dzh(j) + dzf(j) * C(j) / (alpha * dt)
      matrix(j,3) =  0.0_wp
      matrix(j,4) =  Q(j)          + dzf(j) * S(j) + dzf(j) * C(j) * X_old(j) / (alpha * dt) - bot_bound
    ELSE
      ! Middle soil layers
      matrix(j,1) = -P(j) / dzh(j)
      matrix(j,2) =  P(j) / dzh(j) +  P(j+1) / dzh(j+1) + dzf(j) * C(j) / (alpha * dt)
      matrix(j,3) = -P(j+1) / dzh(j+1)
      matrix(j,4) =  Q(j) - Q(j+1) + dzf(j) * S(j) + dzf(j) * C(j) * X_old(j) / (alpha * dt)
    END IF
  END DO

  ! This solver uses LU-decomposition to solve the tridiagonal matrix.
  ! It is based on the routine tridiag (page 43) in the
  ! Numerical Recipes in Fortran. 2nd Edition, 1992, ISBN 0-521-43064-X

  IF( matrix(1,2) .EQ. 0.0) THEN
    WRITE (message_text,*) 'something went terribly wrong ... rewrite equations!', C(1), dzf(1)
    CALL finish('calc_vertical_transport', message_text)
  END IF
  bet  = matrix(1,2)
  X_new(1) = matrix(1,4)/bet
  DO j = 2, n
    gam(j) = matrix(j-1,3) /bet
    bet    = matrix(j,2) - matrix(j,1) * gam(j)
    IF( bet .EQ. 0.0) THEN
      WRITE (message_text,*) 'tridiag: matrix couldn t be solved !'
      CALL finish('calc_vertical_transport', message_text)
    END IF
    X_new(j)   = (matrix(j,4) - matrix(j,1) * X_new(j-1)) / bet
  END DO
  DO j  = n-1, 1, -1
    X_new(j)= X_new(j) - gam(j+1) * X_new(j+1)
  END DO

  END FUNCTION calc_vertical_transport

  SUBROUTINE get_soilhyd_cond_properties(           &
            & soilhydmodel, interpol_mean,          &
            & dt,                                   &
            & nc, nsoil, nsoil_act, dsoil,          &
            & ws, is, wscools,                      &
            & wsat, fc, pwp,                        &
            & k_sat, mpot_sat, bclapp, mvgn,        &
            & ice_impedance, K, D, K_inter, D_inter &
            )

    USE mo_hydro_constants,    ONLY: BrooksCorey_, ClappHornberger_, VanGenuchten_, &
      &                              Upstream_, Arithmetic_

    INTEGER,  INTENT(in)            :: soilhydmodel       !> Approaches to determine conductivites & diffusivities (KD)
    INTEGER,  INTENT(in)            :: interpol_mean      !> and how to derive values at layer interface (AVG scheme)
    REAL(wp), INTENT(in)            :: dt                 !> time step length
    INTEGER,  INTENT(in)            :: nc                 !> vector length
    INTEGER,  INTENT(in)            :: nsoil              !> number of below ground layers
    INTEGER,  INTENT(in)            :: nsoil_act(:)       !> index of deepest soil layers (above bedrock)
    REAL(wp), INTENT(in)            :: dsoil(:,:)         !> soil depth (until bedrock) per layer
    REAL(wp), INTENT(in)            :: ws(:,:)            !> water content
    REAL(wp), INTENT(in)            :: is(:,:)            !> ice content
    REAL(wp), INTENT(in)            :: wscools(:,:)       !> supercooled water
    REAL(wp), INTENT(in)            :: wsat(:,:)          !> soil depth (until bedrock) per layer
    REAL(wp), INTENT(in)            :: fc(:,:)            !> field capacity [m]
    REAL(wp), INTENT(in)            :: pwp(:,:)           !> residual water (until bedrock) per layer
    REAL(wp), INTENT(in)            :: k_sat(:,:)         !> hydraaulic conductivity
    REAL(wp), INTENT(in)            :: mpot_sat(:,:)      !> matric potential
    REAL(wp), INTENT(in)            :: bclapp(:,:)        !> clapp & Hornberger exponent
    REAL(wp), INTENT(in)            :: mvgn(:,:)          !> Van Genuchten m
    REAL(wp), INTENT(out), OPTIONAL :: ice_impedance(:,:) !> Impedance factor to account for ice blocking flowpaths
    REAL(wp), INTENT(out), OPTIONAL :: K(:,:)             !> hydraulic condctivity of the soil layer [m]
    REAL(wp), INTENT(out), OPTIONAL :: D(:,:)             !> diffusivity of the soil layer [m]
    REAL(wp), INTENT(out), OPTIONAL :: K_inter(:,:)       !> hydraulic condctivity at (upper) layer interface [m]
    REAL(wp), INTENT(out), OPTIONAL :: D_inter(:,:)       !> diffusivity at (upper) layer interface [m]

    ! local variables
    INTEGER  :: i,j ! cell index and soil layer index
    REAL(wp) :: ws_vol(nc,nsoil), ice_vol(nc,nsoil), scool_vol(nc,nsoil), &
             &  ws_rel_vol(nc,nsoil),  ws_range_vol(nc,nsoil),  &
             &  wsat_vol(nc,nsoil), fc_vol(nc,nsoil), pwp_vol(nc,nsoil), &
             &  ws_free_vol(nc,nsoil), ice_imp(nc,nsoil)


    ! Set defaults
    ws_rel_vol(:,:) = 0._wp
    ice_imp(:,:)    = 1._wp

    IF (PRESENT(K))       K(:,:)       = 0._wp
    IF (PRESENT(D))       D(:,:)       = 0._wp
    IF (PRESENT(K_inter)) K_inter(:,:) = 0._wp
    IF (PRESENT(D_inter)) D_inter(:,:) = 0._wp

    ! Convert to volumetric quantities
    WHERE(dsoil(:,:)> 0._wp)
      ws_vol(:,:)   =  ws(:,:)      / dsoil(:,:)
      ice_vol(:,:)  =  is(:,:)      / dsoil(:,:)
      scool_vol     =  wscools(:,:) / dsoil(:,:)
      wsat_vol(:,:) =  wsat(:,:)    / dsoil(:,:)
      fc_vol(:,:)   =  fc(:,:)      / dsoil(:,:)
      pwp_vol(:,:)  =  pwp(:,:)     / dsoil(:,:)
    ELSEWHERE
      ws_vol(:,:)   =  0._wp
      ice_vol(:,:)  =  0._wp
      scool_vol     =  0._wp
      wsat_vol(:,:) =  0._wp
      fc_vol(:,:)   =  0._wp
      pwp_vol(:,:)  =  0._wp
    END WHERE

    ! Compute impedance factor due to soil ice
    DO i = 1, nc
      DO j = 1, nsoil
        IF (j <= nsoil_act(i)) THEN
          IF (fc_vol(i,j) + ice_vol(i,j) > 1.0e-10_wp) THEN
            ice_imp(i,j)  = 10._wp**(-6._wp * (ice_vol(i,j) + scool_vol(i,j) ) &
              & / (fc_vol(i,j) + ice_vol(i,j)))
          ELSE
            ice_imp(i,j) = 1._wp
          END IF
        END IF
      END DO
    END DO

    IF (PRESENT(ice_impedance)) THEN
      ice_impedance(:,:) = ice_imp(:,:)
    END IF

    ! Compute relative soil moisture depending on soil hydrological model
    IF (soilhydmodel == ClappHornberger_) THEN
      ws_free_vol(:,:)   = MAX(0._wp, ws_vol(:,:) - scool_vol(:,:))
      ws_range_vol(:,:)  = wsat_vol(:,:)
    ELSE  ! Brooks & Corey; Van Genuchten
      ws_free_vol(:,:)   = MAX(0._wp, ws_vol(:,:) - MAX(scool_vol(:,:), pwp_vol(:,:)))
      ws_range_vol(:,:)  = fc_vol(:,:) - pwp_vol(:,:)
    END IF
    WHERE (ws_range_vol(:,:) > 0._wp)
      ws_rel_vol(:,:)    = ws_free_vol(:,:) / ws_range_vol(:,:)
    END WHERE
    ! Making sure ws_rel is in the physical range
    ! This is especially important when infiltration is added to the first layer
    IF (soilhydmodel == VanGenuchten_) THEN
      ws_rel_vol(:,:) = MIN(0.9999_wp, MAX(0._wp, ws_rel_vol(:,:))) ! doesn't work for wrel == 1; (1-ws_rel(:,:)**(1/mvgn(:,:)))**(-mvgn)
    ELSE  ! Brooks & Corey; Van Genuchten
      ws_rel_vol(:,:) = MIN(1._wp,     MAX(0._wp, ws_rel_vol(:,:)))
    END IF

    ! Compute hydrological conductivity (K) for every layer
    IF (PRESENT(K)) THEN

      IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == ClappHornberger_) THEN
        K(:,:) =  k_sat(:,:)                                                   &
               & * ws_rel_vol(:,:)**(2._wp * bclapp(:,:) + 3._wp)
      ELSE IF (soilhydmodel == VanGenuchten_) THEN
        WHERE(mvgn(:,:) > 0._wp)
          K(:,:) = k_sat(:,:) * ws_rel_vol(:,:)**(0.5_wp) &
                 & *(1._wp-(1._wp-ws_rel_vol(:,:)**(1._wp/mvgn(:,:)))**mvgn(:,:))**2._wp
        ELSEWHERE
          K(:,:) = 0._wp
        ENDWHERE
      END IF

      ! Reduce transport velocity in presence of ice
      K(:,:) = K(:,:) * ice_imp(:,:)

      ! Make sure percolation from top layer is limited to water content !!!
      K(:,1) = MIN(K(:,1), MAX(0._wp,(ws(:,1)-pwp(:,1)) / dt))

    END IF

    ! Compute hydrological diffusivity (D) for every layer
    IF (PRESENT(D)) THEN

      IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == ClappHornberger_) THEN
        WHERE(ws_range_vol(:,:) > 0._wp)
          D(:,:) = -k_sat(:,:) * mpot_sat(:,:) * bclapp(:,:) / ws_range_vol(:,:) &
                 & * ws_rel_vol(:,:)**(1._wp * bclapp(:,:) + 2._wp)
        ELSEWHERE
          D(:,:) = 0._wp
        END WHERE
      ELSE IF (soilhydmodel == VanGenuchten_) THEN
        IF (.NOT. PRESENT(K)) THEN
          WRITE (message_text,*) 'VanGenuchten diffusivity computation requires hydrological conductivity'
          CALL finish ('get_soilhyd_cond_properties', message_text)
        END IF
        WHERE(ws_range_vol(:,:) > 0._wp)
          D(:,:) = -K(:,:)* mpot_sat(:,:) * bclapp(:,:) / ws_range_vol(:,:) &
                 & * ws_rel_vol(:,:)**(1._wp/mvgn(:,:)) * (1._wp-ws_rel_vol(:,:)**(1._wp/mvgn(:,:)))**(-mvgn(:,:))
        ELSEWHERE
          D(:,:) = 0._wp
        END WHERE
      END IF

      ! Reduce transport velocity in presence of ice
      ! However, ice seems to attract water --> is this the right approach for diffusivity?
      D(:,:) = D(:,:) * ice_imp(:,:)

    END IF

    ! Determining K & D on layer interface
    ! Possibility
    ! 1) Upstream mean [K_lev; Dmax(x_lev, x_lev-1)}
    ! 2) Adjusted arithmetic mean A (of k, D)
    ! Note that in the vertical transport scheme the notation refers to the upper interface of a level.
    ! Hence the loop goes from nsoil_act to 2

    IF (PRESENT(K_inter)) THEN

      IF (.NOT. PRESENT(K)) THEN
        WRITE (message_text,*) 'K at layer interface cannot be computed without computing K on layers'
        CALL finish ('get_soilhyd_cond_properties', message_text)
      END IF

      DO i = 1, nc
        DO j = nsoil, 2, -1
          IF (j <= nsoil_act(i)) THEN
            IF (interpol_mean == Upstream_) THEN
              K_inter(i,j)  =   K(i,j-1)
            ELSE IF (interpol_mean == Arithmetic_) THEN
              IF (K(i,j-1) > K(i,j)) THEN ! MOVEMENT OF WETTING FRONTS LIMITED BY SATURATION OF SUBJACENT LAYER. ELSE ...
                K_inter(i,j)  =   (K(i,j)*dsoil(i,j)+K(i,j-1)*dsoil(i,j-1))/(dsoil(i,j)+dsoil(i,j-1))
              ELSE ! ... PREVENT FLUX FROM COMPLETLY DRYING OUT UPPER LAYER
                K_inter(i,j)  = K(i,j-1)
              END IF
            END IF
          END IF
        END DO
      END DO

    END IF

    IF (PRESENT(D_inter)) THEN

      IF (.NOT. PRESENT(D)) THEN
        WRITE (message_text,*) 'D at layer interface cannot be computed without computing D on layers'
        CALL finish ('get_soilhyd_cond_properties', message_text)
      END IF
      DO i = 1, nc
        DO j = nsoil, 2, -1
          IF (j <= nsoil_act(i)) THEN
            IF (interpol_mean == Upstream_) THEN
              D_inter(i,j)  = MAX(D(i,j),D(i,j-1))
            ELSE IF (interpol_mean == Arithmetic_) THEN
              D_inter(i,j)  = (D(i,j)*dsoil(i,j)+D(i,j-1)*dsoil(i,j-1))/(dsoil(i,j)+dsoil(i,j-1))
            END IF
          END IF
        END DO
      END DO

    END IF

  END SUBROUTINE get_soilhyd_cond_properties

#endif
END MODULE mo_hydro_process
