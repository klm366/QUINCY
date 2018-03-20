!>
!! @brief Contains the routines for the hydro processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!
MODULE mo_hydro_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message, message_text, em_warn

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_surface_hydrology_land, calc_surface_hydrology_glacier, calc_soil_hydrology, &
    & calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare, &
    & get_canopy_conductance, get_water_stress_factor

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  REAL(wp), PARAMETER ::  zwdmin = 0.05_wp
  REAL(wp), PARAMETER ::  zwdcri = 0.9_wp
  REAL(wp), PARAMETER ::  zdrmin = 0.001_wp / (3600._wp*1000._wp)  ! factor for minimum drainage (~ 10 hours)
  REAL(wp), PARAMETER ::  zdrmax = 0.1_wp   / (3600._wp*1000._wp)  ! factor for maximum drainage (~ 4 days)
  REAL(wp), PARAMETER ::  zdrexp = 1.5_wp
  REAL(wp), PARAMETER ::  zfcmin = 1.e-10_wp              ! minimum field capacity (or rather epsilon?)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'

CONTAINS
  !
  ! ================================================================================================================================
  !>
  !! @brief   Compute surface hydrology on glacier-free land
  !!
  !! @par Revision History
  !!   Routine based on the former JSBACH3 routine update_surf_down
  !!   Rewritten in Mai 2016 by V. Gayler
  !!   Adapted for JSBACH4 in February 2017 by R. Schnur
  !!
  !! @param [in]     lstart                T: Start of experiment
  !! @param [in]     delta_time            Time step
  !! @param [in]     l_dynsnow             T: Compute density of snow dynamically
  !! @param [in]     t_srf_unfilt          surface tempemperature [K] (unfiltered)
  !! @param [in]     wind_10m              wind speed at 10m height [m/s]
  !! @param [in]     t_air                 lowest layer atmosphere temperature [K]
  !! @param [in]     skinres_canopy_max    capacity of the canopy skin reservoir [m]
  !! @param [in]     skinres_max           
  !! @param [in]     sfract_srf            surface snow fraction []
  !! @param [in]     wfract_srf            skin reservoir fraction []
  !! @param [in]     hcap_grnd             heat capacity of the uppermost soil layer [J m-2 K-1]
  !! @param [in]     evapotrans            evapotranspiration (including sublimation) [kg m-2 s-1]
  !! @param [in]     evapopot              potential evaporation [kg m-2 s-1]
  !! @param [in]     trans                 transpiration [kg m-2 s-1]
  !! @param [in]     rain                  liquid precipitation [kg m-2 s-1]
  !! @param [in]     snow                  solid precipitation [kg m-2 s-1]
  !! @param [in,out] skinres               water content of the skin reservoir (vegetation and bare soil) [m]
  !! @param [in,out] snow_soil             snow depth at the ground [m water equivalent]
  !! @param [in,out] canopy_snow           snow depth on the canopy [m water equivalent]
  !! @param [in,out] snow_soil_dens        snow density on soil at non-glacier points [m water equivalent]
  !! @param [out]    q_snocpymlt           Heating due to snow melt on canopy [W m-2]
  !! @param [out]    snow_accum            snow budget at non-glacier points [m water equivalent]
  !! @param [out]    snow_melt_soil        snow/ice melt at land points (excluding canopy) [kg m-2 s-1]
  !! @param [out]    evapotrans_soil       evapotranspiration from soil w/o skin and snow reservoirs [kg m-2 s-1]
  !! @param [out]    evapo_skin             evaporation from skin reservoir [kg m-2 s-1]
  !! @param [out]    evapo_snow             evaporation from snow [kg m-2 s-1]
  !! @param [out]    water_to_soil         Water available for infiltration into the soil [m water equivalent]
  !! @param [out]    evapo_deficit          Evaporation deficit flux due to inconsistent treatment of snow evap. [m]
  !
  PURE ELEMENTAL SUBROUTINE calc_surface_hydrology_land (          &
    & lstart, delta_time, l_dynsnow,                               &
    & t_srf_unfilt, wind_10m, t_air,                               &
    & skinres_canopy_max, skinres_max, sfract_srf, wfract_srf,     &
    & hcap_grnd, evapotrans, evapopot, trans,                      &
    & rain, snow, skinres, snow_soil, canopy_snow, snow_soil_dens, &
    & q_snocpymlt, snow_accum, snow_melt_soil,                     &
    & evapotrans_soil, evapo_skin, evapo_snow, water_to_soil, evapo_deficit)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf, dens_snow_min, dens_snow_max, dens_snow
    USE mo_hydro_constants,        ONLY: InterceptionEfficiency
    USE mo_sse_constants,          ONLY: snow_depth_min

    LOGICAL,  INTENT(in) :: lstart, l_dynsnow
    REAL(wp), INTENT(in) :: &
      & delta_time, t_srf_unfilt, wind_10m, t_air, skinres_canopy_max, skinres_max, &
      & sfract_srf, wfract_srf, hcap_grnd, evapotrans, evapopot, trans, rain, snow
    REAL(wp), INTENT(inout) :: &
      & skinres, snow_soil, canopy_snow, snow_soil_dens
    REAL(wp), INTENT(out)   :: &
      & q_snocpymlt, snow_accum, snow_melt_soil, &
      & evapotrans_soil, evapo_skin, evapo_snow, water_to_soil, evapo_deficit
    !
    !  local variables
    !
    REAL(wp) :: &
      & rain_in_m,            & !< amount of rainfall within time step [m] 
      & new_snow,             & !< amount of snowfall within time step [m]
      & snow_soil_old,        & !< amount of snow prior to update [m]
      & evapotrans_in_m,      & !< amount of evapotranspiration within time step [m]
      & evapotrans_soil_in_m, & !< amount of evapotranspiration from soil (w/o snow and skin) within time step [m]
      & evapo_soil_in_m,       & !< amount of ground evaporation within time step [m]
      & evapo_skin_in_m,       & !< amount of evaporation from skin reservoir within time step [m]
      & trans_in_m,           & !< amount of transpiration within time step [m]
      & evapo_snow_in_m,       & !< amount of sublimation from snow within time step [m]
      & snow_melt_can,        & !< snow melt from canopy [m]
      & snow_melt_pot_in_m,   & !< potential amount of snow melt [m]
      & new_snow_canopy,      & !< amount of snowfall on canopy within time step [m]
      & new_snow_soil,        & !< amount of snowfall to soil within time step [m]
      & exp_t, exp_w,         & !< exponents needed for unloading of snow due to temperature / wind
      & snow_blown,           & !< amount of snow blown from canopy to the ground
      & canopy_pre_snow,      & !< snow depth on the canopy prior to its update
      & evapotrans_no_snow,   & !< evapotranspiration without snow evaporation [m]
      & evapo_snow_pot_in_m,   & !< potential snow evaporation
      & evapo_skin_pot_in_m,   & !< potential evaporation from wet surface (skin reservoir)
      & rain_to_skinres,      & ! amount of rain going into the skin reservoir
      & skinres_pre_rain        ! amount of water in the skin reservoir prior to the update [m]


    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp !< [K]
    REAL(wp), PARAMETER :: zc2=1.87E5_wp       !< [Ks]
    REAL(wp), PARAMETER :: zc3=1.56E5_wp       !< [m]

    snow_accum = 0._wp
    snow_soil_old = snow_soil
    snow_melt_soil = 0._wp
    evapotrans_soil = 0._wp
    evapo_skin = 0._wp
    evapo_snow = 0._wp
    water_to_soil = 0._wp
    evapo_deficit = 0._wp

    !-----------------------------------------------------------------
    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m          = rain       * delta_time/rhoh2o
    new_snow           = snow       * delta_time/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m    = evapotrans * delta_time/rhoh2o
    trans_in_m         = trans      * delta_time/rhoh2o
    evapo_snow_pot_in_m = sfract_srf * evapopot * delta_time/rhoh2o
    evapotrans_no_snow = evapotrans_in_m - evapo_snow_pot_in_m
    evapo_skin_pot_in_m = (1._wp - sfract_srf) * wfract_srf * evapopot * delta_time/rhoh2o

    !---------------------------------------------
    !  Budgets of snow (canopy, ground)
    !---------------------------------------------

    ! amount of snow remaining on the canopy
    new_snow_canopy = MIN(new_snow * InterceptionEfficiency, MAX(skinres_canopy_max - canopy_snow, 0._wp))

    ! remaining snow falls on the ground 
    new_snow_soil = new_snow - new_snow_canopy

    ! update snow on the canopy
    ! Note: evaporation happens from canopy, as long as there is enough snow
    canopy_pre_snow = canopy_snow
    canopy_snow = MIN(MAX(0._wp, canopy_pre_snow + new_snow_canopy + evapo_snow_pot_in_m), skinres_canopy_max)
    evapo_snow_in_m = evapo_snow_pot_in_m - (canopy_snow - canopy_pre_snow - new_snow_canopy)
    ! unloading of snow from the canopy due to melting
    exp_t = MAX(0._wp, t_air-zc1)/zc2 * delta_time
    snow_melt_can = canopy_snow * (1._wp-EXP(-exp_t))
    canopy_snow = canopy_snow - snow_melt_can

    ! unloading of snow from the canopy due to wind
    exp_w = wind_10m/zc3 * delta_time
    snow_blown = canopy_snow * (1._wp-EXP(-exp_w))
    canopy_snow   = canopy_snow   - snow_blown
    new_snow_soil = new_snow_soil + snow_blown

    ! Heating due to snow melt on canopy
    ! @Todo Why is heating only computed for canopy? 
    q_snocpymlt = snow_melt_can * rhoh2o * alf / delta_time

    !  Snowfall and sublimation
    !-----------------------------------------------
    snow_soil = snow_soil + new_snow_soil + evapo_snow_in_m

    ! Correction if there was too much snow evaporation
    IF (snow_soil < 0._wp) THEN
      evapotrans_no_snow = evapotrans_no_snow + snow_soil
      evapo_deficit = snow_soil
      snow_soil = 0._wp
    ELSE
      evapo_deficit = 0._wp
    END IF

    !  Snow melt
    !------------------------

    IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
      IF (t_srf_unfilt > tmelt) THEN
        snow_melt_pot_in_m = hcap_grnd * (t_srf_unfilt-tmelt)/(alf*rhoh2o) !  potential snow and ice melt [m]
        IF (snow_soil > 0._wp) THEN
          snow_melt_soil = MAX(MIN(snow_melt_pot_in_m, snow_soil), 0._wp)  ! snow melt limited by actual snow depth [m]
          snow_soil = snow_soil - snow_melt_soil                           ! reduce snow depth according to melting
        ! surface temperature correction for snow/ice melt is applied in task "snowmelt_correction"
        END IF
      END IF
    END IF

    !  Snow budget and meltwater
    !-----------------------------------------------------
    skinres = skinres + snow_melt_can                                    ! Add melt water from canopy to skin reservoir
    water_to_soil = snow_melt_soil + MAX(0._wp, skinres-skinres_max)     ! Excess water enters the soil, and
    skinres = MIN(skinres_max, skinres)                                  ! skin reservoir is limited to the maximum
    snow_accum = new_snow + evapo_snow_pot_in_m - snow_melt_soil - snow_melt_can

    !-----------------------------------------------------------
    !  Budget of water in skin reservoir (canopy and ground)
    !-----------------------------------------------------------
    ! Interception of rain
    rain_to_skinres = MIN(rain_in_m * InterceptionEfficiency, MAX(skinres_max - skinres, 0._wp)) ! Part of rain goes to skin reservoir
    ! Note: the MAX in following line is only to account for apparent precision error (small negative values)
    water_to_soil = water_to_soil + MAX(rain_in_m - rain_to_skinres, 0._wp)        ! Rest of rain goes to soil
    skinres_pre_rain = skinres
    skinres = MIN(MAX(0._wp, skinres_pre_rain + rain_to_skinres + evapo_skin_pot_in_m), skinres_max)

    ! note: at this point ignoring the amount of dew that exceeds the skinreservoir
    !       as dew is calculated later based on evapo_soil (f(evapotrans_soil))
    evapo_skin_in_m = skinres - (skinres_pre_rain + rain_to_skinres)
    evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m
    IF (evapo_skin_pot_in_m < 0._wp) THEN
      evapo_deficit = evapo_deficit + evapo_skin_pot_in_m - evapo_skin_in_m
    END IF
    evapo_soil_in_m = evapotrans_soil_in_m - trans_in_m  ! Evaporation from ground

    ! positive values of evaporation and transpiration (dew) are added to water_to_soil
    ! negative fluxes change soil moisture later in calc_soil_hydrology
    IF (evapo_soil_in_m > 0._wp) water_to_soil = water_to_soil + evapo_soil_in_m
    IF (trans_in_m > 0._wp)     water_to_soil = water_to_soil + trans_in_m

    ! Compute snow density
    IF (l_dynsnow) THEN
      IF (snow_soil_old * rhoh2o / snow_soil_dens > snow_depth_min) THEN
        ! Update density of old snow following eq. 31 of Verseghy, 1991
        snow_soil_dens = dens_snow_max - (dens_snow_max - snow_soil_dens) * EXP(-0.01_wp * delta_time/3600._wp)
        ! Calculate weighted mean density from old and fresh snow
        IF (snow_soil > snow_soil_old) &
          & snow_soil_dens = (snow_soil_dens * snow_soil_old + dens_snow_min * (snow_soil - snow_soil_old)) / snow_soil
      ELSE
        snow_soil_dens = dens_snow_min  ! fresh snow
      END IF
    ELSE
      snow_soil_dens = dens_snow
    END IF

    ! Transform fluxes in output from [m] to [kg m-2 s-1]
    snow_melt_soil = snow_melt_soil * rhoh2o / delta_time

    evapotrans_soil = evapotrans_soil_in_m * rhoh2o / delta_time
    evapo_skin = evapo_skin_in_m  * rhoh2o / delta_time
    evapo_snow = sfract_srf * evapopot

    RETURN

  END SUBROUTINE calc_surface_hydrology_land

  ! ================================================================================================================================
  !>
  !! @brief   calc surface hydrology on glaciers
  !!
  !! @par Revision History
  !!   Routine based on the former JSBACH3 routine update_surf_down
  !!   Rewritten in Mai 2016 by V. Gayler
  !!   Adapted for JSBACH4 in February 2017 by R. Schnur
  !!
  !! @param [in]     lstart                T: Start of experiment
  !! @param [in]     delta_time            Time step
  !! @param [in]     t_srf_unfilt          surface tempemperature [K] (unfiltered)
  !! @param [in]     sfract_srf            surface snow fraction []
  !! @param [in]     hcap_grnd             heat capacity of the uppermost soil layer [J m-2 K-1]
  !! @param [in]     evapotrans            evapotranspiration (including sublimation) [kg m-2 s-1]
  !! @param [in]     evapopot              potential evaporation [kg m-2 s-1]
  !! @param [in]     rain                  liquid precipitation [kg m-2 s-1]
  !! @param [in]     snow                  solid precipitation [kg m-2 s-1]
  !! @param [in,out] glacier_depth         glacier depth (snow and ice) [m water equivalent]
  !! @param [out]    q_snocpymlt           Heating due to snow melt on canopy [W m-2]
  !! @param [out]    glacier_melt          Snow/ice melt at glacier points [kg m-2 s-1]
  !! @param [out]    glacier_runoff        glacier runoff (rain+snow/ice melt, no calving) [kg m-2 s-1]
  !! @param [out]    pme_glacier           precipitation minus sublimation on glacier [kg m-2 s-1]
  !
  PURE ELEMENTAL SUBROUTINE calc_surface_hydrology_glacier (   &
    & lstart, delta_time,                                         &
    & t_srf_unfilt,                                               &
    & sfract_srf, hcap_grnd,     &
    & evapotrans, evapopot, rain, snow,                           &
    & glacier_depth, q_snocpymlt,                                 &
    & glacier_melt, glacier_runoff,                               &
    & pme_glacier)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf
    USE mo_hydro_constants,        ONLY: InterceptionEfficiency

    LOGICAL,  INTENT(in) :: lstart
    REAL(wp), INTENT(in) :: &
      & delta_time, t_srf_unfilt, &
      & sfract_srf, hcap_grnd, evapotrans, evapopot, rain, snow
    REAL(wp), INTENT(inout) :: &
      & glacier_depth
    REAL(wp), INTENT(out)   :: &
      & q_snocpymlt, glacier_melt, glacier_runoff, pme_glacier
    !
    !  local variables
    !
    REAL(wp) :: &
      & rain_in_m,           & !< amount of rainfall within time step [m] 
      & new_snow,            & !< amount of snowfall within time step [m]
      & evapotrans_in_m,     & !< amount of evapotranspiration within time step [m]
      & evapo_snow_in_m,      & !< amount of sublimation from snow within time step [m]
      & snow_melt_pot_in_m,  & !< potential amount of snow melt [m]
      & new_snow_soil,       & !< amount of snowfall to soil within time step [m]
      & exp_t, exp_w,        & !< exponents needed for unloading of snow due to temperature / wind
      & evapo_snow_pot_in_m,  & !< potential snow evaporation
      & pme_glacier_in_m,    & !< P-E in m
      & glacier_runoff_in_m, & !< Runoff in m
      & glacier_melt_in_m      !< Snow/ice melt in m

    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp !< [K]
    REAL(wp), PARAMETER :: zc2=1.87E5_wp       !< [Ks]
    REAL(wp), PARAMETER :: zc3=1.56E5_wp       !< [m]

    !----------------------------------------------------------------------------------------------

    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m          = rain       * delta_time/rhoh2o
    new_snow           = snow       * delta_time/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m    = evapotrans * delta_time/rhoh2o
    evapo_snow_pot_in_m = sfract_srf * evapopot * delta_time/rhoh2o

    !  Snowfall and sublimation on glaciers
    !---------------------------------------

    glacier_depth       = glacier_depth + new_snow + evapo_snow_pot_in_m   ! glacier depth [m]
    pme_glacier_in_m    = rain_in_m     + new_snow + evapotrans_in_m      ! P-E on glaciers [m]
    glacier_runoff_in_m = rain_in_m                                       ! there is no infiltration on glaciers


    !  Snow and glacier melt
    !------------------------

    glacier_melt_in_m = 0._wp
    IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
      IF (t_srf_unfilt > tmelt) THEN
        snow_melt_pot_in_m  = hcap_grnd * (t_srf_unfilt-tmelt)/(alf*rhoh2o) !  potential snow and ice melt [m]
        glacier_melt_in_m   = snow_melt_pot_in_m                       ! there is an unlimited amount of snow
        glacier_depth       = glacier_depth - glacier_melt_in_m        ! reduce glacier depth according to melting
        glacier_runoff_in_m = glacier_runoff_in_m + glacier_melt_in_m  ! add melt water to the runoff
        ! surface temperature correction for snow/ice melt is applied in task "snowmelt_correction"
      END IF
    END IF

    q_snocpymlt = 0._wp

    glacier_melt   = glacier_melt_in_m   * rhoh2o / delta_time
    glacier_runoff = glacier_runoff_in_m * rhoh2o / delta_time
    pme_glacier    = pme_glacier_in_m    * rhoh2o / delta_time

    RETURN

  END SUBROUTINE calc_surface_hydrology_glacier

  !
  ! @brief      Calculates runoff, drainage and soil moisture
  !
  ! @param [in]      nidx               Vector lenth
  ! @param [in]      nsoil              Number of soil layers
  ! @param [in]      dz                 (Fixed) Thicknesses of soil layers [m]
  ! @param [in]      nlat               (Effective) number of latitudes
  ! @param [in]      delta_time         Time step [s]
  ! @param [in]      soil_depth_l       Thicknesses of soil layers until bedrocj [m]
  ! @param [in]      root_depth_l       Thicknesses of soil layers until rooting depth [m]
  ! @param [in]      fract_org_l        Fractions per soil layer of organic comonent []
  ! @param [in]      hyd_cond_sat_l     Hydraulic conductivity of saturated soil [m/s]
  ! @param [in]      matrix_pot_l       Soil matrix potential [m]
  ! @param [in]      bclapp_l           Exponent B in Clapp and Hornberger
  ! @param [in]      pore_size_index_l  Soil pore size distribution index used in Van Genuchten method
  ! @param [in]      vol_porosity_l     Soil porosity [m/m]
  ! @param [in]      vol_field_cap_l    Volumetric field capacity [m/m]
  ! @param [in]      vol_p_wilt_l       Volumetric permanent wilting point [m/m]
  ! @param [in]      sfract_srf         Snow cover fraction []
  ! @param [in]      wfract_srf         Wet surface fraction (skin reservoir) []
  ! @param [in]      ws_max_mineral     Maximum water holding capacity of non-organic soil [m]
  ! @param [in]      oro_stddev         Sub-grid standard deviation of orography [m]
  ! @param [in]      t_soil             Temperature of the upper soil layer [K]
  ! @param [in]      water_to_soil      Amount of water going into the soil [m]
  ! @param [in]      evapotrans_soil    Evapotranspiration from the soil (w/o snow or canopy) [kg m-2 s-1]
  ! @param [in]      trans              Evapotranspiration [kg m-2 s-1]
  ! @param [in,out]  ws                 Total water content of soil column in rooting zone [m]
  ! @param [in,out]  ws_l               Content un-frozen water of soil layers in rooting zone [m]
  ! @param [in,out]  soil_ice_l         Content of frozen water os soil layers in rooting zone [m]
  ! @param [out]     wsat_l             Maximum water storage in each soil layers [m]
  ! @param [out]     field_cap_l        Water content at field capacity in soil layers [m]
  ! @param [out]     p_wilt_l           Water content at permanent wilting point in soil layers [m]
  ! @param [out]     runoff             Runoff 
  ! @param [out]     drainage           Drainage
  !
  SUBROUTINE calc_soil_hydrology( &
    ! in
    & nidx, nsoil, dz, nlat, delta_time, ltpe_closed, ltpe_open,                                                   &
    & soil_depth_l, root_depth_l, fract_org_l,                                                                     &
    & hyd_cond_sat_l, matrix_pot_l, bclapp_l, pore_size_index_l, vol_porosity_l, vol_field_cap_l, vol_p_wilt_l,    &
    & sfract_srf, wfract_srf, ws_max_mineral, oro_stddev, t_soil, water_to_soil, evapotrans_soil, trans,           &
    ! inout
    & ws, ws_l, soil_ice_l, w_soil_overflow,                                                                       &
    ! out
    & wsat_l, field_cap_l, p_wilt_l, runoff, drainage                                                              &
    & )

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o
    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max
    USE mo_hydro_util,             ONLY: get_water_in_root_zone

    INTEGER,  INTENT(in) :: &
      & nidx, nsoil, nlat
    REAL(wp), INTENT(in) :: &
      & dz(nsoil)
    REAL(wp), INTENT(in) :: &
      & delta_time
    LOGICAL, INTENT(in) :: ltpe_closed
    LOGICAL, INTENT(in) :: ltpe_open
    REAL(wp), INTENT(in), DIMENSION(nidx, nsoil) :: &
      & soil_depth_l, root_depth_l, fract_org_l, &
      & hyd_cond_sat_l, matrix_pot_l, bclapp_l, pore_size_index_l, vol_porosity_l, vol_field_cap_l, vol_p_wilt_l
    REAL(wp), INTENT(in), DIMENSION(nidx) :: &
      & sfract_srf, wfract_srf, ws_max_mineral, oro_stddev, t_soil, water_to_soil, &
      & evapotrans_soil, trans
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: &
      & ws, w_soil_overflow
    REAL(wp), INTENT(inout), DIMENSION(nidx, nsoil) :: &
      & ws_l, soil_ice_l
    REAL(wp), INTENT(out), DIMENSION(nidx) :: &
      & runoff, drainage
    REAL(wp), INTENT(out), DIMENSION(nidx, nsoil) :: &
      & wsat_l, field_cap_l, p_wilt_l

    ! Local variables
    REAL(wp) ::    &   
      !< variables for infiltration calculation
      & sigma_0,              & !< minimum value
      & sigma_max,            & !< resolution-dependent maximum
      & steepness,            & !< sub-grid scale steepness of the orography
      & zbws, zb1, zbm, zconw1, zvol, &
      ! Other variables
      & evapo_soil_in_m(nidx), & !< evaporation from the soil (without transpiration) [m]
      & trans_in_m(nidx),     & !< transpiration within time step [m]
      & ws_rel,               & !< relative soil moisture (0: wilting point, 1: field capacity)
      & infilt(nidx)            !< infiltration [m]

    REAL(wp) :: &
      & ws_max(nidx),        & ! water holding capacity of the soil [m]
      & ws_crit(nidx),       & ! soil moisture above which drainage strongly increases
      & ws_min_drain(nidx),  & ! minimum soil moisture for drainage
      & ws_incl_ice(nidx),   & ! soil moisture including ice [m]
      & ice_rootzone(nidx)     ! amount of ice in the rootzone [m]

    INTEGER :: jl

    INTEGER, PARAMETER :: ilog = 0 ! Switch for debugging output
    INTEGER :: jllog = 0

    sigma_0   = oro_var_min
    sigma_max = oro_var_max * 64._wp / REAL(nlat, wp) ! @todo: site level

    trans_in_m(:) = trans(:)  * delta_time / rhoh2o

    wsat_l     (:,:) = MAX(vol_porosity_l (:,:) * soil_depth_l(:,:) - soil_ice_l(:,:), 0._wp)
    field_cap_l(:,:) = MAX(vol_field_cap_l(:,:) * soil_depth_l(:,:) - soil_ice_l(:,:), 0._wp)
    WHERE (vol_porosity_l(:,:) > 0._wp)
      p_wilt_l(:,:) = wsat_l(:,:) * vol_p_wilt_l(:,:) / vol_porosity_l(:,:)
    ELSE WHERE
      p_wilt_l(:,:) = 0._wp
    END WHERE

    ws_max(:) = get_water_in_root_zone(field_cap_l(:,:), soil_depth_l(:,:), root_depth_l(:,:))
    WHERE (SUM(fract_org_l(:,:), DIM=2) <= EPSILON(1._wp))
      ws_max(:) = ws_max_mineral(:)
    END WHERE

    ice_rootzone(:) = get_water_in_root_zone(soil_ice_l(:,:), soil_depth_l(:,:), root_depth_l(:,:))
    ws          (:) = get_water_in_root_zone(ws_l      (:,:), soil_depth_l(:,:), root_depth_l(:,:))

    !--------------------------------------------------------------------------------------------
    !  water budgets
    !--------------------------------------------------------------------------------------------

    !  soil water reservoir
    ! ----------------------

    ! note: at this point ignoring the amount of dew that exceeds the skinreservoir
    !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
    evapo_soil_in_m(:) = evapotrans_soil(:) * delta_time / rhoh2o - trans_in_m(:)

    ! infiltration
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

    DO jl = 1, nidx

      ws_min_drain(jl) = zwdmin * ws_max(jl)            ! minimum amount of soil water for drainage

      steepness = MAX(0._wp, oro_stddev(jl)-sigma_0) / (oro_stddev(jl)+sigma_max)

      zbws   = MAX(MIN(steepness,0.5_wp), 0.01_wp)   ! wie passt das zur Doku???
      zb1    = 1._wp+zbws
      zbm    = 1._wp/zb1
      zconw1 = ws_max(jl)*zb1

      !  surface runoff and infiltration
      ! -----------------------------
      IF (ltpe_closed) THEN
         runoff(jl) = 0._wp
         infilt(jl) = water_to_soil(jl)
      ELSE
        IF (t_soil(jl) < tmelt) THEN
          runoff(jl) = water_to_soil(jl)       !  no infiltration as uppermost soil layer is frozen -> runoff
          infilt(jl) = 0._wp
        ELSE
          ws_incl_ice(jl) = ws(jl) + ice_rootzone(jl)                    ! soil water including ice
          IF (water_to_soil(jl) > 0._wp .AND. ws_incl_ice(jl) > ws_min_drain(jl)) THEN  ! soil water content above limit for drainage

            ws_rel = MIN(1._wp, ws_incl_ice(jl) / ws_max(jl))            !  relative soil moisture
            zvol   = (1._wp - ws_rel)**zbm - water_to_soil(jl)/zconw1          ! ???

            runoff(jl) = water_to_soil(jl) - (ws_max(jl)-ws_incl_ice(jl))    ! if >0: amount of water exceeding water holding capacity
            IF (zvol > 0._wp) runoff(jl) = runoff(jl) + ws_max(jl)*zvol**zb1
            runoff(jl) = MAX(MIN(runoff(jl), water_to_soil(jl)), 0._wp)      ! no negative runoff, runoff limited by available water
            infilt(jl) = water_to_soil(jl) - runoff(jl)                      ! water not going into runoff is infiltrated
          ELSE
            runoff(jl) = 0._wp                                           ! no runoff as soil is too dry
            infilt(jl) = water_to_soil(jl)                                   ! all water infiltrated
          END IF
        END IF
      END IF
    ENDDO

    !---------------------
    ! drainage and runoff
    ! ----------------------------------------------------------------------------------

    ! initilalization of drainage 
    drainage(:) = 0._wp

    ! Changes by evaporation fluxes
    CALL digest_evapotrans(nidx, nsoil, root_depth_l, soil_depth_l, &
         field_cap_l, trans_in_m, evapo_soil_in_m, ws_l, soil_ice_l)

    ! For numerical reasons, half of the water is infiltrated 
    ! before soilhyd and the other half of it afterwards
    infilt(:) = infilt(:) / 2._wp

    ! update soil moisture according to infiltration
    !    - first half of the water
    CALL digest_infilt(nidx, nsoil, ltpe_closed, field_cap_l, infilt, &
         ws_l, drainage, w_soil_overflow)

    ! Drainage calculation with multi-layer soil scheme
    CALL soilhyd(nidx, nsoil, delta_time, ilog, jllog, ws_l, soil_depth_l, &
         wsat_l, field_cap_l, hyd_cond_sat_l, vol_porosity_l, bclapp_l,   &
         matrix_pot_l, drainage, pore_size_index_l, p_wilt_l, dz,       &
         soil_ice_l, w_soil_overflow, ltpe_closed, ltpe_open)

    ! update soil moisture according to infiltration
    !    - second half of the water
    CALL digest_infilt(nidx, nsoil, ltpe_closed, field_cap_l, infilt, &
         ws_l, drainage, w_soil_overflow)

    ! recalculate soil moisture from the soil moisture of the individual levels
    ! Note: ws only includes the soil moisture within the root zone
    ws(:) = get_water_in_root_zone(ws_l(:,:), soil_depth_l(:,:), root_depth_l(:,:))

    runoff(:)   = runoff  (:) * rhoh2o / delta_time
    drainage(:) = drainage(:) * rhoh2o / delta_time

  END SUBROUTINE calc_soil_hydrology

  !------------------------------------------------------------------------------------------------
  SUBROUTINE soilhyd(nidx, nsoil, delta_time, ilog, jllog, ws_l, dsoil,      &
       wsat_l, field_cap_l, hyd_cond_sat_l, vol_porosity_l, bclapp_l, matrix_pot_l,  &
       drainage, pore_size_index_l, p_wilt_l, dz, soil_ice, w_soil_overflow, ltpe_closed, ltpe_open)

    !----------------------------------------------------------------------------------------------
    !
    !  SOIL HYDROLOGY - Routine to be build into REMO to calculate soil hydrology for 5 soil layers
    !
    !  Calculation of percolation and vertical diffusion of water
    !
    !  Programming and development: Stefan Hagemann (MPI-M)
    !----------------------------------------------------------------------------------------------

    INTEGER,  INTENT(IN)            :: nidx                    ! vector length
    INTEGER,  INTENT(IN)            :: nsoil                   ! number of soil layers
    REAL(wp), INTENT(IN)            :: delta_time              ! model time step legth [s]
    INTEGER,  INTENT(IN)            :: ilog                    !
    INTEGER,  INTENT(IN)            :: jllog                   !
    REAL(wp), INTENT(IN)            :: dsoil(nidx, nsoil)      ! soil depth until bedrock within each layer [m]
    REAL(wp), INTENT(IN)            :: wsat_l(nidx,nsoil)      ! Maximum storage in each soil layer (porosity) [m]
    REAL(wp), INTENT(IN)            :: field_cap_l(nidx,nsoil) ! field capacity of the soil layer [m]
    REAL(wp), INTENT(IN)            :: hyd_cond_sat_l(nidx,nsoil) ! hydraulic conductivity of saturated soil [m/s]
    REAL(wp), INTENT(IN)            :: vol_porosity_l(nidx,nsoil) ! volumetric soil porosity [m/m]
    REAL(wp), INTENT(IN)            :: bclapp_l(nidx,nsoil)       ! exponent B in Clapp and Hornberger
    REAL(wp), INTENT(IN)            :: matrix_pot_l(nidx,nsoil)   ! matrix potential [m]
    REAL(wp), INTENT(IN)            :: pore_size_index_l(nidx,nsoil) ! soil pore size distribution index used in Van Genuchten
    REAL(wp), INTENT(IN)            :: p_wilt_l(nidx,nsoil )   ! wilting point of the soil layer
    REAL(wp), INTENT(IN)            :: dz(nsoil)               ! vertical depth of each soil layer
    REAL(wp), INTENT(INOUT)         :: ws_l(nidx, nsoil)       ! soil moisture of each layer [m]
    REAL(wp), INTENT(INOUT)         :: drainage(nidx)          ! amount of drainage within timestep [m]
    REAL(wp), INTENT(INOUT)         :: w_soil_overflow(nidx)
    REAL(wp), INTENT(IN), OPTIONAL  :: soil_ice(nidx, nsoil)   ! amount of ice in the soil [m]

    LOGICAL, INTENT(IN) :: ltpe_closed
    LOGICAL, INTENT(IN) :: ltpe_open

    INTEGER,  PARAMETER :: iredu  = 0         ! reduction of conductivity with depth yes/no (1/0)
    REAL(wp), PARAMETER :: credu  = 1.0_wp

    INTEGER  :: i                             ! looping index
    REAL(wp) :: percolation(nidx,0:nsoil)     ! downward flux from one soil layer to the layer below [m/s]
                                              !     and later amount of perculation per time step [m]
    REAL(wp) :: diffus_l(nidx,nsoil)          ! diffusivity at mid-depth of the layer
    REAL(wp) :: ws_min_drain(nidx)            ! minimum soil moisture for drainage
    REAL(wp) :: ws_rel(nidx)                  ! normalized (~relative) soil mositure
    REAL(wp) :: ws_crit(nidx)                 ! soil moisture above which drainage strongly increases [m]
    REAL(wp) :: drain_crit(nidx)              ! drainage above a soil moisture of ws_crit [m/s]
    REAL(wp) :: ws_tmp(nidx)                  ! temporary array for soil moisture
    REAL(wp) :: drain_l(nidx)                 ! drainage from the current layer
    LOGICAL  :: soil_above(nidx)              ! is there a soil layer above? 
    LOGICAL  :: soil_below(nidx)              ! does the layer below have a soil fraction?
    REAL(wp) :: zda(nidx,nsoil)               ! diffusion coefficients
    REAL(wp) :: zdb(nidx,nsoil)               !     ''
    REAL(wp) :: zdc(nidx,nsoil)               !     ''
    REAL(wp) :: ztri(nidx,nsoil)              !     ''
    REAL(wp) :: diffusivity(nidx,nsoil-1)     ! diffusivity between two layers, i.e. at the bottom of each layer
    REAL(wp) :: zkredu(nsoil)                 ! reduction factor of hydraulic conductivity with depth for percolation
    REAL(wp) :: zdiflog(10), wslog(nsoil)     ! needed for debugging output
    REAL(wp) :: remoist(nidx)                 ! for TPE

    REAL(wp) :: zmvg_l(nidx,nsoil)            ! exponent M in Van Genuchten scheme = pore_size_index / (pore_size_index + 1)

    INTEGER :: jl

    ! define reduction factor of hydraulic conductivity with depth (for percolation)
    IF (iredu == 1) THEN
      zkredu(1) = (1._wp + EXP(-credu*dz(1)) )/2._wp
      DO i = 2, nsoil
        zkredu(i) = (EXP(-credu*SUM(dz(1:i-1))) + EXP(-credu*SUM(dz(1:i)))) / 2._wp
      END DO
    ELSE
      ! conductivity is constant with depth
      zkredu(1:nsoil) = 1._wp
    END IF

    zmvg_l(:,:)           = pore_size_index_l(:,:) / (pore_size_index_l(:,:)+1._wp)

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
    !  theta:   current soil moisture        -> ws_l
    !  theta_s: soil moisture at stauration  -> field_cap
    !  theta_r: residual soil moisture       -> p_wilt

    percolation(:,nsoil) = 0._wp
    DO i = 1, nsoil-1
      WHERE (field_cap_l(:,i) > zfcmin .AND. field_cap_l(:,i) > p_wilt_l(:,i))
        ws_rel(:) = MAX(ws_l(:,i) - p_wilt_l(:,i), 0._wp) / (field_cap_l(:,i) - p_wilt_l(:,i))
        ws_rel(:) = MIN(ws_rel(:), 1._wp)
        percolation(:,i) =  hyd_cond_sat_l(:,i)*zkredu(i) * SQRT(ws_rel(:)) &
                           * (1._wp - (1._wp - ws_rel(:)**(1._wp/zmvg_l(:,i)))**zmvg_l(:,i))**2._wp &
                           * delta_time     ! [m/s] -> [m]
      ELSEWHERE
        percolation(:,i) = 0._wp
      END WHERE
    END DO


    ! drainage of the lowest layer above bedrock replaces the percolation calculated above

    DO i = 1, nsoil

      ! is there a layer above or a layer with soil below?
      soil_above(:) = (i > 1)
      IF (i == nsoil) THEN
        soil_below(:) = .false.
      ELSE
        soil_below(:) = dsoil(:,i+1) > 0._wp
      END IF

      ws_crit(:) = field_cap_l(:,i) * zwdcri     ! soil moisture above which drainage strongly increases

      ! calculate drainage of the lowest layer above bedrock
      WHERE (dsoil(:,i) > 0._wp .AND. .NOT. soil_below(:))          ! lowest layer above bedrock

        WHERE (ws_l(:,i) <= p_wilt_l(:,i))
          drain_l(:) = 0._wp              ! no drainage below the wilting point
        ELSEWHERE
          WHERE (field_cap_l(:,i) > p_wilt_l(:,i))
            ! drain_min = minimum factor * normalized soil moisture
            drain_l(:) = zdrmin * (ws_l(:,i) - p_wilt_l(:,i)) / (field_cap_l(:,i) - p_wilt_l(:,i))
          ELSEWHERE
            drain_l(:) = 0._wp
          END WHERE
          WHERE (ws_l(:,i) > ws_crit(:))
            drain_crit(:) = (zdrmax-zdrmin) * ((ws_l(:,i)-ws_crit(:)) / (wsat_l(:,i) - ws_crit(:)))**zdrexp
            drain_l(:) = drain_l(:) + MIN(drain_crit(:), (ws_l(:,i)-ws_crit(:))/delta_time)
          END WHERE
          ! drainage is limited by the available water of the layer
          drain_l(:) = MIN(drain_l(:), (ws_l(:,i)-p_wilt_l(:,i))/delta_time)
        END WHERE
        percolation(:,i) = drain_l(:) * delta_time  ! [m/s] -> [m]

      ELSEWHERE (dsoil(:,i) > 0._wp)

        ! Drainage is calculated with percolation removed from WS_L
        !   de Vrese --- WS_L corrected for incomming percolation

        ! Remove percolating water from the corresponding soil layer and add percolation from the above layer
        WHERE (.NOT. soil_above(:))       ! surface layer
          ws_tmp(:) = ws_l(:,i) - percolation(:,i)
        ELSEWHERE
          ws_tmp(:) = ws_l(:,i) - percolation(:,i) + percolation(:,i-1)
        END WHERE

        WHERE (ws_tmp(:) > p_wilt_l(:,i))          ! drainage occurs only above wilting point
          WHERE (field_cap_l(:,i) > p_wilt_l(:,i))
          ! drain_min = minimum factor * normalized soil moisture
            drain_l(:) = zdrmin * (ws_tmp(:) - p_wilt_l(:,i)) / (field_cap_l(:,i)-p_wilt_l(:,i))
          ELSEWHERE
            drain_l(:) = 0._wp
          END WHERE
          WHERE (ws_tmp(:) > ws_crit(:))
            drain_crit(:) = (zdrmax-zdrmin) * ((ws_tmp(:) - ws_crit(:)) / (wsat_l(:,i) - ws_crit(:)))**zdrexp
            drain_l(:) = drain_l(:) + MIN(drain_crit(:), (ws_tmp(:) - ws_crit(:))/delta_time)
          END WHERE
          ! drainage is limited by the soil moisture within the layer
          drain_l(:) = MIN(drain_l(:) * delta_time, (ws_tmp(:) - p_wilt_l(:,i)))
        END WHERE
      END WHERE

      IF (.NOT. ltpe_closed) THEN
        WHERE (dsoil(:,i) > 0._wp .AND. soil_below(:) .AND. ws_tmp(:) > p_wilt_l(:,i))
          ws_l(:,i) = ws_l(:,i) - drain_l(:)       ! remove drainage of the layer from soil moisture
          drainage(:) = drainage(:) + drain_l(:)   ! and add it to the total drainage
        END WHERE
      END IF                                       ! no drainage for TPE closed case

    END DO  ! soil level

    !---------------------------------------------------------------------------------------------- 
    !  water diffusion
    ! -----------------

    ! numerical recipes/Richtmyer & -Morton diffusion

    ! calculating the diffusivity at mid depth of layer i
    DO i = 1, nsoil

      ws_rel(:) = 0._wp
        WHERE (ws_l(:,i) > 0._wp)
          ws_rel(:) =  ws_l(:,i) / dsoil(:,i)  &
                      / (vol_porosity_l(:,i) - (soil_ice(:,i) / dsoil(:,i)))
        END WHERE

      diffus_l(:,i) = 0._wp
      WHERE (ws_rel(:) > zfcmin)
        diffus_l(:,i) =  bclapp_l(:,i) * hyd_cond_sat_l(:,i) * matrix_pot_l(:,i)  &
                        / (ws_l(:,i) / dsoil(:,i)) *  ws_rel(:) ** (bclapp_l(:,i)+3._wp)
      END WHERE
    END DO

    ! calculating the diffusivity at the bottom of layer i
    DO i = 1, nsoil-1
      WHERE (diffus_l(:,i+1) == 0._wp .OR. diffus_l(:,i) == 0._wp)
        diffusivity(:,i) = 0._wp
      ELSEWHERE
        diffusivity(:,i) = diffus_l(:,i+1) * diffus_l(:,i) / &
                        ((diffus_l(:,i+1)*dsoil(:,i) + diffus_l(:,i)*dsoil(:,i+1))/(dsoil(:,i+1)+dsoil(:,i)))
      END WHERE
    END DO

    ! calculation of diffusion coefficients

    ! deepest layer
    zda(:,nsoil) = 0._wp
    WHERE (dsoil(:,nsoil) > 0._wp)
      zdb(:,nsoil) = diffusivity(:,nsoil-1) * delta_time / dsoil(:,nsoil) / (dsoil(:,nsoil)+dsoil(:,nsoil-1)) * 2._wp
    ELSEWHERE
      zdb(:,nsoil) = 0._wp
    END WHERE

    ! middle layers
    DO i = nsoil-1, 2, -1
      WHERE (dsoil(:,i) > 0._wp)
        WHERE (dsoil(:,i+1) > 0._wp)
          zda(:,i) = diffusivity(:,i) * delta_time / dsoil(:,i) / (dsoil(:,i)+dsoil(:,i+1)) * 2._wp
        ELSEWHERE
          zda(:,i) = 0._wp
        END WHERE
        zdb(:,i) = diffusivity(:,i-1) * delta_time / dsoil(:,i) / (dsoil(:,i)+dsoil(:,i-1)) * 2._wp
      ELSEWHERE
        zda(:,i) = 0._wp
        zdb(:,i) = 0._wp
      END WHERE
    END DO

    ! uppermost layer
    zdb(:, 1) = 0._wp
    WHERE (dsoil(:,1) > 0)
      WHERE (dsoil(:,2) > 0)
        zda(:,1) = diffusivity(:,1) * delta_time / dsoil(:,1) / (dsoil(:,1)+dsoil(:,2)) * 2._wp
      ELSEWHERE
        zda(:,1) = 0._wp
      END WHERE
    ELSEWHERE
      zda(:,1) = 0._wp
    END WHERE

    ! temporary storage for diagnostic printing
    IF (ilog == 1 .OR. ilog == 3) wslog(1:nsoil) = ws_l(jllog, 1:nsoil)


    ! ---------------------
    !  update soil wetness  - with respect to diffusion
    ! ---------------------
    ! routine TRIDIAG from Numerical Recipes, p. 43
    !   -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, WS_L(T)=RI, WS_L(T+1)=UI

    WHERE (dsoil(:,1) > 0._wp)
      zdc(:,1) = zda(:,1) + zdb(:,1) + 1._wp
      ws_l(:,1) = ws_l(:,1) / zdc(:,1)
    END WHERE
 
    ! decomposition and forward substitution
    DO i = 2, nsoil
      WHERE (dsoil(:,i) > 0._wp)
        ztri(:,i) = -zda(:,i-1) / zdc(:,i-1)
        zdc(:,i)  = zda(:,i) + zdb(:,i) + 1._wp + zdb(:,i)*ztri(:,i)
        ws_l(:,i) = (ws_l(:,i) / dsoil(:,i) + zdb(:,i) * ws_l(:,i-1)/dsoil(:,i-1)) / zdc(:,i) * dsoil(:,i)
      END WHERE
    END DO

    ! backsubstitution
    DO i = nsoil-1, 1, -1
      WHERE (dsoil(:,i+1) > 0._wp)
        ws_l(:,i) = ws_l(:,i) - ztri(:,i+1) * ws_l(:,i+1) * (dsoil(:,i)/dsoil(:,i+1))
      END WHERE
    END DO

    ! temporary storage for diagnostic printing
    IF (ilog == 1 .OR. ilog == 3) THEN
      zdiflog(1:nsoil) = wslog(1:nsoil) - ws_l(jllog,1:nsoil)
    END IF


    ! treatment of layer overflow due to diffusion
    IF (.NOT. ltpe_closed) THEN
      ! limit soil moisture to field capacity and add overflow to drainage
      DO i=1,nsoil
        WHERE (dsoil(:,i) > 0._wp .AND. ws_l(:,i) > field_cap_l(:,i))
          drainage(:) = drainage(:) + ws_l(:,i)-field_cap_l(:,i)
          ws_l(:,i) = field_cap_l(:,i)
        END WHERE
      END DO
    ELSE
      ! Modification for TPE (closed)
      DO i=1,nsoil
        WHERE (dsoil(:,i) > 0._wp)
          WHERE (ws_l(:,i) > field_cap_l(:,i))
            ! redirect overflow from drainage to overflow pool
            w_soil_overflow(:) = w_soil_overflow(:) + ws_l(:,i) - field_cap_l(:,i)
            ws_l(:,i) = field_cap_l(:,i)
          ELSEWHERE ! ws_l <= field_cap_l
            ! refill soil moisture from overflow pool
            remoist(:) = MIN(w_soil_overflow(:), field_cap_l(:,i) - ws_l(:,i))
            ws_l(:,i) = ws_l(:,i) + remoist(:)
            w_soil_overflow(:) = w_soil_overflow(:) - remoist(:)
          END WHERE
        END WHERE
      END DO
    END IF

    ! Modification for TPE (open)
    IF (ltpe_open) THEN
      DO i = 1, nsoil
        !IF (i == nsoil - 2 .OR. i == nsoil - 1 .OR. i == nsoil) THEN   ! special adaptation of Shirisha
        IF (i == nsoil-1 .OR. i == nsoil) THEN
          WHERE (dsoil(:,i) > 0._wp)
             ! keep the two lowest soil layers always very wet
             ws_l(:,i) = MAX(ws_l(:,i),field_cap_l(:,i) * 0.9_wp)
          ENDWHERE
        END IF
      END DO
    END IF

    !----------------------------------------------------------------------------------------------
    ! debugging output for grid cell jllog
    IF (ilog == 2 .OR. ilog == 3) THEN
      write(*,'(A10,1X, 11(E16.9E2, 1X) )')  &
        "SOILHYD1: ", hyd_cond_sat_l(jllog,1:nsoil)*1000._wp, ws_l(jllog,1:nsoil)*1000._wp,  &
        field_cap_l(jllog,1:nsoil)*1000._wp
    END IF
    !----------------------------------------------------------------------------------------------


    ! ---------------------
    !  update soil wetness  - with respect to percolation
    ! ---------------------

    ! percolation per time step must be smaller than the usable soil layer water content
    DO i = 1, nsoil
      WHERE (dsoil(:,i) > 0._wp .AND. percolation(:,i) > 0._wp)
        ws_min_drain(:) = zwdmin * field_cap_l(:,i)              ! minimum soil moisture for drainage
        WHERE (ws_l(:,i)-percolation(:,i) < ws_min_drain(:))     ! more percolation than water available 
          WHERE (ws_l(:,i) > ws_min_drain(:))                    !   there is water for perculation
            percolation(:,i) = ws_l(:,i) - ws_min_drain(:)       !   all available water perculates
            ws_l(:,i) = ws_min_drain(:)                          !   soil moisture is limited to minimum
          ELSEWHERE
            percolation(:,i) = 0._wp
          END WHERE
        ELSEWHERE                                                ! there is enough soil water
          ws_l(:,i) = ws_l(:,i) - percolation(:,i)
        END WHERE
      END WHERE
    END DO

    ! Percolation is added to downward layer or, if lowest layer or bedrock is reached, to
    ! drainage (not TPE closed) resp. overflow pool (TPE closed)
    DO i = 1, nsoil-1
      WHERE (dsoil(:,i+1) > 0._wp)
        ws_l(:,i+1) = ws_l(:,i+1) + percolation(:,i)
      END WHERE
      IF (.NOT. ltpe_closed) THEN
        WHERE (dsoil(:,i+1) <= 0._wp)       ! bedrock below
          drainage(:) = drainage(:) + percolation(:,i)
        END WHERE
      ELSE  ! TPE closed case
        WHERE (dsoil(:,i+1) <= 0._wp)       ! bedrock below
          w_soil_overflow(:) = w_soil_overflow(:) + percolation(:,i)
        END WHERE
      END IF
      percolation(:,i) = 0._wp
    END DO

    ! Percolation overflow is not allowed from layers below the surface layer.
    ! Instead it is assumed that water piles upwards from saturated layers.
    DO i = nsoil, 2, -1
      WHERE (dsoil(:,i) > 0._wp .AND. ws_l(:,i) > field_cap_l(:,i))
        ! soil moisture above field capacity:
        !   reduce percolation from the above layer by the amount that does not fit into the below layer
        percolation(:,i-1) = MAX(percolation(:,i-1) - ws_l(:,i) + field_cap_l(:,i), 0._wp)
        ws_l(:,i-1) = ws_l(:,i-1) + ws_l(:,i) - field_cap_l(:,i)
        ws_l(:,i) = field_cap_l(:,i)
      END WHERE
    END DO

    IF (.NOT. ltpe_closed) THEN
      drainage(:) = drainage(:) + percolation(:,nsoil)
    ELSE  ! TPE closed case
      ! redirect percolation from drainage to overflow pool
      w_soil_overflow(:) = w_soil_overflow(:) +  percolation(:,nsoil)
    END IF

    ! Overflow from the surface layer enters drainage, resp. the overflow pool for TPE closed. 
    ! This should not happen, this is just for security.
    IF (.NOT. ltpe_closed) THEN
      WHERE (ws_l(:,1) > field_cap_l(:,1))
        drainage(:) = drainage(:) + ws_l(:,1) - field_cap_l(:,1)
        ws_l(:,1) = field_cap_l(:,1)
      END WHERE
    ELSE      ! TPE closed case
      WHERE (ws_l(:,1) > field_cap_l(:,1))
        ! redirect overflow to overflow pool
        w_soil_overflow(:) = w_soil_overflow(:) + ws_l(:,1) - field_cap_l(:,1)
        ws_l(:,1) = field_cap_l(:,1)
      END WHERE
    END IF

    !----------------------------------------------------------------------------------------------
    ! debugging output for grid cell jllog
    IF (ilog == 1 .OR. ilog == 3) THEN
      DO i=1, nsoil
        zdiflog(nsoil+i) = percolation(jllog,i)
      END DO
      write(*,'(A10,1X,11(E15.6E3,1X))') "SOILFLUX",zdiflog(1:10)*1000._wp, drainage(jllog)*1000._wp
      write(*,'(A10,1X,10(E16.9E2,1X))') "SOILHYD2: ", ws_l(jllog,1:nsoil)*1000._wp, dsoil(jllog,1:nsoil)*1000._wp
    END IF

  END SUBROUTINE soilhyd

  SUBROUTINE digest_evapotrans(nidx, nsoil, droot, dsoil, &
       field_cap_l, trans_in_m, evapo_soil, ws_l, soil_ice)
    !----------------------------------------------------------------------------------------------
    !  The routine is based on routine soilchange by Stefan Hagemann. It replaces the part for
    !  evapotranspiration changes (isch=1)
    !
    !  Update soil moisture due to evaporation and transpiration fluxes. The fluxes had been
    !  calculated with the old bucket scheme, and are not necessarily consistent to soil moisture
    !  of the multi-layer scheme.
    !  Therefore evaporation might need to be reduced (reduce_evapo). 
    !  @todo: However this term is only
    !  exported for diagnostics and the water deficit is distributed between the soil layers.
    !----------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nidx                    ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    REAL(wp), INTENT(in)     :: droot(nidx, nsoil)      ! root depth within layer [m]
    REAL(wp), INTENT(in)     :: dsoil(nidx, nsoil)      ! soil depth until bedrock within layer [m]
    REAL(wp), INTENT(in)     :: field_cap_l(nidx,nsoil) ! field capacity of the layer
    REAL(wp), INTENT(in)     :: trans_in_m(nidx)        ! transpiration [m]
    REAL(wp), INTENT(in)     :: evapo_soil(nidx)         ! bare soil evaporation [m]

    REAL(wp), INTENT(inout)  :: ws_l(nidx,nsoil)        ! water content of the soil layer [m]
    REAL(wp), INTENT(inout)  :: soil_ice(nidx,nsoil)    ! ice content of the soil layer [m]

    ! local variables

    INTEGER  :: jk
    REAL(wp) :: deficit_evapo(nidx)   ! water deficit due to evaporation
    REAL(wp) :: deficit_trans(nidx)  ! water deficit due to tranpiration
    REAL(wp) :: deficit_l(nidx)      ! water deficit within a soil layer
    REAL(wp) :: remaining(nidx)      ! water remaining in the soil layer
    REAL(wp) :: rootfract(nidx)      ! fraction of the soil layer within the root zone
    REAL(wp) :: fixed(nidx)          ! amount of water below the wilting point, not available for plant
    REAL(wp) :: trans_layer(nidx)    ! transpiration from the layer
    REAL(wp) :: ws(nidx)             ! water within all soil layers
    REAL(wp) :: root_depth(nidx)     ! Total depth of rootzone
    REAL(wp) :: reduce_evapo(nidx)

    INTEGER :: max2d(2)

    ! parameters
    ! @todo - take from mo_hydro_constants.f90
    REAL(wp), PARAMETER :: zwilt  = 0.35_wp 

    !----------------------------------------------------------------------------------------------
    !  update soil moisture due to evaporation and transpiration fluxes
    !----------------------------------------------------------------------------------------------

    !  bare soil evaporation
    ! -----------------------

    deficit_evapo(:) = 0._wp   ! water missing for evaporation due to limited soil water
    WHERE (evapo_soil(:) < 0._wp)            ! pos. evaporation flux (not dew)
      ws_l(:,1) = ws_l(:,1) + evapo_soil(:)  ! take water from uppermost soil layer 
      WHERE (ws_l(:,1) < 0._wp)             ! if there is not enough soil water
        deficit_evapo(:) = ws_l(:,1)         !    this is the amount of water missing
        ws_l(:,1) = 0._wp                   !    soil moisture cannot be negative
      END WHERE
    END WHERE

    !  transpiration
    ! ---------------

    root_depth(:) = SUM(droot(:,:), DIM=2)

    ! sealed grid cells (rooting_depth zero) and positive transpiration
    !    no water reachable for transpiration, add flux to deficit term
    ! @todo: what does rooting depth=zero for bare soil tile?
    deficit_trans(:) = MERGE(trans_in_m(:), 0._wp, root_depth(:) <= 0._wp .AND. trans_in_m(:) < 0._wp)

    DO jk =1, nsoil

      deficit_l(:) = 0._wp

      ! grid cell not sealed and positive transpiration (not dew)
      !    take water equally from all over the root zone
      WHERE (root_depth(:) > 0._wp .AND. trans_in_m(:) < 0._wp)
          
        ! water needed for transpiration from the current layer 
        trans_layer(:) = trans_in_m(:) * droot(:,jk)/root_depth(:)

        WHERE (droot(:,jk) > 0._wp)

          rootfract(:) = droot(:,jk) / dsoil(:,jk)             ! fraction of the soil layer that is within the root zone
          fixed(:) = field_cap_l(:,jk) * zwilt * rootfract(:)  ! amount of water in the root zone, not available for plants

          WHERE (ws_l(:,jk) * rootfract(:) >= fixed(:))        ! there is water in the root zone available for plants
            remaining(:) = ws_l(:,jk) * rootfract(:) + trans_layer(:)  ! water that would remain after transpiration
            WHERE (remaining(:) < fixed(:))                                ! transpiration exceeds available water
              deficit_l(:) = remaining(:) - fixed(:)
              ws_l(:,jk) = fixed(:)     + ws_l(:,jk) * (1._wp-rootfract(:))  
            ELSEWHERE                                                        ! enough soil water
              deficit_l(:) = 0._wp
              ws_l(:,jk) = remaining(:) + ws_l(:,jk) * (1._wp-rootfract(:))
            END WHERE
          ELSEWHERE                                            ! too dry: no water available for transpiration
            deficit_l(:) = trans_layer(:)
          END WHERE

          ! sum up deficits from the different layers
          deficit_trans(:) = deficit_trans(:) + deficit_l(:)

        END WHERE
      END WHERE

    END DO

    ! -----------------------------------------
    !  reduction needed for evapotranspiration
    ! -----------------------------------------
      reduce_evapo(:) = deficit_evapo(:) + deficit_trans(:)

    ! --------------------------
    !  soil moisture correction
    ! --------------------------
    ! reduce_evapo is only used for diagnostics. If the therm is not zero, soil moisture must be
    ! changed to close the land surface water balance.

    ws(:) = SUM(ws_l(:,:), DIM=2) + SUM(soil_ice(:,:), DIM=2)  ! water and ice  within the column

    DO jk=1, nsoil

      WHERE (reduce_evapo(:) < 0._wp .AND. ws(:) > 0._wp)
        ! reduce soil moisture and ice relative to the water content of each layer
        ws_l(:,jk)     = MAX(0._wp, ws_l(:,jk)     + reduce_evapo(:) * ws_l(:,jk)/ws(:))
        soil_ice(:,jk) = MAX(0._wp, soil_ice(:,jk) + reduce_evapo(:) * soil_ice(:,jk)/ws(:))
      END WHERE
    END DO

    ! make sure there is no negative soil ice
    IF (ANY(soil_ice(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil ice: ', MINVAL(soil_ice), ' at ', MINLOC(soil_ice)
      CALL finish ('digest_evapotrans', message_text)
      ! CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      ! soil_ice(:,:) = MAX(soil_ice(:,:), 0._wp)
    END IF

    ! make sure there is no negative soil moisture and soil moisture does not exceed field capacity
    IF (ANY(ws_l(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture: ', MINVAL(ws_l), ' at ', MINLOC(ws_l)
      CALL finish ('digest_evapotrans', message_text)
      ! CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      ! ws_l(:,:) = MAX(ws_l(:,:), 0._wp)
    END IF
    IF (ANY(ws_l(:,:) > field_cap_l(:,:) + zfcmin)) THEN
      max2d(:) = MAXLOC(ws_l(:,:) - field_cap_l(:,:))
      WRITE (message_text,*) 'liquid water content above field capacity: at ',max2d(:) ,': ws_l: ', ws_l(max2d(1),max2d(2)), &
           ' field_cap_l: ', field_cap_l(max2d(1),max2d(2))
      CALL finish ('digest_evapotrans', message_text)
      ! CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      ! ws_l(:,:) = MIN(ws_l(:,:), field_cap_l(:,:))
    END IF 

  END SUBROUTINE digest_evapotrans
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE digest_infilt (nidx, nsoil, ltpe_closed, field_cap_l, infilt, ws_l, drainage, w_soil_overflow)
  !------------------------------------------------------------------------------------------------
  !  The routine is based on routine soilchange by Stefan Hagemann. It replaces the part for 
  !  infiltration changes (isch=2)   
  !------------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nidx                    ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    LOGICAL,  INTENT(in)     :: ltpe_closed             ! TPE closed case
    REAL(wp), INTENT(in)     :: field_cap_l(nidx,nsoil) ! field capacity of the layer
    REAL(wp), INTENT(in)     :: infilt(nidx)            ! infiltration in [m/s]

    REAL(wp), INTENT(inout)  :: ws_l(nidx,nsoil)        ! water content of the soil layer [m]
    REAL(wp), INTENT(inout)  :: drainage(nidx)          ! drainage [m]
    REAL(wp), INTENT(inout)  :: w_soil_overflow(nidx)   ! For TPE

    ! local variables

    INTEGER  :: jk
    REAL(wp) :: excess(nidx)         ! excess water after infiltration

    !------------------------------------------------------------------------------------------------
    !  update soil moisture according to infiltration
    !------------------------------------------------------------------------------------------------

    excess(:) = infilt(:)                   ! amount of water for infiltration at the surface [m]
    DO  jk = 1, nsoil
      ws_l(:,jk) = ws_l(:,jk) + excess(:)
      IF (.NOT. ltpe_closed) THEN
        WHERE (ws_l(:,jk) > field_cap_l(:,jk))             ! soil moisture exceeds field capacity
          excess(:) = ws_l(:,jk) - field_cap_l(:,jk)       !   excess water -> goes into level below
          ws_l(:,jk) = field_cap_l(:,jk)                   !   limit soil moisture to field capacity
        ELSEWHERE
          excess(:) = 0._wp
        END WHERE
      ELSE  ! TPE closed case
        WHERE (ws_l(:,jk) > field_cap_l(:,jk))
          w_soil_overflow(:) = w_soil_overflow(:) + ws_l(:,jk) - field_cap_l(:,jk)
          ws_l(:,jk) = field_cap_l(:,jk)
        END WHERE
        excess(:) = 0._wp
      END IF
    END DO
  
    ! if all soil layers are at field capacity and there is still excess water, it is added to drainage
    WHERE (excess(:) > 0._wp)
      drainage(:) = drainage(:) + excess(:)
    END WHERE

    ! make sure there is no negative soil moisture
    IF (ANY(ws_l(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture: ', MINVAL(ws_l), ' at ', MINLOC(ws_l)
      CALL finish ('digest_infilt', message_text)
      ! CALL message ('digest_infilt', message_text, all_print=.TRUE., level=em_warn)
      ! ws_l(:,:) = MAX(ws_l(:,:), 0._wp)
    END IF

  END SUBROUTINE digest_infilt

  !--------------  !!
  ! @brief      Calculates the snow fraction on lake ice.
  !
  ! @param [in]     w_snow_lice      Snow depth on lake ice [m water equivalent]
  ! @param [out]    fract_snow_lice  Fraction of snow on lake ice
  !!
  ELEMENTAL SUBROUTINE calc_wskin_fractions_lice( &
    & w_snow_lice,                                & ! in
    & fract_snow_lice                             & ! out
    & )

    REAL(wp),  INTENT(in)    :: w_snow_lice
    REAL(wp),  INTENT(out)   :: fract_snow_lice

    fract_snow_lice = TANH(w_snow_lice * 100._wp)

  END SUBROUTINE calc_wskin_fractions_lice
      
  ELEMENTAL SUBROUTINE calc_wskin_fractions_veg( &
    & dtime,          & ! in
    & steplen,        & ! in
    & wsr_max,        & ! in
    & oro_stddev,     & ! in
    & t_srf_old,      & ! in
    & press_srf,      & ! in
    & drag_srf,       & ! in
    & q_air,          & ! in
    & w_skin,         & ! in
    & w_snow_soil,    & ! in
    & w_snow_can,     & ! in
    & fract_snow_can, & ! inout
    & fract_water    ,    & ! out
    & fract_snow_soil     & ! out
    & )

    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in) :: &
      & dtime,              &
      & steplen,            &
      & oro_stddev,         &
      & t_srf_old,          &
      & press_srf,          &
      & drag_srf,           &
      & q_air,              &
      & wsr_max,            &
      & w_skin,             &
      & w_snow_soil,        &
      & w_snow_can

    REAL(wp), INTENT(out) :: &
      & fract_snow_can,      &
      & fract_water    ,     &
      & fract_snow_soil

    REAL(wp) ::       &
      & qsat_srf_old, &
      & evapo_pot
    
    ! Snow fraction on canopy
    fract_snow_can = MERGE(MIN(1._wp, w_snow_can / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    ! Wet skin reservoir fraction on soil and canopy (between 0 and 1)
    fract_water = MERGE(MIN(1._wp, w_skin / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    ! Snow fraction on soil
    fract_snow_soil = Get_snow_fract_noforest(w_snow_soil, oro_stddev) ! snow on soil below forest
    fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, &
                               fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))

    
    ! Potential evaporation using old values of air and surface humidity
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evapo_pot =  heat_transfer_coef(drag_srf, steplen) * (q_air - qsat_srf_old) * -1._wp ! Positive upwards
    
    ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy; same for skin reservoir
    IF (fract_snow_soil > 0._wp) THEN
      ! @todo Shouldn't one take rhoice here insteady of rhoh2o?
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evapo_pot * dtime &
        &                                            / (rhoh2o * (w_snow_soil + w_snow_can)) &
        &                                    )
    END IF
    IF (fract_water > 0._wp ) THEN
      fract_water = fract_water    / MAX(1._wp, (1._wp - fract_snow_soil) * evapo_pot * dtime &
        &                                       / (rhoh2o * MAX(EPSILON(1._wp), w_skin)) &
        &                                    )
    END IF

  END SUBROUTINE calc_wskin_fractions_veg

  ELEMENTAL SUBROUTINE calc_wskin_fractions_bare( &
    & dtime,          & ! in
    & steplen,        & ! in
    & wsr_max,        & ! in
    & oro_stddev,     & ! in
    & t_srf,          & ! in
    & press_srf,      & ! in
    & drag_srf,       & ! in
    & q_air,          & ! in
    & w_skin,         & ! in
    & w_snow_soil,    & ! in
    & fract_water,    & ! out
    & fract_snow_soil    & ! out
    & )

    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp),  INTENT(in)    ::                &
                                dtime,         &
                                steplen,       &
                                oro_stddev,    &
                                t_srf,         &
                                press_srf,     &
                                drag_srf,      &
                                q_air,         &
                                wsr_max,       &
                                w_skin,        &
                                w_snow_soil

    REAL(wp),  INTENT(out)   ::                  &
                               fract_water,      &
                               fract_snow_soil

    REAL(wp) ::              &
               t_srf_old,    &
               qsat_srf_old, &
               evapo_pot

    fract_water = MERGE(MIN(1._wp, w_skin / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    ! Snow fraction on soil
    fract_snow_soil = Get_snow_fract_noforest(w_snow_soil, oro_stddev) ! snow on soil below forest

    ! Potential evaporation using old values of air and surface humidity
    t_srf_old = t_srf
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evapo_pot =  heat_transfer_coef(drag_srf, steplen) * (qsat_srf_old - q_air) ! Positive upwards

    ! Reduce snow cover on soil if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy together
    ! R: In the following lines I just kicked out the contribution of w_skin_can as compared to JSBACH3. For the intension of
    !    the parametrization this can be done so easily for the bare type and it is consistent with JSBACH3.
    IF (fract_snow_soil > 0._wp .AND. w_snow_soil > EPSILON(1._wp)) THEN
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evapo_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp),(w_snow_soil))))
    END IF
    ! Same for water cover on soil
    IF (fract_water > 0._wp .AND. w_skin > EPSILON(1._wp)) THEN
      fract_water = fract_water / MAX(1._wp, (1._wp - fract_snow_soil) * evapo_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp), w_skin)))
    END IF

  END SUBROUTINE calc_wskin_fractions_bare

  !< Roesch et al. 2002, Climate Dynamics
  !
  REAL(wp) ELEMENTAL FUNCTION Get_snow_fract_noforest(snow, orodev)

    USE mo_hydro_constants, ONLY: wsn2fract_const, wsn2fract_eps, wsn2fract_sigfac

    REAL(wp), INTENT(in) :: &
      & snow,   &
      & orodev

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_canopy_cond_unstressed_simple'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_canopy_cond_stressed_simple'

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

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
