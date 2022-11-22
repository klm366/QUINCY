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
    & get_canopy_conductance, get_water_stress_factor,                                        &
    & vert_transport, get_conductivies_diffusivities_column, get_conductivies_diffusivities

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  REAL(wp), PARAMETER ::  zwdmin = 0.05_wp
  REAL(wp), PARAMETER ::  zwdcri = 0.9_wp
  REAL(wp), PARAMETER ::  zwdfre = 1.0_wp
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
    & skinres_canopy_max, skinres_max, weq_pond_max,                &
    & fract_snow, fract_wet,                                        &
    & hcap_grnd, evapotrans, evapopot, transpiration, rain, snow,   &
    & wtr_skin, weq_snow_soil, weq_snow_can, snow_soil_dens,        &
    & fract_pond, weq_pond, wtr_pond, ice_pond,                     &
    & q_snocpymlt, snow_accum, snowmelt_soil, pond_freeze, pond_melt, &
    & evapotrans_soil, evapo_skin, evapo_snow, evapo_pond, water_to_soil, evapo_deficit)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf, dens_snow_min, dens_snow_max, dens_snow
    USE mo_hydro_constants,        ONLY: InterceptionEfficiency
    USE mo_sse_constants,          ONLY: snow_depth_min

    LOGICAL,  INTENT(in) :: lstart, l_dynsnow
    REAL(wp), INTENT(in) :: dtime, snow_depth_max
    !$claw model-data
    REAL(wp), INTENT(in) ::                                                         &
      & t_unfilt, wind_10m, t_air, skinres_canopy_max, skinres_max, weq_pond_max,   &
      & fract_snow, fract_wet, hcap_grnd, evapotrans, evapopot, transpiration, rain, snow
    REAL(wp), INTENT(inout) :: &
      & wtr_skin, weq_snow_soil, weq_snow_can, snow_soil_dens, &
      & fract_pond, weq_pond, wtr_pond, ice_pond
    REAL(wp), INTENT(out)   ::                   &
      & q_snocpymlt, snow_accum, snowmelt_soil, pond_freeze, pond_melt, &
      & evapotrans_soil, evapo_skin, evapo_snow, evapo_pond, water_to_soil, evapo_deficit
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
      & evapo_pond_in_m,      & !< amount of evaporation from pond reservoir within time step [m]
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
      & evapo_pond_pot_in_m,  & !< potential evaporation from pond reservoir
      & rain_to_skinres,      & !< amount of rain going into the skin reservoir
      & wtr_skin_pre_rain       !< amount of water in the skin reservoir prior to the update [m]

! pdv >
    REAL(wp) ::               &
      & t_pond,               &
      & pond_change_pot
!pdv

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

    pond_freeze = 0._wp
    pond_melt   = 0._wp

    !-----------------------------------------------------------------
    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m           = rain          * dtime/rhoh2o
    new_snow            = snow          * dtime/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m     = evapotrans    * dtime/rhoh2o
    transpiration_in_m  = transpiration * dtime/rhoh2o
    evapo_snow_pot_in_m = fract_snow    * evapopot * dtime/rhoh2o
    evapotrans_no_snow  = evapotrans_in_m - evapo_snow_pot_in_m
    evapo_skin_pot_in_m = (1._wp - fract_snow) * fract_wet  * evapopot * dtime/rhoh2o
    evapo_pond_pot_in_m = (1._wp - fract_snow) * fract_pond * evapopot * dtime/rhoh2o

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

    !pdv >
    ! Freezing and melting of pond reservoir
    ! 1.) Potential melting:
    t_pond = t_unfilt - snowmelt_soil / hcap_grnd * (alf*rhoh2o)
    IF (t_pond > tmelt .AND. ice_pond > 0._wp) THEN
      pond_change_pot = hcap_grnd * (t_pond-tmelt)/(alf*rhoh2o)
      pond_melt       = MAX(MIN(pond_change_pot, ice_pond), 0._wp)
      ice_pond        = ice_pond - pond_melt
      wtr_pond        = wtr_pond + pond_melt
    ELSE IF (t_pond < tmelt .AND. wtr_pond > 0._wp) THEN
      pond_change_pot = hcap_grnd * (tmelt-t_pond)/(alf*rhoh2o)
      pond_freeze     = MAX(MIN(pond_change_pot, wtr_pond), 0._wp)
      wtr_pond        = wtr_pond - pond_freeze
      ice_pond        = ice_pond + pond_freeze
    END IF
    !pdv <

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

    evapo_pond_in_m = MAX(-(wtr_pond + ice_pond), evapo_pond_pot_in_m)

    evapo_pond_in_m = MAX(evapo_pond_in_m,          &
                    & MIN(0._wp, evapotrans_no_snow - evapo_skin_in_m - transpiration_in_m))


    IF (-evapo_pond_in_m > ice_pond ) THEN
      wtr_pond = wtr_pond + ice_pond + evapo_pond_in_m
      ice_pond = 0._wp
    ELSE
      ice_pond = ice_pond + evapo_pond_in_m
    END IF
    weq_pond   = ice_pond + wtr_pond
    IF (weq_pond_max > 1.E-20_wp .AND. weq_pond > 1.E-20_wp) THEN
      fract_pond = MIN(1._wp, (weq_pond / weq_pond_max)**0.5_wp)
    ELSE
      fract_pond = 0._wp
    END IF

    evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m - evapo_pond_in_m

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
    pond_freeze   = pond_freeze   * rhoh2o / dtime
    pond_melt     = pond_melt     * rhoh2o / dtime

    evapotrans_soil = evapotrans_soil_in_m * rhoh2o / dtime
    evapo_skin = evapo_skin_in_m  * rhoh2o / dtime
    evapo_pond = evapo_pond_in_m  * rhoh2o / dtime
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
  ! @param [in]      oro_stddev          Sub-grid standard deviation of orography [m]
  ! @param [in]      t_soil_sl1          Temperature of the upper soil layer [K]
  ! @param [in]      water_to_soil       Amount of water going into the soil [m]
  ! @param [in]      evapotrans_soil     Evapotranspiration from the soil (w/o snow or canopy) [kg m-2 s-1]
  ! @param [in]      transpiration       Evapotranspiration [kg m-2 s-1]
  ! @param [in]      wtr_soil_pot_scool_sl Amount of potentially supercooled water [m]
  ! @param [in,out]  wtr_soil_sl         Content un-frozen water of soil layers in rooting zone [m]
  ! @param [in,out]  ice_soil_sl         Content of frozen water of soil layers in rooting zone [m]
  ! @param [out]     wsat_sl             Maximum possible amount of liquid water in soil layer
  !                                        - reduced by amount of ice [m]
  ! @param [out]     field_cap_sl        Liquid water content at field capacity in soil layer [m]
  ! @param [out]     p_wilt_sl           Water content at permanent wilting point in soil layers
  !                                        - reduced relative to ice fraction) [m]
  ! @param [out]     runoff              Runoff
  ! @param [out]     drainage            Drainage
  !
  SUBROUTINE calc_soil_hydrology(                                                           &
    ! in
    & nc, nsoil, nlat, dtime, ltpe_closed, ltpe_open,                                       &
    & l_latdrain,    &
    & manning_n, dt_sat, max_slope,                                         &
    & i_shm, i_interpol, i_infil, i_testing1, i_testing2,                                   &
    & l_pf_soil,                                                                            &
    & soil_depth_sl, root_depth_sl, fract_org_sl,                                           &
    & hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, vol_porosity_sl,       &
    & vol_field_cap_sl, vol_p_wilt_sl, wtr_soil_pot_scool_sl,                               &
    & slope, oro_stddev, t_soil_sl1, wtr_pond_max,                                               &
    & water_to_soil, evapotrans_soil, transpiration,                                        &
    ! inout
    & wtr_soil_sl, ice_soil_sl, wtr_pond, tpe_overflow,                                     &
    !pdv >
    & weq_soilhyd_err, pdv_1d_out, &
    !pdv <
    ! out
    & wsat_sl, field_cap_sl, p_wilt_sl, &
    & infilt, runoff, runoff_horton, runoff_dunne, drainage, drain_rock, drain_lat_sl       &
    & )

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o
    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max
    USE mo_hydro_util,             ONLY: get_amount_in_rootzone

    INTEGER,  INTENT(in) :: &
      & nc, nsoil, nlat
    REAL(wp), INTENT(in) :: &
      & dtime
    LOGICAL, INTENT(in) :: ltpe_closed
    LOGICAL, INTENT(in) :: ltpe_open
    LOGICAL, INTENT(in) :: l_latdrain

    REAL(wp), INTENT(in) :: &
      & manning_n, dt_sat, max_slope

    INTEGER, INTENT(in) :: &
      & i_shm, i_interpol, i_infil, i_testing1, i_testing2 

    LOGICAL, INTENT(in), DIMENSION(:) :: l_pf_soil

    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & soil_depth_sl, root_depth_sl, fract_org_sl, hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, &
      & vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl, wtr_soil_pot_scool_sl
    REAL(wp), INTENT(in), DIMENSION(:) ::               &
      & slope, oro_stddev, t_soil_sl1, wtr_pond_max, water_to_soil, &
      & evapotrans_soil, transpiration
    REAL(wp), INTENT(inout), DIMENSION(:) :: &
      & wtr_pond, tpe_overflow
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & wtr_soil_sl, ice_soil_sl, weq_soilhyd_err
    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & infilt, runoff, runoff_horton,     &
      & runoff_dunne, drainage, drain_rock
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & drain_lat_sl
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & wsat_sl, field_cap_sl, p_wilt_sl

    REAL(wp), INTENT(inout), DIMENSION(:) :: &
      & pdv_1d_out

    ! Local variables
    REAL(wp) ::                       &
      & evapo_soil_in_m(nc),          & !< evaporation from the soil (without transpiration) [m]
      & transpiration_in_m(nc)          !< transpiration within time step [m]      

    REAL(wp) ::              &
      & ws_max(nc),          & ! water holding capacity of the soil [m]
      & wtr_rootzone(nc),    & ! amount of liquid water in the rootzone [m]
      & ice_rootzone(nc),    & ! amount of ice in the rootzone [m]
      & manning(nc)            ! Manning roughness coefficient (0.15 for flowplain with trees) 

    REAL(wp) ::                 &
      & wsat_sl_wi(nc,nsoil),      &
      & field_cap_sl_wi(nc,nsoil), &
      & p_wilt_sl_wi(nc,nsoil)

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
    wsat_sl_wi     (:,:) = vol_porosity_sl (:,:) * soil_depth_sl(:,:)  ! -> diffusion
    field_cap_sl_wi(:,:) = vol_field_cap_sl(:,:) * soil_depth_sl(:,:)  ! -> percolation
    p_wilt_sl_wi   (:,:) = vol_p_wilt_sl(:,:)    * soil_depth_sl(:,:)

    DO i = 1, nc
      !$claw acc loop seq
      DO jl = 1, nsoil
        ! Note: here we are calculating soil saturation (this does not correspond to weq_rootzone_max, which was
        !    calculated from field capacity)
        wsat_sl     (i,jl) = MAX(vol_porosity_sl (i,jl) * soil_depth_sl(i,jl) - ice_soil_sl(i,jl), 0._wp)  ! -> diffusion
        field_cap_sl(i,jl) = MAX(vol_field_cap_sl(i,jl) * soil_depth_sl(i,jl) - ice_soil_sl(i,jl), 0._wp)  ! -> percolation
        ! @todo: to be discussed: Does it make sense, to reduce water content at wilting point relative to the ice fraction?
        IF (vol_porosity_sl(i,jl) > 0._wp) THEN
          p_wilt_sl(i,jl) = wsat_sl(i,jl) * vol_p_wilt_sl(i,jl) / vol_porosity_sl(i,jl)
        ELSE
          p_wilt_sl(i,jl) = 0._wp
        END IF
      END DO
    END DO

    !$claw acc end parallel
    !$claw acc end data

    manning(:) = manning_n

    CALL get_amount_in_rootzone(wsat_sl_wi(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), ws_max(:))
    CALL get_amount_in_rootzone(ice_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))


    !
    !> 1.2 Update soil water with respect to evapotransration fluxes
    !
    ! Note: It is important at which point evapotrans is removed from soil !!!
    !$claw acc data present(evapo_soil_in_m,evapotrans_soil,transpiration_in_m,transpiration)
    !$claw acc parallel
    !$claw acc loop gang
    DO i = 1, nc
      transpiration_in_m(i) = transpiration(i)  * dtime / rhoh2o
      ! Note: at this point ignoring the amount of dew that exceeds the skin reservoir
      !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
      evapo_soil_in_m(i) = evapotrans_soil(i) * dtime / rhoh2o - transpiration_in_m(i)
    END DO


    !$claw acc present(ws_min_drain,ws_max,oro_stddev,runoff,infilt)    &
    !$claw acc end parallel

    runoff_dunne(:)     = 0._wp
    runoff_horton(:)    = 0._wp
    runoff(:)           = 0._wp

    IF (ltpe_closed) THEN
      ! Terra planet setup without runoff
      runoff(:) = 0._wp
      infilt(:) = water_to_soil(:)
    ELSE
      CALL infiltration_new(nc,nlat, dt_sat, dtime, max_slope, wtr_pond_max(:), i_infil, i_shm, &
                          & oro_stddev(:), manning(:), soil_depth_sl(:,1),                      &
                          & pore_size_index_sl(:,1), bclapp_sl(:,1), wsat_sl_wi(:,1),           &
                          & field_cap_sl_wi(:,1), p_wilt_sl_wi(:,1), hyd_cond_sat_sl(:,1),      &
                          & wtr_soil_sl(:,1), ice_soil_sl(:,1), wtr_soil_pot_scool_sl(:,1),     &
                          & water_to_soil(:), wtr_pond(:),                                      &
                          & runoff_horton(:),infilt(:),pdv_1d_out(:))
    END IF  

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

! pdv >
    drain_rock(:)       = 0._wp
    drain_lat_sl(:,:)   = 0._wp

    CALL soilhyd_new(nc, nsoil, dtime, i_shm, i_interpol,                                     &
      &     l_latdrain, l_pf_soil,                                     &
      &     wtr_pond_max, slope, soil_depth_sl, wsat_sl_wi, field_cap_sl_wi, hyd_cond_sat_sl, &
      &     bclapp_sl, matrix_pot_sl, pore_size_index_sl, p_wilt_sl_wi,                       &
      &     wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl,  wtr_pond,                       &
      &     infilt,runoff_dunne, drainage, drain_rock, drain_lat_sl,                          &
      &     tpe_overflow, &
      &     weq_soilhyd_err, &
      &     ltpe_closed, ltpe_open)

! pdv <

    ! recalculate soil moisture from the soil moisture of the individual levels
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:), soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    !$claw acc data present(runoff,drainage)
    !$claw acc parallel loop gang vector
!pdv > 
    runoff(:)            = runoff_horton(:) + runoff_dunne(:)
!pdv <


    DO i = 1, nc
      infilt(i)          = infilt     (i)     * rhoh2o / dtime
      runoff(i)          = runoff     (i)     * rhoh2o / dtime
      runoff_dunne(i)    = runoff_dunne(i)    * rhoh2o / dtime
      runoff_horton(i)   = runoff_horton(i)   * rhoh2o / dtime
      drainage(i)        = drainage   (i)     * rhoh2o / dtime
      drain_rock(i)      = drain_rock (i)     * rhoh2o / dtime
      drain_lat_sl(i,:)  = drain_lat_sl(i,:)  * rhoh2o / dtime
    END DO
    !$claw acc end data
    !$claw acc end data
  END SUBROUTINE calc_soil_hydrology


  !------------------------------------------------------------------------------------------------

  SUBROUTINE infiltration_new(nc,nlat, dt_sat, dtime, max_slope, wtr_pond_max, i_infil, i_shm,     &
                            & oro_stddev, manning, dsoil,                                &
                            & pore_size_index, bclapp, ws_max, fc, pwp, hyd_cond_sat,    &
                            & ws_soil, ice_soil, wtr_pot_scool, water_to_soil, wtr_pond, &
                            & runoff, infilt, pdv_1d_out)


    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max

    INTEGER,  INTENT(in) :: &
      & nc, nlat

    REAL(wp), INTENT(in) :: dtime, &
      & dt_sat,          &         !< Time it takes hydraulic conductivity to reach max [sec]
      & max_slope

    INTEGER,  INTENT(in) :: &
      & i_infil, i_shm

    REAL(wp), INTENT(in) ::      &
      & wtr_pond_max(:),             &
      & oro_stddev(:),           &
      & pore_size_index(:),      &
      & bclapp(:),               &
      & ws_max(:),               &
      & fc(:),                   &
      & pwp(:),                  &
      & hyd_cond_sat(:),         &
      & dsoil(:),                &
      & manning(:),              & 
      & ws_soil(:),              &
      & wtr_pot_scool(:),        &
      & ice_soil(:),             &
      & water_to_soil(:)

    REAL(wp), INTENT(inout) ::      &
      & wtr_pond(:)

    REAL(wp), INTENT(out) ::     &
      & runoff(:), infilt(:),    &
      pdv_1d_out(:)
    ! Local variables
    INTEGER :: jl, i, varindex

    REAL(wp) ::                       &
      ! variables for infiltration calculation
      & sigma_0,                      & !< minimum value
      & sigma_max,                    & !< resolution-dependent maximum
      & steepness,                    & !< sub-grid scale steepness of the orography
      & zbws

    REAL(wp) ::              &
      & zeros(nc),       &
      & wtr_scool(nc),       &
      & ws_incl_ice(nc),     &
      & ice_impede(nc),      &   ! Inhibition of infiltration due to teh presence of ice
      & k_start(nc),         &   ! hydrualic conductivity at the beginning ...
      & k_end(nc),           &   ! ... and the end of the timestep (assuming all water to infiltrate)
      & delta_k(nc),         &   ! ... additional infiltration capacity 
      & water_excess(nc),    &   ! Water that doesn infiltrte immediatly but ponds
      & pond_pot(nc),        &   ! Ponding depth reached during saturation time
      & slope(nc),           &   ! Slope of the contributing areas ... for now use b
      & v_max(nc),           &   ! Avg. velocity of water; vertical crossection 
      & sro_max(nc),         &   ! velocity relative to m2 gridcell area
      & fract_infilt(nc),    &   ! fraction of ponding water that can infiltrate due to inicreasing k
      & infilt_pot(nc)           ! Additional infitration flux due to increasing k

    !--------------------------------------------------------------------------------------------------------------------------------------

    zeros(:) = 0._wp
    wtr_scool(:) = MIN(wtr_pot_scool(:), ws_soil(:))
    ws_incl_ice(:) = ws_soil(:) + ice_soil(:)
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
      !>  1.1 Surface runoff and infiltration
      !   -----------------------------------

      varindex      = 1
      call get_conductivies_diffusivities(i_shm,varindex,            &
                    & dsoil(i), ws_soil(i), ice_soil(i), wtr_scool(i),   &
                    & ws_max(i), fc(i), pwp(i), hyd_cond_sat(i),     &
                    & zeros(i), bclapp(i), pore_size_index(i),       &
                    & ice_impede(i))


      varindex      = 2
      call get_conductivies_diffusivities(i_shm,varindex,            &
                    & dsoil(i), ws_soil(i), ice_soil(i), wtr_scool(i),   &
                    & ws_max(i), fc(i), pwp(i), hyd_cond_sat(i),     &
                    & zeros(i), bclapp(i), pore_size_index(i),       &
                    & k_start(i))


      IF (water_to_soil(i) <= 0._wp) THEN
        k_end(i)      = k_start(i)
      ELSE
        IF (i_infil == 1) THEN
          !!---Assuming the wetting front moves with the speed of ksat
          k_end(i)        = ice_impede(i) * hyd_cond_sat(i)
        ELSE IF (i_infil == 2) THEN
          !!---Assuming that by the end of the timestep k = precip*
          k_end(i)        = MIN(ice_impede(i) * hyd_cond_sat(i), water_to_soil(i)/dtime)
        ELSE IF (i_infil == 3) THEN
          !Assuming that pondig water is euqally distributed throughout the soil layer 
          varindex      = 2
          call get_conductivies_diffusivities(i_shm,varindex,            &
                    & dsoil(i), (ws_soil(i)+ water_to_soil(i)), ice_soil(i), wtr_scool(i),   &
                    & ws_max(i), fc(i), pwp(i), hyd_cond_sat(i),     &
                    & zeros(i), bclapp(i), pore_size_index(i),       &
                    & k_end(i))
        ELSE IF (i_infil == 4) THEN
          !!---Assuming the wetting front moves with the speed of ksat
          !!---and in frozen soils water can enter until pore space of top layer is saturated
          k_end(i)        = MIN(hyd_cond_sat(i),                                  &
                                MAX(ice_impede(i) * hyd_cond_sat(i),              &
!                                   (ws_max(i) - ws_incl_ice(i))/dtime))
                                   (fc(i) - ws_incl_ice(i))/dtime))
        END IF
      END IF

      delta_k(i) = k_end(i) - k_start(i)
      infilt(i)  = MIN(water_to_soil(i), k_start(i)* dtime)
      water_excess(i) = MAX(0._wp, water_to_soil(i) - infilt(i) - wtr_pond_max(i))
!      water_excess(i) = water_to_soil(i) - infilt(i)
      pond_pot(i) = water_excess(i) * dt_sat / dtime


      IF (dtime < dt_sat .AND. (i_infil == 1 .OR. i_infil == 2 .OR. i_infil == 4 )) THEN
        k_end(i)   = k_start(i) +  delta_k(i) * dtime / dt_sat
        delta_k(i) = k_end(i) - k_start(i)
      END IF

      slope(i) = MIN(ABS(zbws), max_slope)
      !IF (hydrotype == LOWLAND_) slope(i) = 0.01_wp

      pdv_1d_out(i) = 1.0_wp
          
      IF (water_excess(i) > 0._wp) THEN
      !!v_ini          = RH**(2/3)*slope**(1/2)*1./manning_n ! RH roughly equal to ponding depth
        v_max(i)       = pond_pot(i)**(2._wp/3._wp)*slope(i)**(1._wp/2._wp)*1._wp / manning(i)
        sro_max(i)     = pond_pot(i) * v_max(i)  ! get from outflux through crossection (side) to flux per grid-box area

        fract_infilt(i) = MAX(0._wp, MIN(1._wp, delta_k(i) / (delta_k(i) + sro_max(i))))
        infilt_pot(i)   = water_excess(i) * fract_infilt(i)
            
        IF (water_excess(i) > (delta_k(i) + sro_max(i)) * dtime) THEN
          infilt(i) = infilt(i) + delta_k(i) * dtime
        ELSE
          infilt(i) = infilt(i) + infilt_pot(i)
        END IF

      END IF

      pdv_1d_out(i) = wtr_pond_max(i)

      wtr_pond(i) = MIN(wtr_pond_max(i), water_to_soil(i) - infilt(i))
      runoff(i)   = MAX(0._wp,    water_to_soil(i) - infilt(i)  - wtr_pond(i))


      !!! Unbedingt noch dahingehend korrigieren, das unter gepondeten Bereichen ziemlich sicher mit ksat infiltriert
      IF (infilt(i) < ice_impede(i) * hyd_cond_sat(i) * dtime .AND. wtr_pond(i) > 0._wp) THEN
        infilt(i)   = MIN(ice_impede(i) * hyd_cond_sat(i) * dtime, infilt(i) + wtr_pond(i))
        wtr_pond(i) = water_to_soil(i) - infilt(i) - runoff(i)
      END IF

    END DO

  END SUBROUTINE infiltration_new

  !------------------------------------------------------------------------------------------------
  SUBROUTINE soilhyd_new(nc, nsoil, dtime,     &
    & i_shm, i_interpol, l_latdrain,   &
    & l_pf_soil,                                               &
    & wtr_pond_max, slope, dsoil,                              &
    & wsat_sl, field_cap_sl, hyd_cond_sat_sl,                  &
    & bclapp_sl, matrix_pot_sl, pore_size_index_sl, p_wilt_sl, &
    & wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl, wtr_pond, &
    & infilt, runoff_dunne, drainage, drain_rock, drain_lat,   &
    & tpe_overflow,             &
    !pdv >
    & weq_soilhyd_err, &
    !pdv <
    & ltpe_closed, ltpe_open)


    USE mo_jsb_math_constants,     ONLY: pi

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
    INTEGER,  INTENT(in)            :: i_shm
    INTEGER,  INTENT(in)            :: i_interpol
    LOGICAL, INTENT(IN)             :: l_latdrain
    LOGICAL, INTENT(IN)             :: l_pf_soil(:)
    REAL(wp), INTENT(IN)            :: wtr_pond_max(:)
    REAL(wp), INTENT(IN)            :: slope(:)
    REAL(wp), INTENT(IN)            :: dsoil(:,:)                 ! soil depth until bedrock within each layer [m]
    REAL(wp), INTENT(IN)            :: wsat_sl(:,:)               ! Maximum storage in each soil layer (porosity) [m]
    REAL(wp), INTENT(IN)            :: field_cap_sl(:,:)          ! field capacity of the soil layer [m]
    REAL(wp), INTENT(IN)            :: hyd_cond_sat_sl(:,:)       ! hydraulic conductivity of saturated soil [m/s]
    REAL(wp), INTENT(IN)            :: bclapp_sl(:,:)             ! exponent B in Clapp and Hornberger
    REAL(wp), INTENT(IN)            :: matrix_pot_sl(:,:)         ! matrix potential [m]
    REAL(wp), INTENT(IN)            :: pore_size_index_sl(:,:)    ! soil pore size distribution index used in Van Genuchten
    REAL(wp), INTENT(IN)            :: p_wilt_sl(:,:)             ! wilting point [m]
    REAL(wp), INTENT(IN)            :: wtr_soil_pot_scool_sl(:,:) ! potentially supercooled water [m]
    REAL(wp), INTENT(IN)            :: ice_soil_sl(:,:)           ! amount of ice in the soil [m]
    REAL(wp), INTENT(INOUT)         :: wtr_soil_sl(:,:)           ! soil moisture of each layer [m]
    REAL(wp), INTENT(INOUT)         :: wtr_pond(:)
    REAL(wp), INTENT(INOUT)         :: infilt(:)                  ! amount of infiltratioj within timestep [m]
    REAL(wp), INTENT(INOUT)         :: runoff_dunne(:)             ! amount of overflow (only diagnostic) [m]
    REAL(wp), INTENT(INOUT)         :: drainage(:)                ! amount of drainage within timestep [m]
    REAL(wp), INTENT(INOUT)         :: drain_rock(:)              ! Needed later for hydrotiles
    REAL(wp), INTENT(INOUT)         :: drain_lat(:,:)             ! Needed later for hydrotiles
    REAL(wp), INTENT(INOUT)         :: tpe_overflow(:)
    REAL(wp), INTENT(INOUT)         :: weq_soilhyd_err(:,:)

    LOGICAL, INTENT(IN) :: ltpe_closed
    LOGICAL, INTENT(IN) :: ltpe_open

    INTEGER  :: i, j                        ! looping index
    INTEGER  :: nsoil_act(nc)               ! number of soil layers above the bedrock boundary

    REAL(wp) :: dttest = 120._wp
    REAL(wp) :: wtr_soil_scool_sl(nc,nsoil) ! supercooled water [m]
    REAL(wp) :: ws_inter_sl(nc,nsoil)       ! soil moisture with infiltration added to first layer
    REAL(wp) :: ws_vol_sl(nc,nsoil)         ! volumetric soil mositure
    REAL(wp) :: ws_rel_sl(nc,nsoil)         ! normalized (~relative) soil mositure
    REAL(wp) :: ice_impede_sl(nc,nsoil)     ! ice impedance factor on all soil layers
    REAL(wp) :: hyd_cond_sl(nc,nsoil)       ! hydraulic conductivity of a layer [m/s]
    REAL(wp) :: diffus_sl(nc,nsoil)         ! diffusivity at mid-depth of the layer[m2/s]
    REAL(wp) :: hyd_cond_li(nc,nsoil)       ! hydraulic conductivity at upper layer interface
    REAL(wp) :: diffus_li(nc,nsoil)         ! diffusivity at upper layer interface
    REAL(wp) :: ws_min_drain(nc)            ! minimum soil moisture for drainage
    REAL(wp) :: ws_rel(nc,nsoil)            ! normalized (~relative) soil mositure
    REAL(wp) :: ws_crit(nc)                 ! soil moisture above which drainage strongly increases [m]
    REAL(wp) :: ws_free(nc)                 ! soil moisture above which lateral drainage apporaches free flow [m]
    REAL(wp) :: drain_crit(nc)              ! drainage above a soil moisture of ws_crit [m/s]
    REAL(wp) :: drain_sl(nc)                ! drainage from the current layer
    REAL(wp) :: remoist(nc)                 ! for TPE
    REAL(wp) :: ws_rsdl_sl(nc,nsoil)        ! residual liquid water that can not move through the soil (supercooled or < pwilt_sl)
    REAL(wp) :: wsat_sl_ni(nc,nsoil)                ! Maximum storage in each soil layer (porosity) (reduced by ice) [m]
    REAL(wp) :: field_cap_sl_ni(nc,nsoil)           ! field capacity of the soil layer (reduced by ice) [m]
    REAL(wp) :: p_wilt_sl_ni(nc,nsoil)              ! wilting point - reduced relative to soil ice fraction [m]

    !$claw acc data &
    !$claw acc present(dsoil, wsat_sl, field_cap_sl, hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, matrix_pot_sl, &
    !$claw acc         pore_size_index_sl, p_wilt_sl, wtr_soil_sl, drainage, tpe_overflow, ice_soil_sl)          &
    !$claw acc pcreate(diffus_sl, ws_min_drain, ws_rel, ws_crit, drain_crit, ws_tmp, drain_sl, soil_above,       &
    !$claw acc         soil_below, zda, zdb, zdc, ztri, diffusivity, zkredu, zdiflog, remoist, zmvg_sl,          &
    !$claw acc         percolation, wslog)

    !$claw acc parallel loop gang vector collapse(2)
    ! Determine the number of soil layers above bedrock


    IF (ANY(wtr_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture (beginning): ', MINVAL(wtr_soil_sl), ' at ', MINLOC(wtr_soil_sl)
      CALL finish ('soilhyd', message_text)
      ! CALL message ('soilhyd', message_text, all_print=.TRUE., level=em_warn)
      ! wtr_soil_sl(:,:) = MAX(wtr_soil_sl(:,:), 0._wp)
    END IF

    !pdv > just here while removing ice from soil props

    wsat_sl_ni(:,:) = MAX(0._wp, wsat_sl(:,:) - ice_soil_sl(:,:))
    field_cap_sl_ni(:,:) = MAX(0._wp, field_cap_sl(:,:) - ice_soil_sl(:,:))
    WHERE(wsat_sl(:,:) > 0._wp)
      p_wilt_sl_ni(:,:) = p_wilt_sl(:,:) * wsat_sl_ni(:,:) / wsat_sl(:,:)
    ELSEWHERE
      p_wilt_sl_ni(:,:) = 0._wp
    END WHERE
    !pdv <


    DO j=1,nc
      nsoil_act(j) = 0
      DO i=1,nsoil
        IF (dsoil(j,i) > 0._wp) nsoil_act(j) = nsoil_act(j) + 1
        ! pdv: supercooled water in the soil should not be able to move vertically
        ws_rsdl_sl(j,i) = MAX(wtr_soil_pot_scool_sl(j,i), p_wilt_sl(j,i))
        ws_rel(j,i)     = 0._wp
      END DO
    END DO

    wtr_soil_scool_sl(:,:) = MIN(wtr_soil_pot_scool_sl(:,:), wtr_soil_sl(:,:))

    ws_inter_sl(:,:) = wtr_soil_sl(:,:)
    ws_inter_sl(:,1) = wtr_soil_sl(:,1) + infilt(:)

    ! Determining hydraulic condiuctivity and diffuisivity on layers and at layer interfaces
    ! Note that the indexing uses the upper interface -> interface Ii = I(i,i-1)
    ! passing nsoil and nsoil_act for looping only over the layers above bedrock ... because why do more ?
    CALL get_conductivies_diffusivities_column(                            &
          & i_shm, i_interpol,                                             &
          & dtime,                                                         &
          & nc,nsoil,                                                      &
          & nsoil_act(:), dsoil(:,:),                                      &
          & ws_inter_sl(:,:), ice_soil_sl(:,:), wtr_soil_scool_sl(:,:),    &
          & wsat_sl(:,:), field_cap_sl(:,:), p_wilt_sl(:,:),               &
          & hyd_cond_sat_sl(:,:), matrix_pot_sl(:,:),                      &
          & bclapp_sl(:,:), pore_size_index_sl(:,:),                       &
          & ice_impede_sl(:,:),                                            &
          & hyd_cond_sl(:,:), diffus_sl(:,:),                              &
          & hyd_cond_li(:,:), diffus_li(:,:)                               &
          )

    !
    !> Drainage pt.1
    !
    ! Drainage of the lowest layer above bedrock. It replaces percolation calculated above
    ! Only done if pf contained to bedrock layers
    !$claw acc loop seq
    ! Above a critical point, fast drainage occurs.
    DO j=1,nc
      IF (.NOT. l_pf_soil(j)) THEN
        !!! pdv 
           ! Hier könnte man sicher auch einfach:
           ! drain_sl(j) = MIN(hyd_cond_sl(j,i) * dtime, MAX(0._wp, wtr_soil_sl(j,i)-ws_rsdl_sl(j,i)))
           ! nehmen, wie fst alle anderen Modelle auch ...
        !!! pdv
        ws_crit(j) = field_cap_sl(j,nsoil_act(j)) * zwdcri     ! soil moisture above which drainage strongly increases

        IF (wtr_soil_sl(j,nsoil_act(j)) < ws_rsdl_sl(j,nsoil_act(j))) THEN
            drain_sl(j) = 0._wp              ! no drainage without mobile water
        ELSE
              ! Slow drainage of mobile water

          IF (field_cap_sl(j,nsoil_act(j)) > ws_rsdl_sl(j,nsoil_act(j))) THEN
            ! slow drainage = slow drainage factor * relative soil moisture
            drain_sl(j) = zdrmin * (wtr_soil_sl(j,nsoil_act(j)) - ws_rsdl_sl(j,nsoil_act(j))) &
                     &  / (field_cap_sl(j,nsoil_act(j)) - ws_rsdl_sl(j,nsoil_act(j)))
          ELSE
            drain_sl(j) = 0._wp
          END IF

          ! Fast drainage of mobile water

          ! Note: here we use wsat_sl (from porosity, reduced by amount of soil ice) and no longer
          !    field_cap_sl as in the above drainage calculations.

          IF (wtr_soil_sl(j,nsoil_act(j)) - wtr_soil_pot_scool_sl(j,nsoil_act(j)) > ws_crit(j)) THEN  ! enough liquid water for fast drainage
            drain_crit(j) = (zdrmax-zdrmin) * ((wtr_soil_sl(j,nsoil_act(j))-wtr_soil_pot_scool_sl(j,nsoil_act(j))-ws_crit(j)) &
             &                             / (wsat_sl(j,nsoil_act(j)) - ws_crit(j)))**zdrexp
            drain_sl(j) = drain_sl(j) &
                      & + MIN(drain_crit(j), (wtr_soil_sl(j,nsoil_act(j))-wtr_soil_pot_scool_sl(j,nsoil_act(j))-ws_crit(j))/dtime)
          END IF
        END IF
        ! drainage is limited by the available water of the layer and impeded by ice

        drain_sl(j)   = MIN(ice_impede_sl(j,nsoil_act(j))*drain_sl(j), &
                      & MAX(0._wp, (wtr_soil_sl(j,nsoil_act(j))-ws_rsdl_sl(j,nsoil_act(j)))/dtime))

        drain_rock(j) = drain_rock(j) + drain_sl(j) * dtime             ! [m/s] -> [m]
        drainage(j)   = drainage(j) + drain_rock(j)     ! and add it to the total drainage
      END IF 
    END DO  


    DO j=1,nc
      IF (wsat_sl(j,1) > 0._wp) THEN !!! Need to figure out if this is necessary, but not sure that all variables exist everywhere !!!
        ! Determination of Vertical movement of water through the soil
        ! Note that concentrations (volumetric soil moisture) are needed for transport
        ws_vol_sl(j,1:nsoil_act(j)) = wtr_soil_sl(j,1:nsoil_act(j)) / dsoil(j,1:nsoil_act(j))    ! Or should this be relative to wsat / fc

        !CALL vert_transport(nc,n,d,alpha,dzf,top_bound, bot_bound, S,C,P,Q, X_old)
        ! *-bound, S are fluxes i.e. x/sec
        ws_vol_sl(j:j,1:nsoil_act(j))  =                              &
        &      vert_transport( 1,nsoil_act(j),dtime,1._wp,            &
                           & dsoil(j:j,1:nsoil_act(j)),               &  ! Or should this be relative to wsat / fc
                           & infilt(j:j)/dtime,                       &
                           & drain_rock(j:j)/dtime,                   &
                           & SPREAD(SPREAD(0._wp, DIM=1, ncopies=1), DIM=2,ncopies=nsoil_act(j)), & ! Ideally one could use evapotranspiration here and get rid of digest evapotr...
                           & SPREAD(SPREAD(1._wp, DIM=1, ncopies=1), DIM=2,ncopies=nsoil_act(j)), &
                           & diffus_li(j:j,1:nsoil_act(j)),           &
                           & hyd_cond_li(j:j,1:nsoil_act(j)),         &
                           & ws_vol_sl(j:j,1:nsoil_act(j))            &
                           & )
        !Comming back to absolute water content
        wtr_soil_sl(j,1:nsoil_act(j)) = ws_vol_sl(j,1:nsoil_act(j)) * dsoil(j,1:nsoil_act(j))
  !pdv > 
        !!!! Testing if waterbalance is ruined.
        DO i = 1, nsoil_act(j)
          IF (wtr_soil_sl(j,i) < 0._wp) THEN
            weq_soilhyd_err(j,i) = weq_soilhyd_err(j,i)+wtr_soil_sl(j,i)
            wtr_soil_sl(j,i) = 0._wp
          END IF
        END DO
      ELSE
        weq_soilhyd_err(j,:) = weq_soilhyd_err(j,:)
      END IF
    END DO
!pdv < 

!pdv >
! Update soil hyd.-porperties
    CALL get_conductivies_diffusivities_column(                          &
        & i_shm, i_interpol,                                             &
        & dtime,                                                         &
        & nc,nsoil,                                                      &
        & nsoil_act(:), dsoil(:,:),                                      &
        & ws_inter_sl(:,:), ice_soil_sl(:,:), wtr_soil_scool_sl(:,:),    &
        & wsat_sl(:,:), field_cap_sl(:,:), p_wilt_sl(:,:),               &
        & hyd_cond_sat_sl(:,:), matrix_pot_sl(:,:),                      &
        & bclapp_sl(:,:), pore_size_index_sl(:,:),                       &
        & ice_impede_sl(:,:),                                            &
        & hyd_cond_sl(:,:), diffus_sl(:,:),                              &
        & hyd_cond_li(:,:), diffus_li(:,:)                               &
        )
!pdv <

    DO j=1,nc
      !
      !> Drainage pt.2
      !
      ! Lateral Drainage (i.e. drainage from upper soil layers)
      !      The idea is that drainage occurs through cracks from an upper soil layer directly
      !      down to the bedrock layer.
      ! ---- In supersaturated, frozen soils, the water above the field capacity is assumed to 
      !      move lateraly until it reaches the stream (or the neighbouring topo-tile)
      ! @todo: This concept needs to be discussed. Should lateral drainage occur at all?
      !
      !$claw acc loop seq
      DO i = 1, nsoil_act(j) - 1
      !> Calculating lateral drainage from all but the lowest layer above bedrock
        !
        !> Drainage pt.2a --- from all soils (impeded in the presence of ice)
!          ws_free(j) = MAX(field_cap_sl_ni(j,i) * zwdfre, EPSILON(1._wp))  ! soil moisture above which drainage strongly increases
!          ! 1.) Caluclate drainage for all soils
          drain_sl(j) = MIN(hyd_cond_sl(j,i) * SIN(slope(i)*pi/2._wp) * dtime, MAX(0._wp, wtr_soil_sl(j,i)-ws_rsdl_sl(j,i)))
          IF ((.NOT. ltpe_closed) .AND. l_latdrain) THEN
            IF (wtr_soil_sl(j,i) > ws_rsdl_sl(j,i)) THEN
              wtr_soil_sl(j,i) = wtr_soil_sl(j,i) - drain_sl(j) ! remove drainage of the layer from soil moisture
              drain_lat(j,i)   = drain_lat(j,i)  + drain_sl(j)  ! When taking this for hydrotiles make sure its no longer added to drainage !
              drainage(j)      = drainage(j) + drain_sl(j)      ! and add it to the total drainage
            END IF
          END IF
      END DO

      ! 3. Correct percolation if maximum moisture is exceeded, and correct soil moisture of the
      !    above layer accordingly

      ! Overflow is not allowed from layers below the surface layer.
      ! Instead it is assumed that water piles upwards from saturated layers.
      !$claw acc loop seq
      DO i = nsoil_act(j), 2, -1
        IF (.NOT. ltpe_closed) THEN
          IF (wtr_soil_sl(j,i) > wsat_sl_ni(j,i)) THEN
            ! soil moisture above field capacity:
            wtr_soil_sl(j,i-1) = wtr_soil_sl(j,i-1) + wtr_soil_sl(j,i) - wsat_sl_ni(j,i)
            wtr_soil_sl(j,i)   = wsat_sl_ni(j,i)
          END IF
        ELSE
          IF (wtr_soil_sl(j,i) > field_cap_sl_ni(j,i)) THEN
            ! redirect overflow from drainage to overflow pool
            tpe_overflow(j) = tpe_overflow(j) + wtr_soil_sl(j,i) - field_cap_sl_ni(j,i)
            wtr_soil_sl(j,i) = field_cap_sl_ni(j,i)
          ELSE ! wtr_soil_sl <= field_cap_sl_ni
            ! refill soil moisture from overflow pool
            remoist(j) = MIN(tpe_overflow(j), field_cap_sl_ni(j,i) - wtr_soil_sl(j,i))
            wtr_soil_sl(j,i) = wtr_soil_sl(j,i) + remoist(j)
            tpe_overflow(j) = tpe_overflow(j) - remoist(j)
          END IF
        END IF
      END DO

      ! 4. Add overflow of the surface layer to drainage
      ! Overflow from the surface layer enters drainage, resp. the overflow pool for TPE closed.
      ! This should not happen, this is just for security.
      IF (.NOT. ltpe_closed) THEN
        IF (wtr_soil_sl(j,1) > wsat_sl_ni(j,1)) THEN
          wtr_pond(j)    = wtr_pond(j) + wtr_soil_sl(j,1)-wsat_sl_ni(j,1)
          infilt(j)      = infilt(j) - wtr_soil_sl(j,1) + wsat_sl_ni(j,1)
          IF (wtr_pond(j) > wtr_pond_max(j)) THEN
            runoff_dunne(j) = runoff_dunne(j) + wtr_pond(j) - wtr_pond_max(j)
            !drainage(j)    = drainage(j)    + wtr_pond(j) - wtr_pond_max(j)
            wtr_pond(j)    = wtr_pond_max(j)
          END IF
          wtr_soil_sl(j,1) = wsat_sl_ni(j,1)
        END IF
      ELSE      ! TPE closed case
        wtr_soil_sl(j,1) = wtr_soil_sl(j,1) + tpe_overflow(j)
        IF (wtr_soil_sl(j,1) > field_cap_sl_ni(j,1)) THEN
          ! redirect overflow to overflow pool
          tpe_overflow(j)  = wtr_soil_sl(j,1) - field_cap_sl_ni(j,1)
          wtr_soil_sl(j,1) = field_cap_sl_ni(j,1)
        END IF
      END IF

      ! Modification for TPE (open)
      IF (ltpe_open) THEN
        !$claw acc loop seq
        DO i = 1, nsoil
          !IF (i == nsoil - 2 .OR. i == nsoil - 1 .OR. i == nsoil) THEN   ! special adaptation of Shirisha
          IF (i == nsoil-1 .OR. i == nsoil) THEN
            IF (dsoil(j,i) > 0._wp) THEN
             ! keep the two lowest soil layers always very wet
             wtr_soil_sl(j,i) = MAX(wtr_soil_sl(j,i),field_cap_sl_ni(j,i) * 0.9_wp)
            END IF
          END IF
        END DO
      END IF

    END DO

    IF (ANY(wtr_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture (end): ', MINVAL(wtr_soil_sl), ' at ', MINLOC(wtr_soil_sl)
      CALL finish ('soilhyd', message_text)
      ! CALL message ('soilhyd', message_text, all_print=.TRUE., level=em_warn)
      ! wtr_soil_sl(:,:) = MAX(wtr_soil_sl(:,:), 0._wp)
    END IF

    !$claw acc end parallel
    !$claw acc end data

  END SUBROUTINE soilhyd_new

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
    REAL(wp) :: fixed(nc)          ! amount of water below the wilting point, not available for plants
    REAL(wp) :: trans_layer(nc)    ! transpiration from the layer
    REAL(wp) :: ws(nc)             ! water and ice within all soil layers
    REAL(wp) :: root_depth(nc)     ! Total depth of rootzone
    REAL(wp) :: deficit_evapotrans(nc) ! water deficit due to evapotranspiration
    REAL(wp) :: wtr_soil_sl_new(nc)    ! temporary water storage

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


    IF (ANY(wtr_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture (beginning): ', MINVAL(wtr_soil_sl), ' at ', MINLOC(wtr_soil_sl)
      !CALL finish ('digest_evapotrans', message_text)
       CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
       wtr_soil_sl(:,:) = MAX(wtr_soil_sl(:,:), EPSILON(1._wp))
    END IF

    DO i = 1, nc
      deficit_evapo(i) = 0._wp   ! water missing for evaporation due to limited soil water
      IF (evapo_soil_in_m(i) < 0._wp) THEN                       ! pos. evaporation flux (not dew)
        wtr_soil_sl(i,1) = wtr_soil_sl(i,1) + evapo_soil_in_m(i) ! take water from uppermost soil layer
        IF (wtr_soil_sl(i,1) < EPSILON(1._wp)) THEN              ! if there is not enough soil water
          deficit_evapo(i) = wtr_soil_sl(i,1) - EPSILON(1._wp)   ! this is the amount of water missing
          wtr_soil_sl(i,1) = EPSILON(1._wp)                      ! soil moisture cannot be negative
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
            ! @todo: should we use here p_wilt_sl (instead of field_cap_sl*zwilt)?
            fixed(i) = MAX(EPSILON(1._wp), field_cap_sl(i,jk) * zwilt * rootfract(i)) ! amount of water in the root zone, not available for plants

            IF (wtr_soil_sl(i,jk) * rootfract(i) >= fixed(i)) THEN              ! there is water in the root zone available for plants
              remaining(i) = wtr_soil_sl(i,jk) * rootfract(i) + trans_layer(i)  ! water that would remain after transpiration
              IF (remaining(i) < fixed(i)) THEN                                 ! transpiration exceeds available water
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

      ! Deficit should be taken from layers starting at the surface to accelerate the feedback on evapotranspiration
      DO jk=1, nsoil
        ! @todo: same as above: shouldn't we use here p_wilt_l (instead of filed_cap_l*zwilt)?
        fixed(i) = MAX(EPSILON(1._wp), field_cap_sl(i,jk) * zwilt)
        IF (deficit_evapotrans(i) < 0._wp .AND. wtr_soil_sl(i,jk) > fixed(i)) THEN
          ! reduce soil moisture and ice relative to the water content of each layer
          wtr_soil_sl_new(i)    = MAX(wtr_soil_sl(i,jk) + deficit_evapotrans(i), fixed(i))
          deficit_evapotrans(i) = deficit_evapotrans(i) + wtr_soil_sl(i,jk) - wtr_soil_sl_new(i)
          wtr_soil_sl(i,jk)     = wtr_soil_sl_new(i)
        END IF
      END DO

      ! Note: ws here refers to total soil column, and not to root zone (as ws above).
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
      !WRITE (message_text,*) 'negative soil ice: ', MINVAL(ice_soil_sl), ' at ', MINLOC(ice_soil_sl)
      !CALL finish ('digest_evapotrans', message_text)
       CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
       ice_soil_sl(:,:) = MAX(ice_soil_sl(:,:), EPSILON(1._wp))
    END IF

    ! make sure there is no negative soil moisture and soil moisture does not exceed field capacity
    IF (ANY(wtr_soil_sl(:,:) < 0._wp)) THEN
      WRITE (message_text,*) 'negative soil moisture (end): ', MINVAL(wtr_soil_sl), ' at ', MINLOC(wtr_soil_sl)
      !CALL finish ('digest_evapotrans', message_text)
       CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
       wtr_soil_sl(:,:) = MAX(wtr_soil_sl(:,:), EPSILON(1._wp))
    END IF
    ! pdv > Idea is to allow supersaturation (i.e. fc > ws > wsat)
    !IF (ANY(wtr_soil_sl(:,:) > field_cap_sl(:,:) + zfcmin)) THEN
    !  max2d(:) = MAXLOC(wtr_soil_sl(:,:) - field_cap_sl(:,:))
    !  WRITE (message_text,*) 'liquid water content above field capacity at ', max2d(:), &
    !    &   ': wtr_soil_sl: ', wtr_soil_sl(max2d(1),max2d(2)),                          &
    !    &   ' field_cap_sl: ', field_cap_sl(max2d(1),max2d(2))
    !  CALL finish ('digest_evapotrans', message_text)
    !END IF
    !pdv <
#endif

  END SUBROUTINE digest_evapotrans
  !------------------------------------------------------------------------------------------------



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

  FUNCTION vert_transport(nc,n,dt,alpha,dzf,top_bound, bot_bound, S,C,P,Q, X_old) RESULT(X_new)

  ! INPUT
  INTEGER, INTENT(in) :: nc,n                       ! vector length, number of levels
  REAL(wp), Intent(in):: dt, alpha,  &              ! time step
                       & top_bound(nc), bot_bound(nc), &    ! boundary conditions
                       & dzf(nc,n), C(nc,n), S(nc,n), Q(nc,n), P(nc,n), X_old(nc,n)

  ! OUTPUT
  REAL(wp)           :: X_new(nc,n)

  ! LOCAL
  INTEGER, PARAMETER :: NMAX=500
  INTEGER            :: i,j
  REAL(wp)           :: dzh(nc,n), matrix(nc,n,4), &
                      & bet(nc), gam(nc,NMAX)

  dzh(:,:) = 0.0                             ! thickness of half levels
  DO j=2,n
    dzh(:,j) = (dzf(:,j-1) + dzf(:,j)) / 2._wp
  END DO

  ! Note that P & Q refer to the inter face above a given lvl not below.
  ! Hence for j == 1, there is no P(j) and Q(j)
  ! Make sure that this is consistent with the definition of input params.

  DO j=1,n
    IF(j .eq. 1) THEN
      matrix(:,j,1) =  0.0_wp
      matrix(:,j,2) =                   P(:,j+1) / dzh(:,j+1) + dzf(:,j) * C(:,j) / (alpha * dt)
      matrix(:,j,3) = -P(:,j+1) / dzh(:,j+1)
      matrix(:,j,4) =       - Q(:,j+1) + dzf(:,j) * S(:,j) + dzf(:,j) * C(:,j) * X_old(:,j) / (alpha * dt) + top_bound(:)
    ELSE IF((j .NE. 1) .AND. (j .NE. n)) THEN
      matrix(:,j,1) = -P(:,j) / dzh(:,j)
      matrix(:,j,2) =  P(:,j) / dzh(:,j) +  P(:,j+1) / dzh(:,j+1) + dzf(:,j) * C(:,j) / (alpha * dt)
      matrix(:,j,3) = -P(:,j+1) / dzh(:,j+1)
      matrix(:,j,4) =  Q(:,j) - Q(:,j+1) + dzf(:,j) * S(:,j) + dzf(:,j) * C(:,j) * X_old(:,j) / (alpha * dt)
    ELSE
      matrix(:,j,1) = -P(:,j) / dzh(:,j)
      matrix(:,j,2) =  P(:,j) / dzh(:,j)                      + dzf(:,j) * C(:,j) / (alpha * dt)
      matrix(:,j,3) =  0.0_wp
      matrix(:,j,4) =  Q(:,j)          + dzf(:,j) * S(:,j) + dzf(:,j) * C(:,j) * X_old(:,j) / (alpha * dt) - bot_bound(:)
    END IF
  END DO

  IF( ANY(matrix(:,1,2) .EQ. 0.0)) THEN
    WRITE (message_text,*) 'tridiag: something went terribly wrong ... rewrite equations!', MINVAL(C(:,1)), MINVAL(dzf(:,1))
    CALL finish('vert_transport', message_text)
  END IF
  bet(:)  = matrix(:,1,2)
  X_new(:,1) = matrix(:,1,4)/bet(:)
  DO j = 2, n
    gam(:,j) = matrix(:,j-1,3) /bet(:)
    bet(:)    = matrix(:,j,2) - matrix(:,j,1) * gam(:,j)
    IF( ANY(bet(:) .EQ. 0.0)) THEN
      WRITE (message_text,*) 'tridiag: matrix couldn t be solved !'
      CALL finish('vert_transport', message_text)
    END IF
    X_new(:,j)   = (matrix(:,j,4) - matrix(:,j,1) * X_new(:,j-1)) / bet(:)
  END DO
  DO j  = n-1, 1, -1
    X_new(:,j)= X_new(:,j) - gam(:,j+1) * X_new(:,j+1)
  END DO

  END FUNCTION vert_transport

  SUBROUTINE get_conductivies_diffusivities_column(     &
            & i_shm, i_interpol,                 &
            & dt,                                &
            & nc,nsoil,nsoil_act,dsoil,          &
            & ws, is, wscools,                   &
            & wsat, fc, pwp,                     &
            & k_sat, mpot_sat, bclapp, pore_size_index,     &
            & ice_impede, K, D, K_inter, D_inter &
            )

    ! arguments

    INTEGER,  INTENT(in)     :: i_shm, i_interpol   ! Approaches to determine conductivites&difusivities (KD) and how to derive values at layer interface (AVG scheme)
    REAL(wp), INTENT(in)     :: dt                      ! time step length
    INTEGER,  INTENT(in)     :: nc                      ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of below ground layers
    INTEGER,  INTENT(in)     :: nsoil_act(:)            ! number of soil layers (above bedrock)
    REAL(wp), INTENT(in)     :: dsoil(:,:)              ! soil depth (until bedrock) per layer
    REAL(wp), INTENT(in)     :: ws(:,:)                 ! water content
    REAL(wp), INTENT(in)     :: is(:,:)                 ! ice content
    REAL(wp), INTENT(in)     :: wscools(:,:)            ! supercooled water
    REAL(wp), INTENT(in)     :: wsat(:,:)               ! soil depth (until bedrock) per layer
    REAL(wp), INTENT(in)     :: fc(:,:)                 ! field capacity [m]
    REAL(wp), INTENT(in)     :: pwp(:,:)                ! residual water (until bedrock) per layer
    REAL(wp), INTENT(in)     :: k_sat(:,:)              ! hydraaulic conductivity
    REAL(wp), INTENT(in)     :: mpot_sat(:,:)           ! matric potential
    REAL(wp), INTENT(in)     :: bclapp(:,:)             ! clapp & Hornberger exponent
    REAL(wp), INTENT(in)     :: pore_size_index(:,:)    ! Pore size distribution index

    REAL(wp), INTENT(out)    :: ice_impede(nc,nsoil)    ! Impedence factor to account for ice blocking flowpaths
    REAL(wp), INTENT(out)    :: K(nc,nsoil)             ! hydrualic condctivity of the soil layer [m]
    REAL(wp), INTENT(out)    :: D(nc,nsoil)             ! diffusivity of the soil layer [m]
    REAL(wp), INTENT(out)    :: K_inter(nc,nsoil)       ! hydrualic condctivity at (upper) layer interface [m]
    REAL(wp), INTENT(out)    :: D_inter(nc,nsoil)       ! diffusivity at (upper) layer interface [m]

    ! local variables

    INTEGER  :: i,j, varindex

    REAL(wp) :: outvar(nc,nsoil)    ! Impedence factor to account for ice blocking flowpaths

    varindex       = 1
    call get_conductivies_diffusivities(i_shm,varindex,                      &
                    & dsoil(:,:), ws(:,:), is(:,:), wscools(:,:),       &
                    & wsat(:,:), fc(:,:), pwp(:,:), k_sat(:,:),         &
                    & mpot_sat(:,:), bclapp(:,:), pore_size_index(:,:), &
                    & ice_impede(:,:))

    varindex       = 2
    call get_conductivies_diffusivities(i_shm,varindex,    &
                    & dsoil(:,:), ws(:,:), is(:,:), wscools(:,:),       &
                    & wsat(:,:), fc(:,:), pwp(:,:), k_sat(:,:),         &
                    & mpot_sat(:,:), bclapp(:,:), pore_size_index(:,:), &
                    & K(:,:))

    varindex       = 3
    call get_conductivies_diffusivities(i_shm,varindex,    &
                    & dsoil(:,:), ws(:,:), is(:,:), wscools(:,:),       &
                    & wsat(:,:), fc(:,:), pwp(:,:), k_sat(:,:),         &
                    & mpot_sat(:,:), bclapp(:,:), pore_size_index(:,:), &
                    &  D(:,:))

    K_inter(:,:)    = 0._wp
    D_inter(:,:)    = 0._wp

    !!! making sure percolation from top layer is limited to water content !!!
    K(:,1) = MIN(K(:,1), MAX(0._wp,(ws(:,1)-pwp(:,1)) / dt))

    ! Determining K & D on layer interface
    ! Possibility
    ! 1) Upstream mean [K_lev; Dmax(x_lev, x_lev-1)}
    ! 2) Adjusted arithmetic mean A (of k, D)
    ! Note that in the vertical transport scheme the notation refers to the upper interface of a level.
    ! Hence the loop goes from nsoil_act to 2
    DO j = 1, nc
      DO i = nsoil_act(j), 2, -1
        IF (i_interpol == 1) THEN ! upstream
          K_inter(j,i)  =   K(j,i-1)
          D_inter(j,i)  =   MAX(D(j,i),D(j,i-1))
!          D_inter(j,i)  =   MIN(D(j,i),D(j,i-1))
        ELSE IF (i_interpol == 2) THEN !  Adapted arithmetic mean
          IF (K(j,i-1) > K(j,i)) THEN ! MOVEMENT OF WETTING FRONTS LIMITED BY SATURATION OF SUBJACENT LAYER. ELSE ... 
            K_inter(j,i)  =   (K(j,i)*dsoil(j,i)+K(j,i-1)*dsoil(j,i-1))/(dsoil(j,i)+dsoil(j,i-1))
          ELSE ! ... PREVENT FLUX FROM COMPLETLY DRYING OUT UPPER LAYER
            K_inter(j,i)  = K(j,i-1)
          END IF
          D_inter(j,i)  =   (D(j,i)*dsoil(j,i)+D(j,i-1)*dsoil(j,i-1))/(dsoil(j,i)+dsoil(j,i-1))
        END IF
      END DO
    END DO

  END SUBROUTINE get_conductivies_diffusivities_column

  ! ===============================================================================================================================
  PURE ELEMENTAL SUBROUTINE get_conductivies_diffusivities(            &
            & i_shm, i_delieveroo,               &
            & dsoil,                             &
            & ws, is, wscools,                   &
            & wsat, fc, pwp,                     &
            & k_sat, mpot_sat, bclapp, pore_size_index,     &
            & outvar                             &
            )

    ! arguments

    INTEGER,  INTENT(in)     :: i_shm, i_delieveroo! soil-hydrology model used, variable requested
    REAL(wp), INTENT(in)     :: dsoil              ! soil depth (until bedrock) per layer
    REAL(wp), INTENT(in)     :: ws                 ! water content
    REAL(wp), INTENT(in)     :: is                 ! ice content
    REAL(wp), INTENT(in)     :: wscools            ! supercooled water
    REAL(wp), INTENT(in)     :: wsat               ! soil depth (until bedrock) per layer
    REAL(wp), INTENT(in)     :: fc                 ! field capacity [m]
    REAL(wp), INTENT(in)     :: pwp                ! residual water (until bedrock) per layer
    REAL(wp), INTENT(in)     :: k_sat              ! hydraaulic conductivity
    REAL(wp), INTENT(in)     :: mpot_sat           ! matric potential
    REAL(wp), INTENT(in)     :: bclapp             ! clapp & Hornberger exponent
    REAL(wp), INTENT(in)     :: pore_size_index    ! Pore size distribution index


    REAL(wp), INTENT(out)    :: outvar             ! Impedence factor to account for ice blocking flowpaths

    REAL(wp) ::  ws_vol, ice_vol, scool_vol, &
             &   ws_rel_vol,  ws_range_vol,  &
             &   wsat_vol, fc_vol, pwp_vol,  &
             &   ws_free_vol, ice_impede

    REAL(wp) ::  ck                 ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)


    REAL(wp) ::  nvgn, mvgn,       &
             &   mpot_sat_correct, &
             &   K, D, mpot


    ck           = 8.0_wp           ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    ! switching to volumetric quantities
    IF(dsoil> 0._wp) THEN
      ws_vol    =  ws      / dsoil
      ice_vol   =  is      / dsoil
      scool_vol =  wscools / dsoil
      wsat_vol  =  wsat    / dsoil
      fc_vol    =  fc      / dsoil
      pwp_vol   =  pwp     / dsoil
    ELSE
      ws_vol    =  0._wp
      ice_vol   =  0._wp
      scool_vol =  0._wp
      wsat_vol  =  0._wp
      fc_vol    =  0._wp
      pwp_vol   =  0._wp
    END IF


    IF (i_delieveroo == 3 .OR. i_delieveroo == 4) THEN
      IF (mpot_sat > 0._wp) THEN
          mpot_sat_correct = -mpot_sat
      ELSE
          mpot_sat_correct = mpot_sat
      END IF
    END IF

    IF (i_shm == 3) THEN
      nvgn = pore_size_index + 1._wp
      mvgn = pore_size_index / nvgn
    END IF

    IF (i_shm==2) THEN ! Clapp & Horberger (Campbell)
      ws_free_vol   = MAX(0._wp,ws_vol-scool_vol)
      ws_range_vol  = wsat_vol
    ELSE  ! Brooks & Coorey; Van Genuchten
      ws_free_vol   = MAX(0._wp,ws_vol-MAX(scool_vol,pwp_vol))
      ws_range_vol  = fc_vol - pwp_vol
    END IF
    IF (ws_range_vol > 0._wp) THEN
      ws_rel_vol  = ws_free_vol / ws_range_vol
    ELSE
      ws_rel_vol  = 0._wp
    END IF

    ! Making sure ws_rel is in the pyhsical range
    ! This is especially important when infiltration is added to the first layer
    ws_rel_vol = MIN(1._wp, MAX(0._wp, ws_rel_vol))

    ! Reduce percolation velocity in presence of ice
    ! But what about diffucivity .... ice seems to attract water ?
    IF (wsat_vol > 0._wp) THEN
        ice_impede  = 10._wp**(-6._wp * (ice_vol + scool_vol) / wsat_vol)
    ELSE
       ice_impede  = 1._wp
    END IF

    IF (i_delieveroo == 1) THEN
      outvar = ice_impede
    END IF

    IF (i_delieveroo == 2) THEN
      ! Determining conductivit K
      IF (i_shm==1 .OR. i_shm==2) THEN  !1=Brooks & Correy or 2 = Clapp & Horberger (Campbell)
        K =  k_sat * ws_rel_vol**(2._wp * bclapp + 3._wp)

      ELSE IF (i_shm==3) THEN ! Van Genuchten
        ws_rel_vol = MIN(0.9999_wp, ws_rel_vol) ! doesn't work for wrel == 1; (1-ws_rel(:,:)**(1/mvgn(:,:)))**(-mvgn)

        K = k_sat * ws_rel_vol**(0.5_wp) &
        & *(1._wp-(1._wp-ws_rel_vol**(1._wp/mvgn))**mvgn)**2._wp
      END IF
      K = K * ice_impede
      outvar = K
    END IF

    IF (i_delieveroo == 3) THEN
      ! Determining diffusivity D
      IF (i_shm==1 .OR. i_shm==2) THEN  !1=Brooks & Correy or 2 = Clapp & Horberger (Campbell)
        IF(ws_range_vol > 0._wp) THEN
          D = -k_sat * mpot_sat_correct * bclapp / ws_range_vol &
            & * ws_rel_vol**(1._wp * bclapp + 2._wp)
        ELSE
          D = 0._wp
        END IF

      ELSE IF (i_shm==3) THEN ! Van Genuchten
        IF(ws_range_vol > 0._wp) THEN
          ws_rel_vol = MIN(0.9999_wp, ws_rel_vol) ! doesn't work for wrel == 1; (1-ws_rel(:,:)**(1/mvgn(:,:)))**(-mvgn)

          K = k_sat * ws_rel_vol**(0.5_wp) &
            & *(1._wp-(1._wp-ws_rel_vol**(1._wp/mvgn))**mvgn)**2._wp

          D = -K * mpot_sat_correct * bclapp / ws_range_vol &
            & * ws_rel_vol**(1._wp/mvgn) * (1._wp-ws_rel_vol**(1._wp/mvgn))**(-mvgn)
        ELSE
          D = 0._wp
        END IF
      END IF
      D = D * ice_impede
      outvar = D
    END IF

    IF (i_delieveroo == 4) THEN
      ! Determin matric potential
      !!! Formula from Stuurop 2021 (https://doi.org/10.1016/j.coldregions.2021.103456) adjusted from Zhang 2007 (https://doi.org/10.1175/JHM605.1)
      IF (ws_rel_vol > 0._wp) THEN
        IF (i_shm == 1 .OR. i_shm==2) THEN
          mpot = mpot_sat_correct * (ws_rel_vol**(-bclapp))
        ELSE
          mpot = mpot_sat_correct * (ws_rel_vol**(-1._wp/mvgn) - 1._wp)**(1._wp/nvgn)
        END IF
        mpot = mpot * ((1._wp + ck * ice_vol)**2._wp)
        mpot = MAX(-100._wp,  mpot)
      ELSE
        mpot = -100._wp
      END IF
      outvar = mpot
    END IF


  END SUBROUTINE get_conductivies_diffusivities


#endif
END MODULE mo_hydro_process
