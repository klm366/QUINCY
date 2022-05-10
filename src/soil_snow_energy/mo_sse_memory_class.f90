!> Contains the memory class for the soil and snow energy process.
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
MODULE mo_sse_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE, BARE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes SSE_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  ! ======================================================================================================= !
  !>Type definition for soil_snow_energy memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_sse_memory

    !------------------------------------------------------------------------------------------------------ !
    !>1.0 jsbach4
    !>
    
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      & vol_heat_cap,         & !< Volumetric heat capacity of soil solids (from input file) [J m-3 K-1]
      & heat_cond,            & !< Thermal conductivity of soil solids (from input file) [J m-1 s-1 K-1]
      & vol_heat_cap_snow,    & !< (dynamic) Volumetric heat capacity of snow (total depth) [J m-3 K-1]
      & heat_cond_snow,       & !< Thermal conductivity of snow [J m-1 s-1 K-1]
      & thermal_diffusivity     !< Thermal diffusivity of soil [m2 s-1]
    !
    TYPE(t_jsb_var_real3d) :: &
      & vol_heat_cap_sl,      & !< (dynamic) Volumetric heat capacity of dry soil [J m-3 K-1]
      & heat_cond_sl,         & !< (dynamic) Thermal conductivity of soil [J m-1 s-1 K-1]
      & t_soil_sl,            & !< Temperature of soil layers [K]
      & t_soil_acoef,         & !< Soil A coefficient of Richtmeyer and Morton scheme
      & t_soil_bcoef,         & !< Soil B coefficient of Richtmeyer and Morton scheme
      & snow_depth_sl,        & !< Depth of snow layers [m]
      & t_snow,               & !< Temperature of snow layers [K]
      & t_snow_acoef,         & !< Snow A coefficient of Richtmeyer and Morton scheme
      & t_snow_bcoef            !< Snow B coefficient of Richtmeyer and Morton scheme
    TYPE(t_jsb_var_real2d) :: &
      & hcap_grnd,            & !< Heat capacity of the ground (uppermost layer)        [J m-2 K-1]
      & hcap_grnd_old,        & !< Heat capacity of the ground (uppermost_layer), old   [J m-2 K-1]
      & grnd_hflx,            & !< ground heat flux
      & grnd_hflx_old,        & !< ground heat flux (old)
      & thaw_depth              !< Thawing depth [m]


    !------------------------------------------------------------------------------------------------------ !
    !>2.0 quincy
    !>
    ! present in HYDRO_ or not used outside quincy SSE_ tasks, maybe used in SB_
    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth      , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    root_depth      , & !< actual rooting depth [m]
                                    t_srf_old       , & !< previous surface temperature [K]
                                    zril_old            !< previous timestep's Reynold's number

    ! present in HYDRO_ or not used outside quincy SSE_ tasks, maybe used in SB_
    TYPE(t_jsb_var_real3d)          :: &
                                    ! soil_depth_sl, &            !< Width ! of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_width_sl, &        !< width of each soil layer that can be saturated with water, identical with soil_depth_sl (until bedrock) [m]
                                    soil_lay_depth_center_sl, & !< depth at the center of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_ubound_sl, & !< depth at the upper bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_lbound_sl    !< depth at the lower bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]

    ! present in HYDRO_ or not used outside quincy SSE_ tasks, maybe used in SB_
    ! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d)          :: &
                                    interception,     & !< surface interception evaporation  [kg m-2 s-1]
                                    pet,              & !< potential evaporation  [kg m-2 s-1]
                                    evaporation,      & !< surface evaporation  [kg m-2 s-1]
                                    runoff,           & !< surface runoff  [kg m-2 s-1]
                                    drainage,         & !< drainage        [kg m-2 s-1]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage [1 s-1] 

    TYPE(t_jsb_var_real2d)          :: &
                                    ! present in HYDRO_, not used outside quincy SSE_
                                    w_skin,           & !< Water content in skin reservoir                                          [m]
                                    ! present in HYDRO_, used in ASSIMI_ and VEG_
                                    ! w_soil_root,      & !< Water content in root zone of the soil                                   [m]
                                    ! w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
                                    ! w_soil_root_pwp,  & !< Water content at permanent wilting point in root zone of the soil        [m]
                                    ! not present in jsb4, calculated in quincy SSE_, calc & var moved to jsb4 SSE_ "soil_and_snow_temperature" task
                                    t_soil_root         !< Temperature of fine roots in the soil                                    [K]

    TYPE(t_jsb_var_real3d)          :: &
                                    ! used in quincy SSE_ & SB_, not present in jsb4
                                    percolation_sl,   & !< fraction of mass transferred to below layer                              [fraction s-1]
                                    drainage_sl!,      & !< drainage of soil layers                                                  [kg m-2 s-1]
                                    ! present in HYDRO_,  not used outside quincy SSE_
                                    ! w_soil_sl,        & !< Water content in soil layers                                             [m]
                                    ! w_soil_sat_sl,    & !< Water content of soil layers at saturation 
                                    !                     !< (derived from volumetric capacity at saturation)                         [m]
                                    ! w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                                    !                     !< (derived from volumetric field capacity)                                 [m]
                                    ! w_soil_pwp_sl,    & !< Water content of soil layers at permanent wilting point
                                    !                     !< (derived from vol. perm wilt point)                                      [m]
                                    ! not present in HYDRO_, needed in jsb4 HYDRO_ for quincy for calc of w_soil_root_pot
                                    ! saxtonA,          & !< coefficient in moisture-tension relationship
                                    ! saxtonB,          & !< coefficient in moisture-tension relationship
                                    ! saxtonC ,          & !< exponent of moisture-conductivity relationship
                                    ! not present in HYDRO_, needed only in quincy SSE_
                                    ! kdiff_sat_sl        !< saturated hydraulic conductivity                                         [m s-1]

    ! not present in HYDRO_ and not used outside quincy SSE_
    ! heat fluxes
    TYPE(t_jsb_var_real2d)          :: &                              
                                    sensible_heat_flx     , &  !< sensible heat flux [W m-2]
                                    latent_heat_flx       , &  !< latent heat flux [W m-2]
                                    ground_heat_flx       , &  !< ground heat flux [W m-2]
                                    ground_heat_flx_old        !< ground heat flux of previous timestep [W m-2]

    ! present in HYDRO_ and used by quincy in VEG_
    ! ! transpiration
    ! TYPE(t_jsb_var_real2d)          :: &
    !                                 transpiration       !< transpiration in [kg m-2 s-1]

    ! moved to HYDRO_, needed outside quincy SSE_
    !Soil Water addons for QUINCY in HYDRO
    ! TYPE(t_jsb_var_real2d)          :: &
    !                                 w_soil_root_theta, & !< root zone integrated plant available water, fraction of maximum 
    !                                 w_soil_root_pot      !< soil water potential [MPa]

    ! not present in HYDRO_, not used outside quincy SSE_
    ! TYPE(t_jsb_var_real3d)          :: &
    !                                 w_soil_pot_sl        !< soil water potential per layer [MPa]

    ! present in jsb4 SSE_ too
    ! soil energy 
    ! TYPE(t_jsb_var_real3d)          :: &
    !                                 t_soil_sl            !< soil temperature [K] 

    ! not present in HYDRO_, not used outside quincy SSE_
    ! soil physical properties
    TYPE(t_jsb_var_real3d)          :: & 
                                    ! heat_capa_sl,         &  !< soil heat capacity [J kg-1]
                                    bulk_dens_sl!,         &  !< soil bulk density [kg m-3]
                                    ! therm_cond_sl            !< soil thermal conductivity [W m-2 K-1]

    TYPE(t_jsb_var_real3d)          :: & 
                                    ! not present in HYDRO_, used outside quincy SSE_
                                    sand_sl             , &  !< soil sand proportion of mineral soil
                                    silt_sl             , &  !< soil silt proportion of mineral soil
                                    clay_sl             , &  !< soil clay proportion of mineral soil(re-calculated as: clay = 1.0 -sand -silt)
                                    ! not present in HYDRO_, used in quincy SSE_ & SB_
                                    volume_min_sl       , &  !< soil mineral volumetric fraction
                                    qmax_org_min_sl     , &  !< maximum OM sorption capacity of mineral soil, [mol C kg-1 mineral soil]
                                    qmax_po4_min_sl     , &  !< maximum po4 sorption capacity of mineral soil, [mol P kg-1 mineral soil]
                                    qmax_nh4_min_sl     , &  !< maximum nh4 sorption capacity of mineral soil, [mol N kg-1 mineral soil]
                                    qmax_po4_om_sl           !< maximum po4 sorption capacity of SOM, [mol P kg-1 SOM]

  CONTAINS
    PROCEDURE :: Init => Init_sse_memory
  END TYPE t_sse_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_memory_class'


CONTAINS


  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: soil_snow_energy
  !>
  SUBROUTINE Init_sse_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_JSBACH, MODEL_JSBACH_QUINCY, MODEL_QUINCY
    USE mo_jsb_class,         ONLY: Get_model
    dsl4jsb_Use_config(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_sse_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix
    CHARACTER(len=*),    INTENT(in)            :: suffix
    INTEGER,             INTENT(in)            :: lct_ids(:)
    INTEGER,             INTENT(in)            :: model_id

    dsl4jsb_Def_config(SSE_)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface, soil_e, snow_e      ! Vertical grids
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_snow_e                 ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (SIZE(lct_ids) < 1) CALL finish(routine, 'lcd_ids not specified')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    dsl4jsb_Get_Config(SSE_)
    table         = tables(1)
    hgrid         => Get_grid(mem%grid_id)
    surface       => Get_vgrid('surface')
    soil_e        => Get_vgrid('soil_depth_energy')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    vgrid_snow_e  => Get_vgrid('snow_depth_energy')
    IF( dsl4jsb_Config(SSE_)%l_snow) THEN
      snow_e  => Get_vgrid('snow_depth_energy')
    END IF

    CALL mem%Add_var( 'vol_heat_cap', mem%vol_heat_cap,                 &
      & hgrid, surface,                                                      &
      & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'vol_heat_cap_sl', mem%vol_heat_cap_sl,           &
      & hgrid, soil_e,                                                       &
      & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'vol_heat_cap_snow', mem%vol_heat_cap_snow,       &
      & hgrid, surface,                                                      &
      & t_cf('vol_heat_capacity_snow', 'J m-3 K-1', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'heat_cond', mem%heat_cond,                       &
      & hgrid, surface,                                                      &
      & t_cf('heat_cond_soil', 'J m-1 s-1 K-1', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'heat_cond_sl', mem%heat_cond_sl,                 &
      & hgrid, soil_e,                                                       &
      & t_cf('heat_cond_soil', 'J m-1 s-1 K-1', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'heat_cond_snow', mem%heat_cond_snow,             &
      & hgrid, surface,                                                      &
      & t_cf('heat_cond_snow', 'J m-1 s-1 K-1', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    IF (      .NOT. dsl4jsb_Config(SSE_)%l_heat_cond_map                     &
      & .AND. .NOT. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      CALL mem%Add_var( 'thermal_diffusivity', mem%thermal_diffusivity,        &
        & hgrid, surface,                                                      &
        & t_cf('thermal_diffusivity_soil', 'm2 s-1', ''),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )
    END IF

    CALL mem%Add_var( 't_soil_sl', mem%t_soil_sl,                            &
      & hgrid, soil_e,                                                       &
      & t_cf('soil_temperature', 'K', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=280.0_wp )

    CALL mem%Add_var( 't_soil_acoef', mem%t_soil_acoef,                      &
      & hgrid, soil_e,                                                       &
      & t_cf('soil_temperature_acoef', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 't_soil_bcoef', mem%t_soil_bcoef,                      &
      & hgrid, soil_e,                                                       &
      & t_cf('soil_temperature_bcoef', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    ! snow
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      CALL mem%Add_var( 'snow_depth_sl', mem%snow_depth_sl,                    &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_depth_sl', 'm', 'Depth of snow layers'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.TRUE.,                                                     &
        & output_level=MEDIUM,                                                 &
        & initval_r=0._wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 't_snow', mem%t_snow,                                  &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_temperature', 'K', ''),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=BASIC,                                                  &
        & initval_r=273.15_wp )

      CALL mem%Add_var( 't_snow_acoef', mem%t_snow_acoef,                      &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_temperature_acoef', '', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 't_snow_bcoef', mem%t_snow_bcoef,                      &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_temperature_bcoef', '', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )
    END IF

    CALL mem%Add_var( 'hcap_grnd', mem%hcap_grnd,                            &
      & hgrid, surface,                                                      &
      & t_cf('heat_capacity_ground', 'J m-2 K-1', 'Ground heat capacity'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=1.3e5_wp )

    CALL mem%Add_var( 'hcap_grnd_old', mem%hcap_grnd_old,                            &
      & hgrid, surface,                                                              &
      & t_cf('heat_capacity_ground_old', 'J m-2 K-1', 'Ground heat capacity (old)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC,                                                          &
      & initval_r=1.3e5_wp )

    CALL mem%Add_var( 'grnd_hflx', mem%grnd_hflx,                            &
      & hgrid, surface,                                                      &
      & t_cf('grnd_hflx', 'J m-2 s-1', 'Ground heat flux'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=10.0_wp )

    CALL mem%Add_var( 'grnd_hflx_old', mem%grnd_hflx_old,                    &
      & hgrid, surface,                                                      &
      & t_cf('grnd_hflx_old', 'J m-2 s-1', 'Ground heat flux (old)'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=10.0_wp )

    CALL mem%Add_var( 'thaw_depth', mem%thaw_depth,                          &
      & hgrid, surface,                                                      &
      & t_cf('thaw_depth', 'm', 'Thawing depth'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=MEDIUM, &
      & initval_r=10.0_wp, l_aggregate_all=.TRUE. )
#endif

  END SUBROUTINE Init_sse_memory

#endif
END MODULE mo_sse_memory_class
