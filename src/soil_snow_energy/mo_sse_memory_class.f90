!>
!! @brief Contains the memory class for the soil and snow energy process.
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

MODULE mo_sse_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes SSE_

  ! Use process configurations
  dsl4jsb_Use_config(SSE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  TYPE, EXTENDS(t_jsb_memory) :: t_sse_memory
    !
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
      & t_soil,               & !< Temperature of soil layers [K]
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

    !! quincy
    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth      , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    root_depth      , & !< actual rooting depth [m]
                                    t_srf_old       , & !< previous surface temperature [K]
                                    zril_old            !< previous timestep's Reynold's number

    ! JSBACH SOIL
    TYPE(t_jsb_var_real3d)          :: &
                                    soil_depth_sl, &            !< Width ! of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_width_sl, &        !< width of each soil layer that can be saturated with water, identical with soil_depth_sl (until bedrock) [m]
                                    soil_lay_depth_center_sl, & !< depth at the center of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_ubound_sl, & !< depth at the upper bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_lbound_sl    !< depth at the lower bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]

    ! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d)          :: &
                                    interception,     & !< surface interception evaporation  [kg m-2 s-1]
                                    pet,              & !< potential evaporation  [kg m-2 s-1]
                                    evaporation,      & !< surface evaporation  [kg m-2 s-1]
                                    runoff,           & !< surface runoff  [kg m-2 s-1]
                                    drainage,         & !< drainage        [kg m-2 s-1]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage [1 s-1] 

    ! Additional variables for sse
    TYPE(t_jsb_var_real2d)          :: &
                                    w_skin,           & !< Water content in skin reservoir                                          [m]
                                    w_soil_root,      & !< Water content in root zone of the soil                                   [m]
                                    w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
                                    w_soil_root_pwp,  & !< Water content at permanent wilting point in root zone of the soil        [m]
                                    t_soil_root         !< Temperature of fine roots in the soil                                    [K]

    TYPE(t_jsb_var_real3d)          :: &
                                    percolation_sl,   & !< fraction of mass transferred to below layer                              [fraction s-1]
                                    drainage_sl,      & !< drainage of soil layers                                                  [kg m-2 s-1]
                                    w_soil_sl,        & !< Water content in soil layers                                             [m]
                                    w_soil_sat_sl,    & !< Water content of soil layers at saturation 
                                                        !! (derived from volumetric capacity at saturation)                         [m]
                                    w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                                                        !! (derived from volumetric field capacity)                                 [m]
                                    w_soil_pwp_sl,    & !< Water content of soil layers at permanent wilting point
                                                        !! (derived from vol. perm wilt point)                                      [m]
                                    saxtonA,          & !< coefficient in moisture-tension relationship
                                    saxtonB,          & !< coefficient in moisture-tension relationship
                                    saxtonC,          & !< exponent of moisture-conductivity relationship
                                    kdiff_sat_sl        !< saturated hydraulic conductivity                                         [m s-1]

    ! heat fluxes
    TYPE(t_jsb_var_real2d)          :: &                              
                                    sensible_heat_flx     , &  !> sensible heat flux [W m-2]
                                    latent_heat_flx       , &  !> latent heat flux [W m-2]
                                    ground_heat_flx       , &  !> ground heat flux [W m-2]
                                    ground_heat_flx_old        !> ground heat flux of previous timestep [W m-2]

    ! transpiration
    TYPE(t_jsb_var_real2d)          :: &
                                    transpiration       !< transpiration in [kg m-2 s-1]

    !Soil Water addons for SLM in HYDRO
    TYPE(t_jsb_var_real2d)          :: &
                                    w_soil_root_theta, & !< root zone integrated plant available water, fraction of maximum 
                                    w_soil_root_pot      !< soil water potential [MPa]

    TYPE(t_jsb_var_real3d)          :: &
                                    w_soil_pot_sl        !< soil water potential per layer [MPa]

    ! soil energy 
    TYPE(t_jsb_var_real3d)          :: &
                                    t_soil_sl            !< soil temperature [K] 

    ! soil physical properties (not sure where they should belong to, but they are needed in the long-term)
    TYPE(t_jsb_var_real3d)          :: & 
                                    heat_capa_sl,         &  !< soil heat capacity [J kg-1]
                                    bulk_dens_sl,         &  !< soil bulk density [kg m-3]
                                    therm_cond_sl            !< soil thermal conductivity [W m-2 K-1]
    TYPE(t_jsb_var_real3d)          :: & 
                                    sand_sl             , &  !< soil sand proportion of mineral soil
                                    silt_sl             , &  !< soil silt proportion of mineral soil
                                    clay_sl             , &  !< soil clay proportion of mineral soil(re-calculated as: clay = 1.0 -sand -silt)
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

  SUBROUTINE Init_sse_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM !, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_sse_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: model_id

    dsl4jsb_Def_config(SSE_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_e, snow_e      ! Vertical grids
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_memory'

    IF (SIZE(lct_ids) < 1) CALL finish(routine, 'lcd_ids not specified')

    model => Get_model(model_id)

    dsl4jsb_Get_Config(SSE_)

    table = tables(1)

    hgrid => Get_grid(mem%grid_id)

    surface       => Get_vgrid('surface')
    soil_e        => Get_vgrid('soil_depth_energy')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    IF( dsl4jsb_Config(SSE_)%l_snow) THEN
      snow_e  => Get_vgrid('snow_depth_energy')
    END IF

    !! !$acc enter data pcreate(mem)

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

    CALL mem%Add_var( 'thermal_diffusivity', mem%thermal_diffusivity,        &
      & hgrid, surface,                                                      &
      & t_cf('thermal_diffusivity_soil', 'm2 s-1', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 't_soil', mem%t_soil,                                  &
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
      & t_cf('grnd_hflx', '', 'Ground heat flux'),                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=10.0_wp )

    CALL mem%Add_var( 'grnd_hflx_old', mem%grnd_hflx_old,                    &
      & hgrid, surface,                                                      &
      & t_cf('grnd_hflx_old', '', 'Ground heat flux (old)'),                 &
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

    !!---------------------------------------------------------------------------------------------
    !! quincy
    IF (model%config%use_quincy) THEN
      IF ( (     One_of(LAND_TYPE, lct_ids(:)) > 0 &
        &   .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0)) THEN
        
        CALL mem%Add_var(mem, 'soil_depth', mem%soil_depth, &
              hgrid, surface, &
              t_cf('soil_depth', 'm', 'Soil depth derived from textures (bedrock), for hydrology'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'root_depth', mem%root_depth, &
              hgrid, surface, &
              t_cf('root_depth', 'm', 'actual rooting depth'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 't_srf_old', mem%t_srf_old, &
              hgrid, surface, &
              t_cf('t_srf_old', 'K', 'previous surface temperature'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'zril_old', mem%zril_old, &
              hgrid, surface, &
              t_cf('zril_old', '', 'previous timesteps Reynolds number'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'soil_depth_sl', mem%soil_depth_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('soil_depth_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'soil_lay_width_sl', mem%soil_lay_width_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('soil_lay_width_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'soil_lay_depth_center_sl', mem%soil_lay_depth_center_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('soil_lay_depth_center_sl', 'm', 'Depth at the center of each soil layer, until bedrock'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'soil_lay_depth_ubound_sl', mem%soil_lay_depth_ubound_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('soil_lay_depth_ubound_sl', 'm', 'Depth at the upper boundary of each soil layer, until bedrock'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'soil_lay_depth_lbound_sl', mem%soil_lay_depth_lbound_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('soil_lay_depth_lbound_sl', 'm', 'Depth at the lower boundary of each soil layer, until bedrock'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'interception', mem%interception, &
              hgrid, surface, &
              t_cf('interception', 'kg m-2 s-1', 'surface interception evaporation'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'pet', mem%pet, &
              hgrid, surface, &
              t_cf('pet', 'kg m-2 s-1', 'potential evaporation'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'evaporation', mem%evaporation, &
              hgrid, surface, &
              t_cf('evaporation', 'kg m-2 s-1', 'surface evaporation'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'runoff', mem%runoff, &
              hgrid, surface, &
              t_cf('runoff', 'kg m-2 s-1', 'surface runoff'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'drainage', mem%drainage, &
              hgrid, surface, &
              t_cf('drainage', 'kg m-2 s-1', 'drainage'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'drainage_fraction', mem%drainage_fraction, &
              hgrid, surface, &
              t_cf('drainage_fraction', '1 s-1', 'water mass fraction lost from the soil column by drainage'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_skin', mem%w_skin, &
              hgrid, surface, &
              t_cf('w_skin', 'm', 'Water content in skin reservoir'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_root', mem%w_soil_root, &
              hgrid, surface, &
              t_cf('w_soil_root', 'm', 'Water content in root zone of the soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_root_fc', mem%w_soil_root_fc, &
              hgrid, surface, &
              t_cf('w_soil_root_fc', 'm', 'Water content at field capacity in root zone of the soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_root_pwp', mem%w_soil_root_pwp, &
              hgrid, surface, &
              t_cf('w_soil_root_pwp', 'm', 'Water content at permanent wilting point in root zone of the soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 't_soil_root', mem%t_soil_root, &
              hgrid, surface, &
              t_cf('t_soil_root', 'K', 'Temperature of fine roots in the soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'percolation_sl', mem%percolation_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('percolation_sl', 'fraction s-1', 'fraction of mass transferred to below layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'drainage_sl', mem%drainage_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('drainage_sl', 'kg m-2 s-1', 'drainage of soil layers'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_sl', mem%w_soil_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('w_soil_sl', 'm', 'Water content in soil layers'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_sat_sl', mem%w_soil_sat_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('w_soil_sat_sl', 'm', 'Water content of soil layers at saturation'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_fc_sl', mem%w_soil_fc_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('w_soil_fc_sl', 'm', 'Water content of soil layers at field capacity'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var(mem, 'w_soil_pwp_sl', mem%w_soil_pwp_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('w_soil_pwp_sl', 'm', 'Water content of soil layers at permanent wilting point'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'saxtonA', mem%saxtonA, &
              hgrid, vgrid_soil_sb, &
              t_cf('saxtonA', '', 'coefficient in moisture-tension relationship'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'saxtonB', mem%saxtonB, &
              hgrid, vgrid_soil_sb, &
              t_cf('saxtonB', '', 'coefficient in moisture-tension relationship'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'saxtonC', mem%saxtonC, &
              hgrid, vgrid_soil_sb, &
              t_cf('saxtonC', '', 'exponent of moisture-conductivity relationship'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'kdiff_sat_sl', mem%kdiff_sat_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('kdiff_sat_sl', 'm s-1', 'saturated hydraulic conductivity'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var( mem, 'sensible_heat_flx', mem%sensible_heat_flx, &
              hgrid, surface, &
              t_cf('sensible_heat_flx', 'W m-2', 'sensible heat flux'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp )
  
        CALL mem%Add_var( mem, 'latent_heat_flx', mem%latent_heat_flx, &
              hgrid, surface, &
              t_cf('latent_heat_flx', 'W m-2', 'latent heat flux'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp )
  
        CALL mem%Add_var( mem, 'ground_heat_flx', mem%ground_heat_flx, &
              hgrid, surface, &
              t_cf('ground_heat_flx', 'W m-2', 'ground heat flux'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp )
  
        CALL mem%Add_var( mem, 'ground_heat_flx_old', mem%ground_heat_flx_old, &
              hgrid, surface, &
              t_cf('ground_heat_flx_old', 'W m-2', 'ground heat flux of previous timestep'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp )
  
        CALL mem%Add_var(mem, 'transpiration', mem%transpiration, &
              hgrid, surface, &
              t_cf('transpiration', 'kg m-2 s-1', 'transpiration in'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'w_soil_root_theta', mem%w_soil_root_theta, &
              hgrid, surface, &
              t_cf('w_soil_root_theta', '', 'root zone integrated plant available water, fraction of maximum '), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=1.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'w_soil_root_pot', mem%w_soil_root_pot, &
              hgrid, surface, &
              t_cf('w_soil_root_pot', 'MPa', 'soil water potential'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'w_soil_pot_sl', mem%w_soil_pot_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('w_soil_pot_sl', 'MPa', 'soil water potential per layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 't_soil_sl', mem%t_soil_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('t_soil_sl', 'K', 'soil temperature'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=273.15_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'heat_capa_sl', mem%heat_capa_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('heat_capa_sl', 'J kg-1', 'soil heat capacity'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'bulk_dens_sl', mem%bulk_dens_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('bulk_dens_sl', 'kg m-3', 'soil bulk density'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'therm_cond_sl', mem%therm_cond_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('therm_cond_sl', 'W m-2 K-1', 'soil thermal conductivity'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'sand_sl', mem%sand_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('sand_sl', '', 'soil sand proportion'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
              
        CALL mem%Add_var(mem, 'silt_sl', mem%silt_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('silt_sl', '', 'soil silt proportion'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
              
        CALL mem%Add_var(mem, 'clay_sl', mem%clay_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('clay_sl', '', 'soil clay proportion'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'volume_min_sl', mem%volume_min_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('volume_min_sl', '', 'soil mineral volumetric proportion'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'qmax_org_min_sl', mem%qmax_org_min_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('qmax_org_min_sl', 'mol C kg-1 mineral soil', 'maximum OM sorption capacity of mineral soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'qmax_po4_min_sl', mem%qmax_po4_min_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('qmax_po4_min_sl', 'mol P kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'qmax_nh4_min_sl', mem%qmax_nh4_min_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('qmax_nh4_min_sl', 'mol N kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var(mem, 'qmax_po4_om_sl', mem%qmax_po4_om_sl, &
              hgrid, vgrid_soil_sb, &
              t_cf('qmax_po4_om_sl', 'mol P kg-1 SOM', 'maximum po4 sorption capacity of SOM'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              lrestart=.FALSE., &
              loutput=.FALSE., &
              initval_r=0.0_wp) ! initvalue

      ENDIF ! tile: land/veg
    ENDIF   ! quincy

    END SUBROUTINE Init_sse_memory

#endif
END MODULE mo_sse_memory_class
