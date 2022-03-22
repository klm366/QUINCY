!>
!! @brief Contains the memory class for the hydrology process.
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

MODULE mo_hydro_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE, LAKE_TYPE, GLACIER_TYPE
  USE mo_jsb_varlist,      ONLY: BASIC !, MEDIUM, FULL
  USE mo_jsb_process_class, ONLY: ON_TILE_, AGGREGATE_

  ! Use of prcesses in this module
  dsl4jsb_Use_processes HYDRO_, SSE_, SEB_, ASSIMI_

  ! Use process configurations
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(ASSIMI_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 90

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_hydro_memory

    TYPE(t_jsb_var_real2d) :: &
      & fract_snow,           & !< Snow fraction of tile                                                []
      & w_snow,               & !< Water content of snow reservoir                                      [m water equivalent]
      & evapo_snow,           & !< Evaporation from snow                                                [kg m-2 s-1]
      & snowmelt,             & !< Snow melt                                                            [kg m-2 s-1]
      & evapopot,             & !< Potential evaporation                                                [kg m-2 s-1]
      & evapotrans,           & !< Evapotranspiration independent of the type of tile                   [kg m-2 s-1]
                                !< e.g. for the lake tile it is evapopot and for the land tile it is evapotrans_lnd.
      & water_flux,           & !< All land water fluxes; needed for water budget check                 [m3/s]
      & w_land,               & !< Total land water content                                             [m3]
      & water_budget            !< Land water balance error within a time step                          [m3]

    TYPE(t_jsb_var_real2d) :: &
      & runoff,               & !< surface runoff
      & drainage,             & !< drainage
      & discharge,            & !< discharge (local)
      & discharge_ocean,      & !< discharge to the ocean
      & internal_drain          !< internal drain

    TYPE(t_jsb_var_real2d) :: &
      & elevation,            & !< Elevation
      & oro_stddev              !< Standard deviation of orography

    ! Additional variables for land type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_soil,      & !< Snow fraction on soil                                               []
      & snow_soil,            & !< Water content of snow reservoir of soil                             [m water equivalent]
      & snow_soil_dens,       & !< Density of snow on soil                                             [kg m-3]
      & evapotrans_lnd,       & !< Evapotranspiration over vegetation and bare soil type               [kg m-2 s-1]
      & q_snocpymlt             !< Heating by melting of snow on canopy                                [W m-2]

    ! Additional variables for soil

    ! Parameters
    TYPE(t_jsb_var_real3d) :: &
      & soil_depth_sl,        & !< Depth of each soil layer that can be saturated with water (until bedrock) [m]
      & fract_org_sl            !< Fractions of organic material in soil layers                        []
    TYPE(t_jsb_var_real2d) :: &
      & soil_depth,           & !< Soil depth derived from textures (bedrock)                          [m]
      & rootzone_max ,        & !< Maximum root zone moisture content                                  [m]
      & vol_field_cap,        & !< Volumetric soil field capacity                                      [m/m]
      & vol_p_wilt,           & !< Volumetric permanent wilting point                                  [m/m]
      & vol_porosity,         & !< Volumetric porosity of mineral soil                                 [m/m]
      & pore_size_index,      & !< Soil pore size distribution index                                   []
      & bclapp,               & !< Exponent B in Clapp and Hornberger
      & matrix_pot,           & !< Matrix potential [m]
      & hyd_cond_sat            !< Saturated hydraulic conductivity [m/s]

    TYPE(t_jsb_var_real2d) :: &
      & fract_water,          & !< Wet (skin reservoir) fraction of tile (soil and canopy)                  []
      & w_skin,               & !< Water content in skin reservoir (soil and canopy)                        [m]
      & snow_accum,           & !< Snow accumulation at non-glacier/non-lake land points                    [m water equivalent]
      & evapotrans_soil,      & !< Evapotranspiration from soil w/o snow and skin evaporation               [kg m-2 s-1]
      !& evapo_soil,           & !< Evaporation from ground                                                  [kg m-2 s-1]
      & evapo_skin,           & !< Evaporation from skin reservoir                                          [kg m-2 s-1]
      & evapo_deficit,        & !< Evaporation deficit flux due to inconsistent treatment of snow evap.     [m water equivalent]
      & water_excess,         & !< Water available for infiltration into the soil                           [m water equivalent]
      & overflow,             & !< Water content of reservoir for soil water overflow (terraplanet)         [m]
      & wat_rootzone,         & !< Liquid water content in the rootzone                                     [m]
      & wat_rootzone_rel        !< Water content of the rootzone relative to the maximum possible amout     [1]

    TYPE(t_jsb_var_real3d) :: &
      & vol_field_cap_sl,     & !< Volumetric soil field capacity                                           [m/m]
      & vol_p_wilt_sl,        & !< Volumetric permanent wilting point                                       [m/m]
      & vol_porosity_sl,      & !< Volumetric porosity of mineral                                           [m/m]
      & pore_size_index_sl,   & !< Soil pore size distribution index                                        []
      & bclapp_sl,            & !< Exponent B in Clapp and Hornberger
      & matrix_pot_sl,        & !< Matrix potential                                                         [m]
      & hyd_cond_sat_sl,      & !< Saturated hydraulic conductivity                                         [m/s]
      & soil_wat_sl,          & !< Water content in soil layers                                             [m]
      & soil_ice_sl,          & !< Ice content in soil layers                                               [m]
      & wat_freeze_sl,        & !< Flux from freezing water in soil layers                                  [kg m-2 s-1]
      & ice_melt_sl,          & !< Flux from melting ice in soil layers                                     [kg m-2 s-1]
      & w_soil_sat_sl,        & !< Water content of soil layers at saturation (derived from soil porosity)  [m]
      & w_soil_fc_sl,         & !< Water content of soil layers at field capacity
                                !! (derived from volumetric field capacity)                                 [m]
      & w_soil_pwp_sl           !< Water content of soil layers at permanent wilting point
                                !! (derived from vol. perm wilt point)                                      [m]

    ! Additional variables for PFT lct_type

    ! Parameters
    TYPE(t_jsb_var_real3d) :: &
      & root_depth_sl           !< Rooted depth per soil layer (until rooting depth) [m]
    TYPE(t_jsb_var_real2d) :: &
      & root_depth,           & !< Rooting depth [m]
      & fract_snow_can,       & !< Snow fraction on canopy                                                  []
      & snow_canopy,          & !< Water content of snow reservoir on canopy                                [m water equivalent]
      & rootzone_fc,          & !< Water content at field capacity in root zone of the soil                 [m]
      & rootzone_pwp,         & !< Water content at permanent wilting point in root zone of the soil        [m]
      & water_stress,         & !< Water stress factor of canopy                        (1: no water stress, 0: infinite stress)
      & canopy_cond_unlimited,& !< Canopy conductance without water limit  [m/s] (mean value)
      & canopy_cond_limited,  & !< Canopy conductance with water limit [m/s] (mean value)
      & transpiration           !< Transpiration

    ! Additional variables for GLACIER lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_glac,      & !< Snow fraction on glacier                                                 []
      & w_glac,               & !< Glacier depth including snow on glacier                                  [m water equivalent]
      & runoff_glac             !< Runoff from glacier (rain+melt, no calving)                              [m water equivalent]

    ! Additional variables if no separate glacier lct_type is used and glaciers are treated as part of SOIL lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_glac          !< Glacier fraction of tile                                                 []

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & evapo_wtr,            & !< Evaporation over lake water
      & evapo_ice,            & !< Evaporation over lake ice
      & fract_snow_lice,      & !< Fraction of swow on lake ice                                             []
      & w_snow_lice             !< Water content of snow on lake ice                                        [m water equivalent]

    !! quincy
    !Soil Water addons for QUINCY in HYDRO
    TYPE(t_jsb_var_real2d) :: &
      & w_soil_root_theta,    & !< root zone integrated plant available water, fraction of maximum 
      & w_soil_root_pot         !< soil water potential [MPa]

    ! Diagnostic global land means/sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) :: &
      trans_gmean,            & !< Global mean transpiration                [kg m-2 s-1]
      evapotrans_gmean,       & !< Global land mean evapotranspiration      [kg m-2 s-1]
      w_land_gsum,            & !< Global land water content                [km3]
      discharge_ocean_gsum,   & !< Global water discharge to the oceans     [Sv]
      wat_rootzone_rel_gmean, & !< Global mean relative rootzone moisture   [1]
      fract_snow_gsum,        & !< Global snow area on non-glacier land     [Mio m2]
      w_snow_gsum               !< Global snow amount on non-glacier land   [Gt]

  CONTAINS
    PROCEDURE :: Init => Init_hydro_memory
  END TYPE t_hydro_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_memory_class'

CONTAINS

  SUBROUTINE Init_hydro_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    & TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    USE mo_jsb_physical_constants, ONLY: &
      & dens_snow,                       & ! Density of snow (l_dynsnow=.FALSE.)
      & dens_snow_min                      ! Density of fresh snow (l_dynsnow=.TRUE.)

    CLASS(t_hydro_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)    :: prefix
    CHARACTER(len=*),      INTENT(in)    :: suffix
    INTEGER,               INTENT(in)    :: lct_ids(:)
    INTEGER,               INTENT(in)    :: model_id

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_config(ASSIMI_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid            ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w  ! Vertical grids
    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    REAL(wp) :: ini_snow_dens

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hydro_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    soil_w       => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(ASSIMI_)

    IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
      ini_snow_dens = dens_snow_min
    ELSE
      ini_snow_dens = dens_snow
    END IF

    !! !$acc enter data pcreate(mem)

    ! Common variables

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,202, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF

    CALL mem%Add_var('fract_snow', mem%fract_snow,                             &
      & hgrid, surface,                                                        &
      & t_cf('fract_snow', '-', 'Snow area fraction'),                         &
      & t_grib1(table, 255, grib_bits), grib2_desc,                            &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,212, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF
    CALL mem%Add_var( 'w_snow', mem%w_snow,                                    &
      & hgrid, surface,                                                        &
      & t_cf('w_snow', 'm (water equivalent)', 'Snow amount'),                 &
      & t_grib1(table, 255, grib_bits), grib2_desc,                            &
      & prefix, suffix,                                                        &
      & output_level=BASIC, initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Total evapotranspiration, including sublimation
    CALL mem%Add_var( 'evapotrans', mem%evapotrans,                                      &
      & hgrid, surface,                                                                  &
      & t_cf('evapotrans', 'kg m-2 s-1', 'Evapotranspiration, including sublimation'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Potential evapotranspiration/sublimation
    CALL mem%Add_var( 'evapopot', mem%evapopot,                                          &
      & hgrid, surface,                                                                  &
      & t_cf('evapopot', 'kg m-2 s-1', 'Potential evaporation (including sublimation)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapo_snow', mem%evapo_snow,                                  &
      & hgrid, surface,                                                              &
      & t_cf('evapo_snow', 'kg m-2 s-1', 'Evaporation from snow'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'snowmelt', mem%snowmelt,                                      &
      & hgrid, surface,                                                              &
      & t_cf('snowmelt', 'kg m-2 s-1', 'Snow melt'),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'runoff', mem%runoff,                                          &
      & hgrid, surface,                                                              &
      & t_cf('runoff', 'kg m-2 s-1', 'Surface runoff'),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & loutput=.TRUE., output_level=BASIC,                                          &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'drainage', mem%drainage,                                      &
      & hgrid, surface,                                                              &
      & t_cf('drainage', 'kg m-2 s-1', 'Drainage'),                                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & loutput=.TRUE., output_level=BASIC,                                          &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'discharge', mem%discharge,                                    &
      & hgrid, surface,                                                              &
      & t_cf('discharge', 'kg m-2 s-1', 'Local discharge'),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'discharge_ocean', mem%discharge_ocean,                        &
      & hgrid, surface,                                                              &
      & t_cf('discharge_ocean', 'm3 s-1', 'Discharge to the ocean'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'internal_drain', mem%internal_drain,                          &
      & hgrid, surface,                                                              &
      & t_cf('internal_drain', 'kg m-2 s-1', 'Internal drainage'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'water_flux', mem%water_flux,                                  &
      & hgrid, surface,                                                              &
      & t_cf('water_flux', 'm3 s-1',                                                 &
      &      'All land water fluxes (for water balance check)'),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'w_land', mem%w_land,                                          &
      & hgrid, surface,                                                              &
      & t_cf('w_land', 'm3', 'Total land water content'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget', mem%water_budget,                              &
      & hgrid, surface,                                                              &
      & t_cf('water_budget', 'm3', 'Land water balance error'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'elevation', mem%elevation,                                    &
      & hgrid, surface,                                                              &
      & t_cf('elevation', 'm', 'Mean orography'),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'oro_stddev', mem%oro_stddev,                                  &
      & hgrid, surface,                                                              &
      & t_cf('oro_stddev', 'm', 'Orographic standard deviation'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    ! Additional variables for land type
    IF ( (     One_of(VEG_TYPE,     lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE,    lct_ids(:)) > 0 &
      &   .OR. One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'fract_snow_soil', mem%fract_snow_soil,                      &
        & hgrid, surface,                                                            &
        & t_cf('fract_snow_soil', '-', 'Snow fraction on soil'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_soil', mem%snow_soil,                                  &
        & hgrid, surface,                                                            &
        & t_cf('snow_soil', 'm (water equivalent)', 'Snow amount on soil'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_soil_dens', mem%snow_soil_dens,                        &
        & hgrid, surface,                                                            &
        & t_cf('snow_soil_dens', 'kg m-3', 'Density of snow on soil'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=ini_snow_dens, l_aggregate_all=.TRUE. )

      ! Total evapotranspiration, including sublimation
      CALL mem%Add_var( 'evapotrans_lnd', mem%evapotrans_lnd,                    &
        & hgrid, surface,                                                        &
        & t_cf('evapotrans_lnd', 'kg m-2 s-1',                                   &
        &      'Evapotranspiration from land surface w/o lakes'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE.,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'q_snocpymlt', mem%q_snocpymlt,                          &
        & hgrid, surface,                                                        &
        & t_cf('q_snocpymlt', 'W m-2', 'Heating due to snow melt on canopy'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE.,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,201, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'fract_water', mem%fract_water,                          &
        & hgrid, surface,                                                        &
        & t_cf('fract_water', '-', 'Wet surface fraction'),                      &
        & t_grib1(table, 255, grib_bits), grib2_desc,                            &
        & prefix, suffix,                                                        &
        & output_level=BASIC,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for soil
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'soil_depth', mem%soil_depth,                          &
        & hgrid, surface,                                                      &
        & t_cf('soil_depth', 'm', 'Depth of soil until bedrock'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'soil_depth_sl', mem%soil_depth_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('soil_depth_sl', 'm', 'Soil depth within each soil layer'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
        CALL mem%Add_var( 'fract_org_sl', mem%fract_org_sl,                          &
          & hgrid, soil_w,                                                           &
          & t_cf('fract_org_sl', '-', 'Fraction of organic material in soil layer'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
          & prefix, suffix,                                                          &
          & lrestart=.TRUE.,                                                         &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

      CALL mem%Add_var( 'rootzone_max', mem%rootzone_max,                      &
        & hgrid, surface,                                                      &
        & t_cf('rootzone_max', 'm', 'Maximum root zone moisture'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_field_cap', mem%vol_field_cap,                    &
        & hgrid, surface,                                                      &
        & t_cf('vol_field_capacity', 'm/m', 'Volumetric soil field capacity'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_field_cap_sl', mem%vol_field_cap_sl,              &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_field_capy_sl', 'm/m', 'Volumetric soil field capacity'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_p_wilt', mem%vol_p_wilt,                          &
        & hgrid, surface,                                                      &
        & t_cf('vol_p_wilt', 'm/m', 'Volumetric permanent wilting point'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_p_wilt_sl', mem%vol_p_wilt_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_p_wilt_sl', 'm/m', 'Volumetric permanent wilting point'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_porosity', mem%vol_porosity,                      &
        & hgrid, surface,                                                      &
        & t_cf('vol_porosity', 'm/m', 'Volumetric soil porosity'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_porosity_sl', mem%vol_porosity_sl,                &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_porosity_sl', 'm/m', 'Volumetric soil porosity'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'pore_size_index', mem%pore_size_index,                &
        & hgrid, surface,                                                      &
        & t_cf('pore_size_index', '', 'Soil pore size distribution index'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'pore_size_index_sl', mem%pore_size_index_sl,          &
        & hgrid, soil_w,                                                       &
        & t_cf('pore_size_index_sl', '', 'Soil pore size distribution index'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'bclapp', mem%bclapp,                                  &
        & hgrid, surface,                                                      &
        & t_cf('bclapp_sl', '', 'Clapp and Hornberger (1978) exponent b'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'bclapp_sl', mem%bclapp_sl,                            &
        & hgrid, soil_w,                                                       &
        & t_cf('bclapp_sl', '', 'Clapp and Hornberger (1978) exponent b'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'matrix_pot', mem%matrix_pot,                          &
        & hgrid, surface,                                                      &
        & t_cf('matrix_pot', 'm', 'Soil saturated matrix potential'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'matrix_pot_sl', mem%matrix_pot_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('matrix_pot_sl', 'm', 'Soil saturated matrix potential'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'hyd_cond_sat', mem%hyd_cond_sat,                           &
        & hgrid, surface,                                                           &
        & t_cf('hyd_cond_sat', 'm s-1', 'Hydraulic conductivity at saturation'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE.,                                                         &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'hyd_cond_sat_sl', mem%hyd_cond_sat_sl,                     &
        & hgrid, soil_w,                                                            &
        & t_cf('hyd_cond_sat_sl', 'm s-1', 'Hydraulic conductivity at saturation'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE.,                                                         &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,211, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'w_skin', mem%w_skin,                                       &
        & hgrid, surface,                                                           &
        & t_cf('w_skin', 'm', 'Water content of surface skin reservoir'),           &
        & t_grib1(table, 255, grib_bits), grib2_desc,                               &
        & prefix, suffix,                                                           &
        & output_level=BASIC,                                                       &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_accum', mem%snow_accum,                               &
        & hgrid, surface,                                                           &
        & t_cf('snow_accum', 'm (water equivalent)',                                &
        &      'Snow budget change within the time step'),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    ! END IF

    ! IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
    !   &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
    !   &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
    !   &  )                                       &
    !   &  .AND. .NOT. ASSOCIATED(mem%fract_water%ptr) ) THEN

      CALL mem%Add_var( 'evapotrans_soil', mem%evapotrans_soil,                &
        & hgrid, surface,                                                      &
        & t_cf('evapotrans_soil', 'kg m-2 s-1',                                &
        &      'Evapotranspiration from soil (w/o snow and skin res.'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'evapo_skin', mem%evapo_skin,                          &
        & hgrid, surface,                                                      &
        & t_cf('evapo_skin', 'kg m-2 s-1', 'Evaporation from skin reservoir'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      !Todo: is this variable of any interest?
      CALL mem%Add_var( 'evapo_deficit', mem%evapo_deficit,                    &
        & hgrid, surface,                                                      &
        & t_cf('evapo_deficit', 'm',   'Water that would have been evaporated  &
              & additionally, if there had been more snow or surface water'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      !Todo: Misleading name: water_excess corresponds to water_to_soil of calc_soil_hydrology
      CALL mem%Add_var( 'water_excess', mem%water_excess,                      &
        & hgrid, surface,                                                      &
        & t_cf('water_excess', 'm', 'Amount of water entering the soil'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Overflow pool is only used in the terra planet setup (TPE closed).
      CALL mem%Add_var( 'overflow', mem%overflow,                              &
        & hgrid, surface,                                                      &
        & t_cf('overflow', 'm', 'Terra planet soil water overflow pool'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,213, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'wat_rootzone', mem%wat_rootzone,                      &
        & hgrid, surface,                                                      &
        & t_cf('wat_rootzone', 'm', 'Liquid water content of the root zone'),  &
        & t_grib1(table, 255, grib_bits), grib2_desc,                          &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'soil_wat_sl', mem%soil_wat_sl,                        &
        & hgrid, soil_w,                                                       &
        & t_cf('soil_wat_sl', 'm', 'Water content of soil layers'),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'soil_ice_sl', mem%soil_ice_sl,                        &
        & hgrid, soil_w,                                                       &
        & t_cf('soil_ice_sl', 'm', 'Ice content of soil layers'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wat_freeze_sl', mem%wat_freeze_sl,                              &
        & hgrid, soil_w,                                                                 &
        & t_cf('wat_freeze_sl', 'kg m-2 s-1', 'Freezing water flux in soil layers'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'ice_melt_sl', mem%ice_melt_sl,                                  &
        & hgrid, soil_w,                                                                 &
        & t_cf('ice_melt_sl', 'kg m-2 s-1', 'Melting ice flux in soil layers'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_sat_sl', mem%w_soil_sat_sl,                              &
        & hgrid, soil_w,                                                                 &
        & t_cf('w_soil_sat_sl', 'm', 'Water content in soil layers at saturation'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_fc_sl', mem%w_soil_fc_sl,                                &
        & hgrid, soil_w,                                                                 &
        & t_cf('w_soil_fc_sl', 'm', 'Water content in soil layers at field capacity'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_pwp_sl', mem%w_soil_pwp_sl,                              &
        & hgrid, soil_w,                                                                 &
        & t_cf('w_soil_pwp_sl', 'm',                                                     &
        &      'Water content in soil layers at permanent wilting point'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('wat_rootzone_rel', mem%wat_rootzone_rel,                         &
        &  hgrid, surface,                                                               &
        &  t_cf('wat_rootzone_rel', '-', 'Relative root zone moisture'),                 &
        &  t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & prefix, suffix,                                                                &
        & initval_r=0.0_wp )

    END IF

    ! Additional variables for vegetation
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'root_depth', mem%root_depth,                                    &
        & hgrid, surface,                                                                &
        & t_cf('root_depth', 'm', 'Rooting depth'),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'root_depth_sl', mem%root_depth_sl,                              &
        & hgrid, soil_w,                                                                 &
        & t_cf('root_depth_sl', 'm', 'Rooting depth within the soil layer'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE.,                                                                &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'fract_snow_can', mem%fract_snow_can,                            &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_can', '-', 'Snow fraction on canopy'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_canopy', mem%snow_canopy,                                  &
        & hgrid, surface,                                                                &
        & t_cf('snow_canopy', 'm (water equivalent)', 'Snow amount on canopy'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! todo: Currently only used with quincy, only in update_water_sress, where it could also be defined
      !    as local variable. Do we need that in the output without quincy? Do we need it at all?
      CALL mem%Add_var( 'rootzone_fc', mem%rootzone_fc,                                  &
        & hgrid, surface,                                                                &
        & t_cf('rootzone_fc', 'm', 'Root zone water content at field capacity'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & loutput=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! todo: Currently only used with quincy, only in update_water_sress, where it could also be defined
      !    as local variable. Do we need that in the output without quincy? Do we need it at all?
      CALL mem%Add_var( 'rootzone_pwp', mem%rootzone_pwp,                                &
        & hgrid, surface,                                                                &
        & t_cf('rootzone_pwp', 'm',                                                      &
        &      'Root zone water content at permanent wilting point'),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & loutput=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'water_stress', mem%water_stress,                                &
        & hgrid, surface,                                                                &
        & t_cf('water_stress', '-', 'Water stress factor'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE., lrestart=.FALSE., initval_r=0.0_wp )

      IF (.NOT. dsl4jsb_Config(ASSIMI_)%active) THEN
        CALL mem%Add_var( 'canopy_cond_unlimited', mem%canopy_cond_unlimited,            &
          & hgrid, surface,                                                              &
          & t_cf('canopy_cond_unlimited', 'm/s',                                         &
          &      'Canopy conductance ignoring water limitation'),                        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & loutput=.TRUE., initval_r=0.0_wp )

        CALL mem%Add_var( 'canopy_cond_limited', mem%canopy_cond_limited,                &
          & hgrid, surface,                                                              &
          & t_cf('canopy_cond_limited', 'm/s',                                           &
          &      'Canopy conductance accounting for water limitation'),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & loutput=.TRUE., initval_r=0.0_wp )
      END IF

      CALL mem%Add_var( 'transpiration', mem%transpiration,                              &
        & hgrid, surface,                                                                &
        & t_cf('transpiration', 'kg m-2 s-1', 'Transpiration'),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE.,                                                              &
        & loutput=.TRUE., output_level=BASIC,                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for GLACIER lct
    IF ( (     One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'fract_snow_glac', mem%fract_snow_glac,                          &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_glac', '-', 'Snow frction on glacier'),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_glac', mem%w_glac,                                            &
        & hgrid, surface,                                                                &
        & t_cf('w_glac', 'm (water equivalent)', 'Glacier depth (including snow)'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=BASIC,                                                            &
        & initval_r=20.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'runoff_glac', mem%runoff_glac,                                  &
        & hgrid, surface,                                                                &
        & t_cf('runoff_glac', 'kg m-2 s-1', 'Runoff from glacier (rain + melt)'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & initval_r=20.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for LAKE lct
    IF (      One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 'evapo_wtr', mem%evapo_wtr,                                &
        & hgrid, surface,                                                          &
        & t_cf('evapo_wtr', 'kg m-2 s-1', 'Evaporation from lake water'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & lrestart=.FALSE.,                                                        &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        CALL mem%Add_var( 'evapo_ice', mem%evapo_ice,                              &
          & hgrid, surface,                                                        &
          & t_cf('evapo_ice', 'kg m-2 s-1', 'Evaporation from lake ice'),          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & lrestart=.FALSE.,                                                      &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        !Fraction of snow on lake ice, rel. to ice fraction
        CALL mem%Add_var( 'fract_snow_lice', mem%fract_snow_lice,                  &
          & hgrid, surface,                                                        &
          & t_cf('fract_snow_lice', '-',                                           &
          &      'Snow fraction on lake ice (rel. to ice fraction)'),              &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & output_level=BASIC,                                                    &
          & initval_r=0.0_wp )

        CALL mem%Add_var( 'w_snow_lice', mem%w_snow_lice,                          &
          & hgrid, surface,                                                        &
          & t_cf('w_snow_lice', 'm water equivalent', 'Snow amount on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & output_level=BASIC,                                                    &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

    END IF

    !! quincy
    IF (model%config%use_quincy) THEN
      CALL mem%Add_var('w_soil_root_theta', mem%w_soil_root_theta,                       &
        & hgrid, surface,                                                                &
        & t_cf('w_soil_root_theta', '',                                                  &
        &      'root zone integrated plant available water, fraction of maximum '),      &
        & t_grib1(table, 255, grib_bits),                                                &
        & t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE., lrestart=.TRUE.,                                               &
        & initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('w_soil_root_pot', mem%w_soil_root_pot,                           &
        & hgrid, surface,                                                                &
        & t_cf('w_soil_root_pot', 'MPa', 'soil water potential'),                        &
        & t_grib1(table, 255, grib_bits),                                                &
        & t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                &
        & loutput=.TRUE., lrestart=.TRUE.,                                               &
        & initval_r=0.0_wp) ! initvalue
    ENDIF

#ifdef __ICON__
    ! Diagnostic 1d global land variables for experiment monitoring
    !       (1d stream variables are not supported with echam)
    !
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var('trans_gmean', mem%trans_gmean,                                   &
        & hgrid, surface,                                                                &
        & t_cf('trans_gmean', 'kg m-2 s-1', 'Global mean transpiration'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('evapotrans_gmean', mem%evapotrans_gmean,                         &
        & hgrid, surface,                                                                &
        & t_cf('evapotrans_gmean', 'kg m-2 s-1', 'Global land evapotranspiration'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('w_land_gsum', mem%w_land_gsum,                                   &
        & hgrid, surface,                                                                &
        & t_cf('w_land_gsum', 'km3', 'Global land water content'),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('discharge_ocean_gsum', mem%discharge_ocean_gsum,                 &
        & hgrid, surface,                                                                &
        & t_cf('discharge_ocean_gsum', 'Sv', 'Global water discharge to the oceans'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('wat_rootzone_rel_gmean', mem%wat_rootzone_rel_gmean,             &
        & hgrid, surface,                                                                &
        & t_cf('wat_rootzone_rel_gmean', '-', 'Global mean relative rootzone moisture'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('fract_snow_gsum', mem%fract_snow_gsum,                           &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_gsum', 'Mio km2', 'Global snow area (including glaciers)'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('w_snow_gsum', mem%w_snow_gsum,                                   &
        & hgrid, surface,                                                                &
        & t_cf('w_snow_gsum', 'Gt', 'Global snow amount on non-glacier land'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
    END IF
#endif

    END SUBROUTINE Init_hydro_memory
#endif
END MODULE mo_hydro_memory_class
